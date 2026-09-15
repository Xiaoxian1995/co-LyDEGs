# ==============================================================================
# Part 0: GMT 基因集文件检查与自动下载
# ==============================================================================

base_dir <- "./"
gmt_file_name <- "c2.cp.kegg_legacy.v2026.1.Hs.symbols.gmt"
gmt_path <- file.path(base_dir, gmt_file_name)

# 官方 MSigDB 下载链接（若链接更新，可手动替换 url 地址）
gmt_url <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/c2.cp.kegg_legacy.v2026.1.Hs.symbols.gmt"

if (!file.exists(gmt_path)) {
  message(">> 本地未检测到 GMT 基因集文件，正在尝试从 MSigDB 自动下载...")
  tryCatch({
    download.file(url = gmt_url, destfile = gmt_path, mode = "wb")
    message(">> GMT 文件下载成功！保存路径为：", gmt_path)
  }, error = function(e) {
    cat("\n⚠️ 自动下载失败！请检查网络连接或手动下载该文件：\n")
    cat("🔗 下载网址 1 (MSigDB 官方): https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp\n")
    cat("🔗 下载网址 2 (直接链接): ", gmt_url, "\n")
    cat("将下载好的文件命名为 '", gmt_file_name, "' 并存放在当前工作目录下：", normalizePath(base_dir), "\n\n")
  })
} else {
  message(">> 已检测到本地 GMT 文件：", gmt_path)
}


# ==============================================================================
# Part 1: 构建 GSE15222 全基因组表达矩阵 (full_gene_expr_matrix)
# ==============================================================================

library(tidyverse)
library(limma)
library(GEOquery)
library(AnnotationDbi)
library(org.Hs.eg.db)

# 1.1 读取 GEO 数据与提取表型矩阵
file_path <- file.path(base_dir, "GSE15222_series_matrix.txt.gz")
gse <- getGEO(filename = file_path, GSEMatrix = TRUE, getGPL = FALSE)

sample_info <- pData(gse)

# 基于 description 列提取样本分组 (Control vs Alzheimer/AD)
sample_info <- sample_info %>%
  dplyr::mutate(Group = case_when(
    str_detect(tolower(description), "control") ~ "Control",
    str_detect(tolower(description), "alzheimer") ~ "AD",
    TRUE ~ NA_character_
  ))

valid_samples <- !is.na(sample_info$Group)
sample_info <- sample_info[valid_samples, ]
expr_raw <- exprs(gse)[, valid_samples]

group_list <- factor(sample_info$Group, levels = c("Control", "AD"))

cat("==================================================\n")
cat("📊 样本分组统计：\n")
print(table(group_list))
cat("==================================================\n\n")

# 1.2 GPL2700 本地注释表生成 (探针 ID -> Gene Symbol)
gpl_file <- file.path(base_dir, "GPL2700.gz")
gpl <- GEOquery::getGEO(filename = gpl_file)
gpl_meta <- Table(gpl)

probe_to_acc <- gpl_meta[, c("ID", "GB_ACC")]
probe_to_acc <- probe_to_acc[probe_to_acc$GB_ACC != "" & !is.na(probe_to_acc$GB_ACC), ]
probe_to_acc$GB_ACC_clean <- sub("\\..*$", "", probe_to_acc$GB_ACC)

mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(probe_to_acc$GB_ACC_clean),
  columns = "SYMBOL",
  keytype = "ACCNUM"
)

probe2symbol <- merge(probe_to_acc, mapping, by.x = "GB_ACC_clean", by.y = "ACCNUM") %>%
  dplyr::select(ID, SYMBOL) %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  dplyr::distinct()

# 1.3 数据预处理：Log2 转换与分位数归一化
if (max(expr_raw, na.rm = TRUE) > 100) {
  expr_raw[expr_raw <= 0] <- 1
  expr_matrix_correct <- log2(expr_raw)
  message(">> 检测到原始数据未取 Log，已自动完成 Log2 转换。")
} else {
  expr_matrix_correct <- expr_raw
}

non_zero_var <- apply(expr_matrix_correct, 1, var, na.rm = TRUE) > 0
expr_matrix_correct <- expr_matrix_correct[non_zero_var, ]
expr_matrix_correct <- normalizeBetweenArrays(expr_matrix_correct)

# 1.4 生成以 Gene Symbol 为行名的全基因组表达矩阵
cat("⚙️ 正在转换探针并做均值合并，构建 full_gene_expr_matrix ...\n")

common_probes <- intersect(rownames(expr_matrix_correct), probe2symbol$ID)
expr_filtered <- expr_matrix_correct[common_probes, ]
matched_symbols <- probe2symbol$SYMBOL[match(common_probes, probe2symbol$ID)]

full_gene_expr_matrix <- aggregate(expr_filtered, list(Gene = matched_symbols), mean)
rownames(full_gene_expr_matrix) <- full_gene_expr_matrix$Gene
full_gene_expr_matrix$Gene <- NULL
full_gene_expr_matrix <- as.matrix(full_gene_expr_matrix)

cat(paste0("✨ 全基因组表达矩阵构建成功！共计 ", nrow(full_gene_expr_matrix), " 个基因，", ncol(full_gene_expr_matrix), " 个样本。\n\n"))


# ==============================================================================
# Part 2: 计算 AD 样本 ssGSEA 并绘制 Pearson 相关性热图
# ==============================================================================

library(GSVA)
library(clusterProfiler)
library(pheatmap)

# 2.1 读取本地 5 基因表达数据并锁定 AD 样本
cat("📂 1. 读取本地 5 基因表达数据...\n")
local_csv <- file.path(base_dir, "AD_5_Genes_Sample_Expression.csv")
df_local <- read.csv(local_csv, header = TRUE, row.names = 1, check.names = FALSE)

if (!("Group" %in% colnames(df_local))) {
  stop("⚠️ 错误：AD_5_Genes_Sample_Expression.csv 中未找到 'Group' 列！")
}

ad_df <- df_local[df_local$Group == "AD", ]
upstream_genes <- c("CLTC", "SLC2A13", "SLC7A14", "ACP5", "RAMP3")

df_upstream_expr <- t(as.matrix(ad_df[, upstream_genes]))

# 2.2 样本交集对齐与背景矩阵准备
ad_samples <- rownames(ad_df)
valid_ad_samples <- intersect(ad_samples, colnames(full_gene_expr_matrix))

if (length(valid_ad_samples) == 0) {
  stop("⚠️ 样本 ID 匹配失败！请检查 CSV 文件中的样本名与 GSE15222 表达矩阵列名是否对应。")
}

full_matrix_ad <- full_gene_expr_matrix[, valid_ad_samples, drop = FALSE]
df_upstream_expr <- df_upstream_expr[, valid_ad_samples, drop = FALSE]

# 2.3 读取 GMT 基因集并计算 ssGSEA 得分
cat("🧬 2. 正在基于全背景矩阵计算 AD 样本 ssGSEA 得分...\n")
gmt_all <- read.gmt(gmt_path)
gmt_all$term <- trimws(as.character(gmt_all$term))

pathway_mapping <- c(
  "KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION" = "Neuroactive ligand-receptor interaction",
  "KEGG_LONG_TERM_POTENTIATION"                   = "Long-term potentiation (Synaptic/Glutamate)",
  "KEGG_LONG_TERM_DEPRESSION"                     = "Long-term depression (Synaptic/GABA)",
  "KEGG_ADHERENS_JUNCTION"                        = "Adherens junction (Cadherin signaling)",
  "KEGG_TYPE_I_DIABETES_MELLITUS"                 = "Type I diabetes mellitus",
  "KEGG_CALCIUM_SIGNALING_PATHWAY"               = "Calcium signaling pathway"
)

gene_sets_official <- list()
for (pid in names(pathway_mapping)) {
  p_show_name <- pathway_mapping[pid]
  gmt_genes <- toupper(trimws(gmt_all$gene[gmt_all$term == pid]))
  valid_genes <- intersect(gmt_genes, rownames(full_matrix_ad))
  if (length(valid_genes) > 5) {
    gene_sets_official[[p_show_name]] <- valid_genes
  }
}

my_param <- ssgseaParam(exprData = full_matrix_ad, geneSets = gene_sets_official, normalize = TRUE)
df_downstream_ssgsea <- gsva(my_param)

# 2.4 计算 Pearson 相关性并导出热图
cat("📊 3. 计算 Pearson 相关性矩阵并绘制热图...\n")

cor_matrix <- cor(t(df_upstream_expr), t(df_downstream_ssgsea), method = "pearson")
cor_matrix_plot <- t(cor_matrix)

cor_matrix_plot[is.na(cor_matrix_plot)] <- 0
cor_matrix_plot[is.nan(cor_matrix_plot)] <- 0

# 锁死横纵轴顺序
cor_matrix_plot <- cor_matrix_plot[, upstream_genes, drop = FALSE]
cor_matrix_plot <- cor_matrix_plot[as.character(pathway_mapping), , drop = FALSE]

output_csv <- file.path(base_dir, "AD_Upstream_vs_Downstream_Pearson_r.csv")
write.csv(as.data.frame(cor_matrix_plot), file = output_csv)

output_heatmap_pdf <- file.path(base_dir, "AD_Upstream_vs_Downstream_Heatmap.pdf")
pdf(output_heatmap_pdf, width = 5.2, height = 3.6)

explicit_breaks <- seq(-1, 1, length.out = 101)
vibrant_colors <- colorRampPalette(c("#1A365D", "#4A5568", "#3182CE", "#FFFFFF", "#E53E3E", "#9B2C2C", "#741B1B"))(100)

pheatmap(
  cor_matrix_plot,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_color = "black",
  fontsize_number = 8.0,
  fontsize = 8.5,
  breaks = explicit_breaks,
  color = vibrant_colors,
  main = "Upstream Genes vs Downstream Pathways (AD Cohort)"
)

dev.off()
cat(paste0("✨ 脚本运行完毕！热图已成功导出至：'", output_heatmap_pdf, "'\n"))