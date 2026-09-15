# ==============================================================================
# GSE68719 帕金森病 (PD) 全基因组表达矩阵构建与 ssGSEA Pearson 相关性热图绘制
# ==============================================================================

# ------------------------------------------------------------------------------
# Part 0: 环境准备与 GMT 基因集自动检查/下载
# ------------------------------------------------------------------------------
library(GEOquery)
library(limma)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)
library(GSVA)
library(clusterProfiler)
library(pheatmap)

# 设置基础相对路径
base_dir <- "./"

# 检查并自动下载 GMT 基因集文件
gmt_file_name <- "c2.cp.kegg_legacy.v2026.1.Hs.symbols.gmt"
gmt_path <- file.path(base_dir, gmt_file_name)
gmt_url <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/c2.cp.kegg_legacy.v2026.1.Hs.symbols.gmt"

if (!file.exists(gmt_path)) {
  message(">> 本地未检测到 GMT 基因集文件，正在尝试从 MSigDB 自动下载...")
  tryCatch({
    download.file(url = gmt_url, destfile = gmt_path, mode = "wb")
    message(">> GMT 文件下载成功！保存路径为：", gmt_path)
  }, error = function(e) {
    cat("\n⚠️ 自动下载失败！请检查网络连接或手动下载该文件：\n")
    cat("🔗 MSigDB 官方下载页面: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp\n")
    cat("🔗 直接下载链接: ", gmt_url, "\n")
    cat("将下载好的文件命名为 '", gmt_file_name, "' 并存放在当前工作目录下：", normalizePath(base_dir), "\n\n")
  })
} else {
  message(">> 已检测到本地 GMT 文件：", gmt_path)
}

# ------------------------------------------------------------------------------
# Part 1: 构建 PD (GSE68719) 全基因组表达矩阵 (full_gene_expr_matrix)
# ------------------------------------------------------------------------------

# 1.1 定义相对路径文件
pd_series_file  <- file.path(base_dir, "GSE68719_series_matrix.txt.gz")
norm_counts_file <- file.path(base_dir, "GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz")

# 1.2 读取临床表型信息与数据清洗
message("📂 1. 读取本地 PD 样本及临床表型信息...")
gse_pd <- getGEO(filename = pd_series_file, GSEMatrix = TRUE, getGPL = FALSE)
sample_info_pd <- pData(gse_pd)

# 提取表型分组 (P 为 PD，C 为 Control)
sample_info_pd_clean <- sample_info_pd %>%
  mutate(
    Group = case_when(
      grepl("P", title, ignore.case = FALSE) ~ "PD",
      grepl("C", title, ignore.case = FALSE) ~ "Control",
      TRUE ~ NA_character_
    ),
    Clean_ID = str_trim(sub("\\[.*\\]", "", title))
  )

valid_pd_samples <- sample_info_pd_clean %>% 
  filter(!is.na(Group))

rownames(valid_pd_samples) <- valid_pd_samples$Clean_ID

# 1.3 读取表达矩阵并对齐样本
message("📂 2. 读取 DESeq2 标准化表达矩阵并对齐样本...")
pd_expr_df <- read.table(gzfile(norm_counts_file), header = TRUE, sep = "\t", check.names = FALSE)

rownames(pd_expr_df) <- pd_expr_df[, 1]
expr_matrix_pd_real <- as.matrix(pd_expr_df[, -1])

# 样本 ID 精准对齐
common_pd_samples <- intersect(colnames(expr_matrix_pd_real), rownames(valid_pd_samples))

if (length(common_pd_samples) == 0) {
  stop("❌ 样本 ID 无法匹配，请检查表达矩阵列名与 Clean_ID 是否一致！")
}

expr_matrix_pd_final <- expr_matrix_pd_real[, common_pd_samples]
sample_info_pd_final <- valid_pd_samples[common_pd_samples, ]

# 1.4 数值转换与 Log2 预处理 + 分位数归一化
expr_numeric <- apply(expr_matrix_pd_final, 2, as.numeric)
rownames(expr_numeric) <- rownames(expr_matrix_pd_final)

if (max(expr_numeric, na.rm = TRUE) > 50) {
  message(">> 检测到数据未取 Log，正在自动进行 log2(x + 1) 转换...")
  expr_matrix_pd_log <- log2(expr_numeric + 1)
} else {
  expr_matrix_pd_log <- expr_numeric
}

# 过滤低方差基因并做归一化
non_zero_var <- apply(expr_matrix_pd_log, 1, var, na.rm = TRUE) > 0
expr_matrix_pd_log <- expr_matrix_pd_log[non_zero_var, ]
expr_matrix_pd_norm <- normalizeBetweenArrays(expr_matrix_pd_log)

# 1.5 Ensembl ID 转换为 Gene Symbol 并做均值合并
message("⚙️ 3. 正在将 Ensembl ID 映射为 Gene Symbol 并构建 full_gene_expr_matrix ...")

ensembl_ids_clean <- sub("\\..*$", "", rownames(expr_matrix_pd_norm))

mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(ensembl_ids_clean),
  columns = "SYMBOL",
  keytype = "ENSEMBL"
)

# 构建探针到 Symbol 的映射表
mapping_df <- data.frame(
  Ensembl_Raw = rownames(expr_matrix_pd_norm),
  Ensembl_Clean = ensembl_ids_clean,
  stringsAsFactors = FALSE
) %>%
  left_join(mapping, by = c("Ensembl_Clean" = "ENSEMBL")) %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  distinct(Ensembl_Raw, SYMBOL, .keep_all = TRUE)

# 过滤共有的基因并聚合均值
common_rows <- intersect(rownames(expr_matrix_pd_norm), mapping_df$Ensembl_Raw)
expr_filtered <- expr_matrix_pd_norm[common_rows, ]
matched_symbols <- mapping_df$SYMBOL[match(common_rows, mapping_df$Ensembl_Raw)]

full_gene_expr_matrix <- aggregate(expr_filtered, list(Gene = matched_symbols), mean)
rownames(full_gene_expr_matrix) <- full_gene_expr_matrix$Gene
full_gene_expr_matrix$Gene <- NULL
full_gene_expr_matrix <- as.matrix(full_gene_expr_matrix)

# 1.6 输出节点数据核对结果
cat("\n====================================================================\n")
cat("📌 【GSE68719 构建全基因组表达矩阵各节点数据核对】\n")
cat("--------------------------------------------------------------------\n")
cat("1. 原始表达矩阵 (`norm_counts`) 行数 (基因数):  ", nrow(pd_expr_df), " 个\n")
cat("2. 原始表达矩阵 (`norm_counts`) 列数 (样本数):  ", ncol(pd_expr_df), " 个\n")
cat("3. 成功对齐临床表型的总样本数:                    ", ncol(expr_matrix_pd_final), " 个\n")
cat("   - 其中 帕金森病 (PD) 患者例数:                  ", sum(sample_info_pd_final$Group == "PD"), " 例\n")
cat("   - 其中 正常对照 (Control) 例数:                ", sum(sample_info_pd_final$Group == "Control"), " 例\n")
cat("4. 映射注释转换后最终的 Gene Symbol 数量:          ", nrow(full_gene_expr_matrix), " 个\n")
cat("5. 构建好的 `full_gene_expr_matrix` 最终维度:    ", nrow(full_gene_expr_matrix), " 行 × ", ncol(full_gene_expr_matrix), " 列\n")
cat("====================================================================\n\n")


# ==============================================================================
# Part 2: 计算 PD 样本 ssGSEA 并绘制 Pearson 相关性热图
# ==============================================================================

# 2.1 读取本地 5 基因表达数据并锁定 PD 样本
local_pd_csv <- file.path(base_dir, "PD_5_Genes_Sample_Expression.csv")
cat("📂 4. 读取本地 5 基因表达数据 (", local_pd_csv, ")...\n", sep = "")
df_local_pd <- read.csv(local_pd_csv, header = TRUE, row.names = 1, check.names = FALSE)

if (!("Group" %in% colnames(df_local_pd))) {
  stop("⚠️ 错误：PD_5_Genes_Sample_Expression.csv 中未找到 'Group' 列！")
}

pd_df <- df_local_pd[df_local_pd$Group == "PD", ]
upstream_genes <- c("CLTC", "SLC2A13", "SLC7A14", "ACP5", "RAMP3")

df_upstream_expr <- t(as.matrix(pd_df[, upstream_genes]))

# 2.2 样本交集对齐与背景矩阵准备
pd_samples <- rownames(pd_df)
valid_pd_samples_ssgsea <- intersect(pd_samples, colnames(full_gene_expr_matrix))

if (length(valid_pd_samples_ssgsea) == 0) {
  stop("⚠️ 样本 ID 匹配失败！请检查 PD CSV 文件中的样本名与 GSE68719 表达矩阵列名是否对应。")
}

full_matrix_pd <- full_gene_expr_matrix[, valid_pd_samples_ssgsea, drop = FALSE]
df_upstream_expr <- df_upstream_expr[, valid_pd_samples_ssgsea, drop = FALSE]

# 2.3 读取 GMT 基因集并计算 ssGSEA 得分
cat("🧬 5. 正在基于全背景矩阵计算 PD 样本 ssGSEA 得分...\n")
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
  valid_genes <- intersect(gmt_genes, rownames(full_matrix_pd))
  if (length(valid_genes) > 5) {
    gene_sets_official[[p_show_name]] <- valid_genes
  }
}

my_param <- ssgseaParam(exprData = full_matrix_pd, geneSets = gene_sets_official, normalize = TRUE)
df_downstream_ssgsea <- gsva(my_param)

# 2.4 计算 Pearson 相关性并导出热图
cat("📊 6. 计算 Pearson 相关性矩阵并绘制 PD 热图...\n")

cor_matrix <- cor(t(df_upstream_expr), t(df_downstream_ssgsea), method = "pearson")
cor_matrix_plot <- t(cor_matrix)

cor_matrix_plot[is.na(cor_matrix_plot)] <- 0
cor_matrix_plot[is.nan(cor_matrix_plot)] <- 0

# 锁死横纵轴顺序
cor_matrix_plot <- cor_matrix_plot[, upstream_genes, drop = FALSE]
cor_matrix_plot <- cor_matrix_plot[as.character(pathway_mapping), , drop = FALSE]

# 保存 Pearson r 值矩阵到 CSV
output_csv <- file.path(base_dir, "PD_Upstream_vs_Downstream_Pearson_r.csv")
write.csv(as.data.frame(cor_matrix_plot), file = output_csv)

# 导出学术级 PDF 热图
output_heatmap_pdf <- file.path(base_dir, "PD_Upstream_vs_Downstream_Heatmap.pdf")
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
  main = "Upstream Genes vs Downstream Pathways (PD Cohort)"
)

dev.off()

cat(paste0("\n✨ GSE68719 (PD) 脚本运行完毕！\n"))
cat(paste0("📊 Pearson r 相关性矩阵已导出至：'", output_csv, "'\n"))
cat(paste0("🖼️ 最终热图已成功导出至：'", output_heatmap_pdf, "'\n"))