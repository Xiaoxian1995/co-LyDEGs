# ==============================================================================
# 🎯 [AD & PD 全流程整合：相对路径通用版 (差异分析 -> 共同基因 -> KEGG 富集与矢量出图)]
# ==============================================================================

# ---- 【第 0 步】网络协议、环境依赖与当前相对路径配置 ----
options(clusterProfiler.download.method = "curl")
options(download.file.extra = "-k") # 规避 rest.kegg.jp 的 SSL 证书校验报错
options(timeout = 300)

# 打印并确认当前工作路径（相对路径基准点）
current_dir <- getwd()
message(sprintf("📂 当前工作目录 (Relative Path Base): %s", current_dir))

message("🧬 正在加载分析所需核心 R 包...")
suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(GEOquery)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(ggrepel)
  library(clusterProfiler)
  library(enrichplot)
})

# ==============================================================================
# STEP 1: 分析 GSE15222 (AD vs Control) - 相对路径读取
# ==============================================================================
message("\n🚀 [1/4] 正在运行 AD (GSE15222) 差异表达分析...")

file_gse15222 <- "./GSE15222_series_matrix.txt.gz"
file_gpl2700   <- "./GPL2700.gz"

if (!file.exists(file_gse15222)) {
  stop(sprintf("⚠️ 找不到文件 '%s'，请确认该文件在当前工作目录下！", file_gse15222))
}

gse_ad <- getGEO(filename = file_gse15222, GSEMatrix = TRUE, getGPL = FALSE)
sample_info_ad <- pData(gse_ad) %>%
  dplyr::mutate(Group = case_when(
    str_detect(tolower(description), "control") ~ "Control",
    str_detect(tolower(description), "alzheimer") ~ "AD",
    TRUE ~ NA_character_
  ))

valid_ad <- !is.na(sample_info_ad$Group)
sample_info_ad <- sample_info_ad[valid_ad, ]
expr_ad_raw <- exprs(gse_ad)[, valid_ad]
group_list_ad <- factor(sample_info_ad$Group, levels = c("Control", "AD"))

# 相对路径读取注释 GPL 平台
if (!file.exists(file_gpl2700)) {
  stop(sprintf("⚠️ 找不到文件 '%s'，请确认该文件在当前工作目录下！", file_gpl2700))
}

gpl2700 <- GEOquery::getGEO(filename = file_gpl2700)
gpl_meta <- Table(gpl2700)
probe_to_acc <- gpl_meta[gpl_meta$GB_ACC != "" & !is.na(gpl_meta$GB_ACC), c("ID", "GB_ACC")]
probe_to_acc$GB_ACC_clean <- sub("\\..*$", "", probe_to_acc$GB_ACC)

mapping_ad <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(probe_to_acc$GB_ACC_clean),
  columns = "SYMBOL",
  keytype = "ACCNUM"
)

probe2symbol_ad <- merge(probe_to_acc, mapping_ad, by.x = "GB_ACC_clean", by.y = "ACCNUM") %>%
  dplyr::select(ID, SYMBOL) %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  dplyr::distinct()

# 标准化与 limma 分析
if (max(expr_ad_raw, na.rm = TRUE) > 100) {
  expr_ad_raw[expr_ad_raw <= 0] <- 1
  expr_ad <- log2(expr_ad_raw)
} else {
  expr_ad <- expr_ad_raw
}
expr_ad <- expr_ad[apply(expr_ad, 1, var, na.rm = TRUE) > 0, ]
expr_ad <- normalizeBetweenArrays(expr_ad)

design_ad <- model.matrix(~group_list_ad)
fit_ad <- eBayes(lmFit(expr_ad, design_ad))
deg_ad <- topTable(fit_ad, coef = 2, number = Inf) %>%
  rownames_to_column("ID") %>%
  inner_join(probe2symbol_ad, by = "ID") %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  dplyr::arrange(P.Value) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

ad_plot_data <- deg_ad %>%
  dplyr::mutate(
    Group = case_when(
      logFC > 0.5 & P.Value < 0.05 ~ "Up",
      logFC < -0.5 & P.Value < 0.05 ~ "Down",
      TRUE ~ "None"
    )
  )

# ==============================================================================
# STEP 2: 分析 GSE68719 (PD vs Control) - 相对路径读取
# ==============================================================================
message("\n🚀 [2/4] 正在运行 PD (GSE68719) 差异表达分析...")

file_gse68719      <- "./GSE68719_series_matrix.txt.gz"
file_gse68719_counts <- "./GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz"

if (!file.exists(file_gse68719) || !file.exists(file_gse68719_counts)) {
  stop("⚠️ 找不到 GSE68719 相关表达矩阵文件，请检查当前文件夹！")
}

gse_pd <- getGEO(filename = file_gse68719, GSEMatrix = TRUE, getGPL = FALSE)
sample_info_pd <- pData(gse_pd) %>%
  dplyr::mutate(
    Group = case_when(
      grepl("P", title, ignore.case = FALSE) ~ "PD",
      grepl("C", title, ignore.case = FALSE) ~ "Control",
      TRUE ~ NA_character_
    ),
    Clean_ID = str_trim(sub("\\[.*\\]", "", title))
  ) %>%
  dplyr::filter(!is.na(Group))

rownames(sample_info_pd) <- sample_info_pd$Clean_ID

pd_expr_df <- read.table(gzfile(file_gse68719_counts), header = TRUE, sep = "\t", check.names = FALSE)
rownames(pd_expr_df) <- pd_expr_df[, 1]
expr_pd_raw <- as.matrix(pd_expr_df[, -1])

common_pd_samples <- intersect(colnames(expr_pd_raw), rownames(sample_info_pd))
expr_pd_final <- expr_pd_raw[, common_pd_samples]
sample_info_pd_final <- sample_info_pd[common_pd_samples, ]

expr_pd_num <- apply(expr_pd_final, 2, as.numeric)
rownames(expr_pd_num) <- rownames(expr_pd_final)
if (max(expr_pd_num, na.rm = TRUE) > 50) {
  expr_pd_log <- log2(expr_pd_num + 1)
} else {
  expr_pd_log <- expr_pd_num
}

group_pd <- factor(sample_info_pd_final$Group, levels = c("Control", "PD"))
design_pd <- model.matrix(~group_pd)
fit_pd <- eBayes(lmFit(expr_pd_log, design_pd))

deg_pd_all <- topTable(fit_pd, coef = 2, number = Inf, adjust.method = "BH") %>%
  rownames_to_column("Ensembl_ID")
deg_pd_all$Ensembl_clean <- sub("\\..*$", "", deg_pd_all$Ensembl_ID)

gene_ids_pd <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = deg_pd_all$Ensembl_clean,
  columns = "SYMBOL",
  keytype = "ENSEMBL"
)

deg_pd_mapped <- dplyr::left_join(deg_pd_all, gene_ids_pd, by = c("Ensembl_clean" = "ENSEMBL")) %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

p_col_pd <- if("P.Value" %in% colnames(deg_pd_mapped)) "P.Value" else "adj.P.Val"

pd_plot_data <- deg_pd_mapped %>%
  dplyr::mutate(
    Group = case_when(
      logFC > 0.5 & .data[[p_col_pd]] < 0.05 ~ "Up",
      logFC < -0.5 & .data[[p_col_pd]] < 0.05 ~ "Down",
      TRUE ~ "None"
    )
  )

# ==============================================================================
# STEP 3: 提取共同差异基因并导出 Shared_All_DEGs.csv (相对路径保存)
# ==============================================================================
message("\n🚀 [3/4] 正在提取 AD & PD 共同差异基因...")

ad_deg_genes <- ad_plot_data %>% 
  dplyr::filter(Group %in% c("Up", "Down")) %>% 
  dplyr::pull(SYMBOL) %>% 
  unique()

pd_deg_genes <- pd_plot_data %>% 
  dplyr::filter(Group %in% c("Up", "Down")) %>% 
  dplyr::pull(SYMBOL) %>% 
  unique()

common_deg_genes <- intersect(ad_deg_genes, pd_deg_genes)
common_gene_count <- length(common_deg_genes)

cat("\n====================================================================\n")
cat(sprintf("🔴 AD 显著差异基因数: %d 个\n", length(ad_deg_genes)))
cat(sprintf("🔵 PD 显著差异基因数: %d 个\n", length(pd_deg_genes)))
cat(sprintf("🔥 【AD & PD 共同差异基因总数】：%d 个\n", common_gene_count))
cat("====================================================================\n\n")

if (common_gene_count == 0) {
  stop("⚠️ 未找到重叠的显著差异基因，请检查差异基因筛选阈值！")
}

# 💾 【保存文件 1】：导出共同基因 CSV (相对路径)
file_out_degs <- "./Shared_All_DEGs.csv"
write.csv(data.frame(Symbol = common_deg_genes), file = file_out_degs, row.names = FALSE, quote = FALSE)
message(sprintf("💾 共同基因已导出至相对路径：`%s`", file_out_degs))

# ==============================================================================
# STEP 4: KEGG 富集计算与 7.33 x 3.7 横向矢量图导出 (相对路径保存)
# ==============================================================================
message("\n🚀 [4/4] 正在运行 KEGG 富集分析并制图...")

# 转换基因名至 Entrez ID
gene_ids <- bitr(common_deg_genes, 
                 fromType = "SYMBOL", 
                 toType   = "ENTREZID", 
                 OrgDb    = org.Hs.eg.db)

message("🔬 后台检索 KEGG 官方数据库...")
kegg_res <- enrichKEGG(gene         = gene_ids$ENTREZID, 
                       organism     = 'hsa', 
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.2)

if (is.null(kegg_res) || nrow(as.data.frame(kegg_res)) == 0) {
  message("⚠️ 提示：标准阈值下未富集到通路，正在放宽标准重新检索...")
  kegg_res <- enrichKEGG(gene         = gene_ids$ENTREZID, 
                         organism     = 'hsa', 
                         pvalueCutoff = 0.2, 
                         qvalueCutoff = 0.5)
}

if (is.null(kegg_res) || nrow(as.data.frame(kegg_res)) == 0) {
  stop("❌ 未能在 KEGG 中匹配到通路，请检查基因格式。")
}

# 将 Entrez ID 转换回 SYMBOL 方便写论文阅读
kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
df_kegg_results <- as.data.frame(kegg_res)

# 💾 【保存文件 2】：导出 KEGG 富集表格 (相对路径)
file_out_kegg_csv <- "./Shared_All_DEGs_KEGG_Results.csv"
write.csv(df_kegg_results, file = file_out_kegg_csv, row.names = FALSE)
message(sprintf("💾 KEGG 富集结果表已保存至相对路径：`%s`", file_out_kegg_csv))

# ---- 准备画图数据 (取前 20 个通路) ----
p_col <- if ("p.adjust" %in% colnames(df_kegg_results)) "p.adjust" else "pvalue"

kegg_plot_data <- df_kegg_results %>%
  head(20) %>%
  mutate(
    p.adjust = .data[[p_col]],
    log10P = -log10(p.adjust),
    Description = factor(Description, levels = rev(Description))
  )

# ---- 构建 ggplot 对象 ----
p_custom_bubble_stretched <- ggplot(kegg_plot_data, aes(x = log10P, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust), shape = 16, stroke = 0, alpha = 0.9) +
  scale_color_gradientn(colors = c("#E64B35FF", "#F39B7FFF", "#4DBBD5FF", "#3C5488FF")) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  scale_size_continuous(range = c(2.8, 6.5)) +
  labs(
    title = sprintf("KEGG Enrichment of Common DEGs (%d Genes)\nRe-analyzed Pathway Top 20", common_gene_count),
    x = "-log10(Adjusted P-value)",
    y = "KEGG Pathway Description",
    size = "Gene Count",
    color = "p.adjust"
  ) +
  theme_bw() +
  theme(
    text               = element_text(family = "Helvetica"),
    plot.title         = element_text(face = "bold", hjust = 0.5, size = 9.5),
    axis.title         = element_text(face = "bold", size = 8.5),
    axis.text.y        = element_text(size = 8, color = "black"),
    axis.text.x        = element_text(size = 8, color = "black"),
    legend.title       = element_text(face = "bold", size = 8),
    legend.text        = element_text(size = 8),
    legend.box.spacing = unit(0.1, "cm"),
    panel.grid.major.x = element_line(color = "gray95"),
    panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
    panel.grid.minor   = element_blank()
  )

# 💾 【保存文件 3】：导出 7.33 x 3.7 横向加宽矢量 PDF (相对路径)
file_out_pdf <- "./ReAnalyzed_KEGG_Bubble_Plot_Stretched.pdf"
output_width  <- 5.5 * 4 / 3  # 7.333 Inches

pdf(
  file        = file_out_pdf,
  width       = output_width,
  height      = 3.7,
  family      = "Helvetica", # AI 活字兼容字体
  useDingbats = FALSE        # 确保圆点和文本在 AI 里均可直接鼠标点击编辑
)
print(p_custom_bubble_stretched)
dev.off()

cat("\n========================================================================\n")
cat("🎉 全流程运行完毕！相对路径下成功生成以下 3 个关键文件：\n")
cat(sprintf("1. 基因列表文件：`%s`\n", file_out_degs))
cat(sprintf("2. 通路结果表格：`%s`\n", file_out_kegg_csv))
cat(sprintf("3. 矢量 PDF 图像：`%s` (7.33 x 3.7 Inches)\n", file_out_pdf))
cat("========================================================================\n")