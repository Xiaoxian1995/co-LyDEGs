# ==============================================================================
# 🎯 [Cross-Disease 5-Gene vs 23 Immune Cells Full Pipeline]
# 包含：
# 1. GSE15222 (AD) & GSE68719 (PD) 表达矩阵预处理与 23 种免疫细胞 ssGSEA 计算
# 2. 导出 ad_ssgsea_23_scores.csv 与 pd_ssgsea_23_scores.csv (存入 output/ 目录)
# 3. 读取 5 个基因表达，按 Group 列精准提取 AD/PD 病人样本 (排除 Control)
# 4. 计算 Spearman 相关性并生成最终 Cross-Disease 免疫热图 (存入 output/ 目录)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. 加载所需的依赖包与创建相对路径目录
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(GEOquery)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GSVA)
  library(GSEABase)
  library(pheatmap)
})

# 创建相对路径目录（数据输入目录与结果输出目录）
dir.create("data", showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)

target_genes <- c("CLTC", "SLC2A13", "SLC7A14", "ACP5", "RAMP3")

immune_markers <- list(
  `Activated B cell`            = c("CD19", "CD22", "MS4A1", "CD79A", "CD79B"),
  `Immature B cell`             = c("CD38", "FAIM3", "IL4R"),
  `Memory B cell`               = c("CD27", "TNFRSF13B", "CHST10"),
  `Activated CD4 T cell`        = c("CD4", "IL2RA", "CD69", "ITGA2"),
  `Activated CD8 T cell`        = c("CD8A", "CD8B", "GZMA", "GZMB", "IFNG"),
  `Eosinophil`                  = c("CCR3", "IL5RA", "CLC", "ALOX5"),
  `Macrophage`                  = c("CD68", "CD163", "MSR1", "CSF1R"),
  `Mast cell`                   = c("TPSAB1", "TPSB2", "CMA1", "FCER1A"),
  `MDSC`                        = c("CD11b", "CD33", "FUT4", "IL4R"),
  `Monocyte`                    = c("CD14", "FCGR1A", "CD64", "ITGAM"),
  `Natural killer cell`         = c("NCAM1", "NCR1", "KLRD1", "CD244"),
  `Neutrophil`                  = c("FCGR3B", "CXCR1", "CXCR2", "FPR1"),
  `Plasmacytoid dendritic cell` = c("IL3RA", "CLEC4C", "NRP1"),
  `T follicular helper cell`    = c("CXCR5", "PDCD1", "BCL6", "ICOS"),
  `Type 1 T helper cell`        = c("TBX21", "IFNG", "CXCR3"),
  `Type 17 T helper cell`       = c("RORC", "IL17A", "IL17F", "CCR6"),
  `Type 2 T helper cell`        = c("GATA3", "IL4", "IL5", "IL13"),
  `Activated dendritic cell`    = c("CD83", "CCR7", "LAMP3"),
  `Regulatory T cell`           = c("FOXP3", "IL2RA", "IKZF2"),
  `NK T cell`                   = c("CD1D", "KLRB1"),
  `Gamma delta T cell`          = c("TRGC1", "TRGC2", "TRDC"),
  `Immature dendritic cell`     = c("CD1A", "CD1C", "ITGAX"),
  `Macrophage M1`               = c("NOS2", "TLR4", "IL12A", "IL12B")
)

gset_obj <- GSEABase::GeneSetCollection(lapply(names(immune_markers), function(x) {
  GSEABase::GeneSet(immune_markers[[x]], setName = x)
}))

# ==============================================================================
# 2. GSE15222 (AD) 数据处理与 ssGSEA 计算
# ==============================================================================
cat("\n📂 [1/4] 开始处理 GSE15222 (AD) 数据...\n")
ad_file <- file.path("data", "GSE15222_series_matrix.txt.gz")
gpl_file <- file.path("data", "GPL2700.gz")

if (!file.exists(ad_file)) stop("❌ 未找到 data/GSE15222_series_matrix.txt.gz 文件！")
if (!file.exists(gpl_file)) stop("❌ 未找到 data/GPL2700.gz 文件！")

gse_ad <- getGEO(filename = ad_file, GSEMatrix = TRUE, getGPL = FALSE)
sample_info_ad <- pData(gse_ad) %>%
  dplyr::mutate(Group = case_when(
    str_detect(tolower(description), "control") ~ "Control",
    str_detect(tolower(description), "alzheimer") ~ "AD",
    TRUE ~ NA_character_
  ))

valid_ad <- !is.na(sample_info_ad$Group)
sample_info_ad <- sample_info_ad[valid_ad, ]
expr_raw_ad <- exprs(gse_ad)[, valid_ad]

# GPL2700 探针注释
gpl <- GEOquery::getGEO(filename = gpl_file)
gpl_meta <- Table(gpl)
probe_to_acc <- gpl_meta[, c("ID", "GB_ACC")]
probe_to_acc <- probe_to_acc[probe_to_acc$GB_ACC != "" & !is.na(probe_to_acc$GB_ACC), ]
probe_to_acc$GB_ACC_clean <- sub("\\..*$", "", probe_to_acc$GB_ACC)

mapping_ad <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(probe_to_acc$GB_ACC_clean),
                                    columns = "SYMBOL", keytype = "ACCNUM")

probe2symbol <- merge(probe_to_acc, mapping_ad, by.x = "GB_ACC_clean", by.y = "ACCNUM") %>%
  dplyr::select(ID, SYMBOL) %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  dplyr::distinct()

if (max(expr_raw_ad, na.rm = TRUE) > 100) {
  expr_raw_ad[expr_raw_ad <= 0] <- 1
  expr_raw_ad <- log2(expr_raw_ad)
}
non_zero_var_ad <- apply(expr_raw_ad, 1, var, na.rm = TRUE) > 0
expr_matrix_ad_norm <- normalizeBetweenArrays(expr_raw_ad[non_zero_var_ad, ])

translated_ad <- as.data.frame(expr_matrix_ad_norm) %>%
  rownames_to_column("ID") %>%
  inner_join(probe2symbol, by = "ID") %>%
  dplyr::select(-ID) %>%
  group_by(SYMBOL) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  ungroup()

ad_large_clean <- as.matrix(translated_ad[, -1])
rownames(ad_large_clean) <- translated_ad$SYMBOL

cat("🔄 计算 AD 样本 23 种免疫细胞 ssGSEA 得分...\n")
ssgsea_param_ad <- ssgseaParam(exprData = ad_large_clean, geneSets = gset_obj)
ad_ssgsea_23_matrix <- gsva(ssgsea_param_ad, verbose = FALSE)

ad_ssgsea_23_scores <- as.data.frame(t(ad_ssgsea_23_matrix)) %>%
  rownames_to_column(var = "Sample_ID") %>%
  left_join(sample_info_ad %>% rownames_to_column(var = "Sample_ID") %>% dplyr::select(Sample_ID, Group), by = "Sample_ID") %>%
  dplyr::select(Sample_ID, Group, everything())

ad_ssgsea_out <- file.path("output", "ad_ssgsea_23_scores.csv")
write.csv(ad_ssgsea_23_scores, file = ad_ssgsea_out, row.names = FALSE)
cat(sprintf("✅ AD ssGSEA 结果已成功保存为：'%s'\n", ad_ssgsea_out))

# ==============================================================================
# 3. GSE68719 (PD) 数据处理与 ssGSEA 计算
# ==============================================================================
cat("\n📂 [2/4] 开始处理 GSE68719 (PD) 数据...\n")
pd_series_file <- file.path("data", "GSE68719_series_matrix.txt.gz")
norm_counts_file <- file.path("data", "GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz")

if (!file.exists(pd_series_file)) stop("❌ 未找到 data/GSE68719_series_matrix.txt.gz 文件！")
if (!file.exists(norm_counts_file)) stop("❌ 未找到 data/GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz 文件！")

gse_pd <- getGEO(filename = pd_series_file, GSEMatrix = TRUE, getGPL = FALSE)
sample_info_pd <- pData(gse_pd) %>%
  mutate(
    Group = case_when(
      grepl("P", title, ignore.case = FALSE) ~ "PD",
      grepl("C", title, ignore.case = FALSE) ~ "Control",
      TRUE ~ NA_character_
    ),
    Clean_ID = str_trim(sub("\\[.*\\]", "", title))
  ) %>%
  filter(!is.na(Group))

rownames(sample_info_pd) <- sample_info_pd$Clean_ID

pd_expr_df <- read.table(gzfile(norm_counts_file), header = TRUE, sep = "\t", check.names = FALSE)
rownames(pd_expr_df) <- pd_expr_df[, 1]
expr_matrix_pd_real <- as.matrix(pd_expr_df[, -1])

common_pd_samples <- intersect(colnames(expr_matrix_pd_real), rownames(sample_info_pd))
expr_matrix_pd_final <- expr_matrix_pd_real[, common_pd_samples]
sample_info_pd_final <- sample_info_pd[common_pd_samples, ]

pd_ensembl_raw <- rownames(expr_matrix_pd_final)
pd_ensembl_clean <- sub("\\..*$", "", pd_ensembl_raw)

gene_ids_pd <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(pd_ensembl_clean),
                                     columns = "SYMBOL", keytype = "ENSEMBL")

pd_map_df <- data.frame(Raw_ID = pd_ensembl_raw, Ensembl_clean = pd_ensembl_clean, stringsAsFactors = FALSE) %>%
  inner_join(gene_ids_pd, by = c("Ensembl_clean" = "ENSEMBL")) %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  distinct(Raw_ID, SYMBOL, .keep_all = TRUE)

pd_sub_mat <- apply(expr_matrix_pd_final[pd_map_df$Raw_ID, ], 2, as.numeric)
if (max(pd_sub_mat, na.rm = TRUE) > 50) {
  pd_sub_mat <- log2(pd_sub_mat + 1)
}

unique_pd_symbols <- unique(pd_map_df$SYMBOL)
pd_large_clean <- matrix(0, nrow = length(unique_pd_symbols), ncol = length(common_pd_samples))
rownames(pd_large_clean) <- unique_pd_symbols
colnames(pd_large_clean) <- common_pd_samples

for (i in seq_along(unique_pd_symbols)) {
  sym <- unique_pd_symbols[i]
  idx <- which(pd_map_df$SYMBOL == sym)
  if (length(idx) == 1) {
    pd_large_clean[i, ] <- pd_sub_mat[idx, ]
  } else {
    pd_large_clean[i, ] <- colMeans(pd_sub_mat[idx, , drop = FALSE], na.rm = TRUE)
  }
}

cat("🔄 计算 PD 样本 23 种免疫细胞 ssGSEA 得分...\n")
ssgsea_param_pd <- ssgseaParam(exprData = pd_large_clean, geneSets = gset_obj)
pd_ssgsea_23_matrix <- gsva(ssgsea_param_pd, verbose = FALSE)

pd_ssgsea_23_scores <- as.data.frame(t(pd_ssgsea_23_matrix)) %>%
  rownames_to_column(var = "Sample_ID") %>%
  left_join(sample_info_pd_final %>% rownames_to_column(var = "Sample_ID") %>% dplyr::select(Sample_ID, Group), by = "Sample_ID") %>%
  dplyr::select(Sample_ID, Group, everything())

pd_ssgsea_out <- file.path("output", "pd_ssgsea_23_scores.csv")
write.csv(pd_ssgsea_23_scores, file = pd_ssgsea_out, row.names = FALSE)
cat(sprintf("✅ PD ssGSEA 结果已成功保存为：'%s'\n", pd_ssgsea_out))

# ==============================================================================
# 4. 精准根据 Group 列提取纯 AD / PD 病人样本并进行相关性计算
# ==============================================================================
cat("\n📊 [3/4] 根据 Group 列精准筛选 AD 与 PD 病人样本...\n")

ad_ssgsea_df <- read.csv(file.path("output", "ad_ssgsea_23_scores.csv"), row.names = 1, check.names = FALSE)
pd_ssgsea_df <- read.csv(file.path("output", "pd_ssgsea_23_scores.csv"), row.names = 1, check.names = FALSE)

ad_expr_file <- file.path("data", "AD_5_Genes_Sample_Expression.csv")
pd_expr_file <- file.path("data", "PD_5_Genes_Sample_Expression.csv")

if (!file.exists(ad_expr_file)) stop("❌ 未找到 data/AD_5_Genes_Sample_Expression.csv 文件！")
if (!file.exists(pd_expr_file)) stop("❌ 未找到 data/PD_5_Genes_Sample_Expression.csv 文件！")

ad_expr_df <- read.csv(ad_expr_file, row.names = 1, check.names = FALSE)
pd_expr_df <- read.csv(pd_expr_file, row.names = 1, check.names = FALSE)

ad_patient_ids <- rownames(ad_ssgsea_df)[toupper(trimws(ad_ssgsea_df$Group)) == "AD"]
pd_patient_ids <- rownames(pd_ssgsea_df)[toupper(trimws(pd_ssgsea_df$Group)) == "PD"]

if (!any(target_genes %in% rownames(ad_expr_df)) && any(target_genes %in% colnames(ad_expr_df))) {
  ad_expr_df <- as.data.frame(t(ad_expr_df))
}
if (!any(target_genes %in% rownames(pd_expr_df)) && any(target_genes %in% colnames(pd_expr_df))) {
  pd_expr_df <- as.data.frame(t(pd_expr_df))
}

ad_ssgsea_mat <- as.data.frame(t(ad_ssgsea_df[, -1]))
pd_ssgsea_mat <- as.data.frame(t(pd_ssgsea_df[, -1]))

ad_samples <- intersect(ad_patient_ids, intersect(colnames(ad_expr_df), colnames(ad_ssgsea_mat)))
pd_samples <- intersect(pd_patient_ids, intersect(colnames(pd_expr_df), colnames(pd_ssgsea_mat)))

cat(paste0("   - 🟣 AD 病人组精准样本数: ", length(ad_samples), " 个 (已排除 Control)\n"))
cat(paste0("   - 🔵 PD 病人组精准样本数: ", length(pd_samples), " 个 (已排除 Control)\n"))

if(length(ad_samples) == 0 || length(pd_samples) == 0) {
  stop("⚠️ 报错：未能在 Group 列中成功匹配到 AD 或 PD 病人样本，请检查输入表格！")
}

actual_cells <- rownames(ad_ssgsea_mat)

correct_col_order <- c()
for(g in target_genes) {
  correct_col_order <- c(correct_col_order, paste0(g, " (AD)"), paste0(g, " (PD)"))
}

plot_mat <- matrix(NA, nrow = length(actual_cells), ncol = length(correct_col_order), dimnames = list(actual_cells, correct_col_order))
pval_mat <- matrix(NA, nrow = length(actual_cells), ncol = length(correct_col_order), dimnames = list(actual_cells, correct_col_order))

safe_cor_test <- function(x, y) {
  valid_idx <- is.finite(x) & is.finite(y)
  x_clean <- x[valid_idx]
  y_clean <- y[valid_idx]
  if (length(x_clean) < 3 || sd(x_clean) == 0 || sd(y_clean) == 0) return(list(estimate = 0, p.value = 1))
  tryCatch({ cor.test(x_clean, y_clean, method = "spearman", exact = FALSE) }, error = function(e) list(estimate = 0, p.value = 1))
}

cat("🧬 正在计算 5 个目标基因与 23 种免疫细胞的 Spearman 相关性...\n")
for(g in target_genes) {
  for(cell in actual_cells) {
    if (g %in% rownames(ad_expr_df) && cell %in% rownames(ad_ssgsea_mat)) {
      res_ad <- safe_cor_test(as.numeric(ad_expr_df[g, ad_samples]), as.numeric(ad_ssgsea_mat[cell, ad_samples]))
      plot_mat[cell, paste0(g, " (AD)")] <- res_ad$estimate
      pval_mat[cell, paste0(g, " (AD)")] <- res_ad$p.value
    }
    if (g %in% rownames(pd_expr_df) && cell %in% rownames(pd_ssgsea_mat)) {
      res_pd <- safe_cor_test(as.numeric(pd_expr_df[g, pd_samples]), as.numeric(pd_ssgsea_mat[cell, pd_samples]))
      plot_mat[cell, paste0(g, " (PD)")] <- res_pd$estimate
      pval_mat[cell, paste0(g, " (PD)")] <- res_pd$p.value
    }
  }
}

plot_mat[is.na(plot_mat)] <- 0
pval_mat[is.na(pval_mat)] <- 1

display_signif <- matrix("", nrow = nrow(pval_mat), ncol = ncol(pval_mat))
display_signif[pval_mat < 0.01] <- "**"
display_signif[pval_mat >= 0.01 & pval_mat < 0.05] <- "*"

# ==============================================================================
# 5. 渲染与导出热图
# ==============================================================================
cat("\n🎨 [4/4] 正在渲染 Cross-Disease 病人组免疫相关性热图...\n")

heatmap_colors <- colorRampPalette(c("#1E466E", "#6BAED6", "#F7FBFF", "#E6550D", "#A63603"))(100)

heatmap_pdf_out <- file.path("output", "Group_Patients_Immune_Correlation_Heatmap.pdf")
pdf(heatmap_pdf_out, width = 10, height = 9)

pheatmap(
  mat               = plot_mat,
  color             = heatmap_colors,
  breaks            = seq(-0.7, 0.7, length.out = 101),
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  display_numbers   = display_signif,
  number_color      = "black",
  fontsize_number   = 12,
  cellwidth         = 45,
  cellheight        = 22,
  border_color      = "#CCCCCC",
  fontsize_row      = 11,
  fontsize_col      = 12,
  angle_col         = 45,
  main              = "Spearman Correlation in Patients Only (Group Classified): 5 Genes vs 23 Immune Cells"
)

dev.off()

cat("\n🎉 【全流程运行完成！】生成结果如下：\n")
cat(sprintf("  👉 1. %s\n", file.path("output", "ad_ssgsea_23_scores.csv")))
cat(sprintf("  👉 2. %s\n", file.path("output", "pd_ssgsea_23_scores.csv")))
cat(sprintf("  👉 3. %s\n", heatmap_pdf_out))