# 🎯 [Cross-Disease 跨疾病合并热图·自动过滤零方差基因版（相对路径）]
# ==============================================================================
library(pheatmap)
library(dplyr)

cat("📊 1. 正在通过相对路径读取 AD 与 PD 数据并自动构建矩阵...\n")

# 🛠️ 设置相对路径（支持直接位于当前工作目录，或放在子目录 data/ 中）
possible_ad_paths <- c(
  "AD_Lysosome_Subtypes_With_Expression_and_PCA.csv",
  "data/AD_Lysosome_Subtypes_With_Expression_and_PCA.csv"
)
possible_pd_paths <- c(
  "PD_Lysosome_Subtypes_With_Expression_and_PCA.csv",
  "data/PD_Lysosome_Subtypes_With_Expression_and_PCA.csv"
)

ad_csv_path <- Find(file.exists, possible_ad_paths)
pd_csv_path <- Find(file.exists, possible_pd_paths)

if (is.null(ad_csv_path) || is.null(pd_csv_path)) {
  stop("⚠️ 错误：未能在当前工作目录或 ./data/ 文件夹下找到 AD/PD 的 CSV 数据文件！请检查文件存放位置。")
}

cat(sprintf("📁 已成功定位 AD 文件：%s\n", ad_csv_path))
cat(sprintf("📁 已成功定位 PD 文件：%s\n", pd_csv_path))

# 1. 载入原始 CSV 数据
df_ad_raw <- read.csv(ad_csv_path, stringsAsFactors = FALSE, header = TRUE)
df_pd_raw <- read.csv(pd_csv_path, stringsAsFactors = FALSE, header = TRUE)

# ---- [数据清洗辅助函数] ----
process_raw_data <- function(df) {
  if ("Subtype" %in% colnames(df)) {
    colnames(df)[colnames(df) == "Subtype"] <- "Lysosome_Subtype"
  }
  df$Lysosome_Subtype <- gsub("Cluster", "Cluster_", df$Lysosome_Subtype)
  df$Lysosome_Subtype <- gsub("Cluster__", "Cluster_", df$Lysosome_Subtype)
  df$Lysosome_Subtype <- ifelse(grepl("1", df$Lysosome_Subtype), "Cluster_1", "Cluster_2")
  
  if (!("Sample_ID" %in% colnames(df))) {
    if ("X" %in% colnames(df)) {
      df$Sample_ID <- as.character(df$X)
    } else {
      df$Sample_ID <- paste0("Sample_", 1:nrow(df))
    }
  }
  return(df)
}

df_ad <- process_raw_data(df_ad_raw)
df_pd <- process_raw_data(df_pd_raw)

# 2. 提取基因列并取交集
meta_cols <- c(
  "Sample_ID", "X", "Final_Sample_ID", "Lysosome_Subtype", "Subtype", "Cluster",
  "Disease", "Center", "PC1", "PC2", "PC3", "PC4", "PC5", "Age", "Gender", "Sex", "RIN", "PMI", "Batch"
)

gene_cols_ad <- setdiff(colnames(df_ad), meta_cols)
gene_cols_pd <- setdiff(colnames(df_pd), meta_cols)
common_genes <- intersect(gene_cols_ad, gene_cols_pd)

if (length(common_genes) == 0) {
  stop("⚠️ 错误：未在 AD 与 PD 表格中匹配到共有的基因列！")
}

# 3. 构建数值型表达量矩阵 (行: 基因, 列: 样本)
df_ad_sorted <- df_ad[order(df_ad$Lysosome_Subtype), ]
mat_ad_clean <- t(apply(df_ad_sorted[, common_genes, drop = FALSE], 2, function(x) as.numeric(as.character(x))))
colnames(mat_ad_clean) <- df_ad_sorted$Sample_ID

df_pd_sorted <- df_pd[order(df_pd$Lysosome_Subtype), ]
mat_pd_clean <- t(apply(df_pd_sorted[, common_genes, drop = FALSE], 2, function(x) as.numeric(as.character(x))))
colnames(mat_pd_clean) <- df_pd_sorted$Sample_ID

# 🔥【核心修复】：清洗合并数据中的零方差 (Zero-Variance) 或全 NA 基因
mat_combined_test <- cbind(mat_ad_clean, mat_pd_clean)
row_vars <- apply(mat_combined_test, 1, function(x) var(x, na.rm = TRUE))

# 找出不可用于聚类的基因 (方差为 0 或计算出 NA)
invalid_genes <- names(row_vars[is.na(row_vars) | row_vars == 0])

if (length(invalid_genes) > 0) {
  cat(sprintf("⚠️ 自动过滤掉 %d 个零方差/无变化的基因，以防止 hclust 聚类报错：\n", length(invalid_genes)))
  cat(paste(invalid_genes, collapse = ", "), "\n\n")
  
  common_genes <- setdiff(common_genes, invalid_genes)
  mat_ad_clean <- mat_ad_clean[common_genes, , drop = FALSE]
  mat_pd_clean <- mat_pd_clean[common_genes, , drop = FALSE]
}

# 4. 按四个亚型切片
ad_c1 <- mat_ad_clean[, df_ad_sorted$Lysosome_Subtype == "Cluster_1", drop = FALSE]
ad_c2 <- mat_ad_clean[, df_ad_sorted$Lysosome_Subtype == "Cluster_2", drop = FALSE]
pd_c1 <- mat_pd_clean[, df_pd_sorted$Lysosome_Subtype == "Cluster_1", drop = FALSE]
pd_c2 <- mat_pd_clean[, df_pd_sorted$Lysosome_Subtype == "Cluster_2", drop = FALSE]

# 5. 输出样本统计信息
cat("====================================================================\n")
cat("📊 【样本数量与亚群分布统计】\n")
cat("====================================================================\n")
cat(sprintf("🔹 AD 总样本数 : %d 个 (Cluster_1: %d, Cluster_2: %d)\n", ncol(ad_c1) + ncol(ad_c2), ncol(ad_c1), ncol(ad_c2)))
cat(sprintf("🔹 PD 总样本数 : %d 个 (Cluster_1: %d, Cluster_2: %d)\n", ncol(pd_c1) + ncol(pd_c2), ncol(pd_c1), ncol(pd_c2)))
cat("====================================================================\n\n")

# 6. 插入【双柱宽】全 NA 虚拟空隙矩阵
blank_cols1 <- matrix(NA, nrow = length(common_genes), ncol = 2, dimnames = list(common_genes, c("gap1_1", "gap1_2")))
blank_cols2 <- matrix(NA, nrow = length(common_genes), ncol = 2, dimnames = list(common_genes, c("gap2_1", "gap2_2")))
blank_cols3 <- matrix(NA, nrow = length(common_genes), ncol = 2, dimnames = list(common_genes, c("gap3_1", "gap3_2")))

heatmap_mat_final <- cbind(ad_c1, blank_cols1, ad_c2, blank_cols2, pd_c1, blank_cols3, pd_c2)

# 7. 构建注释与配色彩盘
subtypes_combined <- c(
  rep("Cluster_1", ncol(ad_c1)), rep("Blank", 2),
  rep("Cluster_2", ncol(ad_c2)), rep("Blank", 2),
  rep("Cluster_1", ncol(pd_c1)), rep("Blank", 2),
  rep("Cluster_2", ncol(pd_c2))
)

diseases_combined <- c(
  rep("AD", ncol(ad_c1)), rep("Blank", 2),
  rep("AD", ncol(ad_c2)), rep("Blank", 2),
  rep("PD", ncol(pd_c1)), rep("Blank", 2),
  rep("PD", ncol(pd_c2))
)

annotation_col <- data.frame(
  Subtype = as.factor(subtypes_combined),
  Disease = as.factor(diseases_combined)
)
rownames(annotation_col) <- colnames(heatmap_mat_final)

ann_colors <- list(
  Subtype = c("Cluster_1" = "#708090", "Cluster_2" = "#FA8072", "Blank" = "#FFFFFF"),
  Disease = c("AD" = "#1F4E79", "PD" = "#D97706", "Blank" = "#FFFFFF")
)

heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)

# ==============================================================================
# 🖼️ 导出 PDF（同样使用相对路径保存在当前工作目录下）
# ==============================================================================
output_pdf <- "Cross_Disease_Lysosome_Combined_PerfectGaps_v2_Mini.pdf"
cat(sprintf("🎨 正在导出跨疾病合并热图 PDF 至：%s\n", output_pdf))

pdf(
  file        = output_pdf,
  width       = 3.9,
  height      = 3.5,
  family      = "Helvetica",
  useDingbats = FALSE
)

pheatmap(
  mat               = heatmap_mat_final,
  color             = heatmap_colors,
  scale             = "row",
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  na_col            = "white",
  treeheight_row    = 8,
  cellwidth         = NA,
  cellheight        = 14,
  border_color      = NA,
  fontsize          = 8,
  fontsize_row      = 8,
  show_colnames     = FALSE,
  annotation_col    = annotation_col,
  annotation_colors = ann_colors,
  main              = "Cross-Disease Lysosome Molecular Subtypes Overview"
)

dev.off()

cat(sprintf("🎉 运行成功！PDF 已保存至相对路径：%s\n", output_pdf))