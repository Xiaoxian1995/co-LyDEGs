# ==============================================================================
# 🎯 [Lysosome 分子分型热图·单栏微型·相对路径版]
# ==============================================================================
library(pheatmap)
library(dplyr)

# ------------------------------------------------------------------------------
# 0. 相对路径与输出文件设置
# ------------------------------------------------------------------------------
# 采用相对路径设置
input_csv  <- file.path(".", "AD_Lysosome_Subtypes_With_Expression_and_PCA.csv")
pdf_output <- file.path(".", "Lysosome_Molecular_Subtypes_Overview_Mini.pdf")

if (!file.exists(input_csv)) {
  stop(sprintf("⚠️ 错误：未在相对路径下找到数据文件：\n%s", input_csv))
}

cat(sprintf("📊 1. 正在从相对路径载入数据:\n %s\n", input_csv))
df_tmp <- read.csv(input_csv, stringsAsFactors = FALSE, header = TRUE)

# ------------------------------------------------------------------------------
# 1. 数据预处理与样本排序（锁定左 Cluster_1 右 Cluster_2 顺序）
# ------------------------------------------------------------------------------
if ("Lysosome_Subtype" %in% colnames(df_tmp)) {
  
  # 兼容 Sample_ID 列名
  if ("Sample_ID" %in% colnames(df_tmp)) {
    df_tmp$Final_Sample_ID <- as.character(df_tmp$Sample_ID)
  } else {
    df_tmp$Final_Sample_ID <- rownames(df_tmp)
  }
  
  # 按 Subtype 严格排序，确保左 1 右 2 绝对对齐
  df_sorted <- df_tmp[order(df_tmp$Lysosome_Subtype), ]
  extracted_subtype <- as.character(df_sorted$Lysosome_Subtype)
  names(extracted_subtype) <- df_sorted$Final_Sample_ID
  
  # 指定 5 个核心溶酶体基因列
  target_genes <- c("SLC7A14", "ACP5", "SLC2A13", "RAMP3", "CLTC")
  gene_cols    <- intersect(target_genes, colnames(df_sorted))
  
  if (length(gene_cols) == 0) {
    stop("⚠️ 错误：数据集中未找到指定的溶酶体核心基因！")
  }
  
  pure_gene_data <- df_sorted[, gene_cols, drop = FALSE]
  
  # 构建热图矩阵（行为基因，列为样本）
  heatmap_mat_final <- t(apply(pure_gene_data, 2, as.numeric))
  colnames(heatmap_mat_final) <- df_sorted$Final_Sample_ID
  rownames(heatmap_mat_final) <- gene_cols
} else {
  stop("⚠️ 错误：未找到 'Lysosome_Subtype' 列！")
}

# ------------------------------------------------------------------------------
# 2. 1:1 纯正原图色彩系统与分组配色
# ------------------------------------------------------------------------------
annotation_col <- data.frame(Subtype = as.factor(extracted_subtype))
rownames(annotation_col) <- colnames(heatmap_mat_final)

# 顶部 Subtype 条的配色（保持原图：Cluster_1 灰色，Cluster_2 珊瑚红）
ann_colors <- list(
  Subtype = c("Cluster_1" = "#708090", "Cluster_2" = "#FA8072")
)

# 1:1 还原原图干净清澈的经典学术红白蓝
heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)

# ------------------------------------------------------------------------------
# 3. 极致压缩排版与绘图（物理宽度削减 1/3）
# ------------------------------------------------------------------------------
cat("🎨 2. 正在导出单栏微型热图 PDF...\n")

# 将宽度缩减至 3.7，高度维持 3.5
pdf(pdf_output, width = 3.7, height = 3.5)

pheatmap(
  mat               = heatmap_mat_final,
  color             = heatmap_colors,
  scale             = "row",            # 按行进行 Z-score 标准化
  cluster_rows      = TRUE,             # 纵向基因保持自动聚类
  cluster_cols      = FALSE,            # 隐藏样本聚类树，锁死左 1 右 2 顺序
  
  # 细节微调：左侧树线极致贴边
  treeheight_row    = 8,                # 缩到 8，贴边不占地方
  
  # 单栏微型尺寸排版
  cellwidth         = NA,               # 样本列在 3.7 宽度内自动挤压
  cellheight        = 14,               # 纵向格子高度 14
  border_color      = NA,               # 关闭格子边框
  
  # 字号死守 8 号
  fontsize          = 8,                
  fontsize_row      = 8,                # 右侧基因字号 8 号
  show_colnames     = FALSE,            # 隐藏样本 ID
  
  annotation_col    = annotation_col,   
  annotation_colors = ann_colors,       
  
  main              = "Lysosome Molecular Subtypes Overview in AD"
)

dev.off()

cat("\n🎉 【大功告成！】\n")
cat(sprintf("💾 最终版 PDF 已成功保存至相对路径：\n %s\n", pdf_output))