# ==============================================================================
# 🎯 [Lysosome 分子分型热图·PD数据集·保留聚类树+指定基因顺序·相对路径版]
# ==============================================================================
library(pheatmap)
library(dplyr)

# ------------------------------------------------------------------------------
# 0. 相对路径与输出文件设置
# ------------------------------------------------------------------------------
# 采用相对路径，直接从当前工作目录读取和保存
input_csv  <- file.path(".", "PD_Lysosome_Subtypes_With_Expression_and_PCA.csv")
pdf_output <- file.path(".", "PD_Lysosome_Molecular_Subtypes_Overview_Mini.pdf")

if (!file.exists(input_csv)) {
  stop(sprintf("⚠️ 错误：未在相对路径下找到 PD 数据文件：\n%s", input_csv))
}

cat(sprintf("📊 1. 正在从相对路径载入 PD 数据:\n %s\n", input_csv))
df_tmp <- read.csv(input_csv, stringsAsFactors = FALSE, header = TRUE)

# ------------------------------------------------------------------------------
# 1. 数据预处理与样本排序
# ------------------------------------------------------------------------------
if ("Lysosome_Subtype" %in% colnames(df_tmp)) {
  
  if ("Sample_ID" %in% colnames(df_tmp)) {
    df_tmp$Final_Sample_ID <- as.character(df_tmp$Sample_ID)
  } else {
    df_tmp$Final_Sample_ID <- rownames(df_tmp)
  }
  
  # 按 Subtype 严格排序，确保左 Cluster_1 右 Cluster_2 绝对对齐
  df_sorted <- df_tmp[order(df_tmp$Lysosome_Subtype), ]
  extracted_subtype <- as.character(df_sorted$Lysosome_Subtype)
  names(extracted_subtype) <- df_sorted$Final_Sample_ID
  
  # 🎯 你的目标基因上下顺序（从上到下）：
  target_genes_order <- c("RAMP3", "ACP5", "SLC7A14", "CLTC", "SLC2A13")
  gene_cols          <- intersect(target_genes_order, colnames(df_sorted))
  
  if (length(gene_cols) == 0) {
    stop("⚠️ 错误：数据集中未找到指定的溶酶体核心基因！")
  }
  
  # 按照指定顺序构建表达数据子集
  pure_gene_data <- df_sorted[, gene_cols, drop = FALSE]
  
  # 构建热图矩阵（行为基因，列为样本）
  heatmap_mat_final <- t(apply(pure_gene_data, 2, as.numeric))
  colnames(heatmap_mat_final) <- df_sorted$Final_Sample_ID
  rownames(heatmap_mat_final) <- gene_cols
} else {
  stop("⚠️ 错误：未找到 'Lysosome_Subtype' 列！")
}

# ------------------------------------------------------------------------------
# 2. 计算 Z-score 并在保留聚类树的前提下强制重排树枝叶节点
# ------------------------------------------------------------------------------
# 按行标准化（与 pheatmap 的 scale="row" 逻辑一致）
mat_scaled <- t(scale(t(heatmap_mat_final)))

# 计算行距离并建立 hclust 对象
row_dist   <- dist(mat_scaled)
row_hclust <- hclust(row_dist, method = "complete")

# 转化为 dendrogram 并按指定的 target_genes_order 进行重新排序
row_dend <- as.dendrogram(row_hclust)

# 匹配并重排 dendrogram 节点顺序
desired_order <- match(target_genes_order, rownames(heatmap_mat_final))
row_dend_reordered <- reorder(row_dend, desired_order, aggregateFun = mean)

# 转换回 hclust 传给 pheatmap
row_hclust_final <- as.hclust(row_dend_reordered)

# ------------------------------------------------------------------------------
# 3. 色彩系统与分组配色
# ------------------------------------------------------------------------------
annotation_col <- data.frame(Subtype = as.factor(extracted_subtype))
rownames(annotation_col) <- colnames(heatmap_mat_final)

ann_colors <- list(
  Subtype = c("Cluster_1" = "#708090", "Cluster_2" = "#FA8072")
)

heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)

# ------------------------------------------------------------------------------
# 4. 绘图导出（聚类树保留，且基因严格按指定顺序排列）
# ------------------------------------------------------------------------------
cat("🎨 2. 正在导出带聚类树且按指定顺序排列的 PD 单栏微型热图 PDF...\n")

pdf(pdf_output, width = 3.7, height = 3.5)

pheatmap(
  mat               = heatmap_mat_final,
  color             = heatmap_colors,
  scale             = "row",            # 按行进行 Z-score 标准化
  cluster_rows      = row_hclust_final, # 🌲 传入重排后的聚类树对象（既保留树结构，又锁定基因从上到下顺序）
  cluster_cols      = FALSE,            # 隐藏样本聚类树，锁死左 1 右 2 顺序
  
  treeheight_row    = 8,                # 左侧树线极致贴边
  cellwidth         = NA,               # 样本列在 3.7 宽度内自动挤压
  cellheight        = 14,               # 纵向格子高度 14
  border_color      = NA,               # 关闭格子边框
  
  fontsize          = 8,                
  fontsize_row      = 8,                # 右侧基因字号 8 号
  show_colnames     = FALSE,            # 隐藏样本 ID
  
  annotation_col    = annotation_col,   
  annotation_colors = ann_colors,       
  
  main              = "Lysosome Molecular Subtypes Overview in PD"
)

dev.off()

cat("\n🎉 【大功告成！】\n")
cat(sprintf("💾 带有聚类树且按指定基因顺序（RAMP3->ACP5->SLC7A14->CLTC->SLC2A13）排列的最终版 PDF 已成功保存至相对路径：\n %s\n", pdf_output))