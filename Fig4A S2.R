library(ConsensusClusterPlus)
library(dplyr)
library(tidyr)
library(ggplot2)

# ==============================================================================
# 0. 相对路径与输出文件夹设置
# ==============================================================================
# 设定相对输出目录（保存在当前工作目录下的 PD_Consensus_Clustering_Result 文件夹）
out_dir <- file.path(".", "PD_Consensus_Clustering_Result")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ==============================================================================
# 1. 数据加载与预处理：严格筛选出 PD 患者样本 (相对路径加载)
# ==============================================================================
csv_file <- file.path(".", "PD_5_Genes_Sample_Expression.csv")

if (!file.exists(csv_file)) {
  stop("❌ 未在当前工作目录下找到基因表达数据文件 'PD_5_Genes_Sample_Expression.csv'，请检查路径！")
}

cat(sprintf("📂 成功找到 PD 输入文件 (相对路径): %s\n", csv_file))
all_data <- read.csv(csv_file, stringsAsFactors = FALSE, header = TRUE)

# 自动识别 Group 和 Sample_ID 列名（不区分大小写）
group_col  <- colnames(all_data)[grep("^group$", colnames(all_data), ignore.case = TRUE)][1]
sample_col <- colnames(all_data)[grep("sample", colnames(all_data), ignore.case = TRUE)][1]

# 过滤并仅保留 PD 分组患者
pd_data <- all_data %>% filter(tolower(.data[[group_col]]) == "pd")

# 提取 5 个核心溶酶体基因
cluster_genes <- c("SLC7A14", "ACP5", "SLC2A13", "RAMP3", "CLTC")

# 构建聚类所需的标准矩阵格式 (行为基因，列为样本 ID)
expr_matrix_pd <- pd_data %>%
  select(all_of(sample_col), all_of(cluster_genes)) %>%
  tibble::column_to_rownames(sample_col) %>%
  as.matrix() %>%
  t()

# 执行行中心化 (Mean Centering)
expr_matrix_centered <- t(scale(t(expr_matrix_pd), scale = FALSE))

# ==============================================================================
# 2. 运行一致性聚类 (Consensus Clustering)
# ==============================================================================
cat("\n====================================================================\n")
cat("🧬 正在运行 ConsensusClusterPlus (hc + pearson) [PD Cohort]... \n")
cat("--------------------------------------------------------------------\n")

set.seed(123456) # 固定随机种子，保证每次运行结果 100% 可重复

results <- ConsensusClusterPlus(
  d = expr_matrix_centered,
  maxK = 6,
  reps = 1000,
  pItem = 0.8,
  pFeature = 1,
  title = out_dir,          # 输出 Consensus 图表到相对目录
  clusterAlg = "hc",        # 层次聚类
  distance = "pearson",    # 皮尔逊相关距离
  seed = 123456,
  plot = "pdf"
)

# 提取 K = 2 时的亚型划分结果
best_k <- 2
subtype_assignments <- results[[best_k]][["consensusClass"]]

# 构建分型结果表
pd_subtypes <- data.frame(
  Sample_ID = names(subtype_assignments),
  Lysosome_Subtype = paste0("Cluster_", subtype_assignments),
  stringsAsFactors = FALSE
)

# 保存单纯亚型结果表 (相对路径)
subtypes_only_csv <- file.path(out_dir, "PD_Patients_Lysosome_Subtypes_Result.csv")
write.csv(pd_subtypes, subtypes_only_csv, row.names = FALSE)

# ==============================================================================
# 3. 执行 PCA 降维并映射 Consensus 亚群
# ==============================================================================
pca_expr_mat <- t(expr_matrix_pd)
pca_res <- prcomp(pca_expr_mat, scale. = TRUE)
pca_coordinates <- as.data.frame(pca_res$x)

pc1_var <- round(summary(pca_res)$importance[2, 1] * 100, 2)
pc2_var <- round(summary(pca_res)$importance[2, 2] * 100, 2)

# 映射聚类结果
subtype_map <- setNames(pd_subtypes$Lysosome_Subtype, pd_subtypes$Sample_ID)
pca_coordinates$Subtype <- factor(subtype_map[rownames(pca_coordinates)])

# 统计样本量
c1_n <- sum(pca_coordinates$Subtype %in% c("Cluster_1", "Cluster 1"))
c2_n <- sum(pca_coordinates$Subtype %in% c("Cluster_2", "Cluster 2"))

cat(sprintf("✅ 聚类及 PCA 匹配完成！样本分布：Sample Size: %d (Cluster_1: %d, Cluster_2: %d)\n", 
            nrow(pca_coordinates), c1_n, c2_n))

# ==============================================================================
# 4. 绘制高颜值 PCA 散点图 + 2 亚型置信椭圆并导出 PDF (相对路径)
# ==============================================================================
p_native_pca <- ggplot(pca_coordinates, aes(x = PC1, y = PC2, color = Subtype, fill = Subtype)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill = Subtype), level = 0.95) +
  scale_color_manual(values = c("Cluster_1" = "#6B7280", "Cluster_2" = "#F87171",
                                "Cluster 1" = "#6B7280", "Cluster 2" = "#F87171")) +
  scale_fill_manual(values = c("Cluster_1" = "#6B7280", "Cluster_2" = "#F87171",
                               "Cluster 1" = "#6B7280", "Cluster 2" = "#F87171")) +
  theme_bw() +
  labs(
    title = "PCA of Lysosome-Related Subtypes (PD Cohort)",
    subtitle = paste0("Sample Size: ", nrow(pca_coordinates), " (Cluster_1: ", c1_n, ", Cluster_2: ", c2_n, ")"),
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    color = "Lysosome Subtype",
    fill = "Lysosome Subtype"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "gray30", size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

pdf_path <- file.path(out_dir, "PCA_PD_Lysosome_Consensus_ExactReproduce.pdf")
pdf(pdf_path, width = 6.5, height = 5.5)
print(p_native_pca)
dev.off()

# ==============================================================================
# 5. 整合样本亚型、基因表达量与 PCA 主成分得分并保存本地 (相对路径)
# ==============================================================================
pca_scores <- as.data.frame(pca_res$x)

final_combined_df <- pd_data %>%
  select(all_of(sample_col), all_of(group_col), all_of(cluster_genes)) %>%
  left_join(pd_subtypes, by = setNames("Sample_ID", sample_col)) %>%
  cbind(pca_scores[match(pd_data[[sample_col]], rownames(pca_scores)), ])

colnames(final_combined_df)[colnames(final_combined_df) == sample_col] <- "Sample_ID"
colnames(final_combined_df)[colnames(final_combined_df) == group_col]  <- "Group"

final_combined_df <- final_combined_df %>%
  select(
    Sample_ID, 
    Group, 
    Lysosome_Subtype, 
    starts_with("PC"), 
    all_of(cluster_genes)
  )

output_csv <- file.path(out_dir, "PD_Lysosome_Subtypes_With_Expression_and_PCA.csv")
write.csv(final_combined_df, output_csv, row.names = FALSE)

cat("\n====================================================================\n")
cat("🎉 【PD 队列分析完成】！全部文件已成功相对保存至目标文件夹：\n")
cat(sprintf(" 1. PDF 矢量图表 : '%s'\n", pdf_path))
cat(sprintf(" 2. 整合数据表格 : '%s'\n", output_csv))
cat("====================================================================\n")

# 预览前 5 行整合结果
print(head(final_combined_df, 5))