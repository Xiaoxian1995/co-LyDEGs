# ==============================================================================
# 🎯 [AD 免疫分析一体化脚本] 热图 + Wilcoxon分面箱线图 + 全差异排序CSV (相对路径版)
# ==============================================================================
library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

cat("📊 1. 正在读取与处理 AD 矩阵数据 (使用相对路径)...\n")

# ------------------------------------------------------------------------------
# 0. 相对路径定义与数据读取
# ------------------------------------------------------------------------------
ad_csv_path     <- "AD_Lysosome_Subtypes_With_Expression_and_PCA.csv"
ssgsea_csv_path <- "ad_ssgsea_23_scores.csv"

if (!file.exists(ad_csv_path)) {
  stop(sprintf("⚠️ 错误：未在当前目录下找到 AD 亚型数据文件：%s", ad_csv_path))
}
if (!file.exists(ssgsea_csv_path)) {
  stop(sprintf("⚠️ 错误：未在当前目录下找到 ssGSEA 免疫得分数据文件：%s", ssgsea_csv_path))
}

df_ad_group_raw  <- read.csv(ad_csv_path, stringsAsFactors = FALSE, header = TRUE)
df_ssgsea_scores <- read.csv(ssgsea_csv_path, stringsAsFactors = FALSE, header = TRUE)

# ---- [第一步：提取并规范样本 ID] ----
df_ad_group <- df_ad_group_raw
if ("Sample_ID" %in% colnames(df_ad_group)) {
  df_ad_group$Final_Sample_ID <- as.character(df_ad_group$Sample_ID)
} else if ("X" %in% colnames(df_ad_group)) {
  df_ad_group$Final_Sample_ID <- as.character(df_ad_group$X)
} else {
  df_ad_group$Final_Sample_ID <- rownames(df_ad_group)
}

# ---- [第二步：清洗 ssGSEA 表并构建纯数值矩阵] ----
df_immune <- df_ssgsea_scores

if (is.character(df_immune[[1]]) || is.factor(df_immune[[1]])) {
  rownames(df_immune) <- df_immune[[1]]
  df_immune <- df_immune[, -1, drop = FALSE]
} else if ("Sample_ID" %in% colnames(df_immune)) {
  rownames(df_immune) <- df_immune$Sample_ID
  df_immune$Sample_ID <- NULL
} else if ("X" %in% colnames(df_immune)) {
  rownames(df_immune) <- df_immune$X
  df_immune$X <- NULL
}

if (any(df_ad_group$Final_Sample_ID %in% rownames(df_immune))) {
  mat_cells_pure <- t(as.matrix(df_immune))
} else if (any(df_ad_group$Final_Sample_ID %in% colnames(df_immune))) {
  mat_cells_pure <- as.matrix(df_immune)
} else {
  mat_cells_pure <- as.matrix(df_immune)
}

class(mat_cells_pure) <- "numeric"

# ---- [第三步：样本对齐与横向排序] ----
common_samples <- intersect(df_ad_group$Final_Sample_ID, colnames(mat_cells_pure))

if (length(common_samples) == 0) {
  mat_cells_pure <- t(mat_cells_pure)
  common_samples <- intersect(df_ad_group$Final_Sample_ID, colnames(mat_cells_pure))
}

if (length(common_samples) == 0) {
  stop("⚠️ 错误：AD 亚型表与 ssGSEA 免疫表的样本 ID 无法匹配！")
}

df_ad_aligned <- df_ad_group[df_ad_group$Final_Sample_ID %in% common_samples, ]
df_ad_aligned <- df_ad_aligned[order(df_ad_aligned$Lysosome_Subtype), ]
extracted_subtype <- as.character(df_ad_aligned$Lysosome_Subtype)

mat_aligned <- mat_cells_pure[, df_ad_aligned$Final_Sample_ID, drop = FALSE]

# ------------------------------------------------------------------------------
# 🎨 1. 生成并导出 AI 可编辑的 AD 免疫热图 PDF
# ------------------------------------------------------------------------------
cat("🎨 2. 正在导出 AD 免疫热图 PDF...\n")

# 数据清洗 (处理 zero variance & NAs)
row_vars <- apply(mat_aligned, 1, var, na.rm = TRUE)
zero_var_rows <- which(row_vars == 0 | is.na(row_vars))
heatmap_mat_final <- if (length(zero_var_rows) > 0) mat_aligned[-zero_var_rows, , drop = FALSE] else mat_aligned

heatmap_mat_scaled <- t(apply(heatmap_mat_final, 1, scale))
colnames(heatmap_mat_scaled) <- colnames(heatmap_mat_final)
rownames(heatmap_mat_scaled) <- rownames(heatmap_mat_final)

if (any(is.na(heatmap_mat_scaled) | is.nan(heatmap_mat_scaled) | is.infinite(heatmap_mat_scaled))) {
  for (i in 1:nrow(heatmap_mat_scaled)) {
    bad_idx <- which(is.na(heatmap_mat_scaled[i,]) | is.nan(heatmap_mat_scaled[i,]) | is.infinite(heatmap_mat_scaled[i,]))
    if (length(bad_idx) > 0) {
      valid_vals <- heatmap_mat_scaled[i, -bad_idx]
      fill_val <- if(length(valid_vals) > 0) median(valid_vals, na.rm = TRUE) else 0
      heatmap_mat_scaled[i, bad_idx] <- fill_val
    }
  }
}

annotation_col <- data.frame(Subtype = as.factor(extracted_subtype))
rownames(annotation_col) <- colnames(heatmap_mat_scaled)

subtype_levels <- levels(annotation_col$Subtype)
color_palette  <- c("#708090", "#FA8072")
ann_colors     <- list(Subtype = setNames(color_palette[1:length(subtype_levels)], subtype_levels))
heatmap_colors <- colorRampPalette(c("#1E466E", "#6BAED6", "#F7FBFF", "#E6550D", "#A63603"))(100)

output_heatmap_pdf <- "AD_Lysosome_Subtypes_Immune_Infiltration_AIVersion.pdf"

pdf(file = output_heatmap_pdf, width = 5.04, height = 2.8, family = "Helvetica", useDingbats = FALSE)
pheatmap(
  mat               = heatmap_mat_scaled,
  color             = heatmap_colors,
  scale             = "none",
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  cellwidth         = NA,
  cellheight        = NA,
  treeheight_row    = 8,
  border_color      = NA,
  fontsize          = 5.5,
  fontsize_row      = 5.5,
  show_colnames     = FALSE,
  annotation_col    = annotation_col,
  annotation_colors = ann_colors,
  annotation_legend = TRUE,
  main              = "Immune Infiltration Landscape (AD Subtypes)"
)
dev.off()

# ------------------------------------------------------------------------------
# 🧪 2. 全免疫细胞 Wilcoxon 检验与差异全排序 CSV 导出
# ------------------------------------------------------------------------------
cat("🧪 3. 正在计算所有免疫细胞的 Wilcoxon 检验差异并导出 CSV...\n")

mat_df <- as.data.frame(t(mat_aligned))
mat_df$Final_Sample_ID <- rownames(mat_df)

merged_df <- inner_join(
  df_ad_aligned %>% select(Final_Sample_ID, Lysosome_Subtype),
  mat_df,
  by = "Final_Sample_ID"
) %>% filter(!is.na(Lysosome_Subtype))

merged_df$Lysosome_Subtype <- as.factor(merged_df$Lysosome_Subtype)
cell_types <- rownames(mat_aligned)

diff_results <- data.frame()

for (cell in cell_types) {
  sub_data <- merged_df[, c("Lysosome_Subtype", cell)]
  colnames(sub_data) <- c("Subtype", "Score")
  sub_data$Score     <- as.numeric(sub_data$Score)
  sub_data           <- sub_data %>% filter(!is.na(Score) & !is.na(Subtype))
  
  active_levels <- unique(as.character(sub_data$Subtype))
  
  if (length(active_levels) == 2) {
    wilcox_res  <- wilcox.test(Score ~ Subtype, data = sub_data)
    p_val       <- wilcox_res$p.value
    mean_group1 <- mean(sub_data$Score[sub_data$Subtype == active_levels[1]], na.rm = TRUE)
    mean_group2 <- mean(sub_data$Score[sub_data$Subtype == active_levels[2]], na.rm = TRUE)
  } else {
    p_val       <- NA
    mean_group1 <- NA
    mean_group2 <- NA
  }
  
  diff_results <- rbind(diff_results, data.frame(
    Cell_Type     = cell,
    P_Value       = p_val,
    Mean_Cluster1 = mean_group1,
    Mean_Cluster2 = mean_group2,
    Score_Diff    = mean_group2 - mean_group1
  ))
}

diff_results <- diff_results %>%
  filter(!is.na(P_Value)) %>%
  mutate(FDR_q_value = p.adjust(P_Value, method = "BH")) %>%
  arrange(P_Value) %>%
  mutate(
    Rank = row_number(),
    Is_Significant_0.05 = ifelse(P_Value < 0.05, "Yes", "No")
  )

output_csv_path <- "AD_All_Immune_Cells_Differential_Ranked.csv"
write.csv(diff_results, file = output_csv_path, row.names = FALSE)

# ------------------------------------------------------------------------------
# 🎨 3. 动态筛选显著细胞并导出 Wilcoxon 分面小提琴箱线图 PDF
# ------------------------------------------------------------------------------
cat("🎨 4. 正在生成与导出 top 显著免疫细胞分面图 PDF...\n")

sig_df <- diff_results %>% filter(P_Value < 0.05)
n_sig  <- nrow(sig_df)

if (n_sig == 0) {
  sig_cells <- diff_results$Cell_Type[1:min(3, nrow(diff_results))]
} else if (n_sig <= 5) {
  sig_cells <- sig_df$Cell_Type
} else {
  sig_cells <- sig_df$Cell_Type[1:5]
}

p_values_map <- setNames(diff_results$P_Value, diff_results$Cell_Type)

df_plot <- merged_df %>%
  select(Lysosome_Subtype, all_of(sig_cells)) %>%
  pivot_longer(cols = -Lysosome_Subtype, names_to = "Immune_Cell", values_to = "Infiltration_Score") %>%
  mutate(Infiltration_Score = as.numeric(Infiltration_Score)) %>%
  filter(!is.na(Infiltration_Score))

cluster_levels <- levels(factor(df_plot$Lysosome_Subtype))
plot_list <- list()

for (i in 1:length(sig_cells)) {
  cell_name <- sig_cells[i]
  df_sub    <- df_plot %>% filter(Immune_Cell == cell_name)
  p_val     <- p_values_map[[cell_name]]
  
  if (p_val < 0.0001) {
    p_formatted <- format(p_val, digits = 2, scientific = TRUE)
  } else {
    p_formatted <- format(p_val, digits = 2, scientific = FALSE)
  }
  
  p_label   <- paste0("P = ", trimws(p_formatted))
  color_map <- setNames(c("#708090", "#FA8072")[1:length(cluster_levels)], cluster_levels)
  fill_map  <- setNames(c("#B2BBC6", "#F1948A")[1:length(cluster_levels)], cluster_levels)
  
  p <- ggplot(df_sub, aes(x = Lysosome_Subtype, y = Infiltration_Score)) +
    geom_violin(aes(fill = Lysosome_Subtype), alpha = 0.15, color = NA, width = 0.8, adjust = 1.2) +
    geom_jitter(aes(color = Lysosome_Subtype), shape = 16, size = 1.2, width = 0.18, alpha = 0.85) +
    geom_boxplot(color = "#2C3E50", fill = NA, width = 0.3, outlier.shape = NA, lwd = 0.8) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    scale_color_manual(values = color_map) +
    scale_fill_manual(values = fill_map) +
    labs(x = NULL, y = if (i == 1) "Infiltration Score" else NULL) +
    annotate("text", x = 1.5, y = max(df_sub$Infiltration_Score, na.rm = TRUE) * 1.05,
             label = p_label, size = 3.8, fontface = "bold", family = "sans", color = "black") +
    facet_grid(~ Immune_Cell) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(colour = "#EAEAEA", size = 0.5),
      panel.border       = element_rect(colour = "black", fill = NA, size = 1.0),
      strip.background   = element_rect(fill = "#D6DBDF", colour = "black", size = 1.0),
      strip.text         = element_text(size = 10, face = "bold", color = "black", family = "sans"),
      axis.text.x        = element_text(size = 11, color = "black", family = "sans"),
      axis.text.y        = element_text(size = 10, color = "black", family = "sans"),
      axis.title.y       = element_text(size = 12, face = "bold", color = "black", family = "sans"),
      legend.position    = "none"
    )
  
  plot_list[[cell_name]] <- p
}

canvas_width    <- length(sig_cells) * 3.0 + 1.0
output_boxplot_pdf <- "Figure_AD_Immune_Infiltration_SolidPoints_AI.pdf"

pdf(output_boxplot_pdf, width = canvas_width, height = 5.0, family = "sans", useDingbats = FALSE)
print(
  wrap_plots(plot_list, nrow = 1) +
    plot_annotation(
      title = "Top Significant Immune Cell Infiltration In AD Subtypes",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5, family = "sans"))
    )
)
dev.off()

cat("\n====================================================================\n")
cat("🎉 【运行完成！所有输出文件已保存在当前工作目录下】\n")
cat(sprintf("1️⃣ 矢量热图 PDF： %s\n", output_heatmap_pdf))
cat(sprintf("2️⃣ Wilcoxon 分面图 PDF：%s\n", output_boxplot_pdf))
cat(sprintf("3️⃣ 差异全排序 CSV：  %s\n", output_csv_path))
cat("====================================================================\n")