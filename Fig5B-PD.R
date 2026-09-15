# ==============================================================================
# 🎯 [PD 双亚型特征评分 ROC 曲线 - 相对路径版]
# ==============================================================================
library(pROC)
library(ggplot2)

cat("📊 1. 正在读取 PD 本地 CSV 数据并计算 Signature ROC...\n")

# 使用纯相对路径获取当前工作目录下的文件
pd_csv_path <- "PD_Lysosome_Subtypes_With_Expression_and_PCA.csv"

if (!file.exists(pd_csv_path)) {
  stop("⚠️ 错误：未在当前工作目录下找到 PD_Lysosome_Subtypes_With_Expression_and_PCA.csv 文件！")
}

# 1. 载入原始数据
pd_roc_data <- read.csv(pd_csv_path, stringsAsFactors = FALSE, header = TRUE)

# 2. 自动构建二分类响应变量 (Group_Binary)
subtype_col <- NULL
for (col in c("Lysosome_Subtype", "Subtype", "Cluster")) {
  if (col %in% colnames(pd_roc_data)) {
    subtype_col <- col
    break
  }
}

if (is.null(subtype_col)) {
  stop("⚠️ 错误：未能在 CSV 中找到亚型分类列 (如 Lysosome_Subtype 或 Subtype)！")
}

raw_subtypes <- as.character(pd_roc_data[[subtype_col]])
pd_roc_data$Group_Binary <- ifelse(grepl("1", raw_subtypes), 1, 0)

cat(sprintf("💡 已成功将 '%s' 列转化为二分类标签 (样本总数: %d)\n", subtype_col, nrow(pd_roc_data)))

# 3. 🔍 定义更新后的基因分组
genes_A <- c("ACP5", "RAMP3")
genes_B <- c("SLC7A14", "CLTC", "SLC2A13")

# 检查基因是否存在于列名中
missing_A <- setdiff(genes_A, colnames(pd_roc_data))
missing_B <- setdiff(genes_B, colnames(pd_roc_data))
if (length(missing_A) > 0 || length(missing_B) > 0) {
  warning("⚠️ 提示：部分指定基因在数据框列名中未完全找到，请核对拼写！")
}

# 4. ⚡️ 自动在数据框中计算 Signature 评分（求行均值）
valid_genes_A <- intersect(genes_A, colnames(pd_roc_data))
valid_genes_B <- intersect(genes_B, colnames(pd_roc_data))

pd_roc_data$PD_Signature_A <- rowMeans(pd_roc_data[, valid_genes_A, drop = FALSE], na.rm = TRUE)
pd_roc_data$PD_Signature_B <- rowMeans(pd_roc_data[, valid_genes_B, drop = FALSE], na.rm = TRUE)

# 5. 📈 计算对应的 ROC 曲线
roc_A_pd <- roc(response = pd_roc_data$Group_Binary, 
                predictor = pd_roc_data$PD_Signature_A, 
                quiet = TRUE)

roc_B_pd <- roc(response = pd_roc_data$Group_Binary, 
                predictor = pd_roc_data$PD_Signature_B, 
                quiet = TRUE)

# 6. 提取 ROC 坐标点并转化为标准 FPR (1 - Specificity)
coords_A <- coords(roc_A_pd, "all", ret = c("specificity", "sensitivity"), transpose = FALSE)
coords_B <- coords(roc_B_pd, "all", ret = c("specificity", "sensitivity"), transpose = FALSE)

df_A <- data.frame(
  FPR = 1 - coords_A$specificity,
  TPR = coords_A$sensitivity,
  Signature = "Signature A (ACP5, RAMP3)"
)
df_A <- rbind(data.frame(FPR = 0, TPR = 0, Signature = df_A$Signature[1]), 
              df_A, 
              data.frame(FPR = 1, TPR = 1, Signature = df_A$Signature[1]))
df_A <- df_A[order(df_A$FPR, df_A$TPR), ]

df_B <- data.frame(
  FPR = 1 - coords_B$specificity,
  TPR = coords_B$sensitivity,
  Signature = "Signature B (SLC7A14, CLTC, SLC2A13)"
)
df_B <- rbind(data.frame(FPR = 0, TPR = 0, Signature = df_B$Signature[1]), 
              df_B, 
              data.frame(FPR = 1, TPR = 1, Signature = df_B$Signature[1]))
df_B <- df_B[order(df_B$FPR, df_B$TPR), ]

df_roc_merged <- rbind(df_A, df_B)

# 获取精确 AUC 值
auc_A_val <- sprintf("%.4f", as.numeric(auc(roc_A_pd))) 
auc_B_val <- sprintf("%.4f", as.numeric(auc(roc_B_pd))) 

cat("🎨 2. 正在生成 PD 学术级 ROC 曲线图...\n")

# 7. ggplot2 绘图
p_roc_pd <- ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               color = "#B2BABB", linetype = "dashed", linewidth = 0.6) +
  geom_step(data = df_roc_merged, aes(x = FPR, y = TPR, color = Signature), 
            linewidth = 1.1, direction = "vh") +
  scale_color_manual(values = c(
    "Signature A (ACP5, RAMP3)" = "#C87A2F",
    "Signature B (SLC7A14, CLTC, SLC2A13)" = "#1F4E79"
  )) +
  scale_x_continuous(limits = c(-0.01, 1.01), expand = c(0, 0), breaks = seq(0, 1, 0.2), labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0")) +
  scale_y_continuous(limits = c(-0.01, 1.01), expand = c(0, 0), breaks = seq(0, 1, 0.2), labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0")) +
  labs(
    title = "PD",
    x = "False Positive Rate (1-Specificity)", 
    y = "Sensitivity"
  ) +
  annotate("text", x = 0.53, y = 0.24, 
           label = paste0("Signature A AUC = ", auc_A_val), 
           color = "#C87A2F", size = 4.2, fontface = "bold", hjust = 0, family = "Helvetica") +
  annotate("text", x = 0.53, y = 0.16, 
           label = paste0("Signature B AUC = ", auc_B_val), 
           color = "#1F4E79", size = 4.2, fontface = "bold", hjust = 0, family = "Helvetica") +
  theme_classic() +
  theme(
    plot.title      = element_text(size = 17, face = "bold", hjust = 0.5, family = "Helvetica", margin = margin(b=10)),
    axis.line       = element_line(color = "black", linewidth = 0.8),
    axis.ticks      = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length = unit(0.18, "cm"),
    axis.text.x     = element_text(size = 14, color = "black", family = "Helvetica", margin = margin(t=5)),
    axis.text.y     = element_text(size = 14, color = "black", family = "Helvetica", margin = margin(r=5)),
    axis.title.x    = element_text(size = 15, face = "bold", color = "black", family = "Helvetica", margin = margin(t=8)),
    axis.title.y    = element_text(size = 15, face = "bold", color = "black", family = "Helvetica", margin = margin(r=8)),
    legend.position = "none"
  )

# ------------------------------------------------------------------------------
# 3. 💾 矢量 PDF 导出（使用相对路径保存到当前工作目录）
# ------------------------------------------------------------------------------
output_pdf <- "Figure_PD_Merged_Signatures_ROC_Editable.pdf"
cat(sprintf("💾 3. 正在生成高保真可编辑的 PD 矢量 PDF：%s\n", output_pdf))

pdf(output_pdf, width = 5.2, height = 5.2, useDingbats = FALSE, family = "Helvetica")
print(p_roc_pd)
dev.off()

cat(sprintf("\n🎉 【PD 双特征 ROC 曲线生成成功！】\n💾 相对路径保存文件：%s\n", output_pdf))