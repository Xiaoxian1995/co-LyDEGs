# ==============================================================================
# 🎯 [GSE33000 5基因随机森林特征权重精修矢量图一键通关脚本]
# ==============================================================================
# 📝 备注：运行本脚本前，请提前下载以下 GEO 数据文件并放置于脚本同级目录（./）：
# 
# 1. GSE33000 表达矩阵文件 (GSE33000_series_matrix.txt.gz):
#    下载链接: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE33nnn/GSE33000/matrix/GSE33000_series_matrix.txt.gz
# 
# 2. GPL4372 平台注释文件 (GPL4372.gz):
#    下载链接: https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL4nnn/GPL4372/suppl/GPL4372.gz
# ==============================================================================

library(GEOquery)
library(Biobase)
library(dplyr)
library(tidyr)
library(randomForest)
library(ggplot2)

base_dir <- "./"
matrix_file <- file.path(base_dir, "GSE33000_series_matrix.txt.gz")
gpl_file    <- file.path(base_dir, "GPL4372.gz")

# ------------------------------------------------------------------------------
# 1. 加载 GSE33000 数据与 GPL4372 注释
# ------------------------------------------------------------------------------
cat("📂 1. 正在读取 GSE33000 数据与 GPL4372 注释...\n")
ext_gse <- getGEO(filename = matrix_file, GSEMatrix = TRUE, getGPL = FALSE)
gpl_data <- getGEO(filename = gpl_file)

# 把 GPL4372 的注释嵌入到 fData 中
fData(ext_gse) <- Table(gpl_data)

probe_annotation  <- fData(ext_gse)
raw_matrix        <- exprs(ext_gse)
ext_clinical_info <- pData(ext_gse)

# ------------------------------------------------------------------------------
# 2. 精准匹配探针并提取 5 个核心基因表达量
# ------------------------------------------------------------------------------
cat("🛠️ 2. 正在匹配探针并提取 5 个核心基因表达量...\n")

target_genes <- c("CLTC", "RAMP3", "SLC7A14", "ACP5", "SLC2A13")

# 寻找注释表中真正存放 Gene Symbol 的列
potential_cols <- colnames(probe_annotation)
best_col <- NULL
for(col in potential_cols) {
  test_vals <- as.character(probe_annotation[[col]])
  if(any(c("CLTC", "RAMP3", "ACP5") %in% test_vals)) {
    best_col <- col
    break
  }
}

val_gene_expr <- matrix(0, nrow = length(target_genes), ncol = ncol(raw_matrix))
rownames(val_gene_expr) <- target_genes
colnames(val_gene_expr) <- colnames(raw_matrix)

if(!is.null(best_col)) {
  message(paste("🎯 找到芯片注释列:", best_col))
  probe_id_col <- colnames(probe_annotation)[1]
  for (gene in target_genes) {
    matched_rows <- which(probe_annotation[[best_col]] == gene)
    matched_probes <- as.character(probe_annotation[[probe_id_col]][matched_rows])
    valid_probes <- intersect(matched_probes, rownames(raw_matrix))
    
    if (length(valid_probes) > 0) {
      val_gene_expr[gene, ] <- colMeans(raw_matrix[valid_probes, , drop = FALSE], na.rm = TRUE)
      cat(sprintf("  ├─ 基因 %-8s : 成功提取到 %d 个探针\n", gene, length(valid_probes)))
    } else {
      cat(sprintf("  ⚠️ 警告: 基因 %s 未找到有效探针\n", gene))
    }
  }
} else {
  message("⚠️ 未找到完全匹配的注释列，启用备用正则匹配方案...")
  for (gene in target_genes) {
    matched_probes <- grep(gene, rownames(raw_matrix), value = TRUE)
    if(length(matched_probes) > 0) {
      val_gene_expr[gene, ] <- colMeans(raw_matrix[matched_probes, , drop = FALSE], na.rm = TRUE)
    }
  }
}

validation_set <- as.data.frame(t(val_gene_expr))

# ------------------------------------------------------------------------------
# 3. 提取临床分组并构造机器学习输入数据集 (ad_data)
# ------------------------------------------------------------------------------
cat("\n🔍 3. 正在提取并清洗 AD 与 Ctrl 临床分组...\n")

disease_column_name <- NULL
for(col in colnames(ext_clinical_info)) {
  unique_vals <- unique(as.character(ext_clinical_info[[col]]))
  if(any(grep("alzheimer|control|normal", unique_vals, ignore.case = TRUE))) {
    disease_column_name <- col
    break
  }
}

if(!is.null(disease_column_name)) {
  raw_labels <- as.character(ext_clinical_info[[disease_column_name]])
  is_ad   <- grepl("alzheimer|case", raw_labels, ignore.case = TRUE)
  is_ctrl <- grepl("control|normal|non-demented", raw_labels, ignore.case = TRUE)
  valid_indices <- which(is_ad | is_ctrl)
  
  cleaned_validation_set <- validation_set[valid_indices, , drop = FALSE]
  cleaned_validation_set$Status <- factor(ifelse(is_ad[valid_indices], "AD", "Ctrl"), levels = c("Ctrl", "AD"))
  ad_data <- cleaned_validation_set
  
  cat("🎉 数据清洗成功！总样本量:", nrow(ad_data), "\n")
  print(table(ad_data$Status))
} else {
  stop("❌ 严重错误：未能识别疾病分类标签列！")
}

# ------------------------------------------------------------------------------
# 4. 构建随机森林模型并评估基因调控方向与重要性
# ------------------------------------------------------------------------------
cat("\n🌲 4. 正在构建 5 基因随机森林 (Random Forest) 分类模型...\n")

set.seed(42) # 固定随机种子以保证结果完全复现
rf_ad <- randomForest(Status ~ ACP5 + RAMP3 + SLC7A14 + SLC2A13 + CLTC, 
                      data = ad_data, importance = TRUE, ntree = 1000)

# 计算各基因在 AD 组 vs Ctrl 组的表达均值，自动判别上调/下调方向
gene_means <- ad_data %>%
  group_by(Status) %>%
  summarise(across(all_of(target_genes), mean, na.rm = TRUE))

up_genes <- c()
down_genes <- c()
for (gene in target_genes) {
  mean_ad   <- gene_means[[gene]][gene_means$Status == "AD"]
  mean_ctrl <- gene_means[[gene]][gene_means$Status == "Ctrl"]
  if (mean_ad > mean_ctrl) {
    up_genes <- c(up_genes, gene)
  } else {
    down_genes <- c(down_genes, gene)
  }
}

# 提取基尼不纯度贡献度 (Mean Decrease Gini)
importance_df <- as.data.frame(importance(rf_ad))
importance_df$Gene <- rownames(importance_df)
importance_df <- importance_df %>%
  select(Gene, Importance = MeanDecreaseGini) %>%
  mutate(Regulation = if_else(Gene %in% up_genes, "Up-regulated (in AD)", "Down-regulated (in AD)")) %>%
  arrange(desc(Importance))

cat("\n📊 随机森林基因贡献度列表 (Importance Score)：\n")
print(importance_df)

# ------------------------------------------------------------------------------
# 5. 精修学术随机森林权重图渲染与 PDF 导出
# ------------------------------------------------------------------------------
cat("\n🎨 5. 正在导出精修随机森林权重矢量 PDF 图表...\n")

p_ad_importance <- ggplot(importance_df, aes(x = reorder(Gene, Importance), y = Importance, fill = Regulation)) +
  geom_bar(stat = "identity", width = 0.55, alpha = 0.9, color = "black", linewidth = 0.4) +
  coord_flip() +
  scale_fill_manual(values = c("Up-regulated (in AD)" = "#DC0000FF", "Down-regulated (in AD)" = "#3C5488FF")) +
  labs(
    title = "Random Forest Classifier Weights\n(Alzheimer's Disease: GSE33000)",
    x = "Shared Target Genes",
    y = "Mean Decrease Gini (Importance Score)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray95"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = ggplot2::margin(b = 10)),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(color = "black", size = 10.5),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 9.5),
    aspect.ratio = 0.75
  )

pdf_file <- file.path(base_dir, "Figure_AD_GSE33000_RF_Importance.pdf")
pdf(pdf_file, width = 6.5, height = 5.5, useDingbats = FALSE)
print(p_ad_importance)
dev.off()

cat("\n✨ 运行成功！随机森林特征权重精修矢量图已保存至本地：\n💾", pdf_file, "\n")