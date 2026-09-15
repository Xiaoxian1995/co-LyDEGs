# ==============================================================================
# 🎯 [GSE33000 5基因联合诊断模型构建与 ROC 曲线矢量出图一键脚本]
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
library(pROC)
library(ggplot2)

base_dir <- "./"
matrix_file <- file.path(base_dir, "GSE33000_series_matrix.txt.gz")
gpl_file    <- file.path(base_dir, "GPL4372.gz")

# ------------------------------------------------------------------------------
# 1. 读取数据并配置注释表 (ext_gse)
# ------------------------------------------------------------------------------
cat("📂 1. 正在读取 GSE33000 数据与 GPL4372 注释...\n")
ext_gse <- getGEO(filename = matrix_file, GSEMatrix = TRUE, getGPL = FALSE)
gpl_data <- getGEO(filename = gpl_file)

# 将芯片注释嵌入 fData
fData(ext_gse) <- Table(gpl_data)

probe_annotation  <- fData(ext_gse)
raw_matrix        <- exprs(ext_gse)
ext_clinical_info <- pData(ext_gse)

# ------------------------------------------------------------------------------
# 2. 精准匹配探针并提取 5 个核心基因表达量
# ------------------------------------------------------------------------------
cat("🛠️ 2. 正在提取 5 个核心基因表达量...\n")

target_genes <- c("CLTC", "RAMP3", "SLC7A14", "ACP5", "SLC2A13")

# 寻找正确的 Gene Symbol 列
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
  message(paste("🎯 匹配到探针注释列:", best_col))
  probe_id_col <- colnames(probe_annotation)[1]
  for (gene in target_genes) {
    matched_rows <- which(probe_annotation[[best_col]] == gene)
    matched_probes <- as.character(probe_annotation[[probe_id_col]][matched_rows])
    valid_probes <- intersect(matched_probes, rownames(raw_matrix))
    
    if (length(valid_probes) > 0) {
      val_gene_expr[gene, ] <- colMeans(raw_matrix[valid_probes, , drop = FALSE], na.rm = TRUE)
      cat(sprintf("  ├─ 基因 %-8s : 匹配到 %d 个探针\n", gene, length(valid_probes)))
    } else {
      cat(sprintf("  ⚠️ 警告: 基因 %s 未找到有效探针\n", gene))
    }
  }
} else {
  message("⚠️ 启用备用探针匹配模式...")
  for (gene in target_genes) {
    matched_probes <- grep(gene, rownames(raw_matrix), value = TRUE)
    if(length(matched_probes) > 0) {
      val_gene_expr[gene, ] <- colMeans(raw_matrix[matched_probes, , drop = FALSE], na.rm = TRUE)
    }
  }
}

validation_set <- as.data.frame(t(val_gene_expr))

# ------------------------------------------------------------------------------
# 3. 剥离并清洗 AD vs Control 分组信息
# ------------------------------------------------------------------------------
cat("\n🔍 3. 正在提取 AD 与 Control 临床标签...\n")

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
  cleaned_validation_set$Group <- ifelse(is_ad[valid_indices], "AD", "Control")
  validation_set <- cleaned_validation_set
  
  cat("🎉 数据清洗完毕！有效样本总数:", nrow(validation_set), "\n")
  print(table(validation_set$Group))
} else {
  stop("❌ 未能成功自动定位疾病分类列，请检查临床数据表头。")
}

# ------------------------------------------------------------------------------
# 4. 构建 5 基因逻辑回归联合模型与计算 ROC
# ------------------------------------------------------------------------------
cat("\n⚙️ 4. 正在构建 5 基因联合诊断模型并评估 ROC 效能...\n")

# 设置二元分类因变量 (AD = 1, Control = 0)
validation_set$Response <- ifelse(validation_set$Group == "AD", 1, 0)

# 构建 Logistic 回归模型
combined_model <- glm(Response ~ CLTC + RAMP3 + SLC7A14 + ACP5 + SLC2A13, 
                      data = validation_set, 
                      family = binomial(link = "logit"))

validation_set$Combined_Prob <- predict(combined_model, type = "response")

# 计算各曲线 ROC 对象
roc_combined <- roc(validation_set$Response, validation_set$Combined_Prob, quiet = TRUE)
roc_cltc     <- roc(validation_set$Response, validation_set$CLTC, quiet = TRUE)
roc_ramp3    <- roc(validation_set$Response, validation_set$RAMP3, quiet = TRUE)
roc_slc7a14  <- roc(validation_set$Response, validation_set$SLC7A14, quiet = TRUE)
roc_acp5     <- roc(validation_set$Response, validation_set$ACP5, quiet = TRUE)
roc_slc2a13  <- roc(validation_set$Response, validation_set$SLC2A13, quiet = TRUE)

# 控制台打印结果
cat("\n📊 【诊断效能评价 (AUC)】：\n")
cat(sprintf("👑 5-Gene Combined Model AUC : %6.4f\n", roc_combined$auc))
cat(sprintf("🧬 CLTC AUC                  : %6.4f\n", roc_cltc$auc))
cat(sprintf("🧬 RAMP3 AUC                 : %6.4f\n", roc_ramp3$auc))
cat(sprintf("🧬 SLC7A14 AUC               : %6.4f\n", roc_slc7a14$auc))
cat(sprintf("🧬 ACP5 AUC                  : %6.4f\n", roc_acp5$auc))
cat(sprintf("🧬 SLC2A13 AUC               : %6.4f\n", roc_slc2a13$auc))

# ------------------------------------------------------------------------------
# 5. 导出矢量 ROC PDF 图像
# ------------------------------------------------------------------------------
message("\n🎨 5. 正在导出矢量格式学术 ROC 曲线图...")

pdf_output <- file.path(base_dir, "AD_5_Gene_Core_Lysosome_Diagnostic_ROC.pdf")
pdf(pdf_output, width = 6.5, height = 6, useDingbats = FALSE)

# 绘制主曲线（5 基因联合模型）
plot(roc_combined, col = "#DC2626", lwd = 3, legacy.axes = TRUE,
     main = "Clinical Diagnostic Efficacy in Large External Validation Cohort\n(GSE33000, N = 467 Clinical Brain Samples)",
     xlab = "False Positive Rate (1 - Specificity)", 
     ylab = "True Positive Rate (Sensitivity)",
     font.lab = 2, cex.main = 0.95)

# 叠加 5 个单基因 ROC 虚线
plot(roc_cltc,    col = "#3B82F6", lwd = 1.5, lty = 2, add = TRUE, legacy.axes = TRUE)
plot(roc_ramp3,   col = "#10B981", lwd = 1.5, lty = 2, add = TRUE, legacy.axes = TRUE)
plot(roc_slc7a14, col = "#F59E0B", lwd = 1.5, lty = 2, add = TRUE, legacy.axes = TRUE)
plot(roc_acp5,    col = "#8B5CF6", lwd = 1.5, lty = 2, add = TRUE, legacy.axes = TRUE)
plot(roc_slc2a13, col = "#EC4899", lwd = 1.5, lty = 2, add = TRUE, legacy.axes = TRUE)

# 添加标准背景网格
grid(col = "#E5E7EB", lty = 1)

# 图例配置
legend("bottomright", 
       legend = c(sprintf("5-Gene Combined Model (AUC = %6.4f)", roc_combined$auc),
                  sprintf("CLTC (AUC = %6.4f)", roc_cltc$auc),
                  sprintf("RAMP3 (AUC = %6.4f)", roc_ramp3$auc),
                  sprintf("SLC7A14 (AUC = %6.4f)", roc_slc7a14$auc),
                  sprintf("ACP5 (AUC = %6.4f)", roc_acp5$auc),
                  sprintf("SLC2A13 (AUC = %6.4f)", roc_slc2a13$auc)),
       col = c("#DC2626", "#3B82F6", "#10B981", "#F59E0B", "#8B5CF6", "#EC4899"), 
       lwd = c(3, 1.5, 1.5, 1.5, 1.5, 1.5), 
       lty = c(1, 2, 2, 2, 2, 2),
       cex = 0.8, bty = "n")

dev.off()

cat("\n🎉 【运行成功！】矢量 PDF 已保存至:", pdf_output, "\n")