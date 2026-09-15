# ==============================================================================
# 🎯 [GSE205450 纯本地提取 + 5-Gene 联合诊断模型 + 动态标题 ROC 矢量图导出]
# ==============================================================================
# 📝 备注：运行本脚本前，请提前下载以下 GEO 数据文件并放置于脚本同级目录（./）：
# 
# 1. GSE205450 lcpm 表达矩阵文件 (GSE205450_lcpm.table.txt.gz):
#    下载链接: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE205450&format=file&file=GSE205450%5Flcpm%2Etable%2Etxt%2Egz
# ==============================================================================

library(dplyr)
library(tidyr)
library(pROC)

# ------------------------------------------------------------------------------
# 1. 检查并读取本地矩阵数据
# ------------------------------------------------------------------------------
base_dir <- "./"
local_file <- file.path(base_dir, "GSE205450_lcpm.table.txt.gz")

if (!file.exists(local_file)) {
  stop("⚠️ 错误：未在当前工作目录下找到 'GSE205450_lcpm.table.txt.gz'，请检查文件路径！")
}

cat("🛠️ 1. 正在读取本地矩阵表头，解析临床分组 (Control/PD & Caudate/Putamen)...\n")

con <- gzfile(local_file, "r")
header <- readLines(con, n = 1)
raw_sample_ids <- strsplit(header, "\t")[[1]][-1]
sample_ids <- gsub('"', '', raw_sample_ids)

# 从 Sample Title 中提取 Group 和 Center 信息
clinical_info <- data.frame(
  title = sample_ids,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Group = case_when(
      grepl("control|ctrl|normal|N_", title, ignore.case = TRUE) ~ "Control",
      grepl("PD|parkinson", title, ignore.case = TRUE) ~ "PD",
      TRUE ~ "PD"
    ),
    Center = case_when(
      grepl("caudate|cau", title, ignore.case = TRUE) ~ "Caudate",
      grepl("putamen|put", title, ignore.case = TRUE) ~ "Putamen",
      TRUE ~ "Caudate"
    )
  )

# 从本地 lcpm 矩阵中提取 5 个核心基因
target_order <- c("CLTC", "SLC7A14", "SLC2A13", "RAMP3", "ACP5")
extracted_rows <- list()

cat("🛠️ 2. 正在从本地矩阵检索 5 个核心基因表达量...\n")
line_count <- 0
pb <- txtProgressBar(min = 0, max = 30000, style = 3)

while (TRUE) {
  line <- readLines(con, n = 1)
  if (length(line) == 0) break
  
  line_count <- line_count + 1
  if (line_count %% 500 == 0) {
    setTxtProgressBar(pb, min(line_count, 30000))
  }
  
  parts <- strsplit(line, "\t")[[1]]
  gene_name <- gsub('"', '', parts[1])
  
  if (gene_name %in% target_order) {
    extracted_rows[[gene_name]] <- as.numeric(parts[-1])
  }
}

close(con)
setTxtProgressBar(pb, 30000)
close(pb)

# 合并表达矩阵与临床数据
expr_5_genes_df <- as.data.frame(extracted_rows)
expr_5_genes_df$Matrix_Title <- sample_ids

final_plot_data <- inner_join(
  clinical_info %>% dplyr::select(title, Group, Center),
  expr_5_genes_df,
  by = c("title" = "Matrix_Title")
)

# ------------------------------------------------------------------------------
# 2. 构建 Caudate & Putamen 脑区逻辑回归模型与 ROC 计算
# ------------------------------------------------------------------------------
cat("\n⚙️ 3. 正在按脑区构建 5-Gene 联合逻辑回归诊断模型...\n")

# --- Caudate 脑区 ---
caudate_data <- final_plot_data %>%
  filter(Center == "Caudate") %>%
  mutate(Status = ifelse(Group == "PD", 1, 0))

fit_caudate_5g <- glm(Status ~ CLTC + SLC7A14 + SLC2A13 + RAMP3 + ACP5, 
                      data = caudate_data, family = binomial(link = "logit"))
caudate_data$Prob_5g <- predict(fit_caudate_5g, type = "response")

roc_c_comb_5g <- roc(caudate_data$Status, caudate_data$Prob_5g, quiet = TRUE)
roc_c_cltc    <- roc(caudate_data$Status, caudate_data$CLTC, quiet = TRUE)
roc_c_slc7a14 <- roc(caudate_data$Status, caudate_data$SLC7A14, quiet = TRUE)
roc_c_slc2a13 <- roc(caudate_data$Status, caudate_data$SLC2A13, quiet = TRUE)
roc_c_ramp3   <- roc(caudate_data$Status, caudate_data$RAMP3, quiet = TRUE)
roc_c_acp5    <- roc(caudate_data$Status, caudate_data$ACP5, quiet = TRUE)

# --- Putamen 脑区 ---
putamen_data <- final_plot_data %>%
  filter(Center == "Putamen") %>%
  mutate(Status = ifelse(Group == "PD", 1, 0))

fit_putamen_5g <- glm(Status ~ CLTC + SLC7A14 + SLC2A13 + RAMP3 + ACP5, 
                      data = putamen_data, family = binomial(link = "logit"))
putamen_data$Prob_5g <- predict(fit_putamen_5g, type = "response")

roc_p_comb_5g <- roc(putamen_data$Status, putamen_data$Prob_5g, quiet = TRUE)
roc_p_cltc    <- roc(putamen_data$Status, putamen_data$CLTC, quiet = TRUE)
roc_p_slc7a14 <- roc(putamen_data$Status, putamen_data$SLC7A14, quiet = TRUE)
roc_p_slc2a13 <- roc(putamen_data$Status, putamen_data$SLC2A13, quiet = TRUE)
roc_p_ramp3   <- roc(putamen_data$Status, putamen_data$RAMP3, quiet = TRUE)
roc_p_acp5    <- roc(putamen_data$Status, putamen_data$ACP5, quiet = TRUE)

# 动态统计样本数
n_c_ctrl <- sum(caudate_data$Status == 0)
n_c_pd   <- sum(caudate_data$Status == 1)
n_p_ctrl <- sum(putamen_data$Status == 0)
n_p_pd   <- sum(putamen_data$Status == 1)

# 控制台汇总打印 AUC 结果
cat("\n==================================================\n")
cat("📊 【Caudate (尾状核) 脑区诊断效能评价 (AUC)】：\n")
cat(sprintf("👑 5-Gene Combined Model AUC : %6.4f\n", roc_c_comb_5g$auc))
cat(sprintf("🧬 CLTC AUC                  : %6.4f\n", roc_c_cltc$auc))
cat(sprintf("🧬 SLC7A14 AUC               : %6.4f\n", roc_c_slc7a14$auc))
cat(sprintf("🧬 SLC2A13 AUC               : %6.4f\n", roc_c_slc2a13$auc))
cat(sprintf("🧬 RAMP3 AUC                 : %6.4f\n", roc_c_ramp3$auc))
cat(sprintf("🧬 ACP5 AUC                  : %6.4f\n", roc_c_acp5$auc))

cat("\n📊 【Putamen (壳核) 脑区诊断效能评价 (AUC)】：\n")
cat(sprintf("👑 5-Gene Combined Model AUC : %6.4f\n", roc_p_comb_5g$auc))
cat(sprintf("🧬 CLTC AUC                  : %6.4f\n", roc_p_cltc$auc))
cat(sprintf("🧬 SLC7A14 AUC               : %6.4f\n", roc_p_slc7a14$auc))
cat(sprintf("🧬 SLC2A13 AUC               : %6.4f\n", roc_p_slc2a13$auc))
cat(sprintf("🧬 RAMP3 AUC                 : %6.4f\n", roc_p_ramp3$auc))
cat(sprintf("🧬 ACP5 AUC                  : %6.4f\n", roc_p_acp5$auc))
cat("==================================================\n\n")

# ------------------------------------------------------------------------------
# 3. 绘制全景双脑区 5-Gene 对齐版 ROC 并保存 PDF
# ------------------------------------------------------------------------------
pdf_output_name <- file.path(base_dir, "PD_5_Gene_Core_Lysosome_Diagnostic_ROC.pdf")

pdf(pdf_output_name, width = 11, height = 5.5, family = "Helvetica", useDingbats = FALSE)

# 设置 1 行 2 列布局与学术边距
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1.5), mgp = c(2.5, 0.7, 0), cex.main = 1.05, cex.lab = 1)

# 学术专属 Palette：联合模型为醒目深红实线，单基因采用柔和虚线
palette_colors <- c("#DC0000FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF")

# ------------------------------------------------------------------------------
# 子图 1：Panel B1 (Caudate Center)
# ------------------------------------------------------------------------------
title_c <- sprintf("Panel B1: Clinical Diagnostic Efficacy\n(Caudate Center, Ctrl:%d vs PD:%d)", n_c_ctrl, n_c_pd)

plot(roc_c_comb_5g, col = palette_colors[1], lwd = 3.5, legacy.axes = TRUE,
     xlab = "False Positive Rate (1 - Specificity)",
     ylab = "True Positive Rate (Sensitivity)",
     main = title_c)

plot(roc_c_cltc,    col = palette_colors[2], lwd = 1.5, lty = 2, add = TRUE, legacy.axes = TRUE)
plot(roc_c_slc7a14, col = palette_colors[3], lwd = 1.5, lty = 2, add = TRUE, legacy.axes = TRUE)
plot(roc_c_slc2a13, col = palette_colors[4], lwd = 1.5, lty = 2, add = TRUE, legacy.axes = TRUE)
plot(roc_c_ramp3,   col = palette_colors[5], lwd = 1.5, lty = 2, add = TRUE, legacy.axes = TRUE)
plot(roc_c_acp5,    col = palette_colors[6], lwd = 1.5, lty = 2, add = TRUE, legacy.axes = TRUE)
abline(a = 0, b = 1, lty = 3, col = "gray60")

legend("bottomright",
       legend = c(
         sprintf("5-Gene Combined Model (AUC = %.4f)", roc_c_comb_5g$auc),
         sprintf("CLTC (AUC = %.4f)", roc_c_cltc$auc),
         sprintf("SLC7A14 (AUC = %.4f)", roc_c_slc7a14$auc),
         sprintf("SLC2A13 (AUC = %.4f)", roc_c_slc2a13$auc),
         sprintf("RAMP3 (AUC = %.4f)", roc_c_ramp3$auc),
         sprintf("ACP5 (AUC = %.4f)", roc_c_acp5$auc)
       ),
       col = palette_colors, lty = c(1, 2, 2, 2, 2, 2), lwd = c(3.5, 1.5, 1.5, 1.5, 1.5, 1.5), bty = "n", cex = 0.8)

# ------------------------------------------------------------------------------
# 子图 2：Panel B2 (Putamen Center)
# ------------------------------------------------------------------------------
title_p <- sprintf("Panel B2: Clinical Diagnostic Efficacy\n(Putamen Center, Ctrl:%d vs PD:%d)", n_p_ctrl, n_p_pd)

plot(roc_p_comb_5g, col = palette_colors[1], lwd = 3.5, legacy.axes = TRUE,
     xlab = "False Positive Rate (1 - Specificity)",
     ylab = "True Positive Rate (Sensitivity)",
     main = title_p)

plot(roc_p_cltc,    col = palette_colors[2], lwd = 1.5, lty = 2, add = TRUE, legacy.axes = TRUE)
plot(roc_p_slc7a14, col = palette_colors[3], lwd = 1.5, lty = 2, add = TRUE, legacy.axes = TRUE)
plot(roc_p_slc2a13, col = palette_colors[4], lwd = 1.5, lty = 2, add = TRUE, legacy.axes = TRUE)
plot(roc_p_ramp3,   col = palette_colors[5], lwd = 1.5, lty = 2, add = TRUE, legacy.axes = TRUE)
plot(roc_p_acp5,    col = palette_colors[6], lwd = 1.5, lty = 2, add = TRUE, legacy.axes = TRUE)
abline(a = 0, b = 1, lty = 3, col = "gray60")

legend("bottomright",
       legend = c(
         sprintf("5-Gene Combined Model (AUC = %.4f)", roc_p_comb_5g$auc),
         sprintf("CLTC (AUC = %.4f)", roc_p_cltc$auc),
         sprintf("SLC7A14 (AUC = %.4f)", roc_p_slc7a14$auc),
         sprintf("SLC2A13 (AUC = %.4f)", roc_p_slc2a13$auc),
         sprintf("RAMP3 (AUC = %.4f)", roc_p_ramp3$auc),
         sprintf("ACP5 (AUC = %.4f)", roc_p_acp5$auc)
       ),
       col = palette_colors, lty = c(1, 2, 2, 2, 2, 2), lwd = c(3.5, 1.5, 1.5, 1.5, 1.5, 1.5), bty = "n", cex = 0.8)

dev.off()
par(mfrow = c(1, 1))

cat(paste0("\n✨ 成功导出高质量矢量 ROC 全景图：'", pdf_output_name, "'\n"))