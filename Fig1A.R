# ==============================================================================
# 项目名称：GSE68719 帕金森病 (PD) RNA-Seq 数据预处理与样本-基因矩阵对齐
# 脚本功能：读取本地表型与 DESeq2 标准化表达矩阵，清洗 ID 并完成双重对齐
# ==============================================================================

# ==============================================================================
# 【本地输入文件说明】
# 请在当前工作目录下创建 `data/` 文件夹，并将以下 3 个文件放入 `./data/` 目录中：
# 1. 表型矩阵文件：GSE68719_series_matrix.txt.gz
# 2. 标准化 count 矩阵文件：GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz
# 3. 官方差异分析结果文件：GSE68719_mlpd_DESeq2_diffexp_pmi_age_rin_default.txt.gz
# ==============================================================================

# 1. 加载所需的依赖包
library(GEOquery)
library(tidyverse)

# 自动创建输出文件夹（如果不存在）
if (!dir.exists("./results")) {
  dir.create("./results")
}

data_dir <- "./data"

# ------------------------------------------------------------------------------
# 2. 读取本地 PD series_matrix 提取临床表型信息 (Phenotype Data)
# ------------------------------------------------------------------------------
pd_series_file <- file.path(data_dir, "GSE68719_series_matrix.txt.gz")

if (!file.exists(pd_series_file)) {
  stop("❌ 未在 ./data/ 目录下找到 GSE68719_series_matrix.txt.gz，请检查文件路径！")
}

message("正在读取本地 PD 样本及临床矩阵信息...")
gse_pd <- getGEO(filename = pd_series_file, GSEMatrix = TRUE, getGPL = FALSE)
sample_info_pd <- pData(gse_pd)

cat("成功读取表型矩阵！包含：\n")
cat("👉 原始样本数 (Rows):", nrow(sample_info_pd), "\n")
cat("👉 临床特征指标数 (Columns):", ncol(sample_info_pd), "\n\n")

# ------------------------------------------------------------------------------
# 3. 表型数据清洗与分组提取 (PD vs Control)
# ------------------------------------------------------------------------------
# 根据 title 列清洗提取分组信息（包含 'P' 为 PD，包含 'C' 为 Control）
sample_info_pd_clean <- sample_info_pd %>%
  mutate(
    Group = case_when(
      grepl("P", title, ignore.case = FALSE) ~ "PD",
      grepl("C", title, ignore.case = FALSE) ~ "Control",
      TRUE ~ NA_character_
    ),
    # 提取纯净样本 ID（例如切除 "[...]" 及多余空格，提取如 C_0002）
    Clean_ID = str_trim(sub("\\[.*\\]", "", title))
  )

# 剔除未成功分类的异常样本
valid_pd_samples <- sample_info_pd_clean %>% 
  filter(!is.na(Group))

# 将行名重置为 Clean_ID 以便后续与表达矩阵的列名精确对齐
rownames(valid_pd_samples) <- valid_pd_samples$Clean_ID

cat("==================================================\n")
cat("📊 【PD 数据集分组清洗完成】\n")
cat("--------------------------------------------------\n")
cat("🔴 帕金森病组 (PD) 样本数:      ", sum(valid_pd_samples$Group == "PD"), " 个\n")
cat("🔵 正常对照组 (Control) 样本数:  ", sum(valid_pd_samples$Group == "Control"), " 个\n")
cat("==================================================\n\n")

# ------------------------------------------------------------------------------
# 4. 读取本地 DESeq2 标准化表达矩阵文件
# ------------------------------------------------------------------------------
norm_counts_file <- file.path(data_dir, "GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz")

if (!file.exists(norm_counts_file)) {
  stop("❌ 未在 ./data/ 目录下找到 GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz！")
}

message("正在读取本地 DESeq2 标准化表达矩阵...")
pd_expr_df <- read.table(gzfile(norm_counts_file), header = TRUE, sep = "\t", check.names = FALSE)

# 将第一列（Gene Symbol / Ensembl ID）设置为行名
rownames(pd_expr_df) <- pd_expr_df[, 1]
expr_matrix_pd_real <- as.matrix(pd_expr_df[, -1])

# ------------------------------------------------------------------------------
# 5. 表达矩阵与样本临床信息交集严格对齐
# ------------------------------------------------------------------------------
common_pd_samples <- intersect(colnames(expr_matrix_pd_real), rownames(valid_pd_samples))

if (length(common_pd_samples) == 0) {
  stop("❌ 样本 ID 无法匹配，请检查表达矩阵列名与 title 提取的 Clean_ID 是否一致！")
}

expr_matrix_pd_final <- expr_matrix_pd_real[, common_pd_samples]
sample_info_pd_final <- valid_pd_samples[common_pd_samples, ]

cat("\n==================================================\n")
cat("🎉 【PD 数据集：基因与样本双重对齐终极大功告成！】\n")
cat("--------------------------------------------------\n")
cat("👉 最终对齐后的基因数 (行):   ", nrow(expr_matrix_pd_final), " 个\n")
cat("👉 最终对齐后的样本数 (列):   ", ncol(expr_matrix_pd_final), " 个\n")
cat("🔴 其中 帕金森病 (PD) 样本:    ", sum(sample_info_pd_final$Group == "PD"), " 个\n")
cat("🔵 其中 正常对照 (Control) 样本:", sum(sample_info_pd_final$Group == "Control"), " 个\n")
cat("==================================================\n\n")

print("📋 基因行名前 5 个预览：")
print(head(rownames(expr_matrix_pd_final), 5))

# ------------------------------------------------------------------------------
# 6. 保存清洗并对齐后的数据到 ./results/ 目录
# ------------------------------------------------------------------------------
write.csv(expr_matrix_pd_final, "./results/GSE68719_expr_matrix_aligned.csv")
write.csv(sample_info_pd_final, "./results/GSE68719_sample_info_aligned.csv")

cat("\n>> 对齐后的矩阵与样本信息已分别保存至 ./results/ 目录下。\n")