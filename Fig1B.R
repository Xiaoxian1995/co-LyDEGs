# ==============================================================================
# 项目名称：GSE68719 帕金森病 (PD) RNA-Seq 差异表达分析与火山图绘制
# 说明：本脚本包含表型清洗、表达矩阵对齐、Limma差异分析、Gene ID映射与可视化。
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. 本地数据准备说明 (Data Preparation Guidelines)
# ------------------------------------------------------------------------------
# 请先从 GEO 数据库 (GSE68719) 下载以下两个数据文件，并在当前工作目录下
# 创建名为 `data` 的文件夹，将文件放入其中：
#
# 1. 临床表型矩阵文件：
#    下载链接: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68719
#    放置路径: ./data/GSE68719_series_matrix.txt.gz
#
# 2. DESeq2 标准化表达计数文件：
#    下载链接: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68719 (补充文件 Supp File)
#    放置路径: ./data/GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz
# ------------------------------------------------------------------------------

# 1. 加载所需的依赖包
library(GEOquery)
library(limma)
library(org.Hs.eg.db)
library(ggplot2)
library(tidyverse)

# ------------------------------------------------------------------------------
# 2. 配置相对路径与输入文件 (Relative Paths & Input Files)
# ------------------------------------------------------------------------------
# 自动检测并创建 data 文件夹（如不存在）
if (!dir.exists("data")) {
  dir.create("data")
}

# 使用相对路径定义输入/输出目录
data_dir <- "./data"

pd_series_file  <- file.path(data_dir, "GSE68719_series_matrix.txt.gz")
norm_counts_file <- file.path(data_dir, "GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz")

# 检查输入文件是否存在
if (!file.exists(pd_series_file) || !file.exists(norm_counts_file)) {
  stop("❌ 未在 ./data/ 目录下找到所需文件！请参照文件开头的【数据准备说明】放置文件。")
}

# ------------------------------------------------------------------------------
# 3. 读取本地 PD series_matrix 并清洗提取表型分组信息
# ------------------------------------------------------------------------------
message("正在读取本地 PD 样本及临床矩阵信息...")
gse_pd <- getGEO(filename = pd_series_file, GSEMatrix = TRUE, getGPL = FALSE)
sample_info_pd <- pData(gse_pd)

cat("成功读取表型矩阵！包含：\n")
cat("👉 原始样本数 (Rows):", nrow(sample_info_pd), "\n")
cat("👉 临床特征指标数 (Columns):", ncol(sample_info_pd), "\n\n")

# 根据 title 列清洗提取分组信息（包含 'P' 为 PD，包含 'C' 为 Control）
sample_info_pd_clean <- sample_info_pd %>%
  mutate(
    Group = case_when(
      grepl("P", title, ignore.case = FALSE) ~ "PD",
      grepl("C", title, ignore.case = FALSE) ~ "Control",
      TRUE ~ NA_character_
    ),
    # 提取纯净样本 ID（切除 "[...]" 及多余空格，提取如 C_0002）
    Clean_ID = str_trim(sub("\\[.*\\]", "", title))
  )

# 剔除未成功分类的样本
valid_pd_samples <- sample_info_pd_clean %>% 
  filter(!is.na(Group))

# 将行名重置为 Clean_ID 以便后续与表达矩阵列名精准对齐
rownames(valid_pd_samples) <- valid_pd_samples$Clean_ID

cat("==================================================\n")
cat("📊 【PD 数据集分组清洗完成】\n")
cat("--------------------------------------------------\n")
cat("🔴 帕金森病组 (PD) 样本数:      ", sum(valid_pd_samples$Group == "PD"), " 个\n")
cat("🔵 正常对照组 (Control) 样本数:  ", sum(valid_pd_samples$Group == "Control"), " 个\n")
cat("==================================================\n\n")

# ------------------------------------------------------------------------------
# 4. 读取 DESeq2 标准化表达矩阵与样本严格对齐
# ------------------------------------------------------------------------------
message("正在读取本地 DESeq2 标准化表达矩阵...")
pd_expr_df <- read.table(gzfile(norm_counts_file), header = TRUE, sep = "\t", check.names = FALSE)

# 将第一列（Ensembl ID）设置为行名
rownames(pd_expr_df) <- pd_expr_df[, 1]
expr_matrix_pd_real <- as.matrix(pd_expr_df[, -1])

# 样本 ID 严格取交集对齐
common_pd_samples <- intersect(colnames(expr_matrix_pd_real), rownames(valid_pd_samples))

if (length(common_pd_samples) == 0) {
  stop("❌ 样本 ID 无法匹配，请检查表达矩阵列名与 title 提取的 Clean_ID 是否一致！")
}

expr_matrix_pd_final <- expr_matrix_pd_real[, common_pd_samples]
sample_info_pd_final <- valid_pd_samples[common_pd_samples, ]

cat("==================================================\n")
cat("🎉 【PD 数据集：基因与样本双重对齐成功！】\n")
cat("--------------------------------------------------\n")
cat("👉 对齐后的基因数 (行):   ", nrow(expr_matrix_pd_final), " 个\n")
cat("👉 对齐后的样本数 (列):   ", ncol(expr_matrix_pd_final), " 个\n")
cat("==================================================\n\n")

# 使用相对路径导出对齐好的文件至 data 目录
write.csv(expr_matrix_pd_final, file.path(data_dir, "GSE68719_expr_matrix_aligned.csv"))
write.csv(sample_info_pd_final, file.path(data_dir, "GSE68719_sample_info_aligned.csv"))

# ------------------------------------------------------------------------------
# 5. 数据转换与 Limma 差异表达分析
# ------------------------------------------------------------------------------
# 强制转换矩阵元素为数值型
expr_numeric <- apply(expr_matrix_pd_final, 2, as.numeric)
rownames(expr_numeric) <- rownames(expr_matrix_pd_final)

# 自动检查并进行 log2(x + 1) 转换
if (max(expr_numeric, na.rm = TRUE) > 50) {
  message("检测到数据为非 log 状态，正在自动进行 log2(x + 1) 转换...")
  expr_matrix_pd_log <- log2(expr_numeric + 1)
} else {
  expr_matrix_pd_log <- expr_numeric
}

# 构建实验设计矩阵 (以 Control 组作为基准 reference)
group_factor <- factor(sample_info_pd_final$Group, levels = c("Control", "PD"))
design <- model.matrix(~ group_factor)
colnames(design) <- c("Intercept", "PD_vs_Control")

# 运行 limma 线性模型
message("正在运行 limma 差异表达分析 (PD vs Control)...")
fit <- lmFit(expr_matrix_pd_log, design)
fit <- eBayes(fit)

# 提取完整差异结果
deg_pd_all <- topTable(fit, coef = "PD_vs_Control", number = Inf, adjust.method = "BH") %>%
  rownames_to_column("Ensembl_ID")

# ------------------------------------------------------------------------------
# 6. 基因 ID 注释转换 (Ensembl ID -> Gene Symbol) 与 0.5 阈值数据绑定
# ------------------------------------------------------------------------------
# 去掉 Ensembl ID 的版本号后缀 (如 ENSG00000137731.8 -> ENSG00000137731)
deg_pd_all$Ensembl_clean <- sub("\\..*$", "", deg_pd_all$Ensembl_ID)

message("正在将 Ensembl ID 映射转换为 Gene Symbol...")
gene_ids <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = deg_pd_all$Ensembl_clean,
  columns = "SYMBOL",
  keytype = "ENSEMBL"
)

# 合并转换结果，过滤未匹配项与重复基因名
deg_pd_mapped <- left_join(deg_pd_all, gene_ids, by = c("Ensembl_clean" = "ENSEMBL")) %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  distinct(SYMBOL, .keep_all = TRUE)

# 保存注释后的完整差异分析结果到 ./data/ 相对路径
write.csv(deg_pd_mapped, file.path(data_dir, "GSE68719_DEG_limma_results.csv"), row.names = FALSE)

# 明确绑定 PD 的核心差异表并按 0.5 阈值分类
pd_raw_data <- as.data.frame(deg_pd_mapped)

# 自动适配列名：不管叫 P.Value 还是 adj.P.Val
p_col_pd <- if("P.Value" %in% colnames(pd_raw_data)) "P.Value" else "adj.P.Val"

# 数据清洗与分类标签（阈值严格同步为 0.5）
pd_plot_data <- pd_raw_data %>%
  mutate(
    Group = case_when(
      logFC > 0.5 & .data[[p_col_pd]] < 0.05 ~ "Up",
      logFC < -0.5 & .data[[p_col_pd]] < 0.05 ~ "Down",
      TRUE ~ "None"
    ),
    LogP = -log10(.data[[p_col_pd]])
  ) %>%
  filter(is.finite(logFC) & is.finite(LogP))

# 📊 【实时统计】输出 PD 在 0.5 阈值下的准确上下调基因数量
cat("\n====================================================================\n")
cat("📊 【PD 差异基因精准计数：阈值 |logFC| > 0.5 & P < 0.05】\n")
cat("--------------------------------------------------------------------\n")
pd_counts_table <- table(pd_plot_data$Group)
print(pd_counts_table)
cat("====================================================================\n\n")

# ------------------------------------------------------------------------------
# 7. 学术级火山图绘制并在 Plots 面板刷新弹图
# ------------------------------------------------------------------------------
volcano_colors <- c("Down" = "#2B4C99", "None" = "#BBBBBB", "Up" = "#E62129")

pd_volcano_plot <- ggplot(pd_plot_data, aes(x = logFC, y = LogP)) +
  geom_point(aes(color = Group), size = 1.3, alpha = 0.6) +
  scale_color_manual(values = volcano_colors) +
  
  # 细黑色边界虚线（设定为 0.5）
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black", size = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.4) +
  
  # 横纵刻度轴展示设置
  scale_x_continuous(limits = c(-2.5, 2.5), breaks = seq(-2, 2, by = 1)) +
  scale_y_continuous(limits = c(0, 10.5), breaks = seq(0, 10, by = 2.5), labels = c("0", "2.5", "5.0", "7.5", "10.0")) +
  
  labs(
    title = "Volcano Plot: PD vs Control",
    x = "Log2 fold change",
    y = "-Log10 P"
  ) +
  
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.6),
    
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "black"),
    
    text = element_text(family = "sans", color = "black"),
    axis.title = element_text(size = 14, face = "plain"),
    axis.text = element_text(size = 12, color = "black"),
    
    axis.ticks.length = unit(-0.15, "cm"),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.text.x = element_text(margin = margin(t = 8)),
    axis.text.y = element_text(margin = margin(r = 8)),
    
    legend.position = "right",
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 12, face = "bold")
  )

# 刷新并输出火山图预览
while (!is.null(dev.list())) { dev.off() }
print(pd_volcano_plot)

cat(">> GSE68719 分析全部完成！中间结果已自动保存至 ./data/ 文件夹。\n")