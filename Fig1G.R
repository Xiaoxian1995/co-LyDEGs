# ==============================================================================
# 项目名称: GSE68719 帕金森病 (PD) 数据分析 & 目标 5 基因学术级散点箱线图绘制
# 数据说明:
# 请先从 NCBI GEO 数据库下载以下两个文件，并放置在当前工作目录的 data/ 文件夹下：
#   - data/GSE68719_series_matrix.txt.gz
#   - data/GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz
# ==============================================================================

# 1. 加载所需的依赖包
library(GEOquery)
library(limma)
library(org.Hs.eg.db)
library(ggplot2)
library(tidyverse)
library(ggpubr)

# 自动创建规范的相对路径目录
if (!dir.exists("data")) dir.create("data")
if (!dir.exists("results")) dir.create("results")

# ------------------------------------------------------------------------------
# 2. 读取 GEO 数据与提取表型信息
# ------------------------------------------------------------------------------
pd_series_file <- file.path("data", "GSE68719_series_matrix.txt.gz")
norm_counts_file <- file.path("data", "GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz")

if (!file.exists(pd_series_file)) {
  stop("❌ 请先从 GEO 下载 GSE68719 数据文件并放置在当前工作目录的 data/ 文件夹下！")
}
if (!file.exists(norm_counts_file)) {
  stop("❌ 请先从 GEO 下载 GSE68719 Counts 文件并放置在当前工作目录的 data/ 文件夹下！")
}

message("正在读取 PD 样本及临床矩阵信息...")
gse_pd <- getGEO(filename = pd_series_file, GSEMatrix = TRUE, getGPL = FALSE)
sample_info_pd <- pData(gse_pd)

# 根据 title 清洗提取分组信息（包含 'P' 为 PD，包含 'C' 为 Control）
sample_info_pd_clean <- sample_info_pd %>%
  mutate(
    Group = case_when(
      grepl("P", title, ignore.case = FALSE) ~ "PD",
      grepl("C", title, ignore.case = FALSE) ~ "Control",
      TRUE ~ NA_character_
    ),
    Clean_ID = str_trim(sub("\\[.*\\]", "", title))
  )

valid_pd_samples <- sample_info_pd_clean %>% filter(!is.na(Group))
rownames(valid_pd_samples) <- valid_pd_samples$Clean_ID

# ------------------------------------------------------------------------------
# 3. 读取 DESeq2 标准化表达矩阵与样本对齐
# ------------------------------------------------------------------------------
message("正在读取 DESeq2 标准化表达矩阵...")
pd_expr_df <- read.table(gzfile(norm_counts_file), header = TRUE, sep = "\t", check.names = FALSE)

# 将第一列（Ensembl ID）设置为行名
rownames(pd_expr_df) <- pd_expr_df[, 1]
expr_matrix_pd_real <- as.matrix(pd_expr_df[, -1])

# 样本 ID 匹配对齐
common_pd_samples <- intersect(colnames(expr_matrix_pd_real), rownames(valid_pd_samples))

if (length(common_pd_samples) == 0) {
  stop("❌ 样本 ID 无法匹配，请检查表达矩阵列名与 Clean_ID 是否一致！")
}

expr_matrix_pd_final <- expr_matrix_pd_real[, common_pd_samples]
sample_info_pd_final <- valid_pd_samples[common_pd_samples, ]

# 转换数值矩阵并做 log2(x + 1) 标准化处理
expr_numeric <- apply(expr_matrix_pd_final, 2, as.numeric)
rownames(expr_numeric) <- rownames(expr_matrix_pd_final)

if (max(expr_numeric, na.rm = TRUE) > 50) {
  expr_matrix_pd_log <- log2(expr_numeric + 1)
} else {
  expr_matrix_pd_log <- expr_numeric
}

# ------------------------------------------------------------------------------
# 4. Ensembl ID 转换为 Gene Symbol 并构建基因表达矩阵
# ------------------------------------------------------------------------------
message("正在将 Ensembl ID 映射转换为 Gene Symbol...")
ensembl_clean <- sub("\\..*$", "", rownames(expr_matrix_pd_log))

gene_ids <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(ensembl_clean),
  columns = "SYMBOL",
  keytype = "ENSEMBL"
)

# 构建包含 Gene Symbol 的表达数据框
expr_mapped_df <- as.data.frame(expr_matrix_pd_log) %>%
  mutate(ENSEMBL = sub("\\..*$", "", rownames(expr_matrix_pd_log))) %>%
  inner_join(gene_ids, by = "ENSEMBL") %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  dplyr::select(-ENSEMBL) %>%
  group_by(SYMBOL) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("SYMBOL")

# ------------------------------------------------------------------------------
# 5. 提取 5 个目标基因数据并导出样本表达量表至 results/
# ------------------------------------------------------------------------------
target_order <- c("CLTC", "SLC2A13", "SLC7A14", "ACP5", "RAMP3")

# 检查是否存在缺失基因
missing_genes <- setdiff(target_order, rownames(expr_mapped_df))
if (length(missing_genes) > 0) {
  warning("⚠️ 警告: 以下基因未在表达矩阵中匹配到: ", paste(missing_genes, collapse = ", "))
  target_order <- intersect(target_order, rownames(expr_mapped_df))
}

# 转置并合并表型数据
target_expr <- as.data.frame(t(expr_mapped_df[target_order, ]))
target_expr$Sample_ID <- rownames(target_expr)
target_expr$Group <- sample_info_pd_final$Group

final_plot_data <- target_expr

# 按正确基因顺序导出每样本表达量至相对路径 results/
pd_sample_expr <- final_plot_data %>%
  dplyr::select(Sample_ID, Group, all_of(target_order))

output_csv <- file.path("results", "PD_5_Genes_Sample_Expression.csv")
write.csv(pd_sample_expr, output_csv, row.names = FALSE)
cat("💾 【已完成】样本数据已成功写入相对路径：'", output_csv, "'\n\n", sep = "")

# ------------------------------------------------------------------------------
# 6. 统计真实样本量
# ------------------------------------------------------------------------------
cat("📊 1. 正在统计当前环境中的真实分组样本量...\n")
real_table <- table(final_plot_data$Group)
print(real_table)

n_control <- as.numeric(real_table["Control"])
n_pd <- as.numeric(real_table["PD"])

# ------------------------------------------------------------------------------
# 7. 重塑长格式数据并绘制精致高分学术散点箱线图
# ------------------------------------------------------------------------------
cat("\n🎨 2. 正在应用高分学术样式（实心提亮散点 + 柔和半透明箱线图）...\n")

plot_data_long <- final_plot_data %>%
  pivot_longer(
    cols = all_of(target_order),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  mutate(
    Group = factor(Group, levels = c("Control", "PD")),
    Gene = factor(Gene, levels = target_order) # 锁死分面排列顺序
  )

# 构建基础 ggplot 对象 (X 轴直接为 Group，Control 在前，PD 在后)
p_box_jitter_pd <- ggplot(plot_data_long, aes(x = Group, y = Expression, color = Group)) +
  
  # ① 散点图：纯色实心点 shape=16, size=1.1, width=0.2
  geom_jitter(
    shape = 16, 
    width = 0.2,
    size = 1.1, 
    show.legend = FALSE
  ) + 
  
  # ② 覆盖柔和半透明箱线图：alpha = 0.55
  geom_boxplot(
    aes(fill = Group), 
    width = 0.4, 
    alpha = 0.55, 
    color = "#222222", 
    linewidth = 0.7, 
    outlier.shape = NA,
    show.legend = FALSE
  ) + 
  
  # ③ 分面与坐标轴
  facet_wrap(~ Gene, nrow = 1, scales = "free_y") + 
  
  # ④ 色彩配色：Control (亮蓝 #6495ED)，PD (珊瑚粉红 #F08080)
  scale_fill_manual(values = c("Control" = "#6495ED", "PD" = "#F08080")) + 
  scale_color_manual(values = c("Control" = "#6495ED", "PD" = "#F08080")) + 
  
  # ⑤ 动态匹配安全样本量的 X 轴标签（换行显示样本数）
  scale_x_discrete(labels = c(
    "Control" = paste0("Control\n(n=", n_control, ")"),
    "PD" = paste0("PD\n(n=", n_pd, ")")
  )) +
  
  labs(y = "Expression Level (lcpm)", x = NULL) + 
  theme_bw(base_size = 12) + 
  theme( 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9), 
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.9), 
    strip.text = element_text(face = "bold", size = 13, colour = "black"), 
    axis.text.x = element_text(face = "bold", color = "black", size = 11), 
    axis.text.y = element_text(color = "black", size = 10), 
    axis.title.y = element_text(face = "bold", size = 13, color = "black", margin = ggplot2::margin(r = 10)), 
    panel.spacing = unit(1.3, "lines"),
    legend.position = "none" # 纯两组展示，隐藏冗余图例
  )

# ⑥ 添加两组间的 Wilcoxon P 值检验
p_box_jitter_pd <- p_box_jitter_pd + 
  ggpubr::stat_compare_means( 
    method = "wilcox.test", 
    label = "p.format", 
    label.x = 1.5,
    hjust = 0.5,
    vjust = -0.4,
    size = 3.8
  )

# ------------------------------------------------------------------------------
# 8. 图形窗口刷新与 PDF 导出至 results/
# ------------------------------------------------------------------------------
# 刷新绘图窗口
while (!is.null(dev.list())) { dev.off() }
print(p_box_jitter_pd)

# 保存为高质量发表级 PDF 文件至 relative path results/
output_pdf <- file.path("results", "PD_5_Gene_Exquisite_Jitter_Boxplot.pdf")
pdf(output_pdf, width = 11, height = 4.5)
print(p_box_jitter_pd)
dev.off()

cat("\n🎉 【运行成功结束！】\n")
cat("💾 PDF 箱线图已成功写入相对路径：'", output_pdf, "'\n", sep = "")