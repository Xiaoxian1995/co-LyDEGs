# ==============================================================================
# 🎯 [GSE205450 纯本地提取 + 表头自解析临床分组 + 指定名称文件导出]
# ==============================================================================
# 📝 备注：运行本脚本前，请提前下载以下 GEO 数据文件并放置于脚本同级目录（./）：
# 
# 1. GSE205450 lcpm 表达矩阵文件 (GSE205450_lcpm.table.txt.gz):
#    下载链接: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE205450&format=file&file=GSE205450%5Flcpm%2Etable%2Etxt%2Egz
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# ------------------------------------------------------------------------------
# 1. 检查本地表达矩阵文件
# ------------------------------------------------------------------------------
base_dir <- "./"
local_file <- file.path(base_dir, "GSE205450_lcpm.table.txt.gz")

if (!file.exists(local_file)) {
  stop("⚠️ 错误：未在当前工作目录下找到 'GSE205450_lcpm.table.txt.gz'，请检查文件路径！")
}

# ------------------------------------------------------------------------------
# 2. 从本地文件表头自动解析 clinical_info（无需联网）
# ------------------------------------------------------------------------------
cat("🛠️ 1. 正在读取本地矩阵表头，自动解析分组信息 (Control/PD & Caudate/Putamen)...\n")

con <- gzfile(local_file, "r")
header <- readLines(con, n = 1)
raw_sample_ids <- strsplit(header, "\t")[[1]][-1] # 扣除第一列 "ID"
sample_ids <- gsub('"', '', raw_sample_ids)       # 清理双引号

# 从 Sample Title 中提取 Group 和 Center 信息
clinical_info <- data.frame(
  title = sample_ids,
  stringsAsFactors = FALSE
) %>%
  mutate(
    # 解析分组: Control vs PD
    Group = case_when(
      grepl("control|ctrl|normal|N_", title, ignore.case = TRUE) ~ "Control",
      grepl("PD|parkinson", title, ignore.case = TRUE) ~ "PD",
      TRUE ~ "PD"
    ),
    # 解析脑区: Caudate vs Putamen
    Center = case_when(
      grepl("caudate|cau", title, ignore.case = TRUE) ~ "Caudate",
      grepl("putamen|put", title, ignore.case = TRUE) ~ "Putamen",
      TRUE ~ "Caudate"
    )
  )

cat("✅ 离线临床表型提取完成！前 5 行预览：\n")
print(head(clinical_info, 5))

# ------------------------------------------------------------------------------
# 3. 从本地 lcpm 矩阵中精准提取 5 个核心基因
# ------------------------------------------------------------------------------
target_order <- c("CLTC", "SLC2A13", "SLC7A14", "ACP5", "RAMP3")
extracted_rows <- list()

cat("\n🛠️ 2. 正在从本地 lcpm 矩阵逐行检索 5 个核心基因...\n")

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
cat("\n✨ 5 个核心基因表达量提取完毕！\n")

# ------------------------------------------------------------------------------
# 4. 组装数据，显示详细样本量，保存为 GSE205450_5_Genes_Sample_Expression.csv
# ------------------------------------------------------------------------------
cat("\n🔍 3. 正在组装表达矩阵并导出指定 CSV...\n")

expr_5_genes_df <- as.data.frame(extracted_rows)
expr_5_genes_df$Matrix_Title <- sample_ids

# 精准合并表达数据与临床表型
final_plot_data <- inner_join(
  clinical_info %>% dplyr::select(title, Group, Center),
  expr_5_genes_df,
  by = c("title" = "Matrix_Title")
)

# 💾 保存指定名称的 CSV 矩阵
csv_output_name <- file.path(base_dir, "GSE205450_5_Genes_Sample_Expression.csv")
write.csv(final_plot_data, csv_output_name, row.names = FALSE)
cat(paste0("💾 表达矩阵已保存为: '", csv_output_name, "'\n"))

# 📌 显式计算并打印详细样本量统计
real_table <- table(final_plot_data$Group, final_plot_data$Center)

cat("\n==================================================\n")
cat("📊 【GSE205450 数据集 - 各脑区与临床分组样本量统计】\n")
cat("--------------------------------------------------\n")
print(real_table)
cat("--------------------------------------------------\n")
cat("  - Caudate (尾状核) : Control =", real_table["Control","Caudate"], "例 | PD =", real_table["PD","Caudate"], "例\n")
cat("  - Putamen (壳核)   : Control =", real_table["Control","Putamen"], "例 | PD =", real_table["PD","Putamen"], "例\n")
cat("  - 总有效分析样本数 :", nrow(final_plot_data), "例\n")
cat("==================================================\n\n")

# ------------------------------------------------------------------------------
# 5. 重构长格式数据并绘制精修学术箱线图
# ------------------------------------------------------------------------------
cat("🎨 4. 正在构建精致高分双组织箱线图...\n")

plot_data_long <- final_plot_data %>%
  pivot_longer(
    cols = all_of(target_order),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  mutate(
    Group = factor(Group, levels = c("Control", "PD")),
    Center = factor(Center, levels = c("Caudate", "Putamen")),
    Gene = factor(Gene, levels = target_order) # 强行锁定 5 基因横向排列顺序
  )

p_box_jitter_pd <- ggplot(plot_data_long, aes(x = Center, y = Expression, color = Group)) +
  # 纯色无遮挡散点
  geom_jitter(
    position = position_dodge(0.75),
    shape = 16,
    size = 1.1,
    show.legend = FALSE
  ) +
  # 半透明箱线图
  geom_boxplot(
    aes(fill = Group),
    width = 0.4,
    alpha = 0.55,
    color = "#222222",
    linewidth = 0.7,
    outlier.shape = NA,
    position = position_dodge(0.75)
  ) +
  # 5 基因分面
  facet_wrap(~ Gene, nrow = 1, scales = "free_y") +
  # 学术配色: Control (#6495ED) vs PD (#F08080)
  scale_fill_manual(values = c("Control" = "#6495ED", "PD" = "#F08080")) +
  scale_color_manual(values = c("Control" = "#6495ED", "PD" = "#F08080")) +
  # X 轴包含样本量计数的标签
  scale_x_discrete(labels = c(
    "Caudate" = paste0("Caudate\n(Ctrl:", real_table["Control","Caudate"], ", PD:", real_table["PD","Caudate"], ")"),
    "Putamen" = paste0("Putamen\n(Ctrl:", real_table["Control","Putamen"], ", PD:", real_table["PD","Putamen"], ")")
  )) +
  labs(y = "Expression Level (lcpm)", x = NULL) +
  theme_bw(base_size = 12) +
  theme(
    text             = element_text(family = "Helvetica"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.9),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.9),
    strip.text       = element_text(face = "bold", size = 13, colour = "black"),
    axis.text.x      = element_text(face = "bold", color = "black", size = 10.5),
    axis.text.y      = element_text(color = "black", size = 10),
    axis.title.y     = element_text(face = "bold", size = 13, color = "black", margin = ggplot2::margin(r = 10)),
    panel.spacing    = unit(1.3, "lines"),
    legend.position  = "top",
    legend.title     = element_blank()
  )

# 添加 Wilcoxon 显著性检验
if (requireNamespace("ggpubr", quietly = TRUE)) {
  p_box_jitter_pd <- p_box_jitter_pd +
    ggpubr::stat_compare_means(
      aes(group = Group),
      method      = "wilcox.test",
      label       = "p.format",
      label.y.npc = 0.93,
      size        = 3.5,
      show.legend = FALSE
    )
}

# ------------------------------------------------------------------------------
# 6. 保存为指定名称的矢量 PDF (GSE205450_5_Genes_Boxplot.pdf)
# ------------------------------------------------------------------------------
pdf_output_name <- file.path(base_dir, "GSE205450_5_Genes_Boxplot.pdf")

pdf(
  file        = pdf_output_name,
  width       = 11,
  height      = 4.17,       # 矢量紧凑学术高度 (AI活字保护)
  family      = "Helvetica",
  useDingbats = FALSE
)

print(p_box_jitter_pd)
dev.off()

cat("\n🎉 【纯本地提取与渲染成功完成！】\n")
cat("💾 1. CSV 文件已导出至 : '", csv_output_name, "'\n", sep = "")
cat("💾 2. PDF 文件已导出至 : '", pdf_output_name, "'\n", sep = "")