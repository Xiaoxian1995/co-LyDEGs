# ==============================================================================
# 项目名称: GSE15222 差异表达分析 & 目标基因散点箱线图 (AD vs Control)
# 数据说明:
#   请先从 GEO 数据库下载数据文件，并放置在当前工作目录的 data/ 文件夹下：
#     - data/GSE15222_series_matrix.txt.gz
#     - data/GPL2700.gz
# ==============================================================================

# 1. 加载所需的依赖包
library(tidyverse)
library(limma)
library(GEOquery)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# 自动建立标准的相对路径目录结构
if (!dir.exists("data")) dir.create("data")
if (!dir.exists("results")) dir.create("results")

# ------------------------------------------------------------------------------
# 2. 读取 GEO 数据与提取表型矩阵
# ------------------------------------------------------------------------------
file_path <- file.path("data", "GSE15222_series_matrix.txt.gz")
if (!file.exists(file_path)) {
  stop("❌ 请先从 GEO 下载 GSE15222 数据文件并放置在当前工作目录的 data/ 文件夹下！")
}

gse <- getGEO(filename = file_path, GSEMatrix = TRUE, getGPL = FALSE)
sample_info <- pData(gse)

# 基于 description 列提取样本分组 (Control vs Alzheimer/AD)
sample_info <- sample_info %>%
  dplyr::mutate(Group = case_when(
    str_detect(tolower(description), "control") ~ "Control",
    str_detect(tolower(description), "alzheimer") ~ "AD",
    TRUE ~ NA_character_
  ))

# 剔除未能分组的样本并对齐表达矩阵
valid_samples <- !is.na(sample_info$Group)
sample_info <- sample_info[valid_samples, ]
expr_raw <- exprs(gse)[, valid_samples]

group_list <- factor(sample_info$Group, levels = c("Control", "AD"))

cat("==================================================\n")
cat("📊 样本分组统计：\n")
print(table(group_list))
cat("==================================================\n\n")

# ------------------------------------------------------------------------------
# 3. GPL2700 本地注释表生成 (探针 ID -> Gene Symbol)
# ------------------------------------------------------------------------------
gpl_file <- file.path("data", "GPL2700.gz")
if (!file.exists(gpl_file)) {
  stop("❌ 请先从 GEO 下载 GPL2700 注释文件并放置在当前工作目录的 data/ 文件夹下！")
}

gpl <- GEOquery::getGEO(filename = gpl_file)
gpl_meta <- Table(gpl)

probe_to_acc <- gpl_meta[, c("ID", "GB_ACC")]
probe_to_acc <- probe_to_acc[probe_to_acc$GB_ACC != "" & !is.na(probe_to_acc$GB_ACC), ]

# 去除登录号版本号后缀 (如 NM_007162.2 -> NM_007162)
probe_to_acc$GB_ACC_clean <- sub("\\..*$", "", probe_to_acc$GB_ACC)

# 从 org.Hs.eg.db 数据库比对获取 Gene Symbol
mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(probe_to_acc$GB_ACC_clean),
  columns = "SYMBOL",
  keytype = "ACCNUM"
)

# 显式调用 dplyr::select 解决函数名冲突
probe2symbol <- merge(probe_to_acc, mapping, by.x = "GB_ACC_clean", by.y = "ACCNUM") %>%
  dplyr::select(ID, SYMBOL) %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  dplyr::distinct()

# ------------------------------------------------------------------------------
# 4. 数据预处理：Log2 缩放与分位数归一化 (Normalize Between Arrays)
# ------------------------------------------------------------------------------
if (max(expr_raw, na.rm = TRUE) > 100) {
  expr_raw[expr_raw <= 0] <- 1
  expr_matrix_correct <- log2(expr_raw)
  message(">> 检测到原始数据未取 Log，已自动完成 Log2 转换。")
} else {
  expr_matrix_correct <- expr_raw
}

# 过滤在所有样本中表达量完全不变化的零方差探针
non_zero_var <- apply(expr_matrix_correct, 1, var, na.rm = TRUE) > 0
expr_matrix_correct <- expr_matrix_correct[non_zero_var, ]

# 使用 limma 进行分位数归一化
expr_matrix_correct <- normalizeBetweenArrays(expr_matrix_correct)

# ------------------------------------------------------------------------------
# 5. 构建基因级别的表达矩阵 (Gene-level Expression Matrix)
# ------------------------------------------------------------------------------
cat(">> 正在整合探针表达式并转换为基因级别表达矩阵...\n")

# 将探针表达矩阵与注释结合，多个探针映射到同一基因时取均值
gene_expr_df <- as.data.frame(expr_matrix_correct) %>%
  rownames_to_column("ID") %>%
  inner_join(probe2symbol, by = "ID") %>%
  dplyr::select(-ID) %>%
  group_by(SYMBOL) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("SYMBOL")

# ------------------------------------------------------------------------------
# 6. 提取目标 5 个基因的表达量并构建箱线图数据集
# ------------------------------------------------------------------------------
target_order <- c("CLTC", "SLC2A13", "SLC7A14", "ACP5", "RAMP3")

# 检查目标基因是否全部在表达矩阵中
missing_genes <- setdiff(target_order, rownames(gene_expr_df))
if (length(missing_genes) > 0) {
  warning("⚠️ 警告: 以下基因未在表达矩阵中匹配到: ", paste(missing_genes, collapse = ", "))
  target_order <- intersect(target_order, rownames(gene_expr_df))
}

# 提取目标基因表达量 (转置为 样本 x 基因 矩阵)
validation_set <- as.data.frame(t(gene_expr_df[target_order, ]))
validation_set$Group <- sample_info$Group

# 锁死分组顺序 (Control 在前，AD 在后)
boxplot_data <- validation_set[, c(target_order, "Group")]
boxplot_data$Group <- factor(boxplot_data$Group, levels = c("Control", "AD"))

# 转换为长矩阵格式
melted_data <- boxplot_data %>%
  pivot_longer(cols = all_of(target_order), names_to = "Gene", values_to = "Expression")

# 锁死基因分面展示顺序
melted_data$Gene <- factor(melted_data$Gene, levels = target_order)

# ------------------------------------------------------------------------------
# 7. 绘制无图例精致散点箱线图 (Jitter Boxplot)
# ------------------------------------------------------------------------------
cat("🎨 正在为您生成精致学术箱线图...\n")

p_box_jitter_ad <- ggplot(melted_data, aes(x = Group, y = Expression, color = Group)) +
  
  # ① 散点图：纯色实心点, shape=16, size=1.1
  geom_jitter(aes(color = Group), shape = 16, width = 0.2, size = 1.1, show.legend = FALSE) + 
  
  # ② 柔和半透明箱线图：alpha = 0.55, 去掉异常值
  geom_boxplot(aes(fill = Group), width = 0.4, alpha = 0.55, color = "#222222", linewidth = 0.7, outlier.shape = NA, show.legend = FALSE) + 
  
  # ③ 分面与坐标轴
  facet_wrap(~ Gene, nrow = 1, scales = "free_y") + 
  
  # ④ 经典学术配色：Control (亮蓝 #6495ED)，AD (珊瑚粉红 #F08080)
  scale_fill_manual(values = c("Control" = "#6495ED", "AD" = "#F08080")) + 
  scale_color_manual(values = c("Control" = "#6495ED", "AD" = "#F08080")) + 
  labs(y = "Expression level", x = NULL) + 
  theme_bw(base_size = 12) + 
  theme( 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9), 
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.9), 
    strip.text = element_text(face = "bold", size = 13, colour = "black"), 
    axis.text.x = element_text(face = "bold", color = "black", size = 11.5), 
    axis.text.y = element_text(color = "black", size = 10), 
    axis.title.y = element_text(face = "bold", size = 13, color = "black", margin = ggplot2::margin(r = 10)), 
    panel.spacing = unit(1.3, "lines"),
    legend.position = "none" # 关掉 Group 图例
  )

# ⑤ 添加 Wilcoxon 检验 P 值
p_box_jitter_ad <- p_box_jitter_ad + 
  ggpubr::stat_compare_means( 
    method = "wilcox.test", 
    label = "p.format", 
    label.x = 1.5, 
    hjust = 0.5, 
    size = 3.8, 
    vjust = -0.4
  )

# ------------------------------------------------------------------------------
# 8. 图形设备刷新与 PDF 导出
# ------------------------------------------------------------------------------
# 清空现有图形窗口并在图形界面输出
while (!is.null(dev.list())) { dev.off() }
print(p_box_jitter_ad)

# 保存高分辨率 PDF 文件至 results/ 文件夹
output_pdf <- file.path("results", "AD_5_Gene_Exquisite_Jitter_Boxplot.pdf")
pdf(output_pdf, width = 11, height = 4.5)
print(p_box_jitter_ad)
dev.off()

cat("\n🎉 【运行成功完成！】\n")
cat("💾 精美散点箱线图已写入相对路径文件：'", output_pdf, "'\n", sep = "")