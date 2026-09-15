# ==============================================================================
# 数据下载与输入说明
# ==============================================================================
# 请先准备好以下数据文件，并放置在当前工作目录的 data/ 文件夹下：
#
# 1. 请先从 GEO 数据库下载 GSE15222 系列矩阵文件：
#    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15222
#    下载文件 GSE15222_series_matrix.txt.gz 并放置在当前工作目录的 data/ 文件夹下。
#
# 2. 请先从 GEO 数据库下载 GPL2700 平台注释文件：
#    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL2700
#    下载文件 GPL2700.gz 并放置在当前工作目录的 data/ 文件夹下。
#
# 3. 请先从 GEO 数据库下载 GSE68719 系列矩阵文件及标准化表达矩阵：
#    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68719
#    下载文件 GSE68719_series_matrix.txt.gz 和 GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz
#    并放置在当前工作目录的 data/ 文件夹下。
# ==============================================================================

# ==============================================================================
# AD (GSE15222) 与 PD (GSE68719) 差异基因分析及 ggvenn 可视化整合脚本
# ==============================================================================

# 1. 加载所需的依赖包
library(tidyverse)
library(limma)
library(GEOquery)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(ggvenn)

# 确保输出目录存在
if (!dir.exists("data")) dir.create("data")
if (!dir.exists("results")) dir.create("results")

# ==============================================================================
# 第一部分：运行 AD (GSE15222) 差异表达分析
# ==============================================================================
message(">>> [1/4] 正在进行 AD (GSE15222) 差异表达分析...")

file_path_ad <- file.path("data", "GSE15222_series_matrix.txt.gz")
if (!file.exists(file_path_ad)) stop("❌ 未找到 data/GSE15222_series_matrix.txt.gz！")

gse_ad <- getGEO(filename = file_path_ad, GSEMatrix = TRUE, getGPL = FALSE)
sample_info_ad <- pData(gse_ad)

sample_info_ad <- sample_info_ad %>%
  dplyr::mutate(Group = case_when(
    str_detect(tolower(description), "control") ~ "Control",
    str_detect(tolower(description), "alzheimer") ~ "AD",
    TRUE ~ NA_character_
  ))

valid_samples_ad <- !is.na(sample_info_ad$Group)
sample_info_ad <- sample_info_ad[valid_samples_ad, ]
expr_raw_ad <- exprs(gse_ad)[, valid_samples_ad]

# GPL2700 探针注释映射
gpl_file_ad <- file.path("data", "GPL2700.gz")
if (!file.exists(gpl_file_ad)) stop("❌ 未找到 data/GPL2700.gz！")

gpl_ad <- GEOquery::getGEO(filename = gpl_file_ad)
gpl_meta_ad <- Table(gpl_ad)

probe_to_acc <- gpl_meta_ad[, c("ID", "GB_ACC")]
probe_to_acc <- probe_to_acc[probe_to_acc$GB_ACC != "" & !is.na(probe_to_acc$GB_ACC), ]
probe_to_acc$GB_ACC_clean <- sub("\\..*$", "", probe_to_acc$GB_ACC)

mapping_ad <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(probe_to_acc$GB_ACC_clean),
  columns = "SYMBOL",
  keytype = "ACCNUM"
)

probe2symbol_ad <- merge(probe_to_acc, mapping_ad, by.x = "GB_ACC_clean", by.y = "ACCNUM") %>%
  dplyr::select(ID, SYMBOL) %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  dplyr::distinct()

# 标准化与预处理
if (max(expr_raw_ad, na.rm = TRUE) > 100) {
  expr_raw_ad[expr_raw_ad <= 0] <- 1
  expr_matrix_ad <- log2(expr_raw_ad)
} else {
  expr_matrix_ad <- expr_raw_ad
}

non_zero_var_ad <- apply(expr_matrix_ad, 1, var, na.rm = TRUE) > 0
expr_matrix_ad <- normalizeBetweenArrays(expr_matrix_ad[non_zero_var_ad, ])

# Limma 差异分析
group_list_ad <- factor(sample_info_ad$Group, levels = c("Control", "AD"))
design_ad <- model.matrix(~group_list_ad)
colnames(design_ad) <- c("Intercept", "AD_vs_Control")

fit_ad <- lmFit(expr_matrix_ad, design_ad)
fit_ad <- eBayes(fit_ad)
deg_ad <- topTable(fit_ad, coef = "AD_vs_Control", number = Inf)

# 注释与分类 (按 |logFC| > 0.5 & P < 0.05)
deg_ad_annotated <- deg_ad %>%
  rownames_to_column(var = "ID") %>%
  inner_join(probe2symbol_ad, by = "ID") %>%
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  dplyr::arrange(P.Value) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

ad_p_col <- if("P.Value" %in% colnames(deg_ad_annotated)) "P.Value" else "adj.P.Val"

ad_deg_classified <- deg_ad_annotated %>%
  dplyr::mutate(
    Group = case_when(
      logFC > 0.5 & .data[[ad_p_col]] < 0.05 ~ "Up",
      logFC < -0.5 & .data[[ad_p_col]] < 0.05 ~ "Down",
      TRUE ~ "None"
    )
  )

# 提取 AD 基因列表
ad_up_genes   <- ad_deg_classified %>% dplyr::filter(Group == "Up") %>% dplyr::pull(SYMBOL) %>% unique()
ad_down_genes <- ad_deg_classified %>% dplyr::filter(Group == "Down") %>% dplyr::pull(SYMBOL) %>% unique()
ad_all_degs   <- unique(c(ad_up_genes, ad_down_genes))


# ==============================================================================
# 第二部分：运行 PD (GSE68719) 差异表达分析
# ==============================================================================
message(">>> [2/4] 正在进行 PD (GSE68719) 差异表达分析...")

pd_series_file <- file.path("data", "GSE68719_series_matrix.txt.gz")
norm_counts_file <- file.path("data", "GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz")

if (!file.exists(pd_series_file) || !file.exists(norm_counts_file)) {
  stop("❌ 未在 data/ 目录下找到 GSE68719 相关输入文件！")
}

gse_pd <- getGEO(filename = pd_series_file, GSEMatrix = TRUE, getGPL = FALSE)
sample_info_pd <- pData(gse_pd)

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

pd_expr_df <- read.table(gzfile(norm_counts_file), header = TRUE, sep = "\t", check.names = FALSE)
rownames(pd_expr_df) <- pd_expr_df[, 1]
expr_matrix_pd_real <- as.matrix(pd_expr_df[, -1])

common_pd_samples <- intersect(colnames(expr_matrix_pd_real), rownames(valid_pd_samples))
expr_matrix_pd_final <- expr_matrix_pd_real[, common_pd_samples]
sample_info_pd_final <- valid_pd_samples[common_pd_samples, ]

expr_numeric <- apply(expr_matrix_pd_final, 2, as.numeric)
rownames(expr_numeric) <- rownames(expr_matrix_pd_final)

if (max(expr_numeric, na.rm = TRUE) > 50) {
  expr_matrix_pd_log <- log2(expr_numeric + 1)
} else {
  expr_matrix_pd_log <- expr_numeric
}

group_factor_pd <- factor(sample_info_pd_final$Group, levels = c("Control", "PD"))
design_pd <- model.matrix(~ group_factor_pd)
colnames(design_pd) <- c("Intercept", "PD_vs_Control")

fit_pd <- lmFit(expr_matrix_pd_log, design_pd)
fit_pd <- eBayes(fit_pd)

deg_pd_all <- topTable(fit_pd, coef = "PD_vs_Control", number = Inf, adjust.method = "BH") %>%
  rownames_to_column("Ensembl_ID")

deg_pd_all$Ensembl_clean <- sub("\\..*$", "", deg_pd_all$Ensembl_ID)

gene_ids_pd <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = deg_pd_all$Ensembl_clean,
  columns = "SYMBOL",
  keytype = "ENSEMBL"
)

deg_pd_mapped <- left_join(deg_pd_all, gene_ids_pd, by = c("Ensembl_clean" = "ENSEMBL")) %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  distinct(SYMBOL, .keep_all = TRUE)

pd_p_col <- if("P.Value" %in% colnames(deg_pd_mapped)) "P.Value" else "adj.P.Val"

pd_deg_classified <- deg_pd_mapped %>%
  mutate(
    Group = case_when(
      logFC > 0.5 & .data[[pd_p_col]] < 0.05 ~ "Up",
      logFC < -0.5 & .data[[pd_p_col]] < 0.05 ~ "Down",
      TRUE ~ "None"
    )
  )

# 提取 PD 基因列表
pd_up_genes   <- pd_deg_classified %>% filter(Group == "Up") %>% pull(SYMBOL) %>% unique()
pd_down_genes <- pd_deg_classified %>% filter(Group == "Down") %>% pull(SYMBOL) %>% unique()
pd_all_degs   <- unique(c(pd_up_genes, pd_down_genes))


# ==============================================================================
# 第三部分：使用 ggvenn 绘制 3 张高颜值 Venn 图
# ==============================================================================
message(">>> [3/4] 正在使用 ggvenn 绘制并保存韦恩图...")

# 通用绘图保存函数
plot_and_save_venn <- function(data_list, title, filename_pdf, fill_colors) {
  p <- ggvenn(
    data_list,
    fill_color = fill_colors,
    stroke_size = 0,
    text_size = 7,
    set_name_size = 7,
    show_percentage = FALSE
  ) +
    coord_fixed() +
    labs(title = title) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 15)),
      plot.margin = margin(20, 20, 20, 20)
    )
  
  print(p)
  ggsave(file.path("results", filename_pdf), plot = p, width = 6, height = 5.5, dpi = 300)
  return(p)
}

# 清空 Plots 窗口
while (!is.null(dev.list())) { dev.off() }

# 1️⃣ 图一：所有变化差异基因 (All DEGs)
list_all_data <- list(AD_All = ad_all_degs, PD_All = pd_all_degs)
p_all <- plot_and_save_venn(
  list_all_data,
  title = "All Changed DEGs: AD vs PD",
  filename_pdf = "Venn_1_All_DEGs_Fixed.pdf",
  fill_colors = c("#CBD5E8", "#FDCDAC")
)

# 2️⃣ 图二：上调差异基因 (Up-regulated DEGs)
list_up_data <- list(AD_Up = ad_up_genes, PD_Up = pd_up_genes)
p_up <- plot_and_save_venn(
  list_up_data,
  title = "Up-regulated DEGs: AD vs PD",
  filename_pdf = "Venn_2_Up_DEGs_Fixed.pdf",
  fill_colors = c("#FBB4AE", "#B3CDE3")
)

# 3️⃣ 图三：下调差异基因 (Down-regulated DEGs)
list_down_data <- list(AD_Down = ad_down_genes, PD_Down = pd_down_genes)
p_down <- plot_and_save_venn(
  list_down_data,
  title = "Down-regulated DEGs: AD vs PD",
  filename_pdf = "Venn_3_Down_DEGs_Fixed.pdf",
  fill_colors = c("#CCEBC5", "#DECBE4")
)


# ==============================================================================
# 第四部分：数据计算与统计打印导出
# ==============================================================================
message(">>> [4/4] 导出分析数据...")

overlap_all  <- intersect(ad_all_degs, pd_all_degs)
overlap_up   <- intersect(ad_up_genes, pd_up_genes)
overlap_down <- intersect(ad_down_genes, pd_down_genes)

cat("\n====================================================================\n")
cat("📊 【AD vs PD 数据交集核对统计】\n")
cat("--------------------------------------------------------------------\n")
cat("● All DEGs 交集基因数:  ", length(overlap_all), " 个\n")
cat("● Up DEGs 交集基因数:   ", length(overlap_up), " 个\n")
cat("● Down DEGs 交集基因数: ", length(overlap_down), " 个\n")
cat("====================================================================\n\n")

# 保存结果 CSV 至 results/ 目录
write.csv(data.frame(Symbol = overlap_all), file.path("results", "Shared_All_DEGs.csv"), row.names = FALSE)
write.csv(data.frame(Symbol = overlap_up), file.path("results", "Shared_Up_Genes.csv"), row.names = FALSE)
write.csv(data.frame(Symbol = overlap_down), file.path("results", "Shared_Down_Genes.csv"), row.names = FALSE)

cat(">> 运行全部完成！PDF 图像与 CSV 交集表已成功存入 results/ 目录下。\n")