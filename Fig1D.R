# ==============================================================================
# 项目名称: AD (GSE15222) vs PD (GSE68719) vs Lysosome 基因集 3维 Venn 分析
# 数据说明:
# 请先从 NCBI GEO 数据库及本地准备好数据文件，并放置在当前工作目录的 data/ 文件夹下：
#   - data/GSE15222_series_matrix.txt.gz
#   - data/GPL2700.gz
#   - data/GSE68719_series_matrix.txt.gz
#   - data/GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz
#   - data/lysosome gene.csv
# ==============================================================================

# 1. 加载所需的依赖包
library(tidyverse)
library(limma)
library(GEOquery)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(ggvenn)

# 自动建立标准的相对路径目录结构
if (!dir.exists("data")) dir.create("data")
if (!dir.exists("results")) dir.create("results")

# ==============================================================================
# 第一部分：运行 AD (GSE15222) 差异表达分析
# ==============================================================================
message(">>> [1/5] 正在进行 AD (GSE15222) 差异表达分析...")

file_path_ad <- file.path("data", "GSE15222_series_matrix.txt.gz")
if (!file.exists(file_path_ad)) stop("❌ 未找到 data/GSE15222_series_matrix.txt.gz 文件！")

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
if (!file.exists(gpl_file_ad)) stop("❌ 未找到 data/GPL2700.gz 文件！")

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

# 注释与分类 (筛选条件: |logFC| > 0.5 & P < 0.05)
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
message(">>> [2/5] 正在进行 PD (GSE68719) 差异表达分析...")

pd_series_file <- file.path("data", "GSE68719_series_matrix.txt.gz")
norm_counts_file <- file.path("data", "GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz")

if (!file.exists(pd_series_file) || !file.exists(norm_counts_file)) {
  stop("❌ 未在 data/ 目录下找到 GSE68719 相关数据文件！")
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
# 第三部分：读取并清洗本地溶酶体基因集 (lysosome gene.csv)
# ==============================================================================
message(">>> [3/5] 正在读取本地溶酶体基因集...")

lyso_file_path <- file.path("data", "lysosome gene.csv")
if (!file.exists(lyso_file_path)) stop("❌ 未找到 data/lysosome gene.csv 文件！")

lyso_df <- read.csv(lyso_file_path, stringsAsFactors = FALSE, check.names = FALSE)

# 自动定位匹配文件内的基因 Symbol 列
symbol_col <- colnames(lyso_df)[grep("symbol|gene", colnames(lyso_df), ignore.case = TRUE)][1]

if (is.na(symbol_col)) {
  lyso_genes <- unique(as.character(lyso_df[[1]]))
} else {
  lyso_genes <- unique(as.character(lyso_df[[symbol_col]]))
}

lyso_genes <- lyso_genes[!is.na(lyso_genes) & lyso_genes != ""]
cat("成功提取溶酶体基因集数量：", length(lyso_genes), " 个\n\n")


# ==============================================================================
# 第四部分：使用 ggvenn 绘制 3 维 (AD, PD, Lysosome) 韦恩图
# ==============================================================================
message(">>> [4/5] 正在绘制并保存 3 维韦恩图...")

# 绘图及文件输出函数
plot_and_save_venn3 <- function(data_list, title, filename_pdf, fill_colors) {
  p <- ggvenn(
    data_list,
    fill_color = fill_colors,
    stroke_size = 0.2,
    stroke_color = "white",
    text_size = 5.5,
    set_name_size = 5.5,
    show_percentage = FALSE
  ) +
    coord_fixed() +
    labs(title = title) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 15, face = "bold", margin = margin(b = 15)),
      plot.margin = margin(20, 20, 20, 20)
    )
  
  print(p)
  ggsave(file.path("results", filename_pdf), plot = p, width = 6.5, height = 6, dpi = 300)
  return(p)
}

# 清理当前图形窗口
while (!is.null(dev.list())) { dev.off() }

# 1. 图一：所有有变化的基因 (AD DEGs vs PD DEGs vs Lysosome)
list_all_3way <- list(
  `AD DEGs` = ad_all_degs,
  `PD DEGs` = pd_all_degs,
  `Lysosome` = lyso_genes
)
p_all_3way <- plot_and_save_venn3(
  list_all_3way,
  title = "All Changed DEGs & Lysosome Genes",
  filename_pdf = "Venn_3Way_1_All_DEGs.pdf",
  fill_colors = c("#8DA0CB", "#FC8D62", "#66C2A5")
)

# 2. 图二：共同上调基因 (AD Up vs PD Up vs Lysosome)
list_up_3way <- list(
  `AD Up` = ad_up_genes,
  `PD Up` = pd_up_genes,
  `Lysosome` = lyso_genes
)
p_up_3way <- plot_and_save_venn3(
  list_up_3way,
  title = "Shared Up-regulated DEGs & Lysosome Genes",
  filename_pdf = "Venn_3Way_2_Up_DEGs.pdf",
  fill_colors = c("#FBB4AE", "#B3CDE3", "#DECBE4")
)

# 3. 图三：共同下调基因 (AD Down vs PD Down vs Lysosome)
list_down_3way <- list(
  `AD Down` = ad_down_genes,
  `PD Down` = pd_down_genes,
  `Lysosome` = lyso_genes
)
p_down_3way <- plot_and_save_venn3(
  list_down_3way,
  title = "Shared Down-regulated DEGs & Lysosome Genes",
  filename_pdf = "Venn_3Way_3_Down_DEGs.pdf",
  fill_colors = c("#CCEBC5", "#FED9A6", "#FFFFCC")
)


# ==============================================================================
# 第五部分：三方交集数据导出 CSV
# ==============================================================================
message(">>> [5/5] 计算三方重叠基因并导出 CSV 文件...")

overlap_3way_all  <- Reduce(intersect, list(ad_all_degs, pd_all_degs, lyso_genes))
overlap_3way_up   <- Reduce(intersect, list(ad_up_genes, pd_up_genes, lyso_genes))
overlap_3way_down <- Reduce(intersect, list(ad_down_genes, pd_down_genes, lyso_genes))

cat("\n====================================================================\n")
cat("📊 【AD x PD x Lysosome 三方交集统计结果】\n")
cat("--------------------------------------------------------------------\n")
cat("1. 共同有变化的溶酶体基因数 (Shared All):  ", length(overlap_3way_all), " 个\n")
cat("2. 共同显著上调的溶酶体基因数 (Shared Up):   ", length(overlap_3way_up), " 个\n")
cat("3. 共同显著下调的溶酶体基因数 (Shared Down): ", length(overlap_3way_down), " 个\n")
cat("====================================================================\n\n")

# 保存结果至 results/ 文件夹
write.csv(data.frame(Symbol = overlap_3way_all), file.path("results", "Shared_3Way_All_Lysosome_DEGs.csv"), row.names = FALSE)
write.csv(data.frame(Symbol = overlap_3way_up), file.path("results", "Shared_3Way_Up_Lysosome_Genes.csv"), row.names = FALSE)
write.csv(data.frame(Symbol = overlap_3way_down), file.path("results", "Shared_3Way_Down_Lysosome_Genes.csv"), row.names = FALSE)

cat(">> 运行成功完成！所有 3-way Venn PDF 图像及对应交集表格已成功存入 results/ 目录中。\n")