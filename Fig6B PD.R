# ==============================================================================
# 🎯 [GSE205450 本地提取 + 5-Gene 随机森林特征重要性 (ACP5/RAMP3 上调红色版)]
# ==============================================================================
# 📝 备注：运行本脚本前，请提前下载以下 GEO 数据文件并放置于脚本同级目录（./）：
# 
# 1. GSE205450 lcpm 表达矩阵文件 (GSE205450_lcpm.table.txt.gz):
#    下载链接: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE205450&format=file&file=GSE205450%5Flcpm%2Etable%2Etxt%2Egz
# ==============================================================================

library(dplyr)
library(tidyr)
library(randomForest)
library(ggplot2)

# ------------------------------------------------------------------------------
# 1. 检查并读取本地矩阵数据
# ------------------------------------------------------------------------------
base_dir <- "./"
local_file <- file.path(base_dir, "GSE205450_lcpm.table.txt.gz")

if (!file.exists(local_file)) {
  stop("⚠️ 错误：未在当前工作目录下找到 'GSE205450_lcpm.table.txt.gz'，请检查文件路径！")
}

cat("🛠️ 1. 正在读取本地矩阵表头，解析临床分组...\n")

con <- gzfile(local_file, "r")
header <- readLines(con, n = 1)
raw_sample_ids <- strsplit(header, "\t")[[1]][-1]
sample_ids <- gsub('"', '', raw_sample_ids)

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

genes_list <- c("CLTC", "RAMP3", "SLC7A14", "ACP5", "SLC2A13")
extracted_rows <- list()

cat("🛠️ 2. 正在从本地矩阵检索 5 个核心基因表达量...\n")
line_count <- 0
pb <- txtProgressBar(min = 0, max = 30000, style = 3)

while (TRUE) {
  line <- readLines(con, n = 1)
  if (length(line) == 0) break
  line_count <- line_count + 1
  if (line_count %% 500 == 0) setTxtProgressBar(pb, min(line_count, 30000))
  
  parts <- strsplit(line, "\t")[[1]]
  gene_name <- gsub('"', '', parts[1])
  if (gene_name %in% genes_list) {
    extracted_rows[[gene_name]] <- as.numeric(parts[-1])
  }
}
close(con)
setTxtProgressBar(pb, 30000)
close(pb)

expr_5_genes_df <- as.data.frame(extracted_rows)
expr_5_genes_df$Matrix_Title <- sample_ids

final_plot_data <- inner_join(
  clinical_info %>% dplyr::select(title, Group, Center),
  expr_5_genes_df,
  by = c("title" = "Matrix_Title")
)

# ------------------------------------------------------------------------------
# 2. 随机森林计算与指定上下调方向 (ACP5 与 RAMP3 明确为 Up-regulated)
# ------------------------------------------------------------------------------
cat("\n🌲 3. 构建随机森林模型...\n")

# ---- Caudate 计算 ----
set.seed(42)
caudate_df <- final_plot_data %>% 
  filter(Center == "Caudate") %>% 
  mutate(Group = factor(Group, levels = c("Control", "PD")))

x_matrix_cau <- as.matrix(caudate_df[, genes_list])
mode(x_matrix_cau) <- "numeric"

rf_cau <- randomForest(x = x_matrix_cau, y = caudate_df$Group, importance = TRUE, ntree = 1000)

imp_cau <- as.data.frame(importance(rf_cau)) %>% 
  mutate(
    Gene = rownames(.), 
    Region = "Caudate", 
    Importance = MeanDecreaseGini, 
    # 纠正点：ACP5 和 RAMP3 为 Up-regulated
    Regulation = if_else(Gene %in% c("ACP5", "RAMP3"), "Up-regulated (in AD/PD)", "Down-regulated (in AD/PD)")
  )

# ---- Putamen 计算 ----
set.seed(42)
putamen_df <- final_plot_data %>% 
  filter(Center == "Putamen") %>% 
  mutate(Group = factor(Group, levels = c("Control", "PD")))

x_matrix_put <- as.matrix(putamen_df[, genes_list])
mode(x_matrix_put) <- "numeric"

rf_put <- randomForest(x = x_matrix_put, y = putamen_df$Group, importance = TRUE, ntree = 1000)

imp_put <- as.data.frame(importance(rf_put)) %>% 
  mutate(
    Gene = rownames(.), 
    Region = "Putamen", 
    Importance = MeanDecreaseGini, 
    # 纠正点：ACP5 和 RAMP3 为 Up-regulated
    Regulation = if_else(Gene %in% c("ACP5", "RAMP3"), "Up-regulated (in AD/PD)", "Down-regulated (in AD/PD)")
  )

# ------------------------------------------------------------------------------
# 3. 绘制 Panel B 特征权重条形图 (上调红色=#DC0000FF，下调蓝色=#3C5488FF)
# ------------------------------------------------------------------------------
importance_pd_all <- bind_rows(imp_cau, imp_put)

p_pd_importance <- ggplot(importance_pd_all, aes(x = reorder(Gene, Importance), y = Importance, fill = Regulation)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.9, color = "black", linewidth = 0.3) +
  coord_flip() +
  facet_wrap(~Region, scales = "free_x") +
  # 映射设置：上调为红色，下调为蓝色
  scale_fill_manual(values = c("Up-regulated (in AD/PD)" = "#DC0000FF", "Down-regulated (in AD/PD)" = "#3C5488FF")) +
  labs(
    title = "Panel B: Random Forest Classifier Weights\n(Parkinson's Disease: GSE205450)",
    x = "Shared Target Genes",
    y = "Mean Decrease Gini (Importance Score)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray95"),
    plot.title = element_text(face = "bold", size = 10),
    axis.title = element_text(face = "bold", size = 10),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 8),
    strip.background = element_rect(fill = "gray93", color = "black"),
    strip.text = element_text(face = "bold", size = 9),
    aspect.ratio = 0.85
  )

# ------------------------------------------------------------------------------
# 4. 导出 PDF 矢量图
# ------------------------------------------------------------------------------
pdf_output_name <- file.path(base_dir, "Figure_PD_GSE205450_RF_Importance.pdf")

ggsave(
  filename = pdf_output_name,
  plot = p_pd_importance,
  width = 6.8,
  height = 4.5,
  device = "pdf"
)

cat(paste0("\n✨ 修正后的特征权重图导出成功：'", pdf_output_name, "'\n"))