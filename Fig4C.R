# ==============================================================================
# 🎯 [PD 溶酶体核心基因亚型差异表达分析·相对路径+精准图示 P 值科学计数法版]
# ==============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------------------------
# 0. 相对路径与数据载入
# ------------------------------------------------------------------------------
# 采用相对路径，直接从当前工作目录读取文件
input_csv <- file.path(".", "PD_Lysosome_Subtypes_With_Expression_and_PCA.csv")

if (!file.exists(input_csv)) {
  stop(sprintf("⚠️ 错误：未在相对路径下找到 PD 数据文件：\n%s", input_csv))
}

cat(sprintf("📊 1. 正在从相对路径读取 PD 数据文件:\n %s\n", input_csv))
df_raw <- read.csv(input_csv, stringsAsFactors = FALSE, header = TRUE)

# ------------------------------------------------------------------------------
# 1. 数据重构与准备
# ------------------------------------------------------------------------------
target_genes    <- c("RAMP3", "ACP5", "SLC7A14", "CLTC", "SLC2A13")
available_genes <- intersect(target_genes, colnames(df_raw))

if (length(available_genes) == 0) {
  stop("⚠️ 错误：数据集中未找到指定的溶酶体核心基因！")
}

# 转换为长数据格式
pd_plot_data <- df_raw %>%
  select(all_of(c("Lysosome_Subtype", available_genes))) %>%
  pivot_longer(
    cols      = all_of(available_genes),
    names_to  = "Gene",
    values_to = "Expression"
  ) %>%
  filter(!is.na(Expression)) %>%
  mutate(
    Lysosome_Subtype = factor(Lysosome_Subtype, levels = c("Cluster_1", "Cluster_2"))
  )

# ------------------------------------------------------------------------------
# 2. P 值 plotmath 格式化函数（含防 0/Inf 及科学计数法）
# ------------------------------------------------------------------------------
format_p_plotmath <- function(p) {
  if (is.na(p)) {
    return("italic(P) == 'NA'")
  }
  
  # 针对 p == 0 或小于 2e-16 的极端情况
  if (p == 0 || p < 2e-16) {
    return("italic(P) < 2 %*% 10^-16")
  } else if (p < 0.001) {
    # 科学计数法格式化 (保留 2 位有效数字)
    p_sig <- signif(p, 2)
    exp_val <- floor(log10(p_sig))
    base_val <- p_sig / (10^exp_val)
    
    # 清理基数显示 (例如 3.0 转为 3)
    base_str <- sprintf("%.2g", base_val)
    
    return(sprintf("italic(P) == %s %%*%% 10^%d", base_str, exp_val))
  } else {
    # 常规数值 (保留 2 位有效数字)
    p_formatted <- formatC(signif(p, 2), digits = 2, format = "fg", flag = "#")
    return(sprintf("italic(P) == '%s'", p_formatted))
  }
}

# ------------------------------------------------------------------------------
# 3. 循环绘制并导出 PDF（相对路径输出）
# ------------------------------------------------------------------------------
cat("\n====================================================================\n")
cat("📊 【正在绘制 PD 两个亚型间的溶酶体基因差异表达图...】\n")
cat("--------------------------------------------------------------------\n")

for (g in available_genes) {
  d_sub <- pd_plot_data %>% filter(Gene == g)
  
  # 计算 Cluster 1 vs Cluster 2 的 p 值
  p_val      <- t.test(Expression ~ Lysosome_Subtype, data = d_sub)$p.value
  p_expr_str <- format_p_plotmath(p_val)
  
  max_y   <- max(d_sub$Expression, na.rm = TRUE)
  min_y   <- min(d_sub$Expression, na.rm = TRUE)
  y_range <- max_y - min_y
  
  # 严格复刻经典极简条带 Prism 风格
  p <- ggplot(d_sub, aes(x = Lysosome_Subtype, y = Expression)) +
    # 背景区域颜色块 (Cluster 1 浅灰，Cluster 2 浅粉)
    geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf), fill = "#E5E7EB", alpha = 0.3) +
    geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf), fill = "#FECACA", alpha = 0.3) +
    # 纯黑抖动散点
    geom_jitter(width = 0.15, size = 1.3, alpha = 0.7, color = "black") +
    # 覆盖其上的极简箱线图
    geom_boxplot(width = 0.3, fill = NA, color = "black", outlier.shape = NA, linewidth = 0.8) +
    # 显著性标签
    annotate("text", x = 1.5, y = max_y + y_range * 0.1, label = p_expr_str, parse = TRUE, size = 5) +
    labs(title = g, x = "", y = "Normalized Expression Value") +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.18))) +
    theme_classic() +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 16, face = "bold.italic"),
      axis.text    = element_text(size = 13, color = "black", face = "bold"),
      axis.title.y = element_text(size = 13, face = "plain"),
      axis.line    = element_line(linewidth = 0.8),
      axis.ticks   = element_line(linewidth = 0.8)
    )
  
  # 导出 PDF 至相对路径 (标注 PD 前缀)
  pdf_path <- file.path(".", paste0("PD_Subtype_Diff_", g, ".pdf"))
  ggsave(pdf_path, p, width = 4, height = 4.8)
  
  cat(sprintf("✅ 已生成 PD 亚型差异图至相对路径: %s (Raw p = %.4e)\n", pdf_path, p_val))
}

cat("\n🎉 【PD 数据集全套差异表达图生成大功告成！】\n")