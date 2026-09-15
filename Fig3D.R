# ==============================================================================
# 🎯 [AD 溶酶体 Signature Score 差异分析·相对路径+科学计数法 P 值版]
# ==============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------------------------
# 0. 相对路径与数据载入
# ------------------------------------------------------------------------------
# 使用相对路径从当前工作目录加载文件
input_csv <- file.path(".", "AD_Lysosome_Subtypes_With_Expression_and_PCA.csv")

if (!file.exists(input_csv)) {
  stop(sprintf("⚠️ 错误：未在相对路径下找到数据文件：\n%s", input_csv))
}

cat(sprintf("📊 1. 正在从相对路径读取数据文件:\n %s\n", input_csv))
df_raw <- read.csv(input_csv, stringsAsFactors = FALSE, header = TRUE)

# ------------------------------------------------------------------------------
# 1. 计算 Signature A 与 Signature B 得分
# ------------------------------------------------------------------------------
# Signature A (CLTC & SLC7A14)
# Signature B (SLC2A13 & ACP5 & RAMP3)
required_genes <- c("CLTC", "SLC7A14", "SLC2A13", "ACP5", "RAMP3")
missing_genes  <- setdiff(required_genes, colnames(df_raw))

if (length(missing_genes) > 0) {
  stop(sprintf("⚠️ 错误：数据集中缺失以下核心基因：%s", paste(missing_genes, collapse = ", ")))
}

ad_score_data <- df_raw %>%
  mutate(
    Signature_Score_A = (CLTC + SLC7A14) / 2,
    Signature_Score_B = (SLC2A13 + ACP5 + RAMP3) / 3,
    Lysosome_Subtype  = factor(Lysosome_Subtype, levels = c("Cluster_1", "Cluster_2"))
  )

# ------------------------------------------------------------------------------
# 2. P 值 plotmath 格式化函数（含防 0/Inf 及科学计数法）
# ------------------------------------------------------------------------------
format_p_plotmath <- function(p) {
  if (is.na(p)) {
    return("italic(P) == 'NA'")
  }
  
  if (p == 0 || p < 2e-16) {
    return("italic(P) < 2 %*% 10^-16")
  } else if (p < 0.001) {
    p_sig <- signif(p, 2)
    exp_val <- floor(log10(p_sig))
    base_val <- p_sig / (10^exp_val)
    base_str <- sprintf("%.2g", base_val)
    
    return(sprintf("italic(P) == %s %%*%% 10^%d", base_str, exp_val))
  } else {
    p_formatted <- formatC(signif(p, 2), digits = 2, format = "fg", flag = "#")
    return(sprintf("italic(P) == '%s'", p_formatted))
  }
}

# ------------------------------------------------------------------------------
# 🎨 3. 绘制 Signature A 得分对比图
# ------------------------------------------------------------------------------
p_val_A      <- t.test(Signature_Score_A ~ Lysosome_Subtype, data = ad_score_data)$p.value
p_expr_str_A <- format_p_plotmath(p_val_A)

max_y_A   <- max(ad_score_data$Signature_Score_A, na.rm = TRUE)
min_y_A   <- min(ad_score_data$Signature_Score_A, na.rm = TRUE)
y_range_A <- max_y_A - min_y_A

p_A <- ggplot(ad_score_data, aes(x = Lysosome_Subtype, y = Signature_Score_A)) +
  geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf), fill = "#E5E7EB", alpha = 0.3) +
  geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf), fill = "#FECACA", alpha = 0.3) +
  geom_jitter(width = 0.15, size = 1.3, alpha = 0.7, color = "black") +
  geom_boxplot(width = 0.3, fill = NA, color = "black", outlier.shape = NA, linewidth = 0.8) +
  annotate("text", x = 1.5, y = max_y_A + y_range_A * 0.1, label = p_expr_str_A, parse = TRUE, size = 5) +
  labs(title = "Lysosome Subtype Signature A\n(CLTC & SLC7A14 enriched)", x = "", y = "Signature Score") +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.18))) +
  theme_classic() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text    = element_text(size = 13, color = "black", face = "bold"),
    axis.title.y = element_text(size = 13, face = "plain"),
    axis.line    = element_line(linewidth = 0.8),
    axis.ticks   = element_line(linewidth = 0.8)
  )

pdf_path_A <- file.path(".", "AD_Subtype_Signature_Score_A.pdf")
ggsave(pdf_path_A, p_A, width = 4, height = 4.8)

# ------------------------------------------------------------------------------
# 🎨 4. 绘制 Signature B 得分对比图
# ------------------------------------------------------------------------------
p_val_B      <- t.test(Signature_Score_B ~ Lysosome_Subtype, data = ad_score_data)$p.value
p_expr_str_B <- format_p_plotmath(p_val_B)

max_y_B   <- max(ad_score_data$Signature_Score_B, na.rm = TRUE)
min_y_B   <- min(ad_score_data$Signature_Score_B, na.rm = TRUE)
y_range_B <- max_y_B - min_y_B

p_B <- ggplot(ad_score_data, aes(x = Lysosome_Subtype, y = Signature_Score_B)) +
  geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf), fill = "#E5E7EB", alpha = 0.3) +
  geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf), fill = "#FECACA", alpha = 0.3) +
  geom_jitter(width = 0.15, size = 1.3, alpha = 0.7, color = "black") +
  geom_boxplot(width = 0.3, fill = NA, color = "black", outlier.shape = NA, linewidth = 0.8) +
  annotate("text", x = 1.5, y = max_y_B + y_range_B * 0.1, label = p_expr_str_B, parse = TRUE, size = 5) +
  labs(title = "Lysosome Subtype Signature B\n(SLC2A13 & ACP5 & RAMP3 enriched)", x = "", y = "Signature Score") +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.18))) +
  theme_classic() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text    = element_text(size = 13, color = "black", face = "bold"),
    axis.title.y = element_text(size = 13, face = "plain"),
    axis.line    = element_line(linewidth = 0.8),
    axis.ticks   = element_line(linewidth = 0.8)
  )

pdf_path_B <- file.path(".", "AD_Subtype_Signature_Score_B.pdf")
ggsave(pdf_path_B, p_B, width = 4, height = 4.8)

# ------------------------------------------------------------------------------
# 输出完成提示信息
# ------------------------------------------------------------------------------
cat("\n====================================================================\n")
cat("🎉 【计算与相对路径出图成功！】：\n")
cat("--------------------------------------------------------------------\n")
cat(sprintf(" 1. Signature A 图形: %s (Raw p = %.4e)\n", pdf_path_A, p_val_A))
cat(sprintf(" 2. Signature B 图形: %s (Raw p = %.4e)\n", pdf_path_B, p_val_B))
cat("====================================================================\n")