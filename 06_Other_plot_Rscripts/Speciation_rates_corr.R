library(ggplot2)
library(ggpubr)
setwd("")
df1 <- read.csv("Bamm_locus_bryophyte_rates.csv", header = TRUE)
df2 <- read.csv("Bamm_nuclear_bryophyte_rates.csv" ,header = TRUE)

combined_df <- inner_join(df1, df2, by = "species", suffix = c("_Bamm_locus", "_Bamm_nuclear"))
print(combined_df)

# 计算相关系数和P值
cor_test <- cor.test(combined_df$net_Bamm_locus, combined_df$net_Bamm_nuclear, method = "pearson")
r_value <- round(cor_test$estimate, 2)
p_value <- ifelse(cor_test$p.value < 0.001, "< 0.001", 
                  sprintf("%.3f", cor_test$p.value))
n <- nrow(combined_df)  # 明确定义样本量
# 创建散点图
p <- ggplot(combined_df, aes(x = net_Bamm_nuclear, y = net_Bamm_locus)) +
  geom_point(color = "grey", size = 0.3, alpha = 1, shape = 19) +
  #geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black", size = 0.15) + # 添加1:1参考线
  labs(x = "Bamm in nuclear tree",
       y = "Bamm in locus tree",
       title = "Net speciation rates correlation") +
  theme_minimal(base_size = 7) +
  theme(
    # 去掉网格线
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 添加轴刻度线
    axis.ticks = element_line(color = "black", size = 0.2),
    # 添加矩形框
    panel.border = element_rect(color = "black", fill = NA, size = 0.2),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, 0.01)) +
  scale_y_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, 0.01)) +
  annotate("text", x = 0.02, y = 0.045, label = paste0("n = ", n, "; r = ", r_value, "; P ", p_value), size = 2, fontface = "bold")
# 显示图形
print(p)
# 保存图形（可选）
ggsave("Bamm_rate_correlation.pdf", plot = p, width = 2, height = 2, dpi = 300)

############################################################################################

df1 <- read.csv("Bamm_locus_bryophyte_rates.csv", header = TRUE)
df2 <- read.csv("Revbayes_locus_bryophyte_rates.csv" ,header = TRUE)

combined_df <- inner_join(df1, df2, by = "species", suffix = c("_Bamm_locus", "_Revbayes_locus"))
print(combined_df)

# 计算相关系数和P值
cor_test <- cor.test(combined_df$net_Bamm_locus, combined_df$net_Revbayes_locus, method = "pearson")
r_value <- round(cor_test$estimate, 2)
p_value <- ifelse(cor_test$p.value < 0.001, "< 0.001", 
                  sprintf("%.3f", cor_test$p.value))
n <- nrow(combined_df)  # 明确定义样本量
# 创建散点图
p <- ggplot(combined_df, aes(x = net_Bamm_locus, y = net_Revbayes_locus)) +
  geom_point(color = "grey", size = 0.3, alpha = 1, shape = 19) +
  #geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black", size = 0.15) + # 添加1:1参考线
  labs(x = "Bamm in locus tree",
       y = "Revbayes in locus tree",
       title = "Net speciation rates correlation") +
  theme_minimal(base_size = 7) +
  theme(
    # 去掉网格线
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 添加轴刻度线
    axis.ticks = element_line(color = "black", size = 0.2),
    # 添加矩形框
    panel.border = element_rect(color = "black", fill = NA, size = 0.2),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, 0.01)) +
  scale_y_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, 0.01)) +
  annotate("text", x = 0.02, y = 0.045, label = paste0("n = ", n, "; r = ", r_value, "; P ", p_value), size = 2, fontface = "bold")
# 显示图形
print(p)
# 保存图形（可选）
ggsave("Bamm_Revbayes_locus_rate_correlation.pdf", plot = p, width = 2, height = 2, dpi = 300)

############################################################################################

df1 <- read.csv("Bamm_locus_bryophyte_rates.csv", header = TRUE)
df2 <- read.csv("Revbayes_nuclear_bryophyte_rates.csv" ,header = TRUE)

combined_df <- inner_join(df1, df2, by = "species", suffix = c("_Bamm_locus", "_Revbayes_nuclear"))
print(combined_df)

# 计算相关系数和P值
cor_test <- cor.test(combined_df$net_Bamm_locus, combined_df$net_Revbayes_nuclear, method = "pearson")
r_value <- round(cor_test$estimate, 2)
p_value <- ifelse(cor_test$p.value < 0.001, "< 0.001", 
                  sprintf("%.3f", cor_test$p.value))
n <- nrow(combined_df)  # 明确定义样本量
# 创建散点图
p <- ggplot(combined_df, aes(x = net_Bamm_locus, y = net_Revbayes_nuclear)) +
  geom_point(color = "grey", size = 0.3, alpha = 1, shape = 19) +
  #geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black", size = 0.15) + # 添加1:1参考线
  labs(x = "Bamm in locus tree",
       y = "Revbayes in nuclear tree",
       title = "Net speciation rates correlation") +
  theme_minimal(base_size = 7) +
  theme(
    # 去掉网格线
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 添加轴刻度线
    axis.ticks = element_line(color = "black", size = 0.2),
    # 添加矩形框
    panel.border = element_rect(color = "black", fill = NA, size = 0.2),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, 0.01)) +
  scale_y_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, 0.01)) +
  annotate("text", x = 0.02, y = 0.045, label = paste0("n = ", n, "; r = ", r_value, "; P ", p_value), size = 2, fontface = "bold")
# 显示图形
print(p)
# 保存图形（可选）
ggsave("Bamm_Revbayes_locus_nuclear_rate_correlation.pdf", plot = p, width = 2, height = 2, dpi = 300)

############################################################################################

df1 <- read.csv("Bamm_nuclear_bryophyte_rates.csv", header = TRUE)
df2 <- read.csv("Revbayes_nuclear_bryophyte_rates.csv" ,header = TRUE)

combined_df <- inner_join(df1, df2, by = "species", suffix = c("_Bamm_nuclear", "_Revbayes_nuclear"))
print(combined_df)

# 计算相关系数和P值
cor_test <- cor.test(combined_df$net_Bamm_nuclear, combined_df$net_Revbayes_nuclear, method = "pearson")
r_value <- round(cor_test$estimate, 2)
p_value <- ifelse(cor_test$p.value < 0.001, "< 0.001", 
                  sprintf("%.3f", cor_test$p.value))
n <- nrow(combined_df)  # 明确定义样本量
# 创建散点图
p <- ggplot(combined_df, aes(x = net_Bamm_nuclear, y = net_Revbayes_nuclear)) +
  geom_point(color = "grey", size = 0.3, alpha = 1, shape = 19) +
  #geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black", size = 0.15) + # 添加1:1参考线
  labs(x = "Bamm in nuclear tree",
       y = "Revbayes in nuclear tree",
       title = "Net speciation rates correlation") +
  theme_minimal(base_size = 7) +
  theme(
    # 去掉网格线
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 添加轴刻度线
    axis.ticks = element_line(color = "black", size = 0.2),
    # 添加矩形框
    panel.border = element_rect(color = "black", fill = NA, size = 0.2),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, 0.01)) +
  scale_y_continuous(limits = c(0, 0.05), breaks = seq(0, 0.05, 0.01)) +
  annotate("text", x = 0.02, y = 0.045, label = paste0("n = ", n, "; r = ", r_value, "; P ", p_value), size = 2, fontface = "bold")
# 显示图形
print(p)
# 保存图形（可选）
ggsave("Bamm_Revbayes_nuclear_rate_correlation.pdf", plot = p, width = 2, height = 2, dpi = 300)

#########################箱图绘制##########################################
library(ggplot2)
library(dplyr)
library(patchwork)
# ============================
# 1. 预留颜色参数（可自行修改）
# ============================
col_genus <- "#56B4E9"   # 自己改颜色
point_col <- "black"

# ---- 数据整合函数 ----
load_pair <- function(nuclear_file, locus_file, group_name) {
  df1 <- read.csv(nuclear_file, header = TRUE)
  df2 <- read.csv(locus_file, header = TRUE)
  
  # --- 检查列名 ---
  required_cols <- c("species", "net")
  
  if (!all(required_cols %in% colnames(df1))) {
    stop(paste("❌", nuclear_file, "必须包含: species, net"))
  }
  if (!all(required_cols %in% colnames(df2))) {
    stop(paste("❌", locus_file, "必须包含: species, net"))
  }
  
  # --- 合并 ---
  merged <- inner_join(df1, df2, by = "species",
                       suffix = c("_nuclear", "_locus"))
  
  if (nrow(merged) == 0) {
    stop(paste0("❌ ", group_name, " 没有共同 species，请检查文件"))
  }
  
  # --- 计算差值 ---
  merged$age_diff <- merged$net_nuclear - merged$net_locus
  merged$Group <- group_name
  return(merged)
}

# ============================
# 读取数据（只有 Genus）
# ============================
df_genus <- load_pair(
  "Bamm_nuclear_bryophyte_rates.csv",
  "Bamm_locus_bryophyte_rates.csv",
  "Genus"
)

df_all <- df_genus

# --- 设定 Group 因子顺序 ----
df_all$Group <- factor(df_all$Group, levels = "Genus")

# ============================
# 4. 绘图
# ============================
p <- ggplot(df_all, aes(x = Group, y = age_diff, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.85, color = NA) +
  geom_jitter(width = 0.12, size = 1.6, color = point_col, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#555555") +
  scale_fill_manual(values = c("Genus" = col_genus)) +
  labs(y = "Rate difference", x = "") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  )

print(p)
# 输出
ggsave("Rate_difference_two_bamm.pdf", plot = p, width = 8, height = 6, dpi = 300)
