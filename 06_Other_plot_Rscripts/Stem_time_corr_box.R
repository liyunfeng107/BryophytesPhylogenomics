library(ggplot2)
library(ggpubr)
setwd("")

#########################################################################################
df1 <- read.csv("fam_nuclear.csv", header = TRUE)
df2 <- read.csv("fam_Bechteler.csv" ,header = TRUE)
combined_df <- inner_join(df1, df2, by = "label", suffix = c("_thisstudy", "_Bechteler"))
print(combined_df)
# 计算相关系数和P值
cor_test <- cor.test(combined_df$stemage_thisstudy, combined_df$stemage_Bechteler, method = "pearson")
r_value <- round(cor_test$estimate, 2)
p_value <- ifelse(cor_test$p.value < 0.001, "< 0.001", 
                  sprintf("%.3f", cor_test$p.value))
n <- nrow(combined_df)  # 明确定义样本量
p <- ggplot(combined_df, aes(x = stemage_thisstudy, y = stemage_Bechteler)) +
  geom_point(color = "grey", size = 2, alpha = 1) +
  #geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black") + # 添加1:1参考线
  labs(x = "this study(Ma)",
       y = "Bechteler tree(Ma)",
       title = "Stem Group Time Correlation") +
  theme_minimal(base_size = 14) +
  theme(
    # 去掉网格线
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 添加轴刻度线
    axis.ticks = element_line(color = "black"),
    # 添加矩形框
    panel.border = element_rect(color = "black", fill = NA),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(0, 450), breaks = seq(0, 450, 100)) +
  scale_y_continuous(limits = c(0, 450), breaks = seq(0, 450, 100)) +
  annotate("text", x = 100, y = 400, label = paste0("n = ", n, "; r = ", r_value, "; P ", p_value), size = 5, fontface = "bold")
# 显示图形
print(p)
ggsave("fam_stem_group_Bechteler_correlation.pdf", plot = p, width = 8, height = 6, dpi = 300)

#########################################################################################
df1 <- read.csv("order_nuclear.csv", header = TRUE)
df2 <- read.csv("order_Bechteler.csv" ,header = TRUE)
combined_df <- inner_join(df1, df2, by = "label", suffix = c("_thisstudy", "_Bechteler"))
print(combined_df)
# 计算相关系数和P值
cor_test <- cor.test(combined_df$stemage_thisstudy, combined_df$stemage_Bechteler, method = "pearson")
r_value <- round(cor_test$estimate, 2)
p_value <- ifelse(cor_test$p.value < 0.001, "< 0.001", 
                  sprintf("%.3f", cor_test$p.value))
n <- nrow(combined_df)  # 明确定义样本量
p <- ggplot(combined_df, aes(x = stemage_thisstudy, y = stemage_Bechteler)) +
  geom_point(color = "grey", size = 2, alpha = 1) +
  #geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black") + # 添加1:1参考线
  labs(x = "this study(Ma)",
       y = "Bechteler tree(Ma)",
       title = "Stem Group Time Correlation") +
  theme_minimal(base_size = 14) +
  theme(
    # 去掉网格线
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 添加轴刻度线
    axis.ticks = element_line(color = "black"),
    # 添加矩形框
    panel.border = element_rect(color = "black", fill = NA),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(0, 450), breaks = seq(0, 450, 100)) +
  scale_y_continuous(limits = c(0, 450), breaks = seq(0, 450, 100)) +
  annotate("text", x = 100, y = 400, label = paste0("n = ", n, "; r = ", r_value, "; P ", p_value), size = 5, fontface = "bold")
# 显示图形
print(p)
ggsave("order_stem_group_Bechteler_correlation.pdf", plot = p, width = 8, height = 6, dpi = 300)

#########################箱图绘制##########################################
library(ggplot2)
library(dplyr)
library(patchwork)
# ============================
# 1. 预留颜色参数（可自行修改）
# ============================
col_order  <- "#42dec8"   # 橙色
col_family <- "#55b4e9"   # 蓝色
point_col  <- "#333333"     # 散点颜色

# ============================
# 2. 定义一个读取+合并+计算差值的函数
# ============================
load_pair <- function(nuclear_file, locus_file, group_name) {
  df1 <- read.csv(nuclear_file, header = TRUE)
  df2 <- read.csv(locus_file, header = TRUE)
  # 检查 label
  if (!all(c("label", "stemage") %in% colnames(df1))) {
    stop(paste("❌", nuclear_file, "不包含 label, stemage 两列"))
  }
  if (!all(c("label", "stemage") %in% colnames(df2))) {
    stop(paste("❌", locus_file, "不包含 label, stemage 两列"))
  }
  # 合并
  merged <- inner_join(df1, df2, by = "label",
                       suffix = c("_nuclear", "_locus"))
  if (nrow(merged) == 0) {
    stop(paste0(
      "❌ ", group_name, " 没有共同 label。\n",
      "请检查两份文件的 label 是否一致。"
    ))
  }
  # 正确的列名 → stemage_nuclear / stemage_locus
  merged$age_diff <- merged$stemage_nuclear - merged$stemage_locus
  merged$Group <- group_name
  return(merged)
}

# ============================
# 3. 读取三类数据
# ============================
df_order  <- load_pair("order_nuclear.csv",  "order_Bechteler.csv",  "Order")
df_family <- load_pair("fam_nuclear.csv",    "fam_Bechteler.csv",    "Family")

# 合并一起，便于统一绘图
df_all <- rbind(df_order, df_family)

# 指定顺序：目 → 科
df_all$Group <- factor(df_all$Group, levels = c("Order", "Family"))

# ============================
# 4. 绘图
# ============================
p <- ggplot(df_all, aes(x = Group, y = age_diff, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.85, color = NA) +
  geom_jitter(width = 0.12, size = 1, color = point_col, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c(
    "Order"  = col_order,
    "Family" = col_family
  )) +
  labs(y = "Age difference (Ma)", x = "") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  )

print(p)
# 输出
ggsave("age_difference_Bechteler.pdf", plot = p, width = 8, height = 6, dpi = 300)


#################################genus###############################################

df1 <- read.csv("genus_HPD_20250604.csv", header = TRUE)
df2 <- read.csv("genus_HPD_locus_20250604.csv" ,header = TRUE)

combined_df <- inner_join(df1, df2, by = "label", suffix = c("_nuclear", "_marker"))

print(combined_df)

# 计算相关系数和P值
cor_test <- cor.test(combined_df$stemage_nuclear, combined_df$stemage_marker, method = "pearson")
r_value <- round(cor_test$estimate, 2)
p_value <- ifelse(cor_test$p.value < 0.001, "< 0.001", 
                  sprintf("%.3f", cor_test$p.value))
n <- nrow(combined_df)  # 明确定义样本量
# 创建散点图
p <- ggplot(combined_df, aes(x = stemage_nuclear, y = stemage_marker)) +
  geom_point(color = "grey", size = 2, alpha = 1) +
  #geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black") + # 添加1:1参考线
  labs(x = "nuclear tree(Ma)",
       y = "locus tree(Ma)",
       title = "Stem Group Time Correlation") +
  theme_minimal(base_size = 14) +
  theme(
    # 去掉网格线
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 添加轴刻度线
    axis.ticks = element_line(color = "black"),
    # 添加矩形框
    panel.border = element_rect(color = "black", fill = NA),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(0, 450), breaks = seq(0, 450, 100)) +
  scale_y_continuous(limits = c(0, 450), breaks = seq(0, 450, 100)) +
  annotate("text", x = 100, y = 400, label = paste0("n = ", n, "; r = ", r_value, "; P ", p_value), size = 5, fontface = "bold")
# 显示图形
print(p)


# 保存图形（可选）
ggsave("genus_stem_group_correlation.pdf", plot = p, width = 8, height = 6, dpi = 300)

#################################fam###############################################
df1 <- read.csv("fam_HPD_20250604.csv", header = TRUE)
df2 <- read.csv("fam_HPD_locus_20250604.csv" ,header = TRUE)

combined_df <- inner_join(df1, df2, by = "label", suffix = c("_nuclear", "_marker"))
print(combined_df)

# 计算相关系数和P值
cor_test <- cor.test(combined_df$stemage_nuclear, combined_df$stemage_marker, method = "pearson")
r_value <- round(cor_test$estimate, 2)
p_value <- ifelse(cor_test$p.value < 0.001, "< 0.001", 
                  sprintf("%.3f", cor_test$p.value))
n <- nrow(combined_df)  # 明确定义样本量
# 创建散点图
p <- ggplot(combined_df, aes(x = stemage_nuclear, y = stemage_marker)) +
  geom_point(color = "grey", size = 2, alpha = 1) +
  #geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black") + # 添加1:1参考线
  labs(x = "nuclear tree(Ma)",
       y = "locus tree(Ma)",
       title = "Stem Group Time Correlation") +
  theme_minimal(base_size = 14) +
  theme(
    # 去掉网格线
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 添加轴刻度线
    axis.ticks = element_line(color = "black"),
    # 添加矩形框
    panel.border = element_rect(color = "black", fill = NA),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(0, 450), breaks = seq(0, 450, 100)) +
  scale_y_continuous(limits = c(0, 450), breaks = seq(0, 450, 100)) +
  annotate("text", x = 100, y = 400, label = paste0("n = ", n, "; r = ", r_value, "; P ", p_value), size = 5, fontface = "bold")
# 显示图形
print(p)
# 保存图形（可选）
ggsave("fam_stem_group_correlation.pdf", plot = p, width = 8, height = 6, dpi = 300)

#################################order###############################################
df1 <- read.csv("order_HPD_20250604.csv", header = TRUE)
df2 <- read.csv("order_HPD_locus_20250604.csv" ,header = TRUE)

combined_df <- inner_join(df1, df2, by = "label", suffix = c("_nuclear", "_marker"))

print(combined_df)

# 计算相关系数和P值
cor_test <- cor.test(combined_df$stemage_nuclear, combined_df$stemage_marker, method = "pearson")
r_value <- round(cor_test$estimate, 2)
p_value <- ifelse(cor_test$p.value < 0.001, "< 0.001", 
                  sprintf("%.3f", cor_test$p.value))
n <- nrow(combined_df)  # 明确定义样本量
# 创建散点图
p <- ggplot(combined_df, aes(x = stemage_nuclear, y = stemage_marker)) +
  geom_point(color = "grey", size = 2, alpha = 1) +
  #geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black") + # 添加1:1参考线
  labs(x = "nuclear tree(Ma)",
       y = "locus tree(Ma)",
       title = "Stem Group Time Correlation") +
  theme_minimal(base_size = 14) +
  theme(
    # 去掉网格线
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 添加轴刻度线
    axis.ticks = element_line(color = "black"),
    # 添加矩形框
    panel.border = element_rect(color = "black", fill = NA),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(0, 450), breaks = seq(0, 450, 100)) +
  scale_y_continuous(limits = c(0, 450), breaks = seq(0, 450, 100)) +
  annotate("text", x = 100, y = 400, label = paste0("n = ", n, "; r = ", r_value, "; P ", p_value), size = 5, fontface = "bold")
# 显示图形
print(p)


# 保存图形（可选）
ggsave("order_stem_group_correlation.pdf", plot = p, width = 8, height = 6, dpi = 300)

#########################箱图绘制##########################################
library(ggplot2)
library(dplyr)
library(patchwork)
# ============================
# 1. 预留颜色参数（可自行修改）
# ============================
col_order  <- "#42dec8"   # 橙色
col_family <- "#55b4e9"   # 蓝色
col_genus  <- "#fba898"   # 绿色
point_col  <- "#333333"     # 散点颜色

# ============================
# 2. 定义一个读取+合并+计算差值的函数
# ============================
load_pair <- function(nuclear_file, locus_file, group_name) {
  df1 <- read.csv(nuclear_file, header = TRUE)
  df2 <- read.csv(locus_file, header = TRUE)
  # 检查 label
  if (!all(c("label", "stemage") %in% colnames(df1))) {
    stop(paste("❌", nuclear_file, "不包含 label, stemage 两列"))
  }
  if (!all(c("label", "stemage") %in% colnames(df2))) {
    stop(paste("❌", locus_file, "不包含 label, stemage 两列"))
  }
  # 合并
  merged <- inner_join(df1, df2, by = "label",
                       suffix = c("_nuclear", "_locus"))
  if (nrow(merged) == 0) {
    stop(paste0(
      "❌ ", group_name, " 没有共同 label。\n",
      "请检查两份文件的 label 是否一致。"
    ))
  }
  # 正确的列名 → stemage_nuclear / stemage_locus
  merged$age_diff <- merged$stemage_nuclear - merged$stemage_locus
  merged$Group <- group_name
  return(merged)
}

# ============================
# 3. 读取三类数据
# ============================
df_order  <- load_pair("order_nuclear.csv",  "order_locus.csv",  "Order")
df_family <- load_pair("fam_nuclear.csv",    "fam_locus.csv",    "Family")
df_genus  <- load_pair("genus_nuclear.csv",  "genus_locus.csv",  "Genus")

# 合并一起，便于统一绘图
df_all <- rbind(df_order, df_family, df_genus)

# 指定顺序：目 → 科 → 属
df_all$Group <- factor(df_all$Group, levels = c("Order", "Family", "Genus"))

# ============================
# 4. 绘图
# ============================
p <- ggplot(df_all, aes(x = Group, y = age_diff, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.85, color = NA) +
  geom_jitter(width = 0.12, size = 1, color = point_col, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c(
    "Order"  = col_order,
    "Family" = col_family,
    "Genus"  = col_genus
  )) +
  labs(y = "Age difference (Ma)", x = "") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  )

print(p)
# 输出
ggsave("age_difference_three_groups.pdf", plot = p, width = 8, height = 6, dpi = 300)
