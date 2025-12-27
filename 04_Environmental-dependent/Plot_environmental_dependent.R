setwd("")
library(RevGadgets)
library(ggplot2)

CO <- c(1, 2.6, 2.6, 2.6, 4.9, 7, 6.5, 8, 7, 8, 9, 9, 9, 9, 10, 7, 6, 7, 10, 10, 17, 20, 20, 24, 28, 25, 10, 6.5, 3, 2.3, 3.5, 2.8, 6, 7, 8, 8, 8, 11, 15, 18, 9, 8, 6, 7, 6, 6, 16, 23, 24, 25, 25, 23, 23, 27, 28)

MAX_VAR_AGE = 540
NUM_INTERVALS = 10
CO_age <- seq(from = 0, to = MAX_VAR_AGE, by = NUM_INTERVALS)
print(Pre_age)

env_data <- data.frame(
  age = rev(CO_age),  # 反转时间轴
  CO2 = rev(CO)        # 反转CO2值
)

env_plot <- ggplot(env_data, aes(x = age, y = CO2)) +
  geom_line(color = "black", linewidth = 0.2) +  # 蓝色曲线
  #geom_point(color = "darkblue", size = 3) +     # 深蓝色点
  scale_x_reverse("Age (Ma)",                   # X轴反转（从大到小）
                  limits = c(max(env_data$age), min(env_data$age)),
                  breaks = seq(0, 540, by = 100)) +  # 设置刻度
  scale_y_continuous("CO2 (ppm)", 
                     #position = "right",
                     limits = c(0, 30)) +      # 设置Y轴范围
  ggtitle("CO2 Concentration Over Time") +      # 标题
  theme_minimal(base_size = 6) +
  theme_void() +                                  # 黑白主题
  theme(
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    plot.title = element_text(hjust = 0.5, size = 6, face = "bold"),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6)
  )

# 打印图形
print(env_plot)

# 保存图形
ggsave("CO2_curve.pdf", env_plot, width = 80, height = 80, units = "mm")
rates <- processDivRates(
  speciation_time_log = "EBD_Corr_speciation_times.log",
  speciation_rate_log = "EBD_Corr_speciation_rates.log",
  extinction_time_log = "EBD_Corr_extinction_times.log",
  extinction_rate_log = "EBD_Corr_extinction_rates.log",
  burnin=0.25)
p <- plotDivRates(rates)
print(p)
ggsave("EBD_Corr_liver_Co.pdf", p, width = 160, height = 160, units = "mm")

rates <- processDivRates(
  speciation_time_log = "EBD_Corr_speciation_times.log",
  speciation_rate_log = "EBD_Corr_speciation_rates.log",
  extinction_time_log = "EBD_Corr_extinction_times.log",
  extinction_rate_log = "EBD_Corr_extinction_rates.log",
  burnin=0.25)
p <- plotDivRates(rates)
print(p)
ggsave("EBD_Corr_moss_Co.pdf", p, width = 160, height = 160, units = "mm")

#################################################################################################
setwd("")
Pre <- c(1048, 1033, 1037, 1044, 1148, 1211, 1194, 1208, 1284, 1276, 1259, 1235, 1212, 1176, 1135, 1127, 1127, 1138, 1087, 1110, 1101, 1152, 1179, 1162, 1162, 1195, 1152, 1090, 1022, 970, 1070, 1014, 1102, 1159, 1200, 1187, 1202, 1270, 1329, 1334, 1258, 1238, 1252, 1233, 1144, 1177, 1286, 1358, 1366, 1346, 1329, 1316, 1324, 1311, 1301)
MAX_VAR_AGE = 540
NUM_INTERVALS = 10
Pre_age <- seq(from = 0, to = MAX_VAR_AGE, by = NUM_INTERVALS)
print(Pre_age)

env_data <- data.frame(
  age = rev(Pre_age),  # 反转时间轴
  Pre = rev(Pre)        # 反转CO2值
)
env_plot <- ggplot(env_data, aes(x = age, y = Pre)) +
  geom_line(color = "black", linewidth = 0.2) +  # 蓝色曲线
  #geom_point(color = "darkblue", size = 3) +     # 深蓝色点
  scale_x_reverse("Age (Ma)",                   # X轴反转（从大到小）
                  limits = c(max(env_data$age), min(env_data$age)),
                  breaks = seq(0, 540, by = 100)) +  # 设置刻度
  scale_y_continuous("Pre (mm)", 
                     limits = c(900, 1400)) +      # 设置Y轴范围
  ggtitle("Pre Concentration Over Time") +      # 标题
  theme_minimal(base_size = 6) +
  theme_void() +                                  # 黑白主题
  theme(
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.1, "cm"),
    plot.title = element_text(hjust = 0.5, size = 6, face = "bold"),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6)
  )

# 打印图形
print(env_plot)

# 保存图形
ggsave("Pre_curve.pdf", env_plot, width = 80, height = 80, units = "mm")

rates <- processDivRates(
  speciation_time_log = "EBD_Corr_speciation_times.log",
  speciation_rate_log = "EBD_Corr_speciation_rates.log",
  extinction_time_log = "EBD_Corr_extinction_times.log",
  extinction_rate_log = "EBD_Corr_extinction_rates.log",
  burnin=0.25)
p <- plotDivRates(rates)
print(p)
ggsave("EBD_Corr_liver_Pre.pdf", p, width = 160, height = 160, units = "mm")

rates <- processDivRates(
  speciation_time_log = "EBD_Corr_speciation_times.log",
  speciation_rate_log = "EBD_Corr_speciation_rates.log",
  extinction_time_log = "EBD_Corr_extinction_times.log",
  extinction_rate_log = "EBD_Corr_extinction_rates.log",
  burnin=0.25)
p <- plotDivRates(rates)
print(p)
ggsave("EBD_Corr_moss_Pre.pdf", p, width = 160, height = 160, units = "mm")

