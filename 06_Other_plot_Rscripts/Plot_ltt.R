# 安装需要的包（首次使用时运行）
install.packages(c("ape", "phytools", "ggplot2"))

# 加载包
library(ape)       # 处理系统发育树的基础包
library(phytools)  # 提供LTT计算函数
library(ggplot2)   # 可视化LTT曲线
setwd("")
######################################################################################################
tree <- read.tree("bry_og_time_genus_rename.tre")
ltt_result <- ltt(tree, plot = FALSE, log.lineages = FALSE)
time_points <- ltt_result$times
lineage_counts <- ltt_result$ltt
T_full <- max(time_points)# 整树“现在”的横坐标
desired_labels <- seq(from = 0, to = T_full, by = 100)  # 步长100，自动终止于≤T_full的最大整数
ticks_new <- T_full - desired_labels  # 对应原始时间坐标（距今时间 = T_full - 原始时间）

pdf("Nuclear_LTT_plot.pdf", width = 8, height = 6)
par(mar = c(5,5,2,2))
plot(
  time_points, lineage_counts,
  type = "l",  # 线条连接
  lwd = 2,     # 线宽
  col = "#999999",  # 线条颜色
  xlab = "Time before present (Ma)",
  ylab = "log(Lineage number)",
  main = "Bryohytes Nuclear Lineage-Through-Time",
  log="y", xaxt="n", 
  cex.lab = 1.2,  # 坐标轴标签大小
  cex.axis = 1    # 坐标轴刻度大小
)
axis(1, at = ticks_new, labels = desired_labels, cex.axis = 1)  
dev.off()  # 关闭绘图设备
# 6. 输出LTT数据（可选，用于进一步分析）
ltt_data <- data.frame(
  时间 = time_points,
  谱系数量 = exp(lineage_counts),  # 转换回原始数量（反对数）
  log_谱系数量 = lineage_counts
)
write.csv(ltt_data, "Nuclear_LTT_data.csv", row.names = FALSE)  # 保存为CSV

######################nuclear liverworts LTT#####################################################
tree <- read.tree("bry_og_time_genus_rename.tre")
liverworts <- getMRCA(tree, tip = c("Calobryales_Haplomitriaceae_Haplomitrium", 
                                    "Sphaerocarpales_Monocarpaceae_Monocarpus"))

clade_tree <- extract.clade(tree, node = liverworts)
clade_ltt <- ltt(clade_tree, plot = FALSE, log.lineages = FALSE)
time_points <- clade_ltt$times  # 时间点（从根到现在，根为0，现在为根年龄）
lineage_counts <- clade_ltt$ltt  # 对应时间点的谱系数量（已取对数）
T_full <- max(time_points)# 整树“现在”的横坐标
desired_labels <- seq(from = 0, to = T_full, by = 100)  # 步长100，自动终止于≤T_full的最大整数
ticks_new <- T_full - desired_labels  # 对应原始时间坐标（距今时间 = T_full - 原始时间）

# 绘制LTT曲线
pdf("Nuclear_liverworts_LTT_plot.pdf", width = 8, height = 6)
par(mar = c(5, 5, 2, 2))
plot(
  time_points, lineage_counts,
  type = "l",  # 线条连接
  lwd = 2,     # 线宽
  col = "#f8951e",  # 线条颜色
  xlab = "Time before present (Ma)",
  ylab = "log(Lineage number)",
  main = "Liverworts Nuclear Lineage-Through-Time",
  log="y", xaxt="n", 
  cex.lab = 1.2,  # 坐标轴标签大小
  cex.axis = 1    # 坐标轴刻度大小
)
axis(1, at = ticks_new, labels = desired_labels, cex.axis = 1)  
dev.off()  # 关闭绘图设备
# 保存目标支系的LTT数据
clade_ltt_data <- data.frame(
  Time_before_present = time_points,
  Lineage_number = exp(lineage_counts),  # 原始数量
  log_Lineage_number = lineage_counts
)
write.csv(clade_ltt_data, "Nuclear_liverworts_LTT_data.csv", row.names = FALSE)

######################nuclear mosses LTT#####################################################
tree <- read.tree("bry_og_time_genus_rename.tre")
mosses <- getMRCA(tree, tip = c("Takakiales_Takakiaceae_Takakia", 
                                    "Tetraphidales_Tetraphidaceae_Tetraphis"))
clade_tree <- extract.clade(tree, node = mosses)
clade_ltt <- ltt(clade_tree, plot = FALSE, log.lineages = FALSE)
time_points <- clade_ltt$times  # 时间点（从根到现在，根为0，现在为根年龄）
lineage_counts <- clade_ltt$ltt  # 对应时间点的谱系数量（已取对数）
T_full <- max(time_points)# 整树“现在”的横坐标
desired_labels <- seq(from = 0, to = T_full, by = 100)  # 步长100，自动终止于≤T_full的最大整数
ticks_new <- T_full - desired_labels  # 对应原始时间坐标（距今时间 = T_full - 原始时间）

# 绘制LTT曲线
pdf("Nuclear_mosses_LTT_plot.pdf", width = 8, height = 6)
par(mar = c(5, 5, 2, 2))
plot(
  time_points, lineage_counts,
  type = "l",  # 线条连接
  lwd = 2,     # 线宽
  col = "#5fbea1",  # 线条颜色
  xlab = "Time before present (Ma)",
  ylab = "log(Lineage number)",
  main = "Mosses Nuclear Lineage-Through-Time",
  log="y", xaxt="n", 
  cex.lab = 1.2,  # 坐标轴标签大小
  cex.axis = 1    # 坐标轴刻度大小
)
axis(1, at = ticks_new, labels = desired_labels, cex.axis = 1)  
dev.off()  # 关闭绘图设备

# 保存目标支系的LTT数据
clade_ltt_data <- data.frame(
  Time_before_present = time_points,
  Lineage_number = exp(lineage_counts),  # 原始数量
  log_Lineage_number = lineage_counts
)
write.csv(clade_ltt_data, "Nuclear_mosses_LTT_data.csv", row.names = FALSE)

######################liverworts LTT#####################################################
tree <- read.tree("liverwort_maker_time_genus_rename.tre")
ltt_result <- ltt(tree, plot = FALSE, log.lineages = FALSE)
time_points <- ltt_result$times
lineage_counts <- ltt_result$ltt
T_full <- max(time_points)# 整树“现在”的横坐标
desired_labels <- seq(from = 0, to = T_full, by = 100)  # 步长100，自动终止于≤T_full的最大整数
ticks_new <- T_full - desired_labels  # 对应原始时间坐标（距今时间 = T_full - 原始时间）

pdf("Liverwort_locus_LTT_plot.pdf", width = 8, height = 6)
plot(
  time_points, lineage_counts,
  type = "l",  # 线条连接
  lwd = 2,     # 线宽
  col = "#f8951e",  # 线条颜色
  xlab = "Time before present (Ma)",
  ylab = "log(Lineage number)",
  main = "Liverwort Locus Lineage-Through-Time",
  log="y", xaxt="n", 
  cex.lab = 1.2,  # 坐标轴标签大小
  cex.axis = 1    # 坐标轴刻度大小
)
axis(1, at = ticks_new, labels = desired_labels, cex.axis = 1)  
dev.off()  # 关闭绘图设备
# 6. 输出LTT数据（可选，用于进一步分析）
ltt_data <- data.frame(
  时间 = time_points,
  谱系数量 = exp(lineage_counts),  # 转换回原始数量（反对数）
  log_谱系数量 = lineage_counts
)
write.csv(ltt_data, "Liverwort_locus_data.csv", row.names = FALSE)  # 保存为CSV

######################mosses LTT#####################################################
tree <- read.tree("moss_maker_time_genus_rename.tre")
ltt_result <- ltt(tree, plot = FALSE, log.lineages = FALSE)
time_points <- ltt_result$times
lineage_counts <- ltt_result$ltt
T_full <- max(time_points)# 整树“现在”的横坐标
desired_labels <- seq(from = 0, to = T_full, by = 100)  # 步长100，自动终止于≤T_full的最大整数
ticks_new <- T_full - desired_labels  # 对应原始时间坐标（距今时间 = T_full - 原始时间）

pdf("Moss_locus_LTT_plot.pdf", width = 8, height = 6)
par(mar = c(5,5,2,2))
plot(
  time_points, lineage_counts,
  type = "l",  # 线条连接
  lwd = 2,     # 线宽
  col = "#5fbea1",  # 线条颜色
  xlab = "Time before present (Ma)",
  ylab = "log(Lineage number)",
  main = "Moss Locus Lineage-Through-Time",
  log="y", xaxt="n", 
  cex.lab = 1.2,  # 坐标轴标签大小
  cex.axis = 1    # 坐标轴刻度大小
)
axis(1, at = ticks_new, labels = desired_labels, cex.axis = 1)  
dev.off()  # 关闭绘图设备
# 6. 输出LTT数据（可选，用于进一步分析）
ltt_data <- data.frame(
  时间 = time_points,
  谱系数量 = exp(lineage_counts),  # 转换回原始数量（反对数）
  log_谱系数量 = lineage_counts
)
write.csv(ltt_data, "Moss_locus_data.csv", row.names = FALSE)  # 保存为CSV
