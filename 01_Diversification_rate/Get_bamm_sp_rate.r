setwd("")
library(ape)
library(phytools)
library(BAMMtools)

data <- read.csv("A_tree_tip_name.csv", header=T, sep = ",", stringsAsFactors = T)
rownames(data) <- data$x

# 创建存储 lambda 和 mu 的数据框  
df.lambda <- as.data.frame(matrix(NA, ncol = 1, nrow = nrow(data)))  
rownames(df.lambda) <- rownames(data)  
  
df.mu <- as.data.frame(matrix(NA, ncol = 1, nrow = nrow(data)))  
rownames(df.mu) <- rownames(data)  

# 读取第一次循环对应的树文件和事件数据  
tree_file <- paste("A.tre", sep = "")  
events_file <- paste("event_data.txt", sep = "")  
  
tree <- read.tree(tree_file)  
events <- read.csv(events_file)  
  
# 使用树文件和事件数据获取事件数据对象  
ed <- getEventData(tree, events, burnin = 0.1)  
  
# 获取尖端速率  
tip.rates <- getTipRates(ed)  
  
# 提取 lambda 和 mu 的平均值，并匹配到 data 的 rownames 上  
lambda <- tip.rates$lambda.avg  
lambda <- lambda[rownames(data)]  
df.lambda[, 1] <- lambda  
  
mu <- tip.rates$mu.avg  
mu <- mu[rownames(data)]  
df.mu[, 1] <- mu  
  
# 计算 lambda 和 mu 的中位数（虽然这里只有一列，但中位数计算仍然适用）  
df.lambda$Median <- apply(df.lambda, 1, median, na.rm = TRUE)  
df.mu$Median <- apply(df.mu, 1, median, na.rm = TRUE)  
  
# 创建包含中位数信息的数据框  
median_df <- data.frame(lambda_median = df.lambda$Median, mu_median = df.mu$Median,  
                        net_div_median = df.lambda$Median - df.mu$Median,  
                        row.names = rownames(data))

write.csv(median_df, "Rates.csv")

