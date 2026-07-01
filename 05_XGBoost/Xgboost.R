#---------------- xgboost -----------------#

library (plyr)  # Version 1.8.4
library (tidyverse)  # Version 1.2.1
library (pdp)  # Version 0.7.0
library (xgboost)  # Version 0.82.1
library (ape)  # Version 5.3
library (phytools)  # Version 0.6.99
library (scales)  # Version 1.0


#-----------------#
### Loading data ##
#-----------------#
setwd("")
#wd <- "~/Desktop"   # Set the working directory
#setwd(wd)

db <- read.csv("Filtered_Data_For_PCA_moss.csv", h = T, stringsAsFactors = F)
db <- read.csv("Filtered_Data_For_PCA_liver.csv", h = T, stringsAsFactors = F)
db <- db[complete.cases(db),]

rownames(db) <- db$Species

db$net_div_median <- ifelse(db$net_div_median < 0, 0, db$net_div_median)
#db$net_div_median <- ifelse(db$lambda_median < 0, 0, db$lambda_median)

#-----------------------#
### Setting the model  ##
#-----------------------#


fmod <- formula (~ bio1 + bio2 + bio12 + bio15 + elev + ph + sand + srad6 + wind + AOO)
#---------------------------------------#
#### Tuning XGBOOST: Systematic way #####
#---------------------------------------#
#options(na.action='na.pass')
modmat <- model.matrix(fmod, db)[, -1]  # 创建模型矩阵，去掉截距列
dtrain <- xgb.DMatrix(data = modmat, label = db$net_div_median, missing = NA)

All_rmse.sist <- c()
Param_group.sist <- list()
Best_iter.sist <- c()

set.seed(123)

pardb <- expand.grid(eta = seq(0.1, 0.4, by = 0.05),
                     gamma = seq(0, 0.2, by = 0.05),
                     max_depth = 2:10,
                     subsample = seq(0.5, 0.9, by = 0.1))

for (j in 1:nrow(pardb)) {
  
  params <- list(booster = 'gbtree',
                 objective = 'reg:gamma',
                 eta = pardb[j, 'eta'],
                 gamma = pardb[j, 'gamma'],
                 max_depth = pardb[j, 'max_depth'],
                 subsample = pardb[j, 'subsample'])
  
  xgb.tune <- xgb.cv(params = params,
                     data = dtrain,
                     nrounds = 150,
                     nfold = 5,
                     metrics = list('rmse'),
                     showsd = T,
                     stratified = T,
                     verbose = F,
                     early_stopping_rounds = 50,
                     maximize = F)
  min_rmse <- min(xgb.tune$evaluation_log$test_rmse_mean)
  All_rmse.sist <- append(All_rmse.sist, min_rmse)
  Param_group.sist[[j]] <- params
  Best_iter.sist <- append(Best_iter.sist, xgb.tune$best_iteration)
  
  cat(paste("Tuning step", j), "\n")
}

(params.sist = Param_group.sist[[which.min(All_rmse.sist)]])
(best.iter <- Best_iter.sist[[which.min(All_rmse.sist)]])
min(All_rmse.sist)

#-----------------------------------#
#### Tuning XGBOOST: Random way #####
#-----------------------------------#

obj_functions <- c("reg:squarederror", "reg:gamma", "reg:tweedie")

All_rmse.rand <- c()
Param_group.rand <- list()
Best_iter.rand <- c()
Objective_function <- c()  # 存储使用的目标函数

set.seed(123)

for (i in 1:1000) {
  current_obj <- sample(obj_functions, 1)
  
  params <- list(
    booster = 'gbtree',
    objective = current_obj,
    eta = runif(1, 0.01, 0.3),      # 扩大eta范围
    gamma = runif(1, 0, 0.5),        # 扩大gamma范围
    max_depth = sample(3:12, 1),     # 扩大深度范围
    subsample = runif(1, 0.6, 1),    # 扩大采样范围
    colsample_bytree = runif(1, 0.6, 1),
    min_child_weight = sample(1:10, 1),
    lambda = runif(1, 0, 5),         # 添加L2正则化
    alpha = runif(1, 0, 1)           # 添加L1正则化
  )
  
  # 如果是tweedie回归，添加额外参数
  if (current_obj == "reg:tweedie") {
    params$tweedie_variance_power <- runif(1, 1.5, 2.0)
  }
  
  # 使用tryCatch捕获并处理异常
  xgb.tune <- tryCatch({
    xgb.cv(params = params,
           data = dtrain,
           nrounds = 200,
           nfold = 5,
           metrics = list('rmse'),
           showsd = T,
           stratified = T,
           verbose = F,
           early_stopping_rounds = 50,
           maximize = F)
  }, error = function(e) {
    cat("Error in CV for params:", unlist(params), "\n")
    cat("Error message:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(xgb.tune)) next
  
  min_rmse <- min(xgb.tune$evaluation_log$test_rmse_mean)
  All_rmse.rand[i] <- min_rmse
  Param_group.rand[[i]] <- params
  Best_iter.rand[i] <- xgb.tune$best_iteration
  Objective_function[i] <- current_obj  # 存储目标函数
  
  cat(paste("Tuning step", i, "| Obj:", current_obj, "| RMSE:", round(min_rmse, 4)), "\n")
}

# 找出最佳参数配置
best_idx <- which.min(All_rmse.rand)
(params.rand = Param_group.rand[[best_idx]])
(best.iter.rand <- Best_iter.rand[[best_idx]])
(best.objective <- Objective_function[best_idx])
min(All_rmse.rand)

# 分析不同目标函数的性能
obj_performance <- data.frame(
  Objective = Objective_function,
  RMSE = All_rmse.rand,
  Iteration = Best_iter.rand
)
# 查看每个目标函数的平均性能
aggregate(RMSE ~ Objective, data = obj_performance, FUN = mean)

# 可视化不同目标函数的性能分布
library(ggplot2)
ggplot(obj_performance, aes(x = Objective, y = RMSE, fill = Objective)) +
  geom_boxplot() +
  labs(title = "XGBoost Objective Function Performance Comparison",
       y = "CV RMSE") +
  theme_minimal()

# 使用最佳配置（无论系统性还是随机性调优）
if (exists("min_rmse_sist") && min_rmse_sist < min(All_rmse.rand)) {
  params <- params.sist
  best.iter <- best.iter.sist
} else {
  params <- params.rand
  best.iter <- best.iter.rand
}

#------------------------------------------------------#
#### Now cross-validating predictions with XGBOOST #####
#------------------------------------------------------#
niter <- 1000
set.seed(123)

meanbias <- vector()
medianbias <- vector()
lowci <- vector()
upci <- vector()
R2 <- vector()

for (i in 1:niter) {
  
  samp <- sample(1:nrow(db), round(0.2 * nrow(db)), replace = FALSE)
  
  db.train <- db[-samp, ]
  db.test <- db[samp, ]
  
  modmat.train <- as.matrix(model.matrix(fmod, db.train)[, -1])  # Ensure it's a matrix
  modmat.test <- as.matrix(model.matrix(fmod, db.test)[, -1])    # Ensure it's a matrix
  
  xgbmod = xgboost(
    data = modmat.train, 
    label = db.train$net_div_median,
    nrounds = best.iter, 
    params = params,  # 包含选定的目标函数
    verbose = 0
  )
  
  pred <- predict(xgbmod, newdata = modmat.test, ntreelimit = best.iter)
  meas <- db.test$net_div_median
  
  bias <- meas - pred
  
  meanbias[i] <- mean(bias)
  medianbias[i] <- median(bias)
  q <- qt(.975, df = length(bias) - 1)
  se <- sd(meas - pred) / sqrt(length(pred))
  lowci[i] <- mean(bias) - q * se
  upci[i] <- mean(bias) + q * se
  R2[i] <- summary(lm(pred ~ meas))$r.squared
  
  cat(paste("Model and prediction", i), "\n")
  
}

XGBoutput <- data.frame(medianbias = medianbias, meanbias = meanbias, lowcibias = lowci, upcibias = upci, R2 = R2)

#-------------------------------------#
#### Checking relevant indicators #####
#-------------------------------------#
ma <- mean(XGBoutput$meanbias)
ma
r <- mean(XGBoutput$R2)
round(r, 2)

grid <- expand.grid(bio1 = mean(db$bio1, na.rm = TRUE), bio2 = mean(db$bio2, na.rm = TRUE), 
                    bio12 = mean(db$bio12, na.rm = TRUE), bio15 = mean(db$bio15, na.rm = TRUE), 
                    srad6 = mean(db$srad6, na.rm = TRUE), sand = mean(db$sand, na.rm = TRUE), 
                    elev = mean(db$elev, na.rm = TRUE), wind = mean(db$wind, na.rm = TRUE), 
                    ph = mean(db$ph, na.rm = TRUE), AOO = mean(db$AOO, na.rm = TRUE))
#----------------------------------#
### Creating grid of predictions ###
#----------------------------------#
modmatgrid <- model.matrix (fmod, data = grid) [, -1]  # Ensure it's a matrix
modmatgrid <- as.matrix(modmatgrid)
# 检查 modmatgrid 的行数，确保它与 pred.grid 的行数匹配
cat("Number of rows in grid:", nrow(grid), "\n")
cat("Number of rows in modmatgrid:", nrow(modmatgrid), "\n")

##----------------------------------##
#### Bootstrapping and predicting ####
##----------------------------------##

modmat <- model.matrix (fmod, db) [,-1]  # Ensure it's a matrix

# 初始化空的 pred.grid 数据框，确保列数是 niter，行数是 modmatgrid 的行数
niter = 10000
set.seed(31)

xgbmod <- list()
rel.imp <- list()

# 初始化 pred.grid，确保它的维度与 modmatgrid 的行数匹配
pred.grid <- as.data.frame(matrix(ncol = niter, nrow = nrow(modmatgrid), 
                                  dimnames = list(NULL, paste0('Boot', 1:niter))))

# 查看 pred.grid 的维度
cat("Number of rows in pred.grid:", nrow(pred.grid), "\n")
cat("Number of columns in pred.grid:", ncol(pred.grid), "\n")


for (i in 1:niter) {
  
  set.seed(i)
  
  # 训练 XGBoost 模型
  xgbmod[[i]] <- xgboost(modmat, label = db$net_div_median,
                         nrounds = best.iter, params = params, verbose = 0, print_every = 1000)
  # 提取特征重要性
  rel.imp[[i]] <- xgb.importance(colnames(modmat), model = xgbmod[[i]])
  # 进行预测，并将结果赋值给 pred.grid
  pred.grid[, i] <- predict(xgbmod[[i]], newdata = modmatgrid)
  cat(paste("Bootstrapping the model, round", i), "\n")
}
####################################上面for运行失败运行下一面for#############################################

for (i in 1:niter) {
  
  set.seed(i)
  xgbmod[[i]] <- xgboost(modmat, label = db$net_div_median,
                         nrounds = best.iter, params = params, verbose = 0, print_every = 1000)
  rel.imp[[i]] <- xgb.importance(colnames(modmat), model = xgbmod[[i]])
  # 将数值向量转换为矩阵
  modmatgrid_matrix <- matrix(modmatgrid, nrow = 1)
  colnames(modmatgrid_matrix) <- names(modmatgrid)
  # 转换为 xgb.DMatrix 类型
  dtest <- xgb.DMatrix(data = modmatgrid_matrix)
  pred.grid[, i] <- predict(xgbmod[[i]], newdata = dtest)
  cat(paste("Bootstrapping the model, round", i), "\n")
}

# 检查 pred.grid 是否已正确填充
cat("Finished bootstrapping. Checking the first few rows of pred.grid:\n")
head(pred.grid)

#--------------------------------------------------------#
### Exploring the relative importance of the variables ###
#--------------------------------------------------------#
# 加载所需包
library(plyr)

# 假设 rel.imp 是已加载的数据列表
relimpdf <- do.call(rbind, rel.imp) %>%
  mutate(FeatureProc = substr(Feature, 1, 5))

# 使用 ddply() 进行分组和汇总
relimpdf_grouped <- ddply(relimpdf, .(Feature, FeatureProc), summarise,
                          meanGain = mean(Gain, na.rm = TRUE),
                          uppGain = quantile(Gain, 0.75, na.rm = TRUE),
                          lowGain = quantile(Gain, 0.25, na.rm = TRUE))

# 按 FeatureProc 分组并计算总和
relimpdf_grouped_final <- ddply(relimpdf_grouped, .(FeatureProc), summarise,
                                meanGain = sum(meanGain, na.rm = TRUE),
                                uppGain = sum(uppGain, na.rm = TRUE),
                                lowGain = sum(lowGain, na.rm = TRUE))

# 按 meanGain 排序
relimpdf_grouped_final <- relimpdf_grouped_final %>%
  arrange(desc(meanGain)) %>%
  as.data.frame()

# 查看最终结果
head(relimpdf_grouped_final)

write.csv(relimpdf_grouped_final, "XGboost_lambda_20250706.csv")
#---------------------------#
### Plotting predictions  ###
#---------------------------#

#只展示前十个重要变量
# 增加图形边距以便为 y 轴标签留出足够空间
par(mar = c(5, 15, 4, 2))  # 上、右、下、左边距
# 选择前10个 meanGain 最大的 FeatureProc
relimpdf_grouped_top10 <- head(relimpdf_grouped_final[order(relimpdf_grouped_final$meanGain, decreasing = TRUE), ], 10)

# 调整 barplot 代码：画前十个 FeatureProc 的条形图
names.arg = as.character(relimpdf_grouped_top10$FeatureProc)  # 修改为 top10 的 FeatureProc
col.bars = c(rep("gainsboro", 5), "deepskyblue3", "deepskyblue3","deepskyblue3", "deepskyblue3", "deepskyblue3")  # 前 8 个是默认颜色，最后两个是深蓝色

# 创建条形图
x <- barplot(rev(relimpdf_grouped_top10$meanGain), names.arg = rev(names.arg), xlab = "", xaxt="n", lwd = 1,
             col = col.bars, border = NA, las = 1, horiz = T, xlim = c(0, 0.30), cex.names = 1.5)

# 设置 x 轴标签
axis(side = 1, at = seq(0, 0.25, 0.05), lwd = 1.7, cex.axis = 1.2, tcl = -0.3, mgp = c(3, 0.5, 0), labels = c("0","5","10","15","20","25"))
mtext(side = 1, text = "Relative importance (%)", line = 2, cex = 1)

# 添加一条虚线，表示某个特定的阈值
lines(c(0.1, 0.1), c(-1,13), col = "grey4", lwd = 2, lty = 2)

# 用 for 循环绘制每个条形图的上下置信区间线和点
for(i in 1:nrow(relimpdf_grouped_top10)) {
  # 绘制上下置信区间的线
  lines(c(relimpdf_grouped_top10$uppGain[i], relimpdf_grouped_top10$lowGain[i]),
        c(rev(x)[i], rev(x)[i]), col = "grey8", lwd = 3)
  
  # 在条形图上绘制 meanGain 的点，使用 FeatureProc 的颜色
  points(x = relimpdf_grouped_top10$meanGain[i], rev(x)[i], cex = 1.5, bg = rev(col.bars)[i], pch = 21)
}

