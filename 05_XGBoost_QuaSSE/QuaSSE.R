#!/usr/bin/env Rscript
# ================================================================
# QuaSSE Pipeline
# ================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(diversitree)
  library(phytools)
  library(ape)
  library(geiger)
  library(nlme)
  library(parallel)
  library(ggplot2)
})

# ---------- 命令行参数 ----------
option_list <- list(
  make_option("--tree",  type="character", help="输入树文件 (.tre)"),
  make_option("--trait", type="character", help="trait 数据表 (csv)"),
  make_option("--col",   type="integer",   default=1, help="trait 数据列号 (默认1)"),
  make_option("--out",   type="character", default="QuaSSE_out", help="输出前缀"),
  make_option("--cores", type="integer",   default=4, help="并行线程数")
)
opt <- parse_args(OptionParser(option_list=option_list))

cat("📂 输入树文件:", opt$tree, "\n")
cat("📂 输入特征文件:", opt$trait, "\n")
cat("🧠 使用线程:", opt$cores, "\n\n")

# ---------- 读取数据 ----------
mytree <- read.tree(opt$tree)
mydata <- read.csv(opt$trait, row.names=1, header=TRUE)

# ---------- 检查名称匹配 ----------
comparison <- name.check(phy=mytree, data=mydata)
if (is.list(comparison)) {
  if (length(comparison$tree_not_data) > 0)
    mytree <- drop.tip(mytree, comparison$tree_not_data)
}
name.check(phy=mytree, data=mydata)

# ---------- 创建特征向量 ----------
states <- as.numeric(mydata[, opt$col])
names(states) <- row.names(mydata)
states.sd <- sd(states)

# ---------- 起始参数 ----------
p <- starting.point.quasse(mytree, states)
xr <- range(states) + c(-1, 1) * 20 * p["diffusion"]
linear.x <- make.linear.x(xr[1], xr[2])

make.models <- function(lambda, mu)
  make.quasse(mytree, states, states.sd, lambda, mu)

nodrift <- function(f) constrain(f, drift ~ 0)

# ---------- 模型定义 ----------
f.c <- make.models(constant.x, constant.x)
f.l <- make.models(linear.x, constant.x)
f.s <- make.models(sigmoid.x, constant.x)
f.h <- make.models(noroptimal.x, constant.x)

# ---------- 控制参数 ----------
control <- list(parscale = 0.1, reltol = 0.001)

# ---------- 常数模型拟合 ----------
mle.c <- find.mle(nodrift(f.c), p, lower=0, control=control, verbose=0)
p.c <- mle.c$par

# ---------- 起始参数 ----------
p.l <- c(p.c[1], l.m=0, p.c[2:3])
p.s <- c(p.c[1], p.c[1], mean(xr), 1, p.c[2:3]); names(p.s) <- argnames(nodrift(f.s))
p.h <- c(p.c[1], mean(xr), 1, 1, p.c[2:3]); names(p.h) <- argnames(nodrift(f.h))

# ---------- 并行设置 ----------
n_cores <- opt$cores
cl <- makeCluster(n_cores)
clusterEvalQ(cl, {
  library(diversitree)
  library(phytools)
  library(ape)
  library(geiger)
  library(nlme)
})

# ---------- 安全 MLE 封装 ----------
safe_mle <- function(f, par, ...) {
  tryCatch(find.mle(f, par, ...),
           error = function(e) {
             cat("❌ find.mle 出错:", conditionMessage(e), "\n")
             return(NULL)
           })
}

# ✅ 一定要导出 safe_mle !!!
clusterExport(cl, c("mytree", "states", "states.sd", "p", "p.l", "p.s", "p.h", 
                    "xr", "linear.x", "make.models", "nodrift", "control",
                    "constant.x", "sigmoid.x", "noroptimal.x",
                    "safe_mle"), 
              envir = environment())


# ---------- 安全 MLE 封装 ----------
safe_mle <- function(f, par, ...) {
  tryCatch(find.mle(f, par, ...),
           error = function(e) {
             cat("❌ find.mle 出错:", conditionMessage(e), "\n")
             return(NULL)
           })
}

# ---------- 并行拟合 ----------
results <- parLapply(cl, 1:3, function(i) {
  xr <- range(states) + c(-1,1) * 20 * p["diffusion"]
  linear.x <- make.linear.x(xr[1], xr[2])
  make.models <- function(lambda, mu) make.quasse(mytree, states, states.sd, lambda, mu)
  nodrift <- function(f) constrain(f, drift ~ 0)
  f.c <- make.models(constant.x, constant.x)
  f.l <- make.models(linear.x, constant.x)
  f.s <- make.models(sigmoid.x, constant.x)
  f.h <- make.models(noroptimal.x, constant.x)

  if (i == 1) {
    mle.l <- safe_mle(nodrift(f.l), p.l, control=control)
    if (is.null(mle.l)) return(NULL)
    mle.d.l <- safe_mle(f.l, coef(mle.l, TRUE), control=control)
    return(list(mle.l=mle.l, mle.d.l=mle.d.l))
  } else if (i == 2) {
    mle.s <- safe_mle(nodrift(f.s), p.s, control=control)
    if (is.null(mle.s)) return(NULL)
    mle.d.s <- safe_mle(f.s, coef(mle.s, TRUE), control=control)
    return(list(mle.s=mle.s, mle.d.s=mle.d.s))
  } else {
    mle.h <- safe_mle(nodrift(f.h), p.h, control=control)
    if (is.null(mle.h)) return(NULL)
    mle.d.h <- safe_mle(f.h, coef(mle.h, TRUE), control=control)
    return(list(mle.h=mle.h, mle.d.h=mle.d.h))
  }
})
stopCluster(cl)

# ---------- 提取结果 ----------
mle.l   <- if (!is.null(results[[1]])) results[[1]]$mle.l else NULL
mle.d.l <- if (!is.null(results[[1]])) results[[1]]$mle.d.l else NULL
mle.s   <- if (!is.null(results[[2]])) results[[2]]$mle.s else NULL
mle.d.s <- if (!is.null(results[[2]])) results[[2]]$mle.d.s else NULL
mle.h   <- if (!is.null(results[[3]])) results[[3]]$mle.h else NULL
mle.d.h <- if (!is.null(results[[3]])) results[[3]]$mle.d.h else NULL

# ---------- 模型比较 ----------
models <- list(constant=mle.c, linear=mle.l, sigmoidal=mle.s, hump=mle.h,
               drift.linear=mle.d.l, drift.sigmoidal=mle.d.s, drift.hump=mle.d.h)
models_ok <- models[!vapply(models, is.null, logical(1))]

if (length(models_ok) >= 2) {
  first_model <- models_ok[[1]]
  rest_models <- models_ok[-1]
  model.test <- do.call(anova, c(list(first_model), rest_models))
} else {
  one <- models_ok[[1]]
  model.test <- data.frame(
    Model = names(models_ok)[1],
    lnLik = one$lnLik,
    AIC   = one$aic
  )
}

write.csv(model.test, paste0(opt$out, "_modeltest.csv"))
save.image(file = paste0(opt$out, ".RData"))

# ---------- 选择最佳模型 ----------
best_model_idx <- which.min(model.test$AIC)
best_model_name <- rownames(model.test)[best_model_idx]
cat("\n🏆 最佳模型:", best_model_name, "\n")

# ---------- 绘图 ----------
load(paste0(opt$out, ".RData"))

if (best_model_name %in% c("linear", "drift.linear")) {
  params <- if (best_model_name == "linear") coef(mle.l) else coef(mle.d.l)
  lambda_fun <- function(x) params["l.c"] + params["l.m"] * x
  mu_fun     <- function(x) params["m.c"]

} else if (best_model_name %in% c("sigmoidal", "drift.sigmoidal")) {
  params <- if (best_model_name == "sigmoidal") coef(mle.s) else coef(mle.d.s)
  sigmoid_func <- function(x, l.y0, l.y1, l.xmid, l.r)
    l.y0 + (l.y1 - l.y0)/(1 + exp(-l.r * (x - l.xmid)))
  lambda_fun <- function(x) sigmoid_func(x, params["l.y0"], params["l.y1"], params["l.xmid"], params["l.r"])
  mu_fun     <- function(x) params["m.c"]

} else if (best_model_name %in% c("hump", "drift.hump")) {
  params <- if (best_model_name == "hump") coef(mle.h) else coef(mle.d.h)
  hump_func <- function(x, l.y0, l.y1, l.xmid, l.s2)
    l.y0 + l.y1 * exp(-(x - l.xmid)^2 / (2 * l.s2^2))
  lambda_fun <- function(x) hump_func(x, params["l.y0"], params["l.y1"], params["l.xmid"], params["l.s2"])
  mu_fun     <- function(x) params["m.c"]

} else {
  params <- coef(mle.c)
  lambda_fun <- function(x) rep(params["l.c"], length(x))
  mu_fun     <- function(x) params["m.c"]
}

x_range <- range(states)
x_values <- seq(from = x_range[1] - 1, to = x_range[2] + 1, length.out = 200)
lambda_vals <- pmax(1e-10, lambda_fun(x_values))
mu_vals <- pmax(1e-10, mu_fun(x_values))
delta_vals <- lambda_vals - mu_vals

plot_data <- data.frame(
  Trait_Value = x_values,
  Lambda = lambda_vals,
  Mu = mu_vals,
  Delta = delta_vals
)

p <- ggplot(plot_data, aes(x = Trait_Value)) +
  geom_line(aes(y = Lambda), color = "black", linewidth = 0.25) +
  geom_line(aes(y = Delta), color = "red", linewidth = 0.2, linetype = "dashed") +
  labs(x = opt$out, y = "Rate") +
  theme_minimal(base_size = 6) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.2),
        panel.grid = element_blank(),
        axis.text = element_text(size = 5))

ggsave(paste0(opt$out, "_fit_", best_model_name, ".pdf"), p, width = 2, height = 2, dpi = 300)

cat("\n✅ 完成。结果输出至:\n",
    paste0(opt$out, "_modeltest.csv\n"),
    paste0(opt$out, "_fit_", best_model_name, ".pdf\n"))
