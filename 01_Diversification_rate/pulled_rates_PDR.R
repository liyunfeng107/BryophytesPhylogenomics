#!/usr/bin/env Rscript
################################################################################
# Single-run PDR with AUTO or FIXED grid + Nbootstraps + spline smoothing
# Author: 十七 + ChatGPT
# Ref: Louca & Pennell 2020 Nature; Louca 2021/2022 (castor::fit_hbd_pdr_on_grid)
#
# Usage:
# Rscript pulled_rates_single_clade_PDR_boot_autogrid_spline.R \
#   <tree_file> <output_prefix> [threads=24] [runtime=900] [Nboot=200] \
#   [splines_degree=1] [auto_grid=TRUE] [spline_df=NA] [spline_spar=NA] [fixed_grid_points=30]
################################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript pulled_rates_single_clade_PDR_boot_autogrid_spline.R <tree_file> <output_prefix> [threads=24] [runtime=900] [Nboot=200] [splines_degree=1] [auto_grid=TRUE/FALSE] [spline_df=NA] [spline_spar=NA] [fixed_grid_points=30]")
}

tree_file       <- args[1]
output_prefix   <- args[2]
threads         <- ifelse(length(args)>=3,  as.integer(args[3]),  24)
max_runtime_sec <- ifelse(length(args)>=4,  as.integer(args[4]),  900)
Nboot           <- ifelse(length(args)>=5,  as.integer(args[5]),  200)
spl_deg         <- ifelse(length(args)>=6,  as.integer(args[6]),  1)
auto_grid       <- ifelse(length(args)>=7,  tolower(args[7]) %in% c("true","t","1","yes","y"), TRUE)
spline_df_in    <- ifelse(length(args)>=8 && args[8]!="NA", as.numeric(args[8]), NA_real_)
spline_spar_in  <- ifelse(length(args)>=9 && args[9]!="NA", as.numeric(args[9]), NA_real_)
fixed_grid_pts  <- ifelse(length(args)>=10 && args[10]!="NA", as.integer(args[10]), 30)

suppressMessages({
  if (!requireNamespace("castor", quietly=TRUE)) install.packages("castor", repos="https://cloud.r-project.org")
  if (!requireNamespace("ape", quietly=TRUE)) install.packages("ape", repos="https://cloud.r-project.org")
  library(castor); library(ape)
})

cat("💻 threads=", threads,
    " | runtime=", max_runtime_sec, "s",
    " | Nbootstraps=", Nboot,
    " | splines_degree=", spl_deg,
    " | auto_grid=", auto_grid,
    " | fixed_grid_points=", fixed_grid_pts,
    " | spline_df=", ifelse(is.na(spline_df_in),"auto",spline_df_in),
    " | spline_spar=", ifelse(is.na(spline_spar_in),"auto",spline_spar_in),
    "\n", sep="")

Sys.setenv(OMP_NUM_THREADS = threads)
Sys.setenv(CASTOR_NUM_THREADS = threads)

# ==== Load tree ====
tree <- read.tree(tree_file)
if (!is.ultrametric(tree)) stop("❌ Input tree is not ultrametric")
NTIPS <- length(tree$tip.label)
tree_age <- max(node.depth.edgelength(tree)) * 1.00001
cat("📂 Tree loaded: ", NTIPS, " tips | age=", round(tree_age,2), " Ma\n", sep="")

# ==== Candidate grids for auto selection ====
make_grid <- function(n) seq(0, tree_age, length.out = n)
grid_candidates <- list(
  gp3  = make_grid(3),
  gp4  = make_grid(4),
  gp5  = make_grid(5),
  gp8  = make_grid(8),
  gp10 = make_grid(10),
  gp12 = make_grid(12),
  gp14 = make_grid(14),
  gp15 = make_grid(15)
)

# ==== Helper: safe smoothing ====
smooth_series <- function(x, y, spline_df = NA_real_, spline_spar = NA_real_){
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 5L) return(y)
  ss_ok <- FALSE; yhat <- y
  try({
    if (!is.na(spline_df)) {
      ss <- stats::smooth.spline(x, y, df = spline_df)
    } else if (!is.na(spline_spar)) {
      ss <- stats::smooth.spline(x, y, spar = spline_spar)
    } else {
      ss <- stats::smooth.spline(x, y)
    }
    yhat <- stats::predict(ss, x = x)$y
    ss_ok <- TRUE
  }, silent = TRUE)
  if (ss_ok) return(yhat)
  k <- max(3L, as.integer(0.1 * length(y)) | 1L)
  return(stats::runmed(y, k = k, endrule = "median"))
}

# ==== Phase 1: Auto grid or fixed grid ====
if (auto_grid) {
  cat("🔎 Auto-selecting grid by BIC (no bootstrap)...\n")
  results <- list()
  for (gname in names(grid_candidates)) {
    age_grid <- grid_candidates[[gname]]
    cat("▶ Testing grid ", gname, " (", length(age_grid), " pts)... ", sep="")
    fit_try <- tryCatch(
      fit_hbd_pdr_on_grid(
        tree = tree,
        age_grid = age_grid,
        Ntrials = 10,
        Nthreads = threads,
        max_model_runtime = max_runtime_sec,
        splines_degree = spl_deg,
        condition = "auto",
        verbose = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(fit_try) || is.null(fit_try$fitted_PDR) || !is.finite(fit_try$BIC)) {
      cat("❌ failed\n")
    } else {
      cat("✅ AIC=", round(fit_try$AIC,1), " BIC=", round(fit_try$BIC,1), "\n", sep="")
      results[[gname]] <- list(grid = age_grid, fit = fit_try)
    }
  }
  if (length(results) == 0) stop("❌ All grid tests failed")
  best_name <- names(results)[ which.min(sapply(results, function(z) z$fit$BIC)) ]
  age_grid  <- results[[best_name]]$grid
  cat("🏆 Best grid: ", best_name, " | ", length(age_grid), " points | BIC=",
      round(results[[best_name]]$fit$BIC,1), "\n", sep="")
  grid_note <- paste0("best=", best_name, "(", length(age_grid), " pts)")
} else {
  age_grid <- seq(0, tree_age, length.out = fixed_grid_pts)
  grid_note <- paste0("fixed=", fixed_grid_pts, " pts")
  cat("ℹ️ Auto-grid disabled; using fixed grid of ", fixed_grid_pts, " points.\n", sep="")
}

# ==== Phase 2: final fit WITH bootstraps ====
cat("🚀 Final fit on grid (", grid_note, ") with Nbootstraps=", Nboot, " ...\n", sep="")
fit <- tryCatch(
  fit_hbd_pdr_on_grid(
    tree = tree,
    age_grid = age_grid,
    Nbootstraps = Nboot,
    Ntrials = 5,
    Nthreads = threads,
    max_model_runtime = max_runtime_sec,
    splines_degree = spl_deg,
    condition = "auto",
    verbose = FALSE
  ),
  error=function(e) NULL
)
if (is.null(fit) || is.null(fit$fitted_PDR)) stop("❌ Final fit failed")

# ==== Summarize bootstrap ====
n_pts <- length(age_grid)
Raw_PDR <- if (length(fit$fitted_PDR) == n_pts) fit$fitted_PDR else rep(NA_real_, n_pts)
median_pdr <- Raw_PDR
q05 <- q95 <- rep(NA_real_, n_pts)

if (!is.null(fit$bootstrap_estimates$PDR)) {
  boot_mat <- fit$bootstrap_estimates$PDR
  if (is.matrix(boot_mat) && nrow(boot_mat) > 1) {
    cat("✅ Extracted ", nrow(boot_mat), " bootstrap replicates\n", sep="")
    median_pdr <- apply(boot_mat, 2, median, na.rm = TRUE)
    q05 <- apply(boot_mat, 2, quantile, probs = 0.05, na.rm = TRUE)
    q95 <- apply(boot_mat, 2, quantile, probs = 0.95, na.rm = TRUE)
  }
}

cat("🎨 Smoothing (smooth.spline -> runmed fallback)...\n")
Smoothed <- smooth_series(age_grid, median_pdr, spline_df_in, spline_spar_in)
Lower_s  <- smooth_series(age_grid, q05, spline_df_in, spline_spar_in)
Upper_s  <- smooth_series(age_grid, q95, spline_df_in, spline_spar_in)

out_df <- data.frame(Time_Ma = age_grid, Raw_PDR, Median_PDR=median_pdr, Q05=q05, Q95=q95,
                     Smoothed, Lower_s, Upper_s)
out_csv  <- paste0(output_prefix, "_PDR_spline_summary.csv")
out_pdf  <- paste0(output_prefix, "_PDR_spline_plot.pdf")
out_rdat <- paste0(output_prefix, "_PDR_spline_workspace.RData")
write.csv(out_df, out_csv, row.names = FALSE)
save(fit, out_df, tree, age_grid, grid_note, file = out_rdat)
cat("💾 Saved: ", out_csv, " , ", out_rdat, "\n", sep="")

pdf(out_pdf, width = 7, height = 5)
ylim <- range(c(out_df$Lower_s, out_df$Upper_s, out_df$Smoothed), na.rm = TRUE)
plot(out_df$Time_Ma, out_df$Smoothed, type="n", xlim = rev(range(out_df$Time_Ma)),
     ylim = ylim, xlab = "Time before present (Ma)",
     ylab = expression("Pulled Diversification Rate (" * r[p] * ")"),
     main = paste0("PDR (", NTIPS, " tips; ", grid_note, ")"))
# 平滑置信带（深）
polygon(c(out_df$Time_Ma, rev(out_df$Time_Ma)),
        c(out_df$Lower_s, rev(out_df$Upper_s)),
        col = adjustcolor("deepskyblue3", 0.25), border = NA)
lines(out_df$Time_Ma, out_df$Smoothed, col = "blue", lwd = 2)
#grid()
dev.off()

cat("📊 Plot saved: ", out_pdf, "\n✅ Done (Nbootstraps=", Nboot, ")\n", sep="")
################################################################################

