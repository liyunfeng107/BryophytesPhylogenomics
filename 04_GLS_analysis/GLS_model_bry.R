# -----------------------------
# 0. Packages
# -----------------------------
need_pkgs <- c(
  "dplyr", "tidyr", "purrr", "ggplot2", "nlme", "readr"
)

for (p in need_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(nlme)
library(readr)

# -----------------------------
# 1. User settings
# -----------------------------
setwd("")

################################################################################
# Angiosperm Vs Bryophyte
################################################################################
infile <- "macroevo_timeseries_bry.csv"
time_col <- "time_ma" 

response_vars_user <- c(
  "Moss_lct"  = "moss_Count",
  "Liver_lct"   = "liver_Count"
)
# Main predictor
predictor_var_user   <- "angiosperm_count"
predictor_label_user <- "Angiosperm_count"

# ---------------------------
# Environmental covariates
# ---------------------------
env_vars_final <- c("temperature", "vascular_plant_count")

# Pulse / dummy variables (not z-scored)
pulse_vars <- c("kpg")

# ---------------------------
# Null test settings
# Circular Shift with minimum shift distance
# ---------------------------
N_NULL <- 1000
SEED_NULL <- 142
MIN_SHIFT <- 10

# ---------------------------
# Output
# ---------------------------
STORE_FITS   <- FALSE
output_prefix <- "Angio_Bryo_LCT_sync"
source("source_GLS_model_bry.R")

################################################################################
# Angiosperm Vs Bryophyte order
################################################################################
infile <- "macroevo_timeseries_bry.csv"
time_col <- "time_ma"

response_vars_user <- c(
  "Moss_Hypnales_ltt"  = "Hypnales_count",
  "Liver_Lejeuneales_ltt"   = "Lejeuneales_count"
)

# Main predictor
predictor_var_user   <- "angiosperm_count"
predictor_label_user <- "Angiosperm_count"

# Environmental covariates, z-scored
env_vars_final <- c("temperature", "vascular_plant_count")

# Pulse / dummy variables (not z-scored)
pulse_vars <- c("kpg")

# ---------------------------
# Null test settings
# Circular Shift with minimum shift distance
# ---------------------------
N_NULL <- 1000
SEED_NULL <- 142
MIN_SHIFT <- 10

# ---------------------------
# Output
# ---------------------------
STORE_FITS   <- FALSE
output_prefix <- "Angio_Bryo_ord_LCT_sync"
source("source_GLS_model_bry.R")

