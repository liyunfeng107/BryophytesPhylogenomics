setwd("")
###############################绘图功能函数##########################################################
library(ggplot2)
library(dplyr)
library(ggforce)
library(RevGadgets)

plot_log10_density <- function(
    df,
    rate_col = "rate",
    state_col = "observed_state",
    value_col = "value",
    smooth_adjust = 2,
    title = NULL
){
  # ==== 把列名转成 symbol，后面统一用 !! 引用 ====
  rate_sym  <- rlang::sym(rate_col)
  state_sym <- rlang::sym(state_col)
  value_sym <- rlang::sym(value_col)
  
  # 状态列转成 factor（这里直接用 [[ ]]，不涉及 .data）
  df[[state_col]] <- factor(df[[state_col]])
  
  # 1%–99% trimming
  df_trim <- df %>% 
    dplyr::group_by(!!rate_sym, !!state_sym) %>%
    dplyr::mutate(
      lo = stats::quantile(!!value_sym, 0.01, na.rm = TRUE),
      hi = stats::quantile(!!value_sym, 0.99, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter( (!!value_sym) >= lo, (!!value_sym) <= hi )
  
  # 每个 rate 面板的 x 范围
  panel_ranges <- df_trim %>%
    dplyr::group_by(!!rate_sym) %>%
    dplyr::summarise(
      xmin = min(!!value_sym, na.rm = TRUE) / 3,
      xmax = max(!!value_sym, na.rm = TRUE) * 3,
      .groups = "drop"
    )
  
  df_joined <- dplyr::left_join(df_trim, panel_ranges, by = rate_col)
  
  # 绘图
  p <- ggplot2::ggplot(
    df_joined,
    ggplot2::aes(x = !!value_sym, fill = !!state_sym)
  ) +
    ggplot2::geom_density(alpha = 0.6, adjust = smooth_adjust, color = NA) +
    
    ggforce::facet_row(
      facets = ggplot2::vars(!!rate_sym),
      scales = "free_x"
    ) +
    
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      x = "Rate (log10 scale)",
      y = "Density",
      title = title,
      fill = "State"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      legend.position   = c(0.9, 0.92),
      strip.background  = ggplot2::element_blank(),
      strip.text        = ggplot2::element_text(size = 10, face = "bold"),
      plot.title        = ggplot2::element_text(hjust = 0.5)
    )
  
  return(p)
}

#################################################################################
folder_name <- "Trait A" 
HiSSE_file <- paste0(folder_name, "/HiSSE_2A.log")
df <- processSSE(HiSSE_file)

p <- plot_log10_density(
  df,
  title = "Trait A"
)
print(p)
ggsave(paste0(folder_name, "_HiSSE_2A.pdf"),p, width=5, height=5)

HiSSE_file <- paste0(folder_name, "/HiSSE_2B.log")
df <- processSSE(HiSSE_file)
p <- plot_log10_density(
  df,
  title = "Trait A"
)
print(p)
ggsave(paste0(folder_name, "_HiSSE_2B.pdf"),p, width=5, height=5)

