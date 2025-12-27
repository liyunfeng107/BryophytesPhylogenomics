library(ape)
library(RevGadgets)
setwd("J:/苔藓性状重建/second/trait_recon/")

################################################################################
# Part 1. 读取结果数据
################################################################################
CHARACTER <- "Habitats"
fold <- "Habitats_mosses_ARD_output"
out_path <- CHARACTER
in_path <- fold
STATE_LABELS <- c("0" = "0", "1" = "1", "2" = "2")
tree_file <- paste0(CHARACTER,"_ARD_output/ARD.tree")

CHARACTER <- "Niche state"
fold <- "Niche_state_mosses_ARD_output"
out_path <- CHARACTER
in_path <- fold
STATE_LABELS <- c("0" = "0", "1" = "1", "2" = "2")
tree_file <- paste0(CHARACTER,"_ARD_output/ARD.tree")

CHARACTER <- "Capsule position"
fold <- "Capsule_position_HRM4_ARD_output"
out_path <- CHARACTER
in_path <- fold
STATE_LABELS <- c("0" = "0", "1" = "1", "2" = "2", "3" = "3", "4" = "4", "5" = "5",
                  "6" = "6", "7" = "7", "8" = "8", "9" = "9", "10" = "10", "11" = "11")
tree_file <- paste0(CHARACTER,"_HRM4_ARD_output/HRM4_ARD.tree")

CHARACTER <- "Capsule_position"
fold <- "Capsule_position_HRM2_ARD_output"
out_path <- CHARACTER
in_path <- fold
STATE_LABELS <- c("0" = "0", "1" = "1", "2" = "2", "3" = "3", "4" = "4", "5" = "5")
tree_file <- paste0(CHARACTER,"_HRM2_ARD_output/HRM2_ARD.tree")

CHARACTER <- "Habitats"
fold <- "Habitats_liverworts_ARD_output"
out_path <- CHARACTER
in_path <- fold
STATE_LABELS <- c("0" = "0", "1" = "1", "2" = "2")
tree_file <- paste0(CHARACTER,"_ARD_output/ARD.tree")

CHARACTER <- "Niche state"
fold <- "Niche_state_liverworts_ER_output"
out_path <- CHARACTER
in_path <- fold
STATE_LABELS <- c("0" = "0", "1" = "1", "2" = "2")
tree_file <- paste0(CHARACTER,"_ER_output/ER.tree")

CHARACTER <- "Fungal symbiosis status"
fold <- "Fungal_HRM2_ER_output"
out_path <- CHARACTER
in_path <- fold
STATE_LABELS <- c("0" = "0", "1" = "1", "2" = "2","3" = "3")
tree_file <- paste0(CHARACTER,"_HRM2_ER_output/HRM2_ER.tree")

CHARACTER <- "Sexual_systems"
fold <- "Sexual_systems_ARD_output"
out_path <- CHARACTER
in_path <- fold
STATE_LABELS <- c("0" = "0", "1" = "1")
tree_file <- paste0(CHARACTER,"_ARD_output/ARD.tree")

ase <- processAncStates(tree_file, state_labels = STATE_LABELS)
df <- dplyr::as_tibble(ase)
new_df <- df[, c(2, 5, 6, 7, 8, 9, 10, 11, 12)]
new_df1 <- new_df %>%
  filter(.[[1]] != 1007)

write.csv(new_df1, paste0(fold,"/summary_states.csv"), row.names = FALSE)
write.csv(new_df, paste0(fold,"/summary_states.csv"), row.names = FALSE)

################################################################################
# Part 2. 提取节点number以及height时间
################################################################################
# 1) 用 ER.tree 生成“完全一致”的 phylo（node 编号一致）
#er_tree_file <- tree_file
er_phy <- read.nexus(tree_file)   # RevBayes 的 tree TREE1 = ... 这种格式，ape::read.nexus 能读
# 2) 生成 Node_ancestor（边表）
node_re <- as.data.frame(er_phy$edge)
colnames(node_re) <- c("Ancestor_node", "Node")
# 3) 生成 Node_height（从根到现今的“高度/年龄”）
depths <- node.depth.edgelength(er_phy)
max_depth <- max(depths)
heights <- max_depth - depths

node_height <- data.frame(
  Node = seq_along(heights),
  Height = heights
)
# 4) 如你没有 HPD，就直接复用同一份（或干脆不用 node_hpd）
node_hpd <- node_height

################################################################################
# Part 3. 计算树的枝上状态信息
################################################################################

#处理祖先状态重建结果
find_MAP_df <- function(csv_file) {
  df <- read.csv(csv_file)
  df_out <- df %>%
    transmute(Node = node,
              State = anc_state_1,
              Max_prob = anc_state_1_pp)
  return(df_out)
}
#合并数据并识别状态转变事件
merge_ase_geo <- function(df_re, df_hpd, df_height, reconstruct_file, out_file) {
  node_state_rates <- find_MAP_df(reconstruct_file)
  
  # 合并当前节点的状态
  state <- left_join(df_re, node_state_rates, by = "Node")
  
  # 合并祖先节点的状态
  state_a <- left_join(
    state,
    node_state_rates %>%
      rename(
        Ancestor_node   = Node,
        Ancestor_max_prob = Max_prob,
        Ancestor_state = State
      ),
    by = "Ancestor_node"
  )
  # 合并当前节点的高度
  height <- left_join(state_a, df_height, by = "Node")
  
  # 合并祖先节点的高度
  height_a <- left_join(
    height,
    df_height %>%
      rename(Ancestor_node = Node, Ancestor_height = Height),
    by = "Ancestor_node"
  )
  
  # 添加状态转变事件
  append_geo_events(height_a, out_file)
  
  return(height_a)
}
#识别状态转变事件
bio_events <- function(state, ancestor_state) {
  s <- as.character(state)
  f <- as.character(ancestor_state)
  
  insitu <- ifelse(s == f, s, NA)
  shift_to <- ifelse(s != f, s, NA)
  shift_from <- ifelse(s != f, f, NA)
  
  return(list(insitu = insitu, shift_to = shift_to, shift_from = shift_from))
}
#将识别的事件信息添加到数据框并写入 CSV 文件
append_geo_events <- function(df, out_file) {
  events <- mapply(bio_events, df$State, df$Ancestor_state, SIMPLIFY = FALSE)
  events_df <- do.call(rbind, lapply(events, as.data.frame))
  colnames(events_df) <- c("No_shift", "Shift_to", "Shift_from")
  
  combined_df <- cbind(df, events_df)
  write.csv(combined_df, out_file, row.names = FALSE)
}

files <- file.path(in_path, "summary_states.csv")

for (file in files) {
  out_file <- paste0(out_path, "_", basename(file))
  merge_ase_geo(node_re, node_hpd, node_height, file, out_file)
}

########################## 线ages状态占比 + 转换事件曲线（修正版） ###########################
library(data.table)
library(dplyr)
library(ggplot2)

## 0. 读入数据（保持你原来的方式）
## 注意：这里假定 *_summary_states.csv 已由 merge_ase_geo() 写出
##       里面包含 Node, Ancestor_node, State, Ancestor_state, Height, Ancestor_height, Shift_to, Shift_from 等列
dt <- fread(paste0(CHARACTER, "_summary_states.csv"))

## 基础清洗：去掉高度缺失或祖先状态缺失的边
dt <- dt[!is.na(Ancestor_height) & !is.na(Height) & !is.na(Ancestor_state)]

## 确保 Ancestor_height > Height（老 → 新），如果搞反了就调换
dt[, `:=`(
  Height_min  = pmin(Ancestor_height, Height),
  Height_max  = pmax(Ancestor_height, Height)
)]
# 我们约定：Height_max = 较老的时间（靠近根），Height_min = 较新的时间（靠近现今）

################################################################################
# Part 1. 基于“边状态”的 LTT-style 状态占比曲线（lineage proportion）
################################################################################

# ---- 1. 构造边区间：每条边一段 [Height_min, Height_max] ----
# 我们把“边的状态”定义为祖先节点状态 Ancestor_state（即进入该边时的状态）
# 这样更接近连续过程：状态在进入分支前已确定
make_lineage_blocks <- function(df) {
  # df 需要包含：Ancestor_state, Height_min, Height_max
  lapply(seq_len(nrow(df)), function(i) {
    list(
      start = df$Height_max[i],  # 老的时间
      end   = df$Height_min[i]   # 新的时间
    )
  })
}

# ---- 2. 扫描算法：计算某个状态在任意时间 t 的“活跃边占比” ----
events_through_time_ratio <- function(blocks_num, blocks_den) {
  make_edges <- function(blks, tag) {
    rbindlist(lapply(blks, function(b) {
      data.table(
        Time     = c(b$start, b$end),
        Happened = c(paste0("begin_", tag), paste0("end_", tag))
      )
    }))
  }
  edges <- rbind(
    make_edges(blocks_num, "num"),
    make_edges(blocks_den, "den")
  )
  # 去掉 NA，按时间从老到新（大→小）排序
  edges <- edges[!is.na(Time)][order(-Time)]
  
  e_num <- 0L
  e_den <- 0L
  last_t <- NA_real_
  res <- data.table(Time = numeric(), Ratio = numeric())
  
  for (i in seq_len(nrow(edges))) {
    r <- edges[i]
    
    # 在时间发生改变时，记录上一时间点的比例
    if (!is.na(last_t) && r$Time != last_t) {
      ratio <- if (e_den == 0) {
        if (e_num == 0) 0 else 1
      } else {
        e_num / e_den
      }
      res <- rbind(res, data.table(Time = last_t, Ratio = ratio))
    }
    
    # 更新计数
    if (r$Happened == "begin_num") e_num <- e_num + 1L
    if (r$Happened == "end_num")   e_num <- e_num - 1L
    if (r$Happened == "begin_den") e_den <- e_den + 1L
    if (r$Happened == "end_den")   e_den <- e_den - 1L
    
    last_t <- r$Time
  }
  
  # 最后一个时间点
  ratio <- if (e_den == 0) {
    if (e_num == 0) 0 else 1
  } else {
    e_num / e_den
  }
  res <- rbind(res, data.table(Time = last_t, Ratio = ratio))
  
  res[]
}

# ---- 3. 针对单个状态计算 LTT 状态占比 ----
compute_lineage_ratio <- function(df, target_state) {
  
  # 分母：所有边
  den_df <- df[!is.na(Merged_Ancestor_state)]
  
  # 分子：合并后的祖先状态 == target_state
  num_df <- df[Merged_Ancestor_state == target_state]
  
  if (nrow(num_df) == 0L) {
    warning("State ", target_state, " 在合并后不存在")
    return(NULL)
  }
  
  blocks_den <- make_lineage_blocks(den_df)
  blocks_num <- make_lineage_blocks(num_df)
  
  dt_ratio <- events_through_time_ratio(blocks_num, blocks_den)
  dt_ratio[, State := factor(target_state)]
  dt_ratio
}


# ---- 4. 多状态叠加：得到整体 LTT 状态占比曲线 ----
compute_lineage_ratios_multi <- function(df, states) {
  res_list <- lapply(states, function(s) compute_lineage_ratio(df, s))
  res_list <- res_list[!sapply(res_list, is.null)]
  if (length(res_list) == 0L) return(NULL)
  rbindlist(res_list)
}

################################################################################
# Part 2. 状态转换事件累积曲线（shift-through-time）
################################################################################

# 思想：
#   每一条边中，如果 Ancestor_state != State，则视为该边上发生了一个“状态转变事件”
#   转变时间采用该边中点：(Height_max + Height_min) / 2
#   然后计算：从老到新的过程中，每个状态“被转入”的累积比例
#   即： Cumulative_shifts_to_state(t) / Total_shifts_to_state

compute_shift_ratio <- function(df, target_state) {
  
  shift_df <- df[
    !is.na(Merged_State) & !is.na(Merged_Ancestor_state) &
      (Merged_State != Merged_Ancestor_state) &
      (Merged_State == target_state)
  ]
  
  if (nrow(shift_df) == 0L) {
    warning("State ", target_state, " 没有转变事件")
    return(NULL)
  }
  
  shift_df[, Event_time := (Height_max + Height_min) / 2]
  shift_df <- shift_df[order(-Event_time)]
  
  total <- nrow(shift_df)
  data.table(
    Time  = shift_df$Event_time,
    Ratio = seq_len(total) / total,
    State = factor(target_state)
  )
}

compute_shift_ratios_multi <- function(df, states) {
  res_list <- lapply(states, function(s) compute_shift_ratio(df, s))
  res_list <- res_list[!sapply(res_list, is.null)]
  if (length(res_list) == 0L) return(NULL)
  rbindlist(res_list)
}

################################################################################
# Part 3. 绘图封装（LTT 状态占比 + 转换累积曲线）
################################################################################

# ---- 1. 颜色向量（你可以沿用自己设定的 colors_vec） ----
# 示例：
# colors_vec <- c('#fa9a7c','#82dee6','#4083d5')

# ---- 2. 绘制 LTT 状态占比曲线 ----
plot_lineage_LTT <- function(all_ratios, colors_vec, title = CHARACTER) {
  if (is.null(all_ratios) || nrow(all_ratios) == 0L) {
    stop("LTT 状态占比结果为空。")
  }
  
  ggplot(subset(all_ratios, Ratio != 0),
         aes(x = Time, y = Ratio, color = State)) +
    geom_line(size = 0.5) +
    scale_x_reverse() +
    scale_color_manual(
      name   = "State",
      values = colors_vec
    ) +
    labs(x = "Time (Ma)", y = "Lineage proportion", title = title) +
    theme_minimal() +
    theme(
      axis.ticks        = element_line(color = "black", size = 0.4),
      axis.ticks.length = unit(0.3, "cm"),
      axis.ticks.y.right= element_line(color = "black", size = 0.4),
      axis.line         = element_line(color = "black", size = 0.4)
    )
}

# ---- 3. 绘制状态转换累积曲线 ----
plot_shift_curve <- function(all_shift, colors_vec, title = CHARACTER) {
  if (is.null(all_shift) || nrow(all_shift) == 0L) {
    stop("Shift 累积曲线结果为空。")
  }
  
  ggplot(all_shift,
         aes(x = Time, y = Ratio, color = State)) +
    geom_line(size = 0.5) +
    scale_x_reverse() +
    scale_color_manual(
      name   = "State",
      values = colors_vec
    ) +
    labs(x = "Time (Ma)", y = "Cumulative fraction of shifts", title = title) +
    theme_minimal() +
    theme(
      axis.ticks        = element_line(color = "black", size = 0.4),
      axis.ticks.length = unit(0.3, "cm"),
      axis.ticks.y.right= element_line(color = "black", size = 0.4),
      axis.line         = element_line(color = "black", size = 0.4)
    )
}

################################################################################
# Part 4. 实际调用示例（你可以像以前一样按性状换 states / colors）
################################################################################
#Capsule_position_HRM4_ARD_output
merge_map <- c(
  "0"="0", "3"="0", "6"="0", "9"="0",
  "1"="1", "4"="1", "7"="1", "10"="1",
  "2"="2", "5"="2", "8"="2", "11"="2"
)
dt[, Merged_State := merge_map[ as.character(State) ]]
dt[, Merged_Ancestor_state := merge_map[ as.character(Ancestor_state) ]]

# 两个状态合并
merge_map <- c("0" = "0", "1" = "1", 
               "2" = "0", "3" = "1")
dt[, Merged_State := merge_map[ as.character(State) ]]
dt[, Merged_Ancestor_state := merge_map[ as.character(Ancestor_state) ]]

# 三个状态合并
merge_map <- c("0" = "0", "1" = "1", "2" = "2", 
               "3" = "0", "4" = "1", "5" = "2")
dt[, Merged_State := merge_map[ as.character(State) ]]
dt[, Merged_Ancestor_state := merge_map[ as.character(Ancestor_state) ]]

# 如果你不想合并，只需设置：
merge_map <- setNames(as.character(unique(dt$State)), as.character(unique(dt$State)))
# 在 dt 中添加合并后的状态
dt[, Merged_State := merge_map[ as.character(State) ]]
dt[, Merged_Ancestor_state := merge_map[ as.character(Ancestor_state) ]]

#Fungal
colors_vec <- c("0" = "#4083d5", "1" = "#fa9a7c")

#bry_habitat_liverworts_ARD_output
colors_vec <- c("0" = "#4083d5", "1" = "#fa9a7c",  
                "2" = "#82dee6")

#Niche_state_liverworts_ER_output
colors_vec <- c("0" = "#fa9a7c", "1" = "#4083d5", 
                "2" = "#82dee6")

#Sexual_systems_ARD_output
colors_vec <- c("0" = "#fa9a7c", "1" = "#4083d5" )

#Capsule_position_HRM4_ARD_output
colors_vec <- c("0" = "#fa9a7c", "1" = "#82dee6", "2" = "#4083d5")

#Capsule_position
colors_vec <- c("0" = "#fa9a7c", "1" = "#4083d5")
## 举例：三态性状（0,1,2），画 LTT 状态占比
states_LTT <- c(0, 1)
states_LTT <- c(0, 1, 2)
lineage_ratios <- compute_lineage_ratios_multi(dt, states = states_LTT)

p_LTT <- plot_lineage_LTT(lineage_ratios, colors_vec, title = CHARACTER)
ggsave(paste0(CHARACTER, "_lineages_LTT.pdf"),
       p_LTT, width = 6, height = 4, limitsize = FALSE)

