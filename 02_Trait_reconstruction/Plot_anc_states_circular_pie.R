################################################################################
library(RevGadgets)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
setwd("")

classification_df <- read.csv("mosses_class_clean.csv",stringsAsFactors = FALSE)
CHARACTER <- "Peristome type"
custom_colors <- c("Slow gain of Haplolepideous" = "#3A73C1",  
                   "Slow gain of Diplolepideous" = "#F28A70", 
                   "Slow gain of Absent" = "#9C84C7",
                   "Fast gain of Haplolepideous" = "#c5daf2",
                   "Fast gain of Diplolepideous" = "#fccbbc",
                   "Fast gain of Absent" = "#C0A8D8",
                   "Fast gain of Haplolepideous" = "gray50",
                   "Fast gain of Haplolepideous" = "gray50",
                   "Fast gain of Haplolepideous" = "gray50",
                   "Fast gain of Haplolepideous" = "gray50",
                   "Fast gain of Haplolepideous" = "gray50",
                   "Fast gain of Haplolepideous" = "gray50",
                   "Fast gain of Haplolepideous" = "gray50"
)   

NUM_STATES <- 6
STATE_LABELS <- c("0" = "Slow gain of Haplolepideous", "1" = "Slow gain of Diplolepideous",
                  "2" = "Slow gain of Absent", "3" = "Fast gain of Haplolepideous",
                  "4" = "Fast gain of Diplolepideous", "5" = "Fast gain of Absent")

tree_file <- paste0(CHARACTER,"_HRM2_ARD_output/HRM2_ARD.tree")
ase <- processAncStates(tree_file, state_labels = STATE_LABELS)
  p <- plotAncStatesPie(t = ase,
                        pie_colors = custom_colors,
                        layout= "fan",
                        #layout= "circular",
                        state_transparency = 1.0,
                        open.angle  = 10,
                        ladderize = FALSE,
                        tip_labels=FALSE,
                        node_pie_size = 0.4,
                        tip_pie_size = 0.1,
                        size = 0.1, # 控制枝条粗细
                        colour = "gray50"  # 控制枝条颜色
                        ) +
  theme(legend.position = c(0.92,0.92), legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.key.size= unit(0.2, "cm"),legend.key = element_blank())
  p2 <- rotate_tree(p, 10)
  p3 <- p2 +
    geom_fruit(
      data     = classification_df,
      geom     = geom_tile,
      mapping  = aes(y= label, x = 1, fill = classify),
      offset   = 2,
      pwidth   = 1.05,
      position = "identity",
      color    = NA,
      width = 6
    ) +
    scale_fill_manual(values = c("Bryidae" = "#5fbea1", 
                                 "Dicranidae" = "#4DAF4A",
                                 "Other mosses" = "#5c6e83"
                                 
    ))
ggsave(paste0(CHARACTER,"_ER_Pie.pdf"), p3, width = 4, height = 6, limitsize = FALSE)

#######################################################################################################  
classification_df <- read.csv("liverworts_class.csv",stringsAsFactors = FALSE)
CHARACTER <- "Fungal symbiosis status"
custom_colors <- c("Slow gain of Symbiotic" = "#3A73C1", "Slow gain of None" = "#F28A70",
                   "Fast gain of Symbiotic" = "#c5daf2", "Fast gain of None" = "#fccbbc",
                   "Fast gain of multiple capsule shapes1" = "gray50",
                   "Fast gain of multiple capsule shapes2" = "gray50",
                   "Fast gain of multiple capsule shapes3" = "gray50",
                   "Fast gain of multiple capsule shapes4" = "gray50",
                   "Fast gain of multiple capsule shapes5" = "gray50"
) 

NUM_STATES <- 4
STATE_LABELS <- c("0" = "Slow gain of Symbiotic", "1" = "Slow gain of None",
                  "2" = "Fast gain of Symbiotic", "3" = "Fast gain of None")

tree_file <- paste0(CHARACTER,"_HRM2_ER_output/HRM2_ER.tree")
ase <- processAncStates(tree_file, state_labels = STATE_LABELS)
p <- plotAncStatesPie(t = ase,
                      pie_colors = custom_colors,
                      layout= "fan",
                      #layout= "circular",
                      state_transparency = 1.0,
                      open.angle  = 10,
                      ladderize = FALSE,
                      tip_labels=FALSE,
                      node_pie_size = 0.4,
                      tip_pie_size = 0.1,
                      size = 0.1, # 控制枝条粗细
                      colour = "gray50"  # 控制枝条颜色
) +
  # modify legend location using ggplot2
  theme(legend.position = c(0.92,0.92), legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.key.size= unit(0.2, "cm"),legend.key = element_blank())
p2 <- rotate_tree(p, 10)
p3 <- p2 +
  geom_fruit(
    data     = classification_df,
    geom     = geom_tile,
    mapping  = aes(y= label, x = 1, fill = classify),
    offset   = 2,
    pwidth   = 1.05,
    position = "identity",
    color    = NA,
    width = 6
  ) +
  scale_fill_manual(values = c("Jungermanniidae" = "#f8951e",
                               "Other liverworts" = "#cc0000"
  ))

ggsave(paste0(CHARACTER,"_HRM2_ER_Pie.pdf"), p3, width = 4, height = 6, limitsize = FALSE)

#######################################################################################################  
classification_df <- read.csv("liverworts_class.csv",stringsAsFactors = FALSE)
CHARACTER <- "Sexual systems"
custom_colors <- c("Bisexual" = "#F28A70", "Unisexual" = "#3A73C1" )

NUM_STATES <- 2
STATE_LABELS <- c("0" = "Bisexual", "1" = "Unisexual")

tree_file <- paste0(CHARACTER,"_ARD_output/ARD.tree")
ase <- processAncStates(tree_file, state_labels = STATE_LABELS)
p <- plotAncStatesPie(t = ase,
                      pie_colors = custom_colors,
                      layout= "fan",
                      #layout= "circular",
                      state_transparency = 1.0,
                      open.angle  = 10,
                      ladderize = FALSE,
                      tip_labels=FALSE,
                      node_pie_size = 0.4,
                      tip_pie_size = 0.1,
                      size = 0.1, # 控制枝条粗细
                      colour = "gray50"  # 控制枝条颜色
) +
  # modify legend location using ggplot2
  theme(legend.position = c(0.92,0.92), legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.key.size= unit(0.2, "cm"),legend.key = element_blank())
p2 <- rotate_tree(p, 10)
p3 <- p2 +
  geom_fruit(
    data     = classification_df,
    geom     = geom_tile,
    mapping  = aes(y= label, x = 1, fill = classify),
    offset   = 2,
    pwidth   = 1.05,
    position = "identity",
    color    = NA,
    width = 6
  ) +
  scale_fill_manual(values = c("Jungermanniidae" = "#f8951e",
                               "Other liverworts" = "#cc0000"
  ))

ggsave(paste0(CHARACTER,"_ARD_Pie.pdf"), p3, width = 4, height = 6, limitsize = FALSE)

#######################################################################################################  
classification_df <- read.csv("mosses_class.csv",stringsAsFactors = FALSE)
CHARACTER <- "Capsule position"
custom_colors <- c("VSlow gain of Pleurocarpous" = "#F28A70", "VSlow gain of Acrocarpous" = "#3A73C1",
                   "VSlow gain of Cladocarpous" = "#9C84C7", 
                   "Slow gain of Pleurocarpous" = "#fccbbc", 
                   "Slow gain of Acrocarpous" = "#c5daf2", "Slow gain of Cladocarpous" = "#C0A8D8",
                    "Fast gain of Pleurocarpous" = "#fccbbc", "Fast gain of Acrocarpous" = "#c5daf2",
                   "Fast gain of Cladocarpous" =  "#C0A8D8", 
                   "VFast gain of Pleurocarpous" = "#fccbbc", 
                    "VFast gain of Acrocarpous" = "#c5daf2", "VFast gain of Cladocarpous" = "#C0A8D8",
                   "Fast gain of multiple capsule shapes1" = "gray50",
                   "Fast gain of multiple capsule shapes2" = "gray50",
                   "Fast gain of multiple capsule shapes3" = "gray50",
                   "Fast gain of multiple capsule shapes4" = "gray50",
                   "Fast gain of multiple capsule shapes3" = "gray50",
                   "Fast gain of multiple capsule shapes4" = "gray50",
                   "Fast gain of multiple capsule shapes5" = "gray50"
) 

NUM_STATES <- 12
STATE_LABELS <- c("0" = "VSlow gain of Pleurocarpous", "1" = "VSlow gain of Acrocarpous",
                  "2" = "VSlow gain of Cladocarpous", "3" = "Slow gain of Pleurocarpous",
                  "4" = "Slow gain of Acrocarpous", "5" = "Slow gain of Cladocarpous",
                  "6" = "Fast gain of Pleurocarpous","7" = "Fast gain of Acrocarpous", 
                  "8" = "Fast gain of Cladocarpous", "9" = "VFast gain of Pleurocarpous",
                  "10" = "VFast gain of Acrocarpous", "11" = "VFast gain of Cladocarpous")

tree_file <- paste0(CHARACTER,"_HRM4_ARD_output/HRM4_ARD.tree")
ase <- processAncStates(tree_file, state_labels = STATE_LABELS)
  p <- plotAncStatesPie(t = ase,
                        pie_colors = custom_colors,
                        layout= "fan",
                        #layout= "circular",
                        state_transparency = 1.0,
                        open.angle  = 10,
                        ladderize = FALSE,
                        tip_labels=FALSE,
                        node_pie_size = 0.4,
                        tip_pie_size = 0.1,
                        size = 0.1, # 控制枝条粗细
                        colour = "gray50"  # 控制枝条颜色
  ) +
    # modify legend location using ggplot2
    theme(legend.position = c(0.92,0.92), legend.text = element_text(size = 5),
          legend.title = element_text(size = 5),
          legend.key.size= unit(0.2, "cm"),legend.key = element_blank())
  p2 <- rotate_tree(p, 10)
  p3 <- p2 +
    geom_fruit(
      data     = classification_df,
      geom     = geom_tile,
      mapping  = aes(y= label, x = 1, fill = classify),
      offset   = 2,
      pwidth   = 1.05,
      position = "identity",
      color    = NA,
      width = 6
    ) +
    scale_fill_manual(values = c("Bryidae" = "#5fbea1", 
                                 "Dicranidae" = "#4DAF4A",
                                 "Other mosses" = "#5c6e83"
                                 
    ))
  
ggsave(paste0(CHARACTER,"_HRM4_ARD_Pie.pdf"), p3, width = 4, height = 6, limitsize = FALSE)

#######################################################################################################  
classification_df <- read.csv("liverworts_class.csv",stringsAsFactors = FALSE)
CHARACTER <- "Habitats"
custom_colors <- c("Terrestrial" = "#3A73C1", "Epiphytic" = "#F28A70",  
                   "Aquatic" = "#9C84C7")

NUM_STATES <- 3
STATE_LABELS <- c("0" = "Terrestrial", "1" = "Epiphytic",
                 "2" = "Aquatic")

tree_file <- paste0(CHARACTER,"_ARD_output/ARD.tree")
ase <- processAncStates(tree_file, state_labels = STATE_LABELS)
p <- plotAncStatesPie(t = ase,
                      pie_colors = custom_colors,
                      layout= "fan",
                      #layout= "circular",
                      state_transparency = 1.0,
                      open.angle  = 10,
                      ladderize = FALSE,
                      tip_labels=FALSE,
                      node_pie_size = 0.4,
                      tip_pie_size = 0.1,
                      size = 0.1, # 控制枝条粗细
                      colour = "gray50"  # 控制枝条颜色
) +
  # modify legend location using ggplot2
  theme(legend.position = c(0.92,0.92), legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.key.size= unit(0.2, "cm"),legend.key = element_blank())
p2 <- rotate_tree(p, 10)
p3 <- p2 +
  geom_fruit(
    data     = classification_df,
    geom     = geom_tile,
    mapping  = aes(y= label, x = 1, fill = classify),
    offset   = 2,
    pwidth   = 1.05,
    position = "identity",
    color    = NA,
    width = 6
  ) +
  scale_fill_manual(values = c("Jungermanniidae" = "#f8951e",
                               "Other liverworts" = "#cc0000"
  ))

ggsave(paste0(CHARACTER,"_ARD_Pie.pdf"), p3, width = 4, height = 6, limitsize = FALSE)

#######################################################################################################  
classification_df <- read.csv("mosses_class.csv",stringsAsFactors = FALSE)
CHARACTER <- "Habitats"
custom_colors <- c("Terrestrial" = "#3A73C1", "Epiphytic" = "#F28A70",  
                   "Aquatic" = "#9C84C7")

NUM_STATES <- 3
STATE_LABELS <- c("0" = "Terrestrial", "1" = "Epiphytic",
                 "2" = "Aquatic")

tree_file <- paste0(CHARACTER,"_ARD_output/ARD.tree")
ase <- processAncStates(tree_file, state_labels = STATE_LABELS)
p <- plotAncStatesPie(t = ase,
                      pie_colors = custom_colors,
                      layout= "fan",
                      #layout= "circular",
                      state_transparency = 1.0,
                      open.angle  = 10,
                      ladderize = FALSE,
                      tip_labels=FALSE,
                      node_pie_size = 0.4,
                      tip_pie_size = 0.1,
                      size = 0.1, # 控制枝条粗细
                      colour = "gray50"  # 控制枝条颜色
) +
  # modify legend location using ggplot2
  theme(legend.position = c(0.92,0.92), legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.key.size= unit(0.2, "cm"),legend.key = element_blank())
p2 <- rotate_tree(p, 10)
p3 <- p2 +
  geom_fruit(
    data     = classification_df,
    geom     = geom_tile,
    mapping  = aes(y= label, x = 1, fill = classify),
    offset   = 2,
    pwidth   = 1.05,
    position = "identity",
    color    = NA,
    width = 6
  ) +
  scale_fill_manual(values = c("Bryidae" = "#5fbea1", 
                               "Dicranidae" = "#4DAF4A",
                               "Other mosses" = "#5c6e83"
                               
  ))

ggsave(paste0(CHARACTER,"_ARD_Pie.pdf"), p3, width = 4, height = 6, limitsize = FALSE)

#######################################################################################################  
classification_df <- read.csv("liverworts_class.csv",stringsAsFactors = FALSE)
CHARACTER <- "Niche state"
custom_colors <- c("Tropical and Subtropical Forests" = "#9C84C7", "Temperate Forests" = "#F28A70", 
                    "Savannas Deserts and Shrublands" = "#3A73C1")

NUM_STATES <- 3
STATE_LABELS <- c("0" = "Tropical and Subtropical Forests", "1" = "Temperate Forests",
                 "2" = "Savannas Deserts and Shrublands")

tree_file <- paste0(CHARACTER,"_ER_output/ER.tree")
ase <- processAncStates(tree_file, state_labels = STATE_LABELS)
p <- plotAncStatesPie(t = ase,
                      pie_colors = custom_colors,
                      layout= "fan",
                      #layout= "circular",
                      state_transparency = 1.0,
                      open.angle  = 10,
                      ladderize = FALSE,
                      tip_labels=FALSE,
                      node_pie_size = 0.4,
                      tip_pie_size = 0.1,
                      size = 0.1, # 控制枝条粗细
                      colour = "gray50"  # 控制枝条颜色
) +
  # modify legend location using ggplot2
  theme(legend.position = c(0.92,0.92), legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.key.size= unit(0.2, "cm"),legend.key = element_blank())
p2 <- rotate_tree(p, 10)
p3 <- p2 +
  geom_fruit(
    data     = classification_df,
    geom     = geom_tile,
    mapping  = aes(y= label, x = 1, fill = classify),
    offset   = 2,
    pwidth   = 1.05,
    position = "identity",
    color    = NA,
    width = 6
  ) +
  scale_fill_manual(values = c("Jungermanniidae" = "#f8951e",
                               "Other liverworts" = "#cc0000"
  ))
ggsave(paste0(CHARACTER,"_ER_Pie.pdf"), p3, width = 4, height = 6, limitsize = FALSE)

#######################################################################################################  
classification_df <- read.csv("mosses_class.csv",stringsAsFactors = FALSE)
CHARACTER <- "Niche state"
custom_colors <- c("Tropical and Subtropical Forests" = "#9C84C7", "Temperate Forests" = "#F28A70", 
                   "Savannas Deserts and Shrublands" = "#3A73C1")
NUM_STATES <- 3
STATE_LABELS <- c("0" = "Tropical and Subtropical Forests", "1" = "Temperate Forests",
                 "2" = "Savannas Deserts and Shrublands")

tree_file <- paste0(CHARACTER,"_ARD_output/ARD.tree")
ase <- processAncStates(tree_file, state_labels = STATE_LABELS)
p <- plotAncStatesPie(t = ase,
                      pie_colors = custom_colors,
                      layout= "fan",
                      #layout= "circular",
                      state_transparency = 1.0,
                      open.angle  = 10,
                      ladderize = FALSE,
                      tip_labels=FALSE,
                      node_pie_size = 0.4,
                      tip_pie_size = 0.1,
                      size = 0.1, # 控制枝条粗细
                      colour = "gray50"  # 控制枝条颜色
) +
  # modify legend location using ggplot2
  theme(legend.position = c(0.92,0.92), legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.key.size= unit(0.2, "cm"),legend.key = element_blank())
p2 <- rotate_tree(p, 10)
p3 <- p2 +
  geom_fruit(
    data     = classification_df,
    geom     = geom_tile,
    mapping  = aes(y= label, x = 1, fill = classify),
    offset   = 2,
    pwidth   = 1.05,
    position = "identity",
    color    = NA,
    width = 6
  ) +
  scale_fill_manual(values = c("Bryidae" = "#5fbea1", 
                               "Dicranidae" = "#4DAF4A",
                               "Other mosses" = "#5c6e83"
                               
  ))
ggsave(paste0(CHARACTER,"_ARD_Pie.pdf"), p3, width = 4, height = 6, limitsize = FALSE)

