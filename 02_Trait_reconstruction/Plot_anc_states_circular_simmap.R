library(plotrix)
library(phytools)
library(RevGadgets)
setwd("")

#######################################################################################################
CHARACTER <- "Peristome type"
character_file <- paste0(CHARACTER,"_HRM2_ARD_output/HRM2_ARD_stoch_map_character.tree")
simmap = read.simmap(file=character_file, format="phylip")

custom_colors <- c("0" = "#3A73C1", "1" = "#F28A70", "2"= "#9C84C7",
                   "3" = "#c5daf2","4" = "#fccbbc","5" = "#c3eff3")
states <- unique(unlist(lapply(simmap$maps, names)))
cols <- custom_colors[states]

pdf( paste0(CHARACTER,"_HRM2_ARD_simmap_circular.pdf") )
p <- plotSimmap(simmap, cols, fsize=0.1, lwd=1, split.vertical=TRUE,type = "fan", 
                offset = 20,part = (360 - 10)/360)
#pdf( paste0(CHARACTER,"_ER_simmap.pdf") )
#Ntip <- length(simmap$tip.label)
#p <- plotSimmap(simmap, cols, fsize=0.1, lwd=1, split.vertical=TRUE, offset = 20, ylim  = c(Ntip + 1, 0))
# add legend
leg_all <- c("0" = "Slow gain of Haplolepideous", "1" = "Slow gain of Diplolepideous",
             "2" = "low gain of Absent", "3" = "Fast gain of Haplolepideous",
             "4" = "Fast gain of Diplolepideous", "5" = "Fast gain of Absent")
states <- unique(unlist(lapply(simmap$maps, names)))
leg    <- leg_all[states]
legend("topleft",                            # 从 topright 换成 topleft
       legend    = leg,                        # 标签文字
       pch       = 21,                         # 圆点
       pt.bg     = custom_colors[names(leg)],  # 填充色
       pt.cex    = 1.5,                        # 点大小
       bty       = "n",                        # 无边框
       cex       = 1.15)                       # 字体大小
dev.off()

#######################################################################################################
CHARACTER <- "Fungal symbiosis status"
character_file <- paste0(CHARACTER,"_HRM2_ER_output/HRM2_ER_stoch_map_character.tree")
simmap = read.simmap(file=character_file, format="phylip")

custom_colors <- c("0" = "#3A73C1", "1" = "#F28A70",  
                   "2" = "#c5daf2", "3" = "#fccbbc" )
states <- unique(unlist(lapply(simmap$maps, names)))
cols <- custom_colors[states]

pdf( paste0(CHARACTER,"_HRM2_ER_simmap_circular.pdf") )
p <- plotSimmap(simmap, cols, fsize=0.1, lwd=1, split.vertical=TRUE,type = "fan", 
                offset = 20,part = (360 - 10)/360)
# add legend
leg_all <- c("0" = "Slow gain of Symbiotic", "1" = "Slow gain of None",
             "2" = "Fast gain of Symbiotic", "3" = "Fast gain of None")
states <- unique(unlist(lapply(simmap$maps, names)))
leg    <- leg_all[states]
legend("topleft",                            # 从 topright 换成 topleft
       legend    = leg,                        # 标签文字
       pch       = 21,                         # 圆点
       pt.bg     = custom_colors[names(leg)],  # 填充色
       pt.cex    = 1.5,                        # 点大小
       bty       = "n",                        # 无边框
       cex       = 1.15)                       # 字体大小
dev.off()

#######################################################################################################
CHARACTER <- "Sexual systems"
character_file <- paste0(CHARACTER,"_ARD_output/ARD_stoch_map_character.tree")
simmap = read.simmap(file=character_file, format="phylip")

custom_colors <- c("0" = "#F28A70", "1" = "#3A73C1" )
states <- unique(unlist(lapply(simmap$maps, names)))
cols <- custom_colors[states]

pdf( paste0(CHARACTER,"_ARD_simmap_circular1.pdf") )
p <- plotSimmap(simmap, cols, fsize=0.1, lwd=1, split.vertical=TRUE,type = "fan", 
                offset = 20,part = (360 - 10)/360)
# add legend
leg_all <- c("0" = "Bisexual", "1" = "Unisexual")

states <- unique(unlist(lapply(simmap$maps, names)))
leg    <- leg_all[states]
legend("topleft",                            # 从 topright 换成 topleft
       legend    = leg,                        # 标签文字
       pch       = 21,                         # 圆点
       pt.bg     = custom_colors[names(leg)],  # 填充色
       pt.cex    = 1.5,                        # 点大小
       bty       = "n",                        # 无边框
       cex       = 1.15)                       # 字体大小
dev.off()

#######################################################################################################
CHARACTER <- "Capsule position"
character_file <- paste0(CHARACTER,"_HRM4_ARD_output/HRM4_ARD_stoch_map_character.tree")
simmap = read.simmap(file=character_file, format="phylip")

custom_colors <- c("0" = "#F28A70", "1" = "#3A73C1", "2" = "#9C84C7", 
                   "3" = "#fccbbc", "4" = "#c5daf2","5" = "#c3eff3", 
                   "6" = "#fccbbc", "7" = "#c5daf2","8" = "#c3eff3", 
                   "9" = "#fccbbc", "10" = "#c5daf2", "11" = "#c3eff3")

states <- unique(unlist(lapply(simmap$maps, names)))
cols <- custom_colors[states]

pdf( paste0(CHARACTER,"_HRM4_ARD_simmap_circular.pdf") )
p <- plotSimmap(simmap, cols, fsize=0.1, lwd=1, split.vertical=TRUE,type = "fan", 
                offset = 20,part = (360 - 10)/360)
# add legend
leg_all <- c("0" = "VSlow gain of Pleurocarpous", "1" = "VSlow gain of Acrocarpous",
             "2" = "VSlow gain of Cladocarpous", "3" = "Slow gain of Pleurocarpous",
             "4" = "Slow gain of Acrocarpous", "5" = "Slow gain of Cladocarpous",
             "6" = "Fast gain of Pleurocarpous","7" = "Fast gain of Acrocarpous", 
             "8" = "Fast gain of Cladocarpous", "9" = "VFast gain of Pleurocarpous",
             "10" = "VFast gain of Acrocarpous", "11" = "VFast gain of Cladocarpous")

states <- unique(unlist(lapply(simmap$maps, names)))
leg    <- leg_all[states]
legend("topleft",                            # 从 topright 换成 topleft
       legend    = leg,                        # 标签文字
       pch       = 21,                         # 圆点
       pt.bg     = custom_colors[names(leg)],  # 填充色
       pt.cex    = 1.5,                        # 点大小
       bty       = "n",                        # 无边框
       cex       = 1.15)                       # 字体大小
dev.off()

#######################################################################################################
CHARACTER <- "Habitats"
character_file <- paste0(CHARACTER,"_liverworts_ARD_output/ARD_stoch_map_character.tree")
simmap = read.simmap(file=character_file, format="phylip")

custom_colors <- c("0" = "#3A73C1", "1" = "#F28A70",  
                   "2" = "#9C84C7")
states <- unique(unlist(lapply(simmap$maps, names)))
cols <- custom_colors[states]

pdf( paste0(CHARACTER,"_liverworts_ARD_simmap_circular.pdf") )
p <- plotSimmap(simmap, cols, fsize=0.1, lwd=1, split.vertical=TRUE,type = "fan", 
                offset = 20,part = (360 - 10)/360)
# add legend
leg_all <- c("0" = "Terrestrial", "1" = "Epiphytic",
             "2" = "Aquatic")
states <- unique(unlist(lapply(simmap$maps, names)))
leg    <- leg_all[states]
legend("topleft",                            # 从 topright 换成 topleft
       legend    = leg,                        # 标签文字
       pch       = 21,                         # 圆点
       pt.bg     = custom_colors[names(leg)],  # 填充色
       pt.cex    = 1.5,                        # 点大小
       bty       = "n",                        # 无边框
       cex       = 1.15)                       # 字体大小
dev.off()

#######################################################################################################
CHARACTER <- "Habitats"
character_file <- paste0(CHARACTER,"_mosses_ARD_output/ARD_stoch_map_character.tree")
simmap = read.simmap(file=character_file, format="phylip")

custom_colors <- c("0" = "#3A73C1", "1" = "#F28A70",  
                   "2" = "#9C84C7")
states <- unique(unlist(lapply(simmap$maps, names)))
cols <- custom_colors[states]

pdf( paste0(CHARACTER,"_moss_ARD_simmap_circular.pdf") )
p <- plotSimmap(simmap, cols, fsize=0.1, lwd=1, split.vertical=TRUE,type = "fan", 
                offset = 20,part = (360 - 10)/360)
# add legend
leg_all <- c("0" = "Terrestrial", "1" = "Epiphytic",
             "2" = "Aquatic")
states <- unique(unlist(lapply(simmap$maps, names)))
leg    <- leg_all[states]
legend("topleft",                            # 从 topright 换成 topleft
       legend    = leg,                        # 标签文字
       pch       = 21,                         # 圆点
       pt.bg     = custom_colors[names(leg)],  # 填充色
       pt.cex    = 1.5,                        # 点大小
       bty       = "n",                        # 无边框
       cex       = 1.15)                       # 字体大小
dev.off()

#######################################################################################################
CHARACTER <- "Niche states"
character_file <- paste0(CHARACTER,"_liverworts_ER_output/ER_stoch_map_character.tree")
simmap = read.simmap(file=character_file, format="phylip")

custom_colors <- c("0" = "#9C84C7", "1" = "#F28A70", 
                   "2" = "#3A73C1")
states <- unique(unlist(lapply(simmap$maps, names)))
cols <- custom_colors[states]

pdf( paste0(CHARACTER,"_liverworts_ER_simmap_circular.pdf") )
p <- plotSimmap(simmap, cols, fsize=0.1, lwd=1, split.vertical=TRUE,type = "fan", 
                offset = 20,part = (360 - 10)/360)
# add legend
leg_all <- c("0" = "Tropical and Subtropical Forests", "1" = "Temperate Forests",
             "2" = "Savannas Deserts and Shrublands")
states <- unique(unlist(lapply(simmap$maps, names)))
leg    <- leg_all[states]
legend("topleft",                            # 从 topright 换成 topleft
       legend    = leg,                        # 标签文字
       pch       = 21,                         # 圆点
       pt.bg     = custom_colors[names(leg)],  # 填充色
       pt.cex    = 1.5,                        # 点大小
       bty       = "n",                        # 无边框
       cex       = 1.15)                       # 字体大小
dev.off()

#######################################################################################################
CHARACTER <- "Niche states"
character_file <- paste0(CHARACTER,"_mosses_ARD_output/ARD_stoch_map_character.tree")
simmap = read.simmap(file=character_file, format="phylip")

custom_colors <- c("0" = "#9C84C7", "1" = "#F28A70", 
                   "2" = "#3A73C1")
states <- unique(unlist(lapply(simmap$maps, names)))
cols <- custom_colors[states]

pdf( paste0(CHARACTER,"_moss_ARD_simmap_circular.pdf") )
p <- plotSimmap(simmap, cols, fsize=0.1, lwd=1, split.vertical=TRUE,type = "fan", 
                offset = 20,part = (360 - 10)/360)
# add legend
leg_all <- c("0" = "Tropical and Subtropical Forests", "1" = "Temperate Forests",
             "2" = "Savannas Deserts and Shrublands")
states <- unique(unlist(lapply(simmap$maps, names)))
leg    <- leg_all[states]
legend("topleft",                            # 从 topright 换成 topleft
       legend    = leg,                        # 标签文字
       pch       = 21,                         # 圆点
       pt.bg     = custom_colors[names(leg)],  # 填充色
       pt.cex    = 1.5,                        # 点大小
       bty       = "n",                        # 无边框
       cex       = 1.15)                       # 字体大小
dev.off()

