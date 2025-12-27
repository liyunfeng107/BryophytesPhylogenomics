setwd("G:\\基因组\\tree\\species_tree\\20251001\\geneshoping\\revbayes")
library(RevGadgets)
library(phytools)
library(tibble)
library(ggtree)
library(treeio)
library(ggplot2)
library(RColorBrewer)
source(file = "plot_branch_rates_tree2.R")

# load the files:
my_tree_file = "bry_og_time_genus_rename.tre"
my_branch_rates_file = "monodoreae3_BDSE_10RateCat_4shift_rates.log"

# set the colors:
Colors <- colorRampPalette(rev(c('darkred',brewer.pal(n = 8, name = "Spectral"),'darkblue')))(100)
#plot_branch_rates_tree2(tree_file=my_tree_file,  branch_rates_file=my_branch_rates_file, parameter_name = "lambda", trans = "identity", colors = Colors) + geom_tiplab(color = "black", size = 2) + scale_color_gradientn("Speciation rate", colors = Colors, trans = "identity") + theme(legend.position=c(0.2,0.85))+ scale_x_continuous(limits = c(0,35))
plot_branch_rates_tree2(tree_file=my_tree_file,  branch_rates_file=my_branch_rates_file, parameter_name = "lambda", trans = "identity", colors = Colors) + scale_color_gradientn("Speciation rate", colors = Colors, trans = "identity") + theme(legend.position=c(0.2,0.85))
ggsave("RevBayes_BDSE_nuclear_Speciation_rate.pdf", width=15, height=15, units="cm")

plot_branch_rates_tree2(tree_file=my_tree_file,  branch_rates_file=my_branch_rates_file, parameter_name = "net_div", trans = "identity", colors = Colors) + scale_color_gradientn("Net speciation rate", colors = Colors, trans = "identity") + theme(legend.position=c(0.2,0.85))
ggsave("RevBayes_BDSE_nuclear_net_Speciation_rate.pdf", width=15, height=15, units="cm")

#获取lambda
lambda_raw = plot_branch_rates_tree2(tree_file=my_tree_file, branch_rates_file=my_branch_rates_file)
lambda <- lambda_raw$data
lambda_clean <- lambda[, c("label", "rates")]
data_clean <- lambda_clean[!is.na(lambda_clean$label), ]
write.csv(data_clean, "Revbayes_rates_nuclear_net.csv")

#获取net
net_raw = plot_branch_rates_tree2(tree_file=my_tree_file, branch_rates_file=my_branch_rates_file, parameter_name = "net_div")
net <- net_raw$data
net_clean <- net[, c("label", "rates")]
data_clean <- net_clean[!is.na(net_clean$label), ]
write.csv(data_clean, "Revbayes_rates_nuclear_net.csv")
