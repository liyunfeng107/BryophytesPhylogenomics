library(phytools)
library(ape)
library(ggplot2)
setwd("")
######################################################################################
A<-read.tree("A.tre")
B<-read.tree("B.tre")
A$edge.length <- rep(1, nrow(A$edge))
B$edge.length <- rep(1, nrow(B$edge))
cophylotr <- cophylo(A,B,rorate = T)

#par(family = "Arial", cex = 5 / 12)
pdf("Plot_tree_conflict.pdf", width=7, height=7)
plot(cophylotr,link.type="curved", link.lwd=5, link.lty="solid",  
     lwd = 1, tip.lwd = 1, pts = FALSE, tip.lty = 1, fsize = 0.5,
     link.col=make.transparent("red",0.4))
# ==== 加节点支持率 ====
lsup <- cophylotr$trees[[1]]$node.label
rsup <- cophylotr$trees[[2]]$node.label
# 左树支持率
nodelabels.cophylo(lsup, which="left",
                   frame="none", cex=0.45, adj=c(1.2,-0.4))
# 右树支持率
nodelabels.cophylo(rsup, which="right",
                   frame="none", cex=0.45, adj=c(-0.2,-0.4))
dev.off()

