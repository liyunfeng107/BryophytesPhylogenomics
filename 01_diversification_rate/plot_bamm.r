setwd("")
library(BAMMtools)
##########################################################
# plot bamm results
##########################################################
tree <- read.tree('A.tre')

###1. Assessing MCMC convergence
mcmcout<-read.csv("mcmc_out.txt",header=T) #读取mcmcout
plot(mcmcout$logLik ~ mcmcout$generation) #plot the log-likelihood trace of your MCMC output file.
##PDF: This can give you a ballpark idea of whether your run has converged
burnstart <- floor(0.1 * nrow(mcmcout)) #discard some as burnin (10%)
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
#ESS检测
library(coda) #Use coda library to check the effective sample sizes of the log-likelihood and the number of shift events present in each sample
effectiveSize(postburn$N_shifts) #at least 200
effectiveSize(postburn$logLik) #at least 200

###2. The bammdata object
library(BAMMtools)
tre <- read.tree("dated_Vitaceae_505_1-e12_dropout.tre")
tre <- ladderize(tre, right = FALSE)
edata<-getEventData(tre,eventdata="event_data.txt",burnin=0.1) #读取edata

##3. Overall best shift configuration
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=1, threshold=5)

pdffn = "MAP_best_config.pdf"
pdf(pdffn, width=6, height=8)
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=1, threshold=5)
q <- plot(best, breaksmethod='jenks', lwd=1)
addBAMMshifts(best, cex=1)
addBAMMlegend(q, location='topleft')
dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

pdffn = "MAP_best_config_tip.pdf"
pdf(pdffn, width=20, height=110)
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=1, threshold=5)
q <- plot(best, labels = TRUE, breaksmethod='jenks', lwd=1.25, cex.lab = 1)
addBAMMshifts(best, cex=6)
addBAMMlegend(q, location='topleft')
dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

##4. Rate-through-time analysis

pdffn = "RTTplot_Bryophytes.pdf"
pdf(pdffn, width=8, height=20)
par(mfrow=c(3,1))
plotRateThroughTime(edata, ratetype="speciation", intervalCol="blue", avgCol="blue", ylim=c(0,0.05), cex.axis=2)
plotRateThroughTime(edata, ratetype="extinction", intervalCol="red", avgCol="red", ylim=c(0,0.01), cex.axis=2)
plotRateThroughTime(edata, ratetype="netdiv", intervalCol="green", avgCol="green", ylim=c(0,0.04), cex.axis=2)
dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

moss <- getMRCA(tre, tip = c("Takakiales_Takakiaceae_Takakia", "Sphagnales_Ambuchananiaceae_Eosphagnum"))
pdffn = "RTTplot_Bryophytes_moss.pdf"
pdf(pdffn, width=8, height=20)
par(mfrow=c(3,1))
plotRateThroughTime(edata, node=moss, ratetype="speciation", intervalCol="blue", avgCol="blue", ylim=c(0,0.06), cex.axis=2)
plotRateThroughTime(edata, node=moss, ratetype="extinction", intervalCol="red", avgCol="red", ylim=c(0,0.01), cex.axis=2)
plotRateThroughTime(edata, node=moss, ratetype="netdiv", intervalCol="green", avgCol="green", ylim=c(0,0.06), cex.axis=2)
dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

liverwort <- getMRCA(tre, tip = c("Calobryales_Haplomitriaceae_Haplomitrium", "Lunulariales_Lunulariaceae_Lunularia"))
pdffn = "RTTplot_Bryophytes_liver.pdf"
pdf(pdffn, width=8, height=20)
par(mfrow=c(3,1))
plotRateThroughTime(edata, node=liverwort, ratetype="speciation", intervalCol="blue", avgCol="blue", ylim=c(0,0.04), cex.axis=2)
plotRateThroughTime(edata, node=liverwort, ratetype="extinction", intervalCol="red", avgCol="red", ylim=c(0,0.01), cex.axis=2)
plotRateThroughTime(edata, node=liverwort, ratetype="netdiv", intervalCol="green", avgCol="green", ylim=c(0,0.04), cex.axis=2)
dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

