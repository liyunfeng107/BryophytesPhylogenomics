library(RevGadgets)
library(ggplot2)
setwd("")
# specify the output files
speciation_time_file <- "EBD_speciation_times.log"
speciation_rate_file <- "EBD_speciation_rates.log"
extinction_time_file <- "EBD_extinction_times.log"
extinction_rate_file <- "EBD_extinction_rates.log"

# read in and process rates
rates <- processDivRates(speciation_time_log = speciation_time_file,
                         speciation_rate_log = speciation_rate_file,
                         extinction_time_log = extinction_time_file,
                         extinction_rate_log = extinction_rate_file,
                         burnin = 0.25,
                         summary = "median")

# plot rates through time
p <- plotDivRates(rates = rates) +
  xlab("Millions of years ago") +
  ylab("Rate per million years")

ggsave("EBD_nuclear.pdf", p, width = 160, height = 160, units = "mm")

#############################################################################################
library(RevGadgets)
setwd("")
# specify the output files
speciation_time_file <- "EBD_speciation_times.log"
speciation_rate_file <- "EBD_speciation_rates.log"
extinction_time_file <- "EBD_extinction_times.log"
extinction_rate_file <- "EBD_extinction_rates.log"

# read in and process rates
rates <- processDivRates(speciation_time_log = speciation_time_file,
                         speciation_rate_log = speciation_rate_file,
                         extinction_time_log = extinction_time_file,
                         extinction_rate_log = extinction_rate_file,
                         burnin = 0.25,
                         summary = "median")

# plot rates through time
p <- plotDivRates(rates = rates) +
  xlab("Millions of years ago") +
  ylab("Rate per million years")

ggsave("EBD_locus_liver.pdf", p, width = 160, height = 160, units = "mm")

#############################################################################################
library(RevGadgets)
setwd("")
# specify the output files
speciation_time_file <- "speciation_times.log"
speciation_rate_file <- "speciation_rates.log"
extinction_time_file <- "extinction_times.log"
extinction_rate_file <- "extinction_rates.log"

# read in and process rates
rates <- processDivRates(speciation_time_log = speciation_time_file,
                         speciation_rate_log = speciation_rate_file,
                         extinction_time_log = extinction_time_file,
                         extinction_rate_log = extinction_rate_file,
                         burnin = 0.25,
                         summary = "median")

# plot rates through time
p <- plotDivRates(rates = rates) +
  xlab("Millions of years ago") +
  ylab("Rate per million years")

ggsave("EBD_locus_moss.pdf", p, width = 160, height = 160, units = "mm")
