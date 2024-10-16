############################################################
# For penn genomes, 06.21
# Figure S window sites
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(tidyverse)
library(cowplot)
library(ggbeeswarm)
library(here)
source(here("figs", "scripts", "lib", "design.r"))

############################################################
cat("----------\n")

window_size = 10
# Window size in kb

marker_window_size = 5
# Marker window size in Mb

au_flag = FALSE
# Set to filter out windows that don't pass the AU test

skip_one = FALSE
# Set to only do one test chromosome

read_data = T
# Whether or not to re-read the input data

save_fig = T
# Whether or not to save the fig

aln_stats_file = here("summary-data", "02-genomic-windows", "phykit-out.csv")

# Input options
######################

if(read_data){
  cat(as.character(Sys.time()), " | Reading window data: ", aln_stats_file, "\n")
  in_data = read.csv(aln_stats_file, header=T)
}

cat(as.character(Sys.time()), " | Generating figure\n")

var_sites = select(in_data, var.sites)
names(var_sites)[1] = "num.sites"
var_sites$label = "Variable sites"

inf_sites = select(in_data, informative.sites)
names(inf_sites)[1] = "num.sites"
inf_sites$label = "Informative sites"

sites = rbind(var_sites, inf_sites)
sites$label = factor(sites$label, levels=c("Variable sites", "Informative sites"))

fig = ggplot(sites, aes(x=label, y=num.sites, group=label)) +
  geom_boxplot(alpha=0.7, width=0.5) +
  xlab("") +
  ylab("# of sites\nper 10kb window") +
  bartheme() +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        axis.text.x=element_text(angle=25, hjust=1))
print(fig)

if(save_fig){
  figfile = here("figs", "supp", "figS4.png")
  cat(as.character(Sys.time()), " | Saving figure: ", figfile, "\n")
  ggsave(filename=figfile, fig, width=3, height=4, units="in")
}

# fig_a = ggplot(in_data, aes(x=1, y=aln.length)) +
#   #geom_quasirandom(size=2, alpha=0.7, width=0.25, color="#d3d3d3") +
#   geom_boxplot(outlier.shape=NA, alpha=0.7, width=0.5) +
#   bartheme()
# print(fig_a)
# 
# fig_b = ggplot(in_data, aes(x=1, y=var.sites)) +
#   #geom_quasirandom(size=2, alpha=0.7, width=0.25, color="#d3d3d3") +
#   geom_boxplot(outlier.shape=NA, alpha=0.7, width=0.5) +
#   bartheme()
# print(fig_b)
# 
# fig_c = ggplot(in_data, aes(x=1, y=informative.sites)) +
#   #geom_quasirandom(size=2, alpha=0.7, width=0.25, color="#d3d3d3") +
#   geom_boxplot(outlier.shape=NA, alpha=0.7, width=0.5) +
#   bartheme()
# print(fig_c)
# 










