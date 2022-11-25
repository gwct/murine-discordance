############################################################
# For penn genomes, 06.21
# Figure 4
# Gregg Thomas
############################################################

library(ggplot2)
library(cowplot)
library(ggbeeswarm)
library(dplyr)
library(RColorBrewer)

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

############################################################
cat("----------\n")

window_size = 10
# Window size in kb

marker_window_size = 1
# Marker window size in Mb

marker_window_size_kb = marker_window_size * 1000
# Marker window size in kb

au_flag = FALSE
# Set to filter out windows that don't pass the AU test

skip_one = FALSE
# Set to only do one test chromosome

read_data = T
# Whether or not to re-read the input data

save_fig = F
# Whether or not to save the fig

datadir = "C:/Users/Gregg/Box Sync/rodents/penn/paper/data/"
tree_file = paste(datadir, window_size, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts-tt.csv", sep="")

if(marker_window_size < 1){
  stop()
}else{
  marker_file = paste(datadir, window_size, "kb-", marker_window_size, "mb-recomb-dists.csv", sep="")
}
# Input options
######################

if(read_data){
  cat(as.character(Sys.time()), " | Reading marker window data: ", tree_file, "\n")
  marker_windows = read.csv(marker_file, header=T)
}

cat(as.character(Sys.time()), " | Fig4: Rendering panel A\n")
fig_4a = ggplot(marker_windows, aes(x=recomb.rate, y=first.last.wrf)) +
  geom_point(size=2, color=corecol(numcol=1, pal="wilke", offset=2), alpha=0.5) +
  geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
  ylab("wRF between\nfirst and last tree") +
  xlab("Recombination rate") +
  bartheme() +
  theme(axis.text.x=element_text(angle=45, hjust=1))
print(fig_4a)

cat(as.character(Sys.time()), " | Fig4: Rendering panel B\n")
fig_4b = ggplot(marker_windows, aes(x=recomb.rate, y=wrf.decay)) +
  geom_point(size=2, color=corecol(numcol=1, pal="wilke", offset=2), alpha=0.5) +
  geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
  ylab("Decay rate\nfrom first tree") +
  xlab("Recombination rate") +
  bartheme() +
  theme(axis.text.x=element_text(angle=45, hjust=1))
print(fig_4b)

######################

cat(as.character(Sys.time()), " | Fig4: Combining panels A and B\n")
fig4 = plot_grid(fig_4a, fig_4b, ncol=2, labels=c("A", "B"), label_size=16, align="vh")

######################

if(save_fig){
  fig4file = "../figs/fig4.png"
  cat(as.character(Sys.time()), " | Fig4: Saving figure:", fig4file, "\n")
  ggsave(filename=fig4file, fig4, width=7.5, height=3.5, units="in")
}

######################

marker_windows_high_wrf = subset(marker_windows, first.last.wrf > 0.75)


fig_4a_high = ggplot(marker_windows_high_wrf, aes(x=marker.slopes, y=first.last.wrf)) +
  geom_point(size=2, color=corecol(numcol=1, pal="wilke", offset=2), alpha=0.5) +
  geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
  ylab("wRF between\nfirst and last tree") +
  xlab("Recombination rate") +
  bartheme() +
  theme(axis.text.x=element_text(angle=45, hjust=1))
#print(fig_4a_high)


fig_4b_high = ggplot(marker_windows_high_wrf, aes(x=marker.slopes, y=decay.to.first)) +
  geom_point(size=2, color=corecol(numcol=1, pal="wilke", offset=2), alpha=0.5) +
  geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
  ylab("Decay rate\nfrom first tree") +
  xlab("Recombination rate") +
  bartheme() +
  theme(axis.text.x=element_text(angle=45, hjust=1))
#print(fig_4b_high)
