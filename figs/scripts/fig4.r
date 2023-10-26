############################################################
# For rodent genomes
# Figure 4
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(tidyverse)
#library(ggplot2)
library(cowplot)
library(ggbeeswarm)
#library(dplyr)
library(RColorBrewer)
library(GenomicRanges)
library(here)
source(here("figs", "scripts", "lib", "design.r"))

############################################################
cat("----------\n")

window_size = 10
# Window size in kb

marker_window_size = 5
# Marker window size in Mb

marker_window_size_kb = marker_window_size * 1000
# Marker window size in kb

read_data = T
# Whether or not to re-read the input data

save_fig = F
# Whether or not to save the fig

datadir = here("summary-data", "02-genomic-windows")

#infile = here(datadir, paste(window_size, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts-tt.csv.gz", sep=""))

if(marker_window_size < 1){
  stop()
}else{
  marker_file = here(datadir, paste(window_size, "kb-", marker_window_size, "mb-recomb-dists-st.csv", sep=""))
  #marker_file = paste(datadir, window_size, "kb-", marker_window_size, "mb-recomb-dists.csv", sep="")
}
# Input options
######################

if(read_data){
  cat(as.character(Sys.time()), " | Reading marker window data: ", marker_file, "\n")
  marker_windows = read.csv(marker_file, header=T)
}

marker_windows = marker_windows %>% filter(recomb.rate < 2e-6)

cat(as.character(Sys.time()), " | Fig4: Rendering panel A\n")
fig_4a = ggplot(marker_windows, aes(x=recomb.rate, y=first.last.wrf)) +
  geom_point(size=2, color=corecol(numcol=1, pal="wilke", offset=2), alpha=0.5) +
  geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
  ylab("wRF between\nfirst and last tree") +
  xlab("") +
  bartheme() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=10),
        plot.margin=unit(c(0.5,0.5,0,0), units="cm"))
print(fig_4a)

cat(as.character(Sys.time()), " | Fig4: Rendering panel B\n")
fig_4b = ggplot(marker_windows, aes(x=recomb.rate, y=wrf.decay)) +
  geom_point(size=2, color=corecol(numcol=1, pal="wilke", offset=2), alpha=0.5) +
  geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
  ylab("Decay rate\nfrom first tree") +
  xlab("Recombination rate") +
  bartheme() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=10),
        plot.margin=unit(c(0.5,0.5,0,0), units="cm"))
print(fig_4b)

cat(as.character(Sys.time()), " | Fig4: Rendering panel C\n")
fig_4c = ggplot(marker_windows, aes(x=recomb.rate, y=wrf.st.avg)) +
  geom_point(size=2, color=corecol(numcol=1, pal="wilke", offset=2), alpha=0.5) +
  geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
  ylab("Mean wRF\nto species tree") +
  xlab("") +
  bartheme() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=10),
        plot.margin=unit(c(0.5,0,0,0), units="cm"))
print(fig_4c)

######################

cat(as.character(Sys.time()), " | Fig4: Combining panels A and B and C\n")
fig4 = plot_grid(fig_4a, fig_4b, fig_4c, ncol=3, labels=c("A", "B", "C"), label_size=16, align="vh")

######################

if(save_fig){
  fig_file = here("figs", "fig4.png")
  cat(as.character(Sys.time()), " | Fig4: Saving figure:", fig_file, "\n")
  ggsave(filename=fig_file, fig4, width=7.5, height=2.5, units="in")
}

######################

# hotspots_file = here("summary-data", "02-genomic-windows", "recombination-markers", "smagulova-hotspots", "Smagulova_2010-hotspots-mm10.bed")
# hotspots = read_tsv(hotspots_file, col_names=c("chr", "start", "end"))
# 
# marker_windows = marker_windows %>% separate(window, into=c("chr", "start", "end"), remove=F)
# 
# mw_range = makeGRangesFromDataFrame(marker_windows)
# # Create range for marker windows
# 
# hs_range = makeGRangesFromDataFrame(hotspots)
# hs_overlaps = countOverlaps(mw_range, hs_range)
# marker_windows$hotspots = hs_overlaps
# # Create range for hotspots and get overlaps with marker windows
# 
# fig_4d = ggplot(marker_windows, aes(x=recomb.rate, y=hotspots)) +
#   geom_point(size=2, color=corecol(numcol=1, pal="wilke", offset=2), alpha=0.5) +
#   geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
#   ylab("# hotspots") +
#   xlab("Recombination rate") +
#   bartheme() +
#   theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
#         axis.text.y=element_text(size=8),
#         axis.title.x=element_text(size=14),
#         axis.title.y=element_text(size=10),
#         plot.margin=unit(c(0.5,0,0,0), units="cm"))
# print(fig_4d)

# marker_windows_high_wrf = subset(marker_windows, first.last.wrf > 0.75)
# 
# 
# fig_4a_high = ggplot(marker_windows_high_wrf, aes(x=marker.slopes, y=first.last.wrf)) +
#   geom_point(size=2, color=corecol(numcol=1, pal="wilke", offset=2), alpha=0.5) +
#   geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
#   ylab("wRF between\nfirst and last tree") +
#   xlab("Recombination rate") +
#   bartheme() +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# #print(fig_4a_high)
# 
# 
# fig_4b_high = ggplot(marker_windows_high_wrf, aes(x=marker.slopes, y=decay.to.first)) +
#   geom_point(size=2, color=corecol(numcol=1, pal="wilke", offset=2), alpha=0.5) +
#   geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
#   ylab("Decay rate\nfrom first tree") +
#   xlab("Recombination rate") +
#   bartheme() +
#   theme(axis.text.x=element_text(angle=45, hjust=1))
# #print(fig_4b_high)
