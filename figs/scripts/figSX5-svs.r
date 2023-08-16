############################################################
# For penn genomes, 06.22
# Gets overlaps/correlations between recombination rate in 5Mb windows
# and other genomic features
# Gregg Thomas
############################################################

library(tidyverse)
library(cowplot)
library(ggbeeswarm)
library(RColorBrewer)
library(GenomicRanges)

############################################################

revOrder <- function(df){
  rev_count = 0
  for(i in 1:nrow(df)){
    if(df[i,]$end < df[i,]$start){
      tmp_start = df[i,]$start
      tmp_end = df[i,]$end
      
      df[i,]$start = tmp_end
      df[i,]$end = tmp_start
      
      rev_count = rev_count + 1
      
    }
  }
  
  print(rev_count)
  return(df)
}

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

save_fig = T
# Whether or not to save the fig

source_pafr = T

infile = here("summary-data", "02-genomic-windows", paste0(window_size, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts.csv"))

aln_stats_file = here("summary-data", "02-genomic-windows", "phykit-out.csv")

syn_file = here("summary-data", "04-mouse-rat-svs", "mouse-rat-syntenic.bed")
inv_file = here("summary-data", "04-mouse-rat-svs", "mouse-rat-inversions.bed")
trans_file = here("summary-data", "04-mouse-rat-svs", "mouse-rat-translocations.bed")
inv_trans_file = here("summary-data", "04-mouse-rat-svs", "mouse-rat-inverted-translocations.bed")

gene_file = here("summary-data", "02-genomic-windows", "feature-beds", "mm10-genes.bed")
uce_file = here("summary-data", "02-genomic-windows", "feature-beds", "mm10-uce.bed")
ps_file = here("summary-data", "02-genomic-windows", "feature-beds", "mm10-ps-genes-all-gt.bed")
hs_file = here("summary-data", "02-genomic-windows", "recombination-markers", "smagulova-hotspots", "Smagulova_2010-hotspots-mm10.bed")

paf_file = here("summary-data", "04-mouse-rat-svs", "rn6-mm10.paf")

marker_file = here("summary-data", "02-genomic-windows", paste0(window_size, "kb-", marker_window_size, "mb-recomb-dists.csv"))

# Input options
######################

if(source_pafr){
  setwd("C:/bin/pafr/R/")
  files.sources = list.files()
  sapply(files.sources, source)
}

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

if(read_data){
  cat(as.character(Sys.time()), " | Reading window data: ", infile, "\n")
  all_windows = read_csv(infile, comment="#")
  names(all_windows) = make.names(names(all_windows))
  # Reading 10kb window tree file
  
  cat(as.character(Sys.time()), " | Reading aln stat data: ", aln_stats_file, "\n")
  aln_stats = read_csv(aln_stats_file, comment="#")
  # Reading 10kb window alignment stat file
  
  all_windows = merge(x=all_windows, y=aln_stats, by="window", all.x=TRUE)
  # Combining tree and alignment data
  
  cat(as.character(Sys.time()), " | Reading marker window data: ", marker_file, "\n")
  marker_windows = read_csv(marker_file, comment="#")
  marker_windows = marker_windows %>% 
    separate(window, c("chr", "coords"), sep=":", remove=F) %>% 
    separate(coords, c("start", "end"), sep="-", remove=T)
  marker_windows$start = as.numeric(marker_windows$start)
  marker_windows$end = as.numeric(marker_windows$end)
  # Reading marker window data

  # cat(as.character(Sys.time()), " | Reading PAF: ", tree_file, "\n")
  # paf = read_paf(paf_file)
  # paf = filter_secondary_alignments(paf)
  # paf = subset(paf, nchar(qname) <= 5 & nchar(tname) <= 5)
  # paf = subset(paf, qname != "chrM" & tname != "chrM")
  # paf = subset(paf, qname != "chrY" & tname != "chrY")
  # paf = subset(paf, alen > 1000 & mapq > 40)
  # Reading the PAF file if I want to parse the alignments myself...
  
  cat(as.character(Sys.time()), " | Reading bed: ", gene_file, "\n")
  genes = read.table(gene_file)
  names(genes) = c("chr", "start", "end", "id", "score", "no")
  genes$start = as.numeric(genes$start)
  genes$end = as.numeric(genes$end)
  genes$len = genes$end - genes$start
  #genes$chr = paste("chr", genes$chr, sep="")
  # Reading gene coords
  
  cat(as.character(Sys.time()), " | Reading bed: ", ps_file, "\n")
  ps_genes = read.table(ps_file)
  names(ps_genes) = c("chr", "start", "end", "id", "score", "no")
  ps_genes$start = as.numeric(ps_genes$start)
  ps_genes$end = as.numeric(ps_genes$end)
  ps_genes$len = ps_genes$end - ps_genes$start
  #ps_genes$chr = paste("chr", ps_genes$chr, sep="")
  # Reading positively selected gene coords
  
  cat(as.character(Sys.time()), " | Reading bed: ", uce_file, "\n")
  uces = read.table(uce_file)
  names(uces) = c("chr", "start", "end", "id", "score", "no", "tmp1", "tmp2", "tmp3")
  uces$start = as.numeric(uces$start)
  uces$end = as.numeric(uces$end)
  uces$len = uces$end - uces$start
  # Reading UCE coords
  
  cat(as.character(Sys.time()), " | Reading bed: ", hs_file, "\n")
  hotspots = read.table(hs_file)
  names(hotspots) = c("chr", "start", "end")
  hotspots$start = as.numeric(hotspots$start)
  hotspots$end = as.numeric(hotspots$end)
  hotspots$len = hotspots$end - hotspots$start
  # Reading hotspot coords
  
  cat(as.character(Sys.time()), " | Reading bed: ", syn_file, "\n")
  mr_syn = read.table(syn_file)
  names(mr_syn) = c("chr", "start", "end", "id", "parent.id", "type", "copy.change")
  mr_syn = subset(mr_syn, type!="HDR")
  mr_syn$start = as.numeric(mr_syn$start)
  mr_syn$end = as.numeric(mr_syn$end)
  mr_syn = revOrder(mr_syn)
  mr_syn$len = mr_syn$end - mr_syn$start
  # Reading syntenic alignments
  
  cat(as.character(Sys.time()), " | Reading bed: ", inv_file, "\n")
  mr_inv = read.table(inv_file)
  names(mr_inv) = c("chr", "start", "end", "id", "parent.id", "type", "copy.change")
  mr_inv = subset(mr_inv, type!="HDR")
  mr_inv$start = as.numeric(mr_inv$start)
  mr_inv$end = as.numeric(mr_inv$end)
  mr_inv = revOrder(mr_inv)
  mr_inv$len = mr_inv$end - mr_inv$start
  # Reading inverted alignments
  
  cat(as.character(Sys.time()), " | Reading bed: ", trans_file, "\n")
  mr_trans = read.table(trans_file)
  names(mr_trans) = c("chr", "start", "end", "id", "parent.id", "type", "copy.change")
  mr_trans = subset(mr_trans, type!="HDR")
  mr_trans$start = as.numeric(mr_trans$start)
  mr_trans$end = as.numeric(mr_trans$end)
  mr_trans = revOrder(mr_trans)
  mr_trans$len = mr_trans$end - mr_trans$start
  # Reading translocated alignments
  
  cat(as.character(Sys.time()), " | Reading bed: ", inv_trans_file, "\n")
  mr_inv_trans = read.table(inv_trans_file)
  names(mr_inv_trans) = c("chr", "start", "end", "id", "parent.id", "type", "copy.change")
  mr_inv_trans = subset(mr_inv_trans, type!="HDR")
  mr_inv_trans$start = as.numeric(mr_inv_trans$start)
  mr_inv_trans$end = as.numeric(mr_inv_trans$end)
  mr_inv_trans = revOrder(mr_inv_trans)
  mr_inv_trans$len = mr_inv_trans$end - mr_inv_trans$start
  # Reading inverted translocated alignments
}
# Reading input data
######################

marker_windows_high_wrf = subset(marker_windows, first.last.wrf > 0.75)
marker_windows$high.wrf = "N"
#marker_windows[marker_windows$first.last.wrf>0.75,]$high.wrf = "Y"
# Parsing high wrf windows
######################

mw_range = makeGRangesFromDataFrame(marker_windows)
# Create range for marker windows

gene_range = makeGRangesFromDataFrame(genes)
gene_overlaps = countOverlaps(mw_range, gene_range)
marker_windows$genes = gene_overlaps
# Create range for genes and get overlaps with marker windows

ps_gene_range = makeGRangesFromDataFrame(ps_genes)
ps_gene_overlaps = countOverlaps(mw_range, ps_gene_range)
marker_windows$ps.genes = ps_gene_overlaps
# Create range for positively selected genes and get overlaps with marker windows

uces_range = makeGRangesFromDataFrame(uces)
uces_overlaps = countOverlaps(mw_range, uces_range)
marker_windows$uces = uces_overlaps
# Create range for UCEs and get overlaps with marker windows

hs_range = makeGRangesFromDataFrame(hotspots)
hs_overlaps = countOverlaps(mw_range, hs_range)
marker_windows$hotspots = hs_overlaps
# Create range for hotspots and get overlaps with marker windows

syn_range = makeGRangesFromDataFrame(mr_syn)
syn_overlaps = countOverlaps(mw_range, syn_range)
marker_windows$syntenic.alns = syn_overlaps
# Create range for syntenic alignments and get overlaps with marker windows

inv_range = makeGRangesFromDataFrame(mr_inv)
inv_overlaps = countOverlaps(mw_range, inv_range)
marker_windows$inversion.alns = inv_overlaps
# Create range for inverted alignments and get overlaps with marker windows

trans_range = makeGRangesFromDataFrame(mr_trans)
trans_overlaps = countOverlaps(mw_range, trans_range)
marker_windows$translocation.alns = trans_overlaps
# Create range for translocated alignments and get overlaps with marker windows

inv_trans_range = makeGRangesFromDataFrame(mr_inv_trans)
inv_trans_overlaps = countOverlaps(mw_range, inv_trans_range)
marker_windows$inverted.translocation.alns = inv_trans_overlaps
# Create range for inverted translocated alignments and get overlaps with marker windows

# Getting overlaps
######################

marker_window_repeats = all_windows %>% 
  group_by(marker.window) %>% 
  summarize(avg.repeat=mean(perc.repeat), 
            avg.perc.inf.sites=mean(perc.informative.sites, na.rm=T), 
            avg.perc.var.sites=mean(perc.var.sites, na.rm=T))

names(marker_window_repeats)[1] = "window"

marker_windows = merge(x=marker_windows, y=marker_window_repeats, by="window", all.x=TRUE)

# Average alignment stats over marker windows
######################


wrf_rec_p = ggplot(marker_windows, aes(x=recomb.rate, y=first.last.wrf)) +
  geom_point(aes(color=high.wrf), size=3, alpha=0.4) +
  geom_smooth(method="lm", se=F, color="#666666", linetype="dashed", size=1.5) +
  xlab("Recombination rate (5Mb windows)") +
  ylab("wRF of first and last 10kb sub-window") +
  scale_color_manual(values=corecol(numcol=2, pal="wilke", offset=4)) +
  bartheme() +
  theme(legend.position="none")
print(wrf_rec_p)


wrf_rec_hist = ggplot(marker_windows, aes(x=recomb.rate)) +
  geom_histogram(aes(fill=high.wrf), bins=40, alpha=0.6, color="#999999") +
  #xlab("wRF of first and last 10kb sub-window") +
  xlab("Recombination rate (5Mb windows)") +
  ylab("# of 5Mb windows") +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_manual(values=corecol(numcol=2, pal="wilke", offset=4)) +
  bartheme() +
  theme(legend.position="none")
print(wrf_rec_hist)

r_plots = plot_grid(wrf_rec_p, wrf_rec_hist, ncol=2)

# What is this group of windows with high wRF?
####################
# Is it correlated with...

plot_list = list()

wrf_dist_p = ggplot(marker_windows, aes(x=first.last.dist, y=recomb.rate)) +
  geom_point(aes(color=high.wrf), size=3, alpha=0.4) +
  geom_smooth(method="lm", se=F, color="#666666", linetype="dashed", size=1.5) +
  xlab("Distance between first and last 10kb sub-window") +
  #ylab("wRF of first and last 10kb sub-window") +
  ylab("Recombination rate (5Mb windows)") +
  scale_color_manual(values=corecol(numcol=2, pal="wilke", offset=4)) +
  bartheme() +
  theme(legend.position="none")
print(wrf_dist_p)
plot_list[["dist"]] = wrf_dist_p

# Distance?
####################

wrf_decay_p = ggplot(marker_windows, aes(x=wrf.decay, y=recomb.rate)) +
  geom_point(aes(color=high.wrf), size=3, alpha=0.4) +
  geom_smooth(method="lm", se=F, color="#666666", linetype="dashed", size=1.5) +
  xlab("Rate of decay of wRF") +
  #ylab("wRF of first and last 10kb sub-window") +
  ylab("Recombination rate (5Mb windows)") +
  scale_color_manual(values=corecol(numcol=2, pal="wilke", offset=4)) +
  bartheme() +
  theme(legend.position="none")
print(wrf_decay_p)
plot_list[["decay"]] = wrf_decay_p

# Decay?
####################

wrf_inv_p = ggplot(marker_windows, aes(x=inversion.alns, y=recomb.rate)) +
  geom_point(aes(color=high.wrf), size=3, alpha=0.4) +
  geom_smooth(method="lm", se=F, color="#666666", linetype="dashed", size=1.5) +
  xlab("Inverted alignments (to rat)") +
  #ylab("wRF of first and last 10kb sub-window") +
  ylab("Recombination rate (5Mb windows)") +
  scale_color_manual(values=corecol(numcol=2, pal="wilke", offset=4)) +
  bartheme() +
  theme(legend.position="none")
print(wrf_inv_p)
plot_list[["inv"]] = wrf_inv_p

# Inversions?
####################

wrf_trans_p = ggplot(marker_windows, aes(x=translocation.alns, y=recomb.rate)) +
  geom_point(aes(color=high.wrf), size=3, alpha=0.4) +
  geom_smooth(method="lm", se=F, color="#666666", linetype="dashed", size=1.5) +
  xlab("Translocated alignments (to rat)") +
  #ylab("wRF of first and last 10kb sub-window") +
  ylab("Recombination rate (5Mb windows)") +
  scale_color_manual(values=corecol(numcol=2, pal="wilke", offset=4)) +
  bartheme() +
  theme(legend.position="none")
print(wrf_trans_p)
plot_list[["trans"]] = wrf_trans_p

# Traslocations?
####################

wrf_trans_inv_p = ggplot(marker_windows, aes(x=inverted.translocation.alns, y=recomb.rate)) +
  geom_point(aes(color=high.wrf), size=3, alpha=0.4) +
  geom_smooth(method="lm", se=F, color="#666666", linetype="dashed", size=1.5) +
  xlab("Inverted translocation alignments (to rat)") +
  #ylab("wRF of first and last 10kb sub-window") +
  ylab("Recombination rate (5Mb windows)") +
  scale_color_manual(values=corecol(numcol=2, pal="wilke", offset=4)) +
  bartheme() +
  theme(legend.position="none")
print(wrf_trans_inv_p)
plot_list[["trans_inv"]] = wrf_trans_inv_p

# Traslocated inversions?
####################

wrf_genes_p = ggplot(marker_windows, aes(x=genes, y=recomb.rate)) +
  geom_point(aes(color=high.wrf), size=3, alpha=0.4) +
  geom_smooth(method="lm", se=F, color="#666666", linetype="dashed", size=1.5) +
  xlab("Genes") +
  #ylab("wRF of first and last 10kb sub-window") +
  ylab("Recombination rate (5Mb windows)") +
  scale_color_manual(values=corecol(numcol=2, pal="wilke", offset=4)) +
  bartheme() +
  theme(legend.position="none")
print(wrf_genes_p)
plot_list[["genes"]] = wrf_genes_p

# Genes?
####################

wrf_psgenes_p = ggplot(marker_windows, aes(x=ps.genes, y=recomb.rate)) +
  geom_point(aes(color=high.wrf), size=3, alpha=0.4) +
  geom_smooth(method="lm", se=F, color="#666666", linetype="dashed", size=1.5) +
  xlab("Positively selected genes") +
  #ylab("wRF of first and last 10kb sub-window") +
  ylab("Recombination rate (5Mb windows)") +
  scale_color_manual(values=corecol(numcol=2, pal="wilke", offset=4)) +
  bartheme() +
  theme(legend.position="none")
print(wrf_psgenes_p)
plot_list[["psgenes"]] = wrf_psgenes_p

# Positively selected genes?
####################

wrf_uces_p = ggplot(marker_windows, aes(x=uces, y=recomb.rate)) +
  geom_point(aes(color=high.wrf), size=3, alpha=0.4) +
  geom_smooth(method="lm", se=F, color="#666666", linetype="dashed", size=1.5) +
  xlab("UCEs") +
  #ylab("wRF of first and last 10kb sub-window") +
  ylab("Recombination rate (5Mb windows)") +
  scale_color_manual(values=corecol(numcol=2, pal="wilke", offset=4)) +
  bartheme() +
  theme(legend.position="none")
print(wrf_uces_p)
plot_list[["uces"]] = wrf_uces_p

# UCEs?
####################

wrf_hs_p = ggplot(marker_windows, aes(x=hotspots, y=recomb.rate)) +
  geom_point(aes(color=high.wrf), size=3, alpha=0.4) +
  geom_smooth(method="lm", se=F, color="#666666", linetype="dashed", size=1.5) +
  xlab("Recombination hotspots") +
  #ylab("wRF of first and last 10kb sub-window") +
  ylab("Recombination rate (5Mb windows)") +
  scale_color_manual(values=corecol(numcol=2, pal="wilke", offset=4)) +
  bartheme() +
  theme(legend.position="none")
print(wrf_hs_p)
plot_list[["hs"]] = wrf_hs_p

# Positively selected genes?
####################

wrf_rep_p = ggplot(marker_windows, aes(x=avg.repeat, y=recomb.rate)) +
  geom_point(aes(color=high.wrf), size=3, alpha=0.4) +
  geom_smooth(method="lm", se=F, color="#666666", linetype="dashed", size=1.5) +
  xlab("Avg. % repeat per 10kb sub-window") +
  #ylab("wRF of first and last 10kb sub-window") +
  ylab("Recombination rate (5Mb windows)") +
  scale_color_manual(values=corecol(numcol=2, pal="wilke", offset=4)) +
  bartheme() +
  theme(legend.position="none")
print(wrf_rep_p)
plot_list[["rep"]] = wrf_rep_p

# Repeats?
####################

wrf_var_p = ggplot(marker_windows, aes(x=avg.perc.var.sites, y=recomb.rate)) +
  geom_point(aes(color=high.wrf), size=3, alpha=0.4) +
  geom_smooth(method="lm", se=F, color="#666666", linetype="dashed", size=1.5) +
  xlab("Avg. % variable sites per 10kb sub-window") +
  #ylab("wRF of first and last 10kb sub-window") +
  ylab("Recombination rate (5Mb windows)") +
  scale_color_manual(values=corecol(numcol=2, pal="wilke", offset=4)) +
  bartheme() +
  theme(legend.position="none")
print(wrf_var_p)
plot_list[["var"]] = wrf_var_p

# Variable sites?
####################

wrf_inf_p = ggplot(marker_windows, aes(x=avg.perc.inf.sites, y=recomb.rate)) +
  geom_point(aes(color=high.wrf), size=3, alpha=0.4) +
  geom_smooth(method="lm", se=F, color="#666666", linetype="dashed", size=1.5) +
  xlab("Avg. % informative sites per 10kb sub-window") +
  #ylab("wRF of first and last 10kb sub-window") +
  ylab("Recombination rate (5Mb windows)") +
  scale_color_manual(values=corecol(numcol=2, pal="wilke", offset=4)) +
  bartheme() +
  theme(legend.position="none")
print(wrf_inf_p)
plot_list[["inf"]] = wrf_inf_p

# Informative sites?
####################

if(save_fig){
  figfile = here("figs", "supp", "figS11.png")
  cat(as.character(Sys.time()), " | FigS-svs: Saving figure:", figfile, "\n")
  q_plots = plot_grid(plotlist=plot_list, nrow=3)
  p = plot_grid(r_plots, q_plots, nrow=2, rel_heights=c(1,2), align='vh')
  ggsave(figfile, p, height=20, width=20, unit="in")
}






