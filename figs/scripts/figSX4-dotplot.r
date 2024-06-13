############################################################
# For rodent genomes, 07.20
# Dot plots
# With code from: https://github.com/tpoorten/dotPlotly
# Gregg Thomas
############################################################

library(tidyverse)
library(cowplot)
library(grid)
library(viridis)
library(circlize)
library(scales)
library(here)

############################################################
# Reading input data

read_data = T
# Whether or not to (re)-read the data

source_pafr = T
# Whether or not to source the pafr package

save_fig = F
# Whether or not to save the figures

if(source_pafr){
  #setwd("C:/Users/grt814/bin/pafr/R/")
  cat(as.character(Sys.time()), " | Sourcing pafr\n")
  setwd("C:/bin/pafr/R/")
  #source("C:/Users/grt814/bin/pafr/")
  files.sources = list.files()
  sapply(files.sources, source)
}
# Source the pafr package if specified

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
# Switch to the current directory

paf_file = here("summary-data", "04-mouse-rat-svs", "rn6-mm10.paf")
#paf_file = paste(datadir, "/rn6-mm10.paf", sep="")
# The alignment file

# Input options
######################

if(read_data){
  cat(as.character(Sys.time()), " | Reading alignment file: ", paf_file, "\n")
  paf = read_paf(paf_file)
  paf = subset(paf, nchar(qname) <= 5 & nchar(tname) <= 5)
  paf = subset(paf, qname != "chrM" & tname != "chrM")
  paf = subset(paf, qname != "chrY" & tname != "chrY")
  paf_long = subset(paf, alen > 10000)
}
# Read the alignment

######################

mouse_chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", 
              "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX")
rat_chr =  c("chr9", "chr13", "chr3", "chr2", "chr5", "chr14", "chr12", "chr4", "chr1", "chr19", 
             "chr16", "chr8", "chr20", "chr10", "chr6", "chr17", "chr15", "chr7", "chr11", "chr18", "chrX")
chr_orders = list(mouse_chr, rat_chr)
# Set the order of the chromosomes


cat(as.character(Sys.time()), " | Rendering dotplot\n")
paf_result = dotplot(paf, order_by="provided", label_seqs=TRUE, xlab="Mouse", ylab="Rat", ordering=chr_orders)
# Process paf with pafr

figS8 = paf_result[[1]] + 
  dottheme() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),)
print(figS8)
# Generate the dotplot

if(save_fig){
  figSfile = here("figs", "supp", "figS6.pdf")
  cat(as.character(Sys.time()), " | Saving supplemental figure:", figSfile, "\n")
  ggsave(filename=figSfile, figS8, width=8, height=8, unit="in")
}

######################

cat(as.character(Sys.time()), " | Counting alignment chunk sizes\n")

paf_bins = data.frame("bin"=c(10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000), "count"=c(0,0,0,0,0,0,0,0,0), "size"=c(0,0,0,0,0,0,0,0,0),
                      "bin.labels"=c("1bp-10bp", "10bp-100bp", "100bp-1kb", "1kb-10kb", "10kb-100kb", "100kb-1Mb", "1Mb-10Mb", "10Mb-100Mb", "100Mb-1Gb"))
# Setting bins for aligned chunks as Ensembl does
# https://useast.ensembl.org/info/genome/compara/mlss.html?mlss=1981

for(i in 1:nrow(paf)){
# Go through every alignment block
  
  cur_len = paf[i,]$alen
  # Get the size of the block
  
  cur_bin = 10
  for(bin in paf_bins$bin){
    if(bin > cur_len){
      break
    }else{
      cur_bin = bin
    }
  }
  # Loop through the bins to find which one the current block belongs in
  
  paf_bins[paf_bins$bin==cur_bin,]$count = paf_bins[paf_bins$bin==cur_bin,]$count + 1
  paf_bins[paf_bins$bin==cur_bin,]$size = paf_bins[paf_bins$bin==cur_bin,]$size + cur_len
  # Add one to the count for the bin and add the length of the block to the size of the bin
}

cat(as.character(Sys.time()), " | Rendering alignment chunk distributions\n")

paf_bins$bin.labels = factor(paf_bins$bin.labels, levels=rev(c("1bp-10bp", "10bp-100bp", "100bp-1kb", "1kb-10kb", "10kb-100kb", "100kb-1Mb", "1Mb-10Mb", "10Mb-100Mb", "100Mb-1Gb")))
# Re-order the bin labels

figS9_a = ggplot(paf_bins, aes(x=bin.labels, y=count)) +
  geom_bar(stat="identity", fill=corecol(pal="wilke", numcol=1, "offset"=0)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab("Size range") +
  ylab("# blocks\n") + 
  bartheme() +
  theme(axis.text.x=element_text(angle=20, hjust=1, size=8),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        plot.margin=margin(0.5,0.5,0.5,0.5, unit="cm"),
        legend.position="none") +
  coord_flip()
#print(p1)
# Plot the distribution of bin counts

paf_bins$size.mb = paf_bins$size / 1000000

figS9_b = ggplot(paf_bins, aes(x=bin.labels, y=size.mb)) +
  geom_bar(stat="identity", fill=corecol(pal="wilke", numcol=1, "offset"=1)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab("") +
  ylab("Total bp in Mb\n(including gaps)") + 
  bartheme() +
  theme(axis.text.x=element_text(angle=20, hjust=1, size=8),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        plot.margin=margin(0.5,0.5,0.5,0.5, unit="cm"),
        legend.position="none") +
  coord_flip()
#print(p2)
# Plot the distribution of bin sizes

figS9_ab = plot_grid(figS9_a, figS9_b, ncol=2, labels=c("A", "B"), label_size=12)
# Combine the two block panels

######################

cat(as.character(Sys.time()), " | Counting distances between aligned chunks\n")

paf_sort = paf[order(paf$qname, paf$qstart), ]
# Make sure the data frame is sorted by query start

q_bed = paf_sort %>% select(qname, qstart, qend)
q_inter_dist = rainfallTransform(q_bed, mode="left")
# Query (mouse) inter distances

paf_sort = paf[order(paf$qname, paf$tstart), ]
# Make sure the data frame is sorted by target start

t_bed = paf_sort %>% select(tname, tstart, tend)
t_inter_dist = rainfallTransform(t_bed, mode="left")
# Target (rat) inter distances

# Transform the query and target segments into bed format and get
# the interval distance between aligned chunks
######################

cat(as.character(Sys.time()), " | Rendering inter-dist distributions\n")

names(q_inter_dist) = c("name", "start", "end", "dist")
q_inter_dist$label = "Mouse"
names(t_inter_dist) = c("name", "start", "end", "dist")
t_inter_dist$label = "Rat"
# Add names and labels to both inter dist dfs

inter_dists = rbind(q_inter_dist, t_inter_dist)
# Combine the inter dist dfs

inter_dists = subset(inter_dists, dist > 0)
# Only get chunks that don't overlap

avg_mouse_dist = mean(subset(q_inter_dist, dist<20000 & dist>0)$dist, na.rm=T)
avg_rat_dist = mean(subset(t_inter_dist, dist<20000 & dist>0)$dist, na.rm=T)
# Calculate averages

figS9_c = ggplot(subset(inter_dists, dist<20000), aes(x=label, y=dist, group=label, fill=label)) +
  geom_violin() +
  geom_boxplot(width=0.05, outlier.shape=NA, fill="#ffffff") +
  geom_hline(yintercept=avg_rat_dist, linetype=1, size=1, color="#ffffff") +
  geom_hline(yintercept=avg_rat_dist, linetype=2, size=1, color=corecol(numcol=1, pal="trek", offset=2)) +
  annotate("text", x=2.4, y=avg_rat_dist+650, label=signif(avg_rat_dist,5), size=3) +
  geom_hline(yintercept=avg_mouse_dist, linetype=1, size=1, color="#ffffff") +
  geom_hline(yintercept=avg_mouse_dist, linetype=2, size=1, color=corecol(numcol=1, pal="trek", offset=1)) +
  annotate("text", x=2.4, y=avg_mouse_dist+650, label=signif(avg_mouse_dist,5), size=3) +
  xlab("") +
  ylab("Distance between\naligned segments (bp)") +
  scale_fill_manual(values=corecol(numcol=2, pal="trek", offset=1)) +
  bartheme() +
  theme(axis.text.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        plot.margin=margin(1,0.5,0,0.5, unit="cm"),
        legend.position="none")
print(figS9_c)
# Render the inter dist figure

######################

cat(as.character(Sys.time()), " | Combining panels\n")

figS9_bot = plot_grid(NULL, figS9_c, NULL, ncol=3, rel_widths=c(0.25, 1, 0.25), labels=c("", "C", ""))

figS9 = plot_grid(figS9_ab, figS9_bot, nrow=2)

if(save_fig){
  figSfile = here("figs", "supp", "figS7.png")
  cat(as.character(Sys.time()), " | Saving supplemental figure:", figSfile, "\n")
  ggsave(filename=figSfile, figS9, width=6, height=8, unit="in")
}

######################

