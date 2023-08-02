############################################################
# For penn windows, 10.20
# Displays histograms for recombination rate for each chromosome.
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(ggplot2)


############################################################
cat("----------\n")

window_size = 10
# Window size in kb

infile = paste("../../data-2/", window_size, "kb-0.5-0.5-topo-counts-tt.csv", sep="")
# Input options
######################
if(read_data){
  win_data = read.csv(infile, header=T)
  # Tree data
  
  #chrome_info = read.csv("../../data/recombination-markers/chrome-stats.csv", comment.char="#", header=T)
  marker_data = subset(win_data, select=c(chr, marker.window, num.markers, marker.slope))
  marker_data = unique(marker_data)
  marker_data$chr = gsub("chr", "", marker_data$chr)
  marker_data$chr = factor(marker_data$chr, levels=c(as.character(sort(as.numeric(levels(as.factor(marker_data$chr))))), "X"))
  # Subset to marker data
  
  # marker_medians = aggregate(marker.slope ~ chr, marker_data, median)
  # marker_means = aggregate(marker.slope ~ chr, marker_data, mean)
  # marker_vars = aggregate(marker.slope ~ chr, marker_data, var)
  # marker_sds = aggregate(marker.slope ~ chr, marker_data, sd)
  # 
  # marker_stats = cbind(marker_medians, marker_means, marker_vars, marker_sds)
  # names(marker_stats) = c("chr", "median", "chr2", "mean", "chr3", "var", "chr4", "sd")
  # marker_stats = subset(marker_stats, select=c(chr, median, mean, var, sd))
}
# Read  and filter the data
######################

marker_hists = ggplot(marker_data, aes(x=marker.slope, fill=chr)) +
  geom_histogram(color="#333333") +
  facet_wrap(vars(chr)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab("Recombination rate") +
  ylab("# 5Mb windows") +
  bartheme() +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1))

print(marker_hists)
outfile = "../../figs-2/recombination-rates.png"
ggsave(outfile, marker_hists, height=8, width=12, units="in")
  

# Recombination dists.