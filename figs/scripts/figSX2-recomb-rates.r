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

# datadir = "C:/Users/Gregg/Box Sync/rodents/penn/paper/data/"
# infile = paste(datadir, window_size, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts-tt.csv", sep="")

infile = here("summary-data", "02-genomic-windows", paste0(window_size_kb, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts.csv"))

# Input options
######################

if(read_data){
  cat(as.character(Sys.time()), " | Reading window data: ", infile, "\n")
  win_data = read.csv(infile, comment.char="#", header=T)
}

#chrome_info = read.csv("../../data/recombination-markers/chrome-stats.csv", comment.char="#", header=T)
marker_data = subset(win_data, select=c(chr, marker.window, num.markers, marker.slope))
marker_data = unique(marker_data)
marker_data$chr = gsub("chr", "", marker_data$chr)
marker_data$chr = factor(marker_data$chr, levels=c(as.character(sort(as.numeric(levels(as.factor(marker_data$chr))))), "X"))
# Subset to marker data

# Read  and filter the data
######################

fig = ggplot(marker_data, aes(x=marker.slope, fill=chr)) +
  geom_histogram(color="#333333") +
  facet_wrap(vars(chr)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab("Recombination rate") +
  ylab(paste("# ", marker_window_size, "Mb windows", sep="")) +
  bartheme() +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1, size=10))
print(fig)

if(save_fig){
  figfile = here("figs", "supp", paste0("figSX-", marker_window_size, "mb-rec-rates.png"))
  #figfile = paste("../figs/figSX-", marker_window_size, "mb-rec-rates.png", sep="")
  cat(as.character(Sys.time()), " | Saving figure: ", figfile, "\n")
  ggsave(filename=figfile, fig, width=6, height=6, units="in")
}

