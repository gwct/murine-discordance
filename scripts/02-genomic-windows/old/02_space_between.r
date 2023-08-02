############################################################
# For penn genomes, 04.20
# Plots distance between windows with trees.
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(ggplot2)
library(cowplot)
library(ggbeeswarm)

############################################################
cat("----------\n")

window_size = 10
# Window size in kb

au_flag = FALSE
# Set to filter out windows that don't pass the AU test

skip_one = FALSE
# Set to only do one test chromosome

infile = paste("../../data-2/", window_size, "kb-0.5-0.5-topo-counts-tt.csv", sep="")
# Input options
######################


all_windows = read.csv(infile, header=T)
all_windows_f = subset(all_windows, repeat.filter=="PASS" & missing.filter=="PASS")
if(au_flag){
  all_windows_f = subset(all_windows_f, AU.test=="PASS")
}
# Read and filter the data
######################

space_between = data.frame("chrome"=c(), "dist"=c())

for(chrome in levels(all_windows$chr)){
  # Process every chromosome
  if(skip_one && chrome!="chr19"){
    next
  }
  if(nchar(chrome) == 4 && chrome != "chrX"){
    out_chrome = gsub("chr", "chr0", chrome)
  }else{
    out_chrome = chrome
  }
  # The string for the current chromosome.
  
  cat(as.character(Sys.time()), " | Starting chromosome ", out_chrome, "\n")
  
  chrdata = subset(all_windows, chr==chrome)
  total_windows = length(chrdata[,1])
  chrdata_f = subset(chrdata, repeat.filter=="PASS" & missing.filter=="PASS")
  if(au_flag){
    chrdata_f = subset(chrdata_f, AU.test=="PASS")
  }
  used_windows = length(chrdata_f[,1])
  # Subset and filter data for current chromosome.
  
  first = TRUE
  for(i in 1:nrow(chrdata_f)) {
    window <- chrdata_f[i,]
    if(first){
      prev_start = window$start
      first = FALSE
      next
    }
    
    cur_dist = window$start - prev_start
    space_between = rbind(space_between, data.frame("chrome"=chrome, "dist"=cur_dist))
    
    prev_start = window$start
  }
}



cat(as.character(Sys.time()), " | ", "Plotting figure...\n")
space_between = subset(space_between, dist < 50000)
space_between$chrome = gsub("chr", "", space_between$chrome)
space_between$chrome = factor(space_between$chrome, levels=c(as.character(sort(as.numeric(levels(as.factor(space_between$chrome))))), "X"))


space_p = ggplot(space_between, aes(x=dist, fill=chrome)) +
  #geom_quasirandom(size=2, alpha=0.7, width=0.25, color="#d3d3d3") +
  #geom_boxplot(outlier.shape=NA, alpha=0.7, width=0.5) +
  geom_histogram(stat="count", color="#333333") +
  facet_wrap(vars(chrome)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab("Distance (bp) between adjacent windows") + 
  ylab("# windows") +
  bartheme() +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 45, hjust = 1))
print(space_p)

outfile = paste("../../figs-2/", window_size, "kb-space-between", sep="")
if(au_flag){
  outfile = paste(outfile, "-au-filter", sep="")
}
outfile = paste(outfile, ".png", sep="")
# Output information for this chrome

cat(as.character(Sys.time()), " | ", out_chrome, "Saving figure: ", outfile, "\n")
ggsave(filename=outfile, space_p, width=12, height=8, units="in")
cat("----------\n")
