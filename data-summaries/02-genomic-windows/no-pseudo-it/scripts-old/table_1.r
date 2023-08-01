############################################################
# For penn genomes, 12.20
# Counts for tables, etc.
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(dplyr)
#library(plyr)
library(stringr)

############################################################
cat("----------\n")

window_size_str = 10
# Window size in kb

window_size = window_size_str * 1000
# Window size in bp

au_flag = F
# Set to filter out windows that don't pass the AU test

read_data = T
# Set to read and filter data or not (if its already stored in mem)

infile = paste("../../windows/data-2/", window_size_str, "kb-0.5-0.5-topo-counts-tt.csv", sep="")
outfile = "../tables/table1.csv"
# Input options
######################

if(read_data){
  cat(as.character(Sys.time()), "| Reading input data:", infile, "\n")
  all_windows = read.csv(infile, header=T)
  all_windows_f = subset(all_windows, repeat.filter=="PASS" & missing.filter=="PASS")
  if(au_flag){
    all_windows_f = subset(all_windows_f, AU.test=="PASS")
  }
}
# Read and filter the data
######################

cat("Total windows:            ", nrow(all_windows), "\n")
cat("Total windows with trees: ", nrow(all_windows_f), "\n")

chrome_topos = subset(all_windows_f, select=c(chr, topo.num.chrome))

chrome_topos = chrome_topos %>%
  group_by(chr) %>%
  summarize(num.topos=max(topo.num.chrome))


top_topos = subset(all_windows_f, topo.rank.overall <= 8)
top_topos = unique(subset(top_topos, select=c(topo.rank.overall, topo, topo.count.overall)))

top_topos$prop = signif(top_topos$topo.count.overall / nrow(all_windows_f), 3)
top_topos$topo = str_replace_all(top_topos$topo, regex("<[0-9]+>"), "")
top_topos$topo = str_replace_all(top_topos$topo, regex("RhySor"), "RS")
top_topos$topo = str_replace_all(top_topos$topo, regex("GraDol"), "GD")
top_topos$topo = str_replace_all(top_topos$topo, regex("RhaDil"), "RD")
top_topos$topo = str_replace_all(top_topos$topo, regex("mm10"), "MM")
top_topos$topo = str_replace_all(top_topos$topo, regex("HylAll"), "HA")
top_topos$topo = str_replace_all(top_topos$topo, regex("MasNat"), "MN")
top_topos$topo = str_replace_all(top_topos$topo, regex("PraDel"), "PD")

top_topos = arrange(top_topos, topo.rank.overall)
write.csv(top_topos, outfile, row.names=F)
## TABLE 1
######################