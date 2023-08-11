############################################################
# For penn genomes, 12.20
# Counts for tables, etc.
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(tidyverse)

############################################################
cat("----------\n")

window_size = 10
# Window size in kb

au_flag = F
# Set to filter out windows that don't pass the AU test

read_data = T
# Set to read and filter data or not (if its already stored in mem)

datadir = here("summary-data", "02-genomic-windows")

infile = here(datadir, paste(window_size, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts.csv", sep=""))

outfile1 = here("tables", "table1.csv")
outfile2 = here("tables", "table2.csv")

# Input options
######################

if(read_data){
  cat(as.character(Sys.time()), "| Reading input data:", infile, "\n")
  all_windows = read_csv(infile, comment="#")
  names(all_windows) = make.names(names(all_windows))
  all_windows_f = subset(all_windows, repeat.filter=="PASS" & missing.filter=="PASS")
  if(au_flag){
    all_windows_f = subset(all_windows_f, AU.test=="PASS")
  }
}
# Read and filter the data
######################

cat("Total windows:            ", nrow(all_windows), "\n")
cat("Total windows with trees: ", nrow(all_windows_f), "\n")

top_topos = subset(all_windows_f, topo.rank.overall <= 8)
top_topos = unique(subset(top_topos, select=c(topo.rank.overall, topo, topo.count.overall)))

top_topos$prop = signif(top_topos$topo.count.overall / nrow(all_windows_f), 3)
top_topos$topo = str_replace_all(top_topos$topo, regex("<[0-9]+>"), "")
top_topos$topo = str_replace_all(top_topos$topo, regex("rsor"), "RS")
top_topos$topo = str_replace_all(top_topos$topo, regex("gdol"), "GD")
top_topos$topo = str_replace_all(top_topos$topo, regex("rdil"), "RD")
top_topos$topo = str_replace_all(top_topos$topo, regex("mm10"), "MM")
top_topos$topo = str_replace_all(top_topos$topo, regex("hall"), "HA")
top_topos$topo = str_replace_all(top_topos$topo, regex("mnat"), "MN")
top_topos$topo = str_replace_all(top_topos$topo, regex("pdel"), "PD")

top_topos = arrange(top_topos, topo.rank.overall)
write.csv(top_topos, outfile1, row.names=F)
## TABLE 1
######################

chrome_topos = all_windows_f %>%
  select(chr, topo.num.chrome) %>%
  group_by(chr) %>%
  summarize(num.topos=max(topo.num.chrome)) %>%
  mutate(chr=gsub("chr", "", chr))

top_chrome_topos = all_windows_f %>%
  group_by(chr) %>%
  filter(topo.rank.chrome <= 3) %>%
  select(topo.rank.chrome, topo.rank.overall) %>%
  unique() %>%
  arrange(topo.rank.overall) %>%
  mutate(rank.str = paste0(topo.rank.overall, collapse=",")) %>%
  select(chr, rank.str) %>%
  unique() %>%
  mutate(chr=gsub("chr", "", chr))

x_chrome_topos = top_chrome_topos %>% filter(chr=="X")
a_chrome_topos = top_chrome_topos %>% 
  filter(chr!="X") %>%
  mutate(chr = as.numeric(chr)) %>%
  arrange(chr) %>%
  mutate(chr = as.character(chr))

top_chrome_topos = rbind(a_chrome_topos, x_chrome_topos)

chrome_topos = left_join(top_chrome_topos, chrome_topos, by="chr")
write.csv(chrome_topos, outfile2, row.names=F)
## TABLE 2
######################