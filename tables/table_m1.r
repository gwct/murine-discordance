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

window_size = 10
# Window size in kb

au_flag = F
# Set to filter out windows that don't pass the AU test

read_data = T
# Set to read and filter data or not (if its already stored in mem)

datadir = here("summary-data", "02-genomic-windows")

infile = here(datadir, paste(window_size, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts.csv", sep=""))

outfile = here("tables", "tableM1.csv")

# Input options
######################

if(read_data){
  cat(as.character(Sys.time()), "| Reading input data:", infile, "\n")
  all_windows = read.csv(infile, comment.char="#", header=T)
  all_windows_f = subset(all_windows, repeat.filter=="PASS" & missing.filter=="PASS")
  if(au_flag){
    all_windows_f = subset(all_windows_f, AU.test=="PASS")
  }
}
# Read and filter the data
######################

cat("Total windows:            ", nrow(all_windows), "\n")
cat("Total windows with trees: ", nrow(all_windows_f), "\n")

all_windows$chr = str_replace_all(all_windows$chr, regex("chr"), "")

windows_rep_f = subset(all_windows, repeat.filter=="FILTER")
windows_mis_f = subset(all_windows, repeat.filter=="PASS" & missing.filter=="FILTER")

all_windows_chr = all_windows %>% group_by(chr, chr.len) %>% summarize(Total.windows=n())
names(all_windows_chr)[2] = "Chromosome.length"
rep_windows_chr = windows_rep_f %>% group_by(chr) %>% summarize(Repeat.windows=n())
mis_windows_chr = windows_mis_f %>% group_by(chr) %>% summarize(Missing.windows=n())

all_rep_windows_chr = merge(all_windows_chr, rep_windows_chr, by="chr")
all_rep_mis_windows_chr = merge(all_rep_windows_chr, mis_windows_chr, by="chr")
all_rep_mis_windows_chr$Final.windows = all_rep_mis_windows_chr$Total.windows - all_rep_mis_windows_chr$Repeat.windows - all_rep_mis_windows_chr$Missing.windows
#all_rep_mis_windows_chr$chr = as.character(all_rep_mis_windows_chr$chr)
all_rep_mis_windows_chr = all_rep_mis_windows_chr[order(all_rep_mis_windows_chr$chr),]



totals = data.frame("chr"="Total", "Chromosome.length"=sum(all_rep_mis_windows_chr$Chromosome.length),
                    "Total.windows"=sum(all_rep_mis_windows_chr$Total.windows),
                    "Repeat.windows"=sum(all_rep_mis_windows_chr$Repeat.windows),
                    "Missing.windows"=sum(all_rep_mis_windows_chr$Missing.windows),
                    "Final.windows"=sum(all_rep_mis_windows_chr$Final.windows))



all_rep_mis_windows_chr = rbind(all_rep_mis_windows_chr, totals)


write.csv(all_rep_mis_windows_chr, outfile, row.names=F)
## TABLE M1
######################