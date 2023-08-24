############################################################
# For rodent genomes
# Stat calculation for last section of Results
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(tidyverse)
library(ape)

############################################################

window_size = 10
# Window size in kb

chrome_info_file = here("summary-data", "02-genomic-windows", "recombination-markers", "chrome-stats.csv")

transcript_windows_file = here("summary-data", "03-selection-tests", "mm10-cds-windows.csv")
longest_transcripts_file = here("summary-data", "03-selection-tests", "mm10.ensGene.chromes.longest.bed")

cat(as.character(Sys.time()), " | Reading gene trees\n")
transcript_windows = read_csv(transcript_windows_file, comment="#")
names(transcript_windows) = make.names(names(transcript_windows))
# The 10kb windows that overlap with mouse transcripts

longest_transcripts = read_tsv(longest_transcripts_file, comment="#", 
                               col_names=c("chr", "transcript.start", "transcript.end", "transcript.id", "score", 
                                           "strand", "cds.start", "cds.end", "color", "num.exons", "exon.lens", "exon.starts"))
# The transcripts used for gene trees and dN/dS

longest_transcripts = longest_transcripts %>% mutate(cds.len = cds.end - cds.start)

transcript_windows_conc = transcript_windows %>% filter(gene.tree.match.window.tree == TRUE)
prop_conc_cds_windows = nrow(transcript_windows_conc) / nrow(transcript_windows)

cat(" ---------- \n")
cat("# transcripts:                                                               ", nrow(longest_transcripts), "\n")
cat("# of times a coding sequences overlaps with a 10kb window:                   ", nrow(transcript_windows), "\n")
cat("Avg. distance from start to end of coding sequence:                          ", mean(longest_transcripts$cds.len), "\n")
cat("Proportion of coding windows where gene tree is concordant with window tree: ", prop_conc_cds_windows, "\n")
cat(" ---------- \n")

############################################################

distdir = "D:/data/rodent-genomes/dists/"

cat(as.character(Sys.time()), " | Reading chrome info: ", chrome_info_file, "\n")
chrome_info = read.csv(chrome_info_file, header=T, comment.char="#")

avg_cds_adj = round( mean(longest_transcripts$cds.len) / (window_size * 1000) )

avg_cds_dists = c()

for(chrome in chrome_info$chr){
  
  if(chrome == "chrX"){
    next
  }
  
  print(chrome)
  chr_cds_dists_file = paste0(distdir, "all-dists/", chrome, "-", avg_cds_adj, ".csv")
  chr_cds_dists = read_csv(chr_cds_dists_file)
  
  avg_cds_dists = c(avg_cds_dists, chr_cds_dists$wrf)
  
}

cat("wRF decay over ", avg_cds_adj, " windows:                                   ", mean(avg_cds_dists))

############################################################

window_size = 10
marker_window_size = 5

window_file = here("summary-data", "02-genomic-windows", paste0(window_size, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts.csv"))
cat(as.character(Sys.time()), " | Reading window data: ", window_file, "\n")
all_windows = read_csv(infile, comment="#")
names(all_windows) = make.names(names(all_windows))
all_windows_f = subset(all_windows, repeat.filter=="PASS" & missing.filter=="PASS")

window_trees = read.tree(text=all_windows_f$unparsed.tree)
#total_bl = lapply(window_trees, FUN = function(x) sum(x$edge.length))
#bls = sapply(window_trees, get, x="edge.length")

bls = c()

## Takes a while
#for(i in 1:length(window_trees)){
#  bls = c(bls, window_trees[[i]]$edge.length)
#}
## Takes a while


