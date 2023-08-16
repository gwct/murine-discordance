# For rodent genomes
# Window/feature overlaps for Discussion
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(tidyverse)
library(ape)
library(phangorn)

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

do_rf = F
# Whether or not to calculate RFs to the species tre

if(do_rf){
  infile = here("summary-data", "02-genomic-windows", paste0(window_size, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts.csv"))
}else{
  infile = here("summary-data", "02-genomic-windows", paste0(window_size, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts-st-rf.csv"))
}
# The window file, depending on if the RFs need to be done

gene_file = here("summary-data", "02-genomic-windows", "feature-beds", "mm10-genes.bed")
uce_file = here("summary-data", "02-genomic-windows", "feature-beds", "mm10-uce.bed")
# Feature beds

concat_tree_file = here("summary-data", "03-selection-tests", "concat.cf.rooted.tree")
# Species tree

# Input options
######################

if(read_data){
  cat(as.character(Sys.time()), " | Reading window data: ", infile, "\n")
  all_windows = read_csv(infile, comment="#")
  names(all_windows) = make.names(names(all_windows))
  all_windows$midpoint = (all_windows$start + all_windows$end) / 2
  all_windows_f = all_windows %>%
    filter(repeat.filter=="PASS", missing.filter=="PASS")
  # Reading 10kb window tree file
  
  cat(as.character(Sys.time()), " | Reading bed: ", gene_file, "\n")
  genes = read.table(gene_file)
  names(genes) = c("chr", "start", "end", "id", "score", "no")
  genes$start = as.numeric(genes$start)
  genes$end = as.numeric(genes$end)
  genes$len = genes$end - genes$start
  genes$midpoint = (genes$start + genes$end) / 2
  genes$start = genes$midpoint
  genes$end = genes$midpoint
  # Reading gene coords
  
  cat(as.character(Sys.time()), " | Reading bed: ", uce_file, "\n")
  uces = read.table(uce_file)
  names(uces) = c("chr", "start", "end", "id", "score", "no", "tmp1", "tmp2", "tmp3")
  uces$start = as.numeric(uces$start)
  uces$end = as.numeric(uces$end)
  uces$len = uces$end - uces$start
  uces$midpoint = (uces$start + uces$end) / 2
  uces$start = uces$midpoint
  uces$end = uces$midpoint
  # Reading UCE coords
  
  concat_tree = read.tree(file=concat_tree_file)
  # Read species tree
}

######################

if(do_rf){
  all_windows_f$st.rf = NA
  for(i in 1:nrow(all_windows_f)){
    tree = read.tree(text=all_windows_f[i,]$unparsed.tree)
    all_windows_f[i,]$st.rf = RF.dist(tree, concat_tree)
  }
  # Calculate RF from each window tree to the species tree
  
  window_rf_file = here("summary-data", "02-genomic-windows", paste0(window_size, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts-st-rf.csv"))
  write.csv(all_windows_f, file=window_rf_file, row.names=F)
  # Write out the RF data
}

######################

cat("----------\n")
window_range = makeGRangesFromDataFrame(all_windows_f)
# Create range for marker windows

gene_range = makeGRangesFromDataFrame(genes)
gene_overlaps = countOverlaps(window_range, gene_range)
all_windows_f$gene.overlap = gene_overlaps
# Create range for genes and get overlaps with marker windows

gene_windows = all_windows_f %>% filter(gene.overlap > 0)
gene_concord_windows = gene_windows %>% filter(concat.chrome.topo == TRUE)
perc_gene_concord = nrow(gene_concord_windows) / nrow(gene_windows)
cat("% window trees with genes concordant with species tree: ", perc_gene_concord, "\n")

######################

uce_range = makeGRangesFromDataFrame(uces)
uce_overlaps = countOverlaps(window_range, uce_range)
all_windows_f$uce.overlap = uce_overlaps
# Create range for uces and get overlaps with marker windows

uce_windows = all_windows_f %>% filter(uce.overlap > 0)
uce_concord_windows = uce_windows %>% filter(st.rf == 0)
perc_uce_concord = nrow(uce_concord_windows) / nrow(uce_windows)
cat("% window trees with UCEs concordant with species tree:  ", perc_uce_concord, "\n")

######################

all_concord_windows = all_windows_f %>% filter(st.rf == 0)
perc_all_concord = nrow(all_concord_windows) / nrow(all_windows_f)
cat("% window trees concordant with species tree:            ", perc_all_concord, "\n")
cat("----------\n")

######################





