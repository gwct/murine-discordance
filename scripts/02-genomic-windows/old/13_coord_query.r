############################################################
# For penn genomes, 02.2021
# Given a set of coordinates, calculates RF from a given 
# distance around those coordinates and compares to random 
# coordinates (or a given set of null coordinates).
# Gregg Thomas
############################################################

library(ape)
library(phangorn)
library(dplyr)
library(tidyr)
library(doParallel)

#############################################################
# Functions

isFiltered <- function(w){
  filter_flag = FALSE
  filter_str = ""
  
  if(w$repeat.filter!="PASS"){
    filter_str = paste(filter_str, "REPEAT;")
    filter_flag = TRUE
  }
  if(w$missing.filter!="PASS"){
    filter_str = paste(filter_str, "MISSING;")
    filter_flag = TRUE
  }
  if(au_flag && w$AU.test != "PASS"){
    filter_str = paste(filter_str, "AU;")
    filter_flag = TRUE
  }  
  
  return(filter_flag)
}

######################

checkBounds <- function(direction, adj_start, c_len){
  bounds_flag = TRUE
  if(direction=="backward" && adj_start < 0){
    bounds_flag = FALSE
  }
  
  if(direction=="forward" && adj_start > c_len){
    bounds_flag = FALSE
  }
  
  return(bounds_flag)
}

######################

getRandWindow <- function(window_df){
  rand_window = data.frame("window"="", "chr"="", "start"="", "tree"="")
  
  random_window = as.character(sample(window_df$window, 1))
  rand_window$window = random_window
  
  window_data = window_df[window_df$window==random_window,]
  rand_window$chr = window_data$chr
  rand_window$start = window_data$start
  rand_window$tree = as.character(window_data$unparsed.tree)
  
  return(rand_window)
}

######################

getDists <- function(direction, adj, target_window, target_tree, gid, int_start){
  result = data.frame("chr"=target_window$chr, "start"=target_window$start, "end"=target_window$end, 
                      "gene.id"=gid, "int.start"=int_start, 
                      "window"=NA, "adj"=adj, "filt"=NA,
                      "rf"=NA, "wrf"=NA, "spr"=NA, "kf"=NA, "path"=NA)
  
  
  adj_start = target_window$start + (adj * window_size)
  
  #print(paste("        ADJACENCY:", adj))

  #print("ENTER checkBounds")
  within_bounds = checkBounds(direction, adj_start, target_window$chr.len)
  #print(win_filtered)
  # Will be TRUE if window is beyond bounds of chromosome length.
  
  if(within_bounds){
    query_window = subset(windows, chr==target_window$chr & start==adj_start)
    result$window = query_window$window
    #print(query_window)
    # Get the next window based on the current adjaceny step and window size
    
    #print("ENTER checkFilter")
    win_filtered = isFiltered(query_window)
    #print(win_filtered)
    # Will be TRUE if any of the filters of the query window are unpassed.
    
    #print(win_filtered)

    if(!win_filtered){
      query_tree = read.tree(text=as.character(query_window$unparsed.tree))
      # Get the tree for the adjacent window
      
      result$rf = RF.dist(target_tree, query_tree)
      result$wrf = wRF.dist(target_tree, query_tree)
      spr = SPR.dist(target_tree, query_tree)
      names(spr) = NULL
      result$spr = spr
      result$kf = KF.dist(target_tree, query_tree)
      result$path = path.dist(target_tree, query_tree)
      # Calculate the wRF for the two trees    
    }
  }
  
  #print(paste("        RESULT", do.call(paste, c(result, sep="    "))))
  #print(paste("        RESULT DIMENSIONS", dim(result)))
  # print(result)
  # print(rf)
  # print(wrf)
  # print(spr)
  # print(kf)
  # print(path)
  # print("----") 
  return(result)
}

######################

parseCoords <- function(f_gene){

  #gene_dists = data.frame("chr"="NA", "start"="NA", "end"="NA", "gene.id"="NA", "int.start"="NA", 
  #                        "window"="NA", "adj"="NA", "filt"="NA",
  #                        "rf"="NA", "wrf"="NA", "spr"="NA", "kf"="NA", "path"="NA")
  gene_dists = NULL

  #cur_window = subset(windows, chr==f_gene$chr & start <= f_gene$start & end >= f_gene$start)
  cur_window = subset(windows, chr==f_gene$chr & start <= f_gene$int.start & end >= f_gene$int.start)

  print(paste("    START WINDOW", f_gene$gid, do.call(paste, c(cur_window, sep="    "))))
  #print(paste("    START FILTER:", isFiltered(cur_window)))
  
  if(!isFiltered(cur_window)){
    cur_tree = read.tree(text=as.character(cur_window$unparsed.tree))
    # Get the tree for the current window

    cur_adj = 1
    while(cur_adj <= windows_per_interval){
      #print(paste(f_gene$gid, "-", cur_adj, sep=""))
      
      for_df = getDists("forward", cur_adj, cur_window, cur_tree, f_gene$gid, f_gene$int.start)
      gene_dists = rbind(gene_dists, for_df)
      
      back_df = getDists("backward", (-1*cur_adj), cur_window, cur_tree, f_gene$gid, f_gene$int.start)
      gene_dists = rbind(gene_dists, back_df)
      cur_adj = cur_adj + 1
    }
  }
  #print(paste("    DIST DIMENSIONS:", dim(gene_dists)))
  return(gene_dists)
}


############################################################
# Options

server = T
# Set if running on server

window_size_str = 10
# Window size in kb

window_size = window_size_str * 1000
# Window size in bp

marker_window_size = 5
# Marker window size in Mb

num_cores = 46
# Number of cores

au_flag = FALSE
# Set to filter out windows that don't pass the AU test

gen_random = F
# Set to generate the random windows

do_random = F
# Set to do the random windows

read_data = T
# Set depending on whether or not data has been read into environment

feature = "coding-genes"
# One of genes, coding-genes, ps, uces, hotspots

flank_interval_mb = 5
flank_interval = flank_interval_mb * 1000000
# The distance to calculate outward from both sides of a given coordinate in Mb

windows_per_interval = flank_interval / window_size

if(!server){
  this.dir <- dirname(parent.frame(2)$ofile)
  setwd(this.dir)
  num_cores = 1
}

############################################################
# Input

cat(as.character(Sys.time()), " | Feature:", feature, "\n")
if(feature=="genes"){
  infile = "../../../cds/data/mm10-genes.bed"
  logfile = "coord-query-genes-"
  cols = c("chr", "start", "end", "ids", "score", "strand")
}else if(feature=="coding-genes"){
  infile = "../../../cds/data/mm10-coding-genes.bed"
  logfile = "coord-query-coding-genes-"
  #logfile = "test-"
  cols = c("chr", "start", "end", "ids", "score", "strand")
}else if(feature=="ps"){
  infile = "../../../cds/data/mm10-ps-genes-all-gt.bed"
  logfile = "coord-query-ps-genes-all-gt-"
  #infile = "../../../cds/data/test.bed"
  #logfile = "test-"
  cols = c("chr", "start", "end", "ids", "score", "strand")
}else if(feature=="uces"){
  infile = "../../../uces/mm10-uce/mm10-uce.bed"
  logfile = "coord-query-uces-"
  cols = c("chr", "start", "end", "gid", "score", "strand")
}else if(feature=="hotspots"){
  infile = "../../data/recombination-markers/smagulova-hotspots/Smagulova_2010-hotspots-mm10.bed"
  logfile = "coord-query-hotspots"
  cols = c("chr", "start", "end")
}
logfile = paste(logfile, flank_interval_mb, "mb.log", sep="")


if(read_data){
  cat(as.character(Sys.time()), " | Reading input data\n")
  windows = read.csv(paste("../../data-3/", window_size_str, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts-tt.csv", sep=""), header=T)
  chromes = unique(select(windows, chr, chr.len))
  # The window data with trees.
  
  coords = read.table(infile)
  names(coords) = cols

  #print(nrow(coords))
  #print(head(coords))
  #print(levels(coords$chr))

  rm_chr = c("GL456210.1", "GL456211.1", "GL456212.1", "GL456216.1", "GL456219.1", "GL456221.1", "GL456233.1", "GL456350.1", "GL456354.1", "JH584292.1", "JH584293.1", "JH584294.1", "JH584295.1", "JH584296.1", "JH584297.1", "JH584298.1", "JH584299.1", "JH584303.1", "JH584304.1", "MT", "Y", "chrY", "chr4_GL456216_random", "chrUn_GL456370")

  coords = coords[ ! coords$chr %in% rm_chr, ]

  # coords = subset(coords, chr!="chrY")
  # coords = subset(coords, chr!="chr4_GL456216_random")
  # coords = subset(coords, chr!="chrUn_GL456370")

  coords$chr = factor(coords$chr)

  if(feature=="genes" || feature=="ps" || feature=="coding-genes"){
    coords = coords %>% separate(ids, c("tid", "gid"))
    coords$chr = paste("chr", coords$chr, sep="")
  }else if(feature=="hotspots"){
    coords$gid = paste(coords$chr, ":", coords$start, "-", coords$end, sep="")
  }

  coords = coords[coords$chr %in% windows$chr,]
  coords$int.start = round((coords$start + coords$end) / 2)
  # The coordinates to test in BED format.

  #print(nrow(coords))
  #print(head(coords))
  #print(levels(coords$chr))
  #print(levels(chromes$chr))
  #stop()
}

############################################################
# Main

split_coords = split(coords, f=coords$chr)
# Split by chromosome.

#print(split_coords)

for(chr_coords in split_coords){
  chrome = chr_coords[1,]$chr
  chrome_len = chromes[chromes$chr==chrome,]$chr.len
  
  # if(!chrome == "chr1"){
  #   next
  # }

  # chr_coords = subset(chr_coords, gid=="ENSMUSG00000104010")

  cat(as.character(Sys.time()), " | Starting", chrome, "with", nrow(chr_coords), "features\n")
  
  #outfile = paste("../../data-3/coord-query/", basename(infile), ".", chrome, ".dists", sep="")
  #outfile_random = paste("../../data-3/coord-query/", basename(infile), ".", chrome, ".random", sep="")

  outfile = paste("../../data-3/coord-query/mm10-", feature, "-", flank_interval_mb, "-Mb.bed.", chrome, ".dists", sep="")
  #outfile = "test.dist"
  outfile_random = paste("../../data-3/coord-query/mm10-", feature, "-", flank_interval_mb, "-Mb.bed.", chrome, ".random", sep="")
  
  if(gen_random){
    cat(as.character(Sys.time()), " | ", chrome, " | Selecting random windows\n")
    null_coords = data.frame("chr"=c(), "start"=c(), "end"=c(), "tid"=c(), "gid"=c(), "score"=c(), "strand"=c())
    for(j in 1:nrow(chr_coords)){
      #rand_chr = chromes[sample(1:nrow(chromes), 1), ]
      rand_pos = sample(1:chrome_len, 1)
      null_coords = rbind(null_coords, data.frame("chr"=chrome, "start"=rand_pos, "end"=rand_pos, 
                                                  "tid"=NA, "gid"=NA, "score"=NA, "strand"=NA))
    }
    null_coords$chr = as.character(null_coords$chr)
    null_coords$int.start = null_coords$start
  }
  # The null/random coordinates to compare against. Set FALSE to randomly select 10kb windows

  cat(as.character(Sys.time()), " | ", chrome, " | Calculating tree dists\n")
  #results <- lapply(1:20, parseCoords)
  #registerDoSEQ()
  cl <- makeCluster(num_cores, outfile=logfile)
  registerDoParallel(cl)
  results <- foreach(i=1:nrow(chr_coords), .packages=c("ape", "phangorn")) %dopar% {
    gene = chr_coords[i,]
    print(paste("GENE", i, do.call(paste, c(gene, sep="    "))))
    parseCoords(gene)
  } 
  results = do.call(rbind, results)
  cat(as.character(Sys.time()), " | ", chrome, " | Writing tree dists to file:", outfile, "\n")
  write.csv(results, file=outfile, row.names=F)

  #stopImplicitCluster()
  
  if(do_random){
    cat(as.character(Sys.time()), " | ", chrome, " | Calculating null dists\n")
    null_results <- foreach(i=1:nrow(null_coords), .packages=c("ape", "phangorn")) %dopar% {
      gene = null_coords[i,]
      print(paste("GENE", i, do.call(paste, c(gene, sep="    "))))
      parseCoords(gene)
    } 
    null_results = do.call(rbind, null_results)
    cat(as.character(Sys.time()), " | ", chrome, " | Writing null dists to file:", outfile_random, "\n")
    write.csv(null_results, file=outfile_random, row.names=F)
  }
  stopCluster(cl)
}






