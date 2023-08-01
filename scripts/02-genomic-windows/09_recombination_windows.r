############################################################
# For penn windows, 06.20
# Correlates recombination rate with measures of phylogenetic
# concordance in 5mb windows.
# Gregg Thomas
############################################################

library(ape)
library(phangorn)
#library(ggplot2)
#library(cowplot)
#library(ggsignif)
#library("ggtree")
library(parallel)

############################################################
# Functions

calcRFs <- function(chrome_win_data){

  chrome = as.character(chrome_win_data$chr[1])
  if(is.na(chrome)){
    return(NULL)
  }
  # cat(chrome)
  # cat(nchar(chrome))
  if(nchar(chrome) == 4 && chrome != "chrX"){
    out_chrome = gsub("chr", "chr0", chrome)
  }else{
    out_chrome = chrome
  }
  # The string for the current chromosome.
  
  cat(as.character(Sys.time()), " | Starting chromosome:", out_chrome, "\n")
  
  windows = data.frame("window"=c(), "num.trees"=c(), "num.uniq.trees"=c(), "recomb.rate"=c(), "first.last.dist"=c(), 
                       "first.last.rf"=c(), "rf.decay"=c(), 
                       "first.last.wrf"=c(), "wrf.decay"=c(), 
                       "first.last.spr"=c(), "spr.decay"=c(), 
                       "first.last.kf"=c(), "kf.decay"=c(), 
                       "first.last.path"=c(), "path.decay"=c()
  )
  # Initialize an empty data frame
  
  marker_windows = unique(chrome_win_data$marker.window)
  # Get the marker windows from the set of tree windows
  
  for(mwin in marker_windows){
    # Loop through every marker window
    
    cat(as.character(Sys.time()), "| ->", mwin, "\n")

    cur_win = subset(chrome_win_data, marker.window==mwin)
    # Get the current window data
    
    num_trees = 0
    
    if(cur_win[1,]$marker.uniq.trees < 10){
      next
    }
    # Skip windows that don't have at least 10 trees
    
    cur_win_data = data.frame("dist.to.first"=c(), "rf.to.first"=c(), "wrf.to.first"=c(),
                              "spr.to.first"=c(), "kf.to.first"=c(), "path.to.first"=c())
    # Data frame for current marker window
    
    first = TRUE
    # Flag indicating this the first window
    for(i in 1:nrow(cur_win)){
      # Loop through every 10kb tree window within the current marker window
      
      tree_row1 = cur_win[i,]
      tree_win1 = tree_row1$window
      # Get the current tree window
      
      tree_str1 = as.character(tree_row1$unparsed.tree)
      # Get the current tree string
      
      if(is.na(tree_str1)){
        next
      }
      # If the tree is empty then the window was filtered and we can skip this window
      
      tree1 = read.tree(text=tree_str1)
      # Read the tree as a phylo object
      
      tree_pos = tree_row1$start
      # Get the start position of the current tree to calculate the distance to the
      # first tree
      
      if(first){
        first_tree = tree1
        first_tree_pos = tree_pos
        first = FALSE
        next
      }
      # If this is the first tree of the window, save the tree and the start position
      # for all other comparisons
      
      
      if(!first){
        dist_to_first = tree_pos - first_tree_pos
        # The physical distance to the first tree
        
        rf = RF.dist(tree1, first_tree)
        wrf = wRF.dist(tree1, first_tree)
        spr = SPR.dist(tree1, first_tree)
        names(spr) = NULL
        spr = spr
        kf = KF.dist(tree1, first_tree)
        path = path.dist(tree1, first_tree)
        # The phylogenetic distances to the first tree
        
        cur_win_data = rbind(cur_win_data, data.frame("dist.to.first"=dist_to_first, "rf.to.first"=rf, "wrf.to.first"=wrf,
                                                      "spr.to.first"=spr, "kf.to.first"=kf, "path.to.first"=path))
        # Add the distances to the data frame
      }
      # For all other trees, calculate distances
      
      num_trees = num_trees + 1
      # Add to the number of trees counted
      
      # for(j in 1:nrow(cur_win)){
      #   # For every 100kb window, we check every other 100kb window
      #   tree_row2 = cur_win[j,]
      #   tree_win2 = tree_row2$window
      #   # Get the tree window to check
      #   
      #   if(tree_win1 == tree_win2){
      #     next
      #   }
      #   # If the windows are the same, move one
      #   
      #   tree_str2 = as.character(tree_row2$unparsed.tree)
      #   # Get the current tree string
      #   
      #   if(is.na(tree_str2)){
      #     next
      #   }
      #   # If the tree is empty then the window was filtered and we can skip this window
      #   
      #   tree2 = read.tree(text=tree_str2)
      #   # Read the tree as a phylo object
      #   
      #   pairwise_wrfs = c(pairwise_wrfs, wRF.dist(tree1, tree2, normalize=T))
      #   #pairwise_spr = c(pairwise_spr, SPR.dist(tree1, tree2))
      #   # Calculate tree distance between the two trees and add to vector of distances for this marker window
      # }
      # This was for the pairwise distances, but that takes a really long time.
    }
    
    # last_tree = tree1
    # last_tree_pos = tree_pos
    # 
    # first_last_wrf = wRF.dist(first_tree, last_tree, normalize=T)
    # first_last_dist = last_tree_pos - first_tree_pos
    # # For first-to-last distances
    
    # first_decay_model = lm(wrfs_to_first ~ log10(dists_to_first))
    # decay_to_first = first_decay_model$coefficients[2]
    # For decay from first window
    
    if(num_trees > 1){
      rf_decay = lm(cur_win_data$rf.to.first ~ log10(cur_win_data$dist.to.first))$coefficients[2]
      wrf_decay = lm(cur_win_data$wrf.to.first ~ log10(cur_win_data$dist.to.first))$coefficients[2]
      spr_decay = lm(cur_win_data$spr.to.first ~ log10(cur_win_data$dist.to.first))$coefficients[2]
      kf_decay = lm(cur_win_data$kf.to.first ~ log10(cur_win_data$dist.to.first))$coefficients[2]
      path_decay = lm(cur_win_data$path.to.first ~ log10(cur_win_data$dist.to.first))$coefficients[2]
    }else{
      rf_decay = NA
      wrf_decay = NA
      spr_decay = NA
      kf_decay = NA
      path_decay = NA     
    }
    
    # plot(dists_to_first, wrfs_to_first)
    # lines(dists_to_first,predict(first_decay_model),col='red')
    # summary(first_decay_model)
    # print(first_decay_model$coefficients)
    # 
    # if(is.na(cur_recomb)){
    #   rec_level = NA
    # }else if(cur_recomb >= high_thresh){
    #   rec_level = "high"
    # }else if(cur_recomb <= low_thresh){
    #   rec_level = "low"
    # }
    # Getting whether recombination rate is high or low for this chromosome
    
    windows = rbind(windows, 
                          data.frame("window"=mwin, "num.trees"=num_trees, "num.uniq.trees"=cur_win[1,]$marker.uniq.trees,
                         "recomb.rate"=cur_win[1,]$marker.slope, "first.last.dist"=dist_to_first, 
                         "first.last.rf"=rf, "rf.decay"=rf_decay, 
                         "first.last.wrf"=wrf, "wrf.decay"=wrf_decay, 
                         "first.last.spr"=spr, "spr.decay"=spr_decay, 
                         "first.last.kf"=kf, "kf.decay"=kf_decay, 
                         "first.last.path"=path, "path.decay"=path_decay)
    )
    # Add the data for the current window to the main data frame
  }
  return(windows)
}

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

num_cores = 20
# The number of cores to split chromosomes across

test_run = FALSE
# Only do chromosome 7 as a test (~1hr)

if(test_run){
  # this.dir <- dirname(parent.frame(2)$ofile)
  # setwd(this.dir)
  num_cores = 1
  # Number of cores
}

infile = paste("../data/", window_size, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts.csv", sep="")
outfile = paste("../data/", window_size, "kb-", marker_window_size, "mb-recomb-dists.csv", sep="")
#infile = paste("/mnt/beegfs/gt156213e/penn-genomes/windows/data-2/", window_size, "kb-0.5-0.5-topo-counts-tt.csv", sep="")

# Input options
######################

cat(as.character(Sys.time()), " | Reading data:", infile, "\n")
win_data = read.csv(infile, header=T, comment.char="#")
# Tree data

if(test_run){
  cat(as.character(Sys.time()), " | Subsetting data for test run.\n")
  test_windows = subset(win_data, chr=="chr19")
  outfile = "marker-windows-test.csv"
}
# Get only one chrome for a test run

# Read  and filter the data
######################

if(test_run){
  test_windows = split(test_windows, f=test_windows$chr)
  all_windows = mclapply(test_windows, calcRFs, mc.cores=num_cores)
}else{
  win_data = split(win_data, f=win_data$chr)
  # Split by chromosome.
  
  all_windows = mclapply(win_data, calcRFs, mc.cores=num_cores)
}
out_windows = do.call("rbind", all_windows)
cat(as.character(Sys.time()), " | Writing output:", outfile,"\n")
write.csv(out_windows, outfile, row.names=F)
# Calc dists in parallel for each chromosome

