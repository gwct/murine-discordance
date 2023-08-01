############################################################
# For penn genomes, 04.20
# Checks distributions of genes overlapping CNVs
# Gregg Thomas
############################################################

library(ape)
library(phangorn)
#library(ggplot2)
#library(cowplot)
#library(ggbeeswarm)
library(dplyr)
library(parallel)
library(here)

#############################################################
# Functions

getRandWindow <- function(window_df){
  rand_window = data.frame("window"="", "chr"="", "start"="", "tree"="")
  
  random_window = as.character(sample(window_df$window, 1))
  rand_window$window = random_window
  
  window_data = window_df[window_df$window==random_window,]
  rand_window$chr = window_data$chr
  rand_window$start = window_data$start
  rand_window$tree = as.character(window_data$unparsed.tree)
  #rand_window$tree = read.tree(text=as.character(window_data$unparsed.tree))
  
  return(rand_window)
}

######################
######################

distReplicates <- function(chrdata, all_win_f, m_win_size, rdists, outdir, numreps){
  
  chrome = as.character(chrdata[1,]$chr)
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
  
  cat(as.character(Sys.time()), " | Starting chromosome ", out_chrome, "\n", sep="")
  
  #chrdata = subset(all_windows, chr==chrome)
  total_windows = length(chrdata[,1])
  chrdata_f = subset(chrdata, repeat.filter=="PASS" & missing.filter=="PASS")
  if(au_flag){
    chrdata_f = subset(chrdata_f, AU.test=="PASS")
  }
  used_windows = length(chrdata_f[,1])
  chr_len = chrdata$chr.len[1]
  # Subset data for this chromosome.
  
  #outfile = paste(outdir, out_chrome, sep="")
  outfile = here(outdir, out_chrome)
  #statsfile = paste(outdir, out_chrome, "-stats", sep="") 
  statsfile = here(outdir, paste0(out_chrome, "-stats"))  
  if(au_flag){
    outfile = paste0(outfile, "-au-filter")
    statsfile = paste0(statsfile, "-au-filter")
  }
  outfile = paste0(outfile, ".csv")
  statsfile = paste0(statsfile, ".csv")
  # Set-up output files for this chromosome.
  
  out_data = data.frame("replicate"=c(), "measure"=c(), "nonsig.adj"=c())
  # Output data frame for this chromosome
  
  #stats_data = data.frame("sample"=c(), "adj"=c(),"num.dists"=c(), "mean.wrf"=c(), "median.wrf"=c(), "var.wrf"=c())
  stats_data = data.frame("replicate"=c(), "adj"=c(), "num.dists"=c(), 
                          "mean.rf"=c(), "median.rf"=c(), "var.rf"=c(),
                          "mean.wrf"=c(), "median.wrf"=c(), "var.wrf"=c(),
                          "mean.spr"=c(), "median.spr"=c(), "var.spr"=c(),
                          "mean.kf"=c(), "median.kf"=c(), "var.kf"=c(),
                          "mean.path"=c(), "median.path"=c(), "var.path"=c())
  # Stats per adjacency and replicate for this chromosome.

  
  sig = TRUE
  sig_rf = seq(1:numreps)
  sig_wrf = seq(1:numreps)
  sig_spr = seq(1:numreps)
  sig_kf = seq(1:numreps)
  sig_path = seq(1:numreps)
  cur_adj = 1
    
  while(sig || cur_adj <= 500){
    # Calculate wRF for adjacent windows until the distribution becomes random-like (non-significant for wilcox or KS test)  
    cat(as.character(Sys.time()), " | Chromosome ", out_chrome, " - adjacency ", cur_adj, "\n", sep="")
    
    cur_dists = data.frame("win1"=c(), "win2"=c(), "rf"=c(), "wrf"=c(), "spr"=c(), "kf"=c(), "path"=c())
    # The tree distances for the current adjacency level
    
    dists_file_base = paste0(chrome, "-", cur_adj, ".csv")
    dists_file = here("data", "dists", "all-dists", dists_file_base)
    # The name of the file to output all tree distances for the current adjacency
    
    if(file.exists(dists_file)){
      cur_dists = read.csv(dists_file, header=T)
      # Check if the dists file already exists and read it if it does, otherwise proceed to calculate dists below
    }else{
      for(w in 1:nrow(chrdata)) {
        # For every window on the chromosome
        
        window = chrdata[w,]
        cur_window = as.character(window$window)
        window_start = window$start
        # Get the info for the current window
        
        if(window$repeat.filter!="PASS" || window$missing.filter!="PASS"){
          next
        }
        if(au_flag && window$AU.test != "PASS"){
          next
        }
        # Skip the window if it didn't pass one of the filters.
        
        cur_topo = window$topo.num.chrome
        cur_tree = read.tree(text=as.character(window$unparsed.tree))
        # Get the tree for the current window
        
        cur_adj_start = window_start + (cur_adj * window_size)
        # Sample the window at the current adjacency.
        
        if(cur_adj_start > chr_len){
          next
        }
        # Flag whether the forward start exceeds the chrome length
        
        forward_window = chrdata[chrdata$start==cur_adj_start,]
        # Get the next window based on the current adjaceny step and window size
        
        if(forward_window$repeat.filter!="PASS" || forward_window$missing.filter!="PASS"){
          next
        }
        if(au_flag && forward_window$AU.test != "PASS"){
          next
        }
        # Skip the adjacent window if it didn't pass one of the filters.
        
        forward_tree = read.tree(text=as.character(forward_window$unparsed.tree))
        # Get the tree for the adjacent window
        
        rf = RF.dist(cur_tree, forward_tree)
        wrf = wRF.dist(cur_tree, forward_tree)
        spr = SPR.dist(cur_tree, forward_tree)
        names(spr) = NULL
        spr = spr
        kf = KF.dist(cur_tree, forward_tree)
        path = path.dist(cur_tree, forward_tree)
        # Calculate the wRF for the two trees   
        
        cur_dists = rbind(cur_dists, data.frame("win1"=cur_window, "win2"=as.character(forward_window$window), "rf"=rf, "wrf"=wrf, "spr"=spr, "kf"=kf, "path"=path))
        # Calculate the wRF for the two trees
        # The forward distance.       
      }
      # Tree adjacency loop
    }
    # cur_dists check
    
    for(i in 1:numreps){  
      rdists_rep = subset(rdists, replicate==i)
      # For every replicate, get the random distance distribution to test the current distribution against below

      if(sig_rf[i] != "DONE"){
        wilcox_rf_p = wilcox.test(cur_dists$rf, rdists_rep$rf)$p.value
        if(wilcox_rf_p > 0.01){
          sig_rf[i] = "DONE"
          out_data = rbind(out_data, data.frame("replicate"=i, "measure"="rf", "nonsig.adj"=cur_adj))
        }
      }
      # RF test
      
      if(sig_wrf[i] != "DONE"){
        wilcox_wrf_p = wilcox.test(cur_dists$wrf, rdists_rep$wrf)$p.value
        if(wilcox_wrf_p > 0.01){
          sig_wrf[i] = "DONE"
          out_data = rbind(out_data, data.frame("replicate"=i, "measure"="wrf", "nonsig.adj"=cur_adj))
        }
      }
      # wRF test
      
      if(sig_spr[i] != "DONE"){
        wilcox_spr_p = wilcox.test(cur_dists$spr, rdists_rep$spr)$p.value
        if(wilcox_spr_p > 0.01){
          sig_spr[i] = "DONE"
          out_data = rbind(out_data, data.frame("replicate"=i, "measure"="spr", "nonsig.adj"=cur_adj))
        }
      }
      # SPR test     
      
      if(sig_kf[i] != "DONE"){
        wilcox_kf_p = wilcox.test(cur_dists$kf, rdists_rep$kf)$p.value
        if(wilcox_kf_p > 0.01){
          sig_kf[i] = "DONE"
          out_data = rbind(out_data, data.frame("replicate"=i, "measure"="kf", "nonsig.adj"=cur_adj))
        }
      }
      # KF test
      
      if(sig_path[i] != "DONE"){
        wilcox_path_p = wilcox.test(cur_dists$path, rdists_rep$path)$p.value
        if(wilcox_path_p > 0.01){
          sig_path[i] = "DONE"
          out_data = rbind(out_data, data.frame("replicate"=i, "measure"="path", "nonsig.adj"=cur_adj))
        }
      }
      # path test
    }
    # Replicate loop
    
    if(length(unique(sig_rf)) == 1 && length(unique(sig_wrf)) == 1 && length(unique(sig_spr)) == 1 && length(unique(sig_kf)) == 1 && length(unique(sig_path)) == 1){
      sig = FALSE
    }
    # Exit loop when all measures are not significantly different from random.
    
    write.csv(cur_dists, dists_file, row.names=F)
    # Write the distances for the current adjacency
    
    stats_data = rbind(stats_data, data.frame("adj"=cur_adj, "num.dists"=nrow(cur_dists), 
                                              "mean.rf"=mean(cur_dists$rf), "median.rf"=median(cur_dists$rf), "var.rf"=var(cur_dists$rf),
                                              "mean.wrf"=mean(cur_dists$wrf), "median.wrf"=median(cur_dists$wrf), "var.wrf"=var(cur_dists$wrf),
                                              "mean.spr"=mean(cur_dists$spr), "median.spr"=median(cur_dists$spr), "var.spr"=var(cur_dists$spr),
                                              "mean.kf"=mean(cur_dists$kf), "median.kf"=median(cur_dists$kf), "var.kf"=var(cur_dists$kf),
                                              "mean.path"=mean(cur_dists$path), "median.path"=median(cur_dists$path), "var.path"=var(cur_dists$path))
    )
    # Compile output stats for the current adjacency
    
    cur_adj = cur_adj + 1
    # Iterate adjacency step
  }
  ## Sig loop
  
  cat(as.character(Sys.time()), " | Chromosome ", out_chrome, " - writing data:       ", outfile, "\n", sep="")
  write.csv(out_data, file=outfile, row.names=F)
  cat(as.character(Sys.time()), " | Chromosome ", out_chrome, " - writing stats:      ", statsfile, "\n", sep="")
  write.csv(stats_data, file=statsfile, row.names=F)
  # Output final stats to files
}

############################################################

window_size_kb = 10
# Window size in kb

window_size = window_size_kb * 1000
# Window size in bp

marker_window_size = 5
# Marker window size in Mb

num_cores = 20
# Number of cores

au_flag = FALSE
# Set to filter out windows that don't pass the AU test

skip_one = FALSE
# Set to only do one test chromosome

nrand_win = 12000
# The number of random windows to sample

nreps = 10
# The number of replicates

test_run = TRUE
# Only do chromosome 7 as a test (~1hr)

if(test_run){
  #this.dir <- dirname(parent.frame(2)$ofile)
  #setwd(this.dir)
  num_cores = 1
  # Number of cores
}

#infile = paste("../data/", window_size_kb, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts.csv", sep="")
infile = here("data", paste0(window_size_kb, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts.csv"))

#script_dir = "../data/dists/"
script_dir = here("data", "dists")
#random_file = paste0(script_dir, "random-dists")
random_file = here(script_dir, "random-dists")
if(au_flag){
  random_file = paste0(random_file, "-au")
}
random_file = paste0(random_file, ".csv")
# Input options
######################

cat(as.character(Sys.time()), " | Reading data:      ", infile, "\n", sep="")
all_windows = read.csv(infile, header=T, comment.char="#")
#all_windows = all_windows[1:1000, ]
all_windows_f = subset(all_windows, repeat.filter=="PASS" & missing.filter=="PASS")
if(au_flag){
  all_windows_f = subset(all_windows_f, AU.test=="PASS")
}

if(test_run){
  test_windows = subset(all_windows, chr=="chr19")
  test_windows_f = subset(all_windows, repeat.filter=="PASS" & missing.filter=="PASS")
  if(au_flag){
    test_windows_f = subset(test_windows_f, AU.test=="PASS")
  }
}
# Read and filter the data
######################
cat("----------\n")

random_dists = data.frame("replicate"=c(), "win1"=c(), "win2"=c(), 
                          "rf"=c(), "wrf"=c(), "spr"=c(), "kf"=c(), "path"=c())
# The wRF distances for the random pairs of windows.

if(file.exists(random_file)){
  cat(as.character(Sys.time()), " | Reading random distributions file: ", random_file, "\n", sep="")
  random_dists = read.csv(random_file, header=T)
  # Read the random distributions if they already exist
}else{
  for(i in 1:nreps){
    cat(as.character(Sys.time()), " | Replicate ", i, " - Random sampling\n", sep="")
    # Replicate this procedure N times
    
    for(j in 1:nrand_win){
      # Sample a random pair of windows for each window on the current chromosome that passes all filters.
      
      random_window1 = getRandWindow(all_windows_f)
      random_window2 = random_window1
      while(random_window2$chr == random_window1$chr){
        random_window2 = getRandWindow(all_windows_f)
      }
      # Select the two random windows, making sure they are not on the same chromosome.
      
      random_tree1 = read.tree(text=as.character(random_window1$tree))
      random_tree2 = read.tree(text=as.character(random_window2$tree))
      # Read the trees from the two random windows as phylo objects.
      
      rf = RF.dist(random_tree1, random_tree2)
      wrf = wRF.dist(random_tree1, random_tree2)
      spr = SPR.dist(random_tree1, random_tree2)
      names(spr) = NULL
      spr = spr
      kf = KF.dist(random_tree1, random_tree2)
      path = path.dist(random_tree1, random_tree2)
      # Calculate the wRF for the two trees    
      
      random_dists = rbind(random_dists, data.frame("replicate"=i, "win1"=random_window1$window, "win2"=random_window2$window, 
                                                    "rf"=rf, "wrf"=wrf, "spr"=spr, "kf"=kf, "path"=path))
      # Caclulate the current random wRF.
    }
  }
  cat(as.character(Sys.time()), " | Writing random data:       ", random_file, "\n", sep="")
  write.csv(random_dists, file=random_file, row.names=F)
  # Write out random dists
}
## The random sampling

all_windows = subset(all_windows, chr!="chrX")
# Remove the X because it's weird

if(test_run){
  test_windows = split(test_windows, f=test_windows$chr)
  mclapply(test_windows, distReplicates, all_win_f=test_windows_f, m_win_size=marker_window_size, rdists=random_dists, outdir=script_dir, numreps=nreps, mc.cores=num_cores)
  
}else {
  all_windows = split(all_windows, f=all_windows$chr)
  # Split by chromosome.
  
  mclapply(all_windows, distReplicates, all_win_f=all_windows_f, m_win_size=marker_window_size, rdists=random_dists, outdir=script_dir, numreps=nreps, mc.cores=num_cores)
}

######
