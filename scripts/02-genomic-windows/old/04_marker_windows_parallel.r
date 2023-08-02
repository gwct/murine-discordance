############################################################
# For penn windows, 06.20
# Correlates recombination rate with measures of phylogenetic
# concordance in 5mb windows.
# Gregg Thomas
############################################################

library(ape)
library(phangorn)
library(ggplot2)
library(cowplot)
library(ggsignif)
library("ggtree")
library(parallel)

############################################################
# Functions

calcRFs <- function(chrome_win_data){
  
  windows = data.frame("window"=c(), 
                       "mean.pairwise.wrf"=c(), 
                       "first.last.wrf"=c(),
                       "first.last.dist"=c(),
                       "decay.to.first"=c(),
                       "num.trees"=c(), 
                       "num.uniq.trees"=c(), 
                       "marker.slopes"=c(),
                       "recomb.level"=c()
  )
  # Initialize an empty data frame
  
  marker_windows = unique(chrome_win_data$marker.window)
  # Get the marker windows from the set of tree windows
  
  for(mwin in marker_windows){
    # Loop through every marker window
    print(mwin)
    cur_win = subset(chrome_win_data, marker.window==mwin)
    # Get the current window data
    
    # cur_chr = as.character(unique(cur_win$chr))
    # cur_mean = marker_stats$mean[marker_stats$chr==cur_chr]
    # cur_sd = marker_stats$sd[marker_stats$chr==cur_chr]
    # high_thresh = cur_mean + (2*cur_sd)
    # low_thresh = cur_mean - (2*cur_sd)
    # 
    # cur_median = marker_stats$median[marker_stats$chr==cur_chr]
    # cur_recomb = as.character(unique(cur_win$marker.slope))
    # Recombination stats
    
    num_trees = 0
    
    if(cur_win[1,]$marker.uniq.trees < 10){
      next
    }
    # Skip windows that don't have at least 10 trees
    
    pairwise_wrfs = c()
    wrfs_to_first = c()
    dists_to_first = c()
    first = TRUE
    
    for(i in 1:nrow(cur_win)){
      # Loop through every 100kb tree window within the current marker window
      
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
      
      if(!first){
        dists_to_first = c(dists_to_first, tree_pos-first_tree_pos)
        wrfs_to_first = c(wrfs_to_first, wRF.dist(tree1, first_tree, normalize=T))
      }
      # For decay from first.
      
      
      if(first){
        first_tree = tree1
        first_tree_pos = tree_pos
        first = FALSE
      }
      
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
    }
    
    last_tree = tree1
    last_tree_pos = tree_pos
    
    first_last_wrf = wRF.dist(first_tree, last_tree, normalize=T)
    first_last_dist = last_tree_pos - first_tree_pos
    # For first-to-last distances
    
    first_decay_model = lm(wrfs_to_first ~ log(dists_to_first))
    decay_to_first = first_decay_model$coefficients[2]
    # For decay from first window
    
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
    
    windows = rbind(windows, data.frame("window"=mwin,
                                        "mean.pairwise.wrf"=NA, 
                                        "first.last.wrf"=first_last_wrf,
                                        "first.last.dist"=first_last_dist,
                                        "decay.to.first"=decay_to_first,
                                        "num.trees"=num_trees, 
                                        "num.uniq.trees"=cur_win[1,]$marker.uniq.trees, 
                                        "marker.slopes"=cur_win[1,]$marker.slope,
                                        "recomb.level"=NA))  
  }
  return(windows)
}

############################################################
cat("----------\n")

window_size = 10
# Window size in kb

marker_window_size = 5
# Marker window size in Mb

marker_window_size_kb = marker_window_size * 1000
# Marker window size in kb

au_flag = FALSE
# Set to filter out windows that don't pass the AU test

skip_one = FALSE
# Set to only do one test chromosome

calc_data = T
gen_figs = F
# Which mode to run in

if(gen_figs){
  this.dir <- dirname(parent.frame(2)$ofile)
  setwd(this.dir)
  infile = paste("../../data-3/", window_size, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts-tt.csv", sep="")
  windows = read.csv(paste("../../data-3/", window_size, "kb-", marker_window_size, "mb-dists.csv", sep=""), header=T)
}else{
  infile = paste("../../data-3/", window_size, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts-tt.csv", sep="")
  outfile = paste("../../data-3/", window_size, "kb-", marker_window_size, "mb-dists.csv", sep="")
  #infile = paste("/mnt/beegfs/gt156213e/penn-genomes/windows/data-2/", window_size, "kb-0.5-0.5-topo-counts-tt.csv", sep="")
}
# Input options
######################

cat(as.character(Sys.time()), " | Reading data...\n")
win_data = read.csv(infile, header=T)
# Tree data

win_data = split(win_data, f=win_data$chr)
# Split by chromosome.
  
# Read  and filter the data
######################


if(calc_data){
  #all_windows = lapply(win_data, calcRFs)
  all_windows = mclapply(win_data, calcRFs, mc.cores=20)
  #out_windows = rbind(all_windows)
  out_windows = do.call("rbind", all_windows)
  cat(as.character(Sys.time()), " | Writing:", outfile,"\n")
  write.csv(out_windows, outfile, row.names=F)
}


if(gen_figs){
  cat(as.character(Sys.time()), " | Generating figures...\n")
  
  rf_p = ggplot(windows, aes(x=marker.slopes, y=mean.pairwise.wrf)) +
      geom_point(size=2, color=corecol(numcol=1), alpha=0.5) +
      geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
      ylab("Mean pairwise wRF") +
      xlab("Recombination rate") +
      bartheme()
  print(rf_p)
  
  outfile = paste("../../figs-2/", window_size, "kb-recomb-pairwise-wrf", sep="")
  if(au_flag){
    outfile = paste(outfile, "-au-filter", sep="")
  }
  outfile = paste(outfile, ".png", sep="")
  cat(as.character(Sys.time()), " | Saving figure:", outfile,"\n")
  ggsave(filename=outfile, rf_p, width=5, height=4, units="in")
  # Recombination rate vs pairwise wRF per 5Mb window
  ###########
  
  # norm_rf_p = ggplot(windows, aes(x=marker.slopes, y=mean.pairwise.wrf/num.uniq.trees)) +
  #   geom_point(size=2, color=corecol(numcol=1), alpha=0.5) +
  #   geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
  #   ylab("Mean pairwise wRF\nscaled by # of unique trees") +
  #   xlab("Recombination rate") +
  #   bartheme()
  # 
  # print(norm_rf_p)
  # outfile = "../figs/r-v-pairwise-wrf"
  # if(au_filter){
  #   outfile = paste(outfile, "-au", sep="")
  # }
  # outfile = paste(outfile, "-norm-uniq-trees.png", sep="")
  # if(save_figs){
  #   ggsave(filename=outfile, norm_rf_p, width=5, height=4, units="in")
  # }
  # Recombination rate vs pairwise wRF normalized by number of unique trees per 5Mb window
  ###########
  
  fl_p = ggplot(windows, aes(x=marker.slopes, y=first.last.wrf)) +
    geom_point(size=2, color=corecol(numcol=1), alpha=0.5) +
    geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
    ylab("wRF between first and last tree") +
    xlab("Recombination rate") +
    bartheme()
  print(fl_p)
  
  outfile = paste("../../figs-2/", window_size, "kb-recomb-first-last-wrf", sep="")
  if(au_flag){
    outfile = paste(outfile, "-au-filter", sep="")
  }
  outfile = paste(outfile, ".png", sep="")
  cat(as.character(Sys.time()), " | Saving figure:", outfile,"\n")
  ggsave(filename=outfile, fl_p, width=5, height=4, units="in")
  # Recombination rate vs wRF of first and last tree per 5Mb window
  ###########

  norm_fl_p = ggplot(windows, aes(x=marker.slopes, y=first.last.wrf/first.last.dist)) +
    geom_point(size=2, color=corecol(numcol=1), alpha=0.5) +
    geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
    ylab("wRF between first and last tree\nscaled by distance") +
    xlab("Recombination ratew") +
    bartheme()
  print(norm_fl_p)
  
  outfile = paste("../../figs-2/", window_size, "kb-recomb-first-last-wrf-norm", sep="")
  if(au_flag){
    outfile = paste(outfile, "-au-filter", sep="")
  }
  outfile = paste(outfile, ".png", sep="")
  cat(as.character(Sys.time()), " | Saving figure:", outfile,"\n")
  ggsave(filename=outfile, norm_fl_p, width=5, height=4, units="in")
  # Recombination rate vs wRF of first and last tree scaled by distance per 5Mb window
  ###########
  
  decay_p = ggplot(windows, aes(x=marker.slopes, y=decay.to.first)) +
    geom_point(size=2, color=corecol(numcol=1), alpha=0.5) +
    geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
    ylab("Decay rate from first tree") +
    xlab("Recombination rate") +
    bartheme()
  print(decay_p)
  
  outfile = paste("../../figs-2/", window_size, "kb-recomb-decay", sep="")
  if(au_flag){
    outfile = paste(outfile, "-au-filter", sep="")
  }
  outfile = paste(outfile, ".png", sep="")
  cat(as.character(Sys.time()), " | Saving figure:", outfile,"\n")
  ggsave(filename=outfile, decay_p, width=5, height=4, units="in")
  # Recombination rate vs decay rate from first tree per 5Mb window
  ###########
  stop("OK")
  
  
  
  
  windows_binned = subset(windows, !is.na(recomb.level))
  binned_p = ggplot(windows_binned, aes(x=recomb.level, y=first.last.wrf/first.last.dist, fille=recomb.level)) +
    geom_boxplot(alpha=0.7, width=0.5) +
    #geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
    ylab("wRF between first and last tree\nscaled by distance") +
    xlab("Recombination level") +
    bartheme()
  
  print(binned_p)
  outfile = "../figs/r-v-first-last-binned"
  if(au_filter){
    outfile = paste(outfile, "-au", sep="")
  }
  outfile = paste(outfile, ".png", sep="")
  if(save_figs){
    ggsave(filename=outfile, binned_p, width=5, height=4, units="in")
  }
  # Distance to non-random wRF binned by high or low recombination rate.
  
  
}

# R <- function(a, b, abT) a*(1 - exp(-b*abT))
# form <- Richness ~ R(a,b,Abundance)
# fit <- nls(form, data=d, start=list(a=20,b=0.01))
# plot(d$Abundance,d$Richness, xlab="Abundance", ylab="Richness")
# lines(d$Abundance, predict(fit,list(x=d$Abundance)))


    
#     
#     
#     
#     
#     cur_chrome = as.character(cur_win$chr)[1]
#     if(au_filter){
#       chrdata_f = subset(win_data, chr==cur_chrome & status=="PASS" & AU.topology.test=="PASS")
#     }else{
#       chrdata_f = subset(win_data, chr==cur_chrome & status=="PASS")
#     }
#     # Get data for the current chromosome of the marker window
# 
#     concat_chr_tree = read.tree(text=as.character(chrdata_f$unparsed.tree[chrdata_f$concat.chrome.topo][1]))
#     astral_chr_tree = read.tree(text=as.character(chrdata_f$unparsed.tree[chrdata_f$astral.chrome.topo][1]))
#     # Get the chromosome trees
#     
#     topo_counts = unique(subset(cur_win, select=c(topo.num.chrome, topo.count.chrome, topo.num.overall, topo.rank.overall)))
#     topo_counts = subset(topo_counts, !is.na(topo.num.chrome))
#     # Get the topology counts
#     
#     topo_counts = topo_counts[order(topo_counts$topo.count.chrome, decreasing=TRUE),]
#     # Order the rows by topology count
#     
#     topo_counts$curwin.count = NA
#     topo_counts$chrome.spr = NA
#     # For counting current trees in marker window and calculating SPR to chromosome tree
#     for(i in topo_counts$topo.num.chrome){
#       cur_tree_bool = cur_win$unparsed.tree[cur_win$topo.num.chrome==i]
#       cur_tree_bool = cur_tree_bool[!is.na(cur_tree_bool)];
#       cur_tree = read.tree(text=as.character(cur_tree_bool[1]))
#       #topo_counts$unparsed.tree[topo_counts$topo.num.chrome==i] = cur_tree
#       # Get each tree in the marker window, 
#       
#       cur_spr = SPR.dist(cur_tree, concat_chr_tree)
#       topo_counts$chrome.spr[topo_counts$topo.num.chrome==i] = cur_spr
#       # Calculate the SPR between the current marker tree and the chromosome tree
# 
#       cur_counts = cur_win$topo.num.chrome
#       cur_counts = cur_counts[!is.na(cur_counts)]
#       topo_counts$curwin.count[topo_counts$topo.num.chrome==i] = sum(cur_counts == i)
#       # Count how many times this tree appears in the current marker window
#     }
#     # Add label and color columns to the topo counts
#     
#     num_cur_trees = sum(topo_counts$curwin.count)
#     # Get the total number of trees in the current marker window
#     
#     topo_counts$curwin.frac = topo_counts$curwin.count / num_cur_trees
#     topo_counts$weighted.spr = topo_counts$chrome.spr * topo_counts$curwin.frac
#     # Calculate the SPR for each tree weighted by the fraction of trees that are that topology
#     
#     mean_weighted_spr = mean(topo_counts$weighted.spr)
#     # Take the mean of the weighted SPRs for this marker window
#     
#     pairwise_rfs = c()
#     pairwise_spr = c()
#     num_trees = 0
#     # A vector for pairwise RFs and a tree counter
#     
#     for(i in 1:nrow(cur_win)){
#     # Loop through every 100kb tree window within the current marker window
#       
#       tree_row1 = cur_win[i,]
#       tree_win1 = tree_row1$window
#       # Get the current tree window
#   
#       tree_str1 = as.character(tree_row1$unparsed.tree)
#       # Get the current tree string
#       
#       if(is.na(tree_str1)){
#         next
#       }
#       # If the tree is empty then the window was filtered and we can skip this window
#       
#       tree1 = read.tree(text=tree_str1)
#       # Read the tree as a phylo object
#       
#       num_trees = num_trees + 1
#       # Add to the number of trees counted
#       
#       for(j in 1:nrow(cur_win)){
#       # For every 100kb window, we check every other 100kb window
#         tree_row2 = cur_win[j,]
#         tree_win2 = tree_row2$window
#         # Get the tree window to check
#         
#         if(tree_win1 == tree_win2){
#           next
#         }
#         # If the windows are the same, move one
#         
#         tree_str2 = as.character(tree_row2$unparsed.tree)
#         # Get the current tree string
#         
#         if(is.na(tree_str2)){
#           next
#         }
#         # If the tree is empty then the window was filtered and we can skip this window
#         
#         tree2 = read.tree(text=tree_str2)
#         # Read the tree as a phylo object
#         
#         pairwise_rfs = c(pairwise_rfs, wRF.dist(tree1, tree2))
#         pairwise_spr = c(pairwise_spr, SPR.dist(tree1, tree2))
#         # Calculate tree distance between the two trees and add to vector of distances for this marker window
#       }
#     }
#   
#     if(is.null(pairwise_rfs)){
#       windows = rbind(windows, data.frame("window"=mwin,
#                                           "mean.pairwise.rf"=NA, 
#                                           "mean.pairwise.spr"=NA, 
#                                           "mean.weighted.spr"=NA,
#                                           "num.trees"=num_trees, 
#                                           "num.uniq.trees"=cur_win[1,]$marker.uniq.trees, 
#                                           "marker.slopes"=cur_win[1,]$marker.slope)
#                                           )
#     }else{
#       windows = rbind(windows, data.frame("window"=mwin, 
#                                           "mean.pairwise.rf"=mean(pairwise_rfs), 
#                                           "mean.pairwise.spr"=mean(pairwise_spr),
#                                           "mean.weighted.spr"=mean_weighted_spr,
#                                           "num.trees"=num_trees, 
#                                           "num.uniq.trees"=cur_win[1,]$marker.uniq.trees, 
#                                           "marker.slopes"=cur_win[1,]$marker.slope)
#                                           )  
#     }
#     # Calculate the mean pairwise tree distance for all trees in this window and add to the output data frame
#   }
# 
#   write.csv(windows, file="../data/pairwise-treedists.csv")
#   # Write the pairwise distance data to a file
# }
# 
# ###############
# # RF plots
# rf_p = ggplot(windows, aes(x=marker.slopes, y=mean.pairwise.rf)) +
#   geom_point(size=2, color=corecol(numcol=1), alpha=0.5) +
#   geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
#   ylab("Mean pairwise RF distances between\nall 100kb treesin 5mb window") +
#   xlab("Recombination rate per 5mb window") +
#   bartheme()
# 
# print(rf_p)
# # Plot raw pairwise tree distance vs. recomb rate
# 
# rf_norm_p = ggplot(windows, aes(x=marker.slopes, y=mean.pairwise.rf/num.uniq.trees)) +
#   geom_point(size=2, color=corecol(numcol=1), alpha=0.5) +
#   geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
#   ylab("Mean pairwise RF distances between\nall 100kb trees in 5mb window / the\nnumber of unique 100kb trees in 5mb window") +
#   xlab("Recombination rate per 5mb window") +
#   bartheme()
# 
# print(rf_norm_p)
# # Plot normalized (by number of unique trees in marker windows) pairwise tree distance vs. recomb rate
# 
# rf_pcombo = plot_grid(rf_p, rf_norm_p, nrow=1, labels=c("A","B"), label_size=20, align="v")
# # Combine the tree distance plots
# 
# 
# ###############
# # SPR plots
# 
# spr_p = ggplot(windows, aes(x=marker.slopes, y=mean.pairwise.spr)) +
#   geom_point(size=2, color=corecol(numcol=1), alpha=0.5) +
#   geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
#   ylab("Mean pairwise SPR distances between\nall 100kb trees in 5mb window") +
#   xlab("Recombination rate per 5mb window") +
#   bartheme()
# 
# print(spr_p)
# # Plot raw pairwise tree distance vs. recomb rate
# 
# spr_norm_p = ggplot(windows, aes(x=marker.slopes, y=mean.pairwise.spr/num.uniq.trees)) +
#   geom_point(size=2, color=corecol(numcol=1), alpha=0.5) +
#   geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
#   ylab("Mean pairwise SPR distances between\nall 100kb trees in 5mb window / the\nnumber of unique 100kb trees in 5mb window") +
#   xlab("Recombination rate per 5mb window") +
#   bartheme()
# 
# print(spr_norm_p)
# # Plot normalized (by number of unique trees in marker windows) pairwise tree distance vs. recomb rate
# 
# mw_spr_p = ggplot(windows, aes(x=marker.slopes, y=mean.weighted.spr)) +
#   geom_point(size=2, color=corecol(numcol=1), alpha=0.5) +
#   geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
#   ylab("Mean weighted SPR distances between 100kb trees\nand concatenated chromosome tree in 5mb window") +
#   xlab("Recombination rate per 5mb window") +
#   bartheme()
# 
# print(mw_spr_p)
# # Plot raw pairwise tree distance vs. recomb rate
# 
# 
# spr_pcombo = plot_grid(spr_p, spr_norm_p, nrow=1, labels=c("C","D","E"), label_size=20, align="v")
# # Combine the tree distance plots
# 
# ###############
# # NUT plots
# 
# nut_p = ggplot(windows, aes(x=marker.slopes, y=num.uniq.trees)) +
#   geom_point(size=2, color=corecol(numcol=1, offset=1), alpha=0.5) +
#   geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
#   ylab("The number of unique 100kb trees in 5mb window") +
#   xlab("Recombination rate per 5mb window") +
#   bartheme()
# print(nut_p)
# # Plot the number of unique trees per marker window by the recomb rate
# 
# nut_norm_p = ggplot(windows, aes(x=marker.slopes, y=num.uniq.trees/num.trees)) +
#   geom_point(size=2, color=corecol(numcol=1, offset=1), alpha=0.5) +
#   geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
#   ylab("The number of unique 100kb trees in 5mb window / the\ntotal number of 100kb trees in 5mb window") +
#   xlab("Recombination rate per 5mb window") +
#   bartheme()
# print(nut_norm_p)
# # Plot the normalized (by total number of trees) number of unique trees per marker window by the recomb rate
# 
# nut_pcombo = plot_grid(nut_p, nut_norm_p, nrow=1, labels=c("F","G"), label_size=20, align="v")
# # Combine the tree count plots
# 
# # gcf_p = ggplot(windows, aes(x=marker.slopes, y=avg.gcf)) +
# #   geom_point(size=2, color=corecol(numcol=1, offset=1), alpha=0.5) +
# #   geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
# #   ylab("Average locus concordance factor to concatenated\nchromosome tree among 100kb trees in 5mb window") +
# #   xlab("Recombination rate per 5mb window") +
# #   bartheme()
# # print(nut_p)
# # 
# # scf_p = ggplot(windows, aes(x=marker.slopes, y=avg.scf)) +
# #   geom_point(size=2, color=corecol(numcol=1, offset=1), alpha=0.5) +
# #   geom_smooth(method="lm", color="#333333", linetype="dashed", se=F) +
# #   ylab("Average site concordance factor to concatenated\nchromosome tree among 100kb trees in 5mb window") +
# #   xlab("Recombination rate per 5mb window") +
# #   bartheme()
# # print(nut_p)
# # 
# # cf_pcombo = plot_grid(gcf_p, scf_p, nrow=1, labels=c("E","F"), label_size=20, align="v")
# 
# ###############
# # Combine and save plots
# 
# p = plot_grid(rf_p, rf_norm_p, spr_p, spr_norm_p, mw_spr_p, nut_p, nut_norm_p, ncol=2, nrow=4, align="v", labes=c("A","B","C","D","E","F","G"), label_size=20)
# 
# # p = plot_grid(rf_pcombo, spr_pcombo, nut_pcombo, mw_spr_p, nrow=4, align="v")
# ggsave(filename="../figs/rcomb-v-treedist-tt.png", p, width=14, height=15, units="in")
# # Combine plots and save
