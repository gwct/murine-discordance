############################################################
# For penn genomes, 01.22
# Generates chromoplots for each chromosome
# Gregg Thomas
############################################################

library(ggplot2)
library(cowplot)
library(ggsignif)
library("ggtree")
library(phangorn)
library(zoo)
library(here)
source(here("lib", "get_tree_info.r"))
source(here("lib", "design.r"))

############################################################
cat("----------\n")

window_size = 10
# Window size in kb

marker_window_size = 1
# Marker window size in Mb

au_flag = F
# Set to filter out windows that don't pass the AU test

read_data = T
# Whether to read the initial data or not

save_fig = F
# Whether or not to save the figure

skip_one = F
# Set to only do one test chromosome

chr_order = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX")
# Correct order of chromosomes

datadir = here("data", "02-Genomic-discordance")

infile = here(datadir, paste(window_size, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts-tt.csv.gz", sep=""))

tree_file = here("data", "03-Selection-tests", "05_penn_7spec_iqtree.cf.rooted.tree")
# Input options
######################

if(read_data){
  cat(as.character(Sys.time()), " | Reading data: ", infile, "\n")
  all_windows = read.csv(gzfile(infile), header=T)
  all_windows_f = subset(all_windows, repeat.filter=="PASS" & missing.filter=="PASS")
  if(au_flag){
    all_windows_f = subset(all_windows_f, AU.test=="PASS")
  }
  # Read the window data
  
  cat(as.character(Sys.time()), " | Reading species tree\n")
  concat_tree = read.tree(tree_file)
  
  tree_to_df_list = treeToDF(concat_tree)
  tree_info_concat = tree_to_df_list[["info"]]
  tree_info_concat = tree_info_concat %>% separate(label, c("astral", "gcf", "scf"), sep="/", remove=F)
  tree_info_concat$bootstrap[tree_info_concat$node.type=="tip"] = NA
  tree_info_concat$bootstrap = as.numeric(tree_info_concat$astral)
  tree_info_concat$gcf = as.numeric(tree_info_concat$gcf)
  tree_info_concat$scf = as.numeric(tree_info_concat$scf)
  # Read astral tree data
  
  tree_info_concat$species = c("Rhynchomys soricoides", "Grammomys dolichurus", "Rhabdomys dilectus", "Hylomyscus alleni", "Mastomys natalensis", "Praomys delectorum", "Mus musculus", NA, NA, NA, NA, NA, NA)
  # Add full species names to tree df
}
# Read and filter the data
######################

cat(as.character(Sys.time()), " | Reading chromosomes\n")

plot_list = list()
wrf_df = data.frame("chr"=c(), "start"=c(), "wrf"=c(), "wrf.chr"=c())
diff_df = data.frame("chr"=c(), "wrf.diff"=c(), "wrf.avg"=c())
i = 1
# Some initializations

for(chrome in levels(as.factor(all_windows$chr))){
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
  
  ##########
  
  cat(as.character(Sys.time()), " | Starting chromosome ", out_chrome, "\n")
  
  chrdata = subset(all_windows, chr==chrome)
  total_windows = length(chrdata[,1])
  chrdata_f = subset(chrdata, repeat.filter=="PASS" & missing.filter=="PASS")
  if(au_flag){
    chrdata_f = subset(chrdata_f, AU.test=="PASS")
  }
  used_windows = length(chrdata_f[,1])
  # Subset and filter data for current chromosome.
  
  ##########
  
  chr_tree_file = here("data", "02-Genomic-discordance", "chrome-trees", "10kb-0.5-0.5", paste0(chrome, "-iqtree"), "concat", paste0(chrome, "_10kb_iqtree.treefile"))
  chr_tree = read.tree(file=chr_tree_file)
  # Read the tree for the current chromosome
  
  ##########
  
  cat(as.character(Sys.time()), " | -----> Calculating wRFs\n")
  chrdata$wrf = ifelse(is.na(chrdata$unparsed.tree), NA, wRF.dist(read.tree(text=chrdata$unparsed.tree), concat_tree))
  chrdata$wrf.chr = ifelse(is.na(chrdata$unparsed.tree), NA, wRF.dist(read.tree(text=chrdata$unparsed.tree), chr_tree))
  wrf_df = rbind(wrf_df, select(chrdata, chr, start, wrf, wrf.chr))
  # Calculate wrf to the species tree and the chromosome tree for every window
  
  ##########
  
  chrdata$wrf.diff = chrdata$wrf - chrdata$wrf.chr
  chrdata$wrf.avg = rowMeans(chrdata[c("wrf", "wrf.chr")], na.rm=T)
  diff_df = rbind(diff_df, select(chrdata, chr, wrf.diff, wrf.avg))
  # Make some comparisons between wrf to the species tree and the chromosome tree
  
  ##########
  
  # cat(as.character(Sys.time()), " | -----> Calculating adjacent wRFs\n")
  # chrdata$wrf.adj = NA
  # for(j in 1:nrow(chrdata)){
  #   if(is.na(chrdata[j,]$unparsed.tree)){
  #     next
  #   }
  #   
  #   cur_tree = read.tree(text=chrdata[j,]$unparsed.tree)
  #   
  #   
  #   if(j==1){
  #     prev_wrf = NA
  #   }else if(!is.na(chrdata[j-1,]$unparsed.tree)){
  #     prev_tree = read.tree(text=chrdata[j-1,]$unparsed.tree)
  #     prev_wrf = wRF.dist(cur_tree, prev_tree)
  #   }else{
  #       prev_wrf = NA
  #   }
  #   
  #   if(j==nrow(chrdata)){
  #     next_wrf = NA
  #   }else if(!is.na(chrdata[j+1,]$unparsed.tree)){
  #       next_tree = read.tree(text=chrdata[j+1,]$unparsed.tree)
  #       next_wrf = wRF.dist(cur_tree, next_tree)
  #   }else{
  #       next_wrf = NA
  #   }
  #   
  #   chrdata[j,]$wrf.adj = mean(c(prev_wrf, next_wrf), na.rm=TRUE)
  # }
  # This block calculates between each window and the windows immediately adjacent to it
  
  ##########
  
  
  p = ggplot(chrdata, aes(x=start, y=wrf)) +
    geom_point(size=1, alpha=0.1) +
    #geom_line(aes(y=rollmean(wrf, 10, fill=NA))) +
    geom_smooth(method="loess", se=F, color="#333333") +
    geom_point(aes(x=start, y=wrf.chr), alpha=0.1, color=corecol(numcol=1)) +
    geom_smooth(aes(x=start, y=wrf.chr), method="loess", se=F, color=corecol(numcol=1)) +
    #geom_point(aes(x=start, y=wrf.adj), alpha=0.1, color=corecol(numcol=1, pal="wilke", offset=1)) +
    #geom_smooth(aes(x=start, y=wrf.adj), method="loess", se=F, color=corecol(numcol=1, pal="wilke", offset=1)) +
    ggtitle(chrome) +
    bartheme()
  print(p)
  # Plot the wrfs of the current chromosome
  
  plot_list[[i]] = p
  i = i + 1
  # Save the current plot
}

######################

wrf_df$chr = factor(wrf_df$chr, levels=chr_order)

man_p = ggplot(wrf_df, aes(x=start, y=wrf)) +
  geom_point(size=1, alpha=0.1, color="#808080") +
  geom_point(aes(x=start, y=wrf.chr), alpha=0.1, color="#ffb366") +
  geom_smooth(aes(color="To species tree"), method="loess", se=F) +
  geom_smooth(aes(x=start, y=wrf.chr, color="To chromosome tree"), method="loess", se=F) +
  geom_smooth(data=subset(wrf_df, wrf > 0.15 & wrf.chr > 0.15), aes(x=start, y=wrf), method="loess", se=F, color="#333333", linetype="dashed") +
  geom_smooth(data=subset(wrf_df, wrf > 0.15 & wrf.chr > 0.15), aes(x=start, y=wrf.chr), method="loess", se=F, color=corecol(numcol=1), linetype="dashed") +
  geom_hline(yintercept=0.15, linetype="dotted", color="#999999") +
  scale_color_manual(values=c("To species tree"="#333333", "To chromosome tree"=corecol(numcol=1))) +
  xlab("") +
  ylab("wRF") +
  facet_grid(~chr, switch = "x", scales = "free_x", space = "free_x") +
  bartheme() +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom")
#print(man_p)
figfile = here("figs", "wrf-manhattan.png")
cat(as.character(Sys.time()), " | Saving figure:", figfile, "\n")

ggsave(figfile, man_p, width=24, height=8, units="in")
# Save the figure



if(save_fig){
  fig = plot_grid(plotlist=plot_list, nrow=4, ncol=5)
  # Combine the plots for each chromosome
  
  figfile = here("figs", "wrf-chromoplots.png")
  cat(as.character(Sys.time()), " | Saving figure:", figfile, "\n")
  
  ggsave(figfile, fig, width=32, height=18, units="in")
  # Save the figure
  
  ##########
  
  diff_df$chr = factor(diff_df$chr, levels=chr_order)
  
  p = ggplot(diff_df, aes(x=chr, y=wrf.diff, group=chr)) +
    #geom_boxplot() +
    geom_violin(aes(fill=chr), alpha=0.2) +
    geom_boxplot(width=0.1, outlier.shape=NA) +
    geom_hline(yintercept=0, linetype="dashed", color="#999999") +
    xlab("") +
    ylab("wRF to species tree - wRF to chromosome tree") +
    bartheme() +
    theme(axis.text.x=element_text(angle=25, hjust=1),
          legend.position="none")
  print(p)
  # Plot the distributions in difference in wrf to the species tree and the chromosome tree for each chromosome
  
  figfile = here("figs", "wrf-chrome-st-diff.png")
  cat(as.character(Sys.time()), " | Saving figure:", figfile, "\n")
  
  ggsave(figfile, p, width=10, height=6, units="in")
  # Save the figure
  
  ##########
  
  # p = ggplot(diff_df, aes(x=wrf.avg, y=wrf.diff, color=chr)) +
  #   #geom_point() +
  #   geom_smooth(method="lm", se=F) +
  #   bartheme()
  # print(p)
}


############################################################




  