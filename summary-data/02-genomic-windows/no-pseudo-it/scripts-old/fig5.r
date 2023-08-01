############################################################
# For penn genomes, 06.21
# Figure 5
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(ggbeeswarm)

############################################################

window_size_kb = 10
# Window size in kb

window_size = window_size_kb * 1000
# Window size in bp

marker_window_size = 5
# Marker window size in Mb

au_flag = FALSE
# Set to filter out windows that don't pass the AU test

read_data = T

gen_figs = T
# Set to generate figures.

skip_one = F
# Set to only do one test chromosome

skip_x = T
# Skip the X chromosome.

featured_chr = "chr7"
# The example chromosome for panels A and B

max_dist_mb = 5
# Distance limit for plots

save_fig = T
# Whether or not to save the figure

gen_supp = F

point_alpha = 0.1

datadir = "C:/Users/Gregg/Box Sync/rodents/penn/paper/data/"
#distdir = "C:/Users/Gregg/Desktop/dists/"
distdir = "D:/data/rodent-genomes/dists/"
tree_file = paste(datadir, window_size_kb, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts-tt.csv", sep="")
marker_file = paste(datadir, window_size_kb, "kb-", marker_window_size, "mb-dists.csv", sep="")

chrome_info_file = paste(datadir, "recombination-markers/chrome-stats.csv", sep="")
# Input options
######################

if(read_data){
  cat(as.character(Sys.time()), " | Reading tree window data: ", tree_file, "\n")
  all_windows = read.csv(tree_file, header=T)
  all_windows_f = subset(all_windows, repeat.filter=="PASS" & missing.filter=="PASS")
  if(au_flag){
    all_windows_f = subset(all_windows_f, AU.test=="PASS")
  }
  
  cat(as.character(Sys.time()), " | Reading marker window data: ", marker_file, "\n")
  marker_windows = read.csv(marker_file, header=T)
  
  #cat(as.character(Sys.time()), " | Reading chrome info: ", chrome_info_file, "\n")
}
# Read and filter the window data
######################

slope_data = data.frame("chr"=c(), "slope.wrf"=c(), "type"=c())
int_data = data.frame("chr"=c(), "int.wrf"=c(), "type"=c())
diff_data = data.frame("chr"=c(), "adj"=c(), "wrf.diff"=c(), "type"=c())
# Initialize data frames
######################

for(chrome in levels(as.factor(all_windows$chr))){
  # Generate figures for each chromosome.
  
  if(skip_one && chrome != "chr7"){
    next
  }
  
  if(skip_x && chrome == "chrX"){
    next
  }
  
  if(nchar(chrome) == 4 && chrome != "chrX"){
    out_chrome = gsub("chr", "chr0", chrome)
  }else{
    out_chrome = chrome
  }
  # The string for the current chromosome.
  
  # Parse chromosome strings and other conditions
  ######################
  
  cat(as.character(Sys.time()), " Starting chromosome | ", out_chrome, "\n", sep="")
  
  chrdata = subset(all_windows, chr==chrome)
  total_windows = length(chrdata[,1])
  chrdata_f = subset(chrdata, repeat.filter=="PASS" & missing.filter=="PASS")
  if(au_flag){
    chrdata_f = subset(chrdata_f, AU.test=="PASS")
  }
  used_windows = length(chrdata_f[,1])
  chr_len = chrdata$chr.len[1]
  # Subset data for this chromosome
  ######################
  
  full_stats_file = paste(distdir, out_chrome, "-stats.csv", sep="")
  gene_file = paste(datadir, "coord-query/mm10-genes-", marker_window_size ,"-Mb.bed.", chrome, ".dists", sep="")
  uce_file = paste(datadir, "coord-query/mm10-uces-", marker_window_size ,"-Mb.bed.", chrome, ".dists", sep="")
  hs_file = paste(datadir, "coord-query/mm10-hotspots-", marker_window_size ,"-Mb.bed.", chrome, ".dists", sep="")
  ps_file = paste(datadir, "coord-query/mm10-ps-", marker_window_size ,"-Mb.bed.", chrome, ".dists", sep="")
  
  # File names for this chromosome
  ######################
  
  cat(as.character(Sys.time()), " | Chromosome ", out_chrome, " - reading data: ", full_stats_file, "\n", sep="")
  full_stats = read.csv(full_stats_file, header=T)
  full_stats$adj = full_stats$adj * window_size / 1000000
  full_stats = subset(full_stats, adj <= max_dist_mb)
  full_stats$log.adj = log10(full_stats$adj)
  # Chromosome-wide data
  ######################
  
  cat(as.character(Sys.time()), " | Chromosome ", out_chrome, " - reading data: ", gene_file, "\n", sep="")
  gene_dists = read.csv(gene_file, header=T)
  gene_dists$adj = gene_dists$adj * window_size / 1000000
  gene_dists$adj[gene_dists$adj < 0] = gene_dists$adj[gene_dists$adj < 0] * -1
  gene_dists = subset(gene_dists, adj <= max_dist_mb)
  gene_dists_avg = gene_dists %>% group_by(adj) %>% summarize(avg.wrf=mean(wrf, na.rm=T))
  gene_dists_avg$log.adj = log10(gene_dists_avg$adj)
  gene_diffs = data.frame("adj"=gene_dists_avg$adj, "wrf.diff" = gene_dists_avg$avg.wrf - full_stats$mean.wrf)
  gene_diffs$chr = out_chrome
  gene_diffs$type = "All genes"
  #gene_diffs$type2 = "Genes"
  # Gene data
  ######################
  
  cat(as.character(Sys.time()), " | Chromosome ", out_chrome, " - reading data: ", ps_file, "\n", sep="")
  ps_dists = read.csv(ps_file, header=T)
  ps_dists$adj = ps_dists$adj * window_size / 1000000
  ps_dists$adj[ps_dists$adj < 0] = ps_dists$adj[ps_dists$adj < 0] * -1
  ps_dists = subset(ps_dists, adj <= max_dist_mb)
  ps_dists_avg = ps_dists %>% group_by(adj) %>% summarize(avg.wrf=mean(wrf, na.rm=T))
  ps_dists_avg$log.adj = log10(ps_dists_avg$adj)
  ps_diffs = data.frame("adj"=ps_dists_avg$adj, "wrf.diff" = ps_dists_avg$avg.wrf - full_stats$mean.wrf)
  ps_diffs$chr = out_chrome
  ps_diffs$type = "Positively selected genes"
  #ps_diffs$type2 = "Positively selected genes"
  # PS gene data
  ######################
  
  cat(as.character(Sys.time()), " | Chromosome ", out_chrome, " - reading data: ", uce_file, "\n", sep="")
  uce_dists = read.csv(uce_file, header=T)
  uce_dists$adj = uce_dists$adj * window_size / 1000000
  uce_dists$adj[uce_dists$adj < 0] = uce_dists$adj[uce_dists$adj < 0] * -1
  uce_dists = subset(uce_dists, adj <= max_dist_mb)
  uce_dists_avg = uce_dists %>% group_by(adj) %>% summarize(avg.wrf=mean(wrf, na.rm=T))
  uce_dists_avg$log.adj = log10(uce_dists_avg$adj)
  uce_diffs = data.frame("adj"=uce_dists_avg$adj, "wrf.diff" = uce_dists_avg$avg.wrf - full_stats$mean.wrf)
  uce_diffs$chr = out_chrome
  uce_diffs$type = "UCEs"
  # UCE data
  ######################
  
  cat(as.character(Sys.time()), " | Chromosome ", out_chrome, " - reading data: ", hs_file, "\n", sep="")
  hs_dists = read.csv(hs_file, header=T)
  hs_dists$adj = hs_dists$adj * window_size / 1000000
  hs_dists$adj[hs_dists$adj < 0] = hs_dists$adj[hs_dists$adj < 0] * -1
  hs_dists = subset(hs_dists, adj <= max_dist_mb)
  hs_dists_avg = hs_dists %>% group_by(adj) %>% summarize(avg.wrf=mean(wrf, na.rm=T))
  hs_dists_avg$log.adj = log10(hs_dists_avg$adj)
  hs_diffs = data.frame("adj"=hs_dists_avg$adj, "wrf.diff" = hs_dists_avg$avg.wrf - full_stats$mean.wrf)
  hs_diffs$chr = out_chrome
  hs_diffs$type = "Hotspots"
  # Set up input file names (output files from gen_data) and read data.
  # Hotspot data
  ######################
  
  fill_fit = lm(full_stats$mean.wrf ~ log10(full_stats$adj))
  full_slope = fill_fit$coefficients[2]
  #full_int = fill_fit$coefficients[1]
  full_int = full_stats$mean.wrf[full_stats$adj==0.01]
  # Chromosome-wide fits
  ######################
  
  gene_fit = lm(gene_dists_avg$avg.wrf ~ log10(gene_dists_avg$adj))
  gene_slope = gene_fit$coefficients[2]
  gene_int = gene_dists_avg$avg.wrf[gene_dists_avg$adj==0.01]
  #gene_int = gene_fit$coefficients[1]
  # Gene fits
  ######################
  
  ps_fit = lm(ps_dists_avg$avg.wrf ~ log10(ps_dists_avg$adj))
  ps_slope = ps_fit$coefficients[2]
  ps_int = ps_dists_avg$avg.wrf[ps_dists_avg$adj==0.01]
  #gene_int = gene_fit$coefficients[1]
  # Positively selected gene fits
  ######################
  
  uce_fit = lm(uce_dists_avg$avg.wrf ~ log10(uce_dists_avg$adj))
  uce_slope = uce_fit$coefficients[2]
  uce_int = uce_dists_avg$avg.wrf[uce_dists_avg$adj==0.01]
  #uce_int = uce_fit$coefficients[1]
  # UCE fits
  ######################
  
  hs_fit = lm(hs_dists_avg$avg.wrf ~ log10(hs_dists_avg$adj))
  hs_slope = hs_fit$coefficients[2]
  hs_int = hs_dists_avg$avg.wrf[hs_dists_avg$adj==0.01]
  #uce_int = uce_fit$coefficients[1]
  # Calculate fits for both models
  # Hotspot fits
  ######################
  
  full_stats$type = "Chromosomes"
  full_stats$type = factor(full_stats$type, levels=c("Chromosomes", "Hotspots", "All genes", "Positively selected genes", "UCEs"))
  # Chromosome-wide labels
  ######################
  
  gene_dists_avg$type = "All genes"
  gene_dists_avg$type = factor(gene_dists_avg$type, levels=c("Chromosomes", "Hotspots", "All genes", "Positively selected genes", "UCEs"))
  #gene_dists_avg$type2 = "All genes"
  #gene_dists_avg$type2 = factor(gene_dists_avg$type2, levels=c("All genes", "Positively selected genes"))
  # Gene labels
  ######################
  
  ps_dists_avg$type = "Positively selected genes"
  ps_dists_avg$type = factor(ps_dists_avg$type, levels=c("Chromosomes", "Hotspots", "All genes", "Positively selected genes", "UCEs"))
  #ps_dists_avg$type2 = "Positively selected genes"
  #ps_dists_avg$type2 = factor(ps_dists_avg$type2, levels=c("All genes", "Positively selected genes"))
  # PS gene labels
  ######################
     
  uce_dists_avg$type = "UCEs"
  uce_dists_avg$type = factor(uce_dists_avg$type, levels=c("Chromosomes", "Hotspots", "All genes", "Positively selected genes", "UCEs"))
  # UCE labels
  ######################
  
  hs_dists_avg$type = "Hotspots"
  hs_dists_avg$type = factor(hs_dists_avg$type, levels=c("Chromosomes", "Hotspots", "All genes", "Positively selected genes", "UCEs"))
  # Hotspot labels
  ######################

  cat(as.character(Sys.time()), " | Fig5: Rendering panel A for ", chrome, "\n")
  log_p = ggplot(full_stats, aes(x=adj, y=mean.wrf, color=type)) +
    geom_point(size=2, alpha=point_alpha) +
    geom_smooth(method="lm", formula=y~log10(x), se=F, linetype="dashed") +
    #geom_line() +
    geom_point(data=gene_dists_avg, aes(x=adj, y=avg.wrf), size=2, alpha=point_alpha) +
    geom_smooth(data=gene_dists_avg, aes(x=adj, y=avg.wrf), method="lm", formula=y~log10(x), se=F, linetype="dashed") +
    geom_point(data=ps_dists_avg, aes(x=adj, y=avg.wrf), size=2, alpha=point_alpha) +
    geom_smooth(data=ps_dists_avg, aes(x=adj, y=avg.wrf), method="lm", formula=y~log10(x), se=F, linetype="dashed") +
    geom_point(data=uce_dists_avg, aes(x=adj, y=avg.wrf), size=2, alpha=point_alpha) +
    geom_smooth(data=uce_dists_avg, aes(x=adj, y=avg.wrf), method="lm", formula=y~log10(x), se=F, linetype="dashed") +
    #geom_line(data=uce_dists_avg, aes(x=adj, y=avg.wrf)) +
    geom_point(data=hs_dists_avg, aes(x=adj, y=avg.wrf), size=2, alpha=point_alpha) +
    geom_smooth(data=hs_dists_avg, aes(x=adj, y=avg.wrf), method="lm", formula=y~log10(x), se=F, linetype="dashed") +
    xlab("Distance between\nwindows (Mb)") +
    ylab("Mean wRF") +
    scale_color_manual(name="", values=c('#006ddb', "#333333", '#db6d00', '#004949', '#920000'), labels=c("Chromosomes", "Hotspots", "All genes", "Positively selected genes", "UCEs"), drop=FALSE) +
    bartheme() + 
    theme(legend.position="right",
          axis.title.x = element_text(margin = unit(c(0, 3, 3, 3), "mm")),
          plot.margin = unit(c(0,0,0.25,1), "cm"))
  print(log_p)
  # Log correlation plot
  ######################
  
  cat(as.character(Sys.time()), " | Fig5: Rendering panel B for ", chrome, "\n")
  lin_p = ggplot(full_stats, aes(x=log10(adj), y=mean.wrf, color=type)) +
    geom_point(size=2, alpha=point_alpha) +
    geom_smooth(method="lm", formula=y~x, se=F, linetype="dashed") +
    #geom_line() +
    geom_point(data=gene_dists_avg, aes(x=log10(adj), y=avg.wrf), size=2, alpha=point_alpha) +
    geom_smooth(data=gene_dists_avg, aes(x=log10(adj), y=avg.wrf), method="lm", formula=y~x, se=F, linetype="dashed") +
    geom_point(data=ps_dists_avg, aes(x=log10(adj), y=avg.wrf), size=2, alpha=point_alpha) +
    geom_smooth(data=ps_dists_avg, aes(x=log10(adj), y=avg.wrf), method="lm", formula=y~x, se=F, linetype="dashed") +
    geom_point(data=uce_dists_avg, aes(x=log10(adj), y=avg.wrf), size=2, alpha=point_alpha) +
    geom_smooth(data=uce_dists_avg, aes(x=log10(adj), y=avg.wrf), method="lm", formula=y~x, se=F, linetype="dashed") +
    #geom_line(data=uce_dists_avg, aes(x=log10(adj), y=avg.wrf)) +
    geom_point(data=hs_dists_avg, aes(x=log10(adj), y=avg.wrf), size=2, alpha=point_alpha) +
    geom_smooth(data=hs_dists_avg, aes(x=log10(adj), y=avg.wrf), method="lm", formula=y~x, se=F, linetype="dashed") +
    xlab("Log distance between\nwindows (Mb)") +
    ylab("Mean wRF") +
    scale_color_manual(name="", values=c('#006ddb', "#333333", '#db6d00', '#004949', '#920000'), labels=c("Chromosomes", "Hotspots", "All genes", "Positively selected genes", "UCEs"), drop=FALSE) +
    bartheme() +
    theme(legend.position="right",
          axis.title.x = element_text(margin = unit(c(0, 3, 3, 3), "mm")),
          plot.margin = unit(c(0,0,0.25,1), "cm"))
  print(lin_p)
  # Linear correlation plot
  ######################
  
  if(chrome == featured_chr){
    fig_5a = log_p
    fig_5b = lin_p
  }
  # Save the plots if this is the chromosome in the main fig
  ######################
  
  if(gen_supp){
    log_p = log_p + theme(legend.position="bottom",
                          axis.text.x=element_text(angle=0, hjust=0.5),
                          axis.text=element_text(size=8),
                          axis.title=element_text(size=10),
                          plot.margin = unit(c(1,1,0,1), "cm"))
    lin_p = lin_p + theme(legend.position="bottom",
                          axis.text.x=element_text(angle=0, hjust=0.5),
                          axis.text=element_text(size=8),
                          axis.title=element_text(size=10),
                          axis.title.y=element_blank(),
                          plot.margin = unit(c(1,1,0,1), "cm"),
                          legend.text=element_text(size=8))
    
    
    supp_file = paste("../figs/supp/supp-feature-plots/", chrome, ".png", sep="")
    cat(as.character(Sys.time()), " | Combining plots for supp fig: ", chrome, "\n")
    supp_p_main = plot_grid(log_p + theme(legend.position="none"),
                       lin_p + theme(legend.position="none"),
                       ncol=2, labels=c("A", "B"), label_size=14)
    
    ptitle = ggdraw() + draw_label(chrome, size=16, fontface="bold", x=0, hjust=0) + theme(plot.margin=margin(0,0,7,7))
    supp_p_title = plot_grid(ptitle, supp_p_main, ncol=1, rel_heights=c(0.1,1))
    
    supp_legend = get_legend(lin_p)
    supp_p = plot_grid(supp_p_title, supp_legend, nrow=2, rel_heights=c(1,0.1))
    
    
    cat(as.character(Sys.time()), " | Saving supp figure:", supp_file, "\n")
    ggsave(filename=supp_file, supp_p, width=6, height=3, units="in")
  }

  # Save the supplemental figure for each chromosome
  ######################
  
  # if(chrome == featured_chr){
  #   cat(as.character(Sys.time()), " | Fig5: Rendering panel A2 for ", chrome, "\n")
  #   fig_5_2a = ggplot(gene_dists_avg, aes(x=adj, y=avg.wrf, color=type2)) +
  #     geom_point(size=2, alpha=0.2) +
  #     geom_smooth(method="lm", formula=y~log10(x), se=F, linetype="dashed") +
  #     #geom_line() +
  #     geom_point(data=ps_dists_avg, aes(x=adj, y=avg.wrf), size=2, alpha=0.2) +
  #     geom_smooth(data=ps_dists_avg, aes(x=adj, y=avg.wrf), method="lm", formula=y~log10(x), se=F, linetype="dashed") +
  #     xlab("Distance between\nwindows (Mb)") +
  #     ylab("Mean wRF") +
  #     scale_color_manual(name="", values=c('#db6d00', corecol(pal="wilke", numcol=1, offset=1)), labels=c("All genes","Positively selected genes"), drop=FALSE) +
  #     bartheme() +
  #     theme(legend.position="right",
  #           axis.text.x=element_text(angle=45, hjust=1),
  #           axis.title.x = element_text(margin = unit(c(0, 3, 3, 3), "mm")),
  #           plot.margin = unit(c(0,0,0.25,1), "cm"))
  #   print(fig_5_2a)
  # 
  #   cat(as.character(Sys.time()), " | Fig5: Rendering panel B2 for  ", chrome, "\n")
  #   fig_5_2b = ggplot(gene_dists_avg, aes(x=log10(adj), y=avg.wrf, color=type2)) +
  #     #geom_point(size=2, alpha=0.2) +
  #     geom_smooth(method="lm", formula=y~x, se=F, linetype="dashed") +
  #     #geom_line() +
  #     geom_point(data=ps_dists_avg, aes(x=log10(adj), y=avg.wrf), size=2, alpha=0.2) +
  #     geom_smooth(data=ps_dists_avg, aes(x=log10(adj), y=avg.wrf), method="lm", formula=y~x, se=F, linetype="dashed") +
  #     xlab("Log distance between\nwindows (Mb)") +
  #     ylab("Mean wRF") +
  #     scale_color_manual(name="", values=c('#db6d00', corecol(pal="wilke", numcol=1, offset=7)), labels=c("All genes","Positively selected genes"), drop=FALSE) +
  #     bartheme() +
  #     theme(legend.position="right",
  #           axis.text.x=element_text(angle=45, hjust=1),
  #           axis.title.x = element_text(margin = unit(c(0, 3, 3, 3), "mm")),
  #           plot.margin = unit(c(0,0,0.25,1), "cm"))
  #   print(fig_5_2b)
  #   stop()
  # }
  # PS vs Gene plots
  ######################
  
  # diff_data_chr = rbind(gene_diffs, uce_diffs, hs_diffs)
  # test_chr = ggplot(uce_diffs, aes(x=adj, y=wrf.diff, group=type, color=type)) +
  #   geom_point(alpha=0.3) +
  #   geom_smooth(method="loess", linetype="dashed", se=F) +
  #   geom_line() +
  #   bartheme()
  # print(test_chr)
  # Test plots
  ######################
  
  if(skip_one){
    stop("skip one ok")
  }
  
  slope_data = rbind(slope_data, data.frame("chr"=c(chrome,chrome,chrome,chrome,chrome), "slope.wrf"=c(full_slope,gene_slope,ps_slope,uce_slope,hs_slope), "type"=c("Chromosomes", "All genes", "Positively selected genes", "UCEs", "Hotspots")))
  int_data = rbind(int_data, data.frame("chr"=c(chrome,chrome,chrome,chrome,chrome), "int.wrf"=c(full_int,gene_int,ps_int,uce_int,hs_int), "type"=c("Chromosomes", "All genes", "Positively selected genes", "UCEs", "Hotspots")))
  diff_data = rbind(diff_data, gene_diffs, ps_diffs, uce_diffs, hs_diffs)
  # Add chromosome data to full data frames
}
## Chrome loop
######################

cat(as.character(Sys.time()), " | Fig5: Rendering panel C\n")
slope_data$type = factor(slope_data$type, levels=c("Chromosomes", "Hotspots", "All genes", "Positively selected genes", "UCEs"))
fig_5c = ggplot(slope_data, aes(x=type, y=slope.wrf, group=type)) +
  geom_point(size=3, alpha=0.6, color=corecol(numcol=1, pal="wilke")) +
  geom_boxplot(outlier.shape=NA, alpha=0.3, width=0.5, fill="transparent", color="#666666") +
  xlab("") +
  ylab("Slope out to 1Mb") +
  bartheme() +
  theme(legend.position="none") +
  bartheme() +
  theme(axis.text.x=element_text(angle=25, hjust=1),
        plot.margin = unit(c(0,0,0,1), "cm"))
print(fig_5c)

######################

# chrome_int = subset(int_data, type=="Chrome")
# genes_int = subset(int_data, type=="All genes")
# ps_int = subset(int_data, type=="Positively selected genes")
# uces_int = subset(int_data, type=="UCEs")
# 
# first_seg = data.frame("x1"=chrome_int$type, "y1"=chrome_int$int.wrf, "x2"=genes_int$type, "y2"=genes_int$int.wrf, "type"=NA)
# second_seg = data.frame("x1"=genes_int$type, "y1"=genes_int$int.wrf, "x2"=uces_int$type, "y2"=uces_int$int.wrf, "type"=NA)

######################

cat(as.character(Sys.time()), " | Fig5: Rendering panel D\n")
int_data$type = factor(int_data$type, levels=c("Chromosomes",  "Hotspots", "All genes", "Positively selected genes", "UCEs"))
fig_5d = ggplot(int_data, aes(x=type, y=int.wrf, group=type)) +
  geom_point(size=3, alpha=0.6, color=corecol(numcol=1, pal="wilke", offset=1)) +
  geom_boxplot(outlier.shape=NA, alpha=0.3, width=0.5, fill="transparent", color="#666666") +
  xlab("") +
  ylab("wRF of tres adjacent\nto feature") +
  bartheme() +
  theme(legend.position="none") +
  bartheme() +
  theme(axis.text.x=element_text(angle=25, hjust=1),
        plot.margin = unit(c(0,0,0,1), "cm"))
print(fig_5d)

######################

cat(as.character(Sys.time()), " | Fig5: Rendering panel E\n")
x = subset(diff_data, type=="Genes")
fig_5e = ggplot(diff_data, aes(x=adj, y=wrf.diff, group=type, color=type)) +
  geom_point(alpha=0.1) +
  geom_smooth(method="lm", se=F, linetype="dashed") +
  xlab("Distance from\nfeature (Mb)") +
  ylab("Difference from\nchromosome-wide wRF") +
  scale_color_manual(name="", values=c("#333333", '#db6d00','#004949','#920000'), breaks=c( "Hotspots", "All genes", "Positively selected genes", "UCEs")) +
  bartheme() +
  theme(plot.margin = unit(c(0,0,0,1), "cm"))
print(fig_5e)

fig_5_legend = get_legend(fig_5a)

fig_5ab = plot_grid(fig_5a + theme(legend.position="none"), 
                    fig_5b + theme(legend.position="none"),
                    ncol=2, labels=c("A", "B"), label_size=16, align='vh')
fig_5cd = plot_grid(fig_5c, fig_5d, labels=c("C", "D"), label_size=16, align='vh')
fig_5eleg = plot_grid(fig_5e + theme(legend.position="none"), fig_5_legend, 
                      labels=c("E", ""), label_size=16, rel_widths=c(0.7,0.3), align='vh')
#fig_5cde = plot_grid(fig_5c, fig_5d, fig_5e + theme(legend.position="none"),
#                     ncol=3, labels=c("C", "D", "E"), label_size=16, align='vh', rel_widths=c(0.2,0.2,0.6))
fig_5abcde = plot_grid(fig_5ab, fig_5cd, fig_5eleg, nrow=3, align='h')

if(save_fig){
  fig5file = "../figs/fig5.png"
  cat(as.character(Sys.time()), " | Fig5: Saving figure:", fig5file, "\n")
  ggsave(filename=fig5file, fig_5abcde, width=7.5, height=8, units="in")
}

############################################################







