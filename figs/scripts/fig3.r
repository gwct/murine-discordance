############################################################
# For rodent genomes
# Figure 3
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(ggplot2)
library(cowplot)
library(ggbeeswarm)
library(dplyr)
library(RColorBrewer)
library(here)
library(readr)
library(archive)
source(here("lib", "design.r"))

############################################################

window_size_kb = 10
# Window size in kb

window_size = window_size_kb * 1000
# Window size in bp

marker_window_size = 1
# Marker window size in Mb

au_flag = FALSE
# Set to filter out windows that don't pass the AU test

read_data = T
# Whether or not to re-read the input data

skip_one = F
# Set to only do one test chromosome

skip_x = T
# Skip the X chromosome

save_fig = T
# Whether or not to save the figure

num_reps = 100
# The number of replicates for random tree dists from dists.r

datadir = here("data", "02-Genomic-discordance")

infile = here(datadir, paste(window_size_kb, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts-tt.csv.gz", sep=""))

chrome_info_file = here(datadir, "recombination-markers", "chrome-stats.csv")

script_dir = "D:/data/rodent-genomes/dists/"
dist_archive = "D:/data/rodent-genomes/dists.tar.gz"

# Input options
######################

if(read_data){
  #window_file = paste(datadir, window_size_kb, "kb-", marker_window_size, "mb-dists.csv", sep="")
  cat(as.character(Sys.time()), " | Reading window data: ", infile, "\n")
  all_windows = read.csv(infile, header=T)
  
  cat(as.character(Sys.time()), " | Reading chrome data: ", chrome_info_file, "\n")
  chrome_info = read.csv(chrome_info_file, comment.char="#", header=T)
  
  random_file = paste0(script_dir, "random-dists")
  #random_file = "dists/random-dists"
  if(au_flag){
    random_file = paste(random_file, "-au", sep="")
  }
  random_file = paste(random_file, ".csv", sep="")
  
  cat(as.character(Sys.time()), " | Reading random tree dists: ", random_file, "\n")
  random_data = read.csv(random_file, header=T)
  #random_data = read.csv(archive_read(dist_archive, file=random_file))
}

# Read input data
######################

my_colors = colorRampPalette(brewer.pal(8, "Set2"))(20)
names(my_colors) = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19")

replicate_data = data.frame("chrome"=c(), "replicate"=c(), "nonsig.adj"=c())
summary_data = data.frame("chrome"=c(), "len"=c(), "map.len"=c(), "num.topos"=c(), 
                          "avg.nonsig.adj"=c(), "med.nonsig.adj"=c(), "avg.decay.rate"=c())
# Summary data frame

blah = data.frame("log.distance"=c(), "avg.wrf"=c())
all_chrome_p = ggplot(blah, aes(x=log.distance, y=avg.wrf))
all_chrome_log_p = ggplot(blah, aes(x=log.distance, y=avg.wrf))

i = 1
all_windows$chr = as.factor(all_windows$chr)
for(chrome in levels(all_windows$chr)){
  # Generate figures for each chromosome.
  
  if(skip_x && chrome == "chrX"){
    next
  }
  
  if(nchar(chrome) == 4 && chrome != "chrX"){
    out_chrome = gsub("chr", "chr0", chrome)
  }else{
    out_chrome = chrome
  }
  # The string for the current chromosome.
  
  cat(as.character(Sys.time()), " | Chromosome ", out_chrome, " - subsetting and filtering...\n", sep="")
  chrdata = subset(all_windows, chr==chrome)
  total_windows = length(chrdata[,1])
  chrdata_f = subset(chrdata, repeat.filter=="PASS" & missing.filter=="PASS")
  if(au_flag){
    chrdata_f = subset(chrdata_f, AU.test=="PASS")
  }
  used_windows = length(chrdata_f[,1])
  chr_len = chrdata$chr.len[1]
  # Subset data for this chromosome
  
  chr_infile = paste(script_dir, out_chrome, sep="")
  chr_statsfile = paste(script_dir, out_chrome, "-stats", sep="")
  #figfile = paste(script_dir, out_chrome, sep="")
  if(au_flag){
    chr_infile = paste(chr_infile, "-au-filter", sep="")
    chr_statsfile = paste(chr_statsfile, "-au-filter", sep="")
    #figfile = paste(figfile, "-au-filter", sep="")
  }
  chr_infile = paste(chr_infile, ".csv", sep="")
  chr_statsfile = paste(chr_statsfile, ".csv", sep="")
  #figfile = paste(figfile, ".png", sep="")

  cat(as.character(Sys.time()), " | Chromosome ", out_chrome, " - reading dist data:       ", chr_infile, "\n", sep="")
  dist_data = read.csv(chr_infile, header=T)
  dist_data = subset(dist_data, measure=="wrf")

  cat(as.character(Sys.time()), " | Chromosome ", out_chrome, " - reading stats data:      ", chr_statsfile, "\n", sep="")
  stats_data = read.csv(chr_statsfile, header=T)
  # Set up input file names (output files from gen_data) and read data.
  
  # chr_infile = paste0("dists/",out_chrome)
  # chr_statsfile = paste0("dists/", out_chrome, "-stats")
  # if(au_flag){
  #   chr_infile = paste0(chr_infile, "-au-filter")
  #   chr_statsfile = paste0(chr_statsfile, "-au-filter")
  # }
  # chr_infile = paste0(chr_infile, ".csv")
  # chr_statsfile = paste0(chr_statsfile, ".csv")
  # 
  # cat(as.character(Sys.time()), " | Chromosome ", out_chrome, " - reading dist data:       ", chr_infile, "\n", sep="")
  # dist_data = read.csv(archive_read(dist_archive, file=chr_infile), header=T)
  # dist_data = subset(dist_data, measure=="wrf")
  # 
  # cat(as.character(Sys.time()), " | Chromosome ", out_chrome, " - reading stats data:      ", chr_statsfile, "\n", sep="")
  # stats_data = read.csv(archive_read(dist_archive, file=chr_statsfile), header=T)
  # # Set up input file names (output files from gen_data) and read data.  
  
  cat(as.character(Sys.time()), " | Chromosome ", out_chrome, " - calculating decay rates (slopes)\n", sep="")
  cur_nonsig_adjs = c()
  for(k in 1:num_reps){
    cur_dists = dist_data[dist_data$replicate==k,]
    cur_nonsig_adjs = c(cur_nonsig_adjs, cur_dists$nonsig.adj)
    cur_dists$nonsig.adj = cur_dists$nonsig.adj * window_size / 1000000
    replicate_data = rbind(replicate_data, data.frame("chrome"=chrome, "replicate"=k, 
                                                      "nonsig.adj"=cur_dists$nonsig.adj))
  }
  
  avg_nonsig_adj = mean(cur_nonsig_adjs)
  med_nonsig_adj = median(cur_nonsig_adjs)
  cur_stats = subset(stats_data, adj <= avg_nonsig_adj)
  
  cur_stats$adj = cur_stats$adj * window_size / 1000000
  cur_fit = lm(cur_stats$mean.wrf ~ cur_stats$adj)
  cur_ln_fit = lm(cur_stats$mean.wrf ~ log10(cur_stats$adj))
  #lines(cur_data$adj,predict(cur_ln_fit),col='red')
  cur_slope = cur_fit$coefficients[2]
  
  cur_color = my_colors[chrome]
  
  all_chrome_p = all_chrome_p +
    geom_smooth(data=cur_stats, aes(x=adj, y=median.wrf), method="lm", formula=y~log10(x), se=F, color=cur_color)
  
  all_chrome_log_p = all_chrome_log_p +
    geom_smooth(data=cur_stats, aes(x=log10(adj), y=median.wrf), method="lm", formula=y~x, se=F, color=cur_color)
  
  summary_data = rbind(summary_data, data.frame("chrome"=chrome, "len"=chr_len, 
                                                "map.len"=chrome_info$max.map.site[chrome_info$chr==chrome], 
                                                "num.topos"=max(chrdata_f$topo.num.chrom), 
                                                "avg.nonsig.adj"=avg_nonsig_adj, "med.nonsig.adj"=med_nonsig_adj,
                                                "avg.decay.rate"=cur_slope))
  
  # cur_slopes$dummy = paste(chrome, sep="")
  # dist_data$dummy = paste(chrome, sep="")
  # # Dummy data for boxplot axis
  # 
  # dist_data$nonsig.adj = dist_data$nonsig.adj * window_size / 1000000
  #out_data$rand.adj = out_data$rand.adj / 10
  # Convert to Mb
  
  # stats_data$adj = stats_data$adj * window_size / 1000000
  # stats_data = subset(stats_data, adj <= dist_data$nonsig.adj)
  # 
  # all_chrome_p = all_chrome_p +
  #   geom_smooth(data=stats_data, aes(x=adj, y=median.wrf), method="lm", formula=y~log10(x), se=F, color=my_colors[i])
  # 
  # all_chrome_log_p = all_chrome_log_p +
  #   geom_smooth(data=stats_data, aes(x=log10(adj), y=median.wrf), method="lm", formula=y~x, se=F, color=my_colors[i])
  # 
  # summary_data = rbind(summary_data, data.frame("chrome"=chrome, "len"=chr_len, "map.len"=chrome_info$max.map.site[chrome_info$chr==chrome], "num.topos"=max(chrdata_f$topo.num.chrom), "avg.nonsig.adj"=mean(dist_data$nonsig.adj), "avg.decay.rate"=mean(cur_slopes$slopes)))
  
  i = i + 1
  if(skip_one){
    stop("skip one ok")
  }
}
## Chrome loop

cat(as.character(Sys.time()), " | Fig3: Generating figure panels A and B\n")
all_chrome_p = all_chrome_p + 
  geom_hline(yintercept=median(random_data$wrf), size=2, linetype="dashed", color="#333333") +
  xlab("Distance between\nwindows (Mb)") +
  ylab("Average wRF") +
  bartheme()
print(all_chrome_p)

all_chrome_log_p = all_chrome_log_p + 
  geom_hline(yintercept=median(random_data$wrf), size=2, linetype="dashed", color="#333333") +
  xlab("Distance between\nwindows (log Mb)") +
  ylab("Average wRF") +
  bartheme()
print(all_chrome_log_p)

######################

cat(as.character(Sys.time()), " | Fig3: Generating figure panel C\n")
#avg_nonsig_chr = signif(mean(summary_data$avg.nonsig.adj), 2)
#med_nonsig_chr = signif(median(summary_data$avg.nonsig.adj), 2)
# Summary stats of chromosome summaries

blah = data.frame()
# Dummy data frame for layering

replicate_data$chrome = factor(replicate_data$chrome, levels=rev(c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19")))

avg_p = ggplot(replicate_data, aes(x=chrome, y=nonsig.adj, group=chrome, color=chrome)) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(alpha=0.2, size=2)

# k = 1
# for(chrome in levels(summary_data$chr)){
#   print(chrome)
#   if(skip_x && chrome == "chrX"){
#     next
#   }
#   avg_p = avg_p + geom_segment(x=k, y=0, xend=k, yend=min(replicate_data$nonsig.adj[replicate_data$chrome==chrome]), color="#666666", linetype="dotted")
#   k = k + 1
# }
# Add the dotted lines for each chrome

avg_p = avg_p + xlab("") +
  ylab("Distance to random-like\ndistribution of wRF (Mb)\n(100 replicates)") +
  scale_y_continuous(expand=c(0,0), limits=c(0,90)) +
  scale_color_manual(values=my_colors) +
  bartheme() +
  theme(legend.position="none") +
  coord_flip()
print(avg_p)
# Distance to non-random trees
######

######################

cat(as.character(Sys.time()), " | Fig3: Generating figure panel D\n")

avg_decay_rate = signif(mean(summary_data$avg.decay.rate), 2)
med_decay_rate = signif(median(summary_data$avg.decay.rate), 2)

summary_data$chrome = factor(summary_data$chrome, levels=rev(c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19")))

dec_p = ggplot(summary_data, aes(x=chrome, y=avg.decay.rate)) +
  geom_hline(yintercept=mean(summary_data$avg.decay.rate), size=2, linetype="dashed", color="#d3d3d3") +
  geom_hline(yintercept=median(summary_data$avg.decay.rate), size=1.5, linetype="dashed", color="#999999")

k = 1
for(chrome in levels(summary_data$chrome)){
  print(chrome)
  dec_p = dec_p + geom_segment(x=k, y=-0.0001, xend=k, yend=summary_data$avg.decay.rate[summary_data$chrome==chrome], color="#666666", linetype="dotted")
  k = k + 1
}

dec_p = dec_p + geom_point(size=4, color="#920000") +
  #geom_text(data=blah, aes(x=18.4, y=avg_decay_rate+0.00004, label=paste("Avg. = ", avg_decay_rate, sep="")),  size=4, color="#333333") +
  xlab("Chromosome") +
  ylab("Average\ndecay rate") +
  scale_y_continuous(expand=c(0,0), limits=c(-0.0001,0.012)) +
  bartheme() +
  coord_flip()
print(dec_p)
# Avg. slope

######################

cat(as.character(Sys.time()), " | Fig3: Combining panels\n")
fig_3ab = plot_grid(all_chrome_p, 
                    all_chrome_log_p,
                    nrow=2, labels=c("A", "B"), label_size=16)

fig3 = plot_grid(fig_3ab, avg_p, dec_p, ncol=3, labels=c("", "C", "D"), label_size=16)

######################

if(save_fig){
  figfile = here("figs", "fig3.png")
  cat(as.character(Sys.time()), " | Fig3: Saving figure:", figfile, "\n")
  ggsave(filename=figfile, fig3, width=12, height=7.5, units="in")
}
