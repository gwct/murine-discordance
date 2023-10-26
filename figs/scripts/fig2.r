############################################################
# For rodent
# Figure 2
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(ggplot2)
library(cowplot)
library(ggsignif)
library(here)
library("ggtree")
source(here("figs", "scripts", "lib", "design.r"))

############################################################
cat("----------\n")

window_size = 10
# Window size in kb

marker_window_size = 5
# Marker window size in Mb

au_flag = F
# Set to filter out windows that don't pass the AU test

read_data = T
# Whether to read the initial data or not

gen_chromo = T
# Whether or not to re-genereate the chromoplot

save_fig = T
# Whether or not to save the figure

skip_one = F
# Set to only do one test chromosome

max_tree_rank = 3
# The number of top trees to print/color

datadir = here("summary-data", "02-genomic-windows")

infile = here(datadir, paste(window_size, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts.csv", sep=""))

chrome_info_file = here(datadir, "recombination-markers", "chrome-stats.csv")
# Input options
######################

if(read_data){
  cat(as.character(Sys.time()), " | Reading data: ", infile, "\n")
  all_windows = read.csv(infile, header=T, comment.char="#")
  all_windows_f = subset(all_windows, repeat.filter=="PASS" & missing.filter=="PASS")
  if(au_flag){
    all_windows_f = subset(all_windows_f, AU.test=="PASS")
  }
  
  cat(as.character(Sys.time()), " | Reading chrome info: ", chrome_info_file, "\n")
  chrome_info = read.csv(chrome_info_file, header=T, comment.char="#")
}
# Read and filter the data
######################

topo_counts = unique(subset(all_windows_f, select=c(topo.num.overall, topo.count.overall, topo.rank.overall, topo.color)))
topo_counts = topo_counts[order(topo_counts$topo.rank.overall),]
# Get the topology counts

col_rank = 1
for(i in 1:nrow(topo_counts)){
  row = topo_counts[i,]
  cur_col = row$topo.color
  if(cur_col != "#999999"){
    new_label = paste("Top ", col_rank, sep="")
    col_rank = col_rank + 1
  }else{
    new_label = "Other"
  }
  topo_counts$col.cat[topo_counts$topo.num.overall==i] = new_label
  
  if(row$topo.rank.overall > max_tree_rank){
    topo_counts[i,]$topo.color = "#999999"
  }
}
# Add label and color columns to the topo counts

topo_counts$outline = "transparent"
#topo_counts$outline[topo_counts$is.astral] = "#24ff24"
#topo_counts$outline[topo_counts$is.concat] = "#000000"
# Add an outline column -- black outline is chrome topo

topo_counts$topo.prop = round(topo_counts$topo.count / sum(topo_counts$topo.count), 2)
# Percentage for each topology.

topo_counts_top = subset(topo_counts, topo.rank.overall <= 20)

bar_fills = as.character(topo_counts_top$topo.color)
names(bar_fills) = topo_counts_top$topo.rank.overall

cat(as.character(Sys.time()), " | Fig2A: topology count distribution\n")
fig_2a = ggplot(topo_counts_top, aes(x=reorder(topo.rank.overall, -topo.count.overall), y=topo.count.overall, fill=as.character(topo.rank.overall), color=outline, label=topo.prop)) +
  geom_bar(stat="identity", size=2) +
  geom_text(position=position_dodge(width=0.9), size=4, color="#333333", hjust=-0.15, vjust=-0.25, angle=45) +
  expand_limits(x = c(1,21)) +
  scale_y_continuous(expand=c(0, 0), limits=c(0,max(topo_counts$topo.count.overall)+7000)) +
  scale_fill_manual(name="", values=bar_fills) +
  scale_color_manual(name="", limits=topo_counts$outline, values=topo_counts$outline) +
  xlab("Topology rank") + 
  ylab("# of topologies") +
  bartheme() +
  theme(legend.position="none")
print(fig_2a)

######################

cat(as.character(Sys.time()), " | Fig2B: Plotting top trees\n")

top_topos = head(as.numeric(topo_counts$topo.rank.overall),max_tree_rank)
# Get the top ranking topologies for this chrome

tree_figs = list()
# A list to hold the three tree figs in

for(j in 1:length(top_topos)){
  # Generate figures for each of the top 3 topos
  cur_topo = top_topos[j]
  cur_tree_raw = as.character(all_windows_f[all_windows_f$topo.rank.overall==cur_topo,]$topo[1])
  
  cur_tree_raw = gsub("rdil", "R. dilectus", cur_tree_raw)
  cur_tree_raw = gsub("gdol", "G. dolicurus", cur_tree_raw)
  cur_tree_raw = gsub("mm10", "M. musculus", cur_tree_raw)
  cur_tree_raw = gsub("pdel", "P. delectorum", cur_tree_raw)
  cur_tree_raw = gsub("mnat", "M. natalensis", cur_tree_raw)
  cur_tree_raw = gsub("hall", "H. alleni", cur_tree_raw)
  cur_tree_raw = gsub("rsor", "R. soricoides", cur_tree_raw)
  
  cur_tree = read.tree(text=cur_tree_raw)
  #cur_color_cat = topo_counts$col.cat[topo_counts$topo.rank.overall==cur_topo]
  #cur_color_ind = which(cur_color_cat == bar_labels)
  cur_color = as.character(topo_counts$topo.color[topo_counts$topo.rank.overall==cur_topo])
  # Get all info for the current topo (color, string, etc.)
  
  #nodecheck(cur_tree, tree_type="object", xmax=10)
  
  cur_fig = ggtree(cur_tree, size=1.5, ladderize=T) + 
    xlim_tree(7) + 
    ggplot2::ylim(0, 8) +
    ggplot2::xlim(0, 18) +
    geom_tiplab(color="#333333", fontface='italic', size=4) +
    theme(plot.margin = unit(c(0,1,0,0), "cm"))
    #theme_tree2() + 
    #theme(panel.border=element_rect(color=cur_color, fill="NA", size=5))
  #print(cur_fig)
  # Generate the tree
  
  # if(j==3){
  #   cur_fig = cur_fig %>% rotate(9)
  # }
  
  tree_figs[[j]] = cur_fig
  # Save the tree fig in th elist
}

cat(as.character(Sys.time()), " | Fig2B: Combining tree figures\n")
fig_2b = plot_grid(plotlist=tree_figs, nrow=1, ncol=length(tree_figs), labels=c("Rank: 1", "2", "3"), label_size=14, hjust=1, align="h")
print(fig_2b)

cat(as.character(Sys.time()), " | Fig2: Combining panels A and B\n")
fig_2ab = plot_grid(fig_2a, 
                 fig_2b + theme(plot.margin = unit(c(0,0,0,0), "cm")), 
                 nrow=2, ncol=1, labels=c("A", "B"), label_size=16, aligh="vh")

######################

if(gen_chromo){
  cat(as.character(Sys.time()), " | Fig2C: Generating chromoplot\n")
  chrome_to_plot = "1"
  chrome_str = paste("chr", chrome_to_plot, sep="")
  chrdata = subset(all_windows, chr==chrome_str)
  chrdata_f = subset(all_windows_f, chr==chrome_str)
  
  chr_len = chrome_info[chrome_info$chr==chrome_str, ]$chr.len
  
  chrdata$col.cat = NA
  chrdata$ystart = NA
  chrdata$yend = NA
  cur_col = "#f2f2f2"
  topo_counts$topo.color = as.character(topo_counts$topo.color)
  for(i in 1:nrow(chrdata)) {
    cur_topo_num = chrdata[i,]$topo.num.overall
    if(is.na(cur_topo_num) || au_flag && chrdata[i,]$AU.test != "PASS"){
      cur_col = "#f2f2f2"
      chrdata[i,]$ystart = 0
      chrdata[i,]$yend = 4
    }else{
      cur_col = topo_counts$topo.color[topo_counts$topo.num.overall==cur_topo_num]
      #print(cur_col)
    }
    chrdata[i,]$col.cat = cur_col
  }
  
  cat(as.character(Sys.time()), " | Fig2C: Assigning coordinates\n")
  chrdata$ystart[chrdata$topo.rank.overall==1] = 4
  chrdata$yend[chrdata$topo.rank.overall==1] = 3
  
  chrdata$ystart[chrdata$topo.rank.overall==2] = 3
  chrdata$yend[chrdata$topo.rank.overall==2] = 2
  
  chrdata$ystart[chrdata$topo.rank.overall==3] = 2
  chrdata$yend[chrdata$topo.rank.overall==3] = 1
  
  chrdata$ystart[chrdata$topo.rank.overall>3] = 1
  chrdata$yend[chrdata$topo.rank.overall>3] = 0
  #chrdata$col.cat[chrdata$topo.rank.chrome>3] = "#999999"
  
  cols_labs = levels(as.factor(chrdata$col.cat))
  names(cols_labs) = levels(as.factor(chrdata$col.cat))
}

cat(as.character(Sys.time()), " | Fig2C: Rendering chromoplot\n")
fig_2c = ggplot(chrdata, aes(x=start, y=ystart, color=col.cat)) +
  geom_rect(aes(ymin=ystart, ymax=yend, xmin=start,xmax=end)) +
  #geom_segment(aes(x=start, y=ystart, xend=start, yend=yend)) +
  scale_color_manual(values=cols_labs) +
  scale_y_continuous(limits=c(0,4), breaks=0.5:3.5, labels=c("Other", "3", "2", "1")) +
  scale_x_continuous(limits=c(0,chr_len), breaks=NULL, expand=c(0,0)) +
  # geom_hline(yintercept=2,color="black") +
  ylab("Rank") +
  xlab("") +
  ggtitle("Chromosome 01") +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=16, angle=0, vjust=0.5),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        legend.position="none",
        legend.key.width = unit(0.75,  unit = "cm"),
        legend.spacing.x = unit(0.25, 'cm'),
        legend.title = element_blank(),
        legend.text=element_text(size=12),
        plot.title = element_text(hjust=0.5, size=16),
        plot.margin = unit(c(0,0,0,0), "cm")
  )
#print(fig_2c)

cat(as.character(Sys.time()), " | Fig2: Combining panels A, B, and C\n")
fig2 = plot_grid(fig_2a + theme(plot.margin = unit(c(0,0,0,0), "cm")), 
                 fig_2b + theme(plot.margin = unit(c(0,0,0,0), "cm")),
                 fig_2c + theme(plot.margin = unit(c(0,0,0,0), "cm")), 
                 nrow=3, labels=c("A", "B", "C"), label_size=16, align="vh")

print(fig2)

######################

if(save_fig){
  fig_file = here("figs", "fig2.png")
  cat(as.character(Sys.time()), " | Fig2: Saving figure:", fig_file, "\n")
  ggsave(filename=fig_file, fig2, width=7.5, height=10, units="in")
}