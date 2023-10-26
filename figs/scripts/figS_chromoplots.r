############################################################
# For penn genomes, 01.22
# Generates chromoplots for each chromosome
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(ggplot2)
library(cowplot)
library(ggsignif)
library("ggtree")
library(here)
source(here("figs", "scripts", "lib", "get_tree_info.r"))
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

save_fig = T
# Whether or not to save the figure

skip_one = F
# Set to only do one test chromosome

max_tree_rank = 3
# The number of top trees to print/color

datadir = here("summary-data", "02-genomic-windows")

infile = here(datadir, paste(window_size, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts.csv", sep=""))
#infile = here(datadir, paste(window_size, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts-tt.csv.gz", sep=""))
#infile = here(datadir, paste0(window_size, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts-chr19-pseudo.csv"))
# Input options
######################

if(read_data){
  cat(as.character(Sys.time()), " | Reading data: ", infile, "\n")
  #all_windows = read.csv(gzfile(infile), header=T)
  all_windows = read.csv(infile, header=T, comment.char="#")
  all_windows_f = subset(all_windows, repeat.filter=="PASS" & missing.filter=="PASS")
  if(au_flag){
    all_windows_f = subset(all_windows_f, AU.test=="PASS")
  }
}
# Read and filter the data
######################

cat(as.character(Sys.time()), " | Generating legend theme\n")
color_temp = data.frame("blah"=c(), "rank"=c(), "col"=c())
for(color in levels(as.factor(all_windows_f$topo.color))){
  cur_rank = all_windows_f$topo.rank.overall[all_windows_f$topo.color==color][1]
  cur_tmp = data.frame("blah"=1, "rank"=cur_rank, "col"=color)
  color_temp = rbind(color_temp, cur_tmp)
}

color_temp$rank[color_temp$col=="#999999"] = "Other topologies"
color_temp = color_temp[order(color_temp$rank),]
color_temp$col = as.character(color_temp$col)
color_temp$rank = as.character(color_temp$rank)
# This generates a data frame with the colors for all topologies that finished top 3 in at least one chromosome

color_p = ggplot(color_temp, aes(x=rank, y=blah, fill=col, label=rank)) +
  geom_bar(stat="identity") +
  scale_fill_manual(name="Overall rank of topologies that are top 3 in at least one chromosome:", 
                    limits=color_temp$rank, 
                    values=color_temp$col) +
  bartheme()
#print(color_p)
# A dummy plot for the legend for fills and ranks

outline_temp = data.frame("blah"=c(1,1),
                          "topo.label"=c("Concatenated", "ASTRAL"),
                          "outline"=c("#000000","#24ff24")
                          )

ol_p = ggplot(outline_temp, aes(x=topo.label, y=blah, color=outline, label=topo.label)) +
  geom_bar(stat="identity", fill="transparent", size=2) +
  scale_color_manual(name="Chromosome level topologies (barplot only):", 
                     labels=outline_temp$topo.label, 
                     values=as.character(outline_temp$outline)) +
  bartheme()
#print(ol_p)
# A dummy plot for the legend for outline colors
# Initialize colors and legends
######################


cat(as.character(Sys.time()), " | Reading chromosomes\n")

max_tree_rank = 3
num_top_trees = 1:max_tree_rank
# The number of top trees we want to track

for(chrome in levels(as.factor(all_windows$chr))){
# Process every chromosome
  if(skip_one && chrome!="chr1"){
    next
  }
  if(nchar(chrome) == 4 && chrome != "chrX"){
    out_chrome = gsub("chr", "chr0", chrome)
  }else{
    out_chrome = chrome
  }
  # The string for the current chromosome.
  
  cat(as.character(Sys.time()), " | Starting chromosome ", out_chrome, "\n")
  
  chrdata = subset(all_windows, chr==chrome)
  total_windows = length(chrdata[,1])
  chrdata_f = subset(chrdata, repeat.filter=="PASS" & missing.filter=="PASS")
  if(au_flag){
    chrdata_f = subset(chrdata_f, AU.test=="PASS")
  }
  used_windows = length(chrdata_f[,1])
  # Subset and filter data for current chromosome.
  
  concat_chr_topo = chrdata_f$topo.num.chrome[chrdata_f$concat.chrome.topo][1]
  astral_chr_topo = chrdata_f$topo.num.chrome[chrdata_f$astral.chrome.topo][1]
  chr_len = chrdata$chr.len[1]
  num_topos = max(chrdata_f$topo.num.chrome)
  # Get the concatenated and astral chromosome topology and length
  
  window_title = paste("Rodent phylogenies: ", window_size, "kb windows on ", chrome, sep="")
  window_subtitle = paste("Chromosome length: ", chr_len, "bp, showing ", used_windows, " of ", total_windows, " windows, ", num_topos, " topologies", sep="")
  
  if(au_flag){
    outfile = here("figs", "supp", "chromoplots", paste(out_chrome, "-au-filter.png", sep=""))
  }else{
    outfile = here("figs", "supp", "chromoplots", paste(out_chrome, ".pdf", sep=""))
  }
  # Output information for this chrome
  
  topo_counts = unique(subset(chrdata_f, select=c(topo.num.chrome, topo.count.chrome, topo.rank.chrome, topo.num.overall, topo.rank.overall)))
  # Get the topology counts
  
  topo_counts$is.concat = FALSE
  topo_counts$is.concat[topo_counts$topo.num.chrome==concat_chr_topo] = TRUE
  
  topo_counts$is.astral = FALSE
  topo_counts$is.astral[topo_counts$topo.num.chrome==astral_chr_topo] = TRUE
  # Add a column to flag whether the topology matches the chromosome level topologies or not
  
  topo_counts = topo_counts[order(topo_counts$topo.count.chrome, decreasing=TRUE),]
  # Order the rows by topology count
  
  topo_counts$col = NA
  col_rank = 1
  for(i in topo_counts$topo.num.chrome){
    #cur_tree_bool = cur_win$unparsed.tree[cur_win$topo.num.chrome==i]
    #cur_tree_bool = cur_tree_bool[!is.na(cur_tree_bool)];
    #cur_tree = as.character(cur_tree_bool[1])
    # ??
    
    cur_col = as.character(chrdata_f$topo.color[chrdata_f$topo.num.chrome==i][1])
    topo_counts$col[topo_counts$topo.num.chrome==i] = cur_col
    if(cur_col != "#999999"){
      new_label = paste("Top ", col_rank, sep="")
      col_rank = col_rank + 1
    }else{
      new_label = "Other"
    }
    topo_counts$col.cat[topo_counts$topo.num.chrome==i] = new_label
  }
  # Add label and color columns to the topo counts
  
  topo_counts$outline = "transparent"
  #topo_counts$outline[topo_counts$is.astral] = "#24ff24"
  #topo_counts$outline[topo_counts$is.concat] = "#000000"
  # Add an outline column -- black outline is chrome topo
  
  bar_fills = unique(topo_counts$col)
  bar_labels = unique(topo_counts$col.cat)
  # Get the unique colors out into separate vectors for the plot
  
  topo_counts$topo.prop = round(topo_counts$topo.count / sum(topo_counts$topo.count), 2)
  # Percentage for each topology.
  
  cat(as.character(Sys.time()), " | ", out_chrome, "Plotting topology counts...\n")
  topo_counts_p = ggplot(subset(topo_counts, topo.rank.chrome<=20), aes(x=reorder(topo.rank.overall, -topo.count.chrome), y=topo.count.chrome, fill=col.cat, color=outline, label=topo.prop)) +
    geom_bar(stat="identity", size=2) +
    #geom_text(position=position_dodge(width=0.9), size=4, color="#333333", vjust=-0.25) +
    geom_text(position=position_dodge(width=0.9), size=4, color="#333333", hjust=-0.15, vjust=-0.25, angle=45) +
    scale_y_continuous(expand=c(0, 0), limits=c(0,max(topo_counts$topo.count.chrome)+700)) +
    scale_fill_manual(name="", limits=bar_labels, values=bar_fills) +
    scale_color_manual(name="", limits=topo_counts$outline, values=topo_counts$outline) +
    xlab("Overall topology rank") + 
    ylab("Topology count") +
    bartheme() +
    theme(legend.position="none")
  if(skip_one){
    print(topo_counts_p)
  }
  #stop()
  # Topology distribution
  ######################
  # Top 3 topologies 
  
  cat(as.character(Sys.time()), " | ", out_chrome, "Plotting top trees...\n")
  
  top_topos = head(as.numeric(topo_counts$topo.num.chrome),max_tree_rank)
  # Get the top ranking topologies for this chrome
  
  tree_figs = list()
  # A list to hold the three tree figs in
  
  for(j in 1:length(top_topos)){
  # Generate figures for each of the top 3 topos
    cur_topo = top_topos[j]
    cur_tree_raw = as.character(chrdata_f[chrdata_f$topo.num.chrome==cur_topo,]$topo[1])
    cur_tree = read.tree(text=cur_tree_raw)
    cur_color_cat = topo_counts$col.cat[topo_counts$topo.num.chrome==cur_topo]
    cur_color_ind = which(cur_color_cat == bar_labels)
    cur_color = topo_counts$col[topo_counts$topo.num.chrome==cur_topo]
    # Get all info for the current topo (color, string, etc.)
    
    cur_fig = ggtree(cur_tree, size=2.25) + 
      xlim_tree(7) + 
      ggplot2::ylim(0, 8) +
      #ggplot2::xlim(0, 4.5) +
      geom_tiplab(color="#333333", fontface='italic', size=6) +
      #theme_tree2() + 
      theme(panel.border=element_rect(color=cur_color, fill="NA", size=5))
    #print(cur_fig)
    # Generate the tree
    
    tree_figs[[j]] = cur_fig
    # Save the tree fig in th elist
  }
  
  tree_combo = plot_grid(plotlist=tree_figs, nrow=1, ncol=length(tree_figs), align="h")
  dist_tree_combo = plot_grid(topo_counts_p, tree_combo, nrow=2, labels=c("A","B"), label_size=20, align="v")
  if(skip_one){
    print(dist_tree_combo)
  }
  # Combine the tree figs, then combine with the bar plot

  # Top 3 topologies
  ######################
  # Chromoplot
  cat(as.character(Sys.time()), " | ", out_chrome, "Generating chromoplot...\n")
  
  chrdata$col.cat = NA
  chrdata$ystart = NA
  chrdata$yend = NA
  for(i in 1:nrow(chrdata)) {
    cur_topo_num = chrdata[i,]$topo.num.overall
    if(is.na(cur_topo_num) || au_flag && chrdata[i,]$AU.test != "PASS"){
      cur_col = "#f2f2f2"
      chrdata[i,]$ystart = 0
      chrdata[i,]$yend = 4
    }else{
      cur_col = topo_counts$col[topo_counts$topo.num.overall==cur_topo_num]
    }
    chrdata[i,]$col.cat = cur_col
  }
  
  chrdata$ystart[chrdata$topo.rank.chrome==1] = 4
  chrdata$yend[chrdata$topo.rank.chrome==1] = 3

  chrdata$ystart[chrdata$topo.rank.chrome==2] = 3
  chrdata$yend[chrdata$topo.rank.chrome==2] = 2

  chrdata$ystart[chrdata$topo.rank.chrome==3] = 2
  chrdata$yend[chrdata$topo.rank.chrome==3] = 1

  chrdata$ystart[chrdata$topo.rank.chrome>3] = 1
  chrdata$yend[chrdata$topo.rank.chrome>3] = 0
  #chrdata$col.cat[chrdata$topo.rank.chrome>3] = "#999999"
  
  cols_labs = levels(as.factor(chrdata$col.cat))
  names(cols_labs) = levels(as.factor(chrdata$col.cat))
  
  chr_p = ggplot(chrdata, aes(x=start, y=ystart, color=col.cat)) +
    geom_rect(aes(ymin=ystart, ymax=yend, xmin=start,xmax=end)) +
    #geom_segment(aes(x=start, y=ystart, xend=start, yend=yend)) +
    scale_color_manual(values=cols_labs) +
    scale_y_continuous(limits=c(0,4), breaks=0.5:3.5, labels=c("Other", "3", "2", "1")) +
    scale_x_continuous(limits=c(0,chr_len), breaks=NULL, expand=c(0,0)) +
    # geom_hline(yintercept=2,color="black") +
    ylab(chrome) +
    xlab("") +
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
  
  if(skip_one){
    print(chr_p)
  }
  # Chromoplot
  ######################
  # Render and combine figures 
  
  cat(as.character(Sys.time()), " | ", out_chrome, "Rendering and combining figures...\n")
  
  full_combo = plot_grid(topo_counts_p + theme(legend.position="none"), 
                         tree_combo, 
                         chr_p,
                         rel_heights=c(0.35,0.4,0.25),
                         nrow=3, labels=c("A", "B", "C"), label_size=20, align="v")
  
  psubtitle = ggdraw() + draw_label(window_subtitle, size=18, fontface="bold", x=0, hjust=0) + theme(plot.margin=margin(0,0,0,7))
  p_subtitle = plot_grid(psubtitle, full_combo, ncol=1, rel_heights=c(0.05,1))
  
  ptitle = ggdraw() + draw_label(window_title, size=24, fontface="bold", x=0, hjust=0) + theme(plot.margin=margin(0,0,0,7))
  p_title = plot_grid(ptitle, p_subtitle, ncol=1, rel_heights=c(0.03,1))
  # Add a subtitle and title
  
  legend_a = get_legend(color_p + theme(legend.direction="horizontal",
                                        legend.justification="center",
                                        legend.box.just="bottom",
                                        legend.title=element_text(size=16),
                                        legend.key = element_rect(size=5, color=NA)
                                        )
                        )
  p = plot_grid(p_title, legend_a, ncol=1, rel_heights=c(1, 0.075))
  
  #legend_b = get_legend(ol_p + theme(legend.direction="horizontal",
  #                                    legend.justification="center",
  #                                    legend.box.just="bottom",
  #                                    legend.title=element_text(size=16),
  #                                    legend.key = element_rect(size=5, color=NA)
  #                                    )
  #                      )
  #p = plot_grid(pa, legend_b, ncol=1, rel_heights=c(1, 0.075))
  # Extract the legends from the dummy color plots
  
  if(save_fig){
    cat(as.character(Sys.time()), " | ", out_chrome, "Saving combined figure: ", outfile, "\n")
    ggsave(filename=outfile, p, width=12, height=12, units="in")
  }
  cat("----------\n")
  if(skip_one){
    stop("skip one ok!")
  }
}