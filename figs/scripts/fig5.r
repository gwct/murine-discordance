############################################################
# For rodent genomes
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
library(ggsignif)
library(phangorn)
library(broom)
library(here)
source(here("figs", "scripts", "lib", "design.r"))

########################################################################################################################

featureCalcs <- function(chrome, dists, window_size, max_dist_mb){
  cat(as.character(Sys.time()), "  | Fitting models\n", sep="")
  dists$adj = dists$adj * window_size / 1000000
  # Adjust the adjacenccy to reflect distance in Mb
  
  dists$adj[dists$adj < 0] = dists$adj[dists$adj < 0] * -1
  # Convert the negative adjacencies to positive
  
  dists = subset(dists, adj <= max_dist_mb)
  # Select only windows within the adjacency threshold we want to check
  
  adj_wrf = select(subset(dists, adj==1), feature.window, wrf)
  names(adj_wrf)[2] = "wrf.adjacent"
  # Get the wrfs from windows immediately adjacent to the features
  print(nrow(dists))
  fitted_models = dists %>% group_by(feature.window) %>% do(model = lm(wrf ~ adj, data = ., na.rm=T))
  # Fit a linear regression to every feature window for adjacencies out to the max distance
  print(nrow(fitted_models))
  fitted_models$int = NA
  fitted_models$slope = NA
  for(i in 1:nrow(fitted_models)){
    fitted_models[i,]$int = fitted_models$model[[i]]$coefficients[1]
    fitted_models[i,]$slope = fitted_models$model[[i]]$coefficients[2]
  }
  # Get the slope and intercept as separate columns
  
  fitted_models = fitted_models %>% select(feature.window, int, slope)
  # Select only the window, intercept, and slope columns from the models
  
  #####
  
  fitted_models_log = dists %>% group_by(feature.window) %>% do(model = lm(wrf ~ log10(adj), data = ., na.rm=T))
  # Fit a linear regression to every feature window for the LOG adjacencies out to the max distance
  
  fitted_models_log$log.int = NA
  fitted_models_log$log.slope = NA
  for(i in 1:nrow(fitted_models_log)){
    fitted_models_log[i,]$log.int = fitted_models_log$model[[i]]$coefficients[1]
    fitted_models_log[i,]$log.slope = fitted_models_log$model[[i]]$coefficients[2]
  }
  # Get the LOG slope and intercept as separate columns
  
  fitted_models_log = fitted_models_log %>% select(feature.window, log.int, log.slope)
  # Select only the window, intercept, and slope columns from the log models
  
  fitted_models = merge(fitted_models, fitted_models_log, by="feature.window")
  fitted_models = merge(fitted_models, adj_wrf, by="feature.window")
  # Merge everything together
  
  fitted_models$chr = chrome
  
  return(fitted_models)
}

featureToSpec <- function(spec_tree, fitted_models, chrdata_f){
  cat(as.character(Sys.time()), "  | Getting distances to species tree\n", sep="")
  
  window_trees = subset(chrdata_f, window %in% fitted_models$feature.window)
  window_trees = window_trees %>% select(window, unparsed.tree)
  # Get the windows with trees from the full data file
  
  names(window_trees)[1] = "feature.window"
  # Rename the window column to merge later
  
  window_trees$wrf.spec = wRF.dist(read.tree(text=window_trees$unparsed.tree), spec_tree)
  # Calculate wRF for all the feature window trees to the species tree
  
  fitted_models = merge(fitted_models, window_trees, by="feature.window")
  # Merge the models with the species tree distances
  
  fitted_models = fitted_models %>% select(!unparsed.tree)
  # Remove the tree from the dataframe

  return(fitted_models)
}


########################################################################################################################

window_size_kb = 10
# Window size in kb

window_size = window_size_kb * 1000
# Window size in bp

marker_window_size = 5
# Marker window size in Mb

au_flag = FALSE
# Set to filter out windows that don't pass the AU test

read_data = T
# Whether or not to read the data

do_calcs = F
# Whether or not to do the calculations (takes about 1hr)

gen_figs = T
# Set to generate figures

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
# Whether or not to save the supplemental figs for each chromosome

coding_genes = T
if(coding_genes){
  gene_str = "All coding genes"
}else{
  gene_str = "All genes"
}
# Whether or not to use coding genes or all genes

point_alpha = 0.1
# Point transparency level

datadir = here("summary-data", "02-genomic-windows")

distdir = "D:/data/rodent-genomes/dists/"

tree_file = here("summary-data", "02-genomic-windows", paste0(window_size_kb, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts.csv"))

species_tree_file = here("summary-data", "03-selection-tests", "concat.cf.rooted.tree")

chrome_info_file = here("summary-data", "02-genomic-windows", "recombination-markers", "chrome-stats.csv")

pre_calc_file = here("summary-data", "02-genomic-windows", paste0("all-feature-stats-", max_dist_mb, "mb.csv"))
# Input options
######################

if(read_data){
  cat(as.character(Sys.time()), " | Reading species tree\n")
  concat_tree = read.tree(species_tree_file)
  
  cat(as.character(Sys.time()), " | Reading tree window data: ", tree_file, "\n")
  all_windows = read.csv(tree_file, comment.char="#", header=T)
  all_windows_f = subset(all_windows, repeat.filter=="PASS" & missing.filter=="PASS")
  if(au_flag){
    all_windows_f = subset(all_windows_f, AU.test=="PASS")
  }
}
# Read and filter the window data
######################

if(do_calcs){
# Whether or not to do calculations
  
  non_features = data.frame("window"=c(), "wrf.adjacent"=c(), "slope"=c(), "int"=c(), "log.slope"=c(), "log.int"=c(), "wrf.spec"=c())
  hs = data.frame("window"=c(), "wrf.adjacent"=c(), "slope"=c(), "int"=c(), "log.slope"=c(), "log.int"=c(), "wrf.spec"=c())
  non_ps_genes = data.frame("window"=c(), "wrf.adjacent"=c(), "slope"=c(), "int"=c(), "log.slope"=c(), "log.int"=c(), "wrf.spec"=c())
  ps_genes = data.frame("window"=c(), "wrf.adjacent"=c(), "slope"=c(), "int"=c(), "log.slope"=c(), "log.int"=c(), "wrf.spec"=c())
  uces = data.frame("window"=c(), "wrf.adjacent"=c(), "slope"=c(), "int"=c(), "log.slope"=c(), "log.int"=c(), "wrf.spec"=c())
  # Initialize data frames for each feature type
  
  for(chrome in levels(as.factor(all_windows$chr))){
  # Do calculations for each chromosome
    if(skip_one && chrome != "chr9"){
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
    if(coding_genes){
      gene_base = paste0("coding-genes-", marker_window_size ,"-Mb.bed.", chrome, ".dists")
      gene_file = here(datadir, "coord-query", gene_base)
    }else{
      gene_base = paste0("genes-", marker_window_size ,"-Mb.bed.", chrome, ".dists")
      gene_file = here(datadir, "coord-query", gene_base)
    }
    uce_base = paste0("uces-", marker_window_size ,"-Mb.bed.", chrome, ".dists")
    uce_file = here(datadir, "coord-query", uce_base)
    
    hs_base = paste0("hotspots-", marker_window_size ,"-Mb.bed.", chrome, ".dists")
    hs_file = here(datadir, "coord-query", hs_base)

    ps_base = paste0("ps-", marker_window_size ,"-Mb.bed.", chrome, ".dists")
    ps_file = here(datadir, "coord-query", ps_base)    
    
    # File names for this chromosome
    ######################
   
    cat(as.character(Sys.time()), "  | Chromosome ", chrome, " - reading data: ", ps_file, "\n", sep="")
    ps_dists = read.csv(ps_file, header=T)
    ps_dists$feature.window = paste(ps_dists$chr, ":", ps_dists$start, "-", ps_dists$end, sep="")
    # Get the feature window from the other columns  
    # Read the data for this chromosome and feature
    
    ps_models = featureCalcs(out_chrome, ps_dists, window_size, max_dist_mb)
    ps_models = featureToSpec(concat_tree, ps_models, chrdata_f)
    ps_genes = rbind(ps_genes, ps_models)
  
    # Positively selected genes
    ######################
    
    cat(as.character(Sys.time()), "  | Chromosome ", chrome, " - reading data: ", gene_file, "\n", sep="")
    gene_dists = read.csv(gene_file, header=T)
    gene_dists$feature.window = paste(gene_dists$chr, ":", gene_dists$start, "-", gene_dists$end, sep="")
    # Get the feature window from the other columns  
    # Read the data for this chromosome and feature
    
    non_ps_models = featureCalcs(out_chrome, gene_dists, window_size, max_dist_mb)
    non_ps_models = featureToSpec(concat_tree, non_ps_models, chrdata_f)
    non_ps_genes = rbind(non_ps_genes, non_ps_models)
    
    non_ps_genes = subset(non_ps_genes, !feature.window %in% ps_genes$feature.window)
    # Non-positively selected genes
    ###################### 
    
    cat(as.character(Sys.time()), "  | Chromosome ", chrome, " - reading data: ", hs_file, "\n", sep="")
    hs_dists = read.csv(hs_file, header=T)
    hs_dists$feature.window = paste(hs_dists$chr, ":", hs_dists$start, "-", hs_dists$end, sep="")
    # Get the feature window from the other columns  
    # Read the data for this chromosome and feature
    
    hs_models = featureCalcs(out_chrome, hs_dists, window_size, max_dist_mb)
    hs_models = featureToSpec(concat_tree, hs_models, chrdata_f)
    hs = rbind(hs, hs_models)
    
    # Hotspots
    ######################
    
    cat(as.character(Sys.time()), "  | Chromosome ", chrome, " - reading data: ", uce_file, "\n", sep="")
    uce_dists = read.csv(uce_file, header=T)
    uce_dists$feature.window = paste(uce_dists$chr, ":", uce_dists$start, "-", uce_dists$end, sep="")
    # Get the feature window from the other columns  
    # Read the data for this chromosome and feature
    
    uce_models = featureCalcs(out_chrome, uce_dists, window_size, max_dist_mb)
    uce_models = featureToSpec(concat_tree, uce_models, chrdata_f)
    uces = rbind(uces, uce_models)
    
    # Hotspots
    ######################
  
    cat(as.character(Sys.time()), "  | Chromosome ", chrome, " - reading all chromosome windows\n", sep="")
    chr_dists = data.frame()
    for(i in 1:(max_dist_mb * 1000000 / window_size)){
      cat(i, " ")
      cur_dist_file = paste(distdir, "/all-dists/", chrome, "-", i, ".csv", sep="")
      cur_dists = read.csv(cur_dist_file, header=T)
      cur_dists$adj = i
      chr_dists = rbind(chr_dists, cur_dists)
    }
    cat("\n")
    # Read a file for each adjacency/dist for this chromosome
    
    all_feature_windows = unique(union(ps_models$feature.window, union(non_ps_models$feature.window, union(hs_models$feature.window, uce_models$feature.window))))
    names(chr_dists)[1] = "feature.window"
    chr_dists = subset(chr_dists, !feature.window %in% all_feature_windows)
    # Get only windows that don't overlap with a feature window
    
    chr_models = featureCalcs(out_chrome, chr_dists, window_size, max_dist_mb)
    chr_models = featureToSpec(concat_tree, chr_models, chrdata_f)
    non_features = rbind(non_features, chr_models)
    
    # Non-feature windows
    ######################
  } ## End Chromosome loop
}else{
  cat(as.character(Sys.time()), " | Reading pre-calculated feature data: ", pre_calc_file, "\n")
  features = read.csv(pre_calc_file, header=T)
} ## End do_calcs block
  
######################

non_lab = "Non-feature windows"
hs_lab = "Hotspot windows"
non_ps_lab = "Non-PS gene windows"
ps_lab = "PS gene windows"
uce_lab = "UCE windows"
# Feature column labels

features$label = factor(features$label, levels=c(non_lab, hs_lab, non_ps_lab, ps_lab, uce_lab))

if(do_calcs){
  cat(as.character(Sys.time()), " | Adding labels to full data\n")
  non_features$label = non_lab
  hs$label = hs_lab
  non_ps_genes$label = non_ps_lab
  ps_genes$label = ps_lab
  uces$label = uce_lab
  # Add labels
  
  features = rbind(non_features, hs, non_ps_genes, ps_genes, uces)
  #features$label = factor(features$label, levels=c(non_lab, hs_lab, non_ps_lab, ps_lab, uce_lab))
  write.csv(features, file=pre_calc_file, row.names=F)
  # Combine feature dfs
}
## If calcs are done, need to combine the data here

#features = subset(features, chr != "chr05" & chr != "chr10")

feature_cols = c('#006ddb', "#333333", '#db6d00', '#004949', '#920000')
feature_labs = c(non_lab, hs_lab, non_ps_lab, ps_lab, uce_lab)
# Colors for the feature distribution plots

sig_comp_lvls = c("Hotspot - Non-feature", "Non-PS gene - Non-feature", "PS gene - Non-feature", "UCE - Non-feature",
                  "Non-PS gene - Hotspot", "PS gene - Hotspot", "UCE - Hotspot",
                  "PS gene - Non-PS gene", "UCE - Non-PS gene",
                  "UCE - PS gene")
# Levels for the feature comparisons for the Tukey test plots

sig_cols = c("N.S."="#333333", 
             "*"=corecol(pal="wilke", numcol=1, offset=4), 
             "**"=corecol(pal="wilke", numcol=1, offset=2), 
             "***"=corecol(pal="wilke", numcol=1, offset=6))
# sig_cols = c("#333333", 
#              corecol(pal="wilke", numcol=1, offset=4), 
#              corecol(pal="wilke", numcol=1, offset=2), 
#              corecol(pal="wilke", numcol=1, offset=6))
# sig_labs = factor(c("N.S", "*", "**", "***"), levels=c("N.S", "*", "**", "***"))
# Colors for the significance levels in the Tukey means diff plots.

# Plot setup
######################

cat(as.character(Sys.time()), " | Rendering adjacent distributions\n")

wrf_adj_avgs = features %>% group_by(label) %>% summarize(avg.wrf.adjacent=mean(wrf.adjacent, na.rm=T))
# Get the average adjacent wrf for each feature

wrf_adj_plot = ggplot(features, aes(x=label, y=wrf.adjacent, group=label)) +
  geom_violin(aes(fill=label), alpha=0.2) +
  geom_boxplot(width=0.1, outlier.shape=NA) +
  geom_point(data=wrf_adj_avgs, aes(x=label, y=avg.wrf.adjacent, color=label)) +
  #geom_hline(data=wrf_adj_avgs, aes(yintercept=avg.wrf.adjacent, color=label)) +
  scale_color_manual(name="", values=feature_cols, labels=feature_labs, drop=FALSE) +
  scale_fill_manual(name="", values=feature_cols, labels=feature_labs, drop=FALSE) +
  xlab("") +
  ylab("wRF from window tree to\nadjacent window tree") +
  bartheme() +
  theme(axis.text.x=element_text(angle=40, hjust=1, size=10),
        axis.title.y=element_text(size=12),
        plot.margin=margin(0.5,0.5,0,0.5, unit="cm"),
        legend.position="none")
print(wrf_adj_plot)
# Plot the distributions of adjacent wrf, colored by feature

#####

cat(as.character(Sys.time()), " | Rendering adjacent mean diffs\n")

adj_wrf_anova = aov(wrf.adjacent ~ label, data=features)
adj_wrf_anova_comp = TukeyHSD(adj_wrf_anova)
# Perform ANOVA and Tukey mean difference test on results

adj_wrf_anova_comp = as.data.frame(adj_wrf_anova_comp[1])
adj_wrf_anova_comp$comp = row.names(adj_wrf_anova_comp)
# Convert Tukey results to df

adj_wrf_anova_comp$comp = gsub(" windows-", " - ", adj_wrf_anova_comp$comp)
adj_wrf_anova_comp$comp = gsub(" windows", "", adj_wrf_anova_comp$comp)
adj_wrf_anova_comp$comp[1] = "Hotspot - Non-feature"
# Remove some strings from the labels

adj_wrf_anova_comp$sig = "N.S."
adj_wrf_anova_comp[adj_wrf_anova_comp$label.p.adj <= 0.05,]$sig = "*"
adj_wrf_anova_comp[adj_wrf_anova_comp$label.p.adj <= 0.01,]$sig = "**"
adj_wrf_anova_comp[adj_wrf_anova_comp$label.p.adj <= 0.001,]$sig = "***"
# Add a significance label column

adj_wrf_anova_comp$comp = factor(adj_wrf_anova_comp$comp, levels=sig_comp_lvls)
# Re-order the comparisons by factoring with pre-ordered levels

adj_wrf_anova_plot = ggplot(adj_wrf_anova_comp, aes(x=comp, y=label.diff, color=sig)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=label.lwr, ymax=label.upr), width=0.25, size=1) +
  geom_hline(yintercept=0, size=1, linetype=2) +
  scale_color_manual(values=sig_cols) +
  #xlab("Feature window comparison\nof wRF to adjacent window tree") +
  xlab("") + 
  ylab("Difference in means") +
  bartheme() +
  theme(panel.grid.major.x=element_line(color = "#d3d3d3", size=0.5, linetype=1),
        axis.text.x=element_text(angle=40, hjust=1, size=10),
        axis.title.y=element_text(size=12),
        plot.margin=margin(0.5,0.5,0,0.5, unit="cm"))
  #coord_flip()
print(adj_wrf_anova_plot)
# Plot the differences in mean, colored by significance level

# Adjacent wrf plots
######################

cat(as.character(Sys.time()), " | Removing duplicate rows\n")

features = features %>% dplyr::select(!(wrf.adjacent))
features = unique(features)
features$label = factor(features$label, levels=c(non_lab, hs_lab, non_ps_lab, ps_lab, uce_lab))

# Removing wrf.adjacent and subsequent duplicate rows (since features have two adjacent windows)
######################

cat(as.character(Sys.time()), " | Rendering slope distributions\n")

slope_avgs = features %>% group_by(label) %>% summarize(avg.slope=mean(log.slope, na.rm=T))
# Get the average slope for each feature window

slope_plot = ggplot(features, aes(x=label, y=log.slope, group=label)) +
  geom_violin(aes(fill=label), alpha=0.1) +
  geom_boxplot(width=0.1, outlier.shape=NA) +
  geom_point(data=slope_avgs, aes(x=label, y=avg.slope, color=label)) +
  #geom_hline(data=slope_avgs, aes(yintercept=avg.slope, color=label)) +
  scale_color_manual(name="", values=feature_cols, labels=feature_labs, drop=FALSE) +
  scale_fill_manual(name="", values=feature_cols, labels=feature_labs, drop=FALSE) +
  xlab("") +
  ylab(paste("Log wRF slope to ", max_dist_mb, "Mb", sep="")) +
  bartheme() +
  theme(axis.text.x=element_text(angle=40, hjust=1, size=10),
        axis.title.y=element_text(size=12),
        plot.margin=margin(0.5,0.5,0,0.5, unit="cm"),
        legend.position="none")
print(slope_plot)
# Plot the distributions of slope, colored by feature

#####

cat(as.character(Sys.time()), " | Rendering slope mean diffs\n")

slope_anova = aov(log.slope ~ label, data=features)
slope_anova_comp = TukeyHSD(slope_anova)
# Perform ANOVA and Tukey mean difference test on results

slope_anova_comp = as.data.frame(slope_anova_comp[1])
slope_anova_comp$comp = row.names(slope_anova_comp)
# Convert Tukey results to df

slope_anova_comp$comp = gsub(" windows-", " - ", slope_anova_comp$comp)
slope_anova_comp$comp = gsub(" windows", "", slope_anova_comp$comp)
# Remove some strings from the labels

slope_anova_comp$sig = "N.S."
slope_anova_comp[slope_anova_comp$label.p.adj <= 0.05,]$sig = "*"
slope_anova_comp[slope_anova_comp$label.p.adj <= 0.01,]$sig = "**"
slope_anova_comp[slope_anova_comp$label.p.adj <= 0.001,]$sig = "***"
# Add a significance label column

slope_anova_comp$comp = factor(slope_anova_comp$comp, levels=sig_comp_lvls)
# Re-order the comparisons by factoring with pre-ordered levels

slope_anova_plot = ggplot(slope_anova_comp, aes(x=comp, y=label.diff, color=sig)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=label.lwr, ymax=label.upr), width=0.25, size=1) +
  geom_hline(yintercept=0, size=1, linetype=2) +
  scale_color_manual(values=sig_cols, drop=FALSE) +
  #xlab("wRF slope from feature window\nto all windows within 5Mb") +
  xlab("") +
  ylab("Difference in means") +
  bartheme() +
  theme(panel.grid.major.x=element_line(color = "#d3d3d3", size=0.5, linetype=1),
        axis.text.x=element_text(angle=40, hjust=1, size=10),
        axis.title.y=element_text(size=12),
        plot.margin=margin(0.5,0.5,0,0.5, unit="cm"))
  #coord_flip()
print(slope_anova_plot)
# Plot the differences in slope, colored by significance level

## Slope plots
######################

cat(as.character(Sys.time()), " | Rendering species tree distributions\n")

wrf_spec_avgs = features %>% group_by(label) %>% summarize(avg.wrf.spec=mean(wrf.spec, na.rm=T))
# Get the average wrf to species tree for each feature window

spec_plot = ggplot(features, aes(x=label, y=wrf.spec, group=label)) +
  geom_violin(aes(fill=label), alpha=0.1) +
  geom_boxplot(width=0.1, outlier.shape=NA) +
  geom_point(data=wrf_spec_avgs, aes(x=label, y=avg.wrf.spec, color=label)) +
  #geom_hline(data=wrf_spec_avgs, aes(yintercept=avg.wrf.spec, color=label)) +
  scale_color_manual(name="", values=feature_cols, labels=feature_labs, drop=FALSE) +
  scale_fill_manual(name="", values=feature_cols, labels=feature_labs, drop=FALSE) +
  xlab("") +
  ylab("wRF from window tree to\nspecies tree") +
  bartheme() +
  theme(axis.text.x=element_text(angle=40, hjust=1, size=10),
        axis.title.y=element_text(size=12),
        plot.margin=margin(0.5,0.5,0,0.5, unit="cm"),
        legend.position="none")
print(spec_plot)
# Plot the distributions of spec wrf, colored by feature

#####

cat(as.character(Sys.time()), " | Rendering species tree mean diffs\n")

spec_anova = aov(wrf.spec ~ label, data=features)
spec_anova_comp = TukeyHSD(spec_anova)
# Perform ANOVA and Tukey mean difference test on results

spec_anova_comp = as.data.frame(spec_anova_comp[1])
spec_anova_comp$comp = row.names(spec_anova_comp)
# Convert Tukey results to df

spec_anova_comp$comp = gsub(" windows-", " - ", spec_anova_comp$comp)
spec_anova_comp$comp = gsub(" windows", "", spec_anova_comp$comp)
# Remove some strings from the labels

spec_anova_comp$sig = "N.S."
spec_anova_comp[spec_anova_comp$label.p.adj <= 0.05,]$sig = "*"
spec_anova_comp[spec_anova_comp$label.p.adj <= 0.01,]$sig = "**"
spec_anova_comp[spec_anova_comp$label.p.adj <= 0.001,]$sig = "***"
# Add a significance label column

spec_anova_comp$comp = factor(spec_anova_comp$comp, levels=sig_comp_lvls)
# Re-order the comparisons by factoring with pre-ordered levels

spec_anova_plot = ggplot(spec_anova_comp, aes(x=comp, y=label.diff, color=sig)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=label.lwr, ymax=label.upr), width=0.25, size=1) +
  geom_hline(yintercept=0, size=1, linetype=2) +
  scale_color_manual(values=sig_cols, drop=FALSE) +
  #xlab("Feature window comparison\nof wRF to species tree") +
  xlab("") +
  ylab("Difference in means") +
  bartheme() +
  theme(panel.grid.major.x=element_line(color = "#d3d3d3", size=0.5, linetype=1),
        axis.text.x=element_text(angle=40, hjust=1, size=10),
        axis.title.y=element_text(size=12),
        plot.margin=margin(0.5,0.5,-0.75,0.5, unit="cm"))
  #coord_flip()
print(spec_anova_plot)
# Plot the differences in spec wrf, colored by significance level

## Spec wRF plots
######################

cat(as.character(Sys.time()), " | Combining panels\n")

comp_leg = get_legend(adj_wrf_anova_plot + theme(legend.position="bottom"))

#fig_5new_top = plot_grid(slope_plot, slope_anova_plot + theme(legend.position="none"), ncol=2, labels=c("A", ""))
fig_5new_mid = plot_grid(wrf_adj_plot, adj_wrf_anova_plot + theme(legend.position="none"), ncol=2, labels=c("A", ""))
fig_5new_bot = plot_grid(spec_plot, spec_anova_plot + theme(legend.position="none"), ncol=2, labels=c("B", ""))
#fig_5new_noleg = plot_grid(fig_5new_top, fig_5new_mid, fig_5new_bot, nrow=3)
fig_5new_noleg = plot_grid(fig_5new_mid, fig_5new_bot, nrow=2)
fig_5new = plot_grid(fig_5new_noleg, comp_leg, nrow=2, rel_heights=c(1,0.1))
# Combine figure panels

if(save_fig){
  figfile = here("figs", "fig5.png")
  cat(as.character(Sys.time()), " | Fig5: Saving figure:", figfile, "\n")
  ggsave(filename=figfile, fig_5new, width=7.5, height=8.5, units="in")
}
# Save figure
######################
######################

# spec_anova_comp = spec_anova_comp %>%
#   mutate(sig2 = ifelse(sig=="***", "Significant", sig))
# 
# sig_shapes = c("N.S."=21,
#              "*"=22,
#              "**"=23,
#              "Significant"=24)
# 
# sig_cols = c("N.S."="#333333", 
#              "*"=corecol(pal="wilke", numcol=1, offset=4), 
#              "**"=corecol(pal="wilke", numcol=1, offset=2), 
#              "Significant"=corecol(pal="wilke", numcol=1, offset=6))
# 
# spec_anova_plot = ggplot(spec_anova_comp, aes(x=comp, y=label.diff, color=sig2, fill=sig2, shape=sig2)) +
#   geom_point(size=3) +
#   geom_errorbar(aes(ymin=label.lwr, ymax=label.upr), width=0.25, size=1) +
#   geom_hline(yintercept=0, size=1, linetype=2) +
#   scale_color_manual(values=sig_cols, drop=FALSE) +
#   scale_fill_manual(values=sig_cols, drop=FALSE) +
#   scale_shape_manual(values=sig_shapes, drop=FALSE) +
#   #xlab("Feature window comparison\nof wRF to species tree") +
#   ggtitle("Comparing distributions of\nphylogenetic distance from trees inferred fromx\ngenomic features to the species tree") +
#   xlab("") +
#   ylab("Difference in means of distance\nfrom locus trees to species tree") +
#   bartheme() +
#   theme(panel.grid.major.x=element_line(color = "#d3d3d3", size=0.5, linetype=1),
#         axis.text.x=element_text(angle=40, hjust=1, size=14),
#         axis.title.y=element_text(size=14),
#         plot.title=element_text(size=16, margin=margin(b=20)),
#         plot.margin=margin(0.5,0.5,0,0.5, unit="cm"),
#         legend.position="bottom",
#         legend.margin = margin(-1, 0, 0, 0, unit="cm"))
# #coord_flip()
# print(spec_anova_plot)
# 
# ggsave(filename=here("figs", "scripts", "fig5c.png"), spec_anova_plot, width=6, height=6, units="in")

######################
######################
## STASH

# ps_genes_file = here("summary-data", "02-genomic-windows", "feature-beds", "ps-transcripts-all-gt.txt")
# #ps_genes_file = paste(ps_datadir, "/ps-genes-all-gt.csv", sep="")
# # PS gene IDs
# 
# cat(as.character(Sys.time()), " | Reading gene IDs for PS genes\n")
# ps_gene_ids = read.csv(ps_genes_file, header=F, stringsAsFactors=F)

# tree_file = paste(datadir, window_size_kb, "kb-0.5-0.5-", marker_window_size, "mb-topo-counts-tt.csv", sep="")
# all_windows = read.csv(tree_file, header=T)
# all_windows_f = subset(all_windows, repeat.filter=="PASS" & missing.filter=="PASS")
# if(au_flag){
#   all_windows_f = subset(all_windows_f, AU.test=="PASS")
# }
# 
# uce_windows = subset(all_windows, window %in% uces$feature.window)
# print(nrow(uce_windows[uce_windows$astral.chrome.topo,]) / nrow(uce_windows))
# 
# non_ps_windows = subset(all_windows, window %in% non_ps_genes$feature.window)
# print(nrow(non_ps_windows[non_ps_windows$astral.chrome.topo,]) / nrow(non_ps_windows))
# 
# ps_windows = subset(all_windows, window %in% ps_genes$feature.window)
# print(nrow(ps_windows[ps_windows$astral.chrome.topo,]) / nrow(ps_windows))
# 
# 
# hs_windows = subset(all_windows, window %in% hs$feature.window)
# print(nrow(hs_windows[hs_windows$astral.chrome.topo,]) / nrow(hs_windows))
# 
# non_windows = subset(all_windows, window %in% non_features$feature.window)
# print(nrow(non_windows[non_windows$astral.chrome.topo,]) / nrow(non_windows))
# 
# ######################
# 
# longest_transcript_windows = subset(transcript_windows, Transcript.ID %in% longest_transcripts$feature.id)
# longest_transcript_windows = select(longest_transcript_windows, window, Gene.ID, Transcript.ID, Gene.tree)
# longest_transcript_windows = unique(longest_transcript_windows)
# # Get the gene trees from the longest transcripts
# 
# #ps_genes_uniq = unique(select(ps_genes, !wrf.adjacent))
# ps_windows = subset(longest_transcript_windows, Gene.ID %in% ps_gene_ids$V1)
# ps_gene_trees = subset(ps_windows, window %in% ps_genes$feature.window)
# #ps_gene_trees = ps_windows %>% group_by(Gene.ID) %>% summarize("gene.tree"=Gene.tree)
# 
# absrel = data.frame("Gene.ID"=c(), "Transcript.ID"=c(), "Gene.tree"=c(), "rf"=c(), "wrf"=c())
# busted = data.frame("Gene.ID"=c(), "Transcript.ID"=c(), "Gene.tree"=c(), "rf"=c(), "wrf"=c())
# paml = data.frame("Gene.ID"=c(), "Transcript.ID"=c(), "Gene.tree"=c(), "rf"=c(), "wrf"=c())
# no_ps = data.frame("Gene.ID"=c(), "Transcript.ID"=c(), "Gene.tree"=c(), "rf"=c(), "wrf"=c())
# 
# longest_transcript_windows$rf = NA
# longest_transcript_windows$wrf = NA
# longest_transcript_windows$absrel = "No"
# longest_transcript_windows$busted = "No"
# longest_transcript_windows$paml = "No"
# # Initialize some new columns
# 
# concat_tree_relabel = concat_tree
# concat_tree_relabel[["tip.label"]] = c("rsor", "gdol", "rdil", "hall", "mnat", "pdel", "mmus")
# # Re-label the species tree to match the gene trees
# 
# for(i in 1:nrow(longest_transcript_windows)){
#   # Go over every transcript
#   
#   cur_tree = read.tree(text=longest_transcript_windows[i,]$Gene.tree)
#   # Get the gene tree
#   
#   longest_transcript_windows[i,]$rf = RF.dist(cur_tree, concat_tree_relabel)
#   longest_transcript_windows[i,]$wrf = wRF.dist(cur_tree, concat_tree_relabel)
#   # Calculate RF and wRF
#   
#   if(longest_transcript_windows[i,]$Gene.ID %in% gt_absrel_ps$id){
#     absrel = rbind(absrel, longest_transcript_windows[i,])
#     longest_transcript_windows[i,]$absrel = "Yes"
#   }
#   
#   if(longest_transcript_windows[i,]$Gene.ID %in% gt_busted_ps$id){
#     busted = rbind(busted, longest_transcript_windows[i,])
#     longest_transcript_windows[i,]$busted = "Yes"
#   }
#   
#   if(longest_transcript_windows[i,]$Gene.ID %in% gt_paml_ps$id){
#     paml = rbind(paml, longest_transcript_windows[i,])
#     longest_transcript_windows[i,]$paml = "Yes"
#   }
#   # Check if the current transcript had evidence for selection with any test
#   
#   if(longest_transcript_windows[i,]$absrel == "No" && longest_transcript_windows[i,]$busted == "No" && longest_transcript_windows[i,]$paml == "No"){
#     no_ps = rbind(no_ps, longest_transcript_windows[i,])
#   }
# }

#x_comps = list(c(non_lab, hs_lab), c(non_lab, non_ps_lab), c(non_lab, ps_lab), c(non_lab, uce_lab),
#               c(hs_lab, non_ps_lab), c(hs_lab, ps_lab), c(hs_lab, uce_lab), 
#               c(non_ps_lab, ps_lab), c(non_ps_lab, uce_lab), 
#               c(ps_lab, uce_lab))

## List of comparisons for the Wilcox test with geom_signif()

#geom_quasirandom(aes(color=ps.label), size=2, width=0.25, alpha=0.25) +
#geom_quasirandom(size=1, width=0.25, alpha=0.1, color="#666666") +
#geom_boxplot(outlier.shape=NA, alpha=0.3, width=0.5, fill="transparent", color="#000000") +
#geom_signif(comparisons=x_comps, test=wilcox.test, map_signif_level=TRUE, textsize=4, size=1, step_increase=0.15, margin_top=0.1) +
#scale_color_manual(values=corecol(pal="wilke", numcol=4, offset=2)) +
#scale_y_continuous(limits=c(0,0.9)) +

## The geom_signif() code for plots with significance brackets

#plot(spec_anova_comp, las=1)

## Base method to plot Tukey results

#spec_kw = kruskal.test(wrf.spec ~ label, data=features)
#spec_pw_wil = pairwise.wilcox.test(features$wrf.spec, features$label, p.adjust.method="BH")

## A median test like ANOVA

######################
######################










