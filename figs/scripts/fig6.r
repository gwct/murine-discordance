############################################################
# For penn genomes, 06.21
# Figure 6
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(shadowtext)
library(ape)
library(eulerr)
library(ggVennDiagram)
library(here)
source(here("figs", "scripts", "lib", "get_tree_info.r"))
source(here("figs", "scripts", "lib", "design.r"))

############################################################

save_fig = F
# Whether or not to save the figure

tree_file = here("summary-data", "03-Selection-tests", "concat.cf.rooted.tree")
transcript_windows_file = here("summary-data", "03-Selection-tests", "mm10-cds-windows.csv")
longest_transcripts_file = here("summary-data", "03-Selection-tests", "mm10.ensGene.chromes.longest.cds.bed")

concat_mg_local_file = here("summary-data", "03-Selection-tests", "mg94-local-st.csv")
gt_mg_local_file = here("summary-data", "03-Selection-tests", "mg94-local-gt.csv")

#concat_mg_local_file = here("summary-data", "03-Selection-tests", "no-pseudo-it", "hyphy", "penn-7spec-mg94-local-concat.csv")
#gt_mg_local_file = here("summary-data", "03-Selection-tests", "no-pseudo-it", "hyphy", "penn-7spec-mg94-local-gt.csv")

concat_m1a_file = here("summary-data", "03-Selection-tests", "m1a-st.csv")
concat_m2a_file = here("summary-data", "03-Selection-tests", "m2a-st.csv")
gt_m1a_file = here("summary-data", "03-Selection-tests", "m1a-gt.csv")
gt_m2a_file = here("summary-data", "03-Selection-tests", "m2a-gt.csv")

concat_absrel_file = here("summary-data", "03-Selection-tests", "absrel-st.csv")
gt_absrel_file = here("summary-data", "03-Selection-tests", "absrel-gt.csv")

concat_busted_file = here("summary-data", "03-Selection-tests", "busted-st.csv")
gt_busted_file = here("summary-data", "03-Selection-tests", "busted-gt.csv")

# Input options and files
######################

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

tree_info_concat$species = c("Rhynchomys soricoides", "Grammomys dolichurus", "Rhabdomys dilectus", "Mus musculus", "Praomys delectorum", 
                             "Mastomys natalensis", "Hylomyscus alleni", NA, NA, NA, NA, NA, NA)
# Add full species names to tree df

cat(as.character(Sys.time()), " | Reading gene trees\n")
#transcript_windows_stream = gzfile(transcript_windows_file, 'rt')  
transcript_windows = read.csv(transcript_windows_file, header=T, comment.char="#", stringsAsFactors=F)
# The 10kb windows that overlap with mouse transcripts

longest_transcripts = read.csv(longest_transcripts_file, sep="\t", header=F, comment.char="#", stringsAsFactors=F)
names(longest_transcripts) = c("chr", "start", "end", "transcript", "score", "strand", "cstart", "cend", "rgb", "num.exons", "exon.lens", "exon.starts")
# The transcripts used for gene trees and dN/dS/

twin = subset(transcript_windows, transcript %in% longest_transcripts$transcript)
# Only take windows that have one of the transcripts we used in them.

## Read tree and gene data
######################

concordant_gt = subset(twin, gene.tree.match.concat.tree==TRUE)
concordant_gt = unique(select(concordant_gt, transcript, chr, cds.start, cds.end))
# Get the concordant gene trees

discordant_gt = subset(twin, gene.tree.match.concat.tree==FALSE)
discordant_gt = unique(select(discordant_gt, transcript, chr, cds.start, cds.end))
# Get the discordant gene trees

perc_discordant_gt = nrow(discordant_gt) / (nrow(concordant_gt) + nrow(discordant_gt))
cat(as.character(Sys.time()), " | % of all gene trees that are discordant: ", perc_discordant_gt, "\n")

## Caclulate % of trees that are discordant with species tree
######################

cat(as.character(Sys.time()), " | Reading MG94 local data for dS filter\n")
concat_mg_local = read.csv(concat_mg_local_file, header=T, comment.char="#", stringsAsFactors=F)
concat_mg_local_gene = group_by(concat_mg_local, id) %>% summarize(dn=mean(dn, na.rm=T), ds=mean(ds, na.rm=T), dn.ds=mean(dn.ds, na.rm=T))
# Average dn and ds across all branches for each gene

gt_mg_local = read.csv(gt_mg_local_file, header=T, comment.char="#", stringsAsFactors=F)
gt_mg_local_gene = group_by(gt_mg_local, id) %>% summarize(dn=mean(dn, na.rm=T), ds=mean(ds, na.rm=T), dn.ds=mean(dn.ds, na.rm=T))
# Average dn and ds across all branches for each gene

concat_ds = ggplot(concat_mg_local_gene, aes(x=ds)) +
  geom_histogram(bins=50, fill=corecol(pal="wilke", numcol=1), color="#666666") +
  scale_y_continuous(expand=c(0,0)) +
  xlab("dS to species tree") +
  ylab("# of genes") +
  bartheme()
print(concat_ds)
# Distribution of dS when using concatenated tree

gt_ds = ggplot(gt_mg_local_gene, aes(x=ds)) +
  geom_histogram(bins=50, fill=corecol(pal="wilke", numcol=1), color="#666666") +
  scale_y_continuous(expand=c(0,0)) +
  xlab("dS to gene trees") +
  ylab("# of genes") +
  bartheme()
print(gt_ds)
# Distribution of dS when using gene trees

print(quantile(concat_mg_local_gene$ds, c(.90, .95, .99)))
#ds_filter_level = 0.0575
ds_filter_level = quantile(concat_mg_local_gene$ds, c(0.95))
cat(as.character(Sys.time()), " | Removing genes with dS above", ds_filter_level, "from subsequent analyses.\n")
# Genes above 95.5 percentile

concat_ds_filter = subset(concat_mg_local_gene, ds > ds_filter_level)
gt_ds_filter = subset(gt_mg_local_gene, ds > ds_filter_level)
ds_filter = intersect(concat_ds_filter$id, gt_ds_filter$id)

# Get a list of genes to filter out in subsequent analyses based on dS
## MG94 local
######################

cat(as.character(Sys.time()), " | Fig6 PAML M1 and M2: Reading, filtering, and doing LRT data\n")
concat_m1a = read.csv(concat_m1a_file, header=T, comment.char="#", stringsAsFactors=F)
concat_m1a = unique(select(concat_m1a, id, chr, start, end, lnl, k))
concat_m1a = concat_m1a[!concat_m1a$id %in% ds_filter,]
concat_m2a = read.csv(concat_m2a_file, header=T, comment.char="#", stringsAsFactors=F)
concat_m2a = unique(select(concat_m2a, id, chr, start, end, lnl, k))
concat_m2a = concat_m2a[!concat_m2a$id %in% ds_filter,]
# Read and filter the data

concat_paml = merge(x=concat_m1a, y=concat_m2a, by="id")
concat_paml$lrt = 2 * (concat_paml$lnl.y - concat_paml$lnl.x)
concat_paml$pval = pchisq(concat_paml$lrt, 2, lower.tail=F)
concat_paml$pval = p.adjust(concat_paml$pval, "fdr")
concat_paml_ps = subset(concat_paml, concat_paml$pval < 0.01)
# Perform the likelihood ratio test and p-value adjustment
## CONCAT DATA
##########

gt_m1a = read.csv(gt_m1a_file, header=T, comment.char="#", stringsAsFactors=F)
gt_m1a = unique(select(gt_m1a, id, chr, start, end, lnl, k))
gt_m1a = gt_m1a[!gt_m1a$id %in% ds_filter,]
gt_m2a = read.csv(gt_m2a_file, header=T, comment.char="#", stringsAsFactors=F)
gt_m2a = unique(select(gt_m2a, id, chr, start, end, lnl, k))
gt_m2a = gt_m2a[!gt_m2a$id %in% ds_filter,]
# Read and filter the data

gt_paml = merge(x=gt_m1a, y=gt_m2a, by="id")
gt_paml$lrt = 2 * (gt_paml$lnl.y - gt_paml$lnl.x)
gt_paml$pval = pchisq(gt_paml$lrt, 2, lower.tail=F)
gt_paml$pval = p.adjust(gt_paml$pval, "fdr")
gt_paml_ps = subset(gt_paml, gt_paml$pval < 0.01)
# Perform the likelihood ratio test and p-value adjustment
## GT DATA
##########

cat(as.character(Sys.time()), " | Fig6 PAML M1 and M2: Counting positively selected genes\n")
paml_gt_uniq = gt_paml_ps[!gt_paml_ps$id %in% concat_paml_ps$id,]
paml_gt_uniq_disco = length(intersect(paml_gt_uniq$id, discordant_gt$Gene.ID))
paml_gt_uniq_disco_perc = paml_gt_uniq_disco / nrow(paml_gt_uniq)
# GT unique
#####

paml_concat_uniq = concat_paml_ps[!concat_paml_ps$id %in% gt_paml_ps$id,]
paml_concat_uniq_disco = length(intersect(paml_concat_uniq$id, discordant_gt$Gene.ID))
paml_concat_uniq_disco_perc = paml_concat_uniq_disco / nrow(paml_concat_uniq)
# Concat unique
#####

paml_shared = concat_paml_ps[concat_paml_ps$id %in% gt_paml_ps$id,]
paml_shared_disco = length(intersect(paml_shared$id, discordant_gt$Gene.ID))
paml_shared_disco_perc = paml_shared_disco / nrow(paml_shared)
# Shared
#####

paml_comp = data.frame("type"=c("Gene trees", "Shared", "Concatenated tree"), 
                       "count"=c(nrow(paml_gt_uniq), nrow(paml_shared), nrow(paml_concat_uniq)),
                       "text"=c(paste(nrow(paml_gt_uniq), " (", signif(paml_gt_uniq_disco_perc, 2) * 100, "%)", sep=""), 
                                paste(nrow(paml_shared), " (", signif(paml_shared_disco_perc, 2) * 100, "%)", sep=""), 
                                paste(nrow(paml_concat_uniq), " (", signif(paml_concat_uniq_disco_perc, 2) * 100, "%)", sep=""))
                       )
# Compile paml data

paml_comp$prop = paml_comp$count / sum(paml_comp$count)
paml_comp$type = factor(paml_comp$type, levels=rev(c("Gene trees", "Shared", "Concatenated tree")))
paml_comp$label="M1a vs M2a"

# Combine the genes from different categories into one table
## M1a-M2a
######################

cat(as.character(Sys.time()), " | Fig6 aBSREL: Reading data\n")
concat_absrel = read.csv(concat_absrel_file, header=T, comment.char="#", stringsAsFactors=F)
concat_absrel = concat_absrel[!concat_absrel$id %in% ds_filter,]

gt_absrel = read.csv(gt_absrel_file, header=T, comment.char="#", stringsAsFactors=F)
gt_absrel = gt_absrel[!gt_absrel$id %in% ds_filter,]
# Read and filter data

cat(as.character(Sys.time()), " | Fig6 aBSREL: Counting positively selected genes\n")
concat_absrel_ps = subset(concat_absrel, num.branches.pval.less.than.alpha > 0)
gt_absrel_ps = subset(gt_absrel, num.branches.pval.less.than.alpha > 0)
# Select only genes with evidence for selection on at least one branch

absrel_gt_uniq = gt_absrel_ps[!gt_absrel_ps$id %in% concat_absrel_ps$id,]
absrel_gt_uniq_disco = length(intersect(absrel_gt_uniq$id, discordant_gt$Gene.ID))
absrel_gt_uniq_disco_perc = absrel_gt_uniq_disco / nrow(absrel_gt_uniq)
# GT unique
#####

absrel_concat_uniq = concat_absrel_ps[!concat_absrel_ps$id %in% gt_absrel_ps$id,]
absrel_concat_uniq_disco = length(intersect(absrel_concat_uniq$id, discordant_gt$Gene.ID))
absrel_concat_uniq_disco_perc = absrel_concat_uniq_disco / nrow(absrel_concat_uniq)
# Concat unique
#####

absrel_shared = concat_absrel_ps[concat_absrel_ps$id %in% gt_absrel_ps$id,]
absrel_shared_disco = length(intersect(absrel_shared$id, discordant_gt$Gene.ID))
absrel_shared_disco_perc = absrel_shared_disco / nrow(absrel_shared)
# Shared
#####

absrel_comp = data.frame("type"=c("Gene trees", "Shared", "Concatenated tree"),
                         "count"=c(nrow(absrel_gt_uniq), nrow(absrel_shared), nrow(absrel_concat_uniq)),
                         "text"=c(paste(nrow(absrel_gt_uniq), " (", signif(absrel_gt_uniq_disco_perc, 2) * 100, "%)", sep=""), 
                                  paste(nrow(absrel_shared), " (", signif(absrel_shared_disco_perc, 2) * 100, "%)", sep=""), 
                                  paste(nrow(absrel_concat_uniq), " (", signif(absrel_concat_uniq_disco_perc, 2) * 100, "%)", sep=""))
                          )
# Compile absrel data

absrel_comp$prop = absrel_comp$count / sum(absrel_comp$count)
absrel_comp$type = factor(absrel_comp$type, levels=rev(c("Gene trees", "Shared", "Concatenated tree")))

absrel_comp$label = "aBSREL"
# Combine the genes from different categories into one table
# aBSREL
######################

cat(as.character(Sys.time()), " | Fig6 BUSTED: Reading data\n")
concat_busted = read.csv(concat_busted_file, header=T, comment.char="#", stringsAsFactors=F)
concat_busted = concat_busted[!concat_busted$id %in% ds_filter,]
gt_busted = read.csv(gt_busted_file, header=T, comment.char="#", stringsAsFactors=F)
gt_busted = gt_busted[!gt_busted$id %in% ds_filter,]
# Read and filter the BUSTED data

cat(as.character(Sys.time()), " | Fig6 BUSTED: Counting positively selected genes\n")
#concat_busted_ps = subset(concat_busted, pval < 1 - ( 1-0.01)^(1/nrow(concat_busted)))
#gt_busted_ps = subset(gt_busted, pval < 1 - ( 1-0.01)^(1/nrow(gt_busted)))
# Dunn-Sidak correction: Genes with p-value less than corrected threshold inferred as PS

concat_busted$pval = p.adjust(concat_busted$pval, "fdr")
gt_busted$pval = p.adjust(gt_busted$pval, "fdr")
# Adjust p-values for multiple testing

concat_busted_ps = subset(concat_busted, pval < 0.01)
gt_busted_ps = subset(gt_busted, pval < 0.01)
# Genes with p-values less than 0.01 inferred as PS

busted_gt_uniq = gt_busted_ps[!gt_busted_ps$id %in% concat_busted_ps$id,]
busted_gt_uniq_disco = length(intersect(busted_gt_uniq$id, discordant_gt$Gene.ID))
busted_gt_uniq_disco_perc = busted_gt_uniq_disco / nrow(busted_gt_uniq)
# GT unique
#####

busted_concat_uniq = concat_busted_ps[!concat_busted_ps$id %in% gt_busted_ps$id,]
busted_concat_uniq_disco = length(intersect(busted_concat_uniq$id, discordant_gt$Gene.ID))
busted_concat_uniq_disco_perc = busted_concat_uniq_disco / nrow(busted_concat_uniq)
# Concat unique
#####

busted_shared = concat_busted_ps[concat_busted_ps$id %in% gt_busted_ps$id,]
busted_shared_disco = length(intersect(busted_shared$id, discordant_gt$Gene.ID))
busted_shared_disco_perc = busted_shared_disco / nrow(busted_shared)
# Shared
#####

busted_comp = data.frame("type"=c("Gene trees", "Shared", "Concatenated tree"), 
                         "count"=c(nrow(busted_gt_uniq), nrow(busted_shared), nrow(busted_concat_uniq)),
                         "text"=c(paste(nrow(busted_gt_uniq), " (", signif(busted_gt_uniq_disco_perc, 2) * 100, "%)", sep=""), 
                                   paste(nrow(busted_shared), " (", signif(busted_shared_disco_perc, 2) * 100, "%)", sep=""), 
                                   paste(nrow(busted_concat_uniq), " (", signif(busted_concat_uniq_disco_perc, 2) * 100, "%)", sep=""))
                                   )
# Compile the BUSTED data

busted_comp$prop = busted_comp$count / sum(busted_comp$count)
busted_comp$type = factor(busted_comp$type, levels=rev(c("Gene trees", "Shared", "Concatenated tree")))

busted_comp$label = "BUSTED"
# Combine the genes from different categories into one table
# BUSTED
######################

cat(as.character(Sys.time()), " | Fig6: Rendering panel A \n")
cols = c(corecol(pal="wilke", numcol=1), "#bda988", corecol(pal="wilke", numcol=1, offset=1))
names(cols) = c("Gene trees", "Shared", "Concatenated tree")
# Pick some colors

paml_comp$text2 = NA
paml_comp$text2 = ifelse(paml_comp$type=="Gene trees", paml_comp$text, NA)
paml_comp$text = ifelse(paml_comp$type=="Gene trees", NA, paml_comp$text)
paml_comp$text3 = NA
paml_comp$text4 = NA

absrel_comp$text2 = NA
absrel_comp$text3 = NA
absrel_comp$text3 = ifelse(absrel_comp$type=="Gene trees", absrel_comp$text, NA)
absrel_comp$text = ifelse(absrel_comp$type=="Gene trees", NA, absrel_comp$text)
absrel_comp$text4 = NA
absrel_comp$text4 = ifelse(absrel_comp$type=="Concatenated tree", absrel_comp$text, NA)
absrel_comp$text = ifelse(absrel_comp$type=="Concatenated tree", NA, absrel_comp$text)

busted_comp$text2 = NA
busted_comp$text3 = NA
busted_comp$text4 = NA
busted_comp$text4 = ifelse(busted_comp$type=="Concatenated tree", busted_comp$text, NA)
busted_comp$text = ifelse(busted_comp$type=="Concatenated tree", NA, busted_comp$text)
# Add the text labels

comps = rbind(paml_comp, absrel_comp, busted_comp)
# Combine the data

fig_6a = ggplot(comps, aes(x=label, y=prop, fill=type, label=text)) +
  geom_bar(stat="identity") +
  geom_text(size=2.3, position=position_stack(vjust=0.5), color="#f2f2f2") +
  geom_text(aes(label=text2), size=2.3, color="#f2f2f2", nudge_y=0.05) +
  geom_text(aes(label=text3), size=2.3, color="#f2f2f2", nudge_y=-0.05) +
  geom_text(aes(label=text4), size=2.3, color="#f2f2f2", nudge_y=0.82) +
  xlab("") +
  ylab("Proportion of positively selected genes") +
  scale_y_continuous(expand=c(0,0), limits=c(0,1)) +
  scale_x_discrete(limits=unique(comps$label)) +
  #scale_x_continuous(limits=c(-0.5,3)) +
  #scale_fill_manual(guide=guide_legend(reverse=T), values=cols) +
  scale_fill_manual(values=cols) +
  #scale_fill_discrete(, values=cols) +
  bartheme() +
  theme(legend.position="bottom",
        axis.title.x=element_text(size=12),
        axis.title.y=element_blank(),
        legend.text = element_text(size=10)) +
  coord_flip()
print(fig_6a)
# Plot the data

######################

if(save_fig){
  figfile = here("figs", "fig6.png")
  cat(as.character(Sys.time()), " | Fig6: Saving figure:", figfile, "\n")
  ggsave(filename=figfile, fig_6a, width=5, height=4, units="in")
}
# Save the figure

############################################################

## Old code for venn diagrams of genes between methods

# cat(as.character(Sys.time()), " | Fig6: Rendering panel B \n")
# gt_comp = list("BUSTED"=gt_busted_ps$id, "aBSREL"=gt_absrel_ps$id, "M1avM2a"=gt_paml_ps$id)
# #shared_gt_genes = Reduce(intersect, gt_comp)
# #write.csv(shared_gt_genes, file=paste(datadir, "shared-ps-genes-all-gt.csv", sep=""), col.names=F, row.names=F)
# #all_gt_genes = Reduce(union, gt_comp)
# #write.csv(all_gt_genes, file=paste(datadir, "ps-genes-all-gt.csv", sep=""), col.names=F, row.names=F)
# fig_6b = plot(euler(gt_comp, shape = "ellipse"), labels=FALSE, quantities=TRUE, fill=corecol(numcol=3, offset=5))
# print(fig_6b)
# 
# title = ggdraw() + draw_label("Gene trees", fontface='bold')
# fig_6b_t = plot_grid(title, fig_6b, ncol=1, rel_heights=c(0.1, 1))
# 
# ######################
# 
# cat(as.character(Sys.time()), " | Fig6: Rendering panel C \n")
# concat_comp = list("BUSTED"=concat_busted_ps$id, "aBSREL"=concat_absrel_ps$id, "M1avM2a"=concat_paml_ps$id)
# fig_6c = plot(euler(concat_comp, shape = "ellipse"), labels=FALSE, quantities=TRUE, fill=corecol(numcol=3, offset=5))
# #par(mar=c(4,4,4,4))
# print(fig_6c)
# 
# title = ggdraw() + draw_label("Concatenated tree", fontface='bold')
# fig_6c_t = plot_grid(title, fig_6c, ncol=1, rel_heights=c(0.1, 1))
# 
# cat(as.character(Sys.time()), " | Fig6: Combining panels B and C\n")
# fig_6bc = plot_grid(fig_6b_t, fig_6c_t, ncol=2, labels=c("B", "C"), label_size=16)
# 
# # gt_busted_list = data.frame("id"=gt_busted_ps$id, "test"="BUSTED", "label"="Gene trees")
# # gt_absrel_list = data.frame("id"=gt_absrel_ps$id, "test"="aBSREL", "label"="Gene trees")
# # gt_paml_list = data.frame("id"=gt_paml_ps$id, "test"="M1a vs. M2a", "label"="Gene trees")
# # 
# # concat_busted_list = data.frame("id"=concat_busted_ps$id, "test"="BUSTED", "label"="Concatenated tree")
# # concat_absrel_list = data.frame("id"=concat_absrel_ps$id, "test"="aBSREL", "label"="Concatenated tree")
# # concat_paml_list = data.frame("id"=concat_paml_ps$id, "test"="M1a vs. M2a", "label"="Concatenated tree")
# # 
# # combined_list = rbind(gt_busted_list, gt_absrel_list, gt_paml_list, concat_busted_list, concat_absrel_list, concat_paml_list)
# # 
# # 
# # plot(euler(combined_list, by = list(label)), legend = TRUE)
# 
# ######################
# 
# cat(as.character(Sys.time()), " | Fig6: Generating legend for panels B and C \n")
# 
# num_busted = data.frame("count"=length(concat_busted_ps$id))
# num_busted$label = "BUSTED"
# 
# num_absrel = data.frame("count"=length(concat_absrel_ps$id))
# num_absrel$label = "aBSREL"
# 
# num_paml = data.frame("count"=length(concat_paml_ps$id))
# num_paml$label = "M1a vs. M2a"
# 
# num_ps_long = rbind(num_busted, num_absrel, num_paml)
# 
# dummy_plot = ggplot(num_ps_long, aes(x=label, y=count, fill=label)) +
#   geom_bar(stat="identity") +
#   scale_fill_manual(values=corecol(numcol=3, offset=5)) +
#   bartheme() +
#   theme(legend.position="bottom")
# print(dummy_plot)
# 
# leg_bc = get_legend(dummy_plot)
# 
# fig_6bc_leg = plot_grid(fig_6bc, leg_bc, nrow=2, rel_heights=c(1,0.1))
# print(fig_6bc_leg)
# 
# ######################
# 
# cat(as.character(Sys.time()), " | Fig6: Combining with panel A\n")
# fig = plot_grid(fig_6a, fig_6bc_leg, nrow=2, labels=c("A", ""), label_size=16, rel_heights=c(1,0.8))
# print(fig_6a)

############################################################