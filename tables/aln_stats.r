############################################################
# For penn genomes, 08.23
#
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(tidyverse)
library(here)
source(here("figs", "scripts", "lib", "get_tree_info.r"))
source(here("figs", "scripts", "lib", "design.r"))

############################################################

aln_stats_file = here("summary-data", "03-selection-tests", "aln-stats-spec4-seq20-site50.log")
aln_stats = read_tsv(aln_stats_file, comment="#")

too_few_spec = aln_stats %>% filter(post.num.species != 7)
print(nrow(too_few_spec))

too_few_uniq = aln_stats %>% filter(post.uniq.seqs < 4)
print(nrow(too_few_uniq))
