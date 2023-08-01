#############################################################################
# Pipeline for read mapping simulations with varying divergence
#############################################################################

import os
import re
from itertools import product

#############################################################################
# Example cmd

# snakemake -p -s 06_trees.smk --configfile rodent-windows-config.yaml --profile profiles/slurm_profile/ --keep-going --dryrun

#############################################################################
# Functions

#############################################################################
# Input and output info

work_dir = os.path.dirname(__file__);
#print("working dir: " + work_dir);

winsize_kb = config["winsize_kb"];
#print("winsize_kb: " + winsize_kb);

repeat_cutoff = config["repeat_cutoff"];
#print("repeat_cutoff: " + repeat_cutoff);

missing_cutoff = config["missing_cutoff"];
#print("missing_cutoff: " + missing_cutoff);

OUTGROUP = config["outgroup"];
INPUT_CHROMES = config["chromes"];

#print("reading window ids...");

TRIM_OUTDIR = config["trim_out_dir"];
TREE_OUTDIR = config["tree_out_dir"];
SUB_DIR = winsize_kb + "-" + repeat_cutoff + "-" + missing_cutoff;
# Dirs

######################

WINDOW_INPUT = os.path.join(config["data_dir"], SUB_DIR + "-windows-filter.tsv");
#print("window file: " + WINDOW_INPUT);

CHROMES, WINDOWS, first = [], [], True;
for line in open(WINDOW_INPUT):
    if line[0] == "#":
        continue;
    if first:
        first = False;
        continue;
    # Skip header lines in the window file

    line = line.strip().split("\t");
    if line[0] in INPUT_CHROMES and all( filters == "PASS" for filters in [line[5], line[7], line[8]] ):
        CHROMES.append(line[0]);
        WINDOWS.append(line[3].replace(":", "-"));
## Read the windows to align given the input chromes and the filters

#print(len(WINDOWS), "window ids read");

######################

#############################################################################
# Final rule - rule that depends on final expected output file and initiates all
# the other rules

localrules: all

rule all:
    input:
        gene_trees_rooted = expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci", "{window}", "{window}-mafft-trimal.treefile.rooted"), zip, chrome=CHROMES, window=WINDOWS),
        astral_tree_rooted = expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "{chrome}-astral.treefile.rooted"), chrome=CHROMES),
        astral_cf_tree_rooted = expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "{chrome}-astral.cf.tree.rooted"), chrome=CHROMES),
        concat_tree_rooted = expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat.treefile.rooted"), chrome=CHROMES),
        concat_cf_tree_rooted = expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat.cf.tree.rooted"), chrome=CHROMES)


        # expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci", "{window}", "{window}-mafft-trimal.treefile"), zip, chrome=CHROMES, window=WINDOWS),
        # expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci.treefile"), chrome=CHROMES),
        # expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "astral.treefile"), chrome=CHROMES),

        # expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat.treefile"), chrome=CHROMES)

#############################################################################
# Pipeline rules

rule iqtree:
    input:
        trim_file = os.path.join(TRIM_OUTDIR, "{chrome}", SUB_DIR, "{window}-mafft-trimal.fa")
    output:
        tree_file = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci", "{window}", "{window}-mafft-trimal.treefile"),
        tree_file_rooted = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci", "{window}", "{window}-mafft-trimal.treefile.rooted")
    log:
        os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci", "logs", "{window}-iqtree.log")
    params:
        bootstrap = "1000",
        prefix = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci", "{window}", "{window}-mafft-trimal"),
        outgroup = OUTGROUP
    resources:
        cpus=1
    shell:
        """
        iqtree -s {input.trim_file} --prefix {params.prefix} -B {params.bootstrap} -T {resources.cpus} &> {log}

        nw_reroot {output.tree_file} {params.outgroup} > {output.tree_file_rooted}
        """
# Run each locus through iqtree

######################

rule combine_gt:
    input:
        tree_file = expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci", "{window}", "{window}-mafft-trimal.treefile"), zip, chrome=CHROMES, window=WINDOWS)
    output:
        gene_trees_file = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci.treefile")
    params:
        loci_dir = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci")
    resources:
        cpus=12
    shell:
        """
        find {params.loci_dir} -name *.treefile -exec cat {{}} \\; > {output.gene_trees_file}
        """
# Combine the gene trees into a single file

######################

rule astral:
    input:
        gene_trees_file = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci.treefile")
    output:
        astral_tree = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "{chrome}-astral.treefile"),
        astral_tree_rooted = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "{chrome}-astral.treefile.rooted")
    log:
        os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "astral.log")
    params:
        outgroup = OUTGROUP
    resources:
        mem="24g"
    shell:
        """
        astral -i {input.gene_trees_file} -o {output.astral_tree} 2> {log}

        nw_reroot {output.astral_tree} {params.outgroup} > {output.astral_tree_rooted}
        """
# Infer species tree with ASTRAL

######################

rule astral_cf:
    input:
        trim_file = expand(os.path.join(TRIM_OUTDIR, "{{chrome}}", SUB_DIR, "{{window}}-mafft-trimal.fa"), zip, chrome=CHROMES, window=WINDOWS),
        trim_dir = os.path.join(TRIM_OUTDIR, "{chrome}", SUB_DIR),    
        astral_tree = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "{chrome}-astral.treefile"),
        gene_trees_file = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci.treefile"),
    output:
        concord_tree = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "{chrome}-astral.cf.tree"),
        concord_tree_rooted = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "{chrome}-astral.cf.tree.rooted")
    log:
        os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "{chrome}-astral-iqtree-cf.log")
    params:
        prefix = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "{chrome}-astral"),
        outgroup = OUTGROUP
    resources:
        cpus=1,
        mem="24g"
    shell:
        """
        iqtree -t {input.astral_tree} --gcf {input.gene_trees_file} -p {input.trim_dir} --scf 100 -T {resources.cpus} --prefix {params.prefix} &> {log}

        nw_reroot {output.concord_tree} {params.outgroup} > {output.concord_tree_rooted}
        """
# Calculate concordance factors on the astral tree

######################

rule concat:
    input:
        trim_file = expand(os.path.join(TRIM_OUTDIR, "{{chrome}}", SUB_DIR, "{{window}}-mafft-trimal.fa"), zip, chrome=CHROMES, window=WINDOWS),
        trim_dir = os.path.join(TRIM_OUTDIR, "{chrome}", SUB_DIR)
    output:
        concat_tree = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat.treefile"),
        concat_tree_rooted = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat.treefile.rooted")
    log:
        os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat-iqtree.log")
    params:
        bootstrap = "1000",
        prefix = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat"),
        outgroup = OUTGROUP
    resources:
        cpus=12,
        time="12:00:00",
        mem="24g"
    shell:
        """
        iqtree -s {input.trim_dir} --prefix {params.prefix} -B {params.bootstrap} -T {resources.cpus} &> {log}

        nw_reroot {output.concat_tree} {params.outgroup} > {output.concat_tree_rooted}
        """
# Run each locus through iqtree with concatenation

######################

rule concat_cf:
    input:
        trim_file = expand(os.path.join(TRIM_OUTDIR, "{{chrome}}", SUB_DIR, "{{window}}-mafft-trimal.fa"), zip, chrome=CHROMES, window=WINDOWS),
        trim_dir = os.path.join(TRIM_OUTDIR, "{chrome}", SUB_DIR),    
        concat_tree = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat.treefile"),
        gene_trees_file = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci.treefile"),
    output:
        concord_tree = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat.cf.tree"),
        concord_tree_rooted = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat.cf.tree.rooted")
    log:
        os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat-iqtree-cf.log")
    params:
        prefix = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat"),
        outgroup = OUTGROUP
    resources:
        cpus=1,
        mem="24g"
    shell:
        """
        iqtree -redo -t {input.concat_tree} --gcf {input.gene_trees_file} -p {input.trim_dir} --scf 100 -T {resources.cpus} --prefix {params.prefix} &> {log}
        # -redo needed because the checkpoint from the main tree search reports that the run has already been done

        nw_reroot {output.concord_tree} {params.outgroup} > {output.concord_tree_rooted}
        """
# Calculate concordance factors on the concatenated tree

######################

