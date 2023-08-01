#############################################################################
# Snakemake rule to run MACSE on an input directory with CDS sequences
# Gregg Thomas, January 2023
#############################################################################

# snakemake -p -s 09_iqtree.smk --configfile 7spec_config.yaml --profile slurm_profile/ --dryrun

#############################################################################

import os

#############################################################################

DATASET = config["dataset"];
TREE_OUTDIR = config["tree_directory"];

ALNFILES_FILE = config["aln_filter_pass_file"];
FILTERFILE = config["aln_filter_file"];
ALNDIR = os.path.join(config["aln_directory"], "04-filter-spec4-seq20-site50", "cds");

aln_files_passed = open(ALNFILES_FILE, "r").read().split("\n");

loci = [];
first = True;
for line in open(FILTERFILE):

    if line.startswith("#"):
        continue;
    if first:
        first = False;
        continue;

    line = line.strip().split("\t");

    if int(line[16]) == 7 and int(line[19]) >= 4 and line[0] in aln_files_passed:
        loci.append(line[0].split("-")[0]);
    # Check to make sure the alignment has all 7 species and at least 4 seqs are unique
# for f in os.listdir(ALNDIR):
#     if not os.path.isfile(os.path.join(ALNDIR, f)) or not f.endswith("pretrim.macse.trim.filter.NT.fa"):
#         continue;

#     locus = f.split("-")[0];
#     loci.append(locus);

#print("# bulding trees for", len(loci), "alignments");

#############################################################################
# Final rule - rule that depends on final expected output file and initiates all
# the other rules

localrules: all

rule all:
    input:
        expand(os.path.join(TREE_OUTDIR, "loci", "{locus}", "{locus}-" + DATASET + "-cds-NT.pretrim.macse.trim.filter.NT.treefile"), locus=loci),
        os.path.join(TREE_OUTDIR, "loci.treefile"),
        os.path.join(TREE_OUTDIR, "astral", "astral.treefile"),
        os.path.join(TREE_OUTDIR, "astral", "astral.cf.tree"),
        os.path.join(TREE_OUTDIR, "concat", "concat.treefile"),
        os.path.join(TREE_OUTDIR, "concat", "concat.cf.tree")

        # gene_trees_rooted = expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci", "{window}", "{window}-mafft-trimal.treefile.rooted"), zip, chrome=CHROMES, window=WINDOWS),
        # astral_tree_rooted = expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "{chrome}-astral.treefile.rooted"), chrome=CHROMES),
        # astral_cf_tree_rooted = expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "astral", "{chrome}-astral.cf.tree.rooted"), chrome=CHROMES),
        # concat_tree_rooted = expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat.treefile.rooted"), chrome=CHROMES),
        # concat_cf_tree_rooted = expand(os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "concat", "{chrome}-concat.cf.tree.rooted"), chrome=CHROMES)

#############################################################################
# Pipeline rules

rule iqtree:
    input:
        trim_file = os.path.join(ALNDIR, "{locus}-" + DATASET + "-cds-NT.pretrim.macse.trim.filter.NT.fa")
    output:
        tree_file = os.path.join(TREE_OUTDIR, "loci", "{locus}", "{locus}-" + DATASET + "-cds-NT.pretrim.macse.trim.filter.NT.treefile"),
        #tree_file_rooted = os.path.join(TREE_OUTDIR, "{chrome}", SUB_DIR, "loci", "{window}", "{window}-mafft-trimal.treefile.rooted")
    log:
        os.path.join(TREE_OUTDIR, "loci", "logs", "{locus}-iqtree.log")
    params:
        bootstrap = "1000",
        prefix = os.path.join(TREE_OUTDIR, "loci", "{locus}", "{locus}-" + DATASET + "-cds-NT.pretrim.macse.trim.filter.NT"),
        #outgroup = OUTGROUP
    resources:
        cpus=1
    shell:
        """
        iqtree -s {input.trim_file} --prefix {params.prefix} -B {params.bootstrap} -T {resources.cpus} &> {log}        
        """
#nw_reroot {output.tree_file} {params.outgroup} > {output.tree_file_rooted}
# Run each locus through iqtree

######################

rule combine_gt:
    input:
        tree_file = expand(os.path.join(TREE_OUTDIR, "loci", "{locus}", "{locus}-" + DATASET + "-cds-NT.pretrim.macse.trim.filter.NT.treefile"), locus=loci)
    output:
        gene_trees_file = os.path.join(TREE_OUTDIR, "loci.treefile")
    params:
        loci_dir = os.path.join(TREE_OUTDIR, "loci")
    resources:
        cpus=1
    shell:
        """
        find {params.loci_dir} -name *.treefile -exec cat {{}} \\; > {output.gene_trees_file}
        """
# Combine the gene trees into a single file

######################

rule astral:
    input:
        gene_trees_file = os.path.join(TREE_OUTDIR, "loci.treefile")
    output:
        astral_tree = os.path.join(TREE_OUTDIR, "astral", "astral.treefile"),
        #astral_tree_rooted = os.path.join(TREE_OUTDIR, "astral", "{chrome}-astral.treefile.rooted")
    log:
        os.path.join(TREE_OUTDIR, "astral", "astral.log")
    # params:
    #     outgroup = OUTGROUP
    resources:
        mem="24g"
    shell:
        """
        astral -i {input.gene_trees_file} -o {output.astral_tree} 2> {log}
        """
#nw_reroot {output.astral_tree} {params.outgroup} > {output.astral_tree_rooted}
# Infer species tree with ASTRAL

######################

rule astral_cf:
    input:
        trim_file = expand(os.path.join(ALNDIR, "{locus}-" + DATASET + "-cds-NT.pretrim.macse.trim.filter.NT.fa"), locus=loci),
        trim_dir = ALNDIR,    
        astral_tree = os.path.join(TREE_OUTDIR, "astral", "astral.treefile"),
        gene_trees_file = os.path.join(TREE_OUTDIR, "loci.treefile"),
    output:
        concord_tree = os.path.join(TREE_OUTDIR, "astral", "astral.cf.tree"),
        #concord_tree_rooted = os.path.join(TREE_OUTDIR, "astral", "{chrome}-astral.cf.tree.rooted")
    log:
        os.path.join(TREE_OUTDIR, "astral", "astral-iqtree-cf.log")
    params:
        prefix = os.path.join(TREE_OUTDIR, "astral", "astral"),
        # outgroup = OUTGROUP
    resources:
        cpus=1,
        mem="24g"
    shell:
        """
        iqtree -t {input.astral_tree} --gcf {input.gene_trees_file} -p {input.trim_dir} --scf 100 -T {resources.cpus} --prefix {params.prefix} &> {log}
        """
#nw_reroot {output.concord_tree} {params.outgroup} > {output.concord_tree_rooted}        
# Calculate concordance factors on the astral tree

######################

rule concat:
    input:
        trim_file = expand(os.path.join(ALNDIR, "{locus}-" + DATASET + "-cds-NT.pretrim.macse.trim.filter.NT.fa"), locus=loci),
        trim_dir = ALNDIR
    output:
        concat_tree = os.path.join(TREE_OUTDIR, "concat", "concat.treefile"),
        #concat_tree_rooted = os.path.join(TREE_OUTDIR, "concat", "{chrome}-concat.treefile.rooted")
    log:
        os.path.join(TREE_OUTDIR, "concat", "concat-iqtree.log")
    params:
        bootstrap = "1000",
        prefix = os.path.join(TREE_OUTDIR, "concat", "concat"),
        #outgroup = OUTGROUP
    resources:
        cpus=12,
        time="12:00:00",
        mem="24g"
    shell:
        """
        iqtree -s {input.trim_dir} --prefix {params.prefix} -B {params.bootstrap} -T {resources.cpus} &> {log}
        """
#nw_reroot {output.concat_tree} {params.outgroup} > {output.concat_tree_rooted}        
# Run each locus through iqtree with concatenation

######################

rule concat_cf:
    input:
        trim_file = expand(os.path.join(ALNDIR, "{locus}-" + DATASET + "-cds-NT.pretrim.macse.trim.filter.NT.fa"), locus=loci),
        trim_dir = ALNDIR,    
        concat_tree = os.path.join(TREE_OUTDIR, "concat", "concat.treefile"),
        gene_trees_file = os.path.join(TREE_OUTDIR, "loci.treefile"),
    output:
        concord_tree = os.path.join(TREE_OUTDIR, "concat", "concat.cf.tree"),
        #concord_tree_rooted = os.path.join(TREE_OUTDIR, "concat", "{chrome}-concat.cf.tree.rooted")
    log:
        os.path.join(TREE_OUTDIR, "concat", "concat-iqtree-cf.log")
    params:
        prefix = os.path.join(TREE_OUTDIR, "concat", "concat"),
        #outgroup = OUTGROUP
    resources:
        cpus=1,
        mem="24g"
    shell:
        """
        iqtree -redo -t {input.concat_tree} --gcf {input.gene_trees_file} -p {input.trim_dir} --scf 100 -T {resources.cpus} --prefix {params.prefix} &> {log}
        # -redo needed because the checkpoint from the main tree search reports that the run has already been done
        """
#nw_reroot {output.concord_tree} {params.outgroup} > {output.concord_tree_rooted}
# Calculate concordance factors on the concatenated tree

######################

