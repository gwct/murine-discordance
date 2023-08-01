#############################################################################
# Snakemake rule to run HyPhy's aBSREL on an input directory with CDS alignments
# Gregg Thomas, January 2023
#############################################################################

# snakemake -p -s 10_selection_tests.smk --configfile 7spec_config.yaml --profile slurm_profile/ --dryrun

#############################################################################

import os

#############################################################################

SNAKEDIR = os.path.realpath(workflow.basedir);

DATASET = config["dataset"]
TREEDIR = config["tree_directory"]

ALNFILES_FILE = config["aln_filter_pass_file"]
FILTERFILE = config["aln_filter_file"];
ALNDIR = os.path.join(config["aln_directory"], "04-filter-spec4-seq20-site50", "cds-header-trimmed")

SPECTREE = config["species_tree_file"]
SELDIR = config["sel_test_outdir"]

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

print("# running selection tests on ", len(loci), " loci");

#############################################################################
# Final rule - rule that depends on final expected output file and initiates all
# the other rules

localrules: all

rule all:
    input:
        expand(os.path.join(SELDIR, "gt", "mg94-local", "{locus}.json"), locus=loci),
        expand(os.path.join(SELDIR, "st", "mg94-local", "{locus}.json"), locus=loci)

        #expand(os.path.join(SELDIR, "st", "m2a", "{locus}", "codeml.ctl"), locus=loci),
        #expand(os.path.join(SELDIR, "st", "m2a", "{locus}", "codeml.out"), locus=loci),
        # Expected output for paml m2a with species tree

        #expand(os.path.join(SELDIR, "st", "m1a", "{locus}", "codeml.ctl"), locus=loci),
        #expand(os.path.join(SELDIR, "st", "m1a", "{locus}", "codeml.out"), locus=loci)
        # Expected output for paml m1a with species tree

        #expand(os.path.join(SELDIR, "gt", "m2a", "{locus}", "codeml.ctl"), locus=loci),
        #expand(os.path.join(SELDIR, "gt", "m2a", "{locus}", "codeml.out"), locus=loci)
        # Expected output for paml m2a with gene trees

        #expand(os.path.join(SELDIR, "gt", "m1a", "{locus}", "codeml.ctl"), locus=loci),
        #expand(os.path.join(SELDIR, "gt", "m1a", "{locus}", "codeml.out"), locus=loci)
        # Expected output for paml m1a with gene trees

        #expand(os.path.join(SELDIR, "gt", "busted", "{locus}.json"), locus=loci),
        #expand(os.path.join(SELDIR, "st", "busted", "{locus}.json"), locus=loci),
        # Expected output from rule busted

        #expand(os.path.join(SELDIR, "gt", "absrel", "{locus}.json"), locus=loci),
        #expand(os.path.join(SELDIR, "st", "absrel", "{locus}.json"), locus=loci),
        # Expected output from rule absrel



#############################################################################
# Pipeline rules

rule mg94_gt:
    input:
        aln = os.path.join(ALNDIR, "{locus}-" + DATASET + "-cds-NT.pretrim.macse.trim.filter.NT.trim.fa"),
        tree = os.path.join(TREEDIR, "loci", "{locus}", "{locus}-" + DATASET + "-cds-NT.pretrim.macse.trim.filter.NT.treefile")
        #tree = SPECTREE
    output:
        os.path.join(SELDIR, "gt", "mg94-local", "{locus}.json")
    params:
        mg94 = config["mg94_path"]
    log:
        os.path.join(SELDIR, "gt", "logs", "mg94-local", "{locus}.log")
    shell:
        """
        hyphy {params.mg94} --alignment {input.aln} --tree {input.tree} --type local --output {output} &> {log}
        """
# Run each locus through mg94 with the gene trees

#############################################################################

rule mg94_st:
    input:
        aln = os.path.join(ALNDIR, "{locus}-" + DATASET + "-cds-NT.pretrim.macse.trim.filter.NT.trim.fa"),
        tree = SPECTREE
    output:
        os.path.join(SELDIR, "st", "mg94-local", "{locus}.json")
    params:
        mg94 = config["mg94_path"]
    log:
        os.path.join(SELDIR, "st", "logs", "mg94-local", "{locus}.log")
    shell:
        """
        hyphy {params.mg94} --alignment {input.aln} --tree {input.tree} --type local --output {output} &> {log}
        """
# Run each locus through mg94 with the species tree

#############################################################################

rule absrel_gt:
    input:
        aln = os.path.join(ALNDIR, "{locus}-" + DATASET + "-cds-NT.pretrim.macse.trim.filter.NT.trim.fa"),
        tree = os.path.join(TREEDIR, "loci", "{locus}", "{locus}-" + DATASET + "-cds-NT.pretrim.macse.trim.filter.NT.treefile")
        #tree = SPECTREE
    output:
        os.path.join(SELDIR, "gt", "absrel", "{locus}.json")
    log:
        os.path.join(SELDIR, "gt", "logs", "absrel", "{locus}.log")
    shell:
        """
        hyphy absrel --alignment {input.aln} --tree {input.tree} --output {output} &> {log}
        """
# Run each locus through absrel with the gene trees

#############################################################################

rule absrel_st:
    input:
        aln = os.path.join(ALNDIR, "{locus}-" + DATASET + "-cds-NT.pretrim.macse.trim.filter.NT.trim.fa"),
        tree = SPECTREE
    output:
        os.path.join(SELDIR, "st", "absrel", "{locus}.json")
    log:
        os.path.join(SELDIR, "st", "logs", "absrel", "{locus}.log")
    shell:
        """
        hyphy absrel --alignment {input.aln} --tree {input.tree} --output {output} &> {log}
        """
# Run each locus through absrel with the species tree

#############################################################################

rule busted_gt:
    input:
        aln = os.path.join(ALNDIR, "{locus}-" + DATASET + "-cds-NT.pretrim.macse.trim.filter.NT.trim.fa"),
        tree = os.path.join(TREEDIR, "loci", "{locus}", "{locus}-" + DATASET + "-cds-NT.pretrim.macse.trim.filter.NT.treefile")
        #tree = SPECTREE
    output:
        os.path.join(SELDIR, "gt", "busted", "{locus}.json")
    log:
        os.path.join(SELDIR, "gt", "logs", "busted", "{locus}.log")
    shell:
        """
        hyphy busted --alignment {input.aln} --tree {input.tree} --output {output} &> {log}
        """
# Run each locus through busted with the gene trees

#############################################################################

rule busted_st:
    input:
        aln = os.path.join(ALNDIR, "{locus}-" + DATASET + "-cds-NT.pretrim.macse.trim.filter.NT.trim.fa"),
        tree = SPECTREE
    output:
        os.path.join(SELDIR, "st", "busted", "{locus}.json")
    log:
        os.path.join(SELDIR, "st", "logs", "busted", "{locus}.log")
    shell:
        """
        hyphy busted --alignment {input.aln} --tree {input.tree} --output {output} &> {log}
        """
# Run each locus through absrel with the species tree

#############################################################################

rule m1a_gt_gen:
    input:
        aln = os.path.join(ALNDIR, "{locus}-" + DATASET + "-cds-NT.pretrim.macse.trim.filter.NT.trim.fa"),
        tree = os.path.join(TREEDIR, "loci", "{locus}", "{locus}-" + DATASET + "-cds-NT.pretrim.macse.trim.filter.NT.treefile"),
        tmplate = os.path.join(SNAKEDIR, "templates", "paml-m1a-template.txt")
    output:
        paml_dir = directory(os.path.join(SELDIR, "gt", "m1a", "{locus}")),
        paml_aln = os.path.join(SELDIR, "gt", "m1a", "{locus}", "codeml.fa"),
        paml_tre = os.path.join(SELDIR, "gt", "m1a", "{locus}", "codeml.tre"),
        paml_ctl = os.path.join(SELDIR, "gt", "m1a", "{locus}", "codeml.ctl")
    params:
        paml_out = os.path.join(SELDIR, "gt", "m1a", "{locus}", "codeml.out")
    log:
        os.path.join(SELDIR, "gt", "logs", "m1a", "{locus}.log")
    run:
        import shutil

        ctl_str = open(input.tmplate, "r").read();
        with open(output.paml_ctl, "w") as ctl_outfile:
            ctl_outfile.write(ctl_str.format(seq=output.paml_aln, tree=output.paml_tre, paml_out=params.paml_out));

        shutil.copy(input.aln, output.paml_aln);
        shutil.copy(input.tree, output.paml_tre);

#############################################################################

rule m1a_gt:
    input:
        paml_dir = os.path.join(SELDIR, "gt", "m1a", "{locus}"),
        paml_aln = os.path.join(SELDIR, "gt", "m1a", "{locus}", "codeml.fa"),
        paml_tre = os.path.join(SELDIR, "gt", "m1a", "{locus}", "codeml.tre"),
        paml_ctl = os.path.join(SELDIR, "gt", "m1a", "{locus}", "codeml.ctl")
    output:
        paml_out = os.path.join(SELDIR, "gt", "m1a", "{locus}", "codeml.out")
    log:
        os.path.join(SELDIR, "gt", "m1a", "{locus}", "codeml.log")
    shell:
        """
        cd {input.paml_dir}
        codeml > {log}
        """
#############################################################################

rule m2a_gt_gen:
    input:
        aln = os.path.join(ALNDIR, "{locus}-" + DATASET + "-cds-NT.pretrim.macse.trim.filter.NT.trim.fa"),
        tree = os.path.join(TREEDIR, "loci", "{locus}", "{locus}-" + DATASET + "-cds-NT.pretrim.macse.trim.filter.NT.treefile"),
        tmplate = os.path.join(SNAKEDIR, "templates", "paml-m2a-template.txt")
    output:
        paml_dir = directory(os.path.join(SELDIR, "gt", "m2a", "{locus}")),
        paml_aln = os.path.join(SELDIR, "gt", "m2a", "{locus}", "codeml.fa"),
        paml_tre = os.path.join(SELDIR, "gt", "m2a", "{locus}", "codeml.tre"),
        paml_ctl = os.path.join(SELDIR, "gt", "m2a", "{locus}", "codeml.ctl")
    params:
        paml_out = os.path.join(SELDIR, "gt", "m2a", "{locus}", "codeml.out")
    log:
        os.path.join(SELDIR, "gt", "logs", "m2a", "{locus}.log")
    run:
        import shutil

        ctl_str = open(input.tmplate, "r").read();
        with open(output.paml_ctl, "w") as ctl_outfile:
            ctl_outfile.write(ctl_str.format(seq=output.paml_aln, tree=output.paml_tre, paml_out=params.paml_out));

        shutil.copy(input.aln, output.paml_aln);
        shutil.copy(input.tree, output.paml_tre);

#############################################################################

rule m2a_gt:
    input:
        paml_dir = os.path.join(SELDIR, "gt", "m2a", "{locus}"),
        paml_aln = os.path.join(SELDIR, "gt", "m2a", "{locus}", "codeml.fa"),
        paml_tre = os.path.join(SELDIR, "gt", "m2a", "{locus}", "codeml.tre"),
        paml_ctl = os.path.join(SELDIR, "gt", "m2a", "{locus}", "codeml.ctl")
    output:
        paml_out = os.path.join(SELDIR, "gt", "m2a", "{locus}", "codeml.out")
    log:
        os.path.join(SELDIR, "gt", "m2a", "{locus}", "codeml.log")
    shell:
        """
        cd {input.paml_dir}
        codeml > {log}
        """
#############################################################################

rule m1a_st_gen:
    input:
        aln = os.path.join(ALNDIR, "{locus}-" + DATASET + "-cds-NT.pretrim.macse.trim.filter.NT.trim.fa"),
        tree = SPECTREE,
        tmplate = os.path.join(SNAKEDIR, "templates", "paml-m1a-template.txt")
    output:
        paml_dir = directory(os.path.join(SELDIR, "st", "m1a", "{locus}")),
        paml_aln = os.path.join(SELDIR, "st", "m1a", "{locus}", "codeml.fa"),
        paml_tre = os.path.join(SELDIR, "st", "m1a", "{locus}", "codeml.tre"),
        paml_ctl = os.path.join(SELDIR, "st", "m1a", "{locus}", "codeml.ctl")
    params:
        paml_out = os.path.join(SELDIR, "st", "m1a", "{locus}", "codeml.out")
    log:
        os.path.join(SELDIR, "st", "logs", "m1a", "{locus}.log")
    run:
        import shutil

        ctl_str = open(input.tmplate, "r").read();
        with open(output.paml_ctl, "w") as ctl_outfile:
            ctl_outfile.write(ctl_str.format(seq=output.paml_aln, tree=output.paml_tre, paml_out=params.paml_out));

        shutil.copy(input.aln, output.paml_aln);
        shutil.copy(input.tree, output.paml_tre);

#############################################################################

rule m1a_st:
    input:
        paml_dir = os.path.join(SELDIR, "st", "m1a", "{locus}"),
        paml_aln = os.path.join(SELDIR, "st", "m1a", "{locus}", "codeml.fa"),
        paml_tre = os.path.join(SELDIR, "st", "m1a", "{locus}", "codeml.tre"),
        paml_ctl = os.path.join(SELDIR, "st", "m1a", "{locus}", "codeml.ctl")
    output:
        paml_out = os.path.join(SELDIR, "st", "m1a", "{locus}", "codeml.out")
    log:
        os.path.join(SELDIR, "st", "m1a", "{locus}", "codeml.log")
    shell:
        """
        cd {input.paml_dir}
        codeml > {log}
        """

#############################################################################

rule m2a_st_gen:
    input:
        aln = os.path.join(ALNDIR, "{locus}-" + DATASET + "-cds-NT.pretrim.macse.trim.filter.NT.trim.fa"),
        tree = SPECTREE,
        tmplate = os.path.join(SNAKEDIR, "templates", "paml-m2a-template.txt")
    output:
        paml_dir = directory(os.path.join(SELDIR, "st", "m2a", "{locus}")),
        paml_aln = os.path.join(SELDIR, "st", "m2a", "{locus}", "codeml.fa"),
        paml_tre = os.path.join(SELDIR, "st", "m2a", "{locus}", "codeml.tre"),
        paml_ctl = os.path.join(SELDIR, "st", "m2a", "{locus}", "codeml.ctl")
    params:
        paml_out = os.path.join(SELDIR, "st", "m2a", "{locus}", "codeml.out")
    log:
        os.path.join(SELDIR, "st", "logs", "m2a", "{locus}.log")
    run:
        import shutil

        ctl_str = open(input.tmplate, "r").read();
        with open(output.paml_ctl, "w") as ctl_outfile:
            ctl_outfile.write(ctl_str.format(seq=output.paml_aln, tree=output.paml_tre, paml_out=params.paml_out));

        shutil.copy(input.aln, output.paml_aln);
        shutil.copy(input.tree, output.paml_tre);

#############################################################################

rule m2a_st:
    input:
        paml_dir = os.path.join(SELDIR, "st", "m2a", "{locus}"),
        paml_aln = os.path.join(SELDIR, "st", "m2a", "{locus}", "codeml.fa"),
        paml_tre = os.path.join(SELDIR, "st", "m2a", "{locus}", "codeml.tre"),
        paml_ctl = os.path.join(SELDIR, "st", "m2a", "{locus}", "codeml.ctl")
    output:
        paml_out = os.path.join(SELDIR, "st", "m2a", "{locus}", "codeml.out")
    log:
        os.path.join(SELDIR, "st", "m2a", "{locus}", "codeml.log")
    shell:
        """
        cd {input.paml_dir}
        codeml > {log}
        """

#############################################################################