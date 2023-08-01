#############################################################################
# Snakemake rule to run MACSE on an input directory with CDS sequences
# Gregg Thomas, January 2023
#############################################################################

# snakemake -p -s 05_aln_cds.smk --configfile 7spec_config.yaml --profile slurm_profile/ --dryrun

#############################################################################

# 2 loci would not run (program exits successfully but no output files are generated?) this step:
# ENSMUST00000092956
# ENSMUST00000064257

#############################################################################

import os

#############################################################################

DATASET = config["dataset"]
INFILE = config["locus_file"]
CDSDIR = config["cds_directory"]
BASE_OUTDIR = config["aln_directory"]

PRE_TRIMDIR = os.path.join(BASE_OUTDIR, "01-pre-trim")
ALNDIR = os.path.join(BASE_OUTDIR, "02-macse")
TRIMDIR = os.path.join(BASE_OUTDIR, "03-trim")

loci = [];
first = True;
for line in open(INFILE):
    if line[0] == "#":
        continue;
    if first:
        first = False;
        continue;

    line = line.strip().split("\t");
    locus = line[0];

    if all(cols == "Y" for cols in line[1:]):
        loci.append(locus);

print("# aligning", len(loci), "orthogroups");

#############################################################################
# Final rule - rule that depends on final expected output file and initiates all
# the other rules

localrules: all

rule all:
    input:
        expand(os.path.join(PRE_TRIMDIR, "nt", "{locus}-7spec-cds-NT.pretrim.fa"), locus=loci),
        expand(os.path.join(PRE_TRIMDIR, "aa", "{locus}-7spec-cds-AA.pretrim.fa"), locus=loci),
        expand(os.path.join(PRE_TRIMDIR, "trim-stats", "{locus}-7spec-cds.pretrim.csv"), locus=loci)
        # Output files from run_guidance

#############################################################################
# Pipeline rules

rule run_macse_pre_trim:
    input:
        os.path.join(CDSDIR, "{locus}-7spec-cds.fa")
    output:
        nt = os.path.join(PRE_TRIMDIR, "nt", "{locus}-7spec-cds-NT.pretrim.fa"),
        aa = os.path.join(PRE_TRIMDIR, "aa", "{locus}-7spec-cds-AA.pretrim.fa"),
        info = os.path.join(PRE_TRIMDIR, "trim-stats", "{locus}-7spec-cds.pretrim.csv"),
    params:
        locus = "{locus}",
    log:
        os.path.join(PRE_TRIMDIR, "logs", "{locus}-7spec-cds.pretrim.log")
    resources:
        cpus = 4,
	    time = "8:00:00"
    shell:
        """
        macse -prog trimNonHomologousFragments -seq {input} -out_NT {output.nt} -out_AA {output.aa} -out_trim_info {output.info} > {log} 2>&1
        """

#############################################################################
