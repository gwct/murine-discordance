#############################################################################
# Snakemake rule to run MACSE on an input directory with CDS sequences
# Gregg Thomas, January 2023
#############################################################################

# snakemake -p -s 06_macse.smk --configfile 7spec_config.yaml --profile slurm_profile/ --dryrun

#############################################################################

import os

#############################################################################

# 5 loci would not run (program exits successfully but no output files are generated?) this step:
# ENSMUST00000164502
# ENSMUST00000036934
# ENSMUST00000143764
# ENSMUST00000168960
# ENSMUST00000097897

#############################################################################

DATASET = config["dataset"]
BASE_OUTDIR = config["aln_directory"]

PRE_TRIMDIR = os.path.join(BASE_OUTDIR, "01-pre-trim", "nt")
ALNDIR = os.path.join(BASE_OUTDIR, "02-macse")

loci = [];
for f in os.listdir(PRE_TRIMDIR):
    if not os.path.isfile(os.path.join(PRE_TRIMDIR, f)) or not f.endswith("pretrim.fa"):
        continue;

    locus = f.split("-")[0];
    loci.append(locus);

print("# aligning", len(loci), "orthogroups");

#############################################################################
# Final rule - rule that depends on final expected output file and initiates all
# the other rules

localrules: all

rule all:
    input:
        expand(os.path.join(ALNDIR, "nt", "{locus}-7spec-cds-NT.pretrim.macse.fa"), locus=loci),
        expand(os.path.join(ALNDIR, "aa", "{locus}-7spec-cds-AA.pretrim.macse.fa"), locus=loci)
        # Output files from run_macse

#############################################################################
# Pipeline rules

rule run_macse:
    input:
        os.path.join(PRE_TRIMDIR, "{locus}-7spec-cds-NT.pretrim.fa")
    output:
        nt = os.path.join(ALNDIR, "nt", "{locus}-7spec-cds-NT.pretrim.macse.fa"),
        aa = os.path.join(ALNDIR, "aa", "{locus}-7spec-cds-AA.pretrim.macse.fa")
    params:
        locus = "{locus}",
    log:
        os.path.join(ALNDIR, "logs", "{locus}-7spec-cds.macse.log")
    resources:
        cpus = 4,
	    time = "8:00:00"
    shell:
        """
        macse -prog alignSequences -seq {input} -out_NT {output.nt} -out_AA {output.aa} > {log} 2>&1
        """

#############################################################################
