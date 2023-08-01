#############################################################################
# Snakemake rule to run MACSE on an input directory with CDS sequences
# Gregg Thomas, January 2023
#############################################################################

# snakemake -p -s 07_trim_aln.smk --configfile 7spec_config.yaml --profile slurm_profile/ --dryrun

#############################################################################

import os

#############################################################################

DATASET = config["dataset"]
BASE_OUTDIR = config["aln_directory"]

ALNDIR = os.path.join(BASE_OUTDIR, "02-macse", "nt")
TRIMDIR = os.path.join(BASE_OUTDIR, "03-trim")

loci = [];
for f in os.listdir(ALNDIR):
    if not os.path.isfile(os.path.join(ALNDIR, f)) or not f.endswith("pretrim.macse.fa"):
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
        expand(os.path.join(TRIMDIR, "nt", "{locus}-7spec-cds-NT.pretrim.macse.trim.fa"), locus=loci),
        expand(os.path.join(TRIMDIR, "info", "{locus}-7spec-cds-NT.pretrim.macse.trim.csv"), locus=loci)
        # Output files from run_macse

#############################################################################
# Pipeline rules

rule run_macse_trim:
    input:
        os.path.join(ALNDIR, "{locus}-7spec-cds-NT.pretrim.macse.fa")
    output:
        nt = os.path.join(TRIMDIR, "nt", "{locus}-7spec-cds-NT.pretrim.macse.trim.fa"),
        info = os.path.join(TRIMDIR, "info", "{locus}-7spec-cds-NT.pretrim.macse.trim.csv")
    params:
        locus = "{locus}",
    log:
        os.path.join(TRIMDIR, "logs", "{locus}-7spec-cds.trim.log")
    resources:
        cpus = 4,
	    time = "8:00:00"
    shell:
        """
        macse -prog trimAlignment -align {input} -min_percent_NT_at_ends 0.5 -out_NT {output.nt} -out_trim_info {output.info} -respect_first_RF_ON > {log} 2>&1
        """

#############################################################################
