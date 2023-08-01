#############################################################################
# Pipeline for read mapping simulations with varying divergence
#############################################################################

import os
import re
from itertools import product

#############################################################################
# Example cmd

# snakemake -p -s 04_alignment.smk --configfile rodent-windows-config.yaml --profile profiles/slurm_profile/ --keep-going --dryrun

#############################################################################
# Functions

#############################################################################
# Input and output info

work_dir = os.path.dirname(__file__);
print("working dir: " + work_dir);

winsize_kb = config["winsize_kb"];
print("winsize_kb: " + winsize_kb);

repeat_cutoff = config["repeat_cutoff"];
print("repeat_cutoff: " + repeat_cutoff);

missing_cutoff = config["missing_cutoff"];
print("missing_cutoff: " + missing_cutoff);

INPUT_CHROMES = config["chromes"];

print("reading window ids...");

SEQ_BASE_DIR = config["seq_dir"];
SUB_DIR = winsize_kb + "-" + repeat_cutoff + "-" + missing_cutoff;
# Dirs

######################

WINDOW_INPUT = os.path.join(config["data_dir"], SUB_DIR + "-windows.tsv");
print("window file: " + WINDOW_INPUT);

CHROMES, WINDOWS, first = [], [], True;
for line in open(WINDOW_INPUT):
    if line[0] == "#":
        continue;
    if first:
        first = False;
        continue;
    # Skip header lines in the window file

    line = line.strip().split("\t");
    if line[0] in INPUT_CHROMES and all( filters == "PASS" for filters in [line[5], line[7]] ):
        CHROMES.append(line[0]);
        WINDOWS.append(line[3].replace(":", "-"));
## Read the windows to align given the input chromes and the filters

print(len(WINDOWS), "window ids read");

######################

ALN_OUTDIR = config["aln_out_dir"];
TRIM_OUTDIR = config["trim_out_dir"];
# Output dirs

#############################################################################
# Final rule - rule that depends on final expected output file and initiates all
# the other rules

localrules: all

rule all:
    input:
        expand(os.path.join(TRIM_OUTDIR, "{chrome}", SUB_DIR, "{win}-mafft-trimal.fa"), zip, chrome=CHROMES, win=WINDOWS),

#############################################################################
# Pipeline rules

rule mafft:
    input:
        seq_file = os.path.join(SEQ_BASE_DIR, "{chrome}", SUB_DIR, "{seq}.fa")
    output:
        aln_file = os.path.join(ALN_OUTDIR, "{chrome}", SUB_DIR, "{seq}-mafft.fa")
    log:
        os.path.join(ALN_OUTDIR, "{chrome}", SUB_DIR, "logs", "{seq}-mafft.log")
    shell:
        """
        mafft --adjustdirection --preservecase {input.seq_file} 2> {log} 1> {output.aln_file}
        """
# Run each locus through mafft

######################

rule trimal:
    input:
        aln_file = os.path.join(ALN_OUTDIR, "{chrome}", SUB_DIR, "{seq}-mafft.fa")
    output:
        trim_file = os.path.join(TRIM_OUTDIR, "{chrome}", SUB_DIR, "{seq}-mafft-trimal.fa")
    log:
        os.path.join(TRIM_OUTDIR, "{chrome}", SUB_DIR, "logs", "{seq}-mafft-trimal.log")
    shell:
        """
        trimal -in {input.aln_file} -out {output.trim_file} -gt 0.5 &> {log}
        """
# Run each locus through trimal

######################
