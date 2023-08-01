#!/usr/bin/python
############################################################
# For rodent genomes, 04.2023
# Checks aligned windows and makes sure at least 4 sequences
# are unique to make trees from
# Adds this info to the windows file
############################################################

import sys
import os
import yaml
from collections import Counter
import datetime

############################################################

def fastaReadSeqs(filename, header_sep=False):
# Read a FASTA formatted sequence file
# Great iterator and groupby code from: https://www.biostars.org/p/710/ 
# Returns dictionary with the key:value format as title:sequence.

    import gzip
    from itertools import groupby

    #compression = core.detectCompression(filename);
    compression = "none";

    if compression == "gz":
        file_stream = gzip.open(filename);
        fa_iter = (x[1] for x in groupby(file_stream, lambda line: line.decode()[0] == ">"));
        readstr = lambda s : s.decode().strip();
    elif compression == "none":
        file_stream = open(filename); 
        fa_iter = (x[1] for x in groupby(file_stream, lambda line: line[0] == ">"));
        readstr = lambda s : s.strip();
    # Read the lines of the file depending on the compression level
    # file_stream opens the file as an iterable
    # groupby takes an iterable (file_stream) and a function that indicates the key of the group. It iterates over
    # the iterable and when it encounters a key, it groups all following items to it (perfect for FASTA files).
    # fa_iter is a generator object holding all the grouped iterators.
    # readstr is a function that changes depending on compression level -- for compressed files we also need to decode
    # each string in the iterators below.

    seqdict = {};
    # A dictionary of sequences:
    # <sequence id/header> : <sequence>

    for header_obj in fa_iter:
        #print(header_obj)
        header = readstr(header_obj.__next__());
        # The header object is an iterator. This gets the string.

        curkey = header[1:];
        # This removes the ">" character from the header string to act as the key in seqdict

        if header_sep:
            curkey = curkey.split(header_sep)[0];

        seq = "".join(readstr(s) for s in fa_iter.__next__());
        # The current header should correspond to the current iterator in fa_iter. This gets all those
        # lines and combines them as a string.

        #print(header, len(seq));

        seqdict[curkey] = seq;
        # Save the sequence in seqdict

    return seqdict;

######################

def containsOnlyMissing(test_str):
    missing_chars = set("NnXx-");
    if set(test_str) <= missing_chars:
        return True;
    else:
        return False;
# Checks if a sequence contains only missing characters

############################################################

config_file = "rodent-windows-config.yaml";
# Config file for project

with open(config_file, "r") as f:
    config = yaml.safe_load(f)

winsize_kb = config["winsize_kb"];
print("winsize_kb: " + winsize_kb);

repeat_cutoff = config["repeat_cutoff"];
print("repeat_cutoff: " + repeat_cutoff);

missing_cutoff = config["missing_cutoff"];
print("missing_cutoff: " + missing_cutoff);

INPUT_CHROMES = config["chromes"];
# Input params

SEQ_BASE_DIR = config["seq_dir"];
SUB_DIR = winsize_kb + "-" + repeat_cutoff + "-" + missing_cutoff;
TRIM_OUTDIR = config["trim_out_dir"];
# Dirs

WINDOW_INPUT = os.path.join(config["data_dir"], SUB_DIR + "-windows.tsv");
WINDOW_OUTPUT = os.path.join(config["data_dir"], SUB_DIR + "-windows-filter.tsv");
print("window file: " + WINDOW_INPUT);
# Files

######################

outlines = [];
first = True;
num_non_unique, num_w_missing, num_filtered = 0, 0, 0;
# The list of output lines to compile and a flag for the header line

for line in open(WINDOW_INPUT):
    if line[0] == "#":
        outlines.append(line);
        continue;
    # Keep the lines generated from 03_get_window_seqs.py

    if first:
        outlines.append("# 05_aln_check.py for chromes " + ",".join(INPUT_CHROMES) + " on " + datetime.datetime.now().strftime("%m/%d/%Y %H:%M:%S") + "\n");
        # Add a line about the current run

        outlines.append(line);
        first = False;
        continue;
    # Add the header line to the output

    line_list = line.strip().split("\t");
    # Parse the current line

    if line_list[0] in INPUT_CHROMES and all( filters == "PASS" for filters in [line_list[5], line_list[7]] ):
    # Check if the current window is in the input chromes and passed the other filters

        chrome = line_list[0];
        window = line_list[3].replace(":", "-");
        trim_file = os.path.join(TRIM_OUTDIR, chrome, SUB_DIR, window + "-mafft-trimal.fa");
        assert os.path.isfile(trim_file), trim_file;
        # Get the trimmed alignment file name and make sure it exists

        cur_seqs = fastaReadSeqs(trim_file);
        seq_counts = dict(Counter(list(cur_seqs.values())));
        # Read the seqs and count unique seqs

        uniq_filter = False;
        if len(seq_counts) < 4:
            uniq_filter = True;
            num_non_unique += 1;
        # Set the unique filter if the alignment contains fewer than 4 unique sequences

        missing_filter = False;
        for seq in cur_seqs.values():
            if containsOnlyMissing(seq):
                missing_filter = True;
                num_w_missing += 1;
                break;
        # Set the missing filter if the alignment contains 1 or more sequences that are made up of only missing data/gaps
        
        if uniq_filter or missing_filter:
            filter_str = "FILTER";
            num_filtered += 1;
        else:
            filter_str = "PASS";
        # Set the final filter string if either of the two filters are True

        line_list[8] = filter_str;
        outlines.append("\t".join(line_list) + "\n");
        # Change the current line and add to list of output lines
    else:
        outlines.append(line);
    # If the window isn't in the current chromes or was filtered out already, just re-write it as is
## End line loop

print("num filtered: " + str(num_filtered));

######################

print("writing output table: " + WINDOW_OUTPUT);
with open(WINDOW_OUTPUT, "w") as outfile:
    for outline in outlines:
        outfile.write(outline);
# Write the output table

######################        