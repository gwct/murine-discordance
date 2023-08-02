#!/usr/bin/python3
############################################################
# For Penn genomes, 12.2020
# This script concatenates all the 10kb alignments within
# a larger window size (i.e. 5MB).
############################################################

import sys, os, copy, core, argparse, csv
from collections import defaultdict

############################################################
# Globals
chrome_sizes = {"1" : 195471971, "2" : 182113224, "3" : 160039680, "4" : 156508116, "5" : 151834684, "6" : 149736546, "7" : 145441459,
                "8" : 129401213, "9" : 124595110, "10" : 130694993, "11" : 122082543, "12" : 120129022, "13" : 120421639, "14" : 124902244,
                "15" : 104043685, "16" : 98207768, "17" : 94987271, "18" : 90702639, "19" : 61431566, "X" : 171031299 }; #, "Y" : 91744698

# Mouse chromosome sizes from /mnt/beegfs/gt156213e/ref-genomes/mm10/mm10.chrome.sizes"

############################################################

parser = argparse.ArgumentParser(description="Concatenate 10kb aligns.");
parser.add_argument("-w", dest="window_size", help="The size of the windows which will contain concatenated 10kb windows in kb.", default=False);
parser.add_argument("--overwrite", dest="overwrite", help="If the output directory already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
args = parser.parse_args();
# IO options

if not args.window_size:
    sys.exit(" * Error 1: Window size in kb (-w) must be defined.");
if int(args.window_size) < 1:
    sys.exit(" * Error 2: Window size in kb (-w) must be a positive integer.");
else:
    wsize = int(args.window_size) * 1000;
    wsize_str = str(wsize);
# Parse the window size

alndir = "../aln-2/10kb-0.5-0.5/";
windows_file = "../data-2/10kb-0.5-0.5-topo-counts-tt-hotspots.csv";

outdir = os.path.join("..", "aln-2", args.window_size + "kb-10kb-concat");
if os.path.isdir(outdir) and not args.overwrite:
    sys.exit( " * Error: Output directory (" + outdir + ") already exists! Explicity specify --overwrite to overwrite it and its contents.");
# Make sure --overwrite is set if the output file already exists
if not os.path.isdir(outdir):
    os.system("mkdir " + outdir);

treedir = os.path.join("..", "tree-2", args.window_size + "kb-10kb-concat");
if not os.path.isdir(treedir):
    os.system("mkdir " + treedir);
# File names

windows = { "chr" + chrome : defaultdict(dict) for chrome in  chrome_sizes };
# Initialize main dict

print("# " + core.getDateTime() + " Reading windows file: " + windows_file);
with open(windows_file) as csvfile:
    first = True;
    num_windows_read = 0;
    windows_reader = csv.reader(csvfile);
    for line in windows_reader:
        # print(line);
        if line[0][0] == "#":
            continue;
        # Skip the log lines from 09_topo_counter.py
        if first:
            headers = line;
            first = False;
            continue;
        # Get the headers and add the additional columns.

        if "FILTER" in [line[6], line[7]]:
            continue;
        # Skip the tree if it was filtered out based on repeats or missing data.

        chrome = line[1];
        window_5mb = line[21];
        window_10kb = line[0];
        tree_10kb = line[20];
        # Get the info for this window.

        windows[chrome][window_5mb][window_10kb] = tree_10kb;
        num_windows_read += 1;
        # Add the info from this 10kb window to the 5mb window entry.
print("# " + core.getDateTime() + " Windows read: " + str(num_windows_read));
print("# ----------------");

windows_concat = 0;
windows_skipped = 0;
files_written = 0;
# Some counting vars.

for chrome in windows:
    print("# " + core.getDateTime() + " Concatenating seqs for chrome: " + chrome);
    # Concatenate alingments for every chrome.

    chrome_alndir = os.path.join(alndir, chrome + "-trimal");
    chrome_outdir = os.path.join(outdir, chrome);
    if not os.path.isdir(chrome_outdir):
        os.system("mkdir " + chrome_outdir);
    chrome_treedir = os.path.join(treedir, chrome + "-iqtree", "10kb-trees");
    if not os.path.isdir(chrome_treedir):
        os.system("mkdir " + chrome_treedir);
    # Chrome files.

    for window_5mb in windows[chrome]:
        first_win = True;
        cur_seqs_5mb = {};
        # Get sequences for every 5Mb window on this chromosome.

        for window_10kb in windows[chrome][window_5mb]:
            infilename = os.path.join(chrome_alndir, window_10kb.replace(":", "-") + "-mafft-trimal.fa");
            # Go through each 10kb sequence for the current 5mb window.

            if not os.path.isfile(infilename):
                windows_skipped += 1;
                continue;
            # Check to make sure the alignment file exists for the current 10kb window.

            cur_seqs_10kb = core.fastaGetDict(infilename);
            # Read the sequences.

            if first_win:
                cur_seqs_5mb = copy.deepcopy(cur_seqs_10kb);
                first_win = False;
            else:
                for title in cur_seqs_5mb:
                    cur_seqs_5mb[title] += cur_seqs_10kb[title];
            windows_concat += 1;
            # If this is the first 10kb window of the current 5mb window, initialize by copying the seq dict.
            # Other wise add each sequence on.

        outfilename = os.path.join(chrome_outdir, window_5mb.replace(":", "-") + ".fa");
        with open(outfilename, "w") as outfile:
            for title in cur_seqs_5mb:
                outfile.write(title + "\n");
                outfile.write(cur_seqs_5mb[title] + "\n");
        # Write the concatenated sequence for the current 5mb window.

        treefilename = os.path.join(chrome_treedir, window_5mb.replace(":", "-") + "-10kb-trees.txt");
        with open(treefilename, "w") as treefile:
            for window_10kb in windows[chrome][window_5mb]:
                treefile.write(windows[chrome][window_5mb][window_10kb] + "\n");
        # Write the 10kb trees to a file for the current 5mb window for gCF.
        files_written += 1;
print("# " + core.getDateTime() + " Windows concatenated: " + str(windows_concat));
print("# " + core.getDateTime() + " Windows skipped:      " + str(windows_skipped));
print("# " + core.getDateTime() + " Files written:        " + str(files_written));
print("# ----------------");