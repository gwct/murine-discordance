#!/usr/bin/python3
############################################################
# For Penn genomes, 12.2020
# Parses concordance analyses
############################################################

import sys, os, core, argparse, treeparse as tp
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

#alndir = "../aln-2/5000kb-10kb-concat/";
#windows_file = "../data-2/10kb-0.5-0.5-topo-counts-tt-hotspots.csv";

treedir = os.path.join("..", "tree-2", args.window_size + "kb-10kb-concat");
if not os.path.isdir(treedir):
    os.system("mkdir " + treedir);
# File names

outfilename = "../data-2/5mb-cf.csv";
if os.path.isfile(outfilename) and not args.overwrite:
    sys.exit( " * Error: Output file (" + outfilename + ") already exists! Explicity specify --overwrite to overwrite it and its contents.");
# Current output file name.

with open(outfilename, "w") as outfile:
    core.runTime("# 5Mb concordance factors from 10kb windows", outfile);
    core.PWS("# Tree directory:       " + treedir, outfile);
    core.PWS("# Output file:          " + outfilename, outfile);
    core.PWS("# ----------------");
    # Input info for log file

    headers = "window,num 10kb trees,ID,gCF,gCF_N,gDF1,gDF1_N,gDF2,gDF2_N,gDFP,gDFP_N,gN,sCF,sCF_N,sDF1,sDF1_N,sDF2,sDF2_N,sN,Label,Length,clade,mouse anc,unparsed tree";
    outfile.write(headers + "\n");
    # Headers for the output file.

    for chrome in chrome_sizes:
        chr_str = "chr" + chrome;
        # Get every window on each chromosome.

        chrome_treedir = os.path.join(treedir, chr_str + "-iqtree");
        chrome_loci_dir = os.path.join(chrome_treedir, "loci");
        chrome_conc_dir = os.path.join(chrome_treedir, "concord");
        # Chrome files.

        for window in os.listdir(chrome_loci_dir):
        # Get info for every 5mb window on each chrome.

            concord_treefile = os.path.join(chrome_conc_dir, window, window + ".cf.branch");
            concord_stats = os.path.join(chrome_conc_dir, window, window + ".cf.stat");
            window_10kb_trees = os.path.join(chrome_treedir, "10kb-trees", window + "-10kb-trees.txt");
            # Window files.

            num_trees = str(len(open(window_10kb_trees, "r").readlines()));
            # Get the number of 10kb trees in the current 5mb window

            concord_tree = open(concord_treefile, "r").read();
            tinfo, tree, root = tp.treeParse(concord_tree);
            # Read and parse the tree with iqtree's branch labels.

            mm10_anc = tinfo['mm10'][1];
            mm10_anc = tinfo[mm10_anc][3];
            # Get the ancestral mouse node for the current tree.

            clades = {};
            for node in tinfo:
                if tinfo[node][2] == 'tip':
                    continue;

                iq_lab = tinfo[node][3];
                cur_clade = tp.getClade(node, tinfo);
                cur_clade = ";".join(sorted(cur_clade));

                clades[iq_lab] = cur_clade;
            # Assign the clade labels for each iqtree branch label.

            first = True;
            for line in open(concord_stats):
            # Read the iqtree stats file.
                if line[0] == "#":
                    continue;
                if first:
                    first = False;
                    continue;
                # Skip comment and header lines.

                line = line.strip().split("\t");
                iq_lab = line[0];
                # Parse the line and get the node label.

                anc_flag = "FALSE";
                if iq_lab == mm10_anc:
                    anc_flag = "TRUE"
                # Flag whether the current label is the mouse ancestor.

                w = window.split("-");
                window_out = w[0] + ":" + w[1] + "-" + w[2];

                outline = [window_out, num_trees] + line + [clades[iq_lab], anc_flag, "\"" + concord_tree.strip() + "\""];
                outfile.write(",".join(outline) + "\n");
                # Combine all info and write the output.



        