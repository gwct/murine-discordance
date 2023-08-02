#!/usr/bin/python3
############################################################
# For Penn genomes, 12.2020
# Runs concordance factor analysis on windows.
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

alndir = "../aln-2/5000kb-10kb-concat/";
#windows_file = "../data-2/10kb-0.5-0.5-topo-counts-tt-hotspots.csv";

treedir = os.path.join("..", "tree-2", args.window_size + "kb-10kb-concat");
if not os.path.isdir(treedir):
    os.system("mkdir " + treedir);
# File names

for chrome in chrome_sizes:
    chr_str = "chr" + chrome;
    print("# " + core.getDateTime() + " Concordance factors for chrome: " + chr_str);
    # Run concordance factors for every chromosome.

    chrome_treedir = os.path.join(treedir, chr_str + "-iqtree");
    chrome_loci_dir = os.path.join(chrome_treedir, "loci");
    chrome_conc_dir = os.path.join(chrome_treedir, "concord");
    # Chrome files.

    for window in os.listdir(chrome_loci_dir):
    # Run concordance factors for every 5mb window on the current chromosome.

        window_tree = os.path.join(chrome_loci_dir, window, window + "-rooted.treefile");
        window_10kb_trees = os.path.join(chrome_treedir, "10kb-trees", window + "-10kb-trees.txt");
        window_aln = os.path.join(alndir, chr_str, window + ".fa");
        window_concord_dir = os.path.join(chrome_conc_dir, window);
        if not os.path.isdir(window_concord_dir):
            os.system("mkdir " + window_concord_dir);
        # Window files.

        prefix = os.path.join(window_concord_dir, window);
        # Window output prefix.

        iqtree_cmd = "iqtree -t " + window_tree + " --gcf " + window_10kb_trees + " -s " + window_aln + " --scf 100 --prefix " + prefix + " -T 1";
        os.system(iqtree_cmd);
        #print(iqtree_cmd);
        #sys.exit();
        # Generate and run the iqtree command for concordance factors.