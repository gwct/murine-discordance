#!/usr/bin/python3
############################################################
# For Penn genomes, 10.2020
# Given the results from 10_topology_test_gen.py, this will
# add columns to the 09_topo_counts.py file about the
# results of each test.
############################################################

import sys, os, core, re, argparse, csv

############################################################
# Globals
chrome_sizes = {"1" : 195471971, "2" : 182113224, "3" : 160039680, "4" : 156508116, "5" : 151834684, "6" : 149736546, "7" : 145441459,
                "8" : 129401213, "9" : 124595110, "10" : 130694993, "11" : 122082543, "12" : 120129022, "13" : 120421639, "14" : 124902244,
                "15" : 104043685, "16" : 98207768, "17" : 94987271, "18" : 90702639, "19" : 61431566, "X" : 171031299 }; #, "Y" : 91744698

# Mouse chromosome sizes from /mnt/beegfs/gt156213e/ref-genomes/mm10/mm10.chrome.sizes"

############################################################
# Options
parser = argparse.ArgumentParser(description="Parse topology test results.");
parser.add_argument("-w", dest="window_size", help="The size of the sliding window in kb.", default=False);
parser.add_argument("-m", dest="marker_window_size", help="The size of the recombination marker windows in Mb.", default=False);
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

if not args.marker_window_size:
    sys.exit(" * Error 1: Marker window size in Mb (-m) must be defined.");
if float(args.marker_window_size) <= 0:
    sys.exit(" * Error 2: Marker window size in Mb (-m) must be positive.");
else:
    marker_wsize_str = str(args.marker_window_size);
# Parse the marker window size

in_prefix = args.window_size + "kb-0.5-0.5";
out_prefix = args.window_size + "kb-0.5-0.5-" + marker_wsize_str + "mb";
tt_dir = "../tree-2/" + in_prefix + "/topology-tests/"

infilename = "../data-3/" + out_prefix + "-topo-counts.csv";
outfilename = "../data-3/" + out_prefix + "-topo-counts-tt.csv";

if os.path.isfile(outfilename) and not args.overwrite:
    sys.exit( " * Error: Output file (" + outfilename + ") already exists! Explicity specify --overwrite to overwrite it.");
# Make sure --overwrite is set if the output file already exists
# File names

with open(outfilename, "w") as outfile, open(infilename) as csvfile:
    core.PWS("# Topology test dir:    " + tt_dir);
    core.PWS("# ----------------");

    au_pass, au_fail, au_none, au_one, au_skip = 0,0,0,0,0;
    # Counts for the topology tests.

    first = True;
    windows_reader = csv.reader(csvfile);
    for line in windows_reader:
        # print(line);
        if line[0][0] == "#":
            continue;
        # Skip the log lines from 09_topo_counter.py
        if first:
            headers = line + ["AU test","Alternate topos"];
            outfile.write(",".join(headers) + "\n");
            first = False;
            continue;
        # Get the headers and add the additional columns.

        au_result, alt_topos = "NA", "NA";
        # Test result variables.

        repeat_filter, missing_filter = line[6], line[7];
        if "FILTER" not in [repeat_filter, missing_filter]:
        # Skip windows that have been filtered out.
            
            window, chrome = line[0].replace(":", "-"), line[1];
            window_tt_dir = os.path.join(tt_dir, chrome, window);
            topology_test_file = os.path.join(window_tt_dir, "topo-test.iqtree");
            tree_file = os.path.join(window_tt_dir, window + "-topos.tre");
            # Get the tree and results file for this current window.
            # print(window);
            trees = [ tree.strip() for tree in open(tree_file, "r").readlines() if tree != "\n" ];
            # Read the input trees.

            # print(trees);

            if not os.path.isfile(topology_test_file):
                print("# CANNOT FIND TOPOLOGY TEST FILE, SETTING TEST TO FAIL: " + topology_test_file);
                au_result = "FAIL";
                au_skip += 1;
            # Some windows may have errored out. Skip them.

            if au_result != "FAIL":
                record = False;
                au_pvals = [];
                for iqtree_line in open(topology_test_file):
                    if "p-AU" in iqtree_line:
                        record = True;
                        continue;

                    if record:
                        if iqtree_line[0] == "-":
                            continue;
                        if iqtree_line == "\n":
                            record = False;
                            break;

                        last_char = iqtree_line.strip()[-1];
                        au_pvals.append(last_char);
                        # Gets the last charcater from each test line (+ or -)

                num_pass = au_pvals.count("+");
                # Counts number of positive test results.

                if num_pass == 1 and au_pvals[0] == "+":
                    au_result = "PASS";
                    au_pass += 1;
                # The window passes the AU test only with the inferred tree.

                if num_pass == 1 and au_pvals[0] == "-":
                    au_result = "FAIL-ORIG";
                    au_one += 1;
                    alt_topos = [];
                    for p in range(len(au_pvals)):
                        if p == 0:
                            continue;
                        if au_pvals[p] == "+":
                            alt_topos.append(trees[p]);
                    alt_topos = "\"" + "/".join(alt_topos) + "\"";
                # The window passes the AU test only with a different tree than the one inferred.

                if num_pass == 0:
                    au_result = "FAIL-NONE";
                    au_none += 1;
                # The window passes the AU test with none of the trees

                if num_pass > 1:
                    au_result = "FAIL-MULT";
                    au_fail += 1;
                    alt_topos = [];
                    for p in range(len(au_pvals)):
                        if p == 0:
                            continue;
                        if au_pvals[p] == "+":
                            alt_topos.append(trees[p]);
                    alt_topos = "\"" + "/".join(alt_topos) + "\"";
                # The window passes the AU test with more than one tree

            # print(au_pvals);
            # print(au_result);
            # print(alt_topos);
            # if window == "chr1-3220001-3230000":
            #     sys.exit();

        line[17] = "\"" + line[17] + "\"";
        line[20] = "\"" + line[20] + "\"";

        outline = line + [au_result, alt_topos];
        outfile.write(",".join(outline) + "\n");
        # Write the output.

    core.PWS("# ----------------");
    core.PWS("# Windows that pass the AU test with the inferred tree:        " + str(au_pass));
    core.PWS("# Windows that pass the AU test with a different tree:         " + str(au_one));
    core.PWS("# Windows that pass the AU test with no trees:                 " + str(au_none));
    core.PWS("# Windows that pass the AU test with more than one tree:       " + str(au_fail));
    core.PWS("# Windows that didn't have an AU test file:                    " + str(au_skip));
    core.PWS("# ----------------");
    ###################
    



