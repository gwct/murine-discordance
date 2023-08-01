#!/usr/bin/python
############################################################
# For Penn genomes, 05.2020
# Counts topologies from window-based trees
############################################################

import sys, os, re, argparse
from collections import defaultdict
sys.path.append("/home/gt156213e/bin/core/python/lib/");
import core
import treeparse as tp

############################################################
# Functions
def treeDictRmBlength(tdict):
    for node in tdict:
        tdict[node][0] = 'NA';
    return tdict;

############################################################
# Globals
chrome_sizes = {"1" : 195471971, "2" : 182113224, "3" : 160039680, "4" : 156508116, "5" : 151834684, "6" : 149736546, "7" : 145441459,
                "8" : 129401213, "9" : 124595110, "10" : 130694993, "11" : 122082543, "12" : 120129022, "13" : 120421639, "14" : 124902244,
                "15" : 104043685, "16" : 98207768, "17" : 94987271, "18" : 90702639, "19" : 61431566, "X" : 171031299 }; #, "Y" : 91744698

#chrome_sizes = { "19" : 61431566 }; #, "Y" : 91744698
# Mouse chromosome sizes from /mnt/beegfs/gt156213e/ref-genomes/mm10/mm10.chrome.sizes"

############################################################

parser = argparse.ArgumentParser(description="Rodent window topology counter");
parser.add_argument("-w", dest="window_size", help="The size of the sliding window in kb.", default=False);
parser.add_argument("-m", dest="marker_window_size", help="The size of the recombination marker windows in Mb.", default=False);
parser.add_argument("--overwrite", dest="overwrite", help="If the output file already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
args = parser.parse_args();
# Options

if not args.window_size:
    sys.exit(" * Error 1: Window size in kb (-w) must be defined.");
if int(args.window_size) < 1:
    sys.exit(" * Error 2: Window size in kb (-w) must be a positive integer.");
else:
    wsize = float(args.window_size) * 1000;
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
windowfile = "../data/" + in_prefix + "-windows-filter.tsv";

if float(args.marker_window_size) >= 1:
    markerfile = "../data/recombination-markers/cox-markers/Revised_HSmap_SNPs-mm10-" + marker_wsize_str + "mb-rates.csv";
else:
    marker_wsize_str = str(int(float(args.marker_window_size) * 1000))
    markerfile = "../data/recombination-markers/brunschwig-rates/mm10-markers-" + marker_wsize_str + "kb.csv";

outfilename = "../data/" + out_prefix + "-topo-counts.csv";
if os.path.isfile(outfilename) and not args.overwrite:
    sys.exit( " * Error: Output file (" + outfilename + ") already exists! Explicity specify --overwrite to overwrite it.");
# Make sure --overwrite is set if the output file already exists
# File names

headers = "window,chr,chr len,start,end,perc repeat,repeat filter, missing filter,gcf,scf,topo num overall,topo count overall,topo rank overall,topo color,topo num chrome,topo count chrome,topo rank chrome,topo,concat chrome topo,astral chrome topo,unparsed tree,marker window,num markers,marker slope,marker uniq trees";
# Column headers for the output file.

with open(outfilename, "w") as outfile:
    core.runTime("# Window-based topology counter", outfile);
    core.PWS("# Window size:          " + wsize_str, outfile);
    core.PWS("# Window file:          " + windowfile, outfile);
    core.PWS("# Marker size:          " + marker_wsize_str, outfile);
    core.PWS("# Marker file:          " + markerfile, outfile);
    core.PWS("# Output file:          " + outfilename, outfile);
    core.PWS("# ----------------");
    # Input info for log file

    core.PWS("# " + core.getDateTime() + " Reading marker windows...", outfile);
    marker_windows = defaultdict(dict);
    first = True;
    for line in open(markerfile):
        if line[0] == "#":
            continue;
        if first:
            first = False;
            continue;
        line = line.strip().split(",");
        #print(line);
        if float(args.marker_window_size) >= 1:
            marker_chrome = line[6];
            if not marker_chrome.startswith("chr"):
                marker_chrome = "chr" + marker_chrome;

            marker_windows[line[8]] = { 'chr' : marker_chrome, 'start' : int(line[9]), 'end' : int(line[10]), 'num-markers' : line[11], 'marker-slope' : line[12] };
        else:
            marker_chrome = line[1];
            marker_windows[line[4]] = { 'chr' : marker_chrome, 'start' : int(line[5]), 'end' : int(line[6]), 'num-markers' : line[7], 'marker-slope' : line[3] };

    core.PWS("# Marker windows read: " + str(len(marker_windows)), outfile);
    core.PWS("# ----------------");
    # Read the marker windows.

    core.PWS("# " + core.getDateTime() + " Reading windows...", outfile);
    windows, marker_win_trees, first = {}, {}, True;
    for line in open(windowfile):
        if line[0] == "#":
            continue;
        if first:
            first = False;
            continue;
        line = line.strip().split("\t");
        #print(line);

        chrome, start, end, window, perc_repeat, repeat_filter, missing_filter, aln_filter = line[0], int(line[1]), int(line[2]), line[3], line[4], line[5], line[7], line[8];

        chr_len = str(chrome_sizes[chrome.replace("chr", "")]);
        # print(line);
        for m in marker_windows:
            # print(chrome, marker_windows[m]['chr'])
            # print(start, marker_windows[m]['start'])
            # print(end, marker_windows[m]['end'])
            #print(m);
            #print(marker_windows[m]);
            # print(chrome);
            # sys.exit();
            if marker_windows[m]['chr'] == chrome and marker_windows[m]['start']-1 <= start and marker_windows[m]['end'] >= end:
                marker_window, num_markers, marker_slope = m, marker_windows[m]['num-markers'], marker_windows[m]['marker-slope'];
        # Find the marker window that includes this tree window.

        windows[window] = { 'chr' : chrome, 'chr-len' : chr_len, 'start' : str(start), 'end' : str(end), 
                             'perc-repeat' : perc_repeat, 'repeat-filter' : repeat_filter, 'missing-filter' : missing_filter,
                                'gcf' : "NA", 'scf' : "NA",
                                'topo-num-all' : "NA", 'topo-count-all' : "NA", 'topo-rank-all' : "NA", 
                                'topo-num-chr' : "NA", 'topo-count-chr' : "NA", 'topo-rank-chr' : "NA", 
                                'topo' : "NA", 'concat-topo-match' : "NA", 'astral-topo-match' : "NA", 'unparsed-tree' : "NA", 'clade-set' : "NA",
                                'marker-win' : marker_window, 'num-markers' : num_markers, 'marker-slope' : marker_slope };
        marker_win_trees[marker_window] = [];
        # Initialize info for this tree window.
    core.PWS("# Windows read:         " + str(len(windows)), outfile);
    core.PWS("# ----------------");
    # Read the tree windows.

    topo_colors, num_colors, colors = {}, 0, ["#db6d00", "#004949", "#006ddb", "#920000", "#490092", "#6cb6ff", "#24ff24", "#fdb4da", "#ffff6d", "#009292"];
    # Variables for assigning colors to each topology.

    uniq_topos_all, topo_counts_all = [], defaultdict(int);
    # Variables for overall counting.

    ##########

    for chromosome in chrome_sizes:
    # Count trees for each chromosome.

        chrstr = "chr" + chromosome
        chr_end = chrome_sizes[chromosome];
        chr_end_str = str(chr_end);
        # Chromosome variables.

        indir = "../tree/" + chrstr + "/" + in_prefix + "/";
        # Tree directory for this chromosomes.
        
        locidir = os.path.join(indir, "loci");
        assert os.path.isdir(locidir), "\nCANNOT FIND TREE DIRECTORY: " + locidir;

        concattreefile = os.path.join(indir, "concat", chrstr + "-concat.treefile.rooted");
        assert os.path.isfile(concattreefile), "\nCANNOT FIND CONCATENATED TREE: " + concattreefile;

        concordfile = os.path.join(indir, "astral", chrstr + "-astral.cf.stat");
        assert os.path.isfile(concordfile), "\nCANNOT FIND CONCORDANCE FILE: " + concordfile;

        astraltreefile = os.path.join(indir, "astral", chrstr + "-astral.treefile.rooted");
        assert os.path.isfile(concattreefile), "\nCANNOT FIND ASTRAL TREE: " + astraltreefile;
        # Tree files for this chromosome.

        core.PWS("# Chromosome:           " + chromosome, outfile);
        core.PWS("# Chromosome length:    " + chr_end_str, outfile);
        core.PWS("# Locus tree directory: " + locidir, outfile);
        core.PWS("# Concat tree file:     " + concattreefile, outfile);
        core.PWS("# ASTRAL tree file:     " + astraltreefile, outfile);
        core.PWS("# Concordance file:     " + concordfile, outfile);
        core.PWS("# ----------------", outfile);
        # Log info for this chromosome.

        core.PWS("# " + core.getDateTime() + " Reading concatenated chromosome tree...", outfile);
        concattree = open(concattreefile, "r").read().strip();
        core.PWS("# Concat chromosome tree:      \"" + concattree + "\"", outfile);
        tinfo, concat_tree, root = tp.treeParse(concattree);
        concat_tree_topo = set([frozenset(tp.getClade(node, tinfo)) for node in tinfo if tinfo[node][2] != 'tip']);
        core.PWS("# Concat chromosome topology:  \"" + concat_tree + "\"", outfile);
        # Read the concatenated tree for this chromosome.

        core.PWS("# " + core.getDateTime() + " Reading ASTRAL chromosome tree...", outfile);
        astraltree = open(astraltreefile, "r").read().strip();
        core.PWS("# ASTRAL chromosome tree:      \"" + astraltree + "\"", outfile);
        tinfo, astral_tree, root = tp.treeParse(astraltree);
        astral_tree_topo = set([frozenset(tp.getClade(node, tinfo)) for node in tinfo if tinfo[node][2] != 'tip']);
        core.PWS("# ASTRAL chromosome topology:  \"" + astral_tree + "\"", outfile);
        # Read the ASTRAL tree for this chromosome.

        astral_concat_match = "FALSE";
        if concat_tree_topo == astral_tree_topo:
            astral_concat_match = "TRUE";
        core.PWS("# ASTRAL-concat match:   " + astral_concat_match, outfile);
        core.PWS("# ----------------");
        # Check if the conactenated and ASTRAL topologies are the same for this chromosome.

        core.PWS("# " + core.getDateTime() + " Calculating average CFs...", outfile);
        gcfs, scfs, first = [], [], True;
        for line in open(concordfile):
            if line[0] == "#":
                continue;
            if first:
                first = False;
                continue;
            line = line.strip().split("\t");
            if line[1] != "NA":
                gcfs.append(float(line[1]));
            if line[10] != "NA":
                scfs.append(float(line[10]));

        avg_gcf = str(core.mean(gcfs))
        avg_scf = str(core.mean(scfs));
        core.PWS("# Average gCF:          " + avg_gcf, outfile);
        core.PWS("# Average sCF:          " + avg_scf, outfile);
        core.PWS("# ----------------");
        # Using the IQ-tree stat file to get average CFs for this chromosome.

        core.PWS("# " + core.getDateTime() + " Reading windows and counting topologies...", outfile);
        uniq_topos, topo_counts = [], defaultdict(int);
        # Lists of unique clade sets for each window and the unique lists of clade sets to count topologies.

        num_trees, chrome_to_all = 0, {};
        # The number of trees read and the conversion between chromosome topology number and overall topology number for the coloring.

        filtered_windows, no_tree_file = 0,0;

        for window in windows:
            if windows[window]['chr'] != chrstr:
                continue;
            # Get the tree for each window on this chromosome.

            # if window != "chr1:85280001-85290000":
            #     continue;
            # print(window);

            windows[window]['gcf'] = avg_gcf;
            windows[window]['scf'] = avg_scf;
            # Add concordance factors.

            if windows[window]['repeat-filter'] == "PASS" and windows[window]['missing-filter'] == "PASS":
            # If the window was previously filtered then skip it.

                window_str = window.replace(":", "-");
                window_name = "chr" + chromosome + "-" + windows[window]['start'] + "-" + windows[window]['end'];
                treefile = os.path.join(locidir, window_name, window_name + "-mafft-trimal.treefile.rooted");
                # Get window name and tree file.
                
                if not os.path.isfile(treefile):
                    print("# CANNOT FIND TREE FILE, SETTING STATUS TO FILTER: " + treefile);
                    no_tree_file += 1;
                    windows[window]['missing-filter'] = "FILTER";
                    continue;
                # If the file isn't found because one of the alignment or tree steps failed, just set it to filter here.

                tree = open(treefile, "r").read().strip();
                # Read the tree from the file.

                tinfo, tree_parsed, root = tp.treeParse(tree);
                #tree_parsed = re.sub("<[\d]+>", "", tree_parsed);
                # Parse the tree.

                clade_topo = set([frozenset(tp.getClade(node, tinfo)) for node in tinfo if tinfo[node][2] != 'tip']);
                # Get the topology as the set of all clades present in the tree.

                if clade_topo not in uniq_topos_all:
                    uniq_topos_all.append(clade_topo);
                # If this topology hasn't been seen before at all, add it to the list of unique topologies.

                if clade_topo not in uniq_topos:
                    uniq_topos.append(clade_topo);
                # If this topology hasn't been seen before on this chromosome, add it to the list of unique topologies.

                marker_window = windows[window]['marker-win'];
                if clade_topo not in marker_win_trees[marker_window]:
                    marker_win_trees[marker_window].append(clade_topo);
                # If this topology hasn't been seen before in this marker window, add it to the list of unique topologies.
                
                concat_match = "FALSE";
                if clade_topo == concat_tree_topo:
                    concat_match = "TRUE";
                astral_match = "FALSE";
                if clade_topo == astral_tree_topo:
                    astral_match = "TRUE";
                # Check if the current topology matches the chromosome level topologies.

                topo_num_all = uniq_topos_all.index(clade_topo)+1;
                topo_counts_all[topo_num_all] += 1;
                topo_num = uniq_topos.index(clade_topo)+1;
                topo_counts[topo_num] += 1;
                chrome_to_all[topo_num] = topo_num_all;
                # The topology ID number is then the index in the unique topology list.

                # print(uniq_topos_all);
                # print(topo_num_all);
                # print(topo_counts_all);

                windows[window]['topo-num-all'] = topo_num_all;
                windows[window]['topo-num-chr'] = topo_num;
                windows[window]['topo'] = tree_parsed + ";";
                windows[window]['concat-topo-match'] = concat_match;
                windows[window]['astral-topo-match'] = astral_match;
                windows[window]['unparsed-tree'] = tree;
                windows[window]['clade-set'] = clade_topo;
                # Save output info for this window.

                # print(windows[window]);
                # sys.exit();

                num_trees += 1;
                # Increment the number of trees read.
            # If this window has a tree.
            else:
                filtered_windows += 1;
        core.PWS("# Total trees read:   " + str(num_trees), outfile);
        core.PWS("# No tree file:       " + str(no_tree_file), outfile);
        core.PWS("# Filtered windows:   " + str(filtered_windows), outfile);
        core.PWS("# Unique topologies:  " + str(len(uniq_topos)), outfile);
        core.PWS("# ----------------", outfile);

        core.PWS("# " + core.getDateTime() + " Counting and ranking topologies for this chromosome...", outfile);
        topo_counts_sorted = sorted(topo_counts.items(), key=lambda x: x[1], reverse=True);
        # Sort the counts
        topo_ranks, j = {}, 1;
        for i in topo_counts_sorted:
            topo_ranks[i[0]] = j;
            j += 1;
        # Rank the sorted counts

        # print(topo_counts);
        # print(topo_counts_sorted);
        # print(topo_ranks);
        # print(windows[window]['chr'], chromosome);

        for topo in topo_ranks:
            if topo_ranks[topo] <= 3:
                topo_num_all = chrome_to_all[topo];
                if topo_num_all not in topo_colors:
                    topo_colors[topo_num_all] = colors[num_colors];
                    num_colors += 1;
        # Assign a color to the top 3 topologies if they haven't been assigned already.

        for window in windows:
            if windows[window]['chr'] != chrstr:
                continue;
            if windows[window]['repeat-filter'] == "PASS" and windows[window]['missing-filter'] == "PASS":
                chrome_topo_count = topo_counts[windows[window]['topo-num-chr']];
                chrome_topo_rank = topo_ranks[windows[window]['topo-num-chr']];
                #print(window, chrome_topo_count, chrome_topo_rank);
                windows[window]['topo-count-chr'] = str(chrome_topo_count);
                windows[window]['topo-rank-chr'] = str(chrome_topo_rank);
            # If this window had a tree made, count the number of times that topology occurred.
        core.PWS("# ----------------", outfile);
    # END CHROMSOME LOOP
    ####################

    core.PWS("# " + core.getDateTime() + " Counting and ranking overall topologies...", outfile);
    topo_counts_sorted = sorted(topo_counts_all.items(), key=lambda x: x[1], reverse=True);
    # Sort the counts
    topo_ranks, j = {}, 1;
    for i in topo_counts_sorted:
        topo_ranks[i[0]] = j;
        j += 1;
    # Rank the sorted counts

    for window in windows:
        if windows[window]['repeat-filter'] == "PASS" and windows[window]['missing-filter'] == "PASS":
            #print(window, windows[window]);
            #print(topo_counts_all);
            topo_count = topo_counts_all[windows[window]['topo-num-all']];
            topo_rank = topo_ranks[windows[window]['topo-num-all']];
            windows[window]['topo-count-all'] = str(topo_count);
            windows[window]['topo-rank-all'] = str(topo_rank);
    # If this window had a tree made, count the number of times that topology occurred.
    core.PWS("# ----------------", outfile);         

    print(topo_colors)
    for topo in topo_ranks:
        if topo not in topo_colors:
            topo_colors[topo] = "#999999";
    # Fill in the remaining colors with grey

    core.PWS("# " + core.getDateTime() + " Writing output: " + outfilename, outfile);
    outfile.write(headers + "\n");

    for window in windows:
        w = windows[window];
        marker_window = w['marker-win'];
        num_marker_trees = str(len(marker_win_trees[marker_window]));
        cur_color = "NA";
        if w['topo-num-all'] != "NA":
            cur_color = topo_colors[w['topo-num-all']];
        outline = [ window, w['chr'], w['chr-len'], w['start'], w['end'], w['perc-repeat'], w['repeat-filter'], w['missing-filter'], w['gcf'], w['scf'], str(w['topo-num-all']), w['topo-count-all'], w['topo-rank-all'], '"' + cur_color + '"', str(w['topo-num-chr']), w['topo-count-chr'], w['topo-rank-chr'], '"' + w['topo'] + '"', w['concat-topo-match'], w['astral-topo-match'], '"' + w['unparsed-tree'] + '"', w['marker-win'], w['num-markers'], w['marker-slope'], num_marker_trees ];
        #headers = "window,chr,chr len,start,end,perc repeat,repeat filter, missing filter,gcf,scf,topo num overall,topo count overall,topo rank overall,topo color,topo num chrome,topo count chrome,topo rank chrome,topo,concat chrome topo,astral chrome topo,unparsed tree,marker window,num markers,marker slope,marker uniq trees";
        outfile.write(",".join(outline) + "\n");
        # Write the current output (minus the clade topology if a tree is present for window);
    core.PWS("# ----------------");
    print("# " + core.getDateTime() + " Done!");
