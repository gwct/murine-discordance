#!/usr/bin/python3
###########################################################
# For Penn genomes, 11.2020
# Gets the trees and windows for each aligned transcript.
###########################################################

# grep -v "GL\|JH\|MT" mm10-cds-coords.bed | sed -e 's/^/chr/' > tmp-chr.bed
# bedtools intersect -a tmp-chr.bed -b mm10-10kb-windows.bed -wao | sort > mm10-cds-window-overlap-sorted.bed

import sys, os, core, treeparse as tp, csv
from collections import defaultdict

############################################################

chrome_sizes = {"1" : 195471971, "2" : 182113224, "3" : 160039680, "4" : 156508116, "5" : 151834684, "6" : 149736546, "7" : 145441459,
                "8" : 129401213, "9" : 124595110, "10" : 130694993, "11" : 122082543, "12" : 120129022, "13" : 120421639, "14" : 124902244,
                "15" : 104043685, "16" : 98207768, "17" : 94987271, "18" : 90702639, "19" : 61431566, "X" : 171031299 }; #, "Y" : 91744698

###########################################################

infile = "../../data/03-selection-tests/bed/mm10-cds-10kb-window-overlap-sorted.bed";
# bedtools intersect -a 03-selection-tests/bed/mm10.ensGene.chromes.longest.cds.bed -b 02-genomic-windows/bed/mm10-10kb-windows.bed -wo -sortout > 03-selection-tests/bed/mm10-cds-10kb-window-overlap-sorted.bed

windowfile = "../../summary-data/02-genomic-windows/10kb-0.5-0.5-5mb-topo-counts.csv";

windows_indir = "../../data/02-genomic-windows/tree/";

concat_tree_file = "../../data/03-selection-tests/tree/concat/concat.cf.rooted.tree"

gene_tree_dir = "../../data/03-selection-tests/tree/loci/";

window_prefix = "10kb-0.5-0.5";
#chrome_tree_outdir = "../tree/penn-7spec/chrome/";
#window_tree_outdir = "../tree/penn-7spec/window/";

outfilename = "../../summary-data/03-selection-tests/mm10-cds-windows.csv";

####################

with open(outfilename, "w") as outfile:
    core.runTime("# Transcript/window overlaps", outfile);
    core.PWS("# Transcript/window overlaps: " + infile, outfile);
    core.PWS("# Window file:                " + windowfile, outfile);
    core.PWS("# Window indir :              " + windows_indir, outfile);
    core.PWS("# Gene trees:                 " + gene_tree_dir, outfile);
    #core.PWS("# Chrome tree outdir:         " + chrome_tree_outdir, outfile);
    #core.PWS("# Window tree outdir:         " + window_tree_outdir, outfile);
    core.PWS("# Output file:                " + outfilename, outfile);
    core.PWS("# ----------------");
    # Input info for log file

    core.PWS("# " + core.getDateTime() + " Reading concatenated tree:" + concat_tree_file, outfile);
    concat_tree = open(concat_tree_file, "r").read().strip();
    catinfo, concat_tree_parsed, root = tp.treeParse(concat_tree);
    #tree_parsed = re.sub("<[\d]+>", "", tree_parsed);
    # Parse the tree.

    concat_clade_topo = set([frozenset(tp.getClade(node, catinfo)) for node in catinfo if catinfo[node][2] != 'tip']);
    # For rooted trees, also get the clades to comapre to gene trees later    

    ####################

    core.PWS("# " + core.getDateTime() + " Getting chromosome trees...", outfile);
    chrome_trees = {};
    for chrome in chrome_sizes:
        chrome = "chr" + chrome;
        chrome_trees[chrome] = { "rooted" : "", "clades" : "", "unrooted" : "" };

        for treetype in ["", ".rooted"]:
            chrome_tree_file = chrome + "-concat.cf.tree" + treetype;
            chrome_tree_path = os.path.join(windows_indir, chrome, window_prefix, "concat", chrome_tree_file);
            chrome_tree = open(chrome_tree_path, "r").read().strip();
            #chrome_tree = chrome_tree.replace("HylAll", "hall").replace("RhySor", "rsor").replace("MasNat", "mnat").replace("PraDel", "pdel").replace("GraDol", "gdol").replace("RhaDil", "rdil").replace("mm10", "mmus");
            # Read the tree from the file.

            if treetype == ".rooted":
                ctinfo, chrome_tree_parsed, root = tp.treeParse(chrome_tree);
                #tree_parsed = re.sub("<[\d]+>", "", tree_parsed);
                # Parse the tree.

                chrome_clade_topo = set([frozenset(tp.getClade(node, ctinfo)) for node in ctinfo if ctinfo[node][2] != 'tip']);
                # For rooted trees, also get the clades to comapre to gene trees later

                chrome_trees[chrome]["rooted"] = chrome_tree;
                chrome_trees[chrome]["clades"] = chrome_clade_topo;
                
            else:
                chrome_trees[chrome]["unrooted"] = chrome_tree;
                # For unrooted trees just save the tree
    core.PWS("# ----------------");
    # Read the chromosome trees.

    ####################

    core.PWS("# " + core.getDateTime() + " Reading windows: " + windowfile, outfile);
    windows = {};
    with open(windowfile) as csvfile:
        first = True;
        windows_reader = csv.reader(csvfile);
        for line in windows_reader:
            # print(line);
            if line[0][0] == "#":
                continue;
            # Skip the log lines from 09_topo_counter.py
            if first:
                headers = line + ["transcript", "cds.start", "cds.end", "overlap", "gene.tree", "gene.tree.match.window.tree", "gene.tree.match.chrome.tree", "gene.tree.match.concat.tree"];
                first = False;
                continue;
            # Get the headers and add the additional columns.

            #line = [ l.replace("HylAll", "hall").replace("RhySor", "rsor").replace("MasNat", "mnat").replace("PraDel", "pdel").replace("GraDol", "gdol").replace("RhaDil", "rdil").replace("mm10", "mmus") for l in line ];
            # Switch all the labels in the line to match the ones I used for CDS stuff...

            line[13] = '"' + line[13] + '"';
            line[17] = '"' + line[17] + '"';
            line[20] = '"' + line[20] + '"';
            #line[26] = '"' + line[26] + '"';

            windows[line[0]] = line;
    core.PWS("# Windows read: " + str(len(windows)), outfile);

    ####################

    core.PWS("# " + core.getDateTime() + " Reading overlaps and copying trees: " + infile, outfile);
    outfile.write(",".join(headers) + "\n");
    window_genes = defaultdict(dict);
    genes_read, genes_skipped = 0, 0;
    for line in open(infile):
        line = line.strip().split("\t");
        window, cds_start, cds_end, tid, overlap, chrome = line[15], line[1], line[2], line[3], line[16], line[0];
        #gid, tid = ids.split(";");
        # Read the transcripts and overlaps with windows.

        gene_tree_file = os.path.join(gene_tree_dir, tid, tid + "-7spec-cds-NT.pretrim.macse.trim.filter.NT.treefile.rooted");
        if chrome == "chrY":
            continue;
        if not os.path.isfile(gene_tree_file):
            genes_skipped += 1;
            continue;

        gt_cat_match = "FALSE";
        gt_ct_match = "FALSE";
        gt_wt_match = "FALSE";
        # Tree match variables

        # chrome_treedir = os.path.join(chrome_tree_outdir, gid);
        # if not os.path.isdir(chrome_treedir):
        #     os.system("mkdir " + chrome_treedir);
        # Make the directory for the chromosome tree for the current gene

        # rooted_chrome_treefile = os.path.join(chrome_treedir, gid + "-" + chrome + "-10kb-rooted.treefile");
        # unrooted_chrome_treefile = os.path.join(chrome_treedir, gid + "-" + chrome + "-10kb.treefile");
        # Chromosome tree output files

        # with open(rooted_chrome_treefile , "w") as cfile:
        #     cfile.write(chrome_trees[chrome]['rooted']);

        # with open(unrooted_chrome_treefile , "w") as cfile:
        #     cfile.write(chrome_trees[chrome]['unrooted']);
        # Write chromosome tree for current gene


        # window_treedir = os.path.join(window_tree_outdir, gid);
        # if not os.path.isdir(window_treedir):
        #     os.system("mkdir " + window_treedir);
        # Make the directory for the window trees for the current gene
        ## CHROMOSOME TREES

        window_clade_topo = False;
        if all(filt != "FILTER" for filt in [ windows[window][6], windows[window][7] ]):

            #for treetype in ["", ".rooted"]:
            for treetype in [".rooted"]:
                window_str = window.replace(":", "-") + "-mafft-trimal";
                window_tree_file_in = window_str + ".treefile" + treetype
                window_tree_path = os.path.join(windows_indir, chrome, window_prefix, "loci", window.replace(":", "-"), window_tree_file_in);

                if not os.path.isfile(window_tree_path):
                    print(window_tree_path);
                else:
                    window_tree = open(window_tree_path, "r").read().strip();
                # window_tree = window_tree.replace("HylAll", "hall").replace("RhySor", "rsor").replace("MasNat", "mnat").replace("PraDel", "pdel").replace("GraDol", "gdol").replace("RhaDil", "rdil").replace("mm10", "mmus");
                # Read the tree from the file.

                #window_treefile_out = os.path.join(window_treedir, gid + "-" + window.replace(":", "-") + treetype + ".treefile");
                # Get the output file name for the window tree for the current gene

                if treetype == ".rooted":
                    wtinfo, window_tree_parsed, root = tp.treeParse(window_tree);
                    #tree_parsed = re.sub("<[\d]+>", "", tree_parsed);
                    # Parse the tree.

                    window_clade_topo = set([frozenset(tp.getClade(node, wtinfo)) for node in wtinfo if wtinfo[node][2] != 'tip']);
                    # For rooted trees, also get the clades to comapre to gene trees later

                # with open(window_treefile_out , "w") as wfile:
                #     wfile.write(window_tree);
                # Write the window tree to the output file
        ## WINDOW TREE

        gene_tree = open(gene_tree_file, "r").read().strip();
        # Read the tree from the file.
        gtinfo, gene_tree_parsed, root = tp.treeParse(gene_tree);
        #tree_parsed = re.sub("<[\d]+>", "", tree_parsed);
        # Parse the tree.
        gene_clade_topo = set([frozenset(tp.getClade(node, gtinfo)) for node in gtinfo if gtinfo[node][2] != 'tip']);
        # Get the topology as the set of all clades present in the tree.
        ### GENE TREE

        if gene_clade_topo == concat_clade_topo:
            gt_cat_match = "TRUE";
        if gene_clade_topo == chrome_trees[chrome]["clades"]:
            gt_ct_match = "TRUE";
        if gene_clade_topo == window_clade_topo:
            gt_wt_match = "TRUE";
        # Check if any of the topologies match

        outline = windows[window] + [tid, cds_start, cds_end, overlap, '"' + gene_tree + '"', gt_wt_match, gt_ct_match, gt_cat_match];
        # Get the output

        outfile.write(",".join(outline) + "\n");
        # Write the output

        genes_read += 1;
        #sys.exit();
    core.PWS("# Genes read:    " + str(genes_read), outfile);    
    core.PWS("# Genes skipped: " + str(genes_skipped), outfile);
    core.PWS("# ----------------");
