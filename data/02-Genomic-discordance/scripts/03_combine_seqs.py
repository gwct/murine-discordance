#!/usr/bin/python
############################################################
# For rodent genomes, 04.2023
#
# Takes the window sequences in the bed files and combines
# them based on window id
# Also generates the windows table
#
############################################################

import sys
import os
import argparse
import core

############################################################

parser = argparse.ArgumentParser(description="Get rodent chromosome window coordinates");
parser.add_argument("-w", dest="window_size", help="The size of the sliding window in kb.", type=int, default=False);
parser.add_argument("-r", dest="repeat_thresh", help="Windows that have a proportion of their sequence overlapping with repeats greater than this will be excluded. Default: 0.5", type=float, default=0.5);
parser.add_argument("-m", dest="missing_thresh", help="Windows that have more than 3 sequences with a proportion of their sequence that is missing data (Ns) greater than this will be excluded. Default: 0.5", type=float, default=0.5);
args = parser.parse_args();

if not args.window_size:
    sys.exit(" * Error 1: Window size in kb (-w) must be defined.");
if args.window_size < 1:
    sys.exit(" * Error 2: Window size in kb (-w) must be a positive integer.");
else:
    wsize_orig = args.window_size * 1000;
    wsize = args.window_size * 1000;

window_size_str = str(args.window_size)
repeat_thresh_str = str(args.repeat_thresh);
missing_thresh_str = str(args.missing_thresh);
# Covnert params to string for logging

## Check input options
###################

beddir = "/mnt/beegfs/gt156213e/penn-genomes/windows/pseudo/bed/";
# Input directory with bed files containing windows for each species

outdir = "/mnt/beegfs/gt156213e/penn-genomes/windows/pseudo/seq/";
if not os.path.isdir(outdir):
    os.makedirs(outdir);
# Output directory

param_str = window_size_str + "kb-" + repeat_thresh_str + "-" + str(args.missing_thresh);
# A string of the parameters used in file and dir names

logfilename = os.path.join("/mnt/beegfs/gt156213e/penn-genomes/windows/pseudo/scripts/logs/", "get-window-seqs-" + param_str + ".log");
# Log file for this run.

outfilename = os.path.join("/mnt/beegfs/gt156213e/penn-genomes/windows/pseudo/data/", param_str + "-windows.tsv");
# Output file for this run

specs = ["mm10", "gdol", "hall", "mnat", "pdel", "rdil", "rsor"];
# The list of species to look up

chromes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX"];
# The list of chromosomes to retrieve

## File paths
###################

pad = 35;
# Spacing pad

windows = {};
# Main output table

with open(logfilename, "w") as logfile:
    core.runTime("# Rodent windows", logfile);
    core.PWS(core.spacedOut("# Window size (-w):", pad) + str(args.window_size) + "kb", logfile);
    core.PWS(core.spacedOut("# Log file:", pad) + logfilename, logfile);
    core.PWS("# ----------", logfile);
    # Run time info for logging.
    ###################

    core.PWS("# " + core.getDateTime() + " Reading windows and checking repeats...", logfile);
    repeat_bed = os.path.join(beddir, "mm10-" + str(args.window_size) + "kb-repeat-coverage.bed");
    # The bed file with repeat coverage

    total_windows, repeat_filtered, repeat_frac_sum = 0, 0, 0;
    for line in open(repeat_bed):
        line = line.strip().split("\t");
        # Parse the line

        if line[0] not in chromes:
            continue;
        # Skip any windows not in the list of chromes we want to retrieve

        window = line[3];
        # Get the window name

        windows[window] = { 'chr' : line[0], 'start' : line[1], 'end' : line[2], 'window' : window,
                             'repeat.frac' : line[7], 'repeat.filter' : "NA", 
                             'num.missing' : "NA", 'missing.filter' : "NA", 
                             'uniq.filter' : "NA" 
                            };
        # Initialize the main dict for this table

        if float(windows[window]['repeat.frac']) >= args.repeat_thresh:
            windows[window]['repeat.filter'] = "FILTER";
            repeat_filtered += 1;
        else:
            windows[window]['repeat.filter'] = "PASS";
        # Mark the repeat filter for this window based on the current threshold

        repeat_frac_sum += float(windows[window]['repeat.frac']);
        # Add the total repeat fraction to get a rough average per window at the end

        total_windows += 1;
    ## End window loop

    core.PWS(core.spacedOut("# Total " + window_size_str + "kb windows:", pad) + str(total_windows), logfile);
    core.PWS(core.spacedOut("# Filtered " + window_size_str + "kb windows:", pad) + str(repeat_filtered), logfile);
    core.PWS(core.spacedOut("# Avg. repeat per " + window_size_str + "kb window:", pad) + str(repeat_frac_sum / float(total_windows)), logfile);
    ###################

    core.PWS("# " + core.getDateTime() + " Checking chrome dirs...", logfile);
    
    for chrome in chromes:
        chrdir = os.path.join(outdir, chrome);
        if not os.path.isdir(chrdir):
            os.makedirs(chrdir);
        chr_subdir = os.path.join(chrdir, param_str);
        if not os.path.isdir(chr_subdir):
            os.makedirs(chr_subdir);
    ## End chrome directory loop

    ###################

    core.PWS("# " + core.getDateTime() + " Reading seq beds...", logfile);

    spec_dict = {};
    for spec in specs:
        core.PWS("# " + core.getDateTime() + " " + spec, logfile);
        num_deleted, num_missing = 0, 0;
        spec_dict[spec] = { window : {'status' : "NA", 'missing.frac' : "NA", 'seq' : "NA", 'chr' : "NA", 'start' : "NA", 'end' : "NA" } for window in windows };
        # Initialize the window dict for each species

        ###################

        bed_seq_file = os.path.join(beddir, spec + "-" + window_size_str + "kb-seqs.bed");
        deleted_bed_file = os.path.join(beddir, spec + "-" + window_size_str + "kb.bed.unmap");
        # Bed files for each species -- species that are liftedOver may have unmapped windows to read

        ###################

        if spec != "mm10":
            for line in open(deleted_bed_file):
                if line[0] == "#" or line == "\n":
                    continue;
                deleted_window = line.strip().split("\t")[3];
                spec_dict[spec][deleted_window]['status'] = "deleted";
                num_deleted += 1;
        core.PWS("# " + core.getDateTime() + " " + str(num_deleted) + " deleted", logfile);
        ## For the non-reference species, look up the deleted windows to mark them as missing later

        ###################

        for line in open(bed_seq_file):
            line = line.strip().split("\t");

            if line[0] not in chromes:
                continue;
            # Skip any windows not in the list of chromes we want to retrieve

            window = line[3];
            spec_dict[spec][window]['chr'] = line[0];
            spec_dict[spec][window]['start'] = line[1];
            spec_dict[spec][window]['end'] = line[2];
            spec_dict[spec][window]['seq'] = line[4];
            spec_dict[spec][window]['missing.frac'] = line[7];
            # Read every line in the bed file for the species and parse out some info, including the sequence

            if float(spec_dict[spec][window]['missing.frac']) > args.missing_thresh:
                spec_dict[spec][window]['status'] = "missing";
                num_missing += 1;
            # Check if the amount of missing sequence is above the current threshold and mark it if so
            
        core.PWS("# " + core.getDateTime() + " " + str(num_missing) + " missing", logfile);
        ###################

    ## End spec loop
    ###################

    core.PWS("# " + core.getDateTime() + " Writing sequences...", logfile);

    total_repeat, total_missing, total_written = 0, 0, 0;
    # Some counters for all windows

    for window in windows:
        num_missing = 0;
        # A counter for the current window

        cur_seqs = {};
        # Initialize seq dict for current window

        for spec in spec_dict:
            if spec_dict[spec][window]['status'] == "deleted":
                header = ">" + spec + " deleted";
                cur_seqs[header] = "N" * wsize;
                num_missing += 1;
            ## If the sequence is deleted in the current species, note that here and count it as missing
            ## Also fill in the sequence with Ns

            else:
                header = ">" + spec + " " + spec_dict[spec][window]['chr'] + ":" + spec_dict[spec][window]['start'] + "-" + spec_dict[spec][window]['end'];
                cur_seqs[header] = spec_dict[spec][window]['seq'];
                # Get the sequence for the current species and window

                if spec_dict[spec][window]['status'] == 'missing':
                    num_missing += 1;
                # Count if it was above the missing threshold
        ## End spec seq loop

        ###################

        windows[window]['num.missing'] = str(num_missing);
        # Add the number missing to the main output table

        if num_missing > 3:
            total_missing += 1;
            windows[window]['missing.filter'] = "FILTER";
        ## Check if more than 3 sequences are over the missing filter and mark the filter if so
        else:
            windows[window]['missing.filter'] = "PASS";
            # Set the missing.filter to PASS here

            if windows[window]['repeat.filter'] == "PASS":
                window_str = "-".join([ windows[window]['chr'], windows[window]['start'], windows[window]['end'] ]);
                cur_seq_out = os.path.join(outdir, windows[window]['chr'], param_str, window_str + ".fa");
                # Get the output filename for the current sequence

                with open(cur_seq_out, "w") as seqout:
                    for seq in cur_seqs:
                        seqout.write(seq + "\n" + cur_seqs[seq] + "\n");
                    total_written += 1;
            ## If the window passed the repeat filter, write the sequences
            else:
                total_repeat += 1;
            ## Otherwise just increment a counter and don't write the sequence
        ## End missing block
    ## End window seq loop

    core.PWS(core.spacedOut("# Total " + window_size_str + "kb windows written:", pad) + str(total_written), logfile);
    core.PWS(core.spacedOut("# Filtered MISSING " + window_size_str + "kb windows:", pad) + str(total_missing), logfile);
    core.PWS(core.spacedOut("# Filtered REPEAT " + window_size_str + "kb windows:", pad) + str(total_repeat), logfile);

    ###################

    core.PWS("# " + core.getDateTime() + " Writing output table...", logfile);

    with open(outfilename, "w") as outfile:

        headers = ["chr","start","end","window","repeat.frac","repeat.filter","num.missing","missing.filter","aln.filter"];
        core.runTime("# Rodent window filtering", outfile, False);
        core.PWS("\t".join(headers), outfile);
        # Write the headers of the output table

        ###################

        for window in windows:
            outline = "\t".join( [ str(windows[window][col]) for col in headers ] );
            outfile.write(outline + "\n");
        # Write the output table

    ## Close output file
    ###################
## Close log file
###################