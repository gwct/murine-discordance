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
sys.path.append("/home/gt156213e/bin/core/python/lib/");
import core

############################################################

beddir = "/home/gt156213e/genomes/cds/bed/";
# Input directory with bed files containing windows for each species

outdir = "/home/gt156213e/genomes/cds/seq/";
if not os.path.isdir(outdir):
    os.makedirs(outdir);
# Output directory

logfilename = os.path.join("/home/gt156213e/genomes/cds/scripts/logs/", "combine-seqs.log");
# Log file for this run.

transcript_file = "/home/gt156213e/genomes/cds/annotations/mm10-degenotate/longest-transcripts.txt";
# The file with the transcript IDs

outfilename = os.path.join("/home/gt156213e/genomes/cds/data/", "combine-seqs.tsv");
# Output file for this run

specs = ["mm10", "gdol", "hall", "mnat", "pdel", "rdil", "rsor"];
# The list of species to look up

## File paths
###################

pad = 35;
# Spacing pad

windows = {};
# Main output table

with open(logfilename, "w") as logfile:
    core.runTime("# Combining CDS seqs", logfile);
    core.PWS("# ----------", logfile);
    # Run time info for logging.
    ###################

    core.PWS("# " + core.getDateTime() + " Reading transcript IDs...", logfile);
    transcripts = open(transcript_file, "r").read().split("\n");
    core.PWS(core.spacedOut("# " +  "Transcript IDs read:", pad) + str(len(transcripts)), logfile);

    ###################

    core.PWS("# " + core.getDateTime() + " Reading seq beds...", logfile);

    spec_dict = {};
    for spec in specs:
        core.PWS("# " + core.getDateTime() + " " + spec, logfile);
        num_deleted = 0;
        spec_dict[spec] = { transcript : {'status' : "NA", 'seq' : "NA", 'coords' : "NA" } for transcript in transcripts };
        # Initialize the window dict for each species

        ###################

        bed_seq_file = os.path.join(beddir, spec + "-cds-seqs.tab");
        if spec != "mm10":
            deleted_bed_file = os.path.join(beddir, spec + "-cds.bed.unmap");
        else:
            deleted_bed_file = "";
        # Bed files for each species -- species that are liftedOver may have unmapped windows to read

        ###################

        if spec != "mm10":
            for line in open(deleted_bed_file):
                if line[0] == "#" or line == "\n":
                    continue;
                deleted_transcript = line.strip().split("\t")[3];
                spec_dict[spec][deleted_transcript]['status'] = "deleted";
                num_deleted += 1;
        core.PWS("# " + core.getDateTime() + " " + str(num_deleted) + " deleted", logfile);
        ## For the non-reference species, look up the deleted windows to mark them as missing later

        ###################

        for line in open(bed_seq_file):
            line = line.strip().split("\t");

            seq_id = line[0];
            transcript, coords = seq_id.split("::");
            spec_dict[spec][transcript]['coords'] = coords;
            spec_dict[spec][transcript]['seq'] = line[1];
            # Read every line in the bed file for the species and parse out some info, including the sequence

        ###################

    ## End spec loop
    ###################

    core.PWS("# " + core.getDateTime() + " Writing sequences...", logfile);

    total_written = 0;
    total = 0;
    non_div_three = 0;
    # Some counters for all transcripts

    for transcript in transcripts:
        num_missing = 0;
        # A counter for the current transcript

        cur_seqs = {};
        # Initialize seq dict for current transcript

        for spec in spec_dict:
            if spec_dict[spec][transcript]['status'] == "deleted":
                num_missing += 1;
            ## If the sequence is deleted in the current species, note that here and count it as missing
            ## Also fill in the sequence with Ns

            else:
                header = ">" + spec + " " + transcript + ":" + spec_dict[spec][transcript]['coords'];
                cur_seqs[header] = spec_dict[spec][transcript]['seq'];
                # Get the sequence for the current species and window
        ## End spec seq loop

        ###################

        if num_missing == 0:
            cur_seq_out = os.path.join(outdir, transcript + "-7spec-cds.fa");
            # Get the output filename for the current sequence

            with open(cur_seq_out, "w") as seqout:
                for seq in cur_seqs:
                    total += 1;
                    if len(cur_seqs[seq]) % 3 != 0:
                        non_div_three += 1
                    seqout.write(seq + "\n" + cur_seqs[seq] + "\n");
                total_written += 1;
        ## If the window passed the repeat filter, write the sequences
    ## End window seq loop

    core.PWS(core.spacedOut("# Total transcripts written:", pad) + str(total_written), logfile);
    print(total);
    print(non_div_three);

    ###################

    core.PWS("# " + core.getDateTime() + " Writing output table...", logfile);

    with open(outfilename, "w") as outfile:

        headers = ["transcript"] + specs;
        core.runTime("# Rodent CDS retrieval", outfile, False);
        core.PWS("\t".join(headers), outfile);
        # Write the headers of the output table

        ###################

        for transcript in transcripts:
            outline = [transcript];
            for spec in specs:
                if spec_dict[spec][transcript]['status'] == "deleted":
                    outline.append("N");
                else:
                    outline.append("Y");

            outline = "\t".join( outline );
            outfile.write(outline + "\n");
        # Write the output table

    ## Close output file
    ###################
## Close log file
###################