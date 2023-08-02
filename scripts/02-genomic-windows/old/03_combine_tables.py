#!/usr/bin/python
############################################################
# For Penn genomes, 06.2021
# Combines the tree and recombination tables
############################################################

import sys, os, core
from scipy import stats

############################################################

window_size = 1;
window_size_str = str(window_size);
wsize = window_size * 1000000;
# Window size in MB

infilename = "../data/recombination-markers/Revised_HSmap_SNPs-mm10-" + window_size_str + "mb.csv";
outfilename = "../data/recombination-markers/Revised_HSmap_SNPs-mm10-" + window_size_str + "mb-rates.csv";
logfilename = "logs/recombination-regression-" + window_size_str + "mb.log"
pad = 25;

with open(logfilename, "w") as logfile:
    core.runTime("# Rodent window recombination rates", logfile);
    core.PWS(core.spacedOut("# Marker file:", pad) + infilename, logfile);
    core.PWS(core.spacedOut("# Output file:", pad) + outfilename, logfile);
    core.PWS("# ----------------", logfile);
    # Run time info for logging.
    ###################

    core.PWS("# " + core.getDateTime() + " Reading markers...", logfile);
    markers, rates = {}, {};
    num_markers, num_windows = 0, 0;
    first = True;
    for line in open(infilename):
        line = line.strip().split(",");
        if first:
            headers = line;
            first = False;
            continue;

        chrome, window, pos = line[6], line[8], line[7];

        if pos != "NA":
            pos = float(pos);
            if chrome == "X":
                cm = float(line[3]);
            else:
                cm = float(line[5]);
        else:
            cm = "NA"

        if window not in markers:
            markers[window] = {'bp' : [], 'cm' : []};
            num_windows += 1;
        if window not in rates:
            rates[window] = "NA";

        markers[window]['bp'].append(pos);
        markers[window]['cm'].append(cm);
        num_markers += 1;

    core.PWS(core.spacedOut("# Markers:", pad) + str(num_markers), logfile); 
    core.PWS(core.spacedOut("# Windows:", pad) + str(num_windows), logfile); 
    core.PWS("# ----------------", logfile);

    core.PWS("# " + core.getDateTime() + " Calculating rates per window (at least 5 markers) via linear regression...", logfile);
    num_rates = 0;
    for window in markers:
        assert len(markers[window]['bp']) == len(markers[window]['cm']), "\nUnequal number of marker coordinates:\n" + window;

        if len(markers[window]['bp']) >= 5:
            slope, intercept, r_value, p_value, std_err = stats.linregress(markers[window]['bp'], markers[window]['cm']);
            rates[window] = slope;
            num_rates += 1;
    core.PWS(core.spacedOut("# Rates caclulated:", pad) + str(num_rates), logfile);
    core.PWS("# ----------------", logfile);

    core.PWS("# " + core.getDateTime() + " Writing output...", logfile);
    with open(outfilename, "w") as outfile:
        headers.append("slope");
        outfile.write(",".join(headers) + "\n");

        markers_written = 0;
        first = True;
        for line in open(infilename):
            if first:
                first = False;
                continue;
            line = line.strip().split(",");
            window = line[8];
            slope = str(rates[window]);
            line.append(slope);

            outfile.write(",".join(line) + "\n");
            markers_written += 1;
    core.PWS(core.spacedOut("# Markers written:", pad) + str(markers_written), logfile);
    print("Done!");
    core.PWS("# ----------------", logfile);