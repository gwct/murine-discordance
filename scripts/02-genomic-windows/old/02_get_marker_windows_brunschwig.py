#!/usr/bin/python3
############################################################
# For Penn genomes, 05.2021
# Parses the recombination markers into windows. Brunschwig
# rates.
############################################################

import sys, os, core, time

############################################################

window_size = 100;
wsize = window_size * 1000;
# Window size in kb

markerfile = "/mnt/beegfs/gt156213e/penn-genomes/windows/data/recombination-markers/brunschwig-rates/mm10-markers.bed";
logfilename = "logs/get-marker-windows-brunscwhig" + str(window_size) + "kb.log";
outfilename = markerfile.replace(".bed", "-" + str(window_size) + "kb.csv");
# File names

chrome_sizes = {"chr1" : 195471971, "chr2" : 182113224, "chr3" : 160039680, "chr4" : 156508116, "chr5" : 151834684, "chr6" : 149736546, "chr7" : 145441459,
                "chr8" : 129401213, "chr9" : 124595110, "chr10" : 130694993, "chr11" : 122082543, "chr12" : 120129022, "chr13" : 120421639, "chr14" : 124902244,
                "chr15" : 104043685, "chr16" : 98207768, "chr17" : 94987271, "chr18" : 90702639, "chr19" : 61431566}#, "chrX" : 171031299 } #, "Y" : 91744698};
# Mouse chromosome sizes from /mnt/beegfs/gt156213e/ref-genomes/mm10/mm10.chrome.sizes"

pad = 20;

with open(logfilename, "w") as logfile:
    core.runTime("# Rodent marker windows", logfile);
    core.PWS(core.spacedOut("# Marker file:", pad) + markerfile, logfile);
    core.PWS(core.spacedOut("# Window size:", pad) + str(window_size) + "kb", logfile);
    core.PWS(core.spacedOut("# Output file:", pad) + outfilename, logfile);
    core.PWS(core.spacedOut("# Log file:", pad) + logfilename, logfile);
    core.PWS("# ----------------", logfile);
    # Run time info for logging.
    ###################

    core.PWS("# " + core.getDateTime() + " Reading markers...", logfile);
    markers = {};
    first = True;
    for line in open(markerfile):
        line = line.strip().split("\t");
        # if first:
        #     headers = line;
        #     first = False;
        #     continue;

        chrome, pos, pos_end, ids, rho, strand = line
        marker_id = chrome + ":" + pos;

        markers[marker_id] = { 'chr' : chrome, 'pos' : int(pos), 'rho' : rho, 'window' : "", 'start' : "", 'end' : "" };

    core.PWS(core.spacedOut("# Markers read:", pad) + str(len(markers)), logfile);
    core.PWS("# ----------------", logfile);

    core.PWS("# " + core.getDateTime() + " Generating windows...", logfile);
    windows = {};

    for chrome in chrome_sizes:
        start = 1;
        cur_wsize = wsize;
        stop = cur_wsize;
        chr_end = chrome_sizes[chrome];
        num_windows = 0;
        
        while start < chr_end:
            if stop > chr_end:
                stop = chr_end;
                cur_wsize = chr_end - start;
            # For the last window;

            window_name = chrome + ":" + str(start) + "-" + str(stop);
            windows[window_name] = { 'chr' : chrome, 'start' : start, 'stop' : stop, 'num-markers' : 0 };
            num_windows += 1;

            start = stop + 1;
            stop = stop + cur_wsize;
            # Increment start and stop of window size.

        core.PWS(core.spacedOut("# " + chrome + " windows:", pad) + str(num_windows), logfile);    
    core.PWS("# ----------------", logfile);

    core.PWS("# " + core.getDateTime() + " Sorting markers...", logfile);
    markers_chr_sorted = [];
    for chrome in chrome_sizes:
        #print(chrome);
        cur_chr_markers = { m : int(markers[m]['pos']) for m in markers if markers[m]['chr'] == chrome };
        markers_chr_sorted += sorted(cur_chr_markers, key=lambda k: (cur_chr_markers[k], k));
    core.PWS("# ----------------", logfile);
    # Sort the markers by chromosome and position so we can go through them in order later and fill in
    # windows with no markers.

    core.PWS("# " + core.getDateTime() + " Placing markers in windows...", logfile);
    markers_placed = 0;
    for window in windows:
        for marker in markers:
            #print(window, marker, markers[marker]['chr'], windows[window]['chr']);
            if markers[marker]['chr'] == windows[window]['chr']:
                if int(markers[marker]['pos']) >= windows[window]['start'] and int(markers[marker]['pos']) <= windows[window]['stop']:
                    markers[marker]['window'] = window;
                    markers[marker]['start'] = windows[window]['start'];
                    markers[marker]['end'] = windows[window]['stop'];

                    windows[window]['num-markers'] += 1;

                    markers_placed += 1;
    core.PWS(core.spacedOut("# Markers placed:", pad) + str(markers_placed), logfile);

    no_marker_windows = [];
    for window in windows:
        if windows[window]['num-markers'] == 0:
            core.PWS(core.spacedOut("# ! No markers found in window:", pad) + window, logfile);
            no_marker_windows.append(window);
    # Get the windows with no markers.
    core.PWS("# ----------------", logfile);

    core.PWS("# " + core.getDateTime() + " Writing markers...", logfile);
    markers_written = 0;
    prev_window, next_window = "chr1:1-" + str(wsize), "chr1:1-" + str(wsize);
    with open(outfilename, "w") as outfile:
        headers = ['marker', 'marker.chr', 'marker.pos', 'rho', 'window', 'win.start', 'win.end', 'num.win.markers'];
        outfile.write(",".join(headers) + "\n");

        no_marker_window, prev_window, prev_window_status, prev_window_start = 1, "first", "last", 0;
        # Keep track of the window of the last marker and the expected next window so we can fill in windows with missing markers.
        for marker in markers_chr_sorted:
            m = markers[marker];

            cur_window = m['window'];
            cur_window_chrome, cur_window_pos = cur_window.split(":");
            cur_window_start, cur_window_end = [ int(pos) for pos in cur_window_pos.split("-") ];
            # Get the current window info.

            if cur_window != prev_window:
            # We only want to check the next expected window if the previous marker was the last one in the window (so if it is not the same as the current window)
                while cur_window != next_window:
                # If the next expected window is not the current window, we run through windows until we match the current window.

                    #print(prev_window, cur_window, next_window);
                    nm_window_chrome, nm_window_pos = next_window.split(":");
                    nm_window_start, nm_window_end = [ int(pos) for pos in nm_window_pos.split("-") ];
                    # Parse the next expected window.

                    core.PWS(core.spacedOut("# ! Inserting no marker window:", pad) + next_window, logfile);
                    outline = [ "no-markers-" + str(no_marker_window), nm_window_chrome, "NA", "NA", next_window, str(nm_window_start), str(nm_window_end), "0" ];

                    outfile.write(",".join(outline) + "\n");                  
                    no_marker_window += 1;
                    # Insert a holder marker for this window with no markers.

                    next_chrome = nm_window_chrome;
                    next_start = nm_window_end + 1;
                    next_end = next_start + wsize - 1;
                    # Get the next expected window after the current expected window.

                    if next_end > chrome_sizes[nm_window_chrome]:
                        next_end = chrome_sizes[nm_window_chrome];
                    # If the end of the next expected window exceeds the current chromosome size, replace the end position with the chromosome size.

                    if next_start > chrome_sizes[nm_window_chrome]:
                    # If the start position of the next expected window exceeds the current chromosome size, we move to the next chromosome.
                        if nm_window_chrome != "chr19":
                        # 19 is last so we don't need to do anything.
                            next_chrome = "chr" + str(int(nm_window_chrome.replace("chr", ""))+1);
                            # For all chromosomes except 19 the next chromosome is simply +1. For 19 it is X.

                            next_start = 1;
                            next_end = next_start + wsize - 1;
                            # The start and end positions of the first window on the next chromosome.

                    next_window = next_chrome + ":" + str(next_start) + "-" + str(next_end);
                    # Re-pack the next window as a string.
                    #print(next_window);
                    #time.sleep(1);
                    #sys.exit();
                # Repeat the next window check until we get to the current window.
            # Only check the next window if the previous marker was the last one in the window.

            outline = [ marker, m['chr'], str(m['pos']), m['rho'], m['window'], str(m['start']), str(m['end']) ];
            num_win_markers = windows[m['window']]['num-markers'];
            outline.append(str(num_win_markers));
            outfile.write(",".join(outline) + "\n");
            markers_written += 1;
            # Write the current marker to the output file.

            next_chrome = cur_window_chrome;
            next_start = cur_window_end + 1;
            next_end = next_start + wsize - 1;
            # Get the next expected window after the current expected window.

            if next_end > chrome_sizes[cur_window_chrome]:
                next_end = chrome_sizes[cur_window_chrome];
            # If the end of the next expected window exceeds the current chromosome size, replace the end position with the chromosome size.

            if next_start > chrome_sizes[cur_window_chrome]:
            # If the start position of the next expected window exceeds the current chromosome size, we move to the next chromosome.
                if cur_window_chrome != "chr19":
                # X is last so we don't need to do anything.
                    next_chrome = "chr" + str(int(cur_window_chrome.replace("chr", ""))+1);
                    # For all chromosomes except 19 the next chromosome is simply +1. For 19 it is X.

                    next_start = 1;
                    next_end = next_start + wsize - 1;
                    # The start and end positions of the first window on the next chromosome.

            next_window = next_chrome + ":" + str(next_start) + "-" + str(next_end);
            # Re-pack the next window as a string.

            prev_window = cur_window;
            # Assign the current window as the next previous window.

        if cur_window != next_window:
            nm_window_chrome, nm_window_pos = next_window.split(":");
            nm_window_start, nm_window_end = [ int(pos) for pos in nm_window_pos.split("-") ];
            # Parse the next expected window.

            core.PWS(core.spacedOut("# ! Inserting no marker window:", pad) + next_window, logfile);
            outline = [ "no-markers-" + str(no_marker_window), nm_window_chrome, "NA", "NA", cur_window, str(nm_window_start), str(nm_window_end), "0" ];
            #outline = [ "no-markers-" + str(no_marker_window), nm_window_chrome, "NA", "NA", "NA", "NA", nm_window_chrome, "NA", next_window, str(nm_window_start), str(nm_window_end), "0" ];
            outfile.write(",".join(outline) + "\n");                  
            no_marker_window += 1;
            # Insert a holder marker for this window with no markers.
        # Add in the last window if it has no markers.

    core.PWS(core.spacedOut("# Markers written:", pad) + str(markers_written), logfile);
    core.PWS("# ----------------", logfile);