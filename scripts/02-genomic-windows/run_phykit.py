#!/usr/bin/python
############################################################
# For Phodopus genomes, 06.2021
# Runs phykit to get alignment stats from windows.
############################################################

import sys, os, argparse, core, subprocess, multiprocessing as mp

############################################################
# Functions
def runCMD(cmd):
# Run a command and deal with the output nicely.
    #core.PWS(core.spacedOut("# --> Executing:", pad) + cmd, std_stream=False, o_stream=logfile);
    cmd_result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE);
    if any(ecode in cmd_result.stderr.decode() for ecode in ['error', 'Error', 'ERROR', 'Exception', 'Could not build fai index', 
                                                                'AssertionError', "Can't read file", "Killed", "No such file or directory", 
                                                                "Failed", "returning empty sequence"]):
        core.PWS("# !! CMD ERROR: The following command returned an error:\n\n" + cmd);
        core.PWS("# !! STDOUT:\n" + cmd_result.stdout.decode() + "\n\n");
        core.PWS("# !! STDERR:\n" + cmd_result.stderr.decode() + "\n\n");
        sys.exit();

    return cmd_result.stdout.decode();

#########################

def runPhykit(file_list, cur_chrome, in_dir):
    outlines = [];

    aln_counter = 1;
    num_alns = len(file_list);
    num_alns_str = str(num_alns);

    for f in file_list:
        if aln_counter == 1 or aln_counter % 100 == 0:
            core.PWS("# " + core.getDateTime() + " | " + cur_chrome + " -> " + str(aln_counter) + " / " + num_alns_str);
        aln_counter += 1;

        window = f.replace("-mafft-trimal.fa", "");
        window = window.split("-");
        window = window[0] + ":" + window[1] + "-" + window[2];

        #infile = os.path.join(indir, cur_chrome + "-trimal", f);

        #print(infile);
        vs_cmd = "phykit vs " + f;
        stats = runCMD(vs_cmd);
        stats = stats.strip().split("\t");
        outline = [window, stats[1], stats[0], stats[2]];

        pis_cmd = "phykit pis " + f;
        stats = runCMD(pis_cmd);
        stats = stats.strip().split("\t");

        outline += [stats[0], stats[2]];

        outlines.append(",".join(outline) + "\n");

    return outlines;

############################################################

indir = "../../data/02-genomic-windows/aln/02-trimal/";
outfilename = "../../summary-data/02-genomic-windows/phykit-out.csv";

#10kb-0.5-0.5

core.PWS("# " + core.getDateTime() + " | Reading file names: " + indir);
file_dict = {};
for d in os.listdir(indir):
    if not d.startswith("chr"):
        continue;

    alndir = os.path.join(indir, d, "10kb-0.5-0.5");
    #chrome = d.replace("-trimal", "");

    file_dict[d] = [ os.path.join(alndir, aln) for aln in os.listdir(alndir) if aln.endswith("-mafft-trimal.fa") ];

core.PWS("# " + core.getDateTime() + " | Starting counts");
outlines_main = [];
with mp.Pool(processes=20) as pool:
    for result in pool.starmap(runPhykit, ((file_dict[chrome], chrome, indir) for chrome in file_dict)):
        outlines_main += result;

core.PWS("# " + core.getDateTime() + " | Writing output: " + outfilename);
with open(outfilename, "w") as outfile:
    headers = "window,aln.length,var.sites,perc.var.sites,informative.sites,perc.informative.sites";
    outfile.write(headers + "\n");
    for outline in outlines_main:
        outfile.write(outline);





