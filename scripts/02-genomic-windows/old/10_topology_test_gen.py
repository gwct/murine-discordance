#!/usr/bin/python3
############################################################
# For Penn genomes, 08.2020
# This script generates the files needed to run the topology
# tests using IQtree
############################################################

import sys, os, core, re, argparse, csv

############################################################
# Options
parser = argparse.ArgumentParser(description="Topology tests with IQtree");
parser.add_argument("-w", dest="window_size", help="The size of the sliding window in kb.", default=False);
parser.add_argument("-n", dest="name", help="A short name for all files associated with this job.", default=False);
parser.add_argument("--overwrite", dest="overwrite", help="If the output directory already exists and you wish to overwrite it, set this option.", action="store_true", default=False);
# IO options
parser.add_argument("-part", dest="part", help="SLURM partition option.", default=False);
parser.add_argument("-nodes", dest="nodes", help="SLURM --nodes option.", type=int, default=1);
parser.add_argument("-tasks", dest="tasks", help="SLURM --ntasks option.", type=int, default=1);
parser.add_argument("-cpus", dest="cpus", help="SLURM --cpus-per-task option.", type=int, default=1);
parser.add_argument("-mem", dest="mem", help="SLURM --mem option.", type=int, default=0);
# SLURM options

args = parser.parse_args();

if not args.name:
    sys.exit(" * Error 1: Please provide a prefix for the job and submit files with -n.");

if not args.window_size:
    sys.exit(" * Error 2: Window size in kb (-w) must be defined.");
if int(args.window_size) < 1:
    sys.exit(" * Error 3: Window size in kb (-w) must be a positive integer.");
else:
    wsize = int(args.window_size) * 1000;
    wsize_str = str(wsize);
# Parse the window size

if not args.part:
    sys.exit( " * Error 4: -part must be defined as a valid node partition on your clutser.");
if args.nodes < 1:
    sys.exit( " * Error 5: -nodes must be a positive integer.");
if args.tasks < 1:
    sys.exit( " * Error 6: -tasks must be a positive integer.");
if args.tasks < 1:
    sys.exit( " * Error 7: -cpus must be a positive integer.");
if args.tasks < 1:
    sys.exit( " * Error 8: -mem must be a positive integer.");
# SLURM option error checking

##########################
# Files

prefix = args.window_size + "kb-0.5-0.5"
windows_file = "/mnt/beegfs/gt156213e/penn-genomes/windows/data-2/" + prefix + "-topo-counts.csv";
base_alndir = "/mnt/beegfs/gt156213e/penn-genomes/windows/aln-2/" + prefix + "/";
base_outdir = "/mnt/beegfs/gt156213e/penn-genomes/windows/tree-2/" + prefix + "/topology-tests/";

log_dir = os.path.abspath(os.path.join("logs", args.name + "-logs"));
job_file = os.path.abspath(os.path.join("jobs", args.name + ".sh"));
submit_file = os.path.join("submit", args.name + "_submit.sh");

if (os.path.isfile(job_file) or os.path.isfile(submit_file)) and not args.overwrite:
    sys.exit( " * Error 7: Output directory (-o) already exists! Explicity specify --overwrite to overwrite it.");
# Get the tree file name and then make sure --overwrite is set if it exists

if not os.path.isdir(log_dir):
    os.system("mkdir " + log_dir);
# Make the log directory

##########################
# Run-time info
with open(job_file, "w") as logfile:
    core.runTime("#!/bin/bash\n# Penn -- generate topology tests", logfile);
    #core.PWS("# Psun reference FASTA: " + sun_fa, logfile);
    core.PWS("# Windows file:       " + windows_file, logfile);
    core.PWS("# Align directory:    " + base_alndir, logfile);
    core.PWS("# Output directory:   " + base_outdir, logfile);
    core.PWS("# Log directory:      " + log_dir, logfile);
    core.PWS("# Job file:           " + job_file, logfile);
    core.PWS("# Submit file:        " + submit_file, logfile);
    core.PWS("# ----------------", logfile);
##########################

    windows = {};
    chromes = {};
    # Main dictionaries.

    core.PWS("# " + core.getDateTime() + " Reading topologies per chromosome...", logfile);
    first = True;
    #for line in open(windows_file):
    with open(windows_file) as csvfile:
        windows_reader = csv.reader(csvfile);
        for line in windows_reader:
            if line[0][0] == "#":
                continue;
            # Skip commented lines

            if first:
                first = False;
                continue;
            # Skip the header line

            #line = line.strip().split(",");
            # Parse line

            if line[6] == "FILTER" or line[7] == "FILTER":
                continue;
            # Skip filtered windows

            window, chrome, start, end = line[0], line[1], line[3], line[4];
            topo, topo_rank, topo_num = line[17], int(line[16]), line[14];
            topo = re.sub("<[\d]+>", "", topo);
            # Window and topo info

            chrome_outdir = os.path.join(base_outdir, chrome);
            if not os.path.isdir(chrome_outdir):
                os.system("mkdir " + chrome_outdir);
            # Create the chromosome outdir

            if chrome not in chromes:
                chromes[chrome] = { 1 : "", 2 : "", 3 : "", "outdir" : chrome_outdir };
            # Initialize chrome dict if it hasn't been created

            if topo_rank <= 3:
                chromes[chrome][topo_rank] = topo;
            # Add topo to chrome dict if it is top 3

            windows[window] = { "chr" : chrome, "start" : start, "end" : end, "topo" : topo, "rank" : topo_rank };
            # Add window info to dict

    for chrome in chromes:
        values = list(chromes[chrome].values());
        assert "" not in values, "\nDidn't get three topologies:\n" + chrome + "\n" + str(chromes[chrome]);
    # Make sure all 3 topologies per chromosome are read

    core.PWS("# Chromosomes read: " + str(len(chromes)), logfile);
    core.PWS("# Windows read:     " + str(len(windows)), logfile);       
    core.PWS("# ----------------", logfile);
    # Read topology info

    ##########################

    core.PWS("# " + core.getDateTime() + " Generating topology tests...", logfile);
    core.PWS("# BEGIN CMDS", logfile);
    for window in windows:

        chrome = windows[window]["chr"];
        outdir = os.path.join(chromes[chrome]["outdir"], window.replace(":", "-"));
        if not os.path.isdir(outdir):
            os.system("mkdir " + outdir);
        # Make output directory

        aln_file = os.path.join(base_alndir, chrome + "-trimal", window.replace(":", "-") + "-mafft-trimal.fa");
        assert os.path.isfile(aln_file), "\nCannot find alignment file: " + aln_file;
        # Get alignment file

        treefile = os.path.join(outdir, window.replace(":", "-") + "-topos.tre");
        with open(treefile, "w") as t:
            t.write(windows[window]["topo"] + "\n");
            for rank in chromes[chrome]:
                if rank == "outdir":
                    continue;
                if rank == windows[window]["rank"]:
                    t.write("\n");
                else:
                    t.write(chromes[chrome][rank] + "\n");
        # Write trees file for current window

        prefix = os.path.join(outdir, "topo-test");
        window_log = os.path.join(log_dir, window.replace(":", "-") + "-topo-test.log");
        iqtree_cmd = "iqtree -s " + aln_file + " --trees " + treefile + " --test 10000 --test-weight --test-au --sitelh --prefix " + prefix + " &> " + window_log;
        core.PWS(iqtree_cmd, logfile);

##########################
# Generating the submit script.

with open(submit_file, "w") as sfile:
    submit = '''#!/bin/bash
#SBATCH --job-name={name}
#SBATCH --output={name}-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gregg.thomas@umontana.edu
#SBATCH --partition={partition}
#SBATCH --nodes={nodes}
#SBATCH --ntasks={tasks}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}

parallel -j {tasks} < {output_file}'''

    sfile.write(submit.format(name=args.name, partition=args.part, nodes=args.nodes, tasks=args.tasks, cpus=args.cpus, mem=args.mem, output_file=job_file));

##########################
