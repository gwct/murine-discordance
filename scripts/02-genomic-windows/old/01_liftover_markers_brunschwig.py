#!/usr/bin/python3
############################################################
# For Penn genomes, 05.2021
# This script moves the coordinates of the markers from mm9
# to mm10. For Brunschwig markers.
############################################################

import sys, os, core

markerfile = "/mnt/beegfs/gt156213e/penn-genomes/windows/data/recombination-markers/brunschwig-rates/Brunschwig_2012-TableS1-rates.csv";
chainfile = "/mnt/beegfs/gt156213e/ref-genomes/mm10/mm9ToMm10.over.chain";
mm9bedfile = "/mnt/beegfs/gt156213e/penn-genomes/windows/data/recombination-markers/brunschwig-rates/mm9-markers.bed";
mm10bedfile = "/mnt/beegfs/gt156213e/penn-genomes/windows/data/recombination-markers/brunschwig-rates/mm10-markers.bed";
mm10bedfile_un = "/mnt/beegfs/gt156213e/penn-genomes/windows/data/recombination-markers/brunschwig-rates/mm10-markers-unmapped.bed";
#new_markerfile = "/mnt/beegfs/gt156213e/penn-genomes/windows/data/recombination-markers/Revised_HSmap_SNPs-mm10.csv";
logfilename = "logs/brunschwig-liftover-markers.log";

print("# " + core.getDateTime() + " Reading markers and writing BED file...");
with open(mm9bedfile, "w") as bedfile:
    markers = {};
    first = True;
    for line in open(markerfile):
        line = line.strip().split(",");
        if first:
            headers = line;
            first = False;
            continue;

        chrome, build37_kb, rho = line;
        build37 = int(float(build37_kb) * 1000)

        bedout = "\t".join([chrome, str(build37), str(build37+1), ".", rho, "+"]);
        bedfile.write(bedout + "\n");
    print("Markers read: " + str(len(markers)));
print("# ----------------");

print("# " + core.getDateTime() + " Lifting over marker coordinates...");
cmd = "liftOver " + mm9bedfile + " " + chainfile + " " + mm10bedfile + " " + mm10bedfile_un;
print("Executing liftOver cmd: " + cmd);
os.system(cmd);
print("# ----------------");

# print("# " + core.getDateTime() + " Writing new markers file...");
# with open(new_markerfile, "w") as outfile:
#     headers += ['chr38', 'build38']
#     outfile.write(",".join(headers) + "\n");
#     markers_written = 0;
#     for line in open(mm10bedfile):
#         line = line.strip().split("\t");
#         new_chr, new_start, new_end, snpid = line;
#         new_chr = new_chr.replace("chr", "");

#         outline = [snpid, markers[snpid]['chr37'], markers[snpid]['coord37'], markers[snpid]['fem_cM'], markers[snpid]['mal_cM'], markers[snpid]['ave_cM'], new_chr, new_start];
#         outfile.write(",".join(outline) + "\n");
#         markers_written += 1;
# print("Markers written: " + str(markers_written));
print("Done!");    
print("# ----------------");