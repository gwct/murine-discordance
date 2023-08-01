#!/usr/bin/bash
#
# Gets window coordinates from the mm10 reference genome
#

ref="/projects/gt156213e/rodent-genomes/genomes/mm10.fa"
ref_rep="/mnt/beegfs/gt156213e/ref-genomes/mm10/mm10-RepeatMasker.gtf"
winsize="10000"
winsize_kb="10"
outbed="/mnt/beegfs/gt156213e/penn-genomes/windows/bed/mm10-${winsize_kb}kb-initial.bed"
outbed_win="/mnt/beegfs/gt156213e/penn-genomes/windows/bed/mm10-${winsize_kb}kb.bed"
outbed_rep="/mnt/beegfs/gt156213e/penn-genomes/windows/bed/mm10-${winsize_kb}kb-repeat-coverage.bed"

bedtools makewindows -g $ref.fai -w $winsize > $outbed
# Get the window coords

awk 'OFS="\t" {print $0, $1":"$2"-"$3}' $outbed > $outbed_win
# Add the window ID column

bedtools coverage -a $outbed_win -b $ref_rep > $outbed_rep
# Get the percent of each window that overlaps with repeat regions