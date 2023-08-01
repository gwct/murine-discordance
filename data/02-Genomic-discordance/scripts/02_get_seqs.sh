#!/usr/bin/bash
#
# Gets sequences in bed format for each window in each genome
# Includes awk program to count the number and fraction of N/ns in each sequence and add them as columns to the bed files
# The final column definitions for the files are:
# chromosome    start   end     window id   sequence    window length   number of N/ns  fraction of N/ns
#

winsize="10000"
winsize_kb="10"
ref_bed="/mnt/beegfs/gt156213e/penn-genomes/windows/pseudo/bed/mm10-${winsize_kb}kb.bed"
ref_bed_named="/mnt/beegfs/gt156213e/penn-genomes/windows/pseudo/bed/mm10-${winsize_kb}kb-windows.bed"
beddir="/mnt/beegfs/gt156213e/penn-genomes/windows/pseudo/bed"


genome="/mnt/beegfs/gt156213e/ref-genomes/mm10/mm10.fa"
spec="mm10"
echo ${spec}
bedtools getfasta -bedOut -fi ${genome} -bed ${beddir}/${spec}-${winsize_kb}kb-windows.bed | awk 'BEGIN{OFS="\t"} {l=$3-$2; x=$NF; n=gsub(/N|n/, "", $NF); print $1, $2, $3, $4, x, l, n, n/l}' > ${beddir}/${spec}-${winsize_kb}kb-windows-seqs.bed


genome="/mnt/beegfs/gt156213e/penn-genomes/pseudo-assemblies/gdol-full/iter-03/fa/iter-03-softmask-final.fa"
spec="gdol"
echo ${spec}
bedtools getfasta -bedOut -fi ${genome} -bed ${beddir}/${spec}-${winsize_kb}kb.bed | awk 'BEGIN{OFS="\t"} {l=$3-$2; x=$NF; n=gsub(/N|n/, "", $NF); print $1, $2, $3, $4, x, l, n, n/l}' > ${beddir}/${spec}-${winsize_kb}kb-seqs.bed


genome="/mnt/beegfs/gt156213e/penn-genomes/pseudo-assemblies/hall/iter-03/fa/iter-03-softmask-final.fa"
spec="hall"
echo ${spec}
bedtools getfasta -bedOut -fi ${genome} -bed ${beddir}/${spec}-${winsize_kb}kb.bed | awk 'BEGIN{OFS="\t"} {l=$3-$2; x=$NF; n=gsub(/N|n/, "", $NF); print $1, $2, $3, $4, x, l, n, n/l}' > ${beddir}/${spec}-${winsize_kb}kb-seqs.bed


genome="/mnt/beegfs/gt156213e/penn-genomes/pseudo-assemblies/mnat-full/iter-03/fa/iter-03-softmask-final.fa"
spec="mnat"
echo ${spec}
bedtools getfasta -bedOut -fi ${genome} -bed ${beddir}/${spec}-${winsize_kb}kb.bed | awk 'BEGIN{OFS="\t"} {l=$3-$2; x=$NF; n=gsub(/N|n/, "", $NF); print $1, $2, $3, $4, x, l, n, n/l}' > ${beddir}/${spec}-${winsize_kb}kb-seqs.bed


genome="/mnt/beegfs/gt156213e/penn-genomes/pseudo-assemblies/pdel-full/iter-03/fa/iter-03-softmask-final.fa"
spec="pdel"
echo ${spec}
bedtools getfasta -bedOut -fi ${genome} -bed ${beddir}/${spec}-${winsize_kb}kb.bed | awk 'BEGIN{OFS="\t"} {l=$3-$2; x=$NF; n=gsub(/N|n/, "", $NF); print $1, $2, $3, $4, x, l, n, n/l}' > ${beddir}/${spec}-${winsize_kb}kb-seqs.bed


genome="/mnt/beegfs/gt156213e/penn-genomes/pseudo-assemblies/rdil/iter-03/fa/iter-03-softmask-final.fa"
spec="rdil"
echo ${spec}
bedtools getfasta -bedOut -fi ${genome} -bed ${beddir}/${spec}-${winsize_kb}kb.bed | awk 'BEGIN{OFS="\t"} {l=$3-$2; x=$NF; n=gsub(/N|n/, "", $NF); print $1, $2, $3, $4, x, l, n, n/l}' > ${beddir}/${spec}-${winsize_kb}kb-seqs.bed


genome="/mnt/beegfs/gt156213e/penn-genomes/pseudo-assemblies/rsor-full/iter-03/fa/iter-03-softmask-final.fa"
spec="rsor"
echo ${spec}
bedtools getfasta -bedOut -fi ${genome} -bed ${beddir}/${spec}-${winsize_kb}kb.bed | awk 'BEGIN{OFS="\t"} {l=$3-$2; x=$NF; n=gsub(/N|n/, "", $NF); print $1, $2, $3, $4, x, l, n, n/l}' > ${beddir}/${spec}-${winsize_kb}kb-seqs.bed

