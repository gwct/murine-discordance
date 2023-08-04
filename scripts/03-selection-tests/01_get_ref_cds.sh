#/bin/bash

samtools faidx mm10.fa chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX > mm10.chromes.fa
samtools faidx mm10.chromes.fa
# Get only the full chromosomes from the mouse reference

grep ">" mm10.chromes.fa
grep -c ">" mm10.chromes.fa
# Checks

##########

zcat mm10.ensGene.gtf.gz | awk '/^#/||$1=="chr1"||$1=="chr2"||$1=="chr3"||$1=="chr4"||$1=="chr5"||$1=="chr6"||$1=="chr7"||$1=="chr8"||$1=="chr9"||$1=="chr10"||$1=="chr11"||$1=="chr12"||$1=="chr13"||$1=="chr14"||$1=="chr15"||$1=="chr16"||$1=="chr17"||$1=="chr18"||$1=="chr19"||$1=="chrX"{print}' > mm10.ensGene.chromes.gtf

#zcat Mus_musculus.GRCm38.99.gtf.gz | awk '/^#/||$1=="1"||$1=="2"||$1=="3"||$1=="4"||$1=="5"||$1=="6"||$1=="7"||$1=="8"||$1=="9"||$1=="10"||$1=="11"||$1=="12"||$1=="13"||$1=="14"||$1=="15"||$1=="16"||$1=="17"||$1=="18"||$1=="19"||$1=="X"{print}' | sed 's/^/^chr/g' > Mus_musculus.GRCm38.99.chromes.gtf
# Manually remove chr from header

# Get only the full chromosomes (and the header) from the mouse reference annotation

awk '{print $1}' mm10.ensGene.chromes.gtf | sort | uniq
awk '{print $1}' mm10.ensGene.chromes.gtf | sort | uniq | c
# Checks

##########

time -p python degenotate/degenotate.py -a ../annotations/mm10.ensGene.chromes.gtf -g ../genomes/mm10.chromes.fa -o ../annotations/mm10-degenotate --overwrite

#time -p degenotate.py -a Mus_musculus.GRCm38.99.chromes.gtf -g ../../00-mm10/mm10.chromes.fa -o mm10-degenotate-ensembl
# degenotate just to get the longest transcript ids

awk '$5=="yes"{print}' transcript-counts.tsv | cut -f 1 > longest-transcripts.txt
# Get only the longest transcript ids from the degenotate output

##########

gtfToGenePred mm10.ensGene.chromes.longest.gtf mm10.ensGene.chromes.longest.gp
genePredToBed mm10.ensGene.chromes.longest.gp ../bed/mm10.ensGene.chromes.longest.bed
cd ../bed/
bedparse cds mm10.ensGene.chromes.longest.bed > mm10.ensGene.chromes.longest.cds.bed
# Get reference gene annotations in bed format

##########

gxf_parse.py -i Mus_musculus.GRCm38.99.chromes.gtf -s ensembl -p ENSMUS --notranscripts --noexons -o mm10-genes.tab
tab_to_bed mm10-genes.tab > ../../../summary-data/02-genomic-windows/feature-beds/mm10-genes.bed
# Get all genes -- had to use a different gtf because the original one didn't have gene coords...

awk '$3=="CDS"{print}' mm10.ensGene.chromes.longest.gtf | cut -d \" -f 2 | sort | uniq > mm10-coding-genes.txt
grep -f mm10-coding-genes.txt ../../../summary-data/02-genomic-windows/feature-beds/mm10-genes.bed > ../../../summary-data/02-genomic-windows/feature-beds/mm10-coding-genes.bed
# Get coding genes

grep -f ../../../summary-data/02-genomic-windows/feature-beds/ps-transcripts-all-gt.txt mm10.ensGene.chromes.longest.gtf | cut -d \" -f 2 | sort | uniq > ps-genes-all-gt.txt
grep -f ps-genes-all-gt.txt ../../../summary-data/02-genomic-windows/feature-beds/mm10-genes.bed > ../../../summary-data/02-genomic-windows/feature-beds/mm10-ps-genes-all-gt.bed
# Get PS genes