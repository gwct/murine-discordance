#!/usr/bin/bash
#
# Get the CDS sequences for each genome in tab format
#

genome_dir="/home/gt156213e/genomes/cds/genomes/"
bed_dir="/home/gt156213e/genomes/cds/bed/"

spec="mm10"
echo ${spec}
bedtools getfasta -fi ${genome_dir}/mm10.chromes.fa -fo ${bed_dir}/${spec}-cds-seqs.tab -bed ${bed_dir}/mm10.ensGene.chromes.longest.cds.bed -name -split -tab -s

spec="gdol"
echo ${spec}
bedtools getfasta -fi ${genome_dir}/${spec}/iter-03-softmask-final.fa -fo ${bed_dir}/${spec}-cds-seqs.tab -bed ${bed_dir}/${spec}-cds.bed -name -split -tab -s

spec="hall"
echo ${spec}
bedtools getfasta -fi ${genome_dir}/${spec}/iter-03-softmask-final.fa -fo ${bed_dir}/${spec}-cds-seqs.tab -bed ${bed_dir}/${spec}-cds.bed -name -split -tab -s

spec="mnat"
echo ${spec}
bedtools getfasta -fi ${genome_dir}/${spec}/iter-03-softmask-final.fa -fo ${bed_dir}/${spec}-cds-seqs.tab -bed ${bed_dir}/${spec}-cds.bed -name -split -tab -s

spec="pdel"
echo ${spec}
bedtools getfasta -fi ${genome_dir}/${spec}/iter-03-softmask-final.fa -fo ${bed_dir}/${spec}-cds-seqs.tab -bed ${bed_dir}/${spec}-cds.bed -name -split -tab -s

spec="rdil"
echo ${spec}
bedtools getfasta -fi ${genome_dir}/${spec}/iter-03-softmask-final.fa -fo ${bed_dir}/${spec}-cds-seqs.tab -bed ${bed_dir}/${spec}-cds.bed -name -split -tab -s

spec="rsor"
echo ${spec}
bedtools getfasta -fi ${genome_dir}/${spec}/iter-03-softmask-final.fa -fo ${bed_dir}/${spec}-cds-seqs.tab -bed ${bed_dir}/${spec}-cds.bed -name -split -tab -s