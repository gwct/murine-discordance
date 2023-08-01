#!/usr/bin/bash
#
# Liftsover mm10 CDS coordinates for each genome
#

genome_dir="/home/gt156213e/genomes/cds/genomes/"
bed_dir="/home/gt156213e/genomes/cds/bed/"

ref="${genome_dir}/genomes/mm10.chromes.fa"
ref_bed="${bed_dir}/mm10.ensGene.chromes.longest.cds.bed"


#awk 'OFS="\t" {print $0, $1":"$2"-"$3}' ${ref_bed} > ${ref_bed_named}
# Add the name column to the mouse bed file

chain="${genome_dir}/gdol/iter-03-softmask-final.chain"
spec="gdol"
echo ${spec}
liftOver ${ref_bed} ${chain} ${bed_dir}/${spec}-cds.bed ${bed_dir}/${spec}-cds.bed.unmap

chain="${genome_dir}/hall/iter-03-softmask-final.chain"
spec="hall"
echo ${spec}
liftOver ${ref_bed} ${chain} ${bed_dir}/${spec}-cds.bed ${bed_dir}/${spec}-cds.bed.unmap

chain="${genome_dir}/mnat/iter-03-softmask-final.chain"
spec="mnat"
echo ${spec}
liftOver ${ref_bed} ${chain} ${bed_dir}/${spec}-cds.bed ${bed_dir}/${spec}-cds.bed.unmap

chain="${genome_dir}/pdel/iter-03-softmask-final.chain"
spec="pdel"
echo ${spec}
liftOver ${ref_bed} ${chain} ${bed_dir}/${spec}-cds.bed ${bed_dir}/${spec}-cds.bed.unmap

chain="${genome_dir}/rdil/iter-03-softmask-final.chain"
spec="rdil"
echo ${spec}
liftOver ${ref_bed} ${chain} ${bed_dir}/${spec}-cds.bed ${bed_dir}/${spec}-cds.bed.unmap

chain="${genome_dir}/rsor/iter-03-softmask-final.chain"
spec="rsor"
echo ${spec}
liftOver ${ref_bed} ${chain} ${bed_dir}/${spec}-cds.bed ${bed_dir}/${spec}-cds.bed.unmap

