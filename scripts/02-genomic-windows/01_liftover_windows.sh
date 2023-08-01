#!/usr/bin/bash
#
# Liftsover mm10 window coordinates for each genome
#

ref="/mnt/beegfs/gt156213e/ref-genomes/mm10/mm10.fa"
ref_rep="/mnt/beegfs/gt156213e/ref-genomes/mm10/mm10-RepeatMasker.gtf"
winsize="10000"
winsize_kb="10"
ref_bed="/mnt/beegfs/gt156213e/penn-genomes/windows/pseudo/bed/mm10-${winsize_kb}kb.bed"
ref_bed_named="/mnt/beegfs/gt156213e/penn-genomes/windows/pseudo/bed/mm10-${winsize_kb}kb-windows.bed"
outdir="/mnt/beegfs/gt156213e/penn-genomes/windows/pseudo/bed"

awk 'OFS="\t" {print $0, $1":"$2"-"$3}' ${ref_bed} > ${ref_bed_named}
# Add the name column to the mouse bed file

chain="/mnt/beegfs/gt156213e/penn-genomes/pseudo-assemblies/gdol-full/iter-03/fa/iter-03-softmask-final.chain"
spec="gdol"
echo ${spec}
liftOver ${ref_bed_named} ${chain} ${outdir}/${spec}-${winsize_kb}kb.bed ${outdir}/${spec}-${winsize_kb}kb.bed.unmap

chain="/mnt/beegfs/gt156213e/penn-genomes/pseudo-assemblies/hall/iter-03/fa/iter-03-softmask-final.chain"
spec="hall"
echo ${spec}
liftOver ${ref_bed_named} ${chain} ${outdir}/${spec}-${winsize_kb}kb.bed ${outdir}/${spec}-${winsize_kb}kb.bed.unmap

chain="/mnt/beegfs/gt156213e/penn-genomes/pseudo-assemblies/mnat-full/iter-03/fa/iter-03-softmask-final.chain"
spec="mnat"
echo ${spec}
liftOver ${ref_bed_named} ${chain} ${outdir}/${spec}-${winsize_kb}kb.bed ${outdir}/${spec}-${winsize_kb}kb.bed.unmap

chain="/mnt/beegfs/gt156213e/penn-genomes/pseudo-assemblies/pdel-full/iter-03/fa/iter-03-softmask-final.chain"
spec="pdel"
echo ${spec}
liftOver ${ref_bed_named} ${chain} ${outdir}/${spec}-${winsize_kb}kb.bed ${outdir}/${spec}-${winsize_kb}kb.bed.unmap

chain="/mnt/beegfs/gt156213e/penn-genomes/pseudo-assemblies/rdil/iter-03/fa/iter-03-softmask-final.chain"
spec="rdil"
echo ${spec}
liftOver ${ref_bed_named} ${chain} ${outdir}/${spec}-${winsize_kb}kb.bed ${outdir}/${spec}-${winsize_kb}kb.bed.unmap

chain="/mnt/beegfs/gt156213e/penn-genomes/pseudo-assemblies/rsor-full/iter-03/fa/iter-03-softmask-final.chain"
spec="rsor"
echo ${spec}
liftOver ${ref_bed_named} ${chain} ${outdir}/${spec}-${winsize_kb}kb.bed ${outdir}/${spec}-${winsize_kb}kb.bed.unmap

