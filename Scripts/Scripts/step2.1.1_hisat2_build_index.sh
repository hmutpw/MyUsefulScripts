#!/bin/bash
mainpath="/mnt/hdisk/backup4T/ESC_RNA_SL2"
indir="/mnt/hdisk/backup4T/ESC_RNA_SL/ref_genome"
outdir="/mnt/hdisk/backup4T/ESC_RNA_SL/index"

# Step1------Preparing splice_sites AND exon file from gtf
extract_splice_sites.py $indir/gencode.v28.annotation.gtf > $outdir/gencode.v28.annotation.hisat2.ss
extract_exons.py $indir/gencode.vM4.annotation.gtf > $outdir/gencode.v28.annotation.hisat2.exon

# Step2------Building hisat2 index before alignment
mkdir -p $outdir/hisat2Index/mouseGencodev28

hisat2-build -p 8 $indir/GRCh38.p12.genome.fa $outdir/hisat2Index/humanGencodev28/humanGencodeIndex

