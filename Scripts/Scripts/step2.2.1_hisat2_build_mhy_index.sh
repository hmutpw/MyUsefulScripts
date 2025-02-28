#!/bin/bash
indir="/mnt/hdisk/backup4T/ESC_RNA_SL/ref_genome"
outdir="/mnt/hdisk/backup4T/ESC_RNA_SL/mhy_index"

# Step1------Preparing splice_sites AND exon file from gtf
extract_splice_sites.py $indir/GCF_000313635.1_ASM31363v1_genomic.gff > $outdir/GCF_000313635.1_ASM31363v1_genomic.hisat2.ss
extract_exons.py $indir/GCF_000313635.1_ASM31363v1_genomic.gff > $outdir/GCF_000313635.1_ASM31363v1_genomic.hisat2.exon

# Step2------Building hisat2 index before alignment
mkdir -p $outdir/hisat2Index/mhyncbi

hisat2-build -p 8 $indir/GCF_000313635.1_ASM31363v1_genomic.fna $outdir/hisat2Index/mhyncbi/mhyncbIndex

