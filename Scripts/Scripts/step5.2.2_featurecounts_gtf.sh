#!/bin/bash

for line in `cat /mnt/hdisk/backup4T/ESC_RNA_SL/sample_infor/all_sample_list.txt`
do
	echo "Runing featurecounts with $line"
	indir="/mnt/hdisk/backup4T/ESC_RNA_SL/alignment_res/$line"
	outdir="/mnt/hdisk/backup4T/ESC_RNA_SL/featurecounts_gtf_res/$line"
	mkdir -p $outdir;
	featureCounts -p -T 8 -B -C -s 2 --primary --ignoreDup -a /mnt/hdisk/backup4T/ESC_RNA_SL/ref_genome/gencode.v28.annotation.gtf -o $outdir/gene_exon_counts.txt $indir/accepted_hits_sorted.bam &> $outdir/"quantify.log"

done
