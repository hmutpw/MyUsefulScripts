#!/bin/bash

for line in `cat /mnt/hdisk/backup4T/ESC_RNA_SL/sample_infor/all_sample_list.txt`
do
	echo "Runing stringtie with $line"
	indir="/mnt/hdisk/backup4T/ESC_RNA_SL/alignment_res/$line"
	outdir="/mnt/hdisk/backup4T/ESC_RNA_SL/stringtie_res/nascent_loc/$line"
	mkdir -p $outdir;
	stringtie $indir/accepted_hits_sorted.bam -o $outdir/outRes.gtf -p 8 -G /mnt/hdisk/backup4T/ESC_RNA_SL/ref_genome/gencode.v28.annotation.gtf --rf -A $outdir/gene_abund.tab -B -e 

done

