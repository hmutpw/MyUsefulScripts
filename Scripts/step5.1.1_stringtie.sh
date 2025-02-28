#!/bin/bash
mainpath="/mnt/hdisk/backup4T/ESC_RNA_SL2"
for line in `cat $mainpath/sample_infor/all_sample_list.txt`
do
	echo "Runing stringtie with $line"
	indir="$mainpath/hisat2_res/alignment_res/$line"
	outdir="$mainpath/hisat2_res/stringtie_res/nascent_loc/$line"
	mkdir -p $outdir;
	stringtie $indir/accepted_hits_sorted.bam -o $outdir/outRes.gtf -p 8 -G $mainpath/ref_genome/gencode.v28.annotation.gtf --rf -A $outdir/gene_abund.tab -B -e 

done

