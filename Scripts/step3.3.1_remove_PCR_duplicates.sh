#!/bin/bash
mainpath="/mnt/hdisk/backup4T/ESC_RNA_SL"
for line in `cat $mainpath/sample_infor/all_sample_list.txt`
do
	date;
	echo "[gatk MarkDuplicates] $line"
	outdir="$mainpath/alignment_res/$line"
	mkdir -p $outdir;
	gatk MarkDuplicates --REMOVE_DUPLICATES true -I $outdir/accepted_hits_sorted.bam -M $outdir/duplicate_metrc.txt -O $outdir/accepted_hits_rm_dup.bam
	samtools sort -@ 8 $outdir/accepted_hits_rm_dup.bam -o $outdir/accepted_hits_sorted_rm_dup.bam
	samtools index -@ 8 $outdir/accepted_hits_sorted_rm_dup.bam
	rm $outdir/accepted_hits_rm_dup.bam
done
