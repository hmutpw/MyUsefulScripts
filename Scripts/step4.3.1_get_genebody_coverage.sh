#!/bin/bash
mainpath="/mnt/hdisk/backup4T/ESC_RNA_SL2"
for line in `cat $mainpath/sample_infor/all_sample_list.txt`
do
	date;
	indir="$mainpath/hisat2_res/alignment_res/$line"
	outdir="$mainpath/hisat2_res/alignment_qc/$line"
	mkdir -p $outdir
	geneBody_coverage.py -i $indir/accepted_hits_sorted.bam -r $mainpath/ref_genome/gencode.v28.annotation.bed -o $outdir/genebody_coverage
done



