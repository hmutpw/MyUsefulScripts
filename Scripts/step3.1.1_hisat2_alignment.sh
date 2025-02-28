#!/bin/bash
mainpath="/mnt/hdisk/backup4T/ESC_RNA_SL"
for line in `cat $mainpath/sample_infor/all_sample_list.txt`
do
	date;
	echo "Runing hisat2 with $line"
	indir="$mainpath/trimmed_data/clean_data/$line"
	outdir="$mainpath/hisat2_res/alignment_res/$line"
	mkdir -p $outdir;
	hisat2 --known-splicesite-infile $mainpath/index/gencode.v28.annotation.hisat2.ss --dta -t -p 2 --rna-strandness RF -x $mainpath/index/hisat2Index/humanGencodev28/humanGencodeIndex -1 $indir/$line"_1.clean.fq.gz" -2 $indir/$line"_2.clean.fq.gz" -S $outdir/accepted_hits.sam
	samtools view -@ 2 -bS $outdir/accepted_hits.sam > $outdir/accepted_hits.bam
	gatk MarkDuplicates --REMOVE_DUPLICATES true -I $outdir/accepted_hits.bam -M $outdir/duplicate_metrc.txt -O $outdir/accepted_hits_rm_dup.bam
	samtools sort -@ 2 $outdir/accepted_hits_rm_dup.bam -o $outdir/accepted_hits_sorted_rm_dup.bam
	samtools index -@ 2 $outdir/accepted_hits_sorted_rm_dup.bam
	rm $outdir/accepted_hits.sam
	rm $outdir/accepted_hits_rm_dup.bam
done
