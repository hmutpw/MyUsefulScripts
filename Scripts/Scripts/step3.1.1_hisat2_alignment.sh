#!/bin/bash
mainpath="/mnt/hdisk/backup4T/ESC_RNA_SL2"
for line in `cat $mainpath/sample_infor/all_sample_list.txt`
do
	echo "Runing hisat2 with $line"
	indir="$mainpath/trimmed_data/clean_data/$line"
	outdir="$mainpath/hisat2_res/alignment_res/$line"
	mkdir -p $outdir;
	hisat2 --known-splicesite-infile $mainpath/index/gencode.v28.annotation.hisat2.ss --dta -t -p 8 --rna-strandness RF -x $mainpath/index/hisat2Index/humanGencodev28/humanGencodeIndex -1 $indir/$line"_1.clean.fq.gz" -2 $indir/$line"_2.clean.fq.gz" -S $outdir/accepted_hits.sam &> $outdir/alignment_summary.txt
	samtools view -@ 8 -F 4 -bS $outdir/accepted_hits.sam > $outdir/accepted_hits.bam
	samtools sort -@ 8 $outdir/accepted_hits.bam -o $outdir/accepted_hits_sorted.bam
	samtools index -@ 8 $outdir/accepted_hits_sorted.bam
	rm $outdir/accepted_hits.sam
	rm $outdir/accepted_hits.bam
done
