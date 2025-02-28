#!/bin/bash

for line in `cat /mnt/hdisk/backup4T/ESC_RNA_SL/sample_infor/all_sample_list.txt`
do
	echo "extract unmapped reads with $line"
	indir="/mnt/hdisk/backup4T/ESC_RNA_SL/alignment_res/$line"
	cd $indir
	samtools view -@ 8 -f 4 $indir/accepted_hits_sorted.bam > $indir/unmapped_reads.sam
	samtools view -@ 8 -H $indir/accepted_hits_sorted.bam > $indir/sam_header.txt
	cat $indir/sam_header.txt $indir/unmapped_reads.sam > $indir/unmapped.sam
	samtools view -@ 8 -bS $indir/unmapped.sam > $indir/unmapped.bam
	samtools sort -@ 8 -o $indir/unmapped_sorted.bam $indir/unmapped.bam
	rm $indir/sam_header.txt
	rm $indir/unmapped.sam
	rm $indir/unmapped.bam
done

