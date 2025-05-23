#!/bin/bash
mainpath="/mnt/hdisk/backup4T/ESC_RNA_SL2"
for line in `cat $mainpath/sample_infor/all_sample_list.txt`
do
	echo "Runing featurecounts with $line"
	indir="$mainpath/hisat2_res/alignment_res/$line"
	outdir="$mainpath/hisat2_res/featurecounts_res/$line"
	mkdir -p $outdir;
#------gene level.
	featureCounts -p -T 8 -B -C -s 2 --primary --ignoreDup -F SAF -a $mainpath/ref_genome/intron_exon/gencode.v28.gene.only.saf -o $outdir/gene_counts.txt $indir/accepted_hits_sorted.bam &> $outdir/"gene.quantify.log"
#------exon only level.
	featureCounts -p -T 8 -B -C -s 2 --primary --ignoreDup -F SAF -a $mainpath/ref_genome/intron_exon/gencode.v28.exon.only.saf -o $outdir/exon_only_counts.txt $indir/accepted_hits_sorted.bam &> $outdir/"exon.only.quantify.log"
#------intron only level.
	featureCounts -p -T 8 -B -C -s 2 --primary --ignoreDup -F SAF -a $mainpath/ref_genome/intron_exon/gencode.v28.intron.only.saf -o $outdir/intron_only_counts.txt $indir/accepted_hits_sorted.bam &> $outdir/"intron.only.quantify.log"
#------exon intron together.
	featureCounts -p -T 8 -B -C -s 2 --primary --ignoreDup -F SAF -a $mainpath/ref_genome/intron_exon/gencode.v28.exon.intron.saf -o $outdir/exon_intron_counts.txt $indir/accepted_hits_sorted.bam &> $outdir/"exon.intron.quantify.log"
done
