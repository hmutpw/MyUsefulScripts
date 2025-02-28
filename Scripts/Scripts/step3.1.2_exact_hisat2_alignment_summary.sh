#!/bin/bash
mainpath="/mnt/hdisk/backup4T/ESC_RNA_SL2"
for line in `cat $mainpath/sample_infor/all_sample_list.txt`
do
	date;
	indir="$mainpath/hisat2_res/alignment_res/$line"
	outdir="$mainpath/hisat2_res/alignment_summary/$line"
	mkdir -p $outdir;
	if [ ! -f "$indir/alignment_summary.txt" ];then
	echo "alignment_summary.txt for sample $line does not exists!";
	else
		echo "summarizing alignment result for sample $line";
		header="item\tcount\trate";
		totalpaired=`cat $indir/alignment_summary.txt | grep 'were paired; of these:$' | sed 's/^[ \t]\+\([0-9]\+\) (\(.*\)) were paired; of these:.*/total_count\t\1\t\2/g'`;
		uniquemap=`cat $indir/alignment_summary.txt | grep 'aligned concordantly exactly 1 time$' | sed 's/^[ \t]\+\([0-9]\+\) (\(.*\)) aligned concordantly exactly 1 time.*/unique_mapped\t\1\t\2/g'`;
		multimap=`cat $indir/alignment_summary.txt | grep 'aligned concordantly >1 times$' | sed 's/^[ \t]\+\([0-9]\+\) (\(.*\)) aligned concordantly >1 times.*/multiple_mapped\t\1\t\2/g'`;
		echo -e "$header\n$totalpaired\n$uniquemap\n$multimap\n" > $outdir/$line"_alignment_summary.tab"
	fi
done

