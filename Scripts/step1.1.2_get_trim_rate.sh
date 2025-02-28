#!/bin/bash
mainpath="/mnt/hdisk/backup4T/ESC_RNA_SL2"
for line in `cat $mainpath/sample_infor/all_sample_list.txt`
do
	date;
	indir="$mainpath/trimmed_data/clean_data/$line"
	outdir="$mainpath/trimmed_data/trim_rate"
	mkdir -p $outdir;
	if [ ! -f "$indir/trimmomatics_summary.txt" ];then
	echo "trimmomatics_summary.txt for sample $line does not exists!";
	else
		echo "summarizing alignment result for sample $line";
		header="total_reads\tclean_reads\ttrim_rate";
		trimrate=`cat $indir/trimmomatics_summary.txt | grep '^Input Read Pairs:' | sed 's/^Input Read Pairs: \([0-9]\+\) Both Surviving: \([0-9]\+\) (\(.*\)) Forward Only Surviving.*/\1\t\2\t\3/g'`;
		echo -e "$header\n$trimrate\n" > $outdir/$line"_trim_summary.tab"
	fi
done

