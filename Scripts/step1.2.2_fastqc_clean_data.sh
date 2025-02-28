#/bin/bash
mainpath="/mnt/hdisk/backup4T/ESC_RNA_SL2"
for line in `cat $mainpath/sample_infor/all_sample_list.txt`
do
	date;
	echo "processing $line"
	indir="$mainpath/trimmed_data/clean_data/$line"
	outdir="$mainpath/trimmed_data/quality_control/clean_data/$line"
	mkdir -p $outdir;
	fastqc -t 8 -o $outdir $indir/*.clean.fq.gz
done




