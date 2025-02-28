#/bin/bash
mainpath="/mnt/hdisk/backup4T/ESC_RNA_SL2"
for line in `cat $mainpath/sample_infor/all_sample_list.txt`
do
	date;
	echo "processing $line"
	indir="$mainpath/raw_data/X101SC19010362-Z01-F001-B1-40_20190314/RawData"
	outdir="$mainpath/trimmed_data/quality_control/raw_data/$line"
	mkdir -p $outdir;
	fastqc -t 8 -o $outdir $indir/$line_*
done




