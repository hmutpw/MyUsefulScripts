#/bin/bash
mainpath="/mnt/hdisk/backup4T/ESC_RNA_SL2"
for line in `cat $mainpath/sample_infor/all_sample_list.txt`
do
	date;
	echo "processing $line"
	indir="$mainpath/raw_data/X101SC19010362-Z01-F001-B1-40_20190314/RawData"
	outdir="$mainpath/trimmed_data/clean_data/$line"
	mkdir -p $outdir;

	java -jar /home/hmutpw/Applications/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 -phred33 $indir/$line"_1.fq.gz" $indir/$line"_2.fq.gz" $outdir/$line"_1.clean.fq.gz" $outdir/$line"_1.unpaired.fq.gz" $outdir/$line"_2.clean.fq.gz" $outdir/$line"_2.unpaired.fq.gz" ILLUMINACLIP:/home/hmutpw/Applications/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:25  SLIDINGWINDOW:4:15 MINLEN:150 &> $outdir/"trimmomatics_summary.txt"

done




