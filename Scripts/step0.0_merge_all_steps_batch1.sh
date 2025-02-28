#!/bin/bash
mainpath="/mnt/hdisk/backup4T/ESC_RNA_SL"

fastqc_dir="$mainpath/trimmed_data/quality_control/raw_data/"
trimmamotics_dir="$mainpath/trimmed_data/clean_data"
hisat2_dir="$mainpath/alignment_res/"
RseQC_dir="$mainpath/alignment_Rseqc/"
featurecounts_dir="$mainpath/featurecounts_gtf_res/"

multiqc_dir="$mainpath/MultiQC_res/"
mkdir -p $multiqc_dir
multiqc -o $multiqc_dir -d $fastqc_dir $trimmamotics_dir $hisat2_dir $RseQC_dir $featurecounts_dir
