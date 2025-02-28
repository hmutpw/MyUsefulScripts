#!/bin/bash
mainpath="/mnt/hdisk/backup4T/ESC_RNA_SL2"

fastqc_dir="$mainpath/trimmed_data/quality_control/raw_data/"
trimmamotics_dir="$mainpath/trimmed_data/clean_data/"
hisat2_dir="$mainpath/hisat2_res/alignment_res/"
RseQC_dir="$mainpath/hisat2_res/alignment_qc/"
featurecounts_dir="$mainpath/hisat2_res/featurecounts_res/"

multiqc_dir="$mainpath/MultiQC_res/"
mkdir -p $multiqc_dir
multiqc -o $multiqc_dir $fastqc_dir $trimmamotics_dir $hisat2_dir $RseQC_dir $featurecounts_dir
