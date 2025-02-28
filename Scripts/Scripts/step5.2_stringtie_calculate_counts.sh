#!/bin/bash
indir="/mnt/hdisk/backup4T/ESC_RNA_SL/stringtie_res/nascent_loc"
outdir="/mnt/hdisk/backup4T/ESC_RNA_SL/stringtie_res/count_res/nascent_loc"
prepDE.py -i $indir -g $outdir/gene_count_matrix.csv -t $outdir/transcript_count_matrix.csv





