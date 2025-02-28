#!/bin/bash
ref_dir="/mnt/hdisk/backup4T/ESC_RNA_SL/ref_genome"
perl $ref_dir/gtf2bed.pl $ref_dir/gencode.v28.annotation.gtf > $ref_dir/gencode.v28.annotation.bed
