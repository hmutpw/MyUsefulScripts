#!/bin/bash
echo "extract trim rate"
bash step1.1.2_get_trim_rate.sh
echo "quality control.";
bash step1.1.2_get_trim_rate.sh
echo "alignment with hisat2";
bash step3.1.1_hisat2_alignment.sh
echo "extract_alignmrnt_infor";
bash step3.1.2_exact_hisat2_alignment_summary.sh
echo "get_read_distribution";
bash step4.2.1_get_read_distribution.sh
echo "get_genebody_coverage";
bash step4.3.1_get_genebody_coverage.sh
echo "featurecounts";
bash step5.2.1_featurecounts_saf.sh
