#!/bin/bash
#SBATCH -A MST109178
#SBATCH -J DEseq2
#SBATCH -p ngs186G
#SBATCH -c 28
#SBATCH --mem=186g
#SBATCH -o Deseq2_out.txt
#SBATCH -e Deseq2_err.txt

./1_summaryCounts.R \
 -c ../counts/ \
 -o ../meta/raw_count_table

./2_DESeq2.R \
 -i ../meta/raw_count_table.csv \
 -a ../meta/Sample_Annotation.csv \
 -O ../result/ \
 -j '~ Annotation' -g 'Annotation' \
 -t NAT10_KO_RNA_seq -c NAT10_control_RNA_seq

