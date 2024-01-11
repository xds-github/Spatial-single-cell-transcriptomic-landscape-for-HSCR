#!/bin/bash
#SBATCH --job-name=cellranger_J20100207-B
#SBATCH --partition compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
cd /share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/
/share/soft/cellranger-6.1.1/cellranger vdj --localmem=64 --localcores=8 --id=J20100207-B --sample=J20100207-B --reference=/share/soft/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0 --fastqs=/data/xudeshu/HSCR_data/scRNA-seq/cohort_data/J20100207-B
