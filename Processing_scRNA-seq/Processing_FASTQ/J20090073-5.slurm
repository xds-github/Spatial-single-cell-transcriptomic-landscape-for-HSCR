#!/bin/bash
#SBATCH --job-name=cellranger_J20090073-5
#SBATCH --partition compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
cd /share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/
/share/soft/cellranger-6.1.1/cellranger count --localmem=128 --localcores=16 --id=J20090073-5 --sample=J20090073-5 --transcriptome=/share/soft/refdata-gex-GRCh38-2020-A --fastqs=/data/xudeshu/HSCR_data/scRNA-seq/cohort_data/J20090073-5
source /share/home/xudeshu/.bashrc
conda activate scve
velocyto run10x -@ 16 /share/home/xudeshu/scanpy_dic/HSCR/cellranger_out/J20090073-5 /share/soft/refdata-gex-GRCh38-2020-A/genes/genes.gtf
