#!/bin/bash
#SBATCH --job-name=maf2vcf_HS589
#SBATCH --partition compute
#SBATCH --nodes=1
#SBATCH --nodelist=n01
#SBATCH --cpus-per-task=6
source /share/home/xudeshu/.bashrc
java -Xmx5g -jar /share/soft/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar SelectVariants --exclude-filtered true --exclude-non-variants true -R /share/soft/human_gatk_resource/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -V /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/gvcf_output/filtered_finial.vcf.gz --tmp-dir /share/home/xudeshu/tmp/ --sample-name HS589 -O /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/splited_vcf_v2/HS589.vcf
conda activate vep104
#mkdir /share/home/xudeshu/tmp/HS589
#cd /share/home/xudeshu/tmp/HS589
vep -i /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/splited_vcf_v2/HS589.vcf -o /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/splited_vcf_v2/HS589.vep.vcf --force_overwrite --offline --cache --dir /share/home/xudeshu/Bio_resource/vep_resource/ --everything --fasta /share/soft/human_gatk_resource/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta --species homo_sapiens --cache_version 104 --vcf --assembly GRCh38 --fork 4 
filter_vep -filter "gnomAD_EAS_AF < 0.01 and gnomAD_AF < 0.01 and AF < 0.01 and EAS_AF < 0.01" -i /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/splited_vcf_v2/HS589.vep.vcf --force_overwrite  -o /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/splited_vcf_v2/HS589.vep_filtered.vcf
perl /share/home/xudeshu/Bio_resource/vcf2maf-1.6.21/vcf2maf.pl --input-vcf /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/splited_vcf_v2/HS589.vep_filtered.vcf --output-maf /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/maf_dic_v2/HS589.vep_filtered.maf --tumor-id HS589 --inhibit-vep --max-subpop-af 0.01 --ncbi-build GRCh38 --ref-fasta /share/soft/human_gatk_resource/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
