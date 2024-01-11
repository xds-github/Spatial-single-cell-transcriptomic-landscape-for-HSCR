#!/bin/bash
#SBATCH --job-name=annovar_HD30
#SBATCH --partition compute
#SBATCH --nodes=1
#SBATCH --nodelist=n01
#SBATCH --cpus-per-task=2
source /share/home/xudeshu/.bashrc
perl /share/home/xudeshu/Bio_resource/annovar_analysis/annovar/table_annovar.pl /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/splited_vcf_v2/HD30.vcf /share/home/xudeshu/Bio_resource/annovar_analysis/annovar/humandb/ -buildver hg38 -out /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/annovar_out/HD30 -remove -protocol refGene,cytoBand,avsnp150,gnomad312_genome,dbnsfp42c,intervar_20180118,clinvar_20221231 -operation gx,r,f,f,f,f,f -nastring . -vcfinput -polish
