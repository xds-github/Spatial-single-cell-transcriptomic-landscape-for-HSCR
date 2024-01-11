#!/bin/bash
#SBATCH --job-name=WES_HS572_L1
#SBATCH --partition compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
/share/soft/samtools-1.13/bin/samtools merge /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic2/HS572_L1.sorted.markdup.bam /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic2/HS572_H73H5DSX2_L4.sorted.markdup.bam /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic2/HS572_HFGVVDSX2_L4.sorted.markdup.bam
/share/soft/samtools-1.13/bin/samtools index /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic2/HS572_L1.sorted.markdup.bam
#########################################################################
#####################       BQSR      ###################################
#########################################################################
java -Xmx20g -jar /share/soft/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar BaseRecalibrator  -R /share/soft/human_gatk_resource/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -I /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic2/HS572_L1.sorted.markdup.bam -L /share/soft/human_gatk_resource/hu38_exon.interval_list --known-sites /share/soft/human_gatk_resource/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz --known-sites /share/soft/human_gatk_resource/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites /share/soft/human_gatk_resource/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /share/soft/human_gatk_resource/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz --known-sites /share/soft/human_gatk_resource/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --tmp-dir /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic2/ -O /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/ready_bam/HS572.recal_data1.table
java -Xmx20g -jar /share/soft/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar ApplyBQSR  -R /share/soft/human_gatk_resource/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -I /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic2/HS572_L1.sorted.markdup.bam --bqsr-recal-file /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/ready_bam/HS572.recal_data1.table -L /share/soft/human_gatk_resource/hu38_exon.interval_list --tmp-dir /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic2/ -O /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/ready_bam/HS572.sorted.markdup.realign.BQSR.bam
/share/soft/samtools-1.13/bin/samtools index /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/ready_bam/HS572.sorted.markdup.realign.BQSR.bam
#########################################################################
##################### HaplotypeCaller ###################################
#########################################################################
java -Xmx20g -jar /share/soft/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar HaplotypeCaller -R /share/soft/human_gatk_resource/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -I /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/ready_bam/HS572.sorted.markdup.realign.BQSR.bam -O /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/gvcf_output/HS572.g.vcf.gz -ERC GVCF -L /share/soft/human_gatk_resource/hu38_exon.interval_list --tmp-dir /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic2/

