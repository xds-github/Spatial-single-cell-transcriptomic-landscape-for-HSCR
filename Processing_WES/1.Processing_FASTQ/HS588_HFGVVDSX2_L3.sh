#!/bin/bash
#SBATCH --job-name=WES_HS588_HFGVVDSX2_L3
#SBATCH --partition compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
/share/soft/fastp/fastp -i /data/xudeshu/HSCR_data/WES_data/HS588/HS588_FKDN210187456-1A_HFGVVDSX2_L3_1.fq.gz -o /share/home/xudeshu/scanpy_dic/HSCR/WES_out/filtered_data/HS588_FKDN210187456-1A_HFGVVDSX2_L3_1.fq.gz -I /data/xudeshu/HSCR_data/WES_data/HS588/HS588_FKDN210187456-1A_HFGVVDSX2_L3_2.fq.gz -O /share/home/xudeshu/scanpy_dic/HSCR/WES_out/filtered_data/HS588_FKDN210187456-1A_HFGVVDSX2_L3_2.fq.gz -j /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/data_QCresult/HS588_HFGVVDSX2_L3_fastp.json -w 4 -q 19 -u 50 -l 36
#This pipline is writen for transforing fastq data into ready_sort_BAM
#xudeshu
#2021/7/11
#v1
#software: GATK4.2.0.0 /bwa 0.7.17-r1188 / samtools 1.7 /Picard
#java -Xmx20g -jar /share/soft/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar IndexFeatureFile -I resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf
#########################################################################
##################### Map to Reference###################################
#########################################################################
/share/soft/bwa-0.7.17/bwa mem -t 6 -R "@RG\tID:HS588\tPL:Illumina\tLB:HFGVVDSX2_L3\tPU:HS588\tSM:HS588"  /share/soft/human_gatk_resource/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta  /share/home/xudeshu/scanpy_dic/HSCR/WES_out/filtered_data/HS588_FKDN210187456-1A_HFGVVDSX2_L3_1.fq.gz  /share/home/xudeshu/scanpy_dic/HSCR/WES_out/filtered_data/HS588_FKDN210187456-1A_HFGVVDSX2_L3_2.fq.gz  | /share/soft/samtools-1.13/bin/samtools view -S -b > /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic/HS588_HFGVVDSX2_L3.bam
#########################################################################
##################### Mark Duplication###################################
#########################################################################
java -Xmx20g -Djava.io.tmpdir=/share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic/ -jar /share/soft/picard/picard.jar SortSam I=/share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic/HS588_HFGVVDSX2_L3.bam O=/share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic/HS588_HFGVVDSX2_L3.sorted.bam SORT_ORDER=coordinate
java -Xmx20g -Djava.io.tmpdir=/share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic/ -jar /share/soft/picard/picard.jar MarkDuplicates I=/share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic/HS588_HFGVVDSX2_L3.sorted.bam M=/share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic/HS588_HFGVVDSX2_L3.marked_dup_metrics.txt O=/share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic/HS588_HFGVVDSX2_L3.sorted.markdup.bam 
/share/soft/samtools-1.13/bin/samtools index /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic/HS588_HFGVVDSX2_L3.sorted.markdup.bam
#########################################################################
#####################       BQSR      ###################################
#########################################################################
java -Xmx20g -jar /share/soft/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar BaseRecalibrator  -R /share/soft/human_gatk_resource/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -I /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic/HS588_HFGVVDSX2_L3.sorted.markdup.bam -L /share/soft/human_gatk_resource/hu38_exon.interval_list --known-sites /share/soft/human_gatk_resource/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz --known-sites /share/soft/human_gatk_resource/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites /share/soft/human_gatk_resource/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /share/soft/human_gatk_resource/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz --known-sites /share/soft/human_gatk_resource/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --tmp-dir /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic/ -O /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/ready_bam/HS588.recal_data1.table
java -Xmx20g -jar /share/soft/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar ApplyBQSR -R /share/soft/human_gatk_resource/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -I /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic/HS588_HFGVVDSX2_L3.sorted.markdup.bam --bqsr-recal-file /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/ready_bam/HS588.recal_data1.table -L /share/soft/human_gatk_resource/hu38_exon.interval_list --tmp-dir /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic/ -O /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/ready_bam/HS588.sorted.markdup.realign.BQSR.bam
/share/soft/samtools-1.13/bin/samtools index /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/ready_bam/HS588.sorted.markdup.realign.BQSR.bam
#########################################################################
##################### HaplotypeCaller ###################################
#########################################################################
java -Xmx20g -jar /share/soft/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar HaplotypeCaller -R /share/soft/human_gatk_resource/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -I /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/ready_bam/HS588.sorted.markdup.realign.BQSR.bam -O /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/gvcf_output/HS588.g.vcf.gz -ERC GVCF -L /share/soft/human_gatk_resource/hu38_exon.interval_list --tmp-dir /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic/
