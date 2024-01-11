#!/bin/bash
#SBATCH --job-name=WES_HS572_H73H5DSX2_L4
#SBATCH --partition compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
/share/soft/fastp/fastp -i /data/xudeshu/HSCR_data/WES_data/HS572/HS572_FKDN210187448-1A_H73H5DSX2_L4_1.fq.gz -o /share/home/xudeshu/scanpy_dic/HSCR/WES_out/filtered_data/HS572_FKDN210187448-1A_H73H5DSX2_L4_1.fq.gz -I /data/xudeshu/HSCR_data/WES_data/HS572/HS572_FKDN210187448-1A_H73H5DSX2_L4_2.fq.gz -O /share/home/xudeshu/scanpy_dic/HSCR/WES_out/filtered_data/HS572_FKDN210187448-1A_H73H5DSX2_L4_2.fq.gz -j /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/data_QCresult/HS572_H73H5DSX2_L4_fastp.json -w 4 -q 19 -u 50 -l 36
#This pipline is writen for transforing fastq data into ready_sort_BAM
#xudeshu
#2021/7/11
#v1
#software: GATK4.2.0.0 /bwa 0.7.17-r1188 / samtools 1.7 /Picard
#java -Xmx20g -jar /share/soft/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar IndexFeatureFile -I resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf
#########################################################################
##################### Map to Reference###################################
#########################################################################
/share/soft/bwa-0.7.17/bwa mem -t 6 -R "@RG\tID:HS572\tPL:Illumina\tLB:H73H5DSX2_L4\tPU:HS572\tSM:HS572"  /share/soft/human_gatk_resource/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta  /share/home/xudeshu/scanpy_dic/HSCR/WES_out/filtered_data/HS572_FKDN210187448-1A_H73H5DSX2_L4_1.fq.gz  /share/home/xudeshu/scanpy_dic/HSCR/WES_out/filtered_data/HS572_FKDN210187448-1A_H73H5DSX2_L4_2.fq.gz  | /share/soft/samtools-1.13/bin/samtools view -S -b > /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic2/HS572_H73H5DSX2_L4.bam
#########################################################################
##################### Mark Duplication###################################
#########################################################################
java -Xmx20g -Djava.io.tmpdir=/share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic2/ -jar /share/soft/picard/picard.jar SortSam I=/share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic2/HS572_H73H5DSX2_L4.bam O=/share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic2/HS572_H73H5DSX2_L4.sorted.bam SORT_ORDER=coordinate
java -Xmx20g -Djava.io.tmpdir=/share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic2/ -jar /share/soft/picard/picard.jar MarkDuplicates I=/share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic2/HS572_H73H5DSX2_L4.sorted.bam M=/share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic2/HS572_H73H5DSX2_L4.marked_dup_metrics.txt O=/share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic2/HS572_H73H5DSX2_L4.sorted.markdup.bam 
/share/soft/samtools-1.13/bin/samtools index /share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic2/HS572_H73H5DSX2_L4.sorted.markdup.bam
