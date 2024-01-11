#!/bin/bash
#SBATCH --job-name=WES_combine
#SBATCH --partition compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#########################################################################
##################### Parameters to be changed ##########################
#########################################################################
temdic=/share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/temp_dic/ # Directory of temp data
vcf_output=/share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/gvcf_output/
combinded_vcf=/share/home/xudeshu/scanpy_dic/HSCR/WES_out/hg38_out/gvcf_output/
#########################################################################
##################### Parameters should not be changed ##################
#########################################################################
GATK4=/share/soft/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar
Referencefile=/share/soft/human_gatk_resource/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
GATKResourceBundle=/share/soft/human_gatk_resource/
interval=/share/soft/human_gatk_resource/hu38_exon.interval_list
java -Xmx20g  -Djava.io.tmpdir=${temdic} -jar ${GATK4} GenotypeGVCFs \
-R ${Referencefile} \
-V ${combinded_vcf}raw_output2.vcf.gz \
-G StandardAnnotation \
-O ${combinded_vcf}raw_output3.vcf.gz
#java -Xmx20g -jar ${GATK4} GenotypeGVCFs \
#-R ${Referencefile} \
#-V gendb://my_database \
#--tmp-dir ${temdic} \
#-O ${combinded_vcf}raw_output.vcf.gz
# Recalibrating SNP
java -Xmx20g -jar ${GATK4} VariantRecalibrator \
-R ${Referencefile} \
-V ${combinded_vcf}raw_output3.vcf.gz \
--resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${GATKResourceBundle}resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz \
--resource:omini,known=false,training=true,truth=false,prior=12.0 ${GATKResourceBundle}resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz \
--resource:1000G,known=false,training=true,truth=false,prior=10.0 ${GATKResourceBundle}resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz \
--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${GATKResourceBundle}resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
-mode SNP \
-O ${combinded_vcf}output1.recal \
--tranches-file ${combinded_vcf}output1.tranches \
--rscript-file ${combinded_vcf}output1.plots.R
java -Xmx20g -jar ${GATK4} ApplyVQSR \
-R ${Referencefile} \
-V ${combinded_vcf}raw_output3.vcf.gz \
-O ${combinded_vcf}filtered1.vcf.gz \
--truth-sensitivity-filter-level 99.0 \
--tranches-file ${combinded_vcf}output1.tranches \
--recal-file ${combinded_vcf}output1.recal \
-mode SNP
# Recalibrating INDEL
java -Xmx20g -Djava.io.tmpdir=${temdic} -jar ${GATK4} VariantRecalibrator \
-R ${Referencefile} \
-V ${combinded_vcf}filtered1.vcf.gz \
--resource:mills,known=true,training=true,truth=true,prior=12.0 ${GATKResourceBundle}resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
-mode INDEL \
-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
--max-gaussians 6 \
-O ${combinded_vcf}output2.recal \
--tranches-file ${combinded_vcf}output2.tranches \
--rscript-file ${combinded_vcf}output2.plots.R 
java -Xmx20g -Djava.io.tmpdir=${temdic} -jar ${GATK4} ApplyVQSR \
-R ${Referencefile} \
-V ${combinded_vcf}filtered1.vcf.gz \
--ts-filter-level 99.0 \
-mode INDEL \
--tranches-file ${combinded_vcf}output2.tranches \
--recal-file ${combinded_vcf}output2.recal \
-O ${combinded_vcf}filtered_finial.vcf.gz 
