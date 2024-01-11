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
#########################################################################
##################### Create script #####################################
#########################################################################
# Combining gvcf
echo "#!/bin/bash
#SBATCH --job-name=WES_${sample_Name}_${Lane}
#SBATCH --partition compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
java -Xmx20g -jar ${GATK4} CombineGVCFs \
-R ${Referencefile} \
-L ${interval} \
--tmp-dir ${temdic} \
-O ${combinded_vcf}raw_output2.vcf.gz \\" > ${vcf_output}merge_gvcf.slurm
gvcf_file=$(ls ${vcf_output} | grep '.g.vcf.gz$')
for i in ${gvcf_file}
do
echo "--variant ${vcf_output}${i} \\" >> ${vcf_output}merge_gvcf.slurm
done
cat ${vcf_output}Combining_gvcf.slurm >> ${vcf_output}merge_gvcf.slurm
sbatch ${vcf_output}merge_gvcf.slurm
