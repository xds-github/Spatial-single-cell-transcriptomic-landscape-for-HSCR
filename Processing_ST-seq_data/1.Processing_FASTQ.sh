cat ST_sample_info.txt | while read id
do
arr=(${id})
sample_Name=${arr[0]}
pic=${arr[4]}
slide_ID=${arr[2]}
area_ID=${arr[3]}
echo "#!/bin/bash
#SBATCH --job-name=space_${sample_Name}
#SBATCH --partition compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
cd /ST_out/
spaceranger-1.3.0/spaceranger count \
--id=${sample_Name} \
--transcriptome=refdata-gex-GRCh38-2020-A \
--fastqs=/data/ \
--sample=${sample_Name} \
--image=/slide_meta/${pic} \
--slide=${slide_ID} \
--slidefile=/slide_meta/ \
--area=${area_ID} \
--localcores=16 \
--localmem=128" > /sh_folder/${sample_Name}.slurm
sbatch /sh_folder/${sample_Name}.slurm
done
