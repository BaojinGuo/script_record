##bam files is the same as SYRI's input bam.

ls 01.minimap/*.bam|cut -f2 -d"/"|while read line
do
filename=$(basename "$line"); prefix=${filename%.Morex.sort.bam}; echo '#!/bin/bash
#SBATCH --job-name=SVIM
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate svim
srun --export=all -n 1 -c 128 svim alignment 02.results/'"$prefix"' 01.minimap/'"$line"' MorexV3.fa  --minimum_depth 5 --min_sv_size 50 --sample '"$prefix"' '>$prefix.svim.sh
done
##don't add the parameter --symbolic_alleles due to SURVIVOR does not support implicit ALT


###rename the result vcf files
for folder in /scratch/pawsey0399/bguo1/Murdoch/03.Austrilian_Barley_assembly/02.SV_calling/03.SVIM/02.results/*
do folder_name=$(basename "$folder")
new_filename="${folder_name}.SVIM.vcf"
mv "$folder/variants.vcf" "$folder/$new_filename"
echo "Renamed $folder/variants.vcf to $folder/$new_filename" 
done
###then, merge all sample vcf. 
##Move all samples vcf file to one file '00.ALL_SVIM_vcf' and run SURVIVOR
SURVIVOR merge 00.All_SVIM_vcf/ 1000 1 1 1 0 50 merge_SVIM.vcf
