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
srun --export=all -n 1 -c 128 svim alignment 02.results/'"$prefix"' 01.minimap/'"$line"' Morex.V3.chr.fasta  --min_sv_size 50 --sample '"$prefix"' '>$prefix.svim.sh
done
##don't add the parameter --symbolic_alleles due to SURVIVOR does not support implicit ALT
