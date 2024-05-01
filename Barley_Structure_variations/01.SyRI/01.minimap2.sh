ls *.fasta|grep -v Morex|while read line; do filename=$(basename "$line");prefix=${filename%.V1.fasta}; echo '#!/bin/bash
#SBATCH --job-name=Mini
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
module load singularity/3.11.4-nompi
module load samtools/1.15--h3843a85_0
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/Singularity_image/minimap2.sif minimap2 -ax asm5 -I100g -f100 -t 128 --eqx Morex.V3.chr.fasta '"$line"' > 01.results/'"$prefix"'.Morex.sam' >$prefix.mini.sh
done
