ls 01.minimap/*.bam|cut -f2 -d "/"|while read line; do filename=$(basename "$line"); prefix=${filename%.Morex.sort.bam}; echo '#!/bin/bash
#SBATCH --job-name=pbsv
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate pbsv
srun --export=all -n 1 -c 128 pbsv discover 01.minimap/'"$line"' 02.result/'"$prefix"'.svsig.gz 
srun --export=all -n 1 -c 128 pbsv call Morex.V3.chr.fasta 02.result/'"$prefix"'.svsig.gz 02.result/'"$prefix"'.pbsv.vcf --min-N-in-gap 1000 --filter-near-reference-gap 0 --ccs' >$prefix.pbsv.sh 
done
