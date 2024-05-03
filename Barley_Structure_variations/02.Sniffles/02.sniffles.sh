ls *.bam|while read line; do srun -n 1 -c 128 samtools index -c -@ 128 $line & 
done ###all bam files should be indexed.

ls 01.minimap/*.bam|cut -f2 -d "/"|while read line; do filename=$(basename "$line"); prefix=${filename%.Morex.sort.bam}; echo '#!/bin/bash
#SBATCH --job-name=sniffles
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate sniffle
srun --export=all -n 1 -c 128 sniffles --input 01.minimap/'"$line"' --vcf 02.sniffles/'"$prefix"'.Morex.vcf --reference Morex.V3.chr.fasta --snf 02.sniffles/'"$prefix"'.Morex.snf --threads 128 --minsupport 1 --long-ins-length 100000000   --long-del-length 100000000' >$prefix.sniffle.sh 
done
