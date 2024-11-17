###first, build index
srun -A pawsey0399 -p highmem -c 8 -n 1 -t 24:00:00 hisat2-build CS2Z559P1.fa CS2Z559P1.fa
samtools faidx CS2Z559P1.fa


####second, align
ls *_1.clean.fq.gz|while read line; do base_name=$(basename "$line" "_1.clean.fq.gz");  echo '#!/bin/bash
#SBATCH --job-name=hisat2
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate hisat
srun --export=all -n 1 -c 32 hisat2 -x 0.Index/CS2Z559P1.fa -1 '${base_name}'_1.clean.fq.gz -2 '${base_name}'_2.clean.fq.gz -S '${base_name}'.sam'  > ${base_name}.hisat.sh;  done




