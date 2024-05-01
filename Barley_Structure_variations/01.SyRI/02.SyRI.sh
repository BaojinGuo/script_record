ls *.sam|cut -f1 -d"."|while read line; do echo '#!/bin/bash
#SBATCH --job-name=syri
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate syri
srun --export=all -n 1 -c 128 syri -c '"$line"'.Morex.sam -r Morex.V3.chr.fasta -q '"$line"'.V1.fasta -k -F S --dir 01.result/'"$line"' --prefix '"$line"'_Morex --nc 7' >$line.syri.sh ; done
