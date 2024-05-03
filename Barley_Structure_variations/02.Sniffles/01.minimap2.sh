ls 00.Hifi-data/*|cut -f2 -d "/"|while read line
do filename=$(basename "$line")
prefix=${filename%.fastq.gz}
echo '#!/bin/bash
#SBATCH --job-name=Mini
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
module load singularity/3.11.4-nompi
module load samtools/1.15--h3843a85_0
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/Singularity_image/minimap2.sif minimap2 -ax map-hifi -R ''"'@RG\tID:"$prefix"\tSM:"$prefix"'"' -I100g -t 128 -Y Morex.V3.chr.fasta 00.Hifi-data/'"$line"' |samtools sort -@ 128 -O BAM -o 01.minimap/'"$prefix"'.Morex.sort.bam' >$prefix.mini.sh
done
