####For illumina reads
##contents;ls
Akshinriki_1.fq.gz  Akshinriki_2.fq.gz	Clipper_1.fq.gz  Clipper_2.fq.gz  Hockett_1.fq.gz  Hockett_2.fq.gz  Stirling_1.fq.gz  Stirling_2.fq.gz	Vlamingh_1.fq.gz  Vlamingh_2.fq.gz

##raw data to clean reads
 ls *_1.fq.gz | cut -f1 -d "_" | while read line; do   echo '#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate fastp
srun --export=all -n 1 -c 128 fastp -i '$line'_1.fq.gz -I '$line'_2.fq.gz -o '$line'_clean_1.fq.gz -O '$line'_clean_2.fq.gz -w 128 '>$line.fastp.sh; done


ls *_1.fq.gz | cut -f1 -d "_" | while read line; do
  echo '#!/bin/bash
#SBATCH --job-name=BWA
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
module load bwa/0.7.17--h7132678_9
module load samtools/1.15--h3843a85_0
srun --export=all -n 1 -c 128 bwa mem -t 128 -R "@RG\tID:'$line'\tPL:illumina\tLB:library\tSM:'$line'" MorexV3.fa '${line}'_clean_1.fq.gz '${line}'_clean_2.fq.gz | samtools sort -@ 128 -o 01.bwa/'${line}'_MorexV3.bwa.sort.bam' > ${line}.bwa.sh
done

####For Hifi reads
##ls
Akashinirki.fastq.gz  Clipper.fastq.gz	Hockett.fastq.gz  Stirling.fastq.gz  Vlamingh.fastq.gz

ls *.fastq.gz|while read line; do filename=$(basename "$line"); prefix=${filename%.fastq.gz}; echo '#!/bin/bash
#SBATCH --job-name=Mini
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
module load singularity/4.1.0-nompi
module load samtools/1.15--h3843a85_0
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/Singularity_image/minimap2.sif minimap2 -ax map-hifi -c --MD -R '"'@RG\tID:"$prefix"\tSM:"$prefix"'"' -I100g -t 128 -Y MorexV3.fa '"$line"' |samtools sort -@ 128 -O BAM -o 01.minimap/'"$prefix"'.Morex.sort.bam' >$prefix.mini.sh; done

###determine raw data depth
samtools index -c ${line}_MorexV3.bwa.sort.bam
samtools depth ${line}_MorexV3.bwa.sort.bam >${line}.depth
awk '{sum += $3; count += 1} END {if (count > 0) print sum / count}' ${line}.depth >${line}.depth.av

###according to depth, split raw fastq to 5X, 10X, 15X depth fastq, -s need to be same.
seqkit sample -p $partition -j 16 -s 11 Vlamingh_clean_1.fq.gz|gzip >1.split/Vlamingh_clean5X_1.fq.gz
seqkit sample -p $partition -j 16 -s 11 Vlamingh_clean_2.fq.gz|gzip >1.split/Vlamingh_clean5X_2.fq.gz











