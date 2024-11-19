###first, build index
srun -A pawsey0399 -p highmem -c 8 -n 1 -t 24:00:00 hisat2-build CS2Z559P1.fa CS2Z559P1.fa
samtools faidx CS2Z559P1.fa


####second, align
ls *_1.clean.fq.gz|while read line; do base_name=$(basename "$line" "_1.clean.fq.gz");  echo '#!/bin/bash
#SBATCH --job-name=hisat2
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=12:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate hisat
srun --export=all -n 1 -c 64 hisat2 --dta -x 0.Index/CS2Z559P1.fa -1 '${base_name}'_1.clean.fq.gz -2 '${base_name}'_2.clean.fq.gz -S '${base_name}'.sam'  > ${base_name}.hisat.sh;  done

####third, sort to bam file
ls *.sam|cut -f1 -d"."|while read line; do  echo '#!/bin/bash
#SBATCH --job-name=sort
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --account=pawsey0399
module load samtools/1.15--h3843a85_0
srun --export=all -n 1 -c 16 samtools sort -@ 16 -o '${line}'.sort.bam '${line}'.sam' >$line.sort.sh ; done

###Forth,featurecounts coount reads

srun -A Pawsey0399 -c 64 -n 1 -p work featureCounts -T 64 -t mRNA -g Parent -p -a 0.Index/CS2.1_Ac.gff3 -o 5113_2305_CS2_Ac.counts *.sort.bam










####anyway,if you mainly focus on special genes, maybe needs unique reads
ls *.sam|cut -f1 -d"."|while read line; do srun -A pawsey0399 -c 32 -n 1 -p work -t 2:00:00 grep "NH:i:1" $line.sam|grep "YT:Z:CP" > unique_reads/$line.unique.sam  & done
ls *.sam|cut -f1 -d"."|while read line; do srun -A pawsey0399 -c 32 -n 1 -p work -t 2:00:00 samtools view -H $line.sort.bam >unique_reads/$line.header; done 
cat $line.header $line.unique.sam |samtools sort -@ 32 -o $line.unique.sort.bam








