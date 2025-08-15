#First step fastp 
#!/bin/bash

mapping_file="mapping.txt"
raw_dir="./Rawdata"
clean_dir="./cleandata"
report_dir="./report"

mkdir -p "$clean_dir" "$report_dir"

while read y_id a_id; do
    echo '#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate fastp
srun --export=all -n 1 -c 128 fastp \
    -i '"$raw_dir/$y_id"'_1.fq.gz \
    -I '"$raw_dir/$y_id"'_2.fq.gz \
    -o '"$clean_dir/$a_id"'_clean_1.fq.gz \
    -O '"$clean_dir/$a_id"'_clean_2.fq.gz \
    -h '"$report_dir/$a_id"'.fastp.html \
    -j '"$report_dir/$a_id"'.fastp.json \
    -w 128' > "$a_id.fastp.slurm"
done < "$mapping_file"


##Second Step Filter reads
ls *.sam|cut -f1 -d"_"|while read line; do echo '#!/bin/bash
#SBATCH --job-name='${line}'_filter
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
module load sambamba/1.0--h98b6b92_0 
srun --export=all -n 1 -c 128 sambamba view -F "mapping_quality >= 30 and proper_pair and [NM] < 2" -S -h -p -t 128 -f bam -o /scratch/pawsey0399/the1/Murdoch/2.Oat/Raw_data/04.Filter-BAM/'${line}'_SANG.bwa.filter.bam '${line}'_SANG.bwa.sam ' >$line.filter.sh; done



####Third Step Sort filter reads
ls *.filter.bam|cut -f1 -d"_"|while read line; do echo '#!/bin/bash
#SBATCH --job-name='${line}'_sort
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
module load sambamba/1.0--h98b6b92_0 
srun --export=all -n 1 -c 16 samtools sort -@ 16 -T /scratch/pawsey0399/the1/TMP -O BAM -o '${line}'_SANG.filter.sort.bam '${line}'_SANG.bwa.filter.bam' >$line.sort.sh
done
















