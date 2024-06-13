####this is for SNV pipeline
###samples:RGT_Planet(cultivar, UK) Igri(cultivar, Germany) OUN333(landrace, Nepal) Vlaminagh(cultivar, Austrilian) Hockett(cultivar, USA), 40X paired-end resequencing data
##first step, run pipeline of each SNV callers (DeepVariants GATK bcftools-mpileup)

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

##BWA align
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

####GATK-including Markduplicate and index, in fact, due to the big reference, and time and memory limited, it will be better to split chromosome to run the pipeline using parameter '-L', using module MergeVcfs to gather split-chromosome-vcf files in the end
ls 01.bwa/*.bam|cut -f2 -d "/"|cut -f1 -d "_"|while read line; do echo '#!/bin/bash
#SBATCH --job-name=gatk
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
module load gatk4/4.2.5.0--hdfd78af_0
module load module load samtools/1.15--h3843a85_0                                                                                                                                           
srun --export=all -n 1 -c 128 gatk --java-options "-XX:ParallelGCThreads=128" MarkDuplicates -I 01.bwa/'$line'_MorexV3.bwa.sort.bam -O 2.GATK/'$line'_MorexV3.bwa.sort.dup.bam -M 2.GATK/'$line'_MorexV3.bwa.sort.dup.metirc --REMOVE_DUPLICATES true                                                                                                                                   
srun --export=all -n 1 -c 128 samtools index -@ 128 -c 2.GATK/'$line'_MorexV3.bwa.sort.dup.bam
srun --export=all -n 1 -c 128 gatk --java-options "-XX:ParallelGCThreads=128" HaplotypeCaller -OVI False -R MorexV3.fa -I 2.GATK/'$line'_MorexV3.bwa.sort.dup.bam -O 2.GATK/'$line'_MorexV3.bwa.vcf --native-pair-hmm-threads 128
' >$line.gatk.sh; done


####DeepVariant
ls 01.bwa/*.bam|cut -f2 -d "/"|cut -f1 -d "_"|while read line; do  echo '#!/bin/bash
#SBATCH --job-name=Deepvariant
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=96:00:00
#SBATCH --account=pawsey0399
module load singularity/4.1.0-nompi
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/Singularity_image/deepvariant_latest.sif run_deepvariant --model_type WGS --ref MorexV3.fa --reads 01.bwa/'$line'_MorexV3.bwa.sort.bam --sample_name '$line' --output_vcf 1.deepvariant/'$line'.bwa.deep.vcf.gz --num_shards 128 --logging_dir 1.deepvariant/'$line'.deep.log
' >$line.deep.sh; done


ls 01.bwa/*.bam|cut -f2 -d "/"|cut -f1 -d "_"|while read line; do  echo '#!/bin/bash
#SBATCH --job-name=bcfcall
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
module load bcftools/1.15--haf5b3da_0
srun --export=all -n 1 -c 128 bcftools mpileup -f MorexV3.fa -Q 20 -q 20 -C 50 -Ou 2.GATK/'$line'_MorexV3.bwa.sort.dup.bam | bcftools call -mv -Oz -o 3.bcftools/'$line'.variants.vcf.gz --threads 128
'>$line.bcfcall.sh; done










