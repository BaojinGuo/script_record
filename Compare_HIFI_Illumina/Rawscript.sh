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


####GATK
ls 01.bwa/*.bam|cut -f2 -d "/"|cut -f1 -d "_"|while read line; do echo '#!/bin/bash
#SBATCH --job-name=gatk
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
module load gatk4/4.2.5.0--hdfd78af_0
module load module load samtools/1.15--h3843a85_0                                                                                                                                           srun --export=all -n 1 -c 128 gatk --java-options "-XX:ParallelGCThreads=128" MarkDuplicates -I 01.bwa/'$line'_MorexV3.bwa.sort.bam -O 2.GATK/'$line'_MorexV3.bwa.sort.dup.bam -M 2.GATK/'$line'_MorexV3.bwa.sort.dup.metirc --REMOVE_DUPLICATES true                                                                                                                                   srun --export=all -n 1 -c 128 samtools index -@ 128 -c 2.GATK/'$line'_MorexV3.bwa.sort.dup.bam
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

###bcftools
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

##cuteSV
ls 01.minimap/*.bam|cut -f2 -d "/"|while read line; do filename=$(basename "$line"); prefix=${filename%.Morex.sort.bam};mkdir -p work_dictionary/$prefix; echo '#!/bin/bash
#SBATCH --job-name cutesv
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate cutesv
srun --export=all -n 1 -c 128 cuteSV 01.minimap/'"$line"' MorexV3.fa 04.cuteSV/'"$prefix"'.cuteSV.vcf work_dictionary/'"$prefix"' -s 5 -S '"$prefix"' -l 50 --genotype --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 -t 128' >$prefix.cuteSV.sh; done



##sniffles
ls 01.minimap/*.bam|cut -f2 -d "/"|while read line; do filename=$(basename "$line"); prefix=${filename%.Morex.sort.bam}; echo '#!/bin/bash
#SBATCH --job-name=sniffles
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=2:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate sniffle
srun --export=all -n 1 -c 128 sniffles --input 01.minimap/'"$line"' --vcf 02.sniffles/'"$prefix"'.Morex.vcf --reference MorexV3.fa --snf 02.sniffles/'"$prefix"'.Morex.snf --threads 128 --minsupport 5 --long-ins-length 100000000 --long-del-length 100000000' >$prefix.sniffle.sh ; done


##svim
ls 01.minimap/*.bam|cut -f2 -d"/"|while read line; do filename=$(basename "$line"); prefix=${filename%.Morex.sort.bam}; echo '#!/bin/bash
#SBATCH --job-name=SVIM
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate svim
srun --export=all -n 1 -c 128 svim alignment 03.SVIM/'"$prefix"' 01.minimap/'"$line"' MorexV3.fa  --minimum_depth 5 --min_sv_size 50 --sample '"$prefix"' '>$prefix.svim.sh; done


###sr minimap
ls 1.split/*|cut -f2 -d "/"|cut -f1 -d "_"|sort|uniq|while read line; do echo '#!/bin/bash
#SBATCH --job-name=Mini
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
module load singularity/4.1.0-nompi
module load samtools/1.15--h3843a85_0
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/Singularity_image/minimap2.sif minimap2 -ax sr -c --MD -R '"'@RG\tID:"$line"\tSM:"$line"'"' -I100g -t 128 -Y MorexV3.fa 1.split/'"$line"'_clean_1.fq.gz 1.split/'"$line"'_clean_2.fq.gz |samtools sort -@ 128 -O BAM -o 01.minimap/'"$line"'.Morex.sr.sort.bam' >$line.srmini.sh; done


##lr minimap

ls 1.split/*.fastq.gz|cut -f2 -d"/"|cut -f1 -d "."|while read line; do echo '#!/bin/bash
#SBATCH --job-name=Mini
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
module load singularity/4.1.0-nompi
module load samtools/1.15--h3843a85_0
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/Singularity_image/minimap2.sif minimap2 -ax map-hifi -c --MD -R '"'@RG\tID:"$line"\tSM:"$line"'"' -I100g -t 128 -Y MorexV3.fa 1.split/'"$line"'.fastq.gz |samtools sort -@ 128 -O BAM -o 01.minimap/'"$line"'.Morex.sort.bam' >$line.mini.sh; done



####sr: minimap + sniffle + deepvariant
ls 01.minimap/*.bam|cut -f2 -d "/"|while read line; do filename=$(basename "$line"); prefix=${filename%.Morex.sr.sort.bam}; echo '#!/bin/bash
#SBATCH --job-name=sniffles
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=2:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate sniffle
srun --export=all -n 1 -c 128 sniffles --input 01.minimap/'"$line"' --vcf 02.sniffles/'"$prefix"'.sr.Morex.vcf --reference MorexV3.fa --snf 02.sniffles/'"$prefix"'.sr.Morex.snf --threads 128 --minsupport 1 --long-ins-length 100000000 --long-del-length 100000000' >$prefix.sr.sniffle.sh ; done









