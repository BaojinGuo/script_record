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

####bcftools
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


######filter raw variants
###deepvariants
ls *.bwa.deep.vcf.gz|tail -n4|cut -f1 -d"."|while read line; do bcftools view -f PASS $line.bwa.deep.vcf.gz -Oz -o $line.bwa.deep.pass.vcf.gz; done 
ls *pass.vcf.gz|while read line; do tabix -C $line; done
ls *.bwa.deep.pass.vcf.gz|cut -f1 -d"."|while read line; do bcftools view -i 'GT="1/1"' -s $line -Ov -o $line.bwa.deep.pass.hom.vcf $line.bwa.deep.pass.vcf.gz; done
bcftools view -i 'DP>10' RGT.bwa.deep.pass.hom.vcf -Ov -o RGT.bwa.deep.PASS.hom.vcf

###GATK
ls *.bwa.vcf |cut -f1 -d "."|cut -f1 -d "_"|while read line; do srun -c 128 -n 1 -p debug -A pawsey0399 gatk SelectVariants -R ../../MorexV3.fa -V ${line}_MorexV3.bwa.vcf --select-type-to-include SNP -O ${line}_MorexV3.bwa.snp.vcf; done
ls *.bwa.vcf |cut -f1 -d "."|cut -f1 -d "_"|while read line; do srun -c 128 -n 1 -p debug -A pawsey0399 gatk SelectVariants -R ../../MorexV3.fa -V ${line}_MorexV3.bwa.vcf --select-type-to-include INDEL -O ${line}_MorexV3.bwa.indel.vcf; done

ls *.indel.vcf|cut -f1 -d"_"|while read line; do gatk VariantFiltration -V ${line}_MorexV3.bwa.indel.vcf --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0||DP<10||QUAL<50" --filter-name "Filter" -O ${line}_MorexV3.bwa.indel.filter.vcf; done
ls *.snp.vcf|cut -f1 -d"_"|while read line; do gatk VariantFiltration -V ${line}_MorexV3.bwa.snp.vcf --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0||DP<10||QUAL<50" --filter-name "Filter" -O ${line}_MorexV3.bwa.snp.filter.vcf; done

ls *.indel.vcf|cut -f1 -d"_"|while read line; do bcftools view -f PASS ${line}_MorexV3.bwa.indel.filter.vcf -Oz -o ${line}_MorexV3.bwa.indel.PASS.vcf.gz; done
ls *.snp.vcf|cut -f1 -d"_"|while read line; do bcftools view -f PASS ${line}_MorexV3.bwa.snp.filter.vcf -Oz -o ${line}_MorexV3.bwa.snp.PASS.vcf.gz; done

ls *PASS.vcf.gz|while read line; do tabix -C $line; done

ls *snp.PASS.vcf.gz|cut -f1 -d"_"|while read line; do gatk MergeVcfs -I ${line}_MorexV3.bwa.snp.PASS.vcf.gz -I ${line}_MorexV3.bwa.indel.PASS.vcf.gz -O ${line}_MorexV3.bwa.all.PASS.vcf; done

ls *_MorexV3.bwa.all.PASS.vcf|cut -f1 -d"_"|while read line; do bcftools view -i 'GT="1/1"' -s $line -Ov -o ${line}_MorexV3.bwa.all.PASS.hom.vcf ${line}_MorexV3.bwa.all.PASS.vcf; done



#####bcftools
snp filter rules:
"QUAL < 50 || DP>2*$DP|| MQBZ < -(3.5+4*DP/QUAL) || RPBZ > (3+3*DP/QUAL) || RPBZ < -(3+3*DP/QUAL) || FORMAT/SP > (40+DP/2) || SCBZ > (2.5+DP/30)"

"QUAL < 50 || DP < 10 || DP > 90 || MQBZ < -4.3 || RPBZ > 3.6 || RPBZ < -3.6 || SCBZ > 2.88 || MQ < 40"
indel filter rules
"QUAL < 50 || DP < 10 || DP > 90 ||IDV < 2 || IMF < 0.02+(($qual+1)/($qual+31))*(($qual+1)/($qual+31))/4 || DP > ($DP/2) * (1.7 + 12/($qual+20)) || MQBZ < -(5+DP/20) || RPBZ+SCBZ > 9"

"QUAL < 50 || DP < 10 || DP > 40 ||IDV < 2 || IMF < 0.12 || MQ < 40 || MQBZ < -5.5 || RPBZ+SCBZ > 9"


ls *.variants.vcf.gz|cut -f1 -d"."|while read line; do bcftools view -i 'TYPE="INDEL"' $line.variants.vcf.gz -Oz -o $line.indel.vcf.gz; done
ls *indel.vcf.gz|while read line; do tabix -C $line; done

ls *.variants.vcf.gz|cut -f1 -d"."|tail -n4|while read line; do bcftools view -e 'TYPE="INDEL"' $line.variants.vcf.gz -Oz -o $line.snp.vcf.gz; done
ls *snp.vcf.gz|while read line; do tabix -C $line; done

ls *indel.vcf.gz|cut -f1 -d "."|while read line; do bcftools view -e "QUAL < 50 || DP < 10 || DP > 40 ||IDV < 2 || IMF < 0.12 || MQ < 40 || MQBZ < -5.5 || RPBZ+SCBZ > 9" -Oz -o $line.indel.PASS.vcf.gz $line.indel.vcf.gz
ls *snp.vcf.gz|cut -f1 -d "."|while read line; do bcftools view -e "QUAL < 50 || DP < 10 || DP > 90 || MQBZ < -4.3 || RPBZ > 3.6 || RPBZ < -3.6 || SCBZ > 2.88 || MQ < 40" -Oz -o $line.snp.PASS.vcf.gz $line.snp.vcf.gz

ls *.PASS.vcf.gz|while read line; do tabix -C $line; done

gatk MergeVcfs -I -I -O

ls *PASS.vcf|cut -f1 -d "."|while read line;do bcftools view -i 'GT="1/1"' -s $line -Ov -o $line.PASS.hom.vcf $line.PASS.vcf; done

##########################
####step2 find the variants exist in 2 of 3 callers
ls *.PASS.hom.vcf|cut -f1 -d"."|sort|uniq|while read line; do bcftools sort -Oz -o $line.bcf.PASS.hom.sort.vcf.gz $line.bcf.PASS.hom.vcf; done 
ls *.PASS.hom.vcf|cut -f1 -d"."|sort|uniq|while read line; do bcftools sort -Oz -o $line.GATK.PASS.hom.sort.vcf.gz $line.GATK.PASS.hom.vcf; done 
ls *.PASS.hom.vcf|cut -f1 -d"."|sort|uniq|while read line; do bcftools sort -Oz -o $line.deep.PASS.hom.sort.vcf.gz $line.deep.PASS.hom.vcf; done 
ls *vcf.gz|while read line; do tabix -C $line;done

ls *.PASS.hom.vcf|cut -f1 -d"."|sort|uniq|while read line; do bcftools merge -m none -Ov -o $line.merge.vcf $line.bcf.PASS.hom.sort.vcf.gz $line.deep.PASS.hom.sort.vcf.gz $line.GATK.PASS.hom.sort.vcf.gz --force-samples

ls *merge.vcf|cut -f1 -d"."|while read line; do bcftools view -i 'COUNT(GT!="./.")>=2' -Ov -o $line.benchmark.vcf $line.merge.vcf ;  done


ls *merge.vcf|cut -f1 -d"."|while read line; do 
awk 'BEGIN {OFS="\t"} !/^#/ {
    split($10, sample1, ":");
    split($11, sample2, ":");
    split($12, sample3, ":");
    print $1, $2, $4, $5, sample1[1], sample2[1], sample3[1];
    }' $line.merge.vcf >$line.merge.simplify.txt

python 03.statistics_illumina_benchmark.py Hockett.merge.simplify.txt Hockett.merge












