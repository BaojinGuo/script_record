#####BWA0.7+GATK4.2.5
###assmbly create index
bwa index S1_hifi.asm.bp.p_ctg.fa
bwa index S2_hifi.asm.bp.p_ctg.fa
samtools faidx S1_hifi.asm.bp.p_ctg.fa
samtools faidx S2_hifi.asm.bp.p_ctg.fa
gatk CreateSequenceDictionary -R S1_hifi.asm.bp.p_ctg.fa -O S1_hifi.asm.bp.p_ctg.dict
gatk CreateSequenceDictionary -R S2_hifi.asm.bp.p_ctg.fa -O S2_hifi.asm.bp.p_ctg.dict
##Mapping & Call variants
#S1
srun --export=all -n 1 -c 128 bwa mem -t 128 -R '@RG\tID:S1\tSM:S1\tPL:PACBIO' /scratch/pawsey0399/bguo1/0.assembly/02.quality_assessmnet/02.remap/00.index/S1_hifi.asm.bp.p_ctg.fa /scratch/pawsey0399/bguo1/0.assembly/Raw_data/s1_hifi.fastq.gz |samtools sort -@ 128 -o S1_remap_sort.bam
srun --export=all -n 1 -c 128 gatk --java-options "-Xmx63g -XX:ParallelGCThreads=128" MarkDuplicates -I S1_remap_sort.bam -O S1_remap_sort_dup.bam -M S1.metric
srun --export=all -n 1 -c 128 gatk  HaplotypeCaller -OVI True -ERC GVCF -R /scratch/pawsey0399/bguo1/0.assembly/02.quality_assessmnet/02.remap/00.index/S1_hifi.asm.bp.p_ctg.fa -I S1_remap_sort_dup.bam -O S1_sort_dup.g.vcf
srun --export=all -n 1 -c 128 gatk GenotypeGVCFs -OVI True -R /scratch/pawsey0399/bguo1/0.assembly/02.quality_assessmnet/02.remap/00.index/S1_hifi.asm.bp.p_ctg.fa -V S1_sort_dup.g.vcf -O S1_sort_dup.vcf
#S2
srun --export=all -n 1 -c 128 bwa mem -t 128 -R '@RG\tID:S2\tSM:S2\tPL:PACBIO' /scratch/pawsey0399/bguo1/0.assembly/02.quality_assessmnet/02.remap/00.index/S2_hifi.asm.bp.p_ctg.fa /scratch/pawsey0399/bguo1/0.assembly/Raw_data/s2_hifi.fastq.gz |samtools sort -@ 128 -o S2_remap_sort.bam
srun --export=all -n 1 -c 128 gatk --java-options "-Xmx63g -XX:ParallelGCThreads=128" MarkDuplicates -I S2_remap_sort.bam -O S2_remap_sort_dup.bam -M S2.metric
srun --export=all -n 1 -c 128 gatk  HaplotypeCaller -OVI True -ERC GVCF -R /scratch/pawsey0399/bguo1/0.assembly/02.quality_assessmnet/02.remap/00.index/S2_hifi.asm.bp.p_ctg.fa -I S2_remap_sort_dup.bam -O S2_sort_dup.g.vcf
srun --export=all -n 1 -c 128 gatk GenotypeGVCFs -OVI True -R /scratch/pawsey0399/bguo1/0.assembly/02.quality_assessmnet/02.remap/00.index/S2_hifi.asm.bp.p_ctg.fa -V S2_sort_dup.g.vcf -O S2_sort_dup.vcf

#####Minimap2+DeepVariants
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/minimap2.sif minimap2 -ax map-pb -t 128 /scratch/pawsey0399/bguo1/0.assembly/02.quality_assessmnet/02.remap/00.index/S1_hifi.asm.bp.p_ctg.fa /scratch/pawsey0399/bguo1/0.assembly/Raw_data/s1_hifi.fastq.gz |samtools sort -@ 128 -o S1/S1_mini2_remap.sort.bam
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/minimap2.sif minimap2 -ax map-pb -t 128 /scratch/pawsey0399/bguo1/0.assembly/02.quality_assessmnet/02.remap/00.index/S2_hifi.asm.bp.p_ctg.fa /scratch/pawsey0399/bguo1/0.assembly/Raw_data/s2_hifi.fastq.gz |samtools sort -@ 128 -o S2/S2_mini2_remap.sort.bam

run_deepvariant.py --model_type=PACBIO --ref $Genome --reads {$SAMPLE}.mini2.remap.sort.bam --output_vcf {$SAMPLE}.deepvariants.vcf.gz --intermediate_results_dir {$SAMPLE}_intermediate_results --num_shards=128 --logging_dir {$SAMPLE}_logging
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/deepvariant_latest.sif run_deepvariant --model_type PACBIO --ref /scratch/pawsey0399/bguo1/0.assembly/02.quality_assessmnet/02.remap/00.index/S1_hifi.asm.bp.p_ctg.fa --reads S1/S1_mini2_remap.sort.bam --sample_name=S1 --output_gvcf S1/S1.mini2.deepvariant.g.vcf.gz --output_vcf S1/S1.mini2.deepvariant.vcf.gz --intermediate_results_dir S1/S1_deepvariant_intermediate --num_shards 128  --logging_dir S1/S1_deepvariant_log
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/deepvariant_latest.sif run_deepvariant --model_type PACBIO --ref /scratch/pawsey0399/bguo1/0.assembly/02.quality_assessmnet/02.remap/00.index/S2_hifi.asm.bp.p_ctg.fa --reads S2/S2_mini2_remap.sort.bam --sample_name=S2 --output_gvcf S2/S2.mini2.deepvariant.g.vcf.gz --output_vcf S2/S2.mini2.deepvariant.vcf.gz --intermediate_results_dir S2/S2_deepvariant_intermediate --num_shards 128  --logging_dir S2/S2_deepvariant_log



