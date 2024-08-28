#!/bin/bash
#SBATCH --job-name=BWA
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
module load singularity/4.1.0-nompi
module load bwa/0.7.17--h7132678_9
module load samtools/1.15--h3843a85_0
srun --export=all -n 1 -c 128 bwa mem -t 128 -R "@RG\tID:H002\tPL:illumina\tLB:library\tSM:H002" Clipper.V1.fasta 00.Raw_data/H002_clean_1.fq.gz 00.Raw_data/H002_clean_2.fq.gz | samtools sort -@ 128 -o 02.bwa/H002_Clipper.bwa.sort.bam -m 100M -T /scratch/pawsey0399/bguo1/TMP
srun -c 128 -n 1 --export=all samtools index -c -@ 128 02.bwa/H002_Clipper.bwa.sort.bam
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/Singularity_image/deepvariant_latest.sif run_deepvariant --model_type WGS --ref Clipper.V1.fasta --regions chr1H --reads 02.bwa/H002_Clipper.bwa.sort.bam --sample_name H002 --output_gvcf 03.Deepvariant/H002_Clipper.1H.bwa.deep.gvcf.gz --output_vcf 03.Deepvariant/H002_Clipper.1H.bwa.deep.vcf.gz --intermediate_results_dir 03.Deepvariant/H002.1H.deep.intermediate --num_shards 128 --logging_dir 03.Deepvariant/H002.1H.deep.log
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/Singularity_image/deepvariant_latest.sif run_deepvariant --model_type WGS --ref Clipper.V1.fasta --regions chr2H --reads 02.bwa/H002_Clipper.bwa.sort.bam --sample_name H002 --output_gvcf 03.Deepvariant/H002_Clipper.2H.bwa.deep.gvcf.gz --output_vcf 03.Deepvariant/H002_Clipper.2H.bwa.deep.vcf.gz --intermediate_results_dir 03.Deepvariant/H002.2H.deep.intermediate --num_shards 128 --logging_dir 03.Deepvariant/H002.2H.deep.log
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/Singularity_image/deepvariant_latest.sif run_deepvariant --model_type WGS --ref Clipper.V1.fasta --regions chr3H --reads 02.bwa/H002_Clipper.bwa.sort.bam --sample_name H002 --output_gvcf 03.Deepvariant/H002_Clipper.3H.bwa.deep.gvcf.gz --output_vcf 03.Deepvariant/H002_Clipper.3H.bwa.deep.vcf.gz --intermediate_results_dir 03.Deepvariant/H002.3H.deep.intermediate --num_shards 128 --logging_dir 03.Deepvariant/H002.3H.deep.log
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/Singularity_image/deepvariant_latest.sif run_deepvariant --model_type WGS --ref Clipper.V1.fasta --regions chr4H --reads 02.bwa/H002_Clipper.bwa.sort.bam --sample_name H002 --output_gvcf 03.Deepvariant/H002_Clipper.4H.bwa.deep.gvcf.gz --output_vcf 03.Deepvariant/H002_Clipper.4H.bwa.deep.vcf.gz --intermediate_results_dir 03.Deepvariant/H002.4H.deep.intermediate --num_shards 128 --logging_dir 03.Deepvariant/H002.4H.deep.log
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/Singularity_image/deepvariant_latest.sif run_deepvariant --model_type WGS --ref Clipper.V1.fasta --regions chr5H --reads 02.bwa/H002_Clipper.bwa.sort.bam --sample_name H002 --output_gvcf 03.Deepvariant/H002_Clipper.5H.bwa.deep.gvcf.gz --output_vcf 03.Deepvariant/H002_Clipper.5H.bwa.deep.vcf.gz --intermediate_results_dir 03.Deepvariant/H002.5H.deep.intermediate --num_shards 128 --logging_dir 03.Deepvariant/H002.5H.deep.log
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/Singularity_image/deepvariant_latest.sif run_deepvariant --model_type WGS --ref Clipper.V1.fasta --regions chr6H --reads 02.bwa/H002_Clipper.bwa.sort.bam --sample_name H002 --output_gvcf 03.Deepvariant/H002_Clipper.6H.bwa.deep.gvcf.gz --output_vcf 03.Deepvariant/H002_Clipper.6H.bwa.deep.vcf.gz --intermediate_results_dir 03.Deepvariant/H002.6H.deep.intermediate --num_shards 128 --logging_dir 03.Deepvariant/H002.6H.deep.log
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/Singularity_image/deepvariant_latest.sif run_deepvariant --model_type WGS --ref Clipper.V1.fasta --regions chr7H --reads 02.bwa/H002_Clipper.bwa.sort.bam --sample_name H002 --output_gvcf 03.Deepvariant/H002_Clipper.7H.bwa.deep.gvcf.gz --output_vcf 03.Deepvariant/H002_Clipper.7H.bwa.deep.vcf.gz --intermediate_results_dir 03.Deepvariant/H002.7H.deep.intermediate --num_shards 128 --logging_dir 03.Deepvariant/H002.7H.deep.log

########it is an example for 1H, actrually we run 1H to 7H seperately.
#!/bin/bash --login
#SBATCH --job-name=gln
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=96:00:00
#SBATCH --account=pawsey0399
#SBATCH --export=NONE
module load bcftools/1.15--haf5b3da_0
module load singularity/4.1.0-nompi
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/Singularity_image/Glnexus.sif glnexus_cli --config DeepVariant --bed 1H --threads 128 /scratch/pawsey0399/bguo1/Murdoch/08.He/03.Deepvariant/*1H*.gvcf.gz |bcftools view --threads 128 - |bgzip -@ 128 -c > WGRS_Clipper.1H.cohort.vcf.gz


######remove non-AUS samples
vcftools --gzvcf 1H/WGRS_Clipper.1H.cohort.vcf.gz --remove remove_no_aus.txt --recode --recode-INFO-all --out WGRS_Clipper.1H.AUS.vcf
######make genotype with 'DP <10' as './.' for each sample
bcftools +setGT WGRS_Clipper.1H.AUS.vcf.recode.vcf -- -t q -i 'DP<10' -n './.' > WGRS_Clipper.1H.AUS.DP10.vcf
#######extract snps only
bcftools view -v snps WGRS_Clipper.1H.AUS.DP10.vcf -o WGRS_Clipper.1H.AUS.DP10.snps.vcf
#######convert format vcf to bcf and index it for software local_pca
bcftools convert -O b WGRS_Clipper.1H.AUS.DP10.snps.vcf >WGRS_Clipper.1H.AUS.DP10.snps.bcf
bcftools index WGRS_Clipper.1H.AUS.DP10.snps.bcf
#######convert format vcf to bcf and index it for software fastEPRR
bcftools convert -O z WGRS_Clipper.1H.AUS.DP10.snps.vcf >WGRS_Clipper.1H.AUS.DP10.snps.vcf.gz



#############local_pca
####install packages
install.packages("devtools")
install.packages("data.table")
devtools::install_github("petrelharp/local_pca/lostruct")
library(lostruct)

snps <- vcf_windower("WGRS_Clipper.1H.AUS.DP10.snps.bcf",size=5000000,type='bp') ######type = 'bp' means split windows with position, and type = "snps" means split it with the number of snps.
pcs <- eigen_windows(snps,k=2)
pcdist <- pc_dist(pcs,npc=2)
write.table(pcs,"WGRS_Clipper.1H.AUS.DP10.snps.PCA_W5M.txt",sep="\t")
write.table(pcdist,"WGRS_Clipper.1H.AUS.DP10.snps.PCA_W5M.pwdist.txt",sep="\t")

fit2d <- cmdscale(pcdist, eig=TRUE, k=2 )

plot(fit2d$points, xlab="Coordinate 1", ylab="Coordinate 2", col=rainbow(1.2*nrow(pcdist)) )

 Rscript run_on_barley.R -t bp -s 1000000 -o barley_pca_results

#############FastEPRR






