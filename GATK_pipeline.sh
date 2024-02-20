####30X resequencing data of RILs parental varieties BASS and Flinders
###Reference Morexv1 and Morexv3, BWA and picard index
##first, Reference Morexv1.0
bwa index 160404_barley_pseudomolecules_masked.fasta
samtools faidx 160404_barley_pseudomolecules_masked.fasta
gatk CreateSequenceDictionary -R 160404_barley_pseudomolecules_masked.fasta -O 160404_barley_pseudomolecules_masked.dict
bwa mem -t 15 -R "@RG\tID:BASS\tPL:illumina\tLB:library\tSM:BASS" /data/baojin/00.RILs_parent_reseq/Reference/160404_barley_pseudomolecules_masked.fasta /data/baojin/00.RILs_parent_reseq/raw_data/Bass/BASS_clean_1.fq.gz /data/baojin/00.RILs_parent_reseq/raw_data/Bass/BASS_clean_2.fq.gz | samtools sort -@ 28 -o BASS_morexv1.sort.bam
bwa mem -t 15 -R "@RG\tID:Flinders\tPL:illumina\tLB:library\tSM:Flinders" /data/baojin/00.RILs_parent_reseq/Reference/160404_barley_pseudomolecules_masked.fasta /data/baojin/00.RILs_parent_reseq/raw_data/Flinders/Flinders_clean_1.fq.gz /data/baojin/00.RILs_parent_reseq/raw_data/Flinders/Flinders_clean_2.fq.gz | samtools sort -@ 15 -o Flinders_morexv1.sort.bam
gatk --java-options "-Xmx50g -XX:ParallelGCThreads=15" MarkDuplicates -I ../02.bwa_align/BASS_morexv1.sort.bam -O BASS_morexv1.markdup.bam -M BASS_morexv1.metric
gatk --java-options "-Xmx50g -XX:ParallelGCThreads=15" MarkDuplicates -I ../02.bwa_align/Flinders_morexv1.sort.bam -O Flinders_morexv1.markdup.bam -M Flinders_morexv1.metric
samtools view -h -b -@ 15 -q 30 BASS_morexv1.markdup.bam |samtools sort -o BASS_morexv1.markdup.q30.bam -@ 15 -
samtools view -h -b -@ 15 -q 30 Flinders_morexv1.markdup.bam |samtools sort -o Flinders_morexv1.markdup.q30.bam -@ 15 -
samtools index -c BASS_morexv1.markdup.q30.bam
samtools index -c Flinders_morexv1.markdup.q30.bam
gatk HaplotypeCaller -OVI False -ERC GVCF -R /data/baojin/00.RILs_parent_reseq/Reference/160404_barley_pseudomolecules_masked.fasta -I BASS_morexv1.markdup.q30.bam -O BASS_morexv1.q30.gvcf
gatk HaplotypeCaller -OVI False -ERC GVCF -R /data/baojin/00.RILs_parent_reseq/Reference/160404_barley_pseudomolecules_masked.fasta -I Flinders_morexv1.markdup.q30.bam -O Flinders_morexv1.q30.gvcf
gatk CombineGVCFs -R ../Reference/160404_barley_pseudomolecules_masked.fasta --variant BASS_morexv1.q30.gvcf --variant Flinders_morexv1.q30.gvcf -O morexv1.q30.com.vcf
gatk GenotypeGVCFs -OVI False -R ../Reference/160404_barley_pseudomolecules_masked.fasta -V morexv1.q30.com.vcf -O morexv1.q30.geno.vcf
gatk SelectVariants -R ../Reference/160404_barley_pseudomolecules_masked.fasta -V morexv1.q30.geno.vcf --select-type-to-include SNP -O morexv1.q30.geno.snp.vcf
gatk VariantFiltration -O morexv1.q30.geno.snp.filt.vcf -V morexv1.q30.geno.snp.vcf --cluster-window-size 10 --filter-expression "DP < 5 " --filter-name "LowCoverage" --filter-expression "QUAL < 30.0 " --filter-name "VeryLowQual" --filter-expression "QUAL > 30.0 && QUAL < 50.0 " --filter-name "LowQual" --filter-expression "QD < 1.5 " --filter-name "LowQD"
gatk SelectVariants -R ../Reference/160404_barley_pseudomolecules_masked.fasta -V morexv1.q30.geno.vcf --select-type-to-include INDEL -O morexv1.q30.geno.indel.vcf
gatk VariantFiltration -O morexv1.q30.geno.indel.filt.vcf -V morexv1.q30.geno.indel.vcf --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0" --filter-name "Filter" 
bcftools view -i 'FILTER=="PASS"' morexv1.q30.geno.indel.filt.vcf >morexv1.q30.geno.indel.PASS.vcf
bcftools view -i 'FILTER=="PASS"' morexv1.q30.geno.snp.filt.vcf >morexv1.q30.geno.snp.PASS.vcf
bgzip morexv1.q30.geno.indel.PASS.vcf
bgzip morexv1.q30.geno.snp.PASS.vcf
tabix -C morexv1.q30.geno.snp.PASS.vcf.gz
tabix -C morexv1.q30.geno.indel.PASS.vcf.gz
bcftools concat -a -O z morexv1.q30.geno.snp.PASS.vcf.gz morexv1.q30.geno.indel.PASS.vcf.gz -o compare-variant/morexv1.q30.variant.PASS.vcf.gz
##bgzip morexv1.q30.variant.PASS.vcf
##tabix -C morexv1.q30.variant.PASS.vcf.gz

##Tassel input morexv1.q30.variant.PASS.vcf.gz, output morexv1.q30.variant.PASS-noinfo.vcf
run_pipeline.pl -Xmx50g -fork1 -vcf BF.morexv3.q30.variants.PASS.vcf.gz -sortPositions -export BF.morexv3.q30.variants.PASS.noinf.vcf -exportType VCF -exportIncludeAnno false -exportIncludeDepth false -runfork1

###R program, print position.txt
data<-read.table("morexv1.q30.variant.PASS-noinfo.vcf",header = F)
head(data)
diff<-subset(data,data$V10 != data$V11)
head(diff)
write.table(diff[,1:2],"diff.txt",sep = "\t",quote=FALSE,row.names = FALSE)
##subset the different variants between parents, remember to modify the chromosome row due to the abbrevation of Tassel output
vcftools --gzvcf morexv1.q30.variant.PASS.vcf.gz --positions diff.txt --recode --recode-INFO-all --out morexv1.q30.variant.PASS-diff.vcf


###snpEFF gene annotation
java -jar snpEff.jar build -gff3 -v Morexv1
java -jar snpEff.jar Morexv1 data/morexv1.q30.variant.PASS-diff.vcf.recode.vcf >morexv1.q30.diff.eff.vcf
###extract snps in the region of QTL
vcftools --vcf morexv1.q30.diff.eff.vcf --chr chr6H --from-bp 384800000 --to-bp 413640000 --recode --recode-INFO-all --out TGW-QTL.morexv1.diff.eff.vcf
###extract the related severe mutants
java -jar /data/baojin/software/snpEff/SnpSift.jar filter "ANN[0].IMPACT='HIGH' || ANN[0].IMPACT='MODERATE'" TGW-QTL.morexv1.diff.eff.vcf.recode.vcf >TGW-QTL.morexv1.diff.eff-severe.vcf


