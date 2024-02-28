#######This pipelines is for genome assembly (based on HIFI reads), quality assessment, TE annotation & gene annotation
#######
#######The installation steps are missing here. If needed, please refer to the specific step-by-step script in this folder
SAMPLE=SAMPLE
HIFIreads=HIFIreads.fastq.gz
HIFIoutput=HIFI_assembly.asm ###this is a prefix, whole name of used file here will be HIFI_assembly.asm.bp.p_ctg.gfa
Genome=HIFI_assembly.fa
Ref=Ref.fa
01.Genome assembly
hifiasm -o $HIFIoutput -t128 -l0 $HIFIreads
gfatools gfa2fa $HIFIoutput > $Genome

02.Genome quality assessment
###the genome quality should be assessed in contiguity, completeness and accuracy, AND general information must be included.
###contiguity is based on N50, L50 and BUSCO, some higher LTR assembly could be based on LAI (LTR index)
###completeness is based on BUSCO
###accuracy is based on reads remapping -- short reads for BWA+GATK &long reads for minimap2+DeepVariant

assembly-stats $Genome >{$Genome}.stats ###general information about whole length, contig number&length, N50 etc.
python quast.py --fragmented -f -b -t 128 $Genome -r $Ref -o $quast_results --contig-thresholds 0,1000000,5000000,10000000,50000000,80000000 ###GC contents, N50, L50 etc. and some useful graphs
busco -i $Genome -m genome -l eukaryota_odb10 --cpu 128 --offline ####eukaryota_odb10 should be downloaded at wget -c https://busco.ezlab.org/datasets/eukaryota_odb10.2024-01-08.tar.gz and uncompressed; busco assess the contiguity and completeness of the assembly.

####by BWA, it is better to use short reads rather than longreads due to BWA is mainly a short reads mapping software, the long reads mapping results is not desirable than Minimap2
bwa index $Genome
samtools faidx $Genome
gatk CreateSequenceDictionary -R $Genome -O $Genome
bwa mem -t 128 -R '@RG\tID:S1\tSM:S1\tPL:PACBIO' $Genome $HIFIreads |samtools sort -@ 128 -o {$SAMPLE}.sort.bam
gatk --java-options "-Xmx63g -XX:ParallelGCThreads=128" MarkDuplicates -I {$SAMPLE}.sort.bam -O {$SAMPLE}.sort.dup.bam -M {$SAMPLE}.metric
gatk  HaplotypeCaller -OVI True -ERC GVCF -R $Genome -I {$SAMPLE}.sort.dup.bam -O {$SAMPLE}.sort.dup.g.vcf
gatk GenotypeGVCFs -OVI True -R $Genome -V {$SAMPLE}.sort.dup.g.vcf -O {$SAMPLE}.sort.dup.vcf

####Minimap2
minimap2 -ax map-pb -t 128 $Genome $HIFIreads|samtools sort -@ 128 -o {$SAMPLE}.mini2.remap.sort.bam
samtools index -@ 128 {$SAMPLE}.mini2.remap.sort.bam
run_deepvariant --model_type=PACBIO --ref $Genome --reads {$SAMPLE}.mini2.remap.sort.bam --output {$SAMPLE}.deepvariants.vcf.gz --intermediate_results_dir {$SAMPLE}_intermediate_results --num_shards 128 --logging_dir {$SAMPLE}_logging


03.TE annotation
####EDTA, 
EDTA.pl --genome $Genome --species others --step all --cds homo.cdhit.cds --sensitive 1 --anno 1 --evaluate 1 -t 128 ##homo.cdhit.cds is the nonredundant homologous genes based on related species and removed duplicated by CD-Hit 

####RepeatMasker
##build library
BuildDatabase --name Ser1 $Genome
RepeatModeler -database Ser1 --threads 128 -LTRStruct
RepeatMasker $Genome -lib Ser1-families.fa -e rmblast -xsmall -s gff -pa 128

04.Gene annotation
####Braker-den ovo prediction, if RNA-Seq data is available, using in this step
braker.pl --genome $Genome --species=Ser1 --prot_seq allhomo.pep.cdhit.fa --useexisting --threads 128 --workingdir={$SAMPLE}_workdiction ##allhomo.pep.cdhit.fa is the nonredundant homologous proteins based on related species and removed duplicated by CD-Hit 
####Homology-based prediction
###Gemoma
java -XX:ParallelGCThreads=128 $GeMoMa-1.9.jar CLI GeMoMaPipeline threads=128 AnnotationFinalizer.r=NO p=false o=true tblastn=true t=$Genome outdir={$SAMPLE}_outdir s=own i=$Related_species_name a=$related_species_annotation.gff3 g=$related_species.fa ###related species can be multiple
###genomeThreader
gth -genomic $Genome -protein allhomo.pep.cdhit.fa -o {$SAMPLE}_gth.gff -gff3out -intermediate
###exonerate
exonerate -q allhomo.pep.cdhit.fa -t $Genome --model pretein2genome --blastn 1 --showtargetgff --showalignment no >{$SAMPLE}_exonerate.gff ####high memmory required, due to memmory limited, using --targetchunkid --targetchunktotal


####Evidence Moduler, Integration of gene prediction results













