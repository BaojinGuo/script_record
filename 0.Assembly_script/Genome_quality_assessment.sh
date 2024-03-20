###Here, using assembly-stat, quast and BUSCO to assess the genome
##Install gfatools assembly-stats quast and BUSCO
conda install -c bioconda gfatools
conda install -c bioconda assembly-stats
singularity pull docker://staphb/quast:latest
singularity pull BUSCO.sif  docker://ezlabgva/busco:v5.6.1_cv1

###transfor *.gfa to *.fa and statistic the genome data
gfatools gfa2fa S1_hifi.asm.bp.p_ctg.gfa >S1_hifi.asm.bp.p_ctg.fa
gfatools gfa2fa S2_hifi.asm.bp.p_ctg.gfa >S2_hifi.asm.bp.p_ctg.fa
samtools faidx S1_hifi.asm.bp.p_ctg.fa
samtools faidx S2_hifi.asm.bp.p_ctg.fa
assembly-stats S1_hifi.asm.bp.p_ctg.fa >S1_hifi.asm.bp.p_ctg.fa.stats
assembly-stats S2_hifi.asm.bp.p_ctg.fa >S2_hifi.asm.bp.p_ctg.fa.stats

##stats for S1_hifi.asm.bp.p_ctg.fa
#sum = 630465039, n = 2866, ave = 219980.82, largest = 93660658
#N50 = 25161254, n = 7
#N60 = 19120655, n = 10
#N70 = 14401794, n = 13
#N80 = 5914009, n = 20
#N90 = 40153, n = 712
#N100 = 12410, n = 2866
#N_count = 0
#Gaps = 0

##stats for S2_hifi.asm.bp.p_ctg.fa
#sum = 644871918, n = 2960, ave = 217862.13, largest = 99614926
#N50 = 39209605, n = 6
#N60 = 28150159, n = 8
#N70 = 19690660, n = 11
#N80 = 1627050, n = 18
#N90 = 40851, n = 969
#N100 = 15279, n = 2960
#N_count = 0
#Gaps = 0

###quast5.2.0 the output will have statistics and some graph; the software is based on minimap2, so if there is a reference, the software could compare your assembles with it and get the uesful information like misassemblies, misassembled contigs, genome fraction, etc..  
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/quast_latest.sif python ../../software/quast-5.1.0rc1/quast.py --fragmented -f -b  -t 128 ../01.hifi_assembly/S1_HIFI_RESULT/S1_hifi.asm.bp.p_ctg.fa -r ../01.hifi_assembly/S2_HIFI_RESULT/S2_hifi.asm.bp.p_ctg.fa -o quast_to_SS --contig-thresholds 0,1000000,5000000,10000000,50000000,80000000
##Options:
#-o  --output-dir  <dirname>       Directory to store all result files [default: quast_results/results_<datetime>]
#-r                <filename>      Reference genome file
#-g  --features [type:]<filename>  File with genomic feature coordinates in the reference (GFF, BED, NCBI or TXT)
#                                  Optional 'type' can be specified for extracting only a specific feature type from GFF
#-m  --min-contig  <int>           Lower threshold for contig length [default: 500]
#-t  --threads     <int>           Maximum number of threads [default: 25% of CPUs]
#-f  --gene-finding                    Predict genes using GeneMarkS (prokaryotes, default) or GeneMark-ES (eukaryotes, use --eukaryote)
#    --mgm                             Use MetaGeneMark for gene prediction (instead of the default finder, see above)
#    --glimmer                         Use GlimmerHMM for gene prediction (instead of the default finder, see above)
#    --gene-thresholds <int,int,...>   Comma-separated list of threshold lengths of genes to search with Gene Finding module
#                                      [default: 0,300,1500,3000]
#-b  --conserved-genes-finding         Count conserved orthologs using BUSCO (only on Linux)
#    --operons  <filename>             File with operon coordinates in the reference (GFF, BED, NCBI or TXT)
#    --est-ref-size <int>              Estimated reference size (for computing NGx metrics without a reference)
#    --contig-thresholds <int,int,...> Comma-separated list of contig length thresholds [default: 0,1000,5000,10000,25000,50000]
#    --x-for-Nx <int>                  Value of 'x' for Nx, Lx, etc metrics reported in addition to N50, L50, etc (0, 100) [default: 90]
#--fragmented                      Reference genome may be fragmented into small pieces (e.g. scaffolded reference) 

###BUSCO
wget -c https://busco.ezlab.org/datasets/eukaryota_odb10.2024-01-08.tar.gz
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/BUSCO.sif busco -i /scratch/pawsey0399/bguo1/0.assembly/01.hifi_assembly/S1_HIFI_RESULT/S1_hifi.asm.bp.p_ctg.fa -o S1/S1_busco -m genome -l /scratch/pawsey0399/bguo1/0.assembly/02.quality_assessmnet/03.BUSCO/eukaryota_odb10 --cpu 128 --offline
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/BUSCO.sif busco -i /scratch/pawsey0399/bguo1/0.assembly/01.hifi_assembly/S2_HIFI_RESULT/S2_hifi.asm.bp.p_ctg.fa -o S2/S2_busco -m genome -l /scratch/pawsey0399/bguo1/0.assembly/02.quality_assessmnet/03.BUSCO/eukaryota_odb10 --cpu 128 --offline

mv $short_summary.specific.fabales_odb10.fabales_S1.txt plot/short_summary.specific.fabales_odb10.Fabales_S1.txt 
mv $short_summary.specific.fabales_odb10.fabales_S2.txt plot/short_summary.specific.fabales_odb10.Fabales_S2.txt 
mv $short_summary.specific.embryophyta_odb10.embryophyta_S1.txt plot/short_summary.specific.embryophyta_odb10.Embryophyta_S1.txt 
mv $short_summary.specific.embryophyta_odb10.embryophyta_S2.txt plot/short_summary.specific.embryophyta_odb10.Embryophyta_S2.txt

generate_plot.py -wd plot/

