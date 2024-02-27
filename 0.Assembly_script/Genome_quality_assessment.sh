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





