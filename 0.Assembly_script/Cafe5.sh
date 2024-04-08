###Cafe5 is a software for gene expansion and constriction, https://github.com/hahnlab/CAFE5
##Installation
git clone https://github.com/hahnlab/CAFE5/releases
cd CAFE5
./configure
make
##CAFE5 needs two input files, gene family counts and tree file with divergence time which came from Orthofinder and r8s or mcmctree, respectively.
#gene family counts could ues Orthogroups.GeneCount.tsv form Orthofinder, then change the format to CAFE5 using
awk 'BEGIN{FS="\t"; OFS="\t"} {NF--; print "(null)", $0}'  Orthogroups.GeneCount.tsv>Gene.counts
#next, delete the orthogroups which have large copy number variations (>100)
awk 'NR==1 || $3<100 && $4<100 && $5<100 && $6<100 && $7<100 && $8<100 && $9<100 && $10<100 && $11<100 && $12<100 && $13<100 && $14<100 && $15<100 && $16<100 && $17<100 && $18<100 && $19<100 && $20<100 {print $0}' Gene.counts >Gene.counts.filter

#tree files preparation
#install r8s
singularity pull r8s.sif docker://reslp/r8s:1.81
#get script
git clone https://github.com/shanedenecke/SLC_ID_SCRIPTS.git
#the tree file is from orthofinder Species tree, then infer ultrametric tree, first step is to write r8s running script
python cafetutorial_prep_r8s.py -i SpeciesTree_rooted.txt -o r8s_ctl_file.txt -s 272 -p 'Arabidopsis_thaliana,Lupinus_angustifolius' -c '108' #-i tree file; -o output file; -s number of sites in alignment that was used to infer species tree; -p list of species name tuples e.g. [('ENSG00','ENSPTR'),('ENSFCA','ENSECA')]; -c divergence tims, list of flats e.g. [6.4,80].
#calculate divergence time of all nodes by r8s
r8s -b -f r8s_ctl_file.txt >r8s_tmp.txt
tail -n 1 r8s_tmp.txt|cut -c 16- >r8s_ultrametric.txt
#running CAFE5
cafe5 -c 128 -i Gene.counts.filter -t r8s_ultrametric.txt -o out #-c cores, CPU number; -i gene family numbers; -t tree files; -o output path; -p possion distribution otherwise uniform distribution

#subset any species significant expansion and constriction genes
cat Base_family_results.txt |grep 'y' |cut -f1|grep -f - Base_change.tab |cut -f1,3|awk '$2>0' >S2.sinnificant.expand.ID
cat Base_family_results.txt |grep 'y' |cut -f1|grep -f - Base_change.tab |cut -f1,3|awk '$2<0' >S2.sinnificant.constrict.ID
cut -f1 S1.sinnificant.expand.ID |grep -f - ../../06.Orthofinder/0.workdic/02.prj/OrthoFinder/Results_Apr04/Orthogroups/Orthogroups.tsv |cut -f14|sed "s/ /\n/g"|sed "s/\t/\n/g"| sed "s/,//g"|sort|uniq >S1.significant.expand.genes
cut -f1 S1.sinnificant.constrict.ID |grep -f - ../../06.Orthofinder/0.workdic/02.prj/OrthoFinder/Results_Apr04/Orthogroups/Orthogroups.tsv |cut -f14|sed "s/ /\n/g"|sed "s/\t/\n/g"| sed "s/,//g"|sort|uniq >S1.significant.constrict.genes

cut -f1 S2.sinnificant.expand.ID |grep -f - ../../06.Orthofinder/0.workdic/02.prj/OrthoFinder/Results_Apr04/Orthogroups/Orthogroups.tsv |cut -f15|sed "s/ /\n/g"|sed "s/\t/\n/g"| sed "s/,//g"|sort|uniq >S2.significant.expand.genes
cut -f1 S2.sinnificant.constrict.ID |grep -f - ../../06.Orthofinder/0.workdic/02.prj/OrthoFinder/Results_Apr04/Orthogroups/Orthogroups.tsv |cut -f15|sed "s/ /\n/g"|sed "s/\t/\n/g"| sed "s/,//g"|sort|uniq >S2.significant.constrict.genes

#final step is make enrichment analysis

