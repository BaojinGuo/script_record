##Ref https://zhuanlan.zhihu.com/p/379464361
########server is nimbus and setonix (perl has some problems in nimbus for this script's last step EvidenceModeler)

####Gene Prediction -- do novo prediction, homology-based prediction, transcriptome-based prediction


####Software preparation
###Install the packages -- AUGUSTUSv3.2.2
##conda install -c bioconda augustus==3.2.2
###error while loading shared libraries: libboost_iostreams.so.1.60.0: solved sudo ln -s libboost_iostreams.so.1.65.1 libboost_iostreams.so.1.60.0 (https://www.biostars.org/p/389360/)
augustus

###Install geneidv1.4
##git clone https://github.com/guigolab/geneid
##make
##vim ~/.bashrc
##export PATH="/data/baojin/software/geneid/bin:$PATH"
##source ~/.bashrc
geneid 

###Install genemarkv4.71
###you should register first at http://exon.gatech.edu/GeneMark/license_download.cgi
##wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_Q5vQH/gmes_linux_64_4.tar.gz
##wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_Q5vQH/gm_key_64.gz  
##tar -xzvf gmes_linux_64_4.tar.gz
##export PATH="/data/baojin/software/gmes_linux_64_4:$PATH" 
##You'll run into the problem of missing a lot of perl modules, and just install whichever one is missing！ You must try first!
##perl gmes_petap.pl
###conda install cpanminus
###cpanm YAML Hash::Merge Logger::Simple Parallel::ForkManager MCE::Mutex
perl gmes_petap.pl

###Install GeMoMa1.6.4
###conda install -c bioconda gemoma -y
##the path is at /data/tools/miniconda3/envs/gene_predict/share/gemoma-1.6.4-1
java -jar /data/tools/miniconda3/envs/gene_predict/share/gemoma-1.6.4-1/GeMoMa-1.6.4.jar

###Install genomeThreaderv1.7.3
##wget https://genomethreader.org/distributions/gth-1.7.3-Linux_x86_64-64bit.tar.gz
##tar -xzvf gth-1.7.3-Linux_x86_64-64bit.tar.gz
##Make sure the bssm folder, gthdata folder, and executable files are in the same directory.
##vim ~/.bashrc
##export PATH="/data/baojin/software/gth-1.7.3-Linux_x86_64-64bit/bin:$PATH"
##export BSSMDIR=/data/baojin/software/gth-1.7.3-Linux_x86_64-64bit/bin/bssm
##export GTHDATADIR=/data/baojin/software/gth-1.7.3-Linux_x86_64-64bit/bin/gthdata
gth

###Install Exonerate v2.4.0
##conda install -c bioconda exonerate -y
exonerate

###Install EvidenceModeler v2.1.0
##git clone https://github.com/EVidenceModeler/EVidenceModeler.git


###Install RepeatMsker v4.1.6
##in setonix, using singularity container, just singularity pull docker://pegi3s/repeat_masker
##conda install -c bioconda trf
##conda install -c bioconda hmmer==3.2.1
##conda install -c bioconda -c conda-forge rmblast

##wget https://repeatmasker.org/RepeatMasker/RepeatMasker-4.1.6.tar.gz
##tar -zxvf RepeatMasker-4.1.6.tar.gz

##wget https://www.biochen.org/public/software/RepBaseRepeatMaskerEdition-20181026.tar.gz 
##wget https://www.dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.0.h5.gz  
###tar -xzvf dfam38_full.0.h5.gz  ###uncompress to the path $RepeatMasker/Libraries/famdb

##conda install h5py
cd RepeatMasker-4.1.6
./configure
##根据提示安装,安装成功后加入环境
##export PATH="/data/baojin/software/RepeatMasker:$PATH"
#RepeatMasker S1_hifi.asm.bp.p_ctg.fa -species "Serradella" -e hmmer -xsmall -s -gff -pa 15  #到这里，发现没有这个物种，所以就要重头预测，应用RepeatModeler
#-e search engine,rmblast hmmer crossmatch abblast, anyone of them
#-species species name, must be in the NCBI Species classification database, Latin scientific names are recommended
#-s  Slow search; 0-5% more sensitive, 2-3 times slower than default
#-q  Quick search; 5-10% less sensitive, 2-5 times faster than default
#-nolow Does not mask low_complexity DNA or simple repeats
#-xsmall Returns repetitive regions in lowercase (rest capitals) rather than masked-----sm, soft-mask
#-pa(rallel) [number] The number of sequence batch jobs [50kb minimum] to run in parallel. RepeatMasker will fork off this number of parallel jobs, each running the search engine specified. For each search engine invocation ( where applicable ) a fixed the number of cores/threads


###Install RepeatModeler 2.0.5 (https://github.com/Dfam-consortium/RepeatModeler)
##RECON 1.0.8
conda install -c bioconda -c conda-forge RECON
##CD-HIT 4.8.1
conda install -c bioconda cd-hit
##RepeatScout
wget http://www.repeatmasker.org/RepeatScout-1.0.6.tar.gz
tar -xzvf RepeatScout-1.0.6.tar.gz 
cd RepeatScout-1.0.6/
make
##LtrHarvest---The LtrHarvest program is part of the GenomeTools suite
wget https://genometools.org/pub/binary_distributions/gt-1.6.2-Linux_x86_64-64bit-complete.tar.gz
tar -xzvf gt-1.6.2-Linux_x86_64-64bit-complete.tar.gz gt-1.6.2-Linux_x86_64-64bit-complete
##Ltr_retriever
conda install -c bioconda -c conda-forge ltr_retriever
git clone https://github.com/oushujun/LTR_retriever.git
##MAFFT - A multiple sequence alignment program.
wget https://mafft.cbrc.jp/alignment/software/mafft-7.505-with-extensions-src.tgz 
tar -xzvf mafft-7.505-with-extensions-src.tgz
vim /data/baojin/software/mafft-7.505-with-extensions/core/Makefile 
modify the first line 'PREFIX = /usr/local' to 'PREFIX = /data/baojin/software/mafft-7.505-with-extensions'  ##install path modify
cd /data/baojin/software/mafft-7.505-with-extensions/core/
make clean
make
make install

##NINJA
wget https://github.com/TravisWheelerLab/NINJA/archive/refs/tags/0.98-cluster_only.tar.gz
tar -zxvf 0.98-cluster_only.tar.gz
cd /data/baojin/software/NINJA-0.98-cluster_only/NINJA
make all
cd ..
cp -r NINJA/ /data/baojin/software

##UCSC TwoBit Tools
mkdir kent
cd kent
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo
chmod +x *

###the last, install RepeatModeler2.0.5
git clone https://github.com/Dfam-consortium/RepeatModeler.git
cd RepeatModeler
##install perl modle the server missing
cpan File::Which
cpan URI
#cpan Devel::Size ###can not install
#cpan LWP::UserAgent ###can not install, so i modified the configure file to delete the two dependence module, these two modules were useless.
根据提示设置上述软件的路径即可

If you use singularity containers, it will be easy, just using following command, this container including RepeatMasker, RepeatModeler, and coseg. 
singularity pull docker://dfam/tetools

####Run RepeatModeler
##First build library
mkdir db
BuildDatabase --name db/Serradalle1 S1_hifi.asm.bp.p_ctg.fa
BuildDatabase --name db/Serradalle2 S2_hifi.asm.bp.p_ctg.fa#--name is the name of library
##Error: Unknown argument: "blastdb_version"; here i have this error, maybe due to the blast version, so i modify the BuilDatabase scipt to remove this augument, then it works! 

##Second run RepeatModeler
RepeatModeler -database -database db/Serradalle1 -threads 128 -LTRStruct
RepeatModeler -database -database db/Serradalle2 -threads 128 -LTRStruct

##Third run RepeatMasker
RepeatMasker S1_hifi.asm.bp.p_ctg.fa -lib db/Serradalle1-families.fa -e hmmer -xsmall -s -gff -pa 128
RepeatMasker S2_hifi.asm.bp.p_ctg.fa -lib db/Serradella2-families.fa -e hmmer -xsmall -s -gff -pa 128



##Files preparation; nimbus server
#First step. make direction and reference
cd /data/baojin/01.Species_genome_assembly/03.gene_predict
mkdir augustus gemoma geneid gmes gth exo evm results
因为只知道该物种为豆科作物，所以我选了进化关系上近的物种豌豆为训练参考模型

#Reference path: ubuntu@node-1:/data/baojin/01.Species_genome_assembly/Reference
Pisum_sativum (pea) (https://plants.ensembl.org/Pisum_sativum/Info/Annotation/#assembly)
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/gff3/pisum_sativum/Pisum_sativum.Pisum_sativum_v1a.58.gff3.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/pisum_sativum/dna/Pisum_sativum.Pisum_sativum_v1a.dna.toplevel.fa.gz

##Medicago truncatula (Barrel Medic) (https://plants.ensembl.org/Medicago_truncatula/Info/Annotation/#assembly)
##wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/medicago_truncatula/dna/Medicago_truncatula.MedtrA17_4.0.dna.toplevel.fa.gz
##wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/gff3/medicago_truncatula/Medicago_truncatula.MedtrA17_4.0.58.chr.gff3.gz

##Glycine max (soybean) https://plants.ensembl.org/Glycine_max/Info/Annotation/#assembly
##wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/glycine_max/dna/Glycine_max.Glycine_max_v2.1.dna.toplevel.fa.gz
##wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/gff3/glycine_max/Glycine_max.Glycine_max_v2.1.58.chr.gff3.gz

#Second Step. mask reference repeat sequences
repeatmask


#get *rm* files

##############Third step. simplify sequence names. I did not use this setp.
#Sequence names that contain multiple characters will result in different software predictions, and we will deal with simplified sequence ids uniformly here
#fasta_id_sim.py

#{#!/usr/bin/env python
import sys
import pysam

with pysam.FastxFile(sys.argv[1]) as fh:
    for r in fh:
        new_name = r.name.split(' ')[0]
        print(">"+new_name)
        print(r.sequence)
}
#commands
ls *.fasta|cut -d "." -f1|while read id;do
fasta_id_sim.py ${id}.fasta >${id}_sim.fasta
done



####Geneome annotations 
##cd-hit 因尚未有转录组数据，我这里先找了四个同源物种的蛋白序列（豌豆，苜蓿，三叶草，大豆）合并为一个文件，用该软件去冗余,用于braker3的近缘蛋白输入文件

#蛋白使用以下参数
cd-hit -i seq.fasta -o seq-out.fasta -c 0.4 -T 4 -n 2
#-i 输入文件，fasta格式
#-o 输出文件
#-c 相似性,0.4代表相似度大于40%的为一类
#-n 两两序列进行序列比对时选择的 word size
#-T 使用的线程数
#Choose of word size:
#-n 5 for thresholds 0.7 ~ 1.0
#-n 4 for thresholds 0.6 ~ 0.7
#-n 3 for thresholds 0.5 ~ 0.6
#-n 2 for thresholds 0.4 ~ 0.5

#核苷酸用以下参数
#cd-hit -i seq.fasta -o seq-out.fasta -c 0.8 -T 4 -n 4
#-i 输入文件，fasta格式
#-o 输出文件
#-c 相似性,0.8代表相似度大于80%的为一类
#-T 使用的线程数
#-n 两两序列进行序列比对时选择的 word size
#Choose of word size:
#-n 10, 11 for thresholds 0.95 ~ 1.0
#-n 8,9 for thresholds 0.90 ~ 0.95
#-n 7 for thresholds 0.88 ~ 0.9
#-n 6 for thresholds 0.85 ~ 0.88
#-n 5 for thresholds 0.80 ~ 0.85
#-n 4 for thresholds 0.75 ~ 0.8






###AUGUSTUS
##

perl /data/tools/miniconda3/pkgs/augustus-3.2.2-0/scripts/autoAugTrain.pl --genome=Pisum_sativum.Pisum_sativum_v1a.dna_rm.toplevel.fa --trainingset=Pisum_sativum.Pisum_sativum_v1a.58.chr.gff3 --species=Pea


##Braker3
braker.pl   Pipeline for predicting genes with GeneMark-EX and AUGUSTUS with RNA-Seq and/or proteins

srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/braker3.sif braker.pl --genome /scratch/pawsey0399/bguo1/0.assembly/01.hifi_assembly/S1_HIFI_RESULT/S1_hifi.asm.bp.p_ctg.fa.masked --species=Ser1 --prot_seq /scratch/pawsey0399/bguo1/0.assembly/03.gene_annotation/Ref_pep/allhomo.pep.cdhit.fa --useexisting --threads 128 \
--workingdir=/scratch/pawsey0399/bguo1/0.assembly/03.gene_annotation/1.Braker/S1/
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/braker3.sif braker.pl --genome /scratch/pawsey0399/bguo1/0.assembly/01.hifi_assembly/S2_HIFI_RESULT/S1_hifi.asm.bp.p_ctg.fa.masked --species=Ser2 --prot_seq /scratch/pawsey0399/bguo1/0.assembly/03.gene_annotation/Ref_pep/allhomo.pep.cdhit.fa --useexisting --threads 128 \
--workingdir=/scratch/pawsey0399/bguo1/0.assembly/03.gene_annotation/1.Braker/S2/

