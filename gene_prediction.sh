########server is nimbus

####Gene Prediction -- do novo prediction, homology-based prediction, transcriptome-based prediction
###Install the packages -- AUGUSTUSv3.2.2
##conda install -c bioconda augustus==3.2.2
###error while loading shared libraries: libboost_iostreams.so.1.60.0: solved sudo ln -s libboost_iostreams.so.1.65.1 libboost_iostreams.so.1.60.0 (https://www.biostars.org/p/389360/)


###Install geneidv1.4
##git clone https://github.com/guigolab/geneid
##make
##vim ~/.bashrc
##export PATH="/data/baojin/software/geneid/bin:$PATH"
##source ~/.bashrc


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


###Install GeMoMa1.6.4
###conda install -c bioconda gemoma -y
##the path is at /data/tools/miniconda3/envs/gene_predict/share/gemoma-1.6.4-1


###Install genomeThreaderv1.7.3
##wget https://genomethreader.org/distributions/gth-1.7.3-Linux_x86_64-64bit.tar.gz
##tar -xzvf gth-1.7.3-Linux_x86_64-64bit.tar.gz
##Make sure the bssm folder, gthdata folder, and executable files are in the same directory.
##vim ~/.bashrc
##export PATH="/data/baojin/software/gth-1.7.3-Linux_x86_64-64bit/bin:$PATH"
##export BSSMDIR=/data/baojin/software/gth-1.7.3-Linux_x86_64-64bit/bin/bssm
##export GTHDATADIR=/data/baojin/software/gth-1.7.3-Linux_x86_64-64bit/bin/gthdata


###Install Exonerate v2.4.0
##conda install -c bioconda exonerate -y


###Install EvidenceModeler
##git clone https://github.com/EVidenceModeler/EVidenceModeler.git






