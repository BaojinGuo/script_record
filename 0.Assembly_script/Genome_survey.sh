####jellyfish+genome scope
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.1/jellyfish-2.3.1.tar.gz
git clone https://github.com/tbenavi1/genomescope2.0.git

jellyfish count -t 128 -C -m 21 -o s1_21mer_out -s 64G s1_hifi.fastq
jellyfish histo -o S1_21mer/s1_21mer_out.histo 21mer_out
jellyfish count -t 128 -C -m 21 -o s2_21mer_out -s 64G s2_hifi.fastq 
jellyfish histo -o S2_21mer/s2_21mer_out.histo s2_21mer_out
jellyfish count -t 128 -C -m 21 -o s1_31mer_out -s 64G s1_hifi.fastq
jellyfish histo -o S1_31mer/s1_31mer_out.histo s1_31mer_out
jellyfish count -t 128 -C -m 21 -o s2_21mer_out -s 64G s2_hifi.fastq
jellyfish histo -o S2_31mer/s2_31mer_out.histo s2_31mer_out

Rscript genomescope.R s2_21mer_out.histo 21 15000 ./
Rscript genomescope.R s2_21mer_out.histo 21 15000 ./
Rscript genomescope.R s1_31mer_out.histo 31 15000 ./
Rscript genomescope.R s2_31mer_out.histo 31 15000 ./

