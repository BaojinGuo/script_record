###Download databases & software
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz ##swissprot
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz ##treml
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz ##NR
##eggnog
wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz
wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz
wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz 
wget http://eggnog5.embl.de/download/emapperdb-5.0.2/pfam.tar.gz
##interproscan
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.66-98.0/interproscan-5.66-98.0-64-bit.tar.gz
##diamond 
singularity pull diamond.sif docker://buchfink/diamond:version2.0.13

###Build diamond index
diamond makedb --in uniprot_sprot.fasta.gz --db uniprot_sprot --threads 128 
diamond makedb --in uniprot_trembl.fasta.gz --db uniprot_thrembl --threads 128 
diamond makedb --in nr.gz --db NR --threads 128

###blastp alignment

diamond blastp -d ../../Gene_Function_annotation_database/uniprot_sprot.dmnd -q ../03.gene_annotation/5.EVM/S1/S1.EVM.pep --evalue 1e-5 -f 6 -p 128 -o 2.Swissprot/S1.pep.swissprot.outfmt6 -k 1 --header
diamond blastp -d ../../Gene_Function_annotation_database/uniprot_thrembl.dmnd -q ../03.gene_annotation/5.EVM/S1/S1.EVM.pep --evalue 1e-5 -f 6 -p 128 -o 3.Trembl/S1.pep.trembl.outfmt6 -k 1 --header
diamond blastp -d ../../Gene_Function_annotation_database/uniprot_sprot.dmnd -q ../03.gene_annotation/5.EVM/S2/S2.EVM.pep --evalue 1e-5 -f 6 -p 128 -o 2.Swissprot/S2.pep.swissprot.outfmt6 -k 1 --header
diamond blastp -d ../../Gene_Function_annotation_database/uniprot_thrembl.dmnd -q ../03.gene_annotation/5.EVM/S2/S2.EVM.pep --evalue 1e-5 -f 6 -p 128 -o 3.Trembl/S2.pep.trembl.outfmt6 -k 1 --header
diamond blastp -d ../../Gene_Function_annotation_database/NR.dmnd -q ../03.gene_annotation/5.EVM/S1/S1.EVM.pep --evalue 1e-5 -f 6 -p 128 -o 1.NR/S1.pep.nr.outfmt6 -k 1 --header
diamond blastp -d ../../Gene_Function_annotation_database/NR.dmnd -q ../03.gene_annotation/5.EVM/S2/S2.EVM.pep --evalue 1e-5 -f 6 -p 128 -o 1.NR/S2.pep.nr.outfmt6 -k 1 --header

##interproscan, modify interproscan.sh java '-XX:ParallelGCThreads=8 -Xms2028M -Xmx9216M' to 'java -XX:ParallelGCThreads=128 -Xms9216M -Xmx512G' based on the server ability
bash /scratch/pawsey0399/bguo1/software/interproscan-5.66-98.0/interproscan.sh -i 5.Interpro/S1.EVM.pep -b 5.Interpro/S1.out.iprscan -goterms -iprlookup -pa -cpu 128 
bash /scratch/pawsey0399/bguo1/software/interproscan-5.66-98.0/interproscan.sh -i 5.Interpro/S2.EVM.pep -b 5.Interpro/S2.out.iprscan -goterms -iprlookup -pa -cpu 128 

