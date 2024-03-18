###Gene annotation assessment, two methods, busco and OMArk
##OMArk, using https://omark.omabrowser.org/ directly or using server
#Installation
conda install -c bioconda omark
singulairty pull omark.sif docker://phuongdoan96/omark:0.3.0

##Download database
wget https://omabrowser.org/All/LUCA.h5
##First: run OMAmer on the proteome FASTA
omamer search -d Ref/LUCA.h5 -q 0.assembly/03.gene_annotation/5.EVM/S1/S1.EVM.pep --nthreads 128 --out 0.assembly/03.gene_annotation/6.OMArk/S1/S1.EVM.omamer
omamer search -d Ref/LUCA.h5 -q 0.assembly/03.gene_annotation/5.EVM/S2/S2.EVM.pep --nthreads 128 --out 0.assembly/03.gene_annotation/6.OMArk/S2/S2.EVM.omamer
##Then, use OMArk
omark -f 0.assembly/03.gene_annotation/6.OMArk/S1/S1.EVM.omamer -d Ref/LUCA.h5 -o 0.assembly/03.gene_annotation/6.OMArk/S1/S1.omark_output
omark -f 0.assembly/03.gene_annotation/6.OMArk/S2/S2.EVM.omamer -d Ref/LUCA.h5 -o 0.assembly/03.gene_annotation/6.OMArk/S2/S2.omark_output

##Visualising multiple OMArk results
#plot_all_results.py -i 0.assembly/03.gene_annotation/6.OMArk/S1/S1.omark_output -o 0.assembly/03.gene_annotation/6.OMArk/S1/S1.omark_output/S1.fig.png
#plot_all_results.py -i 0.assembly/03.gene_annotation/6.OMArk/S2/S2.omark_output -o 0.assembly/03.gene_annotation/6.OMArk/S2/S2.omark_output/S2.fig.png

##BUSCO
busco -i ../5.EVM/S1/S1.EVM.pep -m prot -o S1_annotation_gene_busco_fabales -l ../../02.quality_assessmnet/03.BUSCO/fabales_odb10 --cpu 128 --offline
busco -i ../5.EVM/S2/S2.EVM.pep -m prot -o S2_annotation_gene_busco_fabales -l ../../02.quality_assessmnet/03.BUSCO/fabales_odb10 --cpu 128 --offline
busco -i ../5.EVM/S1/S1.EVM.pep -m prot -o S1_annotation_gene_busco_embryophyta -l ../../02.quality_assessmnet/03.BUSCO/embryophyta_odb10 --cpu 128 --offline
busco -i ../5.EVM/S2/S2.EVM.pep -m prot -o S2_annotation_gene_busco_embryophyta -l ../../02.quality_assessmnet/03.BUSCO/embryophyta_odb10 --cpu 128 --offline
busco -i ../5.EVM/S1/S1.EVM.pep -m prot -o S1_annotation_gene_busco_eukaryota -l ../../02.quality_assessmnet/03.BUSCO/eukaryota_odb10 --cpu 128 --offline
busco -i ../5.EVM/S2/S2.EVM.pep -m prot -o S2_annotation_gene_busco_eukaryota -l ../../02.quality_assessmnet/03.BUSCO/eukaryota_odb10 --cpu 128 --offline
