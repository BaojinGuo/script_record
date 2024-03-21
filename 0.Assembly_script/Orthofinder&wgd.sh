##Orthofinder
singularity  pull docker://davidemms/orthofinder:latest
##simplify all protein fasta file name as *.fa and put all pep file to the same file "pep/", delete charactor '.' and '*' in fasta.
ls | while read file; do
    if [ -f "$file" ]; then
        prefix=$(echo "$file" | cut -d '.' -f 1)
        new_name="${prefix}.fa"
        mv "$file" "$new_name"
        echo "Renamed $file to $new_name"
    fi
done
orthofinder -f pep/ -t 128 -a 128

##wgd, python==3.8
conda install -c bioconda wgd
pip install numpy==1.19.0 
git clone https://github.com/heche-psb/wgd.git
pip install -r requirement.txt
pip install .

wgd dmd Arabidopsis_thaliana.cds.fa Cajanus_cajan.cds.fa Cicer_arietinum.cds.fa Cicer_echinospermum.cds.fa Cicer_reticulatum.cds.fa Glycine_max.cds.fa Glycine_soja.cds.fa Lotus_japonicus.cds.fa Lupinus_angustifolius.cds.fa Medicago_truncatula.cds.fa Phaseolus_vulgaris.cds.fa Pisum_sativum.cds.fa S1.EVM.cds.fa S2.EVM.cds.fa Trifolium_pratense.cds.fa Vigna_angularis.cds.fa Vigna_radiata.cds.fa Vigna_unguiculata.cds.fa -e 1e-10 -n 128 -o ../../0.workdic/fablas.dmd


