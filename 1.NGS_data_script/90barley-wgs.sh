







vcftools --gzvcf barley_1H.cohort.vcf.gz --remove-indels --maf 0.25 --max-missing 0.8 --recode --recode-INFO-all --out barley_1H.snp.filter

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' barley_1H.snp.filter.cohort.vcf.gz.recode.vcf |awk -v OFS="\t" '{het=0;hom=0;for(i=5;i<=NF;i++) if($i=="0/1" || $i=="1/0") het+=1; else hom+=1;} {if(het/(het+hom)<0.1) print $1,$2;}' > low_het_sites.txt

bcftools view -T low_het_sites.txt barley_1H.snp.filter.cohort.vcf.gz.recode.vcf -Oz -o 1H_filtered.vcf.gz
