
#first step, convert heterozygous to missing genotype
python convert_vcf_het_to_miss.py combineSV_Clipper_sniffle.vcf 1.combineSV_Clipper_sniffle_het_to_miss.vcf
#second step, extarct vcf file with parameter "maf 2.5% and missing below 20%"
vcftools 	--vcf 1.combineSV_Clipper_sniffle_het_to_miss.vcf
	--recode-INFO-all
	--maf 0.025
	--max-missing 0.8
	--out 2.combineSV_Clipper_sniffle_het_to_miss80_maf25
	--recode
#third step, divide to each vcf per sample
mkdir 3.split_sample
grep CHROM 2.combineSV_Clipper_sniffle_het_to_miss80_maf25.recode.vcf |cut -f10- |sed 's/\t/\n/g'|while read line
do bcftools view -s $line 2.combineSV_Clipper_sniffle_het_to_miss80_maf25.recode.vcf >3.split_sample/$line.vcf
done
#forth step, extract genotype which is '1/1'
cd 3.split_sample
ls |while read line
do
bcftools view -i 'GT="1/1"' $line >$line.11.vcf
done
#fifth step, subset information for plot per sample
ls *.vcf.11.vcf|cut -f1 -d '.'|while read line
do grep -v "##" $line.vcf.11.vcf |cut -f1-3|grep -v "BND"|sed 's/Sniffles2\.//g'|sed 's/\..*//g'|sed 's/\#CHROM/Chromosome/g'|sed 's/ID/Type/g'|sed 's/POS/Position/g' >$line.Clipper.hte_to_miss80_maf25.SV.txt
done
##sixth step, add accession information
ls *.txt|cut -f1 -d "."|while read line
do sed -i 's/^/'$line'\t/g' $line.Clipper.hte_to_miss80_maf25.SV.txt
sed -i '1s/'$line'/Accession/g' $line.Clipper.hte_to_miss80_maf25.SV.txt
done
##seventh step, combine all files and add header
cat *.txt |grep -v Accession|sed '1i Accession\tChromosome\tPosition\tType' >Barley.hte_to_miss80_maf25.allSV.txt
##eighth step, delete lines Clipper have
library(dplyr)
data<-read.table("Barley.hte_to_miss80_maf25.allSV.txt",header = T)
clipper_data <- filter(data, Accession == "Clipper")
filtered_data <- data %>%
anti_join(clipper_data, by = c("Chromosome", "Position", "Type"))
write.table(filtered_data, "Barley.hte_to_miss80_maf25.allSV_noClipper.txt", sep = "\t", row.names = FALSE, quote = FALSE)
######ninth step, statistic SV length
bcftools query -f '%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\n' 2.combineSV_Clipper_sniffle_het_to_miss80_maf25.recode.vcf |grep -v BND >4.combineSV_Clipper_sniffle_het_to_miss80_maf25.svlength



