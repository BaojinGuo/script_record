#!/bin/bash

###############################################################################
# Script Name: extract_minor_allele_sites_CM72.sh
#
# Description:
#   This script extracts all variant sites in sample CM72 with 
#   minor allele frequency (MAF) < 5%, where the allele carried 
#   by CM72 is different from the major allele. 
#   The analysis is based on the reference genome MorexV3.
#
#   Steps include:
#     1. Filtering low MAF variants
#     2. Calculating allele frequencies
#     3. Extracting CM72 genotypes and filtering missing/heterozygous calls
#     4. Standardizing variant IDs
#     5. Identifying major alleles
#     6. Normalizing file delimiters (space to tab)
#     7. Merging genotype with major allele info
#     8. Determining actual CM72 allele
#     9. Selecting variants where CM72 allele ≠ major allele
#     10. Merging results across chromosomes
#
# Author: Baojin
# Date: 2025-08-04
###############################################################################
# Step 1: Extract variants with MAF < 0.05 from the original dataset
for i in {1..7}H; do vcftools --gzvcf AUS_100+_Morexv3.${i}.cohort.vcf.gz --recode-INFO-all --chr ${i} --max-maf 0.05 --out AUS_MorexV3_${i}.lessMAF05 --recode; done

# Step 2: Calculate allele frequency for each variant
for i in {1..7}H; do vcftools --gzvcf AUS_MorexV3_${i}.lessMAF05.vcf.gz --chr${i} --freq --out ${i}.allele_freq; done

# Step 3: Extract genotype of the target sample (H008), and remove missing and heterozygous genotypes
for i in {1..7}H; do bcftools view -s H008 AUS_MorexV3_${i}.lessMAF05.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ID\t[%GT]\n' | grep -v '0/1' | grep -v '\./\.' > try/H008_${i}.genotype.txt & done

# Step 4: Standardize the variant ID (keep only the part before the semicolon)
for i in {1..7}H; do awk 'BEGIN{OFS="\t"} {split($3,a,";"); $3=a[1]; print}' H008_${i}.genotype.txt  > tmp && mv tmp H008_${i}.genotype.txt ; done

# Step 5: Identify the major allele for each variant
for i in {1..7}H; do awk 'NR==1 {
        print $1, $2, $3, $4, "MAX_ALLELE"; next
    }
    {
        max_freq = -1
        max_allele = ""
        for (i = 5; i <= NF; i++) {
            split($i, a, ":")
            if (a[2] + 0 > max_freq) {
                max_freq = a[2] + 0
                max_allele = a[1]
            }
        }
        print $1, $2, $3, $4, max_allele
    }' ${i}.allele_freq.frq | cut -f1,2,5 -d" " > ${i}.major_allele & done

# Step 6: Normalize file format: replace multiple spaces with tab delimiters
for i in {1..7}H; do sed -i 's/ \+/\t/g' ${i}.major_allele & done

# Step 7: Merge major allele information with sample genotype
for i in {1..7}H; do awk '
        NR==FNR && FNR > 1 { key = $1 FS $2; allele[key] = $3; next }
        {
            key = $1 FS $2
            print $0, (key in allele ? allele[key] : "NA")
        }
    ' ${i}.major_allele H008_${i}.genotype.txt > ${i}.merge.txt & done

# Step 8: Normalize merged file format (space to tab)
for i in {1..7}H; do sed -i 's/ \+/\t/g' ${i}.merge.txt & done

# Step 9: Compare genotype with major allele, and determine actual allele of the sample
for i in {1..7}H; do awk '
    BEGIN { OFS="\t" }
    NR==1 {
        print $0, "ALLELE"; next
    }
    {
        split($3, a, "_");  # e.g., ID = 1H_322_T_TA → a[3]=T, a[4]=TA
        gt = $4
        if (gt == "0/0") allele = a[3]
        else if (gt == "1/1") allele = a[4]
        else if (gt == "0/1" || gt == "1/0") allele = a[3] "/" a[4]
        else allele = "NA"
        print $0, allele
    }
    ' ${i}.merge.txt > ${i}.merge.geno.txt & done

# Step 10: Retain only the sites where the sample's allele is not the major allele
for i in {1..7}H; do awk '$5 != $6' ${i}.merge.geno.txt > ${i}.diff & done

# Final filtering: remove lines with NA alleles and combine all chromosomes into one file
for i in {1..7}H; do grep -v "NA" ${i}.diff > ${i}.diff.txt; done
cat *diff.txt >> CM72.MAF05.diff.site
