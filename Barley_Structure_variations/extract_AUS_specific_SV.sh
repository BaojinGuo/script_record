#This script is used to extract specific information from a VCF file based on the genotypes of samples belonging to two groups ('Aus' and 'Other'). 
#It first filters the VCF file to select samples from these two groups, then calculates genotype counts for each group at each site. 
#It selects sites where the genotype counts meet certain criteria and extracts their positions. 
#Finally, it uses the position information to extract complete data from the original VCF file and saves it to a new VCF file named 'AUS_special_SV.vcf'.



#!/bin/bash

# Define group lists
aus_group=(
    "Compass" "Leabrook" "Fathom" "Clipper" "Prior" "Buloke" "Schooner" "Yambla" 
    "Galleon" "Flagship" "Vlamingh" "Stirling" "Baudin" "LaTrobe" "Maximus" 
    "Buff" "Mundah" "RGT_Planet" "Gairdner" "Halcyon"
)
other_group=(
    "Yeti" "Laperouse" "Eth69" "Beecher" "Sahara3771" "ILAN15" "Betzes" "Harrington"
)

# Convert group lists to comma-separated strings
aus_samples=$(IFS=, ; echo "${aus_group[*]}")
other_samples=$(IFS=, ; echo "${other_group[*]}")

# Define input and output files, modify this
input_vcf="combineSV_Clipper_sniffle.vcf"
output_vcf="AUS_special_SV.vcf"

# Use bcftools to extract sample information for Aus and Other groups, and filter specific sites of AUS group with the standard AUS 1/1 >15 & Other 0/0 >5 or AUS 0/0 >15 & Other 1/1 >5. 
bcftools view -s $aus_samples,$other_samples $input_vcf -Ou | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | awk -v aus_list="$aus_samples" -v other_list="$other_samples" '
BEGIN {
    split(aus_list, aus, ",");
    split(other_list, other, ",");
}
{
    aus_1_1 = 0;
    aus_0_0 = 0;
    other_0_0 = 0;
    other_1_1 = 0;
    for (i = 1; i <= length(aus); i++) {
        if ($(i+4) == "1/1") aus_1_1++;
        if ($(i+4) == "0/0") aus_0_0++;
    }
    for (i = length(aus) + 1; i <= length(aus) + length(other); i++) {
        if ($(i+4) == "0/0") other_0_0++;
        if ($(i+4) == "1/1") other_1_1++;
    }
    if ((aus_1_1 > 15 && other_0_0 > 5) || (aus_0_0 > 15 && other_1_1 > 5)) {
        print $0
    }
}' > inter.vcf

# Extract positions from intermediate VCF file
cut -f1,2 inter.vcf > position.txt

# Step 2: Extract complete data from the original file based on position information and save to the output file
# Use vcftools to extract lines matching the positions from the input VCF file and save to the output VCF file
vcftools --vcf $input_vcf --positions position.txt --recode --recode-INFO-all --stdout > $output_vcf

# Print completion message
echo "Extraction complete. Results saved to $output_vcf"
