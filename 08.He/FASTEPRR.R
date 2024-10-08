#!/programs/bin/R/R

######### R parallelized wrapper for FastEPRR #######################

# You need: vcf files imputed and phased per chromosome
# This wrapper takes 4 arguments:
# 1) vcf_prefix             --> if vcf files look like this chr1.vcf.gz, chr2.vcf.gz...   then vcf_prefix = chr
# 2) windows size           --> rho will be calculated in windows with this size
# 3) number of chromosomes  
#
# Usage example: Rscript FASTEPRR.R WGRS_Clipper. 1000000 7
#
# Script optimized to work with 128 cores to change this modify line 23, 45, and 46
# Three folders will be created: Step1_output, Step2_output, Step3_output
#
#####################################################################

args = commandArgs(trailingOnly=TRUE)

library(FastEPRR)
library(foreach)
library(doMC)
registerDoMC(128)

### STEP 1

phased_prefix = args[1]
window_size = args[2]
nchr = as.numeric(args[3])

# Create the output folder for Step 1
system("mkdir -p Step1_output")

foreach(i=1:nchr) %dopar% {
    path = system2("pwd", stdout=TRUE)
    input = sprintf("%s/%s%dH.AUS.DP10.snps.vcf.gz", path, phased_prefix, i)
    step1_output = sprintf("Step1_output/%dH.coor", i)

    FastEPRR_VCF_step1(vcfFilePath=input, erStart=1, winLength=1000000, srcOutputFilePath=step1_output)
}

### STEP 2
system("mkdir -p Step2_output")

foreach(I=1:128) %dopar% {
    FastEPRR_VCF_step2(srcFolderPath="Step1_output", jobNumber=128, currJob=I, DXOutputFolderPath="Step2_output")
}

### STEP 3
system("mkdir -p Step3_output")
FastEPRR_VCF_step3(srcFolderPath="Step1_output", DXFolderPath="Step2_output", finalOutputFolderPath="Step3_output")
