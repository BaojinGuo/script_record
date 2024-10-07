#####geneHapR is for gene haplotype statistics, phenotype association and visualization.
#####http://127.0.0.1:28913/library/geneHapR/doc/Introduction_of_geneHapR.html

install.packages("geneHapR")

library(geneHapR)

# import vcf file
vcf <- import_vcf("RAP_MAF5_U80.vcf.recode.vcf")

###import annotation files
# import GFFs
gff <- import_gff("your_gff_file_path.gff", format = "GFF")

# import GFFs
bed <- import_bed("RAP.gff.bed")


# import phynotype and accession group information data
pheno <- import_AccINFO("pheno.txt")

head(pheno)
#       Subpopulation Longitude Latitude Grain length
#NLL001  intermediate        30       10            1
#NLL002      japonica       120       50           10
#NLL003      japonica        30       10            5
#NLL004        indica        60       10            7
#NLL005  intermediate        30       10            3
#NLL006        indica        60       30            3

# filter VCF by position and annotation
# gene as a unit
vcf_f2 <- filter_vcf(vcf, mode = "both",
                     Chr = "HiC_scaffold_6",
                     start = 21051814, end = 21054935,
                     gff = bed,
                     type="gene")
# cds as a unit
vcf_f3 <- filter_vcf(vcf, mode = "both",
                     Chr = "HiC_scaffold_6",
                     start = 21051814, end = 21054935,
                     gff = bed,
                     type = "CDS")

# haplotype calculate
hapResult <- vcf2hap(vcf_f3,
                     hapPrefix = "H",
                     hetero_remove = TRUE,
                     na_drop = TRUE)
# set position of ATG as zero in hapResult/hapSummary
newhap <- hapSetATGas0(gff = bed, hap = hapResult,
                       geneID = "OIV89004_R0",
                       Chr = "HiC_scaffold_6",
                       POS = c(21051813,21054935))


# Summary hapResult
hapSummary_raw<-hap_summary(hapResult) #raw result
hapSummary <- hap_summary(newhap) #set position as zero result

# Visualize haplotye as table
pdf(file = "01.hap_as_table.pdf",height = 8,width = 12)
plotHapTable(hapSummary)
dev.off()
# Display variations on gene model.
pdf(file = "02.hap_structure.pdf",height = 4,width = 12)
displayVarOnGeneModel(hapSummary_raw, bed,
                      Chr = "HiC_scaffold_6",
                      startPOS = 21051813, endPOS = 21054935,
                      type = "pin", cex = 0.7,
                      CDS_h = 0.05, fiveUTR_h = 0.02, threeUTR_h = 0.01)
dev.pff()

# hapNet calculation and visualization
hapNet <- get_hapNet(hapSummary,
                     AccINFO = pheno,
                     groupName = "Subpopulation")
# plot haploNet
pdf(file = "03.hap_net_plot.pdf",height = 12,width = 6)

plotHapNet(hapNet,
           size = "freq",                   # circle size
           scale = "log10",                 # scale circle with 'log10(size + 1)'
           cex = 0.5,                       # size of hap symbol
           col.link = 2,                    # link colors
           link.width = 1,                  # link widths
           show.mutation = 2,               # mutation types one of c(0,1,2,3)
           legend = c(-60, 20),        # legend position
           pie.lim = c(10,20),
           labels.cex = 0.8,legend_version = 1,show_color_legend = T,show_size_legend = T)
dev.off()

#Geography distribution of main haplotypes

pdf(file = "04.hap_world_distribution_plot.pdf",height = 6,width = 10)
pheno$Longitude <- as.numeric(pheno$Longitude)
pheno$Latitude <- as.numeric(pheno$Latitude)
hapDistribution(hapResult,
                AccINFO = pheno,
                LON.col = "Longitude",
                LAT.col = "Latitude", 
                hapNames = c("H001", "H002", "H003","H004","H005","H006","H007","H008","H009","H010","H011","H012","H013"), 
                legend = TRUE,label.cex = 1)
dev.off()

##Phenotype association analysis
#merge fig
results <-hapVsPheno(hapResult,
                     hapPrefix = "H",
                     title = "",
                     mergeFigs = TRUE,
                     pheno = pheno,
                     phenoName = "Grain length",
                     minAcc = 15) ##the min accession is set to 15

pdf(file = "05.hap_pheno.pdf",height = 8,width = 12)

plot(results$figs)
dev.off()

#seperate fig
results2 <-hapVsPheno(hapResult,
                     hapPrefix = "H",
                     title = "",
                     mergeFigs = FALSE,
                     pheno = pheno,
                     phenoName = "Grain length",
                     minAcc = 15)

pdf(file = "05-1.hap_pheno_pvalue.pdf",height = 6,width = 6)
plot(results2$fig_pvalue)
dev.off()

pdf(file = "05-2.hap_pheno_violin.pdf",height = 6,width = 6)
plot(results2$fig_Violin)
dev.off()



















