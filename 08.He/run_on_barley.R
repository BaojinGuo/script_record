#!/usr/bin/env Rscript
library(optparse)

invocation <- commandArgs()

usage <- "
Does a local PCA analysis on barley data.
To run an analysis on windows of 1000 SNPs, for instance:
    Rscript run_on_barley.R -t snp -s 1000
while to run it on windows of 1000 bps:
    Rscript run_on_barley.R -t bp -s 1000

Do
    Rscript run_on_barley.R -h
for options.
"

option_list <- list(
    make_option(c("-t", "--type"), type="character", help="Window by SNP or by bp?"),
    make_option(c("-s", "--size"), type="integer", help="Size of the window, in units of type."),
    make_option(c("-k", "--npc"), type="integer", default=2L, help="Number of principal components to compute for each window. [default: %default]"),
    make_option(c("-o", "--outdir"), type="character", help="Directory to save results to. [default: lostruct/results_type_%type_size_%size_jobid_%jobid/]"),
    make_option(c("-j", "--jobid"), type="character", default=formatC(1e6*runif(1), width=6, format="d", flag="0"), help="Unique job id. [default random]")
)

opt <- parse_args(OptionParser(option_list=option_list, description=usage))
if (is.null(opt$outdir)) {
    opt$outdir <- file.path("lostruct_results", 
                            sprintf("type_%s_size_%d_jobid_%s", 
                                    opt$type, 
                                    opt$size, 
                                    opt$jobid))
}

if (is.null(opt$type) || is.null(opt$size)) {
    stop(usage)
}

opt$start.time <- Sys.time()
opt$run.dir <- normalizePath(".")
print(opt)

# Load necessary libraries
library(lostruct)

# Define BCF files
chroms <- paste0(1:7, "H")
bcf.files <- file.path("Local_PCA_input", paste0("WGRS_Clipper.", chroms, ".AUS.DP10.snps.bcf"))
names(bcf.files) <- chroms

# Create output directory
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

# Perform local PCA by chromosome
all.pcas <- numeric(0)
all.lengths <- numeric(0)
all.regions <- data.frame()

for (bcf.file in bcf.files) {
    pca.file <- file.path(opt$outdir, sprintf(gsub(".bcf", ".pca.csv", basename(bcf.file))))
    regions.file <- file.path(opt$outdir, sprintf(gsub(".bcf", ".regions.csv", basename(bcf.file))))
    
    if (file.exists(pca.file)) {
        warning(paste("File", pca.file, "already exists! Not recomputing."))
        pca.stuff <- as.matrix(data.table::fread(pca.file, header=TRUE))
        these.regions <- data.table::fread(regions.file, header=TRUE)
    } else {
        cat("Finding PCs for", bcf.file, "and writing out to", pca.file, "and", regions.file, "\n")
        win.fn <- vcf_windower(bcf.file, size=opt$size, type=tolower(opt$type))
        these.regions <- region(win.fn)()
        system.time(
            pca.stuff <- eigen_windows(win.fn, k=opt$npc)
        )
        write.csv(pca.stuff, file=pca.file, row.names=FALSE)
        write.csv(these.regions, file=regions.file, row.names=FALSE)
    }
    
    all.pcas <- rbind(all.pcas, pca.stuff)
    all.lengths <- c(all.lengths, nrow(pca.stuff))
    all.regions <- rbind(all.regions, these.regions)
}

names(all.lengths) <- chroms

# Calculate distance matrix
cat("Done finding PCs, computing distances.\n")
system.time(pc.distmat <- pc_dist(all.pcas, npc=opt$npc))

# Perform MDS
mds.file <- file.path(opt$outdir, "mds_coords.csv")
cat("Done computing distances, running MDS and writing results to", mds.file, "\n")
na.inds <- is.na(all.pcas[, 1])
mds.coords <- cbind(
    data.frame(
        chrom=rep(chroms, all.lengths),
        window=unlist(lapply(all.lengths, seq_len))
    ),
    cmdscale(pc.distmat[!na.inds, !na.inds], k=opt$npc)[ifelse(na.inds, NA, cumsum(!na.inds)),]
)
colnames(mds.coords)[-(1:2)] <- paste0("MDS", seq_len(opt$npc))
write.csv(mds.coords, mds.file, row.names=FALSE)

cat("All done!\n")
Sys.time()
