conda create -n gapit -c bioconda -c conda-forge -c hcc r-devtools r-lme4 r-base=4.4.1 r-ggplot2 r-dplyr r-data.table r-readr r-tibble r-tidyr
devtools::install_github("jiabowang/GAPIT3", build_vignettes = FALSE)
library(GAPIT)


hmp_folder <- "/Users/baojinguo/Downloads/Generegion/hmp"
pheno_file <- "/Users/baojinguo/Downloads/487Oat/oat_ph.txt"

output_base <- "/Users/baojinguo/Downloads/487Oat/gR"
phenotype <- read.table(pheno_file, header = TRUE)
hmp_files <- list.files(hmp_folder, pattern = "\\.hmp\\.txt$", full.names = TRUE)

for (hmp_path in hmp_files) {
    
    # 提取染色体名称作为输出文件夹，例如 1A1 从 Oat_Adata_Generegion_U10_MAF5_1A1.hmp.txt
    file_name <- basename(hmp_path)
    chr_name <- sub(".*_([0-9A-Za-z]+)\\.hmp\\.txt$", "\\1", file_name)
    
    # 输出路径
    out_dir <- file.path(output_base, chr_name)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    
    # 设置工作目录为当前输出文件夹
    setwd(out_dir)
    
    cat(">> 运行 GAPIT for", chr_name, "\n")
    myG<-read.table(hmp_path)
    # 运行 GAPIT
    myGAPIT <- GAPIT(Y = phenotype,G = myG,model = c("MLM","MLMM","SUPER","FarmCPU","Blink"),PCA.total = 3)
    # 返回主输出路径
    setwd(output_base)
}
