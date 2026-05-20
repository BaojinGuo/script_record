grep -Ei -l "error|fail|exception" slurm-*
comm -23 <(sort ../../mapping.txt ) <(ls *7A*.g.vcf.gz |sed 's#.*/##' | cut -d'.' -f1 | sort -u) >../7Adeep.missing

tree -d
.
├── 01.GWAS
├── 02.KGWAS
├── 03.Depth
├── Clean_data
├── Log
│   └── fastp
├── Raw_data
│   └── shell
└── Ref

REF=/scratch/pawsey0399/bguo1/Murdoch/11.Oat/Pinyan/PY6.fa
##########STEP1#############
while read y_id; do echo '#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
module load singularity/4.1.0-nompi
srun --export=all -n 1 -c 16 singularity exec /scratch/pawsey0399/bguo1/Singularity_image/fastp.sif fastp -i Raw_data/'${y_id}'_f1.fq.gz -I Raw_data/'${y_id}'_r2.fq.gz -o Clean_data/'${y_id}'_clean_1.fq.gz -O Clean_data/'${y_id}'_clean_2.fq.gz -h Log/fastp/'${y_id}'.fastp.html -j Log/fastp/'${y_id}'.fastp.json -w 16' > ${y_id}.fastp.sh
done <mapping.txt 

##########STEP2############
ls *1.fq.gz|cut -f1 -d"_"|while read line; do echo '#!/bin/bash
#SBATCH --job-name='${line}'BWA
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
module load bwa/0.7.17--h7132678_9
srun --export=all -n 1 -c 128 bwa mem -t 128 -R "@RG\tID:'${line}'\tPL:illumina\tLB:library\tSM:'${line}'" /scratch/pawsey0399/bguo1/Murdoch/11.Oat/Pinyan/Ref/PY6.fa '${line}'_clean_1.fq.gz '${line}'_clean_2.fq.gz >/scratch/pawsey0399/bguo1/Murdoch/11.Oat/Pinyan/01.GWAS/01.SAM/'${line}'.PY6.sam 
'>$line.bwa.sh; done

#################STEP3###########################
ls *sam|cut -f1 -d"."|while read line; do echo '#!/bin/bash
#SBATCH --job-name='${line}'_sam
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
module load samtools/1.15--h3843a85_0
mkdir /scratch/pawsey0399/bguo1/TMP/'${line}'
srun --export=all -n 1 -c 48 samtools view -@ 16 -b '${line}'.PY6.sam |samtools sort -@ 32 -m 5G -T /scratch/pawsey0399/bguo1/TMP/'${line}' -o /scratch/pawsey0399/bguo1/Murdoch/11.Oat/Pinyan/01.GWAS/02.BAM/'${line}'.PY6.sort.bam
' >$line.sort.sh; done
samtools index -c $line
#################STEP4###########################
#!/bin/bash

# 染色体表
cat <<EOF > chr_map.txt
1A GWHCBGG00000044
2A GWHCBGG00000010
3A GWHCBGG00000054
4A GWHCBGG00000007
5A GWHCBGG00000075
6A GWHCBGG00000048
7A GWHCBGG00000030
1C GWHCBGG00000087
2C GWHCBGG00000033
3C GWHCBGG00000072
4C GWHCBGG00000029
5C GWHCBGG00000050
6C GWHCBGG00000086
7C GWHCBGG00000011
1D GWHCBGG00000078
2D GWHCBGG00000035
3D GWHCBGG00000067
4D GWHCBGG00000039
5D GWHCBGG00000043
6D GWHCBGG00000080
7D GWHCBGG00000004
EOF

# 样本循环
for bam in *.bam
do
    sample=$(basename "$bam" .PY6.sort.bam)

    for group in A C D
    do
        script=${sample}.${group}.sh

        # 写头
        cat <<EOF > "$script"
#!/bin/bash
#SBATCH --job-name=${sample}_${group}_Deep
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399

module load singularity/4.1.0-nompi

EOF

        while read chr region
        do
            if [[ $chr == *${group} ]]; then
                cat <<EOF >> "$script"
echo "Running ${sample} ${chr} ..."

srun --export=all -n 1 -c 64 singularity exec /scratch/pawsey0399/bguo1/Singularity_image/deepvariant.sif run_deepvariant \\
--model_type WGS \\
--ref /scratch/pawsey0399/bguo1/Murdoch/11.Oat/Pinyan/Ref/PY6.fa \\
--regions ${region} \\
--reads /scratch/pawsey0399/bguo1/Murdoch/11.Oat/Pinyan/01.GWAS/02.BAM/${sample}.PY6.sort.bam \\
--sample_name ${sample} \\
--output_gvcf /scratch/pawsey0399/bguo1/Murdoch/11.Oat/Pinyan/01.GWAS/03.Deep/${sample}.${chr}.deep.g.vcf.gz \\
--output_vcf /scratch/pawsey0399/bguo1/Murdoch/11.Oat/Pinyan/01.GWAS/03.Deep/${sample}.${chr}.deep.vcf.gz \\
--num_shards 64

EOF
            fi
        done < chr_map.txt   

        chmod +x "$script"
    done
done

############################STEP4-2################################
# 先准备 chr_map
declare -A chr_map
while read chr region
do
    chr_map[$chr]=$region
done < chr_map.txt


# 遍历所有 missing 文件
for file in *.deep.missing
do
    chr=$(basename $file .deep.missing)   # 例如 3C

    while read sample
    do
        script=${sample}.${chr}.rerun.sh

        cat <<EOF > "$script"
#!/bin/bash
#SBATCH --job-name=${sample}_${chr}_Deep
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399

module load singularity/4.1.0-nompi

echo "Re-running ${sample} ${chr} ..."

srun --export=all -n 1 -c 64 singularity exec /scratch/pawsey0399/bguo1/Singularity_image/deepvariant.sif run_deepvariant \\
--model_type WGS \\
--ref /scratch/pawsey0399/bguo1/Murdoch/11.Oat/Pinyan/Ref/PY6.fa \\
--regions ${chr_map[$chr]} \\
--reads /scratch/pawsey0399/bguo1/Murdoch/11.Oat/Pinyan/01.GWAS/02.BAM/${sample}.PY6.sort.bam \\
--sample_name ${sample} \\
--output_gvcf /scratch/pawsey0399/bguo1/Murdoch/11.Oat/Pinyan/01.GWAS/03.Deep/${sample}.${chr}.deep.g.vcf.gz \\
--output_vcf /scratch/pawsey0399/bguo1/Murdoch/11.Oat/Pinyan/01.GWAS/03.Deep/${sample}.${chr}.deep.vcf.gz \\
--num_shards 64

EOF

        chmod +x "$script"

    done < $file

done

#######################STEP5##################
for i in {1..7}{A,C,D}; do echo '#!/bin/bash        
#SBATCH --job-name='${i}'_gln
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=96:00:00
#SBATCH --account=pawsey0399
#SBATCH --export=NONE
module load bcftools/1.15--haf5b3da_0
module load singularity/4.1.0-nompi                                                                                                                                         

srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/Singularity_image/Glnexus.sif glnexus_cli --config DeepVariant --bed '${i}' --threads 128 /scratch/pawsey0399/bguo1/Murdoch/11.Oat/Pinyan/01.GWAS/03.Deep/*'${i}'*.g.vcf.gz |bcftools view --threads 128 - |bgzip -@ 128 -c > Oat_hull_PY.'${i}'.cohort.vcf.gz
' >$i.glnexus.sh; done

##################STEP6######################
for i in {1..7}{A,C,D}; do bcftools annotate --threads 6 --rename-chrs chr.rename -Oz -o ${i}.rename.vcf.gz Oat_hull_PY.${i}.maf25.U80.recode.vcf.gz & done
############change chromosome name first###################
for i in {1..7}{A,C,D}
> do
> plink --vcf Oat_hull_PY.${i}.maf25.U80.recode.vcf.gz --make-bed --out ${i} --allow-extra-chr --vcf-half-call m --threads 128
> done


###############STEP7######################
########kinship file is from emma############
for i in {1..7}{A,C,D}; do gemma -bfile ${i} -k ${i}.kinship.txt -lmm 4 -p pheno2.txt -o gwas_${i} & done

################STEP8#########################
# ===============================
# 1. 加载包
# ===============================
library(ggplot2)
library(dplyr)
library(patchwork)

# ===============================
# 2. 参数
# ===============================
chr_order <- c(paste0(1:7,"A"),
               paste0(1:7,"C"),
               paste0(1:7,"D"))

threshold <- 0.05 / 12708035
log_threshold <- -log10(threshold)

# ===============================
# 3. 读取数据
# ===============================
gwas_list <- list()

for(chr in chr_order){

  file <- paste0("gwas_", chr, ".assoc.txt")

  df <- read.table(file, header=TRUE)

  df$chr_label <- chr
  df$pos_mb <- df$ps / 1e6
  df$logp <- -log10(df$p_score)

  # 标记显著
  df$significant <- df$p_score < threshold

  gwas_list[[chr]] <- df
}

# ===============================
# 4. QQ plot（单染色体）
# ===============================
makeQQ <- function(df){

  obs <- -log10(sort(df$p_score))
  exp <- -log10(ppoints(length(obs)))

  qq <- data.frame(exp=exp, obs=obs)

  ggplot(qq, aes(exp, obs)) +
    geom_point(size=0.8) +
    geom_abline(slope=1, intercept=0, color="red") +
    theme_classic() +
    labs(title=unique(df$chr_label),
         x="Expected -log10(P)",
         y="Observed -log10(P)")
}

# ===============================
# 5. Manhattan（单染色体）
# ===============================
makeMan <- function(df){

  ggplot(df, aes(x=pos_mb, y=logp)) +
    geom_point(data=subset(df, significant==FALSE),
               color="grey70", size=0.5) +
    geom_point(data=subset(df, significant==TRUE),
               color="red", size=0.7) +
    geom_hline(yintercept = log_threshold,
               color="red", linetype="dashed") +
    theme_classic() +
    labs(title=unique(df$chr_label),
         x="Position (Mb)",
         y="-log10(P)")
}

# ===============================
# 6. 输出单染色体图（21个）
# ===============================
for(chr in chr_order){

  df <- gwas_list[[chr]]

  # QQ
  p1 <- makeQQ(df)
  ggsave(paste0("QQ_", chr, ".pdf"),
         p1, width=4, height=4)

  # Manhattan
  p2 <- makeMan(df)
  ggsave(paste0("Manhattan_", chr, ".pdf"),
         p2, width=8, height=6)
}

# ===============================
# 7. 无标题 Manhattan（用于拼图）
# ===============================
makeMan_noTitle <- function(df){

  ggplot(df, aes(x=pos_mb, y=logp)) +
    geom_point(data=subset(df, significant==FALSE),
               color="grey70", size=0.5) +
    geom_point(data=subset(df, significant==TRUE),
               color="red", size=0.7) +
    geom_hline(yintercept = log_threshold,
               color="red", linetype="dashed") +
    theme_classic() +
    theme(
      plot.title = element_blank(),
      axis.title = element_blank()
    )
}

# ===============================
# 8. 3×7 Manhattan 拼图
# ===============================
man_list <- list()

for(chr in chr_order){
  man_list[[chr]] <- makeMan_noTitle(gwas_list[[chr]])
}

manhattan_3x7 <- wrap_plots(man_list, ncol=7)

ggsave("Manhattan_3x7.pdf",
       manhattan_3x7,
       width=24,
       height=21)

# ===============================
# 9. （可选）QQ 拼图
# ===============================
makeQQ_noTitle <- function(df){

  obs <- -log10(sort(df$p_score))
  exp <- -log10(ppoints(length(obs)))

  qq <- data.frame(exp=exp, obs=obs)

  ggplot(qq, aes(exp, obs)) +
    geom_point(size=0.6) +
    geom_abline(slope=1, intercept=0) +
    theme_classic() +
    theme(
      plot.title = element_blank(),
      axis.title = element_blank()
    )
}

qq_list <- lapply(chr_order, function(chr){
  makeQQ_noTitle(gwas_list[[chr]])
})

qq_3x7 <- wrap_plots(qq_list, ncol=7)

ggsave("QQ_3x7.pdf",
       qq_3x7,
       width=20,
       height=10)

###############################################################
###############################################################
###############################################################
###############################################################
###############################################################


###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
#Split plot
library(ggplot2)
library(dplyr)

# ===============================
# 2. 参数
# ===============================
threshold <- 0.05 / 12708035
log_threshold <- -log10(threshold)

# 染色体顺序
chr_order <- c(paste0(1:7,"A"),
               paste0(1:7,"C"),
               paste0(1:7,"D"))

# ===============================
# 3. 读取数据
# ===============================
gwas <- read.table("gwas_all.assoc.txt", header=TRUE)

# ===============================
# 4. 染色体转换（关键）
# ===============================
gwas$chr_label <- NA

gwas$chr_label[gwas$chr %in% 1:7] <- paste0(gwas$chr[gwas$chr %in% 1:7], "A")
gwas$chr_label[gwas$chr %in% 8:14] <- paste0(gwas$chr[gwas$chr %in% 8:14]-7, "C")
gwas$chr_label[gwas$chr %in% 15:21] <- paste0(gwas$chr[gwas$chr %in% 15:21]-14, "D")

# ===============================
# 5. 清理数据（非常重要！！！）
# ===============================
gwas <- gwas %>%
  filter(!is.na(p_wald)) %>%
  filter(p_wald > 0) %>%      # 防止 log10(0)
  filter(!is.na(ps))

# 转换
gwas$pos_mb <- gwas$ps / 1e6
gwas$logp <- -log10(gwas$p_wald)
gwas$significant <- gwas$p_wald < threshold

# ===============================
# 6. 分染色体作图
# ===============================
gwas$chr_label <- as.character(gwas$chr_label)

for(chr in chr_order){

  df <- gwas[gwas$chr_label == chr, ]

  cat("Processing:", chr, "N =", nrow(df), "\n")

  if(nrow(df) == 0){
    cat("Skip:", chr, "\n")
    next
  }

  # QQ
  obs <- -log10(sort(df$p_wald))
  exp <- -log10(ppoints(length(obs)))

  qq <- data.frame(exp=exp, obs=obs)

  p_qq <- ggplot(qq, aes(exp, obs)) +
    geom_point(size=0.8) +
    geom_abline(slope=1, intercept=0, color="red") +
    theme_classic() +
    labs(title=chr)

  ggsave(paste0("QQ_", chr, ".pdf"), p_qq, width=4, height=4)

  # Manhattan
  p_man <- ggplot(df, aes(x=pos_mb, y=logp)) +
    geom_point(color="grey70", size=0.5) +
    geom_point(data=df[df$significant, ],
               color="red", size=0.7) +
    geom_hline(yintercept = log_threshold,
               color="red", linetype="dashed") +
    theme_classic() +
    labs(title=chr)

  ggsave(paste0("Manhattan_", chr, ".pdf"), p_man, width=8, height=6)
}
######################################################
#####combined##########

library(patchwork)
makeQQ_noTitle <- function(df){

  # 去掉异常值（防止空图）
  df <- df[!is.na(df$p_wald) & df$p_wald > 0, ]

  obs <- -log10(sort(df$p_wald))
  exp <- -log10(ppoints(length(obs)))

  qq <- data.frame(exp=exp, obs=obs)

  ggplot(qq, aes(exp, obs)) +
    geom_point(size=0.5) +
    geom_abline(slope=1, intercept=0) +
    theme_classic() +
    theme(
      plot.title = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size=6)
    )
}

# ===============================
# 4. 生成21个QQ图
# ===============================
qq_list <- list()

for(chr in chr_order){

  df <- gwas[gwas$chr_label == chr, ]

  cat("QQ:", chr, "N =", nrow(df), "\n")

  if(nrow(df) == 0){
    next
  }

  qq_list[[chr]] <- makeQQ_noTitle(df)
}

# ===============================
# 5. 拼图（3×7）
# ===============================
qq_3x7 <- wrap_plots(qq_list, ncol=7)

# 保存
ggsave("QQ_3x7.pdf",
       qq_3x7,
       width=20,
       height=10)
#########################################


############genotype heatmap-group#############
############################################################
## Genotype heatmap
## - row = samples
## - column = SNPs
## - group order: hulled -> wild -> hulless
## - only cluster samples
## - NO SNP clustering
## - NO sample names
## - NO group label text
## - left dendrogram retained
## - red/blue genotype color
############################################################

library(vcfR)
library(ComplexHeatmap)
library(circlize)
library(grid)

############################################################
# 1. read VCF
############################################################

vcf <- read.vcfR("gene_only.vcf.gz")

gt <- extract.gt(vcf, element = "GT")

############################################################
# 2. GT -> numeric
############################################################

gt[gt %in% c("0/0", "0|0")] <- 0
gt[gt %in% c("0/1", "1/0", "0|1", "1|0")] <- 1
gt[gt %in% c("1/1", "1|1")] <- 2
gt[gt %in% c("./.", ".|.")] <- NA

############################################################
# 3. build matrix
# row = sample
# col = SNP
############################################################

mat <- t(apply(gt, 2, as.numeric))

############################################################
# 4. read group file
############################################################

group <- read.table(
  "group.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

rownames(group) <- group$sample

############################################################
# 5. match sample order
############################################################

group <- group[rownames(mat), , drop = FALSE]

############################################################
# 6. set group order
############################################################

group$type <- factor(
  group$type,
  levels = c("hulled", "wild", "hulless")
)

############################################################
# 7. reorder samples
############################################################

ord <- order(group$type)

mat <- mat[ord, ]
group <- group[ord, , drop = FALSE]

############################################################
# 8. genotype color
# ref = blue
# het = white
# alt = red
############################################################

col_fun <- colorRamp2(
  c(0, 1, 2),
  c("#2166AC", "white", "#B2182B")
)

############################################################
# 9. left group annotation
############################################################

ha <- rowAnnotation(
  type = group$type,

  col = list(
    type = c(
      hulled  = "#1B9E77",
      wild    = "#E6AB02",
      hulless = "#D95F02"
    )
  ),

  show_annotation_name = FALSE,

  width = unit(8, "mm")
)

############################################################
# 10. draw heatmap
############################################################

pdf(
  "Genotype_heatmap.pdf",
  width = 12,
  height = 10
)

Heatmap(
  mat,

  ##########################################################
  # heatmap color
  ##########################################################

  col = col_fun,

  ##########################################################
  # left annotation
  ##########################################################

  left_annotation = ha,

  ##########################################################
  # sample clustering only
  ##########################################################

  cluster_rows = TRUE,

  ##########################################################
  # NO SNP clustering
  ##########################################################

  cluster_columns = FALSE,

  ##########################################################
  # keep group blocks
  ##########################################################

  row_split = group$type,

  ##########################################################
  # remove names
  ##########################################################

  show_row_names = FALSE,
  show_column_names = FALSE,

  ##########################################################
  # dendrogram on left
  ##########################################################

  row_dend_side = "left",

  ##########################################################
  # legend
  ##########################################################

  heatmap_legend_param = list(
    title = "Genotype",
    at = c(0, 1, 2),
    labels = c("Ref", "Het", "Alt")
  ),

  ##########################################################
  # cleaner border
  ##########################################################

  border = FALSE,

  ##########################################################
  # remove column title
  ##########################################################

  column_title = NULL,
  row_title = NULL
)

dev.off()

###########################################################
###################genotype heatmap-ungroup-clustured######
############################################################
## Genotype heatmap
## - NORMAL clustering
## - no forced group order
## - cluster samples by genotype
## - NO SNP clustering
## - NO sample names
## - NO group label text
## - left dendrogram retained
## - red/blue genotype color
############################################################

library(vcfR)
library(ComplexHeatmap)
library(circlize)
library(grid)

############################################################
# 1. read VCF
############################################################

vcf <- read.vcfR("gene_only.vcf.gz")

gt <- extract.gt(vcf, element = "GT")

############################################################
# 2. GT -> numeric
############################################################

gt[gt %in% c("0/0", "0|0")] <- 0
gt[gt %in% c("0/1", "1/0", "0|1", "1|0")] <- 1
gt[gt %in% c("1/1", "1|1")] <- 2
gt[gt %in% c("./.", ".|.")] <- NA

############################################################
# 3. build matrix
# row = sample
# col = SNP
############################################################

mat <- t(apply(gt, 2, as.numeric))

############################################################
# 4. read group file
############################################################

group <- read.table(
  "group.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

rownames(group) <- group$sample

############################################################
# 5. match sample order
############################################################

group <- group[rownames(mat), , drop = FALSE]

############################################################
# 6. left annotation
############################################################

ha <- rowAnnotation(
  type = group$type,

  col = list(
    type = c(
      hulled  = "#1B9E77",
      wild    = "#E6AB02",
      hulless = "#D95F02"
    )
  ),

  show_annotation_name = FALSE,

  width = unit(8, "mm")
)

############################################################
# 7. genotype color
# ref = blue
# het = white
# alt = red
############################################################

col_fun <- colorRamp2(
  c(0, 1, 2),
  c("#2166AC", "white", "#B2182B")
)

############################################################
# 8. draw heatmap
############################################################

pdf(
  "Genotype_heatmap_clustered.pdf",
  width = 12,
  height = 10
)

Heatmap(
  mat,

  ##########################################################
  # heatmap color
  ##########################################################

  col = col_fun,

  ##########################################################
  # left annotation
  ##########################################################

  left_annotation = ha,

  ##########################################################
  # NORMAL clustering
  ##########################################################

  cluster_rows = TRUE,

  ##########################################################
  # NO SNP clustering
  ##########################################################

  cluster_columns = FALSE,

  ##########################################################
  # remove names
  ##########################################################

  show_row_names = FALSE,
  show_column_names = FALSE,

  ##########################################################
  # dendrogram on left
  ##########################################################

  row_dend_side = "left",

  ##########################################################
  # legend
  ##########################################################

  heatmap_legend_param = list(
    title = "Genotype",
    at = c(0, 1, 2),
    labels = c("Ref", "Het", "Alt")
  ),

  ##########################################################
  # cleaner border
  ##########################################################

  border = FALSE,

  ##########################################################
  # remove titles
  ##########################################################

  column_title = NULL,
  row_title = NULL
)

dev.off()






###############Kmer-gwas##########################
###############STEP1#############################
while read id; do echo '#!/bin/bash
#SBATCH --job-name='${id}'
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
module load singularity/4.1.0-nompi
mkdir 02.KGWAS/01_kmc/tmp/'${id}'
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/Singularity_image/kmcv3.2.4.sif kmc -k31 -ci3 -cx100000 -t128 -m200 @Raw_data/'${id}'.list 02.KGWAS/01_kmc/canon/'${id}' 02.KGWAS/01_kmc/tmp/'${id}'  
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/Singularity_image/kmcv3.2.4.sif kmc -k31 -ci0 -b -t128 -m200 @Raw_data/'${id}'.list 02.KGWAS/01_kmc/noncanon/'${id}' 02.KGWAS/01_kmc/tmp/'${id}' '>$id.kmc.sh
done <mapping.txt

#############STEP2##############################
cat mapping.txt |while read id; do echo '#!/bin/bash
#SBATCH --job-name='${id}'_strand
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
module load singularity/4.1.0-nompi
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate kgwas
srun --export=all -n 1 -c 4  /scratch/pawsey0399/bguo1/software/kgwas/bin/kmers_add_strand_information -c 01_kmc/canon/'${id}' -n 01_kmc/noncanon/'${id}' -k 31 -o 02_kmers_strand/'${id}' '>$id.strand.sh; done

###############STEP3###########################
ls $STRAND/* > $LIST/kmers_paths_raw.txt
awk -F'/' '{print $0"\t"$NF}' $LIST/kmers_paths_raw.txt > $LIST/kmers_list_paths.txt


#!/bin/bash
#SBATCH --job-name=kgws-03
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=800G
#SBATCH --cpus-per-task=4
#SBATCH --time=96:00:00
#SBATCH --account=pawsey0399
module load singularity/4.1.0-nompi
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate kgwas
srun --export=all -n 1 -c 4  /scratch/pawsey0399/bguo1/software/kgwas/bin/list_kmers_found_in_multiple_samples -l /scratch/pawsey0399/bguo1/Murdoch/11.Oat/Pinyan/02.KGWAS/03_kmers_list/kmers_list_paths.txt -k 31 --mac 50 -p 0.1 -o /scratch/pawsey0399/bguo1/Murdoch/11.Oat/Pinyan/02.KGWAS/04_kmers_filter/kmers_to_use

################STEP4##################################

#!/bin/bash
#SBATCH --job-name=kgws-04
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=800G
#SBATCH --cpus-per-task=4
#SBATCH --time=96:00:00
#SBATCH --account=pawsey0399
module load singularity/4.1.0-nompi
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate kgwas
srun --export=all -n 1 -c 4  /scratch/pawsey0399/bguo1/software/kgwas/bin/build_kmers_table -l /scratch/pawsey0399/bguo1/Murdoch/11.Oat/Pinyan/04.target-kgwas/03_kmers_list/kmers_list_paths-CRR.txt -k 31 -a /scratch/pawsey0399/bguo1/Murdoch/11.Oat/Pinyan/04.target-kgwas/04_kmers_filter/kmers_to_use-CRR -o /scratch/pawsey0399/bguo1/Murdoch/11.Oat/Pinyan/04.target-kgwas/05_kmers_table/kmers_table-CRR


############STEP6##############################################
/scratch/pawsey0399/bguo1/Murdoch/11.Oat/Pinyan/04.target-kgwas> /scratch/pawsey0399/bguo1/software/kgwas/bin/emma_kinship_kmers -t 05_kmers_table/kmers_table-SRR -k 31 --maf 0.2 >06_kinship/kmers_table-SRR.kinship

######################STEP7####################################
python2.7 /scratch/pawsey0399/bguo1/software/kgwas2/kmers_gwas.py --pheno pheno_ssr.txt --kmers_table 05_kmers_table/kmers_table-SRR -p 128 -l 31 --maf 0.05 --outdir 07_kgwas/kmers-SRR-gwas_out --dont_remove_intermediates


####################STEP8#####################################
awk '{print$2}' phenotype_value.assoc.txt |sort|uniq |awk '{print ">"NR"\n"$1}' >sig_kmers.fa

bowtie2 -x PY6.4D -f ../07_kgwas/kmers-CRR-gwas_out-2/kmers/output/sig_kmers.fa -k 10 --very-sensitive -p 128 -S sig_kmers.sam

samtools sort -o -o sig_kmers.sort.bam sig_kmers.sam
samtools index -c sig_kmers-all.sort.bam

bedtools makewindows -g PY6.fa.fai -w 5000 > windows.5k.bed
bedtools coverage -a windows.5k.bed -b sig_kmers.sort.bam >kmer_density.5k.txt

 



##################################################################
###################################################################
####################run target region Kgwas########################
###################STEP1##########################################
for f in SRR*.sort.bam; do line=$(basename $f .sort.bam)
cat <<EOF > ${line}.reads.sh
#!/bin/bash
#SBATCH --job-name=${line}
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399

module load samtools/1.15--h3843a85_0
srun --export=all -n 1 -c 8 samtools view -@ 8 -bh -f 4 ${line}.sort.bam > ../00.kmer-bam/${line}.unmapped.bam
srun --export=all -n 1 -c 8 samtools view -h ${line}.sort.bam | awk 'BEGIN{OFS="\t"} /^@/ || $6 ~ /S|H/' | samtools view -@ 8 -b - > ../00.kmer-bam/${line}.clipped.tmp.bam

srun --export=all -n 1 -c 8 samtools view -@ 8 -bh -L region.bed -U ../00.kmer-bam/${line}.clipped.bam ../00.kmer-bam/${line}.clipped.tmp.bam
rm ../00.kmer-bam/${line}.clipped.tmp.bam

EOF
 done

ls CRR773*.sort.bam|cut -f1 -d"."|while read line; do samtools view -bh $line.PY6.sort.bam GWHCBGG00000039:400000000-456028319 >../00.kmer-bam/$line.candi-region.bam & done

#############STEP2#########################################
######all the bam file should be sorted by name first, which means samtools sort -n####################

ls CRR*.clipped.bam|cut -f1 -d"."|while read line; do
echo '#!/bin/bash        
#SBATCH --job-name='${line}'-combine
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
module load samtools/1.15--h3843a85_0
srun --export=all -n 1 -c 8 samtools merge -@ 8 -o ../00-1.kmer-combinebam/'${line}'.merge.bam '${line}'.candi-region.bam -h '${line}'.PY6.clipped.bam '${line}'.PY6.unmapped.bam' >$line.merge.sh
done

#############STEP3###################################

ls |cut -f1 -d"."|while read line; do echo '#!/bin/bash        
#SBATCH --job-name='${line}'-fastq
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
module load samtools/1.15--h3843a85_0
srun --export=all -n 1 -c 8 samtools fastq -@ 8 -1 ../00-3.kmer-fastq/'${line}'_1.fq.gz -2 ../00-3.kmer-fastq/'${line}'_2.fq.gz -s ../00-3.kmer-fastq/'${line}'_s.fq.gz '${line}'.merge.bam' >$line.fastq.sh; done


############STEP4####################################

contindued kgwas STEP1



##################kgwas plot#############################
library(data.table)
library(ggplot2)
py_map <- data.frame(
  chr = c("1A","2A","3A","4A","5A","6A","7A",
          "1C","2C","3C","4C","5C","6C","7C",
          "1D","2D","3D","4D","5D","6D","7D"),
  ref = c("GWHCBGG00000044","GWHCBGG00000010","GWHCBGG00000054","GWHCBGG00000007",
          "GWHCBGG00000075","GWHCBGG00000048","GWHCBGG00000030",
          "GWHCBGG00000087","GWHCBGG00000033","GWHCBGG00000072","GWHCBGG00000029",
          "GWHCBGG00000050","GWHCBGG00000086","GWHCBGG00000011",
          "GWHCBGG00000078","GWHCBGG00000035","GWHCBGG00000067","GWHCBGG00000039",
          "GWHCBGG00000043","GWHCBGG00000080","GWHCBGG00000004"),
  stringsAsFactors = FALSE
)
convert_chr <- function(df, ref_name) {
  
  # PY: scaffold → chr
  if (ref_name == "PY") {
    df$RNAME <- sapply(df$RNAME, function(x) {
      hit <- py_map$chr[py_map$ref == x]
      if (length(hit) == 0) return(NA)
      return(hit)
    })
  }
  
  # SFS: chr4D → 4D
  if (ref_name == "SFS") {
    df$RNAME <- gsub("^chr", "", df$RNAME)
  }
  
  return(df)
}
read_sam_safe <- function(f, ref_name) {
  
  sam <- fread(f, header = FALSE, fill = TRUE)
  
  # ⭐只保留SAM标准11列（去掉optional tags）
  sam <- sam[, 1:11]
  
  setnames(sam, c(
    "QNAME","FLAG","RNAME","POS","MAPQ","CIGAR",
    "RNEXT","PNEXT","TLEN","SEQ","QUAL"
  ))
  
  sam$POS <- as.numeric(sam$POS)
  sam$REF <- ref_name
  
  sam <- convert_chr(sam, ref_name)
  
  return(sam)
}
files <- list.files(pattern = "sig_kmers-all-.*\\.sam$")

all_data <- list()

for (f in files) {
  
  ref <- sub("sig_kmers-all-(.*)\\.sam", "\\1", f)
  cat("Reading:", ref, "\n")
  
  sam <- read_sam_safe(f, ref)
  
  all_data[[ref]] <- sam
}

df <- rbindlist(all_data)

# 去掉无法映射染色体
df <- df[!is.na(df$RNAME)]
chr_keep <- c(
  "1A","2A","3A","4A","5A","6A","7A",
  "1C","2C","3C","4C","5C","6C","7C",
  "1D","2D","3D","4D","5D","6D","7D"
)

df <- df[df$RNAME %in% chr_keep]
df$group <- ifelse(
  df$MAPQ == 0, "multi",
  ifelse(df$MAPQ <= 30, "ambiguous",
    ifelse(df$MAPQ == 255, "MAPQ255", "uniq")
  )
)
bin_size <- 100000

df$RNAME <- factor(df$RNAME, levels = chr_order)
df$MB <- df$POS / 1e6
df$BIN <- floor(df$POS / bin_size) * bin_size
plot_df <- df[, .N, by = .(REF, RNAME, BIN, group)]
plot_df$MB_BIN <- plot_df$BIN / 1e6
chr_order <- c(
  paste0(1:7, "A"),
  paste0(1:7, "C"),
  paste0(1:7, "D")
)

plot_df$RNAME <- factor(plot_df$RNAME, levels = chr_order)

p1 <- ggplot(plot_df, aes(x = MB_BIN, y = N, color = group)) +
    
    geom_line(alpha = 0.8, linewidth = 0.4) +
    
    facet_grid(REF ~ RNAME, scales = "free_x") +
    
    theme_bw(base_size = 12) +
    
    labs(
        title = "Genome-wide k-mer distribution across 4 references",
        x = "Genomic position (Mb)",
        y = "k-mer count"
    ) +
    
    scale_color_manual(values = c(
        "multi" = "black",
        "ambiguous" = "grey60",
        "uniq" = "blue",
        "MAPQ255" = "red"
    )) +
    
    theme(
        strip.text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

ggsave("kgwas-4ref_21chr_Mb-order.pdf", p1, width = 18, height = 10, dpi = 300)
###############only 4D zoom in 4 ref####################################
get_4d_terminal_real <- function(df_ref) {
    
    d <- df_ref[df_ref$RNAME == "4D"]
    
    end_pos <- max(d$POS, na.rm = TRUE)
    start_pos <- end_pos - 1e6
    
    sub <- d[d$POS >= start_pos & d$POS <= end_pos]
    
    return(sub)
}
df_list <- split(df, df$REF)

df_4D_win <- rbindlist(lapply(df_list, get_4d_terminal_real))
bin_size <- 50000

df_4D_win$BIN <- floor(df_4D_win$POS / bin_size) * bin_size
df_4D_win$MB_BIN <- df_4D_win$BIN / 1e6
plot_zoom <- df_4D_win[, .N, by = .(REF, MB_BIN, group)]
pB <- ggplot(plot_zoom, aes(x = MB_BIN, y = N, fill = group)) +
    
    geom_bar(stat = "identity", position = "stack") +
    
    facet_wrap(~REF, ncol = 2, scales = "free_x") +
    
    theme_bw(base_size = 12) +
    
    labs(
        title = "Chr4D candidate region ",
        x = "Genomic position (Mb)",
        y = "k-mer count"
    ) +
    
    scale_fill_manual(values = c(
        "multi" = "black",
        "ambiguous" = "grey60",
        "uniq" = "blue",
        "MAPQ255" = "red"
    )) +
    
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10)
    )

ggsave("4D_terminal_realCoord_4ref.pdf", pB, width = 12, height = 8, dpi = 300)












##############zoom in candidate region####################
df_py <- df[df$REF == "PY"]
target_chr <- "4D"
start <- 455700000
end   <- 456000000
sub <- df_py[
  df_py$RNAME == target_chr &
  df_py$POS >= start &
  df_py$POS <= end
]
bin_size <- 100000
sub$BIN <- floor(sub$POS / bin_size) * bin_size
sub$MB_BIN <- sub$BIN / 1e6
plot_df2 <- sub[, .N, by = .(MB_BIN, group)]
p2 <- ggplot(plot_df2, aes(x = MB_BIN, y = N, color = group)) +
    
    geom_line(linewidth = 0.8) +
    
    theme_bw(base_size = 12) +
    
    labs(
        title = "PY chr4D: 450Mb to end",
        x = "Genomic position (Mb)",
        y = "k-mer count"
    ) +
    
    scale_color_manual(values = c(
        "multi" = "black",
        "ambiguous" = "grey60",
        "uniq" = "blue",
        "MAPQ255" = "red"
    ))

ggsave("PY_4D_450-460_zoom.pdf", p2, width = 10, height = 5, dpi = 300)











###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
####################DEPTH##################################
##################STEP1#####################################
ls *.bam|cut -f1 -d"."|while read line; do samtools view -bh -o ../../../../03.Depth/$line.PY6.4D.sort.bam $line.PY6.sort.bam GWHCBGG00000039 & done



################STEP2#####################################

bedtools makewindows -g ../Ref/PY6.fa.fai -w 100000|grep GWHCBGG00000039 > 100Kb-windows/4D.100k.bed
cat name.txt |while read line; do singularity exec /scratch/pawsey0399/bguo1/Singularity_image/mosdepth.sif mosdepth -t 4 -c GWHCBGG00000039 -Q 30 -b 4D.100k.bed ${line} ../${line}.PY6.4D.sort.bam& done

################STEP3#####################################

library(ggplot2)
library(dplyr)
library(stringr)

files <- list.files(pattern = "*.bed.gz")
all_data <- lapply(files, function(f) {
    df <- read.table(f, header = FALSE, sep = "\t")
    colnames(df) <- c("chr", "start", "end", "depth")
    
    # 样本名取文件名第一个点之前
    sample_name <- str_split(basename(f), "\\.")[[1]][1]
    df$sample <- sample_name
    
    # 去除异常值 (>5*mean)
    m <- mean(df$depth)
    df <- df[df$depth <= 5*m, ]
    
    # 归一化
    df$depth_norm <- df$depth/m
    
    # 区间中点作为位置
    df$pos <- (df$start + df$end)/2
    df$pos_Mb <- df$pos / 1e6
    
    return(df)
})
plot_data <- do.call(rbind, all_data)
trait_info <- read.table("trait.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
trait_info <- trait_info %>% filter(Accession %in% plot_data$sample)
trait_info$Type<-factor(trait_info$Type,levels = c("Hulless","Hulled","Wild"))
plot_data <- plot_data %>%
    left_join(trait_info, by = c("sample" = "Accession")) %>%
    arrange(desc(Type), sample)
plot_data$sample <- factor(plot_data$sample, levels = unique(plot_data$sample))
group_boundaries <- plot_data %>%
    group_by(Type) %>%
    summarise(y = max(as.numeric(factor(sample)))) %>%
    arrange(desc(Type))
group_boundaries <- group_boundaries %>%
  mutate(y = cumsum(y))
boundary_lines <- head(group_boundaries$y, -1)+0.5

p_all <- ggplot(plot_data, aes(x = pos_Mb, y = sample, fill = depth_norm)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5) +
    geom_hline(yintercept = boundary_lines, color = "black", size = 0.8) +
    annotate(
        "text",
        x = min(plot_data$pos_Mb) - 20,
        y = group_boundaries$y - diff(c(0, group_boundaries$y))/2,
        label = group_boundaries$Type,
        angle = 90, vjust = 0.5, size = 5
    ) +
    theme_bw() +
    labs(
        x = "Genomic position (Mb)",
        y = "Accession",
        fill = "Normalized Depth",
        title = "Chr4D Depth Heatmap"
    ) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(10, 20, 10, 70)
    )
ggsave("depth_heatmap_grouped_full.pdf", p_all, width = 12, height = 102, dpi = 300,limitsize = F)


zoom_data <- plot_data %>% filter(pos >= 455000000)
zoom_group_boundaries <- zoom_data %>%
    group_by(Type) %>%
    summarise(y = max(as.numeric(factor(sample)))) %>%
    arrange(desc(Type))
zoom_group_boundaries <- zoom_group_boundaries %>%
    mutate(y = cumsum(y))
zoom_boundary_lines <- head(zoom_group_boundaries$y, -1)+0.5
p_zoom <- ggplot(zoom_data, aes(x = pos_Mb, y = sample, fill = depth_norm)) +
    geom_tile() +
    geom_text(aes(label = round(depth_norm, 2)), size = 2) +  # 显示数字
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5) +
    geom_hline(yintercept = zoom_boundary_lines, color = "black", size = 0.8) +
    annotate(
        "text",
        x = min(zoom_data$pos_Mb) - 0.1,
        y = zoom_group_boundaries$y - diff(c(0, zoom_group_boundaries$y))/2,
        label = zoom_group_boundaries$Type,
        angle = 90, vjust = 0.5, size = 5
    ) +
    theme_bw() +
    labs(
        x = "Genomic position (Mb)",
        y = "Accession",
        fill = "Normalized Depth",
        title = "Chr4D:455Mb+ Depth Heatmap "
    ) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(10, 20, 10, 70)
    )
ggsave("depth_heatmap_grouped_zoom_455Mb.pdf", p_zoom, width = 12, height = 102, dpi = 300,limitsize = F)

write.table(zoom_data, "zoom_455Mb_grouped_data.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

##########################boxplot#########################################








focus_data2 <- plot_data %>%    
    filter(       
        pos_Mb >= 455,        
        pos_Mb <= 456,        
        Type %in% c("Hulled", "Hulless")        
    ) %>%    
    mutate(        
        bin = cut(            
            pos_Mb,            
            breaks = seq(455, 456.0001, by = 0.1),            
            include.lowest = TRUE,            
            right = FALSE            
        )        
    )

stats_df <- focus_data2 %>%    
    group_by(bin) %>%    
    summarise(      
        hulled_median =            
            median(depth_norm[Type == "Hulled"], na.rm = TRUE),        
        hulless_median =           
            median(depth_norm[Type == "Hulless"], na.rm = TRUE),        
        diff =           
            hulless_median - hulled_median,      
        fold =          
            hulless_median / hulled_median,      
        p =           
           wilcox.test(depth_norm ~ Type)$p.value,        
        .groups = "drop"        
    ) %>%    
    mutate(        
        padj = p.adjust(p, method = "BH"),        
        label = paste0(          
            "Δ=",           
            round(diff, 2),           
            "\nFDR=",         
            format(padj, scientific = TRUE, digits = 2)          
        )      
    )
p_box <- ggplot(    
    focus_data2,   
    aes(x = bin, y = depth_norm, fill = Type)    
) +    
    geom_boxplot(      
        outlier.shape = NA,       
        alpha = 0.8,      
        width = 0.6,      
        position = position_dodge(width = 0.7)      
    ) + 
    geom_jitter(       
        aes(color = Type),        
        position = position_jitterdodge(          
            jitter.width = 0.2,           
            dodge.width = 0.7          
        ),      
        size = 0.8,      
        alpha = 0.35     
    ) + 
    scale_fill_manual(    
        values = c(        
            "Hulless" = "#E64B35",        
            "Hulled" = "#4DBBD5"       
        )    
    ) +    
    scale_color_manual(      
           values = c(          
            "Hulless" = "#E64B35",         
            "Hulled" = "#4DBBD5"        
        )     
    ) +
    # 显示统计信息
    geom_text(     
        data = stats_df,    
        inherit.aes = FALSE,   
        aes(       
            x = bin,      
            y = 5,      
            label = label     
        ),  
        size = 3,  
        lineheight = 0.9    
    ) +
    theme_bw() + 
    labs(        
        x = "Genomic position (Mb, 100 kb bins)",        
        y = "Normalized depth",        
        fill = "Type",        
        color = "Type",        
        title = "Depth comparison in Chr4D:455–456 Mb"        
    ) +    
    coord_cartesian(ylim = c(0, 5)) +    
    theme(        
        axis.text.x = element_text(            
            angle = 45,            
            hjust = 1            
        ),        
        panel.grid.minor = element_blank(),        
        plot.title = element_text(            
            hjust = 0.5,            
            face = "bold"            
        )        
    )




















