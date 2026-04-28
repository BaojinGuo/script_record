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





















