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
or i in {1..7}{A,C,D}; do echo '#!/bin/bash        
#SBATCH --job-name='${i}'_gln
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
#SBATCH --export=NONE
module load bcftools/1.15--haf5b3da_0
module load singularity/4.1.0-nompi                                                                                                                                         cd /scratch/pawsey0399/bguo1/Murdoch/11.Oat/Pinyan/01.GWAS/04.GLnexus/'${i}'
srun --export=all -n 1 -c 64 singularity exec /scratch/pawsey0399/bguo1/Singularity_image/Glnexus.sif glnexus_cli --config DeepVariant --bed '${i}' --threads 64 /scratch/pawsey0399/bguo1/Murdoch/11.Oat/Pinyan/01.GWAS/03.Deep/*'${i}'*.g.vcf.gz |bcftools view --threads 64 - |bgzip -@ 64 -c > Oat_hull_PY.'${i}'.cohort.vcf.gz
' >$i.glnexus.sh; done























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












