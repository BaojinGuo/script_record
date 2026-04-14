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


