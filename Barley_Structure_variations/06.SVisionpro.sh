ls 02.minimap/*.bam|cut -f2 -d"/"|cut -f1 -d "."|while read line; do echo '#!/bin/bash
#SBATCH --job-name=SVisionpro
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate svisionpro
srun --export=all -n 1 -c 128 SVision-pro --target_path 02.minimap/'${line}'.Morex.sort.bam --genome_path MorexV3.fa --model_path /scratch/pawsey0399/bguo1/software/SVision-pro/src/pre_process/model_liteunet_256_8_16_32_32_32.path --out_path 03.svisionpro/ --detect_mode germline --sample '$line' --sample_name '$line' --preset hifi --process_num 128 --min_supp 3 --max_sv_size 100000' >$line.svisionpro.sh; done
