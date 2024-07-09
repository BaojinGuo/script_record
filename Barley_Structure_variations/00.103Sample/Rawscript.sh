


###SVision-pro1.8
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


###Sniffles2.4
ls 02.minimap/*.bam|cut -f2 -d"/"|while read line; do filename=$(basename "$line"); prefix=${filename%.Morex.sort.bam}; echo '#!/bin/bash
#SBATCH --job-name=sniffles
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate sniffle
srun --export=all -n 1 -c 128 sniffles --input 02.minimap/'"$line"'--vcf 03.sniffles/'"$prefix"'.Morex.vcf --reference MorexV3.fa --snf 03.sniffles/'"$prefix"'.Morex.snf --threads 128 --minsupport 3 --long-ins-length 100000 --long-del-length 100000' >$prefix.sniffle.sh
done

####SVIM2.0.0
ls 02.minimap/*.bam|cut -f2 -d"/"|while read line; do filename=$(basename "$line"); prefix=${filename%.Morex.sort.bam}; echo '#!/bin/bash
#SBATCH --job-name=SVIM
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate svim
srun --export=all -n 1 -c 128 svim alignment 03.SVIM/'"$prefix"' 02.minimap/'"$line"' MorexV3.fa  --minimum_depth 3 --min_sv_size 50 --max_sv_size 100000 --sample '"$prefix"' '>$prefix.svim.sh; done

####cuteSV2.1.1
ls 02.minimap/*.bam|cut -f2 -d "/"|while read line; do filename=$(basename "$line"); prefix=${filename%.Morex.sort.bam};mkdir -p work_dictionary/$prefix; echo '#!/bin/bash
#SBATCH --job-name cutesv
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate cutesv
srun --export=all -n 1 -c 128 cuteSV 02.minimap/'"$line"' MorexV3.fa 03.cuteSV/'"$prefix"'.cuteSV.vcf work_dictionary/'"$prefix"' -s 3 -S '"$prefix"' -l 50 -L 100000 --genotype --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 -t 128' >$prefix.cuteSV.sh; done 

###PBSV2.9.0
ls 02.minimap/*.bam|cut -f2 -d "/"|while read line; do filename=$(basename "$line"); prefix=${filename%.Morex.sort.bam}; echo '#!/bin/bash
#SBATCH --job-name=pbsv
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate pbsv
srun --export=all -n 1 -c 128 pbsv discover --hifi -s '"$prefix"' -y 97 02.minimap/'"$line"' 03.pbsv/'"$prefix"'.svsig.gz
srun --export=all -n 1 -c 128 pbsv call --hifi -m 50 --max-ins-length 100000 --max-dup-length 100000 -O 3 --min-N-in-gap 1000 --filter-near-reference-gap 0 MorexV3.fa 03.pbsv/'"$prefix"'.svsig.gz 03.pbsv/'"$prefix"'.pbsv.vcf'>$prefix.pbsv.sh 
done
