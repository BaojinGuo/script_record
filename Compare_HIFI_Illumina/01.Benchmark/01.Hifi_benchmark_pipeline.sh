###this is for SV benchmark
###samples:RGT_Planet(cultivar, UK) Igri(cultivar, Germany) OUN333(landrace, Nepal) Vlaminagh(cultivar, Austrilian) Hockett(cultivar, USA)
##first step, run pipeline of each SV callers (sniffles SVIM cuteSV)
ls *.fastq.gz|while read line; do filename=$(basename "$line"); prefix=${filename%.fastq.gz}; echo '#!/bin/bash
#SBATCH --job-name=Mini
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
module load singularity/4.1.0-nompi
module load samtools/1.15--h3843a85_0
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/Singularity_image/minimap2.sif minimap2 -ax map-hifi -c --MD -R '"'@RG\tID:"$prefix"\tSM:"$prefix"'"' -I100g -t 128 -Y MorexV3.fa '"$line"' |samtools sort -@ 128 -O BAM -o 01.minimap/'"$prefix"'.Morex.sort.bam' >$prefix.mini.sh; done


##cuteSV
ls 01.minimap/*.bam|cut -f2 -d "/"|while read line; do filename=$(basename "$line"); prefix=${filename%.Morex.sort.bam};mkdir -p work_dictionary/$prefix; echo '#!/bin/bash
#SBATCH --job-name cutesv
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate cutesv
srun --export=all -n 1 -c 128 cuteSV 01.minimap/'"$line"' MorexV3.fa 04.cuteSV/'"$prefix"'.cuteSV.vcf work_dictionary/'"$prefix"' -s 5 -S '"$prefix"' -l 50 --genotype --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 -t 128' >$prefix.cuteSV.sh; done



##sniffles
ls 01.minimap/*.bam|cut -f2 -d "/"|while read line; do filename=$(basename "$line"); prefix=${filename%.Morex.sort.bam}; echo '#!/bin/bash
#SBATCH --job-name=sniffles
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=2:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate sniffle
srun --export=all -n 1 -c 128 sniffles --input 01.minimap/'"$line"' --vcf 02.sniffles/'"$prefix"'.Morex.vcf --reference MorexV3.fa --snf 02.sniffles/'"$prefix"'.Morex.snf --threads 128 --minsupport 5 --long-ins-length 100000 --long-del-length 100000' >$prefix.sniffle.sh ; done


##svim
ls 01.minimap/*.bam|cut -f2 -d"/"|while read line; do filename=$(basename "$line"); prefix=${filename%.Morex.sort.bam}; echo '#!/bin/bash
#SBATCH --job-name=SVIM
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate svim
srun --export=all -n 1 -c 128 svim alignment 03.SVIM/'"$prefix"' 01.minimap/'"$line"' MorexV3.fa  --minimum_depth 5 --min_sv_size 50 --sample '"$prefix"' '>$prefix.svim.sh; done

###second step, retain homozygous sites and only concern about SV types with INS DEL INV and DUP
ls *.sh|cut -f1 -d '.'|sort|uniq|while read line; do awk 'BEGIN {OFS="\t"} /^#/ || ($0 ~ /1\/1/) {print}' 02.sniffles/${line}.Morex.vcf >05.SURVIVOR/${line}/${line}.sniffles.hom.vcf; done
ls *.sh|cut -f1 -d '.'|sort|uniq|while read line; do awk 'BEGIN {OFS="\t"} /^#/ || ($0 ~ /1\/1/) {print}' 04.cuteSV/${line}.cuteSV.vcf >05.SURVIVOR/${line}/${line}.cuteSV.hom.vcf; done
ls *.sh|cut -f1 -d '.'|sort|uniq|while read line; do awk 'BEGIN {OFS="\t"} /^#/ || ($0 ~ /1\/1/) {print}' 03.SVIM/${line}/variants.vcf >05.SURVIVOR/${line}/${line}.SVIM.hom.vcf; done

ls|awk 'NR==1 || NR==2 || NR==3 || NR==4 || NR==9 {print$0}'|while read line; do python 00.SVvcf.filter.py $line/$line.cuteSV.hom.vcf $line/$line.cuteSV.final.vcf; done
ls|awk 'NR==1 || NR==2 || NR==3 || NR==4 || NR==9 {print$0}'|while read line; do python 00.SVvcf.filter.py $line/$line.SVIM.hom.vcf $line/$line.SVIM.final.vcf; done
ls|awk 'NR==1 || NR==2 || NR==3 || NR==4 || NR==9 {print$0}'|while read line; do python 00.SVvcf.filter.py $line/$line.sniffles.hom.vcf $line/$line.sniffles.final.vcf; done

###third step, detect union sites with SURVIVOR with parameters '50 2 1 1 0 50' which means the sites will be exist if 2 of 3 callers get 
ls *.vcf >samplefile
ls |while read line; do cd $line; SURVIVOR merge samplefile 50 2 1 1 0 50 ${line}.SURVIVOR.final.vcf; cd ..; done

###forth step, Summary and statistics 
ls *.vcf|cut -f1 -d'.'|while read line; do python 01.extract_SV_info.py ${line}.SURVIVOR.final.vcf 01.Summary/${line}.summary; done
ls *.vcf|cut -f1 -d'.'|while read line; do python 02.statistics.py 01.Summary/${line}.summary 01.Summary/${line}.total 01.Summary/${line}.chr 01.Summary/${line}.len; done













