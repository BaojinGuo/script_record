ls *.bam|while read line; do srun -n 1 -c 128 samtools index -c -@ 128 $line & 
done ###all bam files should be indexed.

ls 01.minimap/*.bam|cut -f2 -d "/"|while read line; do filename=$(basename "$line"); prefix=${filename%.Morex.sort.bam}; echo '#!/bin/bash
#SBATCH --job-name=sniffles
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate sniffle
srun --export=all -n 1 -c 128 sniffles --input 01.minimap/'"$line"' --vcf 02.sniffles/'"$prefix"'.Morex.vcf --reference Morex.V3.chr.fasta --snf 02.sniffles/'"$prefix"'.Morex.snf --threads 128 --minsupport 1 --long-ins-length 100000000   --long-del-length 100000000' >$prefix.sniffle.sh 
done


###multisamples
srun -c 128 -n 1 sniffles --input list.tsv --vcf 02.sniffles/combineSV_sniffle.vcf --reference Morex.V3.chr.fasta --threads 128 --minsupport 1 --long-ins-length 100000000 --long-del-length 100000000
###list.tsv----all snf files of samples
02.sniffles/Baudin.Morex.snf
02.sniffles/Beecher.Morex.snf
02.sniffles/Betzes.Morex.snf
02.sniffles/Buff.Morex.snf
02.sniffles/Buloke.Morex.snf
02.sniffles/Clipper.Morex.snf
02.sniffles/Compass.Morex.snf
02.sniffles/Eth69.Morex.snf
02.sniffles/Fathom.Morex.snf
02.sniffles/Flagship.Morex.snf
02.sniffles/Gairdner.Morex.snf
02.sniffles/Galleon.Morex.snf
02.sniffles/Halcyon.Morex.snf
02.sniffles/Harrington.Morex.snf
02.sniffles/ILAN15.Morex.snf
02.sniffles/Laperouse.Morex.snf
02.sniffles/LaTrobe.Morex.snf
02.sniffles/Leabrook.Morex.snf
02.sniffles/Maximus.Morex.snf
02.sniffles/Mundah.Morex.snf
02.sniffles/Prior.Morex.snf
02.sniffles/RGT_Planet.Morex.snf
02.sniffles/Sahara3771.Morex.snf
02.sniffles/Schooner.Morex.snf
02.sniffles/Stirling.Morex.snf
02.sniffles/Vlamingh.Morex.snf
02.sniffles/Yambla.Morex.snf
02.sniffles/Yeti.Morex.snf
