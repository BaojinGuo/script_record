ls 01.minimap/*.bam|cut -f2 -d "/"|while read line; do filename=$(basename "$line"); prefix=${filename%.Morex.sort.bam}; echo '#!/bin/bash
#SBATCH --job-name=pbsv
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate pbsv
srun --export=all -n 1 -c 128 pbsv discover 01.minimap/'"$line"' 02.result/'"$prefix"'.svsig.gz 
srun --export=all -n 1 -c 128 pbsv call Morex.V3.chr.fasta 02.result/'"$prefix"'.svsig.gz 02.result/'"$prefix"'.pbsv.vcf --min-N-in-gap 1000 --filter-near-reference-gap 0 --ccs' >$prefix.pbsv.sh 
done


###merge and multisample call

 pbsv call Morex.V3.chr.fasta 02.result/Baudin.svsig.gz 02.result/Beecher.svsig.gz 02.result/Betzes.svsig.gz 02.result/Buff.pbsv.svsig.gz 02.result/Buff.svsig.gz 02.result/Buloke.svsig.gz 02.result/Clipper.svsig.gz 02.result/Compass.svsig.gz 02.result/Eth69.svsig.gz 02.result/Fathom.svsig.gz 02.result/Flagship.svsig.gz 02.result/Gairdner.svsig.gz 02.result/Galleon.svsig.gz 02.result/Halcyon.svsig.gz 02.result/Harrington.svsig.gz 02.result/ILAN15.svsig.gz 02.result/Laperouse.svsig.gz 02.result/LaTrobe.svsig.gz 02.result/Leabrook.svsig.gz 02.result/Maximus.svsig.gz 02.result/Mundah.svsig.gz 02.result/Prior.svsig.gz 02.result/RGT_Planet.svsig.gz 02.result/Sahara3771.svsig.gz 02.result/Schooner.svsig.gz 02.result/Stirling.svsig.gz 02.result/Vlamingh.svsig.gz 02.result/Yambla.svsig.gz 02.result/Yeti.svsig.gz 02.result/merge_pbsv.vcf --min-N-in-gap 1000 --filter-near-reference-gap 0 --ccs
