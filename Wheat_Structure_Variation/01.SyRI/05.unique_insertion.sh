#get all assembly all chromosome lists
ls */syri/*/*.out|grep -v 3504 >list.syri.txt
#get 3504 all chromosome syri out
ls */syri/*/*.out|grep 3504|while read line
  do
  cat $line >>3504_CS21syri.out
  done
#run script unique_ins.py 
python unique_ins.py -F 3504_CS21syri.out -L list.syri.txt -O 3504_CS21syri.unique_insertion.out
#extract desirable length insertion sequences
awk '$8 - $7 > 10000 && $8 - $7 < 50000 {print $0}' 3504_CS21syri.unique_insertion.out > 3504_unique_insertion.10k_to_50k.out
cut -f6,7,8 3504_unique_insertion.10k_to_50k.out >3504_unique_insertion.10k_to_50k.bed
bedtools getfasta -fi ../../../3504_map_field/3504.v1.fa -bed 3504_unique_insertion.10k_to_50k.bed -fo 3504_unique_insertion.10k_to_50k.fa

srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/minimap2.sif minimap2 -x asm5 -I100g --cs=long --secondary=no -t 128 Ac.1gChr.fa 3504_unique_insertion.200_to_1k.fa >3504_1chr.200-1k.asm5.paf
