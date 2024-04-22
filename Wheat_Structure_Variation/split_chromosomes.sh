ls *.fa>>list.txt
#list.txt:
#3504_v1.fa
#Arinalrfor.fa
#CS21.fa
#Fielder.fa
#Jagger.fa
#Julius.fa
#Kariega.fa
#Lancer.fa
#Landmark.fa
#Mace.fa
#Mattis.fa
#Norin61.fa
#Renan.fa
#Stanley.fa

cat list.txt |while read line; 
  do seqkit sort -l -r $line | seqkit head -n 21 > 0.chromo/$line
  samtools faidx $line
  done
###################################################################################################
ls *.fa | while read line
  do     
  filename=$(basename "$line")  # 提取文件名
  prefix=${filename%%.*}  # 提取文件名中第一个"."之前的内容
  for chr in 1A 1B 1D 2A 2B 2D 3A 3B 3D 4A 4B 4D 5A 5B 5D 6A 6B 6D 7A 7B 7D
    do
    mkdir -p "$chr"  # 确保染色体目录存在，如果不存在则创建  
    samtools faidx "$line" "$chr" > "$chr/${prefix}.$chr.fa"
    done 
  done
