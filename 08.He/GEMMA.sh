LD_modify.vcf 为vcf输入文件

header_file 为表型数据表头文件，和表型数据对应

pheno_file 为全部表型数据，只有数值

vcftools --vcf LD_modify.vcf --plink --out LD.geno

plink --file LD.geno --make-bed --out LD.gemma.input

##########run all phenotypes by gemma
#!/bin/bash

# 文件路径
header_file="pheno.header"
pheno_file="gemma.pheno.txt"
input_bfile="LD.gemma.input"

# 读取表型名
read -r -a pheno_names < "$header_file"

# 获取表型数量
num_pheno=${#pheno_names[@]}

# 循环处理每个表型
for i in $(seq 0 $((num_pheno-1))); do
    # 提取当前表型名称
    pheno="${pheno_names[i]}"

    # 输出文件名
    output_file="${pheno}"

    # 运行 GEMMA 的第一个命令
    gemma -bfile "$input_bfile" -gk 2 -p "$pheno_file" -o "$output_file"

    # 运行 GEMMA 的第二个命令
    gemma -bfile "$input_bfile" -k "output/${pheno}.sXX.txt" -lmm 4 -p "$pheno_file" -o "$output_file"
done

#####merge all phenotypes
for file in *assoc.txt; do
    # 获取文件名（第一个 . 前的内容）
    base_name="${file%%.*}"

    # 在每个文件的每一行前面添加文件名，并追加到输出文件中
    awk -v name="$base_name" '{print name, $0}' "$file" >> merge.assoc.gemma.txt
done
