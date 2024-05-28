####INV
grep -v "##" Eth69_Clippersyri.vcf |grep "#"|sed '1s/$/\tGT\tEth69/' >Eth69_Clipper.INV.vcf
grep "<INV>" Eth69_Clippersyri.vcf | sed 's/$/\tGT\t1\/1/' >> Eth69_Clipper.INV.vcf

ls | while read line; do
    echo "$line"
    cd "$line"
    
    # 处理 VCF 文件头部（第一行）
    grep -v "##" "${line}_Clippersyri.vcf" | grep "#" | sed '1s/$/\tFORMAT\t'"$line"'/' > "${line}_Clipper.INV.vcf"
    
    # 处理包含 <INV> 的数据行，并进行列替换
    grep "<INV>" "${line}_Clippersyri.vcf" | sed 's/$/\tGT\t1\/1/' | awk -F'\t' -v OFS='\t' '{ $3="."; $6="30"; $8="."; print }' >> "${line}_Clipper.INV.vcf"
    
    cd ..
done

