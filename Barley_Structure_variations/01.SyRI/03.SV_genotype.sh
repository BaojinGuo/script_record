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

ls *INV.vcf|while read line
do 
    run_pipeline.pl -Xmx50g -fork1 -vcf $line -sortPositions -export 2_$line -exportType VCF -exportIncludeAnno false -exportIncludeDepth false -runfork1
done

ls 2_*|while read line 
do 
    bgzip $line
done

ls 2_*|while read line; do tabix -C $line; done
ls 2_*.gz >list.txt
bcftools merge -l list.txt -o INV.vcf

####INS/DEL
ls | while read line
do    
    echo "$line"
    cd "$line"
    grep -v "##" "${line}_Clippersyri.vcf" | grep "#" | sed '1s/$/\tFORMAT\t'"$line"'/' > "${line}_Clipper.INS.vcf"
    grep -v "##" "${line}_Clippersyri.vcf" | grep "#" | sed '1s/$/\tFORMAT\t'"$line"'/' > "${line}_Clipper.DEL.vcf"
    egrep 'INS|DEL' "${line}_Clippersyri.vcf" |awk 'length($4)-length($5)>50'|sed 's/$/\tGT\t1\/1/'| awk -F'\t' -v OFS='\t' '{ $3="."; $6="30"; $8="."; print }'>>"${line}_Clipper.DEL.vcf"
    egrep 'INS|DEL' "${line}_Clippersyri.vcf" |awk 'length($5)-length($4)>50'|sed 's/$/\tGT\t1\/1/'| awk -F'\t' -v OFS='\t' '{ $3="."; $6="30"; $8="."; print }'>>"${line}_Clipper.INS.vcf"
    cd ..
done 

##INS
ls *INS.vcf|while read line
do 
    run_pipeline.pl -Xmx50g -fork1 -vcf $line -sortPositions -export 2_$line -exportType VCF -exportIncludeAnno false -exportIncludeDepth false -runfork1
done

ls 2_*|while read line 
do 
    bgzip $line
done

ls 2_*|while read line; do tabix -C $line; done

ls 2_*.gz >list.txt
bcftools merge -l list.txt -o INS.vcf

##DEL
ls *DEL.vcf|while read line
do 
    run_pipeline.pl -Xmx50g -fork1 -vcf $line -sortPositions -export 2_$line -exportType VCF -exportIncludeAnno false -exportIncludeDepth false -runfork1
done

ls 2_*|while read line 
do 
    bgzip $line
done

ls 2_*|while read line; do tabix -C $line; done

ls 2_*.gz >list.txt
bcftools merge -l list.txt -o DEL.vcf








