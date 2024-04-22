for i in 1A 1B 1D 2A 2B 2D 3A 3B 3D 4A 4B 4D 5A 5B 5D 6A 6B 6D 7A 7B 7D; do
    cd "$i" || exit 1  # 进入染色体目录，如果失败则退出脚本
    ls *.fa | while read line; do
        filename=$(basename "$line")
        prefix=${filename%.fa}
        mkdir -p "syri/$prefix"  # 创建文件夹，如果不存在则创建
        echo '#!/bin/bash
#SBATCH --job-name=syri
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
source /scratch/pawsey0399/bguo1/software/miniconda/bin/activate syri
srun --export=all -n 1 -c 128 syri -c '"$prefix"'.CS21.sam -r CS21.'"$i"'.fa -q '"$line"' -k -F S --dir syri/'"$prefix"' --prefix '"$prefix"'_CS21 --nc 1' > "$prefix.syri.sh"
    done
    cd ..
done
