for i in 1A 1B 1D 2A 2B 2D 3A 3B 3D 4A 4B 4D 5A 5B 5D 6A 6B 6D 7A 7B 7D; do
    cd "$i" || exit 1  # 进入染色体目录，如果失败则退出脚本
    ls *.fa | while read line; do
        filename=$(basename "$line")
        prefix=${filename%.fa}
        echo '#!/bin/bash
#SBATCH --job-name=Mini
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --account=pawsey0399
module load singularity/3.11.4-nompi
module load samtools/1.15--h3843a85_0
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/minimap2.sif minimap2 -ax asm5 -I100g -f100 -t 128 --eqx CS21.'"$i"'.fa '"$line"' >'"$prefix"'.CS21.sam' >"$prefix.minisam.sh"
    done
    cd ..  # 返回上一级目录
done
