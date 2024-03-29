####EDTA+REPEATMASKER annotate TE
###EDTA2.2.0 INSTALL
singularity pull docker://quay.io/biocontainers/edta:2.2.0--hdfd78af_1
###download congeneric species' cds 
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/pisum_sativum/cds/Pisum_sativum.Pisum_sativum_v1a.cds.all.fa.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/medicago_truncatula/cds/Medicago_truncatula.MedtrA17_4.0.cds.all.fa.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/trifolium_pratense/cds/Trifolium_pratense.Trpr.cds.all.fa.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-58/fasta/glycine_max/cds/Glycine_max.Glycine_max_v2.1.cds.all.fa.gz
gunzip *.gz
cat *.fa >homo.cds
srun --export=all -n 1 -c 128 cd-hit -i homo.cds -o homo.cdhit.cds -c 0.8 -n 4 -T 128 -M 0

###EDTA RUN
srun --export=all -n 1 -c 128 singularity exec /scratch/pawsey0399/bguo1/edta_2.2.0--hdfd78af_1.sif EDTA.pl --genome /scratch/pawsey0399/bguo1/0.assembly/01.hifi_assembly/S1_HIFI_RESULT/S1_hifi.asm.bp.p_ctg.fa --species others --step all --cds /scratch/pawsey0399/bguo1/0.assembly/04.TE_annotation/Ref_cds/homo.cdhit.cds --sensitive 1 --anno 1 --evaluate 1 -t 128
##–genome 基因组文件
##–cds 提供这个种或近缘种的CDS序列（不能包括introns和UTR），用于最终过滤。
##–rmout 提供其他软件做的同源TE注释（比如repeatmasker的.out文件），如果不提供则默认使用EDTA - library用于masking。
##–species [Rice|Maize|others]三种可选
##–step [all|filter|final|anno]
##-sensitive: 是否用RepeatModeler分析剩下的TE，默认是0，不要。RepeatModeler会增加运行时间。
##-anno: 是否在构建TE文库后进行全基因组预测，默认是0.
##-evalues: 默认是0，需要同时设置-anno 1才能使用。能够评估注释质量，但会显著增加分析时间。
##–overwrite默认是0，设定为1会删除已有结果重新运行，建议保持默认，运行中断可以继续运行。

###LAI 
LAI -t 128 -genome S1/S1_hifi.asm.bp.p_ctg.fa -intact S1/S1_hifi.asm.bp.p_ctg.fa.mod.EDTA.raw/LTR/S1_hifi.asm.bp.p_ctg.fa.mod.pass.list -all S1/S1_hifi.asm.bp.p_ctg.fa.mod.EDTA.anno/S1_hifi.asm.bp.p_ctg.fa.mod.out >S1.LAI.out
LAI -t 128 -genome S2/S2_hifi.asm.bp.p_ctg.fa -intact S2/S2_hifi.asm.bp.p_ctg.fa.mod.EDTA.raw/LTR/S2_hifi.asm.bp.p_ctg.fa.mod.pass.list -all S2/S2_hifi.asm.bp.p_ctg.fa.mod.EDTA.anno/S2_hifi.asm.bp.p_ctg.fa.mod.out >S2.LAI.out

###Modify insertion time, because of If the parameter --u is not set in edta, the default insertion time 1.3e-8 for gramineae is used to calculate the insertion time, so we should transform it to the specific crop. 
###git clone https://gitee.com/zheng-tianze/time_convert.py.git
###time_convert.py input.txt 3e-10; Among them, the input. TXT file for input, generally for EDTA output file XX.fa.mod.EDTA.raw/LTR/XX.fa.mod.pass.list. 3e-10 is the mutation rate of the species, here changed to the mutation rate of the species itself.
time_convert.py ../S1/S1_hifi.asm.bp.p_ctg.fa.mod.EDTA.raw/LTR/S1_hifi.asm.bp.p_ctg.fa.mod.pass.list 7e-9
time_convert.py ../S2/S2_hifi.asm.bp.p_ctg.fa.mod.EDTA.raw/LTR/S2_hifi.asm.bp.p_ctg.fa.mod.pass.list 7e-9




