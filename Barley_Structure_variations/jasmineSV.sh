#1.align reads in each sample
#2.Call SVs in each sample with SVcaller
#3.Convert duplications to insertions temporarily for breakpoint refinement and better cross-sample comparison
ls *.vcf|while read line; do jasmine --dup_to_ins --preprocess_only file_list=$line.txt out_file=$line.DtoI.vcf genome_file=Morex.V3.chr.fasta threads=15; done
#4.Refine SVs in each sample with Iris

java -jar /scratch/pawsey0399/bguo1/software/miniconda/pkgs/jasminesv-1.1.5-hdfd78af_0/bin/jasmine.jar --run_iris iris_args="--ngmlr" file_list=irl.list out_file= genome_file=Morex.V3.chr.fasta threads=15 bam_list=bamlist
