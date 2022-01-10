#!/bin/bash

## Modify this job script accordingly and submit with the command:
##    sbatch HPC.sbatch
#SBATCH --nodes=1   # number of nodes
#BATCH --ntasks-per-node=1   # 16 processor core(s) per node
#SBATCH --job-name='171107'
#SBATCH --mem=100000
#SBATCH --partition="all"
#SBATCH --mail-user=weijia.su@duke.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output="13-%j.out"
#SBATCH --error="13-%j.err"
## LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE


reads1="/data/zhanglab/Weijia_Su/Nanopore_Raw_Data/211013-fly-vasKD-F2emb-genomicDNA/barcode13_vas-white.fq";
readName=$(basename $reads1)


TE="/data/zhanglab/Weijia_Su/CommonDataSet/HMS-Beagle.fa"
Reporter="/data/zhanglab/Weijia_Su/CommonDataSet/GFP_seq.fa"
TE2="/data/zhanglab/Weijia_Su/CommonDataSet/HMS-Beagle_GPF.fa"
genome="/data/zhanglab/Weijia_Su/Genomes/Dro/dm6.fa"
AllTE="/data/zhanglab/Weijia_Su/CommonDataSet/TE_full.fa"
TEName=$(basename $TE)

for i in {1..5}
do
reads="/data/zhanglab/Weijia_Su/Nanopore_Raw_Data/210806_vas_Alt-EJ_circle-seq/barcode0"$i".fastq.chop.fastq";
name=$(basename $reads);
for te in HMS-Beagle blood 3S18;
do

python3 /data/zhanglab/Weijia_Su/PythonScrip/TEinsertion_FilReads.py -p $name"_"$te -bam $name"_"$te".bam" -TE "/data/zhanglab/Weijia_Su/CommonDataSet/TE_full/"$te".fasta";

python3 /data/zhanglab/Weijia_Su/PythonScrip/TEinsertion_circle.py -p $name"_"$te -TE "/data/zhanglab/Weijia_Su/CommonDataSet/TE_full/"$te".fasta";

done
done
