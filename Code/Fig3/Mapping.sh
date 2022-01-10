#!/bin/bash

## Modify this job script accordingly and submit with the command:
##    sbatch HPC.sbatch
#SBATCH --nodes=1   # number of nodes
#BATCH --ntasks-per-node=1   # 16 processor core(s) per node
#SBATCH --job-name='200116'
#SBATCH --mem=100000
#SBATCH --partition="all"
#SBATCH --mail-user=weijia.su@duke.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output="Mapping-%j.out"
#SBATCH --error="Mapping-%j.err"
## LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE


reads1="/data/zhanglab/Weijia_Su/Nanopore_Raw_Data/200116_AubAgo3_eccDNA/200116_AubAgo3_eccDNA.fastq.chop.fastq"
reads2="/data/zhanglab/Weijia_Su/Nanopore_Raw_Data/200117_white_eccDNA/200117_white_eccDNA.fastq.chop.fastq"

TE="/data/zhanglab/Weijia_Su/CommonDataSet/HMS-Beagle.fa"
Reporter="/data/zhanglab/Weijia_Su/CommonDataSet/GFP_seq.fa"
TE2="/data/zhanglab/Weijia_Su/CommonDataSet/HMS-Beagle_GPF.fa"
genome="/data/zhanglab/Weijia_Su/Genomes/Dro/dm6.fa"
AllTE="/data/zhanglab/Weijia_Su/CommonDataSet/TE_full.fa"
TEName=$(basename $TE)

for reads in $reads1 $reads2;
do
name=$(basename $reads);
for ref in $AllTE;
do
refName=$(basename $ref);
minimap2 -ax map-ont $ref $reads -Y -t 16 | samtools view -bS | samtools sort > $name"-"$refName".bam";

done
done
