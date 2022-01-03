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
#SBATCH --output="GFP-%j.out"
#SBATCH --error="GFP-%j.err"
## LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

reads1="171107_LW1_aubago_eggs.fastq.chop.fastq";
reads2="20201024_HMS_embryo_0-6h_gDNA.fastq.chop.fastq"

HMS="HMS-Beagle.fa"
Reporter="GFP_seq.fa"
HMSBeagle_GFP="HMS-Beagle_GPF.fa"
genome="dm6.fa"

AllTE="TE_full.fa"

for reads in $reads1 $reads2;
do
name=$(basename $reads);

########################################## Mapping reads to transposon (eg. HMS-Beagle) --> TE reads ##################################################
minimap2 -ax map-ont $HMS $reads -Y -t 16 | samtools view -bS -F 4 | samtools fasta >$name"-"$HMS".sort.fa";

########################################## Mapping TE reads to GFP reporter --> TE+GFP+ reads  ##################################################
minimap2 -ax map-ont $Reporter $name"-"$refName".sort.fa" -Y -t 16 | samtools view -bS -F 4 | samtools fasta > $name"-"$refName".TE+GFP+.fa"

########################################## Remove TE+GFP+ reads from the original chr2L locus ##################################################
minimap2 -ax map-ont "chr2L:6984548-6991487.fa" $name"-"$refName".TE+GFP+.fa" -Y -t 16 | samtools view -bS -f 4 | samtools fasta > $name"-"$refName".TE+GFP+Non2L1.fa";
minimap2 -ax map-ont "chr2L:6998548-7006507.fa" $name"-"$refName".TE+GFP+Non2L1.fa" -Y -t 16 | samtools view -bS -f 4 | samtools fasta > $name"-"$refName".TE+GFP+Non2L2.fa";

########################################### Align reads to HMS-Beagle_GFP consensus  ####################################################################################
minimap2 -ax map-ont $TE2 $name"-"$refName".TE+GFP+Non2L2.fa" -Y -t 16 | samtools view -bS -F 4 | samtools sort > $name"-"$refName".TE+GFP+HMSGFP.bam"
#
########################################### Seperate reads into 3 datasets : 1.short reads; 2.LiTE reads; 3.candidate reads ################################################
python3 /data/zhanglab/Weijia_Su/PythonScrip/TEinsertion_FilReads.py -p $name"-"$refName".TE+GFP+" -bam $name"-"$refName".TE+GFP+HMSGFP.bam" -TE $TE2;
#
########################################### Get circular reads from the candidate reads  ####################################################################################
python3 /data/zhanglab/Weijia_Su/PythonScrip/TEinsertion_circle_lgs.py -p $name"-"$refName".TE+GFP+" -TE $TE2;
#
########################################### Get Insertions from the NC reads ####################################################################################
python3 /data/zhanglab/Weijia_Su/PythonScrip/TEinsertion_Noncircle.py -p $name"-"$refName".TE+GFP+" -fa $name"-"$refName".TE+GFP+Non2L2.fa";
minimap2 -ax map-ont $genome $name"-"$refName".TE+GFP+NC.txt.fa" -Y -t 16 | samtools view -bS | samtools sort > $name"-"$refName".TE+GFP+NC.txt.fa_GENOME.bam"
minimap2 -ax map-ont $TE2 $name"-"$refName".TE+GFP+NC.txt.fa" -Y -t 16 | samtools view -bS | samtools sort > $name"-"$refName".TE+GFP+NC.txt.fa_TE.bam";
python3 /data/zhanglab/Weijia_Su/PythonScrip/TEinsertion_20211117.py -Ta $name"-"$refName".TE+GFP+NC.txt.fa_TE.bam" -Ga $name"-"$refName".TE+GFP+NC.txt.fa_GENOME.bam" -pName $name"-"$refName".TE+GFP+" -TE $TE2
done

