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

reads1="/data/zhanglab/Weijia_Su/Nanopore_Raw_Data/171107_LW1/171107_LW1_aubago_eggs.fastq.chop.fastq";
readName=$(basename $reads1)


#TE="/data/zhanglab/Weijia_Su/CommonDataSet/HMS-Beagle.fa"
Reporter="/data/zhanglab/Weijia_Su/CommonDataSet/GFP_seq.fa"
genome="/data/zhanglab/Weijia_Su/Genomes/Dro/dm6.fa"
TE="/data/zhanglab/Weijia_Su/CommonDataSet/TE_full.fa"

for reads in $reads1;
do
name=$(basename $reads);
for ref in $TE
do
refName=$(basename $ref);
######################################### Mapping reads to transposon (eg. HMS-Beagle) --> TE reads ##################################################
minimap2 -ax map-ont $ref $reads -Y -t 16 | samtools view -bS -F 4 | samtools fasta >$name"-"$refName".sort.fa";
######################################### Mapping TE reads to GFP reporter --> TE+GFP_ reads  ##################################################
minimap2 -ax map-ont $Reporter $name"-"$refName".sort.fa" -Y -t 16 | samtools view -bS -f 4 | samtools fasta > $name"-"$refName".TE+GFP_.fa"
########################################## Align reads to HMS-Beagle consensus  ####################################################################################
minimap2 -ax map-ont $TE $name"-"$refName".TE+GFP_.fa" -Y -t 16 | samtools view -bS -F 4 | samtools sort > $name"-"$refName".TE+GFP_HMS.bam"

########################################## Seperate reads into 3 datasets : 1.short reads; 2.LiTE reads; 3.candidate reads ################################################
python3 /data/zhanglab/Weijia_Su/PythonScrip/TEinsertion_FilReads.py -p $name"-"$refName".TE+GFP_" -bam $name"-"$refName".TE+GFP_HMS.bam" -TE $TE;

########################################## Get circular reads from the candidate reads  ####################################################################################
python3 /data/zhanglab/Weijia_Su/PythonScrip/TEinsertion_circle_lgs.py -p $name"-"$refName".TE+GFP_" -TE $TE;

########################################## Get Insertions from the NC reads ####################################################################################
python3 /data/zhanglab/Weijia_Su/PythonScrip/TEinsertion_Noncircle.py -p $name"-"$refName".TE+GFP_" -fa $name"-"$refName".TE+GFP_.fa";

minimap2 -ax map-ont $genome $name"-"$refName".TE+GFP_NC.txt.fa" -Y -t 16 | samtools view -bS | samtools sort > $name"-"$refName".TE+GFP_NC.txt.fa_GENOME.bam"
minimap2 -ax map-ont $TE $name"-"$refName".TE+GFP_NC.txt.fa" -Y -t 16 | samtools view -bS | samtools sort > $name"-"$refName".TE+GFP+NC.txt.fa_TE.bam";

for te in /data/zhanglab/Weijia_Su/CommonDataSet/TE_full/*;
do
tename=$(basename $te);
echo $tename
python3 /data/zhanglab/Weijia_Su/PythonScrip/TEinsertion_20211203.py -Ta $name"-"$refName".TE+GFP+NC.txt.fa_TE.bam"; -Ga $name"-"$refName".TE+GFP_NC.txt.fa_GENOME.bam" -pName  $name"-"$refName-$tename".TE+GFP_" -TE $te;
done
done
done
