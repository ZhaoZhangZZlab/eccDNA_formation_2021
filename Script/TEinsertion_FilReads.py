import pandas as pd
import os
import argparse
import re
from Bio import SeqIO
import pysam

pd.set_option("display.max_columns",40)

parser=argparse.ArgumentParser()
parser.add_argument("-p","--prefix")
parser.add_argument("-bam","--bamFile")
parser.add_argument("-TE","--TEfile")
parser.add_argument("-FuR","--FuR")
parser.add_argument("-f_name","--Fname")
args=parser.parse_args()

bamFile=args.bamFile
pre=args.prefix
TE=args.TEfile
full_r=args.FuR
p_file=args.Fname

records=list(SeqIO.parse(TE,"fasta"))
d=dict(zip([rec.id for rec in records],[len(str(rec.seq)) for rec in records]))

from TEinsertion_bamConverter import bamConverter
#try:
bamConverter().ConverAlignment(bamFile)
#except:
#	print("no alingments file")


def GetTEreads(File):
	f=pd.read_table(File,header=0,sep="\t")
	f["Refname"]=f["Refname"].apply(str)
	f["Reflen"]=f["Refname"].apply(lambda x: d[x])
	full_reads=f.loc[(f["ReadStart"]<=100) & (f["ReadEnd"]>=f["ReadLen"]-100)]
	TE_reads=full_reads.loc[(full_reads["RefStart"]<=100) & (full_reads["RefEnd"]>=full_reads["Reflen"]-100)]
	TE_reads.to_csv(pre+"_LiTEreads.csv",index=None,sep="\t")
	short=full_reads.loc[~full_reads["Readname"].isin(list(TE_reads["Readname"]))]
	short.to_csv(pre+"_shortReads.csv",index=None,sep="\t")
	print("Total: ", len(set(f["Readname"])))
	print("TEreads: ", len(set(TE_reads["Readname"])))
	print("short reads: ", len(set(short["Readname"])))

	candidateReads=f.loc[f["Readname"].apply(lambda x: x not in list(TE_reads["Readname"]) and x not in list(short["Readname"]))]
	print("reads_left: ", len(set(candidateReads["Readname"])))
	candidateReads.to_csv(pre+".candi.tsv",index=None,sep="\t")

#if full_r==0:
GetTEreads(bamFile+"_AligTable.tsv")
#else:
#GetTEreads(bamFile+"_fuR.tsv")
#GetTEreads(p_file)
#GetTEreads(bamFile+"_AligTable.tsv")
