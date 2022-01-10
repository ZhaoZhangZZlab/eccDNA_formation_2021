import pandas as pd
import os
import os.path
import argparse
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time
from cigar import Cigar
from collections import Counter
from TEinsertion_bamConverter import bamConverter
import numpy as np
import warnings
warnings.filterwarnings('ignore')


pd.set_option("display.max_column",40)

parser=argparse.ArgumentParser()
parser.add_argument("-Bed1","--bed1")
parser.add_argument("-Bed2","--bed2")
parser.add_argument("-bam","--bamFile")
parser.add_argument("-file","--File")
args=parser.parse_args()

b1=args.bed1
b2=args.bed2
bamfile=args.bamFile

ali_file=args.File

def map_ratio(sub_f):
	r_len=int(list(sub_f["ReadLen"])[0])
	l=zip(sub_f["ReadStart"],sub_f["ReadEnd"])
	a=np.array([0]*r_len)
	for i in l:
		a[i[0]:i[1]+1]=1
	return list(a).count(1)
	


def getChimeric_reads(bamfile):
	
	bamConverter().ConverAlignment(bamfile)
	f=pd.read_table(bamfile+"_AligTable.tsv",header=0,sep="\t")
	#f=pd.read_table(bamfile,header=0,sep="\t")
	f=f.sort_values(["Readname","ReadStart","ReadEnd"])
	s1=f.drop_duplicates(["Readname"],keep="first").shape[0]
	f["m"]=0
	for r in set(f["Readname"]):
		sub_f=f.loc[f["Readname"]==r]
		p=map_ratio(sub_f)
		f.loc[f["Readname"]==r,"m"]=p
	f["p"]=f["m"]/f["ReadLen"]
	f=f.loc[f["p"]>=0.9]
	f.to_csv(bamfile+"_fuR.tsv",index=None,sep="\t")
	s2=f.drop_duplicates(["Readname"],keep="first").shape[0]
	print("Total reads %s; Fully aligned reads %s (%s)"%(s1,s2,round(s2/s1*100,2)))
getChimeric_reads(bamfile)
#getChimeric_reads(ali_file)



def getChimeric_two(file1,file2):
	f1=pd.read_table(file1,header=0,sep="\t")
	f2=pd.read_table(file2,header=0,sep="\t")
	f=pd.merge(f1,f2,on=["3"],how="inner")
	#f=f.loc[f["0_y"]!="HMS-Beagle"]
	f=f.loc[f["0_x"]!=f["0_y"]]
	f=f.loc[(f["read_end_x"]<=f["read_start_y"]) | (f["read_start_x"]>=f["read_end_y"])]
	f["Map1_len"]=f["read_end_x"]-f["read_start_x"]+1
	f["Map2_len"]=f["read_end_y"]-f["read_start_y"]+1
	f=f[["3","read_length_x","0_x","1_x","2_x","5_x","read_start_x","read_end_x","Map1_len","0_y","1_y","2_y","5_y","read_start_y","read_end_y","Map2_len"]]
	f.columns=["read","read_length","Ref1","Ref1_s","Ref1_e","Ref1_dir","read_start_ref1","read_end_ref1","Map1_len","Ref2","Ref2_s","Ref2_e","Ref2_dir","read_start_ref2","read_end_ref2","Map2_len"]
	f=f.sort_values(["read","Map1_len","Map2_len"],ascending=[True,False,False])
	r=f.drop_duplicates(["read"],keep="first")
#getChimeric_two(b1,b2)
