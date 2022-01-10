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

import warnings
warnings.filterwarnings('ignore')


pd.set_option("display.max_column",40)

parser=argparse.ArgumentParser()
parser.add_argument("-TEinsert","--TE")
parser.add_argument("-Genomeinsert","--Genome")
parser.add_argument("-pre","--Prefix")
args=parser.parse_args()

TEi=args.TE
Genomeinsert=args.Genome
pre=args.Prefix

def filter(I1,I2):
	f1=pd.read_table(I1,header=0,sep="\t")
	f2=pd.read_table(I2,header=0,sep="\t")
	print(f1)
	print(f1.shape)
	f2=f2.loc[~f2["read"].isin(list(f1["read"]))]
	print(f2)
	print(f2.shape)
	f2.to_csv(I2,index=None,sep="\t")

filter(TEi,Genomeinsert)
