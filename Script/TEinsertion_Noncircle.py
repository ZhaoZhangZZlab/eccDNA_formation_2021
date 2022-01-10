import pandas as pd
import os
import argparse
import re
from Bio import SeqIO


pd.set_option("display.max_columns",40)

parser=argparse.ArgumentParser()
parser.add_argument("-p","--prefix")
parser.add_argument("-fa","--fasta")
args=parser.parse_args()

pre=args.prefix
fa=args.fasta

def GetNonCircle(file):
	f=pd.read_table(file,header=0,sep="\t")
	NC_reads=list(f["Readname"])
	reads=list(SeqIO.parse(fa,"fasta"))
	SeqIO.write((rec for rec in reads if rec.id in NC_reads),file+".fa","fasta")


GetNonCircle(pre+"NC.txt")
