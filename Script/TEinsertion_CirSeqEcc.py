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
import pysam
from TEinsertion_bamConverter import bamConverter
import warnings
warnings.filterwarnings('ignore')


pd.set_option("display.max_column",40)

parser=argparse.ArgumentParser()
parser.add_argument("-f","--Circle_File")
args=parser.parse_args()

cf=args.Circle_File

#bamConverter().ConverAlignment(Ta)
def filtReads(file):
	f=pd.read_table(file,header=0,sep="\t")
	f=f.sort_values(["Readname","ReadStart","ReadEnd"])
	#f=f.drop_duplicates(["Readname","RefStart"],keep=False)
	f=f.drop_duplicates(["Readname"],keep="first")
	print(f.shape)
filtReads(cf)
