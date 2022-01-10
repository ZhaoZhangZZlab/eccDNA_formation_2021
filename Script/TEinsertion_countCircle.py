import pandas as pd
import os
import argparse
import re
from Bio import SeqIO
import pysam

pd.set_option("display.max_columns",40)

def get(file):
	f=pd.read_table(file,header=0,sep="\t")
	r=f.drop_duplicates(["Readname"],keep="first")
	r=r.groupby(["Circle"],as_index=False).count()[["Circle","Readname"]].sort_values(["Circle"])
	print(r)

dir1="/data/zhanglab/Weijia_Su/2021_fly_ecc/Fig4/rep1-210806_vas_Alt-EJ_circle-seq/"
for i in range(1,6):
	get("barcode0%s_circleAnalyze.txt"%(i))
