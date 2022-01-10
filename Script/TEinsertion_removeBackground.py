import pandas as pd
import os
import argparse
import re
from Bio import SeqIO


pd.set_option("display.max_columns",40)

parser=argparse.ArgumentParser()
parser.add_argument("-f1","--File1")
parser.add_argument("-f2","--File2")
args=parser.parse_args()

f1=args.File1
f2=args.File2



f_1=pd.read_table(f1,header=0)
f_1["coor"]=f_1["RefName_x"]+"_"+f_1["cluster"].apply(str)
print(f_1.shape)
print(f_1.drop_duplicates(["coor"],keep="first").shape)


f_2=pd.read_table(f2,header=0)
f_2["coor"]=f_2["RefName_x"]+"_"+f_2["cluster"].apply(str)
print(f_2.shape)
print(f_2.drop_duplicates(["coor"],keep="first").shape)

f_1=f_1.loc[~f_1["coor"].isin(f_2["coor"])]
print(f_1.shape)
print(f_1.drop_duplicates(["coor"],keep="first").shape)

f_1.to_csv(f1+"_rebg.tsv",index=None,sep="\t")
