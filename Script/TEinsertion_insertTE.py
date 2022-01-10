import pandas as pd
import os
import argparse
import re
from Bio import SeqIO


pd.set_option("display.max_columns",40)

parser=argparse.ArgumentParser()
parser.add_argument("-f1","--File1")
parser.add_argument("-f2","--File2")
parser.add_argument("-TE","--TEfile")
args=parser.parse_args()

f1=args.File1
f2=args.File2
TE=args.TEfile


records=list(SeqIO.parse(TE,"fasta"))
d=dict(zip([rec.id for rec in records],[len(str(rec.seq)) for rec in records]))

f_1=pd.read_table(f1,header=0)
f_2=pd.read_table(f2,header=0)
f_2["RefLen"]=f_2["RefName_x"].apply(lambda x: d[str(x)])
remove=f_2.loc[(f_2["Junction_1"]<=10) | (f_2["Junction_1"]>=f_2["RefLen"]-10) | ((f_2["Junction_2"]<=10) & (f_2["Junction_2"]>0)) | (f_2["Junction_2"]>=f_2["RefLen"]-10)]
print(f_2.shape)
print(remove.shape)

print(f_1.shape)
f_1=f_1.loc[~f_1["Readname"].isin(list(remove["Readname"]))]
#print(f_1.shape)
#f_1=f_1.loc[~f_1["Readname"].isin(list(f_2["Readname"]))]
#f_1["TEinsertion"]=f_1["Readname"].apply(lambda x: x in list(f_2["Readname"]))
f_1.to_csv(f1+"_Label.tsv",index=None,sep="\t")
