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
parser.add_argument("-Ta","--TE_bam")
parser.add_argument("-Ga","--Genome_bam")
parser.add_argument("-pName","--Prefix")
parser.add_argument("-TE","--TE_file")
parser.add_argument("-raw","--rawReads")
args=parser.parse_args()

Ta=args.TE_bam
Ga=args.Genome_bam
pName=args.Prefix
TE=args.TE_file
RawReads=args.rawReads

records=list(SeqIO.parse(TE,"fasta"))
TE_l=dict(zip([rec.id for rec in records],[len(str(rec.seq)) for rec in records]))

# get read length from cigar
def __getLenfromCigar__(cigar):
	c=Cigar(cigar)
	return len(c)

# get the first item from cigar such as (10,S)
def __getStartfromCigar__(cigar):
	c=Cigar(cigar)
	return list(c.items())[0]

# get the last item from cigar such as (10,S)
def __getEndfromCigar__(cigar):
	c=Cigar(cigar)
	return list(c.items())[-1]

# convert bam to bed file
def bamTobed(bam):
	print("Converting %s to bed"%(bam))
	bamTobed="bedtools bamtobed -i %s -cigar > %s"%(bam,bam+".bed")
	os.system(bamTobed)

bamTobed(Ta)	
bamTobed(Ga)


def SelectBed(bed):
	print("Modify %s file to get read start and read end of the each alignment"%(bed))

	f=pd.read_table(bed,header=None,sep="\t")
	f["read_length"]=f[6].apply(lambda x: __getLenfromCigar__(x))
	f["cigar_start"]=f[6].apply(lambda x: __getStartfromCigar__(x)[1])
	f["cigar_end"]=f[6].apply(lambda x: __getEndfromCigar__(x)[1])
	f["cigar_position_s"]=f[6].apply(lambda x:int(__getStartfromCigar__(x)[0]))
	f["cigar_position_e"]=f[6].apply(lambda x: int(__getEndfromCigar__(x)[0]))
	f["read_start"]=0
	f["read_end"]=0
	
	f.loc[f["cigar_start"]!="S", "cigar_position_s"]=0
	f.loc[f["cigar_end"]!="S","cigar_position_e"]=0
	

	f.loc[f[5]=="+","read_start"]=f["cigar_position_s"]+1
	f.loc[f[5]=="+","read_end"]=f["read_length"]-f["cigar_position_e"]
	f.loc[f[5]=="-","read_start"]=f["cigar_position_e"]+1
	f.loc[f[5]=="-","read_end"]=f["read_length"]-f["cigar_position_s"]
	
	f=f.to_csv(bed+"_modified.bed",index=None,sep="\t")

SelectBed(Ta+".bed")   
SelectBed(Ga+".bed")

def combineTEandGenome(TEbed,GenomeBed):
	print("Combining TE and Genome mappings")
	TE=pd.read_table(TEbed,header=0,sep="\t")
	Genome=pd.read_table(GenomeBed,header=0,sep="\t")
	combine=pd.merge(TE,Genome,how="inner",on="3")
	combine["0_x"]=combine["0_x"].apply(str)
	combine["TE_length"]=combine["0_x"].apply(lambda x:TE_l[x])
	combine=combine[["3","0_x","1_x","2_x","TE_length","5_x","6_x","0_y","1_y","2_y","5_y","6_y","read_start_x","read_end_x","read_start_y","read_end_y","read_length_y"]]
	combine.columns=["read","TE","TE_start","TE_end","TE_length","TEstrand","TEcigar","Chr","Chr_start","Chr_end","Chrstrand","Chrcigar","read-TE-start","read-TE-end","read-Chr-start","read-Chr-end","read_length"]
	combine.to_csv(pName+"_combine.tsv",index=None,sep="\t")	

combineTEandGenome(Ta+".bed_modified.bed",Ga+".bed_modified.bed")


def Filter(combinedFile):
	print("Filtering")
	f=pd.read_table(combinedFile,header=0,sep="\t")
	#Remove Overlaps
	f["o"]=((f["read-Chr-start"]<f["read-TE-start"]) & (f["read-Chr-end"]>f["read-TE-start"]+100)) | ((f["read-Chr-end"]>f["read-TE-end"]) & (f["read-TE-end"]-f["read-Chr-start"]>100))
	#f["o"]=((f["read-Chr-start"]<=f["read-TE-start"]) & (f["read-Chr-end"]>=f["read-TE-start"]))
	f=f.loc[f["o"]==False]
	
	#Get TE ends
	###f=f.loc[((f["TE_start"]<=10) & (f["TE_end"]>=f["TE_length"]/2)) | ((f["TE_end"]>=f["TE_length"]-10) & (f["TE_start"]<=f["TE_length"]/2))]	
	f=f.loc[((f["TE_start"]<=10) & (f["TE_end"]>=500)) | ((f["TE_end"]>=f["TE_length"]-10) & (f["TE_start"]<=f["TE_length"]-500))] 
	f["dis1"]=abs(f["read-Chr-end"]-f["read-TE-start"])
	f["dis2"]=abs(f["read-Chr-start"]-f["read-TE-end"])
	f=f.loc[(f["dis1"]<=100) | (f["dis2"]<=100)]
	f=f.drop(["dis1","dis2","o"],axis=1)
	f.to_csv(pName+"_filtered.tsv",index=None,sep="\t")
	
Filter(pName+"_combine.tsv")

def get_insertions(filteredFile):
	f=pd.read_table(filteredFile,header=0,sep="\t")
	infor_col=["read","TE","read-TE-start","read-TE-end"]
	f["infor"]=f[infor_col].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
	r=f.groupby(["infor"],as_index=False).count()
	r1=r.loc[r["TE_length"]==1]
	r2=r.loc[r["TE_length"]==2]
	f1=f.loc[f["infor"].isin(list(r1["infor"]))]
	f2=f.loc[f["infor"].isin(list(r2["infor"]))]
	
	f1=f1.drop(["infor"],axis=1)
	f2=f2.drop(["infor"],axis=1)
		
	f2=f2.sort_values(["read","TE","TE_start","TE_end","read-Chr-start","read-Chr-end"])

	f2_1=f2.drop_duplicates(["read","TE","TE_start","TE_end"],keep="first")
	f2_2=f2.drop_duplicates(["read","TE","TE_start","TE_end"],keep="last")
	new_f2=f2_1.merge(f2_2,on=["read","TE","TE_start","TE_end"],how="inner")
	new_f2=new_f2.drop(["TE_length_y","TEstrand_y","TEcigar_y","read-TE-start_y","read-TE-end_y"],axis=1)
	new_f2.to_csv(pName+"_org_double.tsv",index=None,sep="\t")
	f1.to_csv(pName+"_org_single.tsv",index=None,sep="\t")

get_insertions(pName+"_filtered.tsv")

def __Select2__(x,y,z,w):
	l=sorted([x,y,z,w])
	return l[1]

def __Select3__(x,y,z,w):
	l=sorted([x,y,z,w])
	return l[2]


def Filter_double(double_df):
	f=double_df
	print(f[0:10])
	f=f.sort_values(["read","TE","TE_start","TE_end"])
	f=f.loc[f["read_length_y"]>=500]
	f=f.loc[(f["read-TE-end_x"]-f["read-TE-start_x"]>=100) & (f["read-Chr-end_y"]-f["read-Chr-start_y"]>=100)]
	f=f[["read","TE","TE_start","TE_end","TE_length_x","TEstrand_x","Chr_x","Chr_start_x","Chr_end_x","Chrstrand_x","Chr_y","Chr_start_y","Chr_end_y","Chrstrand_y","read-Chr-start_x","read-Chr-end_x","read-TE-start_x","read-TE-end_x","read-Chr-start_y","read-Chr-end_y","read_length_y"]]
	f["dis1"]=abs(f["read-Chr-end_x"]-f["read-TE-start_x"])
	f["dis2"]=abs(f["read-Chr-start_y"]-f["read-TE-end_x"])	
	f=f.loc[(f["dis1"]<=100) & (f["dis2"]<=100)]
	f["sum"]=f["read-Chr-end_x"]-f["read-Chr-start_x"]+1+f["read-TE-end_x"]-f["read-TE-start_x"]+1+f["read-Chr-end_x"]-f["read-Chr-start_x"]+1
	f=f.loc[f["sum"]/f["read_length_y"]>=0.8]
	f=f.loc[(f["read-Chr-start_x"]<=1000) | (f["read-Chr-end_y"]>=f["read_length_y"]-1000)]
	return f
	
double_df=pd.read_table(pName+"_org_double.tsv")
double_df=Filter_double(double_df)
insertion=double_df.loc[double_df["Chr_x"]==double_df["Chr_y"]]
trans=double_df.loc[double_df["Chr_x"]!=double_df["Chr_y"]]
trans.to_csv(pName+"_transloc.tsv",sep="\t",index=None)

def re_countDouble(double_df):
	f=double_df
	f["junctionStart"]=f.apply(lambda x : __Select2__(x["Chr_start_x"],x["Chr_end_x"],x["Chr_start_y"], x["Chr_end_y"]), axis=1)
	f["junctionEnd"]=f.apply(lambda x : __Select3__(x["Chr_start_x"],x["Chr_end_x"],x["Chr_start_y"], x["Chr_end_y"]), axis=1)
	I=f[["Chr_x","junctionStart","junctionEnd","read","TE","TE_length_x","read_length_y","read-Chr-start_x","read-Chr-end_x","read-TE-start_x","read-TE-end_x","read-Chr-start_y","read-Chr-end_y"]]
	deletion=f.loc[f["junctionEnd"]-f["junctionStart"]>100]
	deletion["D"]=deletion["junctionEnd"]-deletion["junctionStart"]
	deletion=deletion.sort_values(["D"],ascending=[False]).drop_duplicates(["read"],keep="first")
	deletion.to_csv(pName+"_deletion.tsv",index=None,sep="\t")
	I["Start_cluster"]=I["junctionStart"].apply(lambda x: round(x,-3))
	I["End_cluster"]=I["junctionEnd"].apply(lambda x: round(x,-3))
	I["j"]="double"
	I=I[["Chr_x","junctionStart","junctionEnd","read","TE","j","Start_cluster","End_cluster","TE_length_x","read_length_y","read-Chr-start_x","read-Chr-end_x","read-TE-start_x","read-TE-end_x","read-Chr-start_y","read-Chr-end_y"]]
	I.columns=["Chr","junctionStart","junctionEnd","read","TE","j","Start_cluster","End_cluster","TE_length","read_length","read-Chr-start","read-Chr-end","read-TE-start","read-TE-end","read-Chr-start_y","read-Chr-end_y"]
	return I

I2=re_countDouble(insertion)


def Re_countSingle(orgFile):
	f=pd.read_table(orgFile,header=0,sep="\t")
	f=f.loc[f["read_length"]>=500]
	f=f.loc[(f["read-TE-end"]-f["read-TE-start"]>=100) & (f["read-Chr-end"]-f["read-Chr-start"]>=100)]
	f["sum"]=f["read-Chr-end"]-f["read-Chr-start"]+1+f["read-TE-end"]-f["read-TE-start"]+1
	f=f.loc[f["sum"]/f["read_length"]>=0.8]
	f["dis1"]=abs(f["read-Chr-end"]-f["read-TE-start"])
	f["dis2"]=abs(f["read-Chr-start"]-f["read-TE-end"])
	f["junctionStart"]=0
	f.loc[(f["Chrstrand"]=="-")&(f["dis1"]<=100),"junctionStart"]=f["Chr_start"]
	f.loc[(f["Chrstrand"]=="-")&(f["dis2"]<=100),"junctionStart"]=f["Chr_end"]
	f.loc[(f["Chrstrand"]=="+")&(f["dis1"]<=100),"junctionStart"]=f["Chr_end"]
	f.loc[(f["Chrstrand"]=="+")&(f["dis2"]<=100),"junctionStart"]=f["Chr_start"]
	f["min"]=0
	f["max"]=0
	f.loc[f["dis1"]<=100,"min"]=f["read-Chr-start"]
	f.loc[f["dis1"]<=100,"max"]=f["read-TE-end"]
	f.loc[f["dis2"]<100,"min"]=f["read-TE-start"]
	f.loc[f["dis2"]<100,"max"]=f["read-Chr-end"]
	f=f.loc[(f["min"]<=1000) | (f["max"]>=f["read_length"]-1000)]
	f=f.drop_duplicates(["read","read-TE-start","read-TE-end"],keep="first")
	f=f.sort_values(["read","read-TE-start","read-TE-end"])
	f["Start_cluster"]=f["junctionStart"].apply(lambda x: round(x,-3))
	f["End_cluster"]=0
	f["junctionEnd"]=0
	f["j"]="single"
	f["read-Chr-start_y"]=0
	f["read-Chr-end_y"]=0
	I=f[["Chr","junctionStart","junctionEnd","read","TE","j","Start_cluster","End_cluster","TE_length","read_length","read-Chr-start","read-Chr-end","read-TE-start","read-TE-end","read-Chr-start_y","read-Chr-end_y"]]
	return I

I1=Re_countSingle(pName+"_org_single.tsv")


final=I2.append(I1)
#final=I1
final=final.drop_duplicates(["read","Start_cluster"],keep="first")
print(final[0:10])
final_inser=final.groupby(["Chr","Start_cluster"],as_index=False).count().sort_values(["read"],ascending=[False])
print(final_inser.shape)
final.to_csv(pName+"_Final_insertion.tsv",index=None,sep="\t")
