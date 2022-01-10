import pandas as pd
import os
import argparse
import re
from Bio import SeqIO


pd.set_option("display.max_columns",40)

parser=argparse.ArgumentParser()
parser.add_argument("-p","--prefix")
parser.add_argument("-TE","--TEfile")
args=parser.parse_args()

pre=args.prefix
TE=args.TEfile


records=list(SeqIO.parse(TE,"fasta"))
d=dict(zip([rec.id for rec in records],[len(str(rec.seq)) for rec in records]))

def Junction(list1,list2):
	n1,n2,n3,n4=list1[0],list1[1],list2[0],list2[1]
	if n3>n2 and n3-n2>100:
		return False
	elif n3<n2 and n3-n2<-500:
		return False
	elif abs(n3-n1)<100 or abs(n2-n4)<100:
		return False
	else:
		return True 


def isCircle(list1,list2):
	
	n1,n2,n3,n4=list1[0],list1[1],list2[0],list2[1]	
	t1,t2,t3,t4=list1[2],list1[3],list2[2],list2[3]
	n_max=max(n1,n2,n3,n4)
	n_min=min(n1,n2,n3,n4)
	d1,d2=list1[4],list2[4]
	t_min=min(t1,t2,t3,t4)
	t_max=max(t1,t2,t3,t4)
	j=Junction(list1,list2)
	readlength=list1[-1]
	if d1!=d2 or j==False:
		return False
	else:
#		if d1=="+" and t_max==t2 and t_min==t3 and n_min<=100 and n_max>=readlength-100:
		if d1=="+" and t_max==t2 and t_min==t3:
			return True
#		if d1=="-" and t_max==t4 and t_min==t1 and n_min<=100 and n_max>=readlength-100:
		if d1=="-" and t_max==t4 and t_min==t1:
			return True
		else:
			return False
	
def CircleType(list1,list2):
	
	n1,n2,n3,n4=list1[0],list1[1],list2[0],list2[1]
	t1,t2,t3,t4=list1[2],list1[3],list2[2],list2[3]
	t_max=max(t1,t2,t3,t4)
	t_min=min(t1,t2,t3,t4)
	n_max=max(n1,n2,n3,n4)
	n_min=min(n1,n2,n3,n4)
	TElength=list1[-2]
	readlength=list1[-1]
	if isCircle(list1,list2) ==False:
		return "NC"
	else:
		if t_max>=TElength-100 and t_min<=100 and n2-n3>=100:
			return "FC_1LTR"
		elif t_max>=TElength-100 and t_min<=100 and n2-n3<100:
			return "FC_2LTR"
		elif t_max<TElength-100 and t_min>100:
			return "PC_nonLTR"
		else:
			return "PC_1LTR"


def getCircle(f):
	d={}
	reads=list(set(list(f["Readname"])))
	for read in reads:
		d[read]="NC"
		sub=f.loc[f["Readname"]==read]
		sub=sub.sort_values(["Refname","ReadStart","ReadEnd"])
		l=list(zip(sub["ReadStart"],sub["ReadEnd"],sub["RefStart"],sub["RefEnd"],sub["Strand"],sub["Reflen"],sub["ReadLen"]))
		i=0
		while i<len(l)-1:
			j=i+1
			while j<len(l):
				list1=l[i]
				list2=l[j]
				cirType=CircleType(list1,list2)
				if cirType!="NC":
					d[read]=cirType
					break
					i=len(l)
				else:
					j=j+1
			i+=1
	f["Circle"]=f["Readname"].apply(lambda x: d[x])
	return f


def GetTEreads(canfile):
	f_c=pd.read_table(canfile,header=0,sep="\t")
	f_c=f_c.sort_values(["Readname","ReadStart","ReadEnd"])
	f_circle=getCircle(f_c)
	f_circle.to_csv(pre+"_circleAnalyze.txt",index=None,sep="\t")
	fc=f_circle.loc[f_circle["Circle"]!="NC"]
	fc.to_csv(pre+"circles.txt",index=None,sep="\t")
	nc=f_circle.loc[f_circle["Circle"]=="NC"]
	nc.to_csv(pre+"NC.txt",index=None,sep="\t")
	print("Circles: ", len(set(fc["Readname"])))
	print("NC: ", len(set(nc["Readname"])))

GetTEreads(pre+".candi.tsv")
