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
import warnings

warnings.filterwarnings('ignore')


pd.set_option("display.max_column",40)

#parser=argparse.ArgumentParser()
#parser.add_argument("-bam","--bam")
#args=parser.parse_args()

#bam_file=args.bam

class bamConverter:

	#Conver bam to bed file
	def ConverAlignment(self,bam):
	  index="samtools index %s"%(bam)
	  os.system(index)
	  samfile = pysam.AlignmentFile(bam, "rc")
	  f=open(bam+"_AligTable.tsv","w")
	  columns=["Refname","RefStart","RefEnd","Readname","ReadLen","ReadStart","ReadEnd","Strand"]
	  f.write("\t".join(columns)+"\n")
	  for read in samfile.fetch():
	    if read.is_unmapped==False:
	      Refname=read.reference_name
	      RefStart=read.reference_start+1
	      RefEnd=read.reference_end
	      Readname=read.query_name
	      ReadLen=read.query_length
	
	      r1=read.query_alignment_start+1
	      r2=read.query_alignment_end
	      if read.is_reverse:
	        Strand="-"
	        ReadStart=ReadLen-r2+1
	        ReadEnd=ReadLen-r1+1
	      else:
	        Strand="+"
	        ReadStart=r1
	        ReadEnd=r2
	      l=[Refname,RefStart,RefEnd,Readname,ReadLen,ReadStart,ReadEnd,Strand]
	      l=[str(i) for i in l]
	      s="\t".join([str(i) for i in l]) 
	      f.write(s+"\n")
	  f.close()

#if __name__=="__main__":
#	myBam = bamConverter(bam_file)
#	myBam.ConverAlignment()
