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
parser.add_argument("-bam","--bam")
args=parser.parse_args()

bam=args.bam

#records=list(SeqIO.parse(TE,"fasta"))
#TE_l=dict(zip([rec.id for rec in records],[len(str(rec.seq)) for rec in records]))


class bamConverter:

	def __init__(self, bam):
		self.bam = bam

	# convert bam to bed file
	def _bamTobed(self):
  	print("Converting %s to bed"%(self.bam))
  	bamTobed="bedtools bamtobed -i %s -cigar > %s"%(self.bam,self.bam+".bed")
  	os.system(bamTobed)
		return self.bam+".bed"

  def selectBed(self):
		bed = self._bamTobed()
		print("Modify %s file to get read start and read end of the each alignment"%(bed))
		f = pd.read_table(bed,header=None,sep="\t")
	  f["read_length"]=f[6].apply(lambda x: bamElement(x)._getLenfromCigar_())
  	f["cigar_start"]=f[6].apply(lambda x: bamElement(x)._getStartfromCigar_()[1])
  	f["cigar_end"]=f[6].apply(lambda x: bamElement(x)._getEndfromCigar_()[1])
  	f["cigar_position_s"]=f[6].apply(lambda x:int(bamElement(x)._getStartfromCigar_()[0]))
  	f["cigar_position_e"]=f[6].apply(lambda x: int(bamElement(x)._getEndfromCigar_()[0]))
  	f["read_start"]=0
  	f["read_end"]=0

  	f.loc[f["cigar_start"]!="S", "cigar_position_s"]=0
  	f.loc[f["cigar_end"]!="S","cigar_position_e"]=0


  	f.loc[f[5]=="+","read_start"]=f["cigar_position_s"]+1
  	f.loc[f[5]=="+","read_end"]=f["read_length"]-f["cigar_position_e"]
  	f.loc[f[5]=="-","read_start"]=f["cigar_position_e"]+1
  	f.loc[f[5]=="-","read_end"]=f["read_length"]-f["cigar_position_s"]

	  f=f.to_csv(bed+"_modified.bed",index=None,sep="\t")

	class bamElement:

		def __init__(self, cigar):
			self.cigar = cigar

		# get read length from cigar
		def _getLenfromCigar_(self):
  		c=Cigar(self.cigar)
  		return len(c)

		# get the first item from cigar such as (10,S)
		def _getStartfromCigar_(self):
  		c=Cigar(self.cigar)
  		return list(c.items())[0]

		# get the last item from cigar such as (10,S)
		def _getEndfromCigar_(self):
  		c=Cigar(self.cigar)
		  return list(c.items())[-1]

if __name__=="__main__":
	myBam = args.bam
	myBam.selectBed(myBam+".bed")
