import pandas as pd
import os
import argparse
import re
from Bio import SeqIO
import pysam

pd.set_option("display.max_columns",40)

parser=argparse.ArgumentParser()
parser.add_argument("-bam","--bamFile")
args=parser.parse_args()

bamFile=args.bamFile


from TEinsertion_bamConverter import bamConverter
bamConverter().ConverAlignment(bamFile)
