# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 15:05:12 2018

@author: Administrator
"""

import argparse
import os
parser= argparse.ArgumentParser(description="""
说明：
filter_fasta,默认args.r为TRUE修改序列名；否则请改成FALSE
""")
parser.add_argument('-i', help="输入文件名", default='merge.list')
parser.add_argument('-m', help="最小重合区域的大小", default=50)
parser.add_argument('-l', help="最小拼接长度", default=400)
parser.add_argument('-r', help="rename", default="TRUE")
args= parser.parse_args()

def readtable(filename):    
    '''读取table文件 .append or .extend '''  
    data = []
    for ln in open(filename,'rt'):
        data.append(ln.strip().split('\t'))
    return data

data1 = readtable(args.i)

for i in data1 :#i=data1[0]
    os.system('''
              usearch11.0.667_i86linux32 -fastq_mergepairs %s -reverse %s -fastq_minovlen %s -fastq_minmergelen %s -fastq_trunctail 2 -fastqout %s
              ''' % (i[1],i[2],args.m,args.l,i[0]+"_merged.fastq"))
    if args.r == "TRUE":
        os.system('''
               usearch11.0.667_i86linux32 -fastq_filter %s -fastaout %s -fastqout %s -relabel %s -fastq_eeout -fastq_maxee 1.0
               ''' % (i[0]+"_merged.fastq",i[0]+"_filtered.fasta",i[0]+"_filtered.fastq",i[0]+"_1000000"))
    else:
        os.system('''
               usearch11.0.667_i86linux32 -fastq_filter %s -fastaout %s -fastqout %s -fastq_eeout -fastq_maxee 1.0
               ''' % (i[0]+"_merged.fastq",i[0]+"_filtered.fasta",i[0]+"_filtered.fastq"))
        
    
#    os.system('''
#usearch11.0.667_i86linux32 -fastq_filter %s -fastqout %s -relabel %s -fastq_eeout -fastq_maxee 1.0
#''' % (i[0]+"_merged.fastq",i[0]+"_filtered.fastq",i[0]+"_1000000"))
    
print("complete!")