# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 17:57:45 2019

@author: Administrator
"""

import pandas as pd
import argparse
parser= argparse.ArgumentParser(description="""
说明：
otu表格不带注释信息
fasta名只包含otu名 
""")
parser.add_argument('-i', help="输入文件名", default="rarefac.otu_table.xls")
parser.add_argument('-m', help="分组文件", default="map-group.txt")
parser.add_argument('-f', help="代表序列", default="otu_reps.raw.fasta")
args= parser.parse_args()

def read_fa(faname,falist): # faname,falist = opts.f,d2["seq_name"].tolist()
    '''读取fasta为字典'''
    out={}
    name = seq = ''
    for ln in open(faname,'rt'):
        if ln[0]=='>':
           name=ln.strip()[1:]
           seq=''
        else:
           seq+=ln.strip()
        out[name]=seq
    # 取出
    data=[]
    for i in falist:
        data.append(['>'+i,out[i]])
    return data

def write_fasta(flist,outname):
    #写入表格
    out=open(outname,'w')
    for i in flist:
        out.write('\n'.join(i)+'\n')
    out.close()

data1 = pd.read_csv(args.i,sep="\t",index_col=0)
data2 = pd.read_csv(args.m,sep="\t")
data3 = data1.loc[data1.apply(lambda x: x.sum(), axis=1) > 0,data2[r"#SampleID"].tolist()]
data3.to_csv("sub." + args.i,sep='\t',na_rep='')
               
write_fasta(read_fa(args.f,data3.index.tolist()),"sub." + args.f)

                 
