# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 10:39:14 2018

@author: Colin Liu
"""

import sys,os,time,shutil
import argparse
import copy

parser= argparse.ArgumentParser(description="""
将含有taxonomy的OTU表格的门纲目科属种拆分                                
""")
parser.add_argument('-i', help="输入文件名", default='rarefac.otu_taxa_table.xls')
parser.add_argument('-o', help="输出文件名", default='rarefac.otu_taxa_split_table.xls')
parser.add_argument('-x', help="输出文件名", default='rarefac.otu_taxa_split_index_table.xls')
args= parser.parse_args()

table = []
for i in open(args.i,'rt'):
    table.append(i.strip('\n').split('\t'))

adict={'d__':'kingdom','p__':'phylum','c__':'class','o__':'order','f__':'family','g__':'genus','s__':'species'}

liaa=['kingdom','phylum','class','order','family','genus','species','unclassified']

table[0].extend(liaa)

data=[]
for i in range(0,len(table)):     # 列表后补齐缺省值
    data.append(table[i]+[""]*(len(table[0])-len(table[i])))

taxwz=data[0].index('taxonomy')

def tax_split(tax):
    data=tax.split(";")
    return data

def handle_str(tax):
    if "__" in tax:
        tax=tax[tax.index("__")+2:]
    if "(" and ")" in tax:
        tax=tax[:tax.index('(')]
    return tax

for i in range(1,len(data)):
    da=data[i]
    ts=tax_split(da[taxwz])
    for j in range(len(ts)):
        if "__" in ts[j]:
            k1=ts[j][:ts[j].index("__")+2]
            if k1 in adict:
                da[data[0].index(adict[k1])]=handle_str(ts[j])
            else:
                da[-1] = handle_str(ts[j])
        elif 'unclassified' not in ts[j]:
            da[-1]==handle_str(ts[j])

datax = copy.deepcopy(data)
for i in range(1,len(datax)):
    da_c=datax[i]
    ts=tax_split(da_c[taxwz])
    for j in range(len(ts)):
        if "__" in ts[j]:
            k1=ts[j][:ts[j].index("__")+2]
            if k1 in adict:
                da_c[datax[0].index(adict[k1])]=ts[j]
            else:
                da_c[-1] = ts[j]
        elif 'unclassified' not in ts[j]:
            da_c[-1]==ts[j]

output=args.o
out=open(output,'w')
for i in range(len(data)):
    for j in range(len(data[i])):
        out.write(data[i][j]+'\t')
    out.write('\n')
out.close()    

output_index=args.x
dex=open(output_index,'w')
for j in range(len(datax)):
    for k in range(len(datax[j])):
        dex.write(datax[j][k]+'\t')
    dex.write('\n')
dex.close()

print("output:%s"%output)               
print("output:%s"%output_index)               

