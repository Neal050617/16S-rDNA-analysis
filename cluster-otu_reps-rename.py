# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 14:45:29 2018

@author: Administrator
"""

import argparse
import os

def option():
    parser= argparse.ArgumentParser(description="""
    说明：
    python cluster-otu_reps-rename.py -i cluster.fasta
    """)
    parser.add_argument('-i', help="代表序列", default='cluster.fasta')
    parser.add_argument('-t', help="代表序列的注释信息", default='cluster_tax_assignments.txt')
    args= parser.parse_args()
    return args

def write_data(flist,outname):
    #写入表格
    out=open(outname,'w')
    for i in flist:
        out.write('\t'.join(i)+'\n')
    out.close()

def read_assign(filename):  # filename = fn1  
    '''读取table文件 .append or .extend '''  
    data = []
    with open(filename,'r') as df:
        for i,j in enumerate([d.strip().split('\t') for d in df]):
            data.append(["OTU"+str(i+1),j[0].split(";")[0],j[1]])
    write_data([[k[0],k[2]] for k in data],"otu_reps_tax_assignments.txt")
    return data

def read_fa_dict(faname): # faname='cluster.fasta'
    '''读取fasta为字典'''
    out = {}
    name1 = seq = ''
    with open(faname,'r') as df:
        for i in [d.strip() for d in df]:
            if i[0]=='>':#ln=gg[0]
                name1 = i[1:].split(";")[0]
                seq=''
            else:
                seq+=i
                out[name1]=[seq]
    return out

def write_out(d1,d2): # d1,d2 = data1,data2
    Data = []
    for i in d1: # i = d1[0]
        Data.append([r">"+i[0]+"_"+i[1],''.join(d2[i[1]])])
    # 写入文件
    out=open("otu_reps.fasta",'w')
    for i in Data:
        out.write('\n'.join(i)+'\n')
    out.close()
       
if __name__=='__main__':
    opts=option()
    data1 = read_assign(opts.t)
    data2 = read_fa_dict(opts.i)
    write_out(data1,data2)
    print('''输出:
        otu_reps.fasta
        otu_reps_tax_assignments.txt
          ''')
    