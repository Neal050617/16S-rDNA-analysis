# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 11:04:18 2018

@author: Colin Liu
"""

import argparse
import itertools
parser= argparse.ArgumentParser(description="""
说明：
 
""")
parser.add_argument('-f', help="输入文件名", default='trimed.fasta')
parser.add_argument('-s', help="from", default='select.list')
parser.add_argument('-o', help="输出文件名", default='')
args= parser.parse_args()

if args.o == '':
    args.o = 'select.'+args.f
    
def read_fa_dict(faname):
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
    return out

def select_fa(fa,lt):
    data=[]
    for i in lt:
        data.append(['>'+i,fa[i]])
    return data

def write_data(flist,outname):
    #写入表格
    out=open(outname,'w')
    for i in flist:
        out.write(''.join(i)+'\n')
    out.close()
    
data = []
with open(args.s,'rt') as df:
    for ln in df:
        data.append(ln.strip())
             
outdata = read_fa_dict(args.f)
select_data = list(itertools.chain.from_iterable(select_fa(outdata,data)))
write_data(select_data,args.o)

