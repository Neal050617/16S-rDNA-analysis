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
parser.add_argument('-o', help="输出文件名", default='trimed.M2O.fasta')
args= parser.parse_args()

if args.o == '':
    args.o = 'select.'+args.f

def readfa(faname):  #fasta文件名
    '''读取fasta'''
    data = []
    name = seq = ''
    for ln in file(faname,'rt'):
        if ln[0]=='>':
           data.append([name,seq])
           name=ln.strip()
           seq=''
        else:
           seq+=ln.strip()
    data.append([name,seq])
    return data[1:]

def write_data(flist,outname):
    #写入表格
    out=open(outname,'w')
    for i in flist:
        out.write('\n'.join(i)+'\n')
    out.close()

outdata = readfa(args.f)
write_data(outdata,args.o)

