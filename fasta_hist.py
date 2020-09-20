# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 10:14:41 2019

@author: Administrator
"""

import argparse
import rpy2.robjects as robjects
from rpy2.robjects import IntVector
from rpy2.robjects import StrVector

parser= argparse.ArgumentParser(description="""
说明：
 
""")
parser.add_argument('-i', help="merge.list", default='meta.fasta')
parser.add_argument('-o', help="输出文件名", default='seq.count.xls')
parser.add_argument('-b', help="分箱数", default=25)
parser.add_argument('-m', help="最大值", default=700)
parser.add_argument('-wd', help="宽度", default=8)
parser.add_argument('-hd', help="高度", default=6)

args= parser.parse_args()


def readfa(faname): #fasta文件名 faname = 'F01A_filtered.fasta'
    '''#读取fasta'''
    data,seq = [],''
    for ln in open(faname,'rt'):
        if ln.startswith(">"):
            data.append(seq)
            seq = ''
            data.append(ln.strip())
        else:
            seq = ln.strip() + seq
    data2 =  [len(data[1:][i*2+1]) for i in range(int(len(data[1:])/2))] # i = 1                 
    return (len(data2),data2,sum(data2),round(sum(data2)/len(data2),2))

def write_data(flist,outname):
    #写入表格
    out=open(outname,'w')
    for i in flist:
        out.write(''.join(i)+'\n')
    out.close()

def draw_hist(length,pdfname='hist.pdf',b=5,m=700,wd=8,hd=6): #length = d5
  
    robjects.globalenv["dd"] = IntVector(length)
    robjects.globalenv["nm"] = StrVector([pdfname])
    robjects.globalenv["b"] = IntVector([b])
    robjects.globalenv["m"] = IntVector([m])
    robjects.globalenv["wd"] = IntVector([wd])
    robjects.globalenv["hd"] = IntVector([hd])
    
    r_script = '''
    library(ape)
    pdf(nm,width=wd,height=hd)
    xcol=seq(0,m,b)

    hist(dd,freq=TRUE,breaks=xcol,col='#228B22',xlab='Sequence Length',ylab='Sequence number',main='Distribution of Sequence Length')
    dev.off()       
    '''
    robjects.r(r_script)

fa1 = readfa(args.i)
draw_hist(fa1[1],b=args.b,m=args.m,wd=args.wd,hd=args.hd)
