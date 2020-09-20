# -*- coding: utf-8 -*-
"""
Created on Sun Dec 16 13:39:45 2018

@author: Administrator
"""

import argparse
import os
import rpy2.robjects as robjects
from rpy2.robjects import IntVector
from rpy2.robjects import StrVector
import pandas as pd
#from numpy import array
#from numpy.random import normal
#from matplotlib import pyplot
#from itertools import groupby

parser= argparse.ArgumentParser(description="""
说明：
 
""")
parser.add_argument('-i', help="merge.list", default='merge.list')
parser.add_argument('-o', help="输出文件名", default='seq.count.xls')
parser.add_argument('-b', help="分箱数", default=25)
parser.add_argument('-m', help="最大值", default=700)
parser.add_argument('-wd', help="宽度", default=8)
parser.add_argument('-hd', help="高度", default=6)

args= parser.parse_args()

def readfq(fqname):  #fastq文件名 fqname = 'F01A_merged.fastq'
    '''#读取fastq'''
    data = []
    for ln in open(fqname,'rt'):        
        data.append(ln.strip())
    data2 = [len(data[i*4+1]) for i in range(int(len(data)/4))] # i = 1
    return (len(data2),data2)

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

def draw_hist(length,pdfname='hist.pdf',b=25,m=700,wd=8,hd=6): #length = d5
#    分区间统计
#    for z1, z2 in groupby(sorted(length), key=lambda x: x//5):
#        print('{}-{}: {}'.format(z1*5, (z1+1)*5-1, len(list(z2))))
#    matplotlib   作图
#    lenths = array(length)
#    pyplot.hist(x=lenths,bins=50)
#    pyplot.xlabel('Sequence Length')
#    pyplot.xlim(400,500)
#    pyplot.ylabel('Sequence Number')
#    pyplot.title('Sequence Length Distribution')
#    pyplot.show()    
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
     
d1 = pd.read_csv(args.i,sep="\t",header=None).sort_values(by=0)
d2 = []
for i in range(d1.shape[0]): # i =0
    d2.append([readfq(d1.iloc[i,1])[0],readfq(d1.iloc[i,0] + '_merged.fastq')[0],readfa(d1.iloc[i,0] + '_filtered.fasta')])

#zz = readfa(d1.iloc[i,0] + '_filtered.fasta')
# 频率分布直方图
draw_hist([d2[j][2][1] for j in range(len(d2))][0],pdfname='hist.pdf',b=args.b,m=args.m,wd=args.wd,hd=args.hd)

d3 = pd.core.frame.DataFrame({'Sample_ID'  : d1.iloc[:,0].tolist(),
                              'Rawdata'    : [d2[j][0] for j in range(len(d2))],
                              'Merge'      : [d2[j][1] for j in range(len(d2))],
                              'Filter'     : [d2[j][2][0] for j in range(len(d2))],
                              'Total_base' : [d2[j][2][2] for j in range(len(d2))],
                              'Average'    : [d2[j][2][3] for j in range(len(d2))]})

d3 = d3.loc[:,['Sample_ID','Rawdata','Merge','Filter','Total_base','Average']]
d3.to_csv(args.o,sep="\t",index=0)

# 文件夹整理
os.system('''
mkdir data rawdata merge filter_fa filter_fq
mv *_merged.fastq merge
mv *_filtered.fasta filter_fa
mv *_filtered.fastq filter_fq
mv *fastq rawdata
mv merge filter_fa filter_fq rawdata data
''')
print("OK!")


