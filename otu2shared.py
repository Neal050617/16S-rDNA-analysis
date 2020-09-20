# -*- coding: utf-8 -*-
"""
Created on Fri May 24 14:55:54 2019

@author: Administrator
"""

import pandas as pd
import argparse
parser= argparse.ArgumentParser(description="""
说明：
 
""")
parser.add_argument('-i', help="输入文件名", default="map.otu_table.xls")
parser.add_argument('-o', help="输出文件名", default="otus.shared")
args= parser.parse_args()

data1 = pd.read_csv(args.i,sep="\t",header=None,index_col=0).T
data1.insert(1,"numOtus",data1.shape[1]-1)
data1.insert(0,"label",0.97)

d = [i.replace('OTU','') for i in data1.columns.tolist()[3::]]
dd = len(str(len(d)))
data1.columns = ['label', 'Group', 'numOtus'] + ["Phylo" + (dd - len(d[i-1]))*"0" + str(i) for i in range(1,len(d)+1)]

#data1[''] = ''
data1.to_csv(args.o,sep='\t',index=0,na_rep='')


