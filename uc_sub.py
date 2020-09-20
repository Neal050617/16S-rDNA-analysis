# -*- coding: utf-8 -*-
"""
Created on Thu Dec 27 18:16:48 2018

@author: Administrator
"""

import argparse
import pandas as pd

def option():
    parser= argparse.ArgumentParser(description="""
    说明：
    """)
    parser.add_argument('-i', help="usearch uparse file", default='map.uc')
    parser.add_argument('-g', help="subsample groups", default='otu.subsample.groups')
    args= parser.parse_args()
    return args

def sub_uc(data,gp): # data,gp = opts.i,opts.g
    uc1 = pd.read_csv(data,sep='\t',encoding='utf8',header=None)
    uc1.index = uc1.iloc[:,8]
    subg = pd.read_csv(gp,sep='\t',encoding='utf8',header=None)
    uc1.loc[subg[0].tolist(),:].to_csv('map.sub.uc',sep='\t',index=0,na_rep='',header=0)

if __name__ == "__main__":
    opts=option()
    sub_uc(opts.i,opts.g)
    print('''生成：
          map.sub.uc
          ''')