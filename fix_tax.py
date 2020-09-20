# -*- coding: utf-8 -*-
"""
Created on Fri Dec 28 13:44:50 2018

@author: Administrator
"""

import argparse
import pandas as pd
import numpy as np
import re
from pandas.core.frame import DataFrame

def option():
    parser= argparse.ArgumentParser(description="""
    将含有taxonomy的OTU表格的门纲目科属种拆分                                
    """)
    parser.add_argument('-i', help="输入文件名", default='otu_select_fix_cluster_tax_assignments.txt')
    args= parser.parse_args()
    return args

def fix_tax(q): # q = opts.i
    taxl = ['d__','p__','c__','o__','f__','g__','s__']
    data1 = pd.read_csv(q,header=None,sep="\t")
    
    # Unknown 替换
    data1[1].replace("__Unknown_\w*","__Unknown",regex=True,inplace = True)
    
    # taxonomy 拆分
    data2 = [i.split(";") for i in data1.iloc[:,1]]
    bb = []
    for j in range(len(data2)):
        bb.append(data2[j]+[""]*(7-len(data2[j])))
    data3 = DataFrame(bb)
    data3.columns = taxl
    regex = re.compile(r'(;s__[a-zA-Z0-9_-]*)',flags=re.IGNORECASE)
    data3['tax'] = data1.iloc[:,1].replace(regex,"")
    data3['OTU_ID'] = data1.iloc[:,0]
    
    # 挑出需要修改的行 
    re1 = re.compile("uncultured|Incertae_Sedis$|Unknown|Subgroup_|norank|Family_",flags=re.IGNORECASE) # 门到属包含指定字符的行
    bl1 = data3['tax'].apply(lambda x:bool(re1.search(x)))
    re2 = re.compile(r'[dpcofgs]__',re.I)
    data4_1 = data3.loc[bl1,:].replace(re2,"").drop(['tax'], axis=1)
    #data4_1.to_csv("test.txt",header=0,index=0)
    for i in range(1,data4_1.shape[1]-2):# i = 4
        for k in range(data4_1.shape[0]):# k = 278
            if bool(re.search(data4_1.iloc[k,i]+"$",data4_1.iloc[k,i-1],re.I)):
                data4_1.iloc[k,i] = data4_1.iloc[k,i-1]
            elif bool(re.search(re1,data4_1.iloc[k,i])):
                data4_1.iloc[k,i] = data4_1.iloc[k,i-1] + "_" + data4_1.iloc[k,i]
                
    for i in range(0,data4_1.shape[1]-1):# i = 4
        for k in range(data4_1.shape[0]):# k = 278
            if data4_1.iloc[k,i] != "":
                data4_1.iloc[k,i] = data4_1.columns[i] + data4_1.iloc[k,i]
    data4_1['tax'] = data4_1.iloc[:,:-1].apply(lambda row: ';'.join(row.values.astype(str)), axis=1)
    # 合并
    data4_2 = data3.loc[~bl1,:]
    regex2 = re.compile(r'(_[0-9]*$)')
    data5 = pd.concat([data4_1[data4_2.columns.tolist()],data4_2]).replace(regex2,"")
    data5.index = data5["OTU_ID"]
    data6 = np.array(data5.loc[data1[0].tolist(),:]).tolist()
    data7 = DataFrame([[data6[pp][-1],";".join(data6[pp][0:7])] for pp in range(len(data6))])
    data7.replace(";+",";",regex=True,inplace = True)
    data7.replace(";$","",regex=True,inplace = True)
    data7.to_csv("fix_cluster_tax_assignments.txt",sep='\t',na_rep='',index=0,header=0)

if __name__ == "__main__":
    opts=option()
    fix_tax(opts.i)
