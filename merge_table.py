# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 16:14:04 2018

@author: Colin Liu
"""
import argparse
import pandas as pd
import copy
parser= argparse.ArgumentParser(description="""
说明：
 
""")
parser.add_argument('-i1', help="输入文件名", default='1genus.xls')
parser.add_argument('-i2', help="输入文件名", default='2genus.xls')
parser.add_argument('-d1', help="第一个样本名后缀", default='1')
parser.add_argument('-d2', help="第二个样本名后缀", default='2')
parser.add_argument('-s', help="挑选的列表", default='sample.txt')
parser.add_argument('-o', help="输出文件名", default='')
args= parser.parse_args()

if args.o == '':
    args.o = args.i1.split('.')[0] + "_" + args.i2.split('.')[0] + '.xls'

def read_tb(file,num):
    data = pd.read_table(file)
    data.columns = [x+'_'+num for x in data.columns]
    dd = ['Taxon'] + copy.deepcopy(data.columns).tolist()[1:]  
    data.columns = dd
    return data

def write_data(flist,outname):
    #写入表格
    out=open(outname,'w')
    for i in flist:
        out.write('\t'.join(i)+'\n')
    out.close()

data1 = read_tb(args.i1,args.d1)
data2 = read_tb(args.i2,args.d2)

if args.s != 'none':
    Data = pd.read_table(args.s)
    Data.iloc[:,0] = Data.iloc[:,0]+"_"+args.d1
    Data.iloc[:,1] = Data.iloc[:,1]+"_"+args.d2
    Data = [Data.iloc[j,:].tolist() for j in range(Data.shape[0])]
    DD = []
    for i in range(len(Data)):
        for j in Data[i]:
            DD.append([j,"S"+str(i)])
        
    write_data(DD,"map-group.txt")
    Data = sum(Data,[])

data3 = pd.merge(data1,data2,how='outer',on='Taxon')
#for i in data3.columns:
#    if data3[i].count != len(data3):#非NA的数目
#        loc = data3[i][data3[i].isnull().values==True].index.tolist()
#        print('列名："{}", 第{}行位置有缺失值'.format(i,loc))
data3 = data3.fillna(0)
data3.index = data3.iloc[:,0]
data3 = data3.drop(['Taxon'],axis=1).applymap(int)

if args.s != 'none':
    data3['Col_sum'] = data3.apply(pd.to_numeric).apply(lambda x: x.sum(), axis=1)
    data3 = data3.loc[data3['Col_sum']>0,:]
    data3 = data3[Data]
    
data3.loc[:,].to_csv(args.o,sep="\t")
