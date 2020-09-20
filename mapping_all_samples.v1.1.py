# -*- coding: utf-8 -*-
"""
Created on Mon Dec 25 10:11:22 2017

@author: Colin Liu
"""

import sys,os,time,shutil
import argparse
import pandas as pd
import re

parser= argparse.ArgumentParser(description="""
说明：这个脚本是。。。。。 
根据mapping不同的分组画含分组的图，并做高级分析
准备一个mapping.txt文件
#SampleID	ABC	AC
A1	A	A
A2	A	A
A3	A	A
B2	B	
B3	B	
B4	B	
C1	C	C
C2	C	C
C3	C	C

""")

parser.add_argument('-m', help="mapping文件", default='20191223.mapping.txt')
parser.add_argument('-i', help="输出文件名", default='')
parser.add_argument('-s', help="是否自动排序True", default=True)
args= parser.parse_args()

re_digits = re.compile(r'(\d+)')  
def embedded_numbers(s): # s= files[0] 
     pieces = re_digits.split(s)               # 切成数字与非数字  
     pieces[1::2] = map(int, pieces[1::2])     # 将数字部分转成整数  
     return pieces
    
def sort_strings_with_embedded_numbers(alist):  
     return sorted(alist, key=embedded_numbers)

def readmapping(filename):#filename = args.m
    '''读取table文件 .append or .extend '''  
    data = [] ##读取数据
    for ln in open(filename,'rt'):
        data.append(ln.strip('\n').split('\t'))
    group = [] ##得到所有的group名称
    for gp in data[0][1:]:
        if gp not in ['#SampleID', 'BarcodeSequence', 'LinkerPrimerSequence']:
            group.append(gp)
    #print len(data)    
    out={}
    for i in group:#i = group[0]
        gpvalu=[]
        for j in data[1:]:
            if j[data[0].index(i)] !='':
               gpvalu.append([j[0],j[data[0].index(i)]])
    return out

def dupList(ll):#ll = sort_strings_with_embedded_numbers(data_.iloc[:,1])
    list2=[]
    for i in ll:
        if i not in list2:
            list2.append(i)
    return list2

def readMap(filename,Sort):
    data = pd.read_csv(filename,header=0,sep='\t')
    for i in range(1,data.shape[1]):#i = 2
        data_ = data.iloc[:,[0,i]].dropna(axis=0)
        if Sort:
            LL = dupList(sort_strings_with_embedded_numbers(data_.iloc[:,1]))
            data_.iloc[:,1] = data_.iloc[:,1].astype("category").cat.set_categories(LL)
            
        data__ = data_.sort_values(by=data_.columns[1])
        data__.to_csv(data__.columns[1]+'group.txt',encoding='utf-8',header=0,index=False,sep="\t")
    return data

def writetable(lista,outname):
    out=open(outname,'w')    
    for i in lista:
        out.write('\t'.join(i)+'\n')
    out.close()
    
### mapping 列包含所有的样本
def plotall(i):
    os.system('''
## 放分析结果的文件######
gn=./按%s分组结果
\\rm -rf $gn
mkdir    $gn
cd $gn

\cp -r ../results/* ./
\cp ../%s map-group.txt

# cat map-group.txt|sort -k 2,2 >map-group1.txt
# mv map-group1.txt map-group.txt

cp ../group.analysis.v2.1.sh ./
sh group.analysis.v2.1.sh
cd ..
    '''%(i,i+'group.txt'))
    print(i)

#args.m = "mapping.txt"    
#mapping = readmapping(args.m)
   
#for i in mapping.keys():
#    writetable(mapping[i][1],i+'group.txt')
#    if mapping[i][0]:
#        print('all')
#        plotall(i)
#    else:
#        print('part')
#        plotall(i)

mapping = readMap(args.m,args.s)

for i in range(1,mapping.shape[1]):
    plotall(mapping.columns[i])
              
print('HAHA ok....')

