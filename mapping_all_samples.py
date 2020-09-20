# -*- coding: utf-8 -*-
"""
Created on Mon Dec 25 10:11:22 2017

@author: Colin Liu
"""

import sys,os,time,shutil
import argparse

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

parser.add_argument('-m', help="mapping文件", default='mapping.txt')
parser.add_argument('-i', help="输出文件名", default='')
args= parser.parse_args()

def readmapping(filename):    
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
    for i in group:
        gpvalu=[]
        for j in data[1:]:
            if j[data[0].index(i)] !='':
               gpvalu.append([j[0],j[data[0].index(i)]])
        #print len(data)-1==len(gpvalu)
        out[i]=[len(data)-1==len(gpvalu),gpvalu]
    return out

def writetable(lista,outname):
    out=open(outname,'w')    
    for i in lista:
        out.write('\t'.join(i)+'\n')
    out.close()
    
#args.m = "mapping.txt"    
mapping = readmapping(args.m)

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

cat map-group.txt|sort -k 2,2 >map-group1.txt
mv map-group1.txt map-group.txt

cp /work/users/huifangan/step3_shell/less_than_100/step3.sh ./
sh step3.sh
cd ..
    '''%(i,i+'group.txt'))
    print(i)

   
for i in mapping.keys():
   writetable(mapping[i][1],i+'group.txt')
   if mapping[i][0]:
       print('all')
       plotall(i)
   else:
       print('part')
       plotall(i)
              
print('HAHA ok....')






