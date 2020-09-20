# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 14:38:14 2018

@author: Colin Liu
"""

import re
import argparse

parser= argparse.ArgumentParser(description="""
说明：拆分RDP数据库，与silva格式一致 
""")

parser.add_argument('-i', help="输入文件名", default='1.txt')
parser.add_argument('-f', help="输出文件名", default='')
parser.add_argument('-t', help="输出文件名", default='')
args= parser.parse_args()

imput=args.i

if args.f == '' or args.t == '':
    output1='fasta_out_'+imput
    output2='tax_out_.'+imput
else:
    output1=args.f
    output2=args.t


#读取fasta_tax文件
data = []
for ln in open(imput,'rt'):
    data.extend(ln.strip().split('\n'))

#序列名所在的行
hs=[] 
for i in range(len(data)):
    if '>' in data[i]:
        hs.append(i)
hs.append(len(data)+1)

def Re(d):
    a=re.compile(r'#|/')
    b=a.sub('_',d)
    b=b.replace(';','')
    b=b.replace('*','_')            
    return b

def Split_tax(tax):#拆分序列名tax=['Bacteria', 'domain', 'unclassified_Bacteria', '']
    List=['domain','phylum','class','order','family','genus']
    D2={}
    if tax[-1]!='':
        num=range(len(tax))
    else:
        num=range(len(tax)-3)
    
    for j in num:
        if j%2==0 :
            D2[tax[j+1]]=tax[j]
            
    tmp = [val for val in D2.keys() if val in List]
    
    Tax=[]
    for i in tmp:
        Tax.append(D2[i])
    Tax1="; ".join(Tax)
    return Tax1

def Split_fa(data,hs):
    #得到序列名和数据库里的编号
    fasta=[]
    fa_name=[]
    for i in range(len(hs)-1):
        fasta.append([data[hs[i]].split(' ')[0]])
        fasta.append([''.join(data[hs[i]+1:hs[i+1]]).upper()])
        #fa_name.append([data[hs[i]]])
    for i in range(len(hs)-1):    
        fa_name.append(data[hs[i]].replace("\"",""))
    ft=[]
    for i in range(len(fa_name)):
        data1=fa_name[i].split(' ')[0][1:]
        data2=Re("_".join(fa_name[i].split('\t')[0].split(' ')[1:]))
        data3=Split_tax(fa_name[i].split('Lineage=Root;rootrank;')[1].split(';'))
        ft.append([data1,data3,data2])
    return fasta,ft


def pretax(data):
    cas=['d__','p__','c__','o__','f__','g__','s__']
    data2=[]
    for i in range(len(data)):
        for j in range(7):
            data2.extend([cas[j]+data[i][j]])
    return data2
 
def write_data(flist,outname):
    #写入表格
    out=open(outname,'w')
    for i in flist:
        out.write('\t'.join(i)+'\n')
    out.close()

DA1=Split_fa(data,hs)
DA11=DA1[0]
DA12=[[DA1[1][i][0]] for i in range(len(DA1[1]))]

DA2=[]
for i in range(len(DA1[1])):
    da2=[]
    da2=DA1[1][i][1].split("; ")
    da2=da2+["unclassified_"+da2[-1]]*(6-len(da2))+[DA1[1][i][2]]
    DA2.append(da2)

DA3=pretax(DA2)
DA3=[DA3[i:i+7] for i in range(0,len(DA3),7)]

DA4=[]
for i in range(len(DA3)):
    DA4.append([DA12[i][0]+"\t"+"; ".join(DA3[i])])

write_data(DA11,output1)
write_data(DA4,output2)

print('成功')

 