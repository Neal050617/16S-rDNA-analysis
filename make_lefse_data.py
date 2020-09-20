# -*- coding: utf-8 -*-
"""
Created on Sat Aug 08 14:12:29 2015

@author: BioLinker05
"""

import getopt,sys,os

opts, args = getopt.getopt(sys.argv[1:], "hi:o:g:", ["help", "imput=","output="])

def usage():
        print '''
        准备画lefse的文件，
        -i   输入文件    rank文件，门纲目科属种分开为列的文件
        -o   输出文件
        -g   group文件
        
        修改记录：
        160629  替换多个|||为一个
        161112  词前面加上前缀
'''

input='rarefac.otu_taxa_split_index_table.xls'
outname='lefse.data.xls'
groupfile='map-group.txt'

for o, a in opts:
        if o in ('-i','--imput'):
                imput = a
        elif o in ('-o','--output'):
                outname = a
        elif o in ('-g','--group'):
                groupfile = a
                
        elif o in ('-h'):
                usage()
                sys.exit()
#


file = open(input,'r')

data1 = []
for ln in file:
    data1.append(ln.strip().split('\t'))
    
data=[]
for i in range(0,len(data1)):     # 列表后补齐缺省值
    data.append(data1[i]+[""]*(len(data1[0])-len(data1[i])))

    
def unique(old_list):
    newList = []
    for x in old_list:
        if x not in newList :
            newList.append(x)
    return newList    

#taxonomy 的位置
OTUwz=data[0].index('taxonomy')

def sumtax(list):  #计算每个具体分类级别的丰度。list就是得到的taxnum
    newlist=[]
    for i in range(1,OTUwz):
        sum=0
        for j in list:
            sum=sum+float(data[j][i])
        newlist.append(sum)
    return newlist

def writedata(data,name):
    out=open(name,'w')
    for i in data:
        out.write('\t'.join(i)+'\n')
    out.close()
    
def readgroup(groupfile):
   #读取group文件为字典   
   group = {}
   with open(groupfile,'r') as df:
        for kv in [d.strip().split('\t') for d in df]:
            group[kv[0]] = kv[1]
   return group
   
def pickdataby_group(data,samplelist):
    #挑选group中含有的列
    #out=[[data[0][0]]+samplelist]
    out=[]
    for j in data[:]:#j=data[0]
        out1=[j[0]]
        classa=['class']
        for i in samplelist:#i='ko142'
            if i in data[0]:
                classa.append(group[i])
                out1.append(j[data[0].index(i)])
            else:
                print('sample '+i+ ' not in data')
        out.append(out1)
    out3=[classa]+out
    return out3        


adict={'kingdom':'d__','phylum':'p__','class':'c__','order':'o__','family':'f__','genus':'g__','species':'s__'}

lineage=["kingdom","phylum","class","order","family","genus","species"]

#lineage=["superkingdom"]
lineage2=lineage
for i in lineage2:
    if i not in data[0]:
       lineage.remove(i)
 
out=open(outname,'w')

#写第一行列标
#out.write('"number"'+'\t'+'"classification"'+'\t'+'"tax"')
out.write('tax')
for i in data[0][1:OTUwz]:
    out.write('\t'+i)
out.write('\n')

#写后面的
for i in lineage:
    #i="phylum"
    id=data[0].index(i)

    tax=[]
    for j in range(len(data)):  #phylum下的所有门
        tax.append(data[j][id])
    
    uniquetax=unique(tax)    #所有的门去除重复
    if '' in uniquetax:
       uniquetax.remove('')
    
    for a in range(1,len(uniquetax)):#a=1
        taxnum=[h for h,v in enumerate(tax) if v==uniquetax[a]]#门所在的行
        #list(enumerate(tax))[0]
        #tax.count(uniquetax[0])
        tax1=[]
        for li in range(lineage.index(i)+1):#li=0
            hang=taxnum[0] #这是在data里的位置
            shu=data[0].index(lineage[li])#22 #该分类级别所在的列
            if data[hang][shu] !='':
                #tax1.append(adict[data[0][shu]]+data[hang][shu])
                tax1.append(data[hang][shu])

                #print adict[data[0][shu]]
                #print adict[data[0][shu]]+data[hang][shu]
                #print hang,shu
                #print tax1
            
        out.write('|'.join(tax1))  #'|'竖线分开每个词          
        for item in sumtax(taxnum):
            out.write("\t"+str(item))   #//list中一项占一行
        out.write('\n')
out.close()          

   
if groupfile!='':
   #读取文件     
    data22 = []
    for ln in open(outname,'rt'):
        data22.append(ln.strip().split('\t'))
        
    group=readgroup(groupfile)
    del group['#SampleID']
    outdata=pickdataby_group(data22,group.keys())
    #print outdata[0:2]
    writedata(outdata,outname)
    

print "output:%s"%outname



