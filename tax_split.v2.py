# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 08:55:05 2018

@author: Colin Liu
"""

import os
import copy
import argparse
import pandas as pd
import numpy as np
from pandas.core.frame import DataFrame
import matplotlib.pyplot as plt
import matplotlib as mpl

def option():
    parser= argparse.ArgumentParser(description="""
    将含有taxonomy的OTU表格的门纲目科属种拆分
    20190114-unroot tree annotation
    20190116-group选项、无根树丰度筛选
    20190325-核心微生物列出
    """)
    parser.add_argument('-i', help="输入文件名", default='rarefac.otu_taxa_table.xls')
    parser.add_argument('-m', help="分组", default='map-group.txt')
    parser.add_argument('-x', help="otu_tax", default='genus')
    parser.add_argument('-r', help="OTUrep,只包含OTU序号和序列:otu_reps.raw.fasta", default='none')
    parser.add_argument('-t', help="无根树筛选OTU丰度", default='0.001')
    args= parser.parse_args()
    return args

def tax_split(tax):
    data=tax.split(";")
    return data

# 取出键 
def handle_str(tax,m):
    if m == "_":
        if "__" in tax:
            tax=tax[tax.index("__")+2:]
        if "(" and ")" in tax:
            tax=tax[:tax.index('(')]
    return tax

# 列表转为字典
def ListoDF(data):
    if isinstance(data, list):    
        Df = DataFrame(data)  # 转为数据框
        Df.columns = Df.iloc[0,:]  # 修改列名
        Df.drop(0,axis=0,inplace=True)   # 删除第一行 
    else:
        Df = data          
    return Df
 
# 表格均一化
def normlize_tb(data): # data = copy.deepcopy(data_) # data = ms4.iloc[:,range(0,ms4.shape[1]-1)]
    # 判断输入文件类型 
    if type(data) == list:
       Df = ListoDF(data) 
    else:
       Df = data
     # 判断是否有物种信息 
    if 'taxonomy' in Df.columns:
        cols = [i for i,x in enumerate(Df.columns) if Df.columns[i]=='taxonomy']  # taxonomy的位置
    else:
        cols = [Df.shape[1]+1]
    Df_1 = Df.iloc[:,1:cols[0]].apply(pd.to_numeric)  # 数据框转数字
    # 均一化为1
    for i in range(Df_1.shape[1]):#i=0
        Df_1.iloc[:,i] = Df_1.iloc[:,i].apply(lambda x: x/Df_1.iloc[:,i].sum())
    # 输出文件  
    Df_2 = pd.concat([Df.iloc[:,0],Df_1,Df.iloc[:,cols[0]:]],axis=1)
    return Df_2

# 分割成门纲目科属种 
def taxon_split(data,name):#name=".percent.xls" 
    datan = ListoDF(data)
    cols=[i for i,x in enumerate(datan.columns) if datan.columns[i]=='taxonomy']  # taxonomy的位置    
    for p in liaa[1:-1]:#p='phylum'
        if any(datan[p] != ''):
            datan_1 = datan[~datan[p].isin([""])] # 去除包含空字符串的行
            datan_2 = datan_1.iloc[:,1:cols[0]].apply(pd.to_numeric).groupby(datan_1[p]).sum()
            datan_2.to_csv('Community/'+p+name,sep='\t',na_rep='')
            datan_2.index.name = "Taxon"

# taxonomy 分割 
def split_tax(datas,m):
    data = copy.deepcopy(datas)
    for i in range(1,len(data)):#i=1
        da=data[i]
        ts=tax_split(da[taxwz])
        for j in range(len(ts)):#j=0            
            if "__" in ts[j]:
                k1=ts[j][:ts[j].index("__")+2]
                if k1 in adict:
                    da[data[0].index(adict[k1])]=handle_str(ts[j],m)
                else:
                    da[-1] = handle_str(ts[j],m)
            elif 'unclassified' not in ts[j]:
                da[-1]==handle_str(ts[j],m)
    return data

def otu_tax(tt,st): # tt,st = copy.deepcopy(data_),"genus"
    tt1 = ListoDF(tt)
    tt1.insert(0,r"#OTU_tax",tt1.iloc[:,0]+r" "+tt1[st])
    tt1.drop(r"#OTU ID",axis=1,inplace=True)
    num = tt1.columns.tolist().index("taxonomy")   
    tt1.iloc[:,0:num].to_csv("otu."+st+".xls",sep='\t',na_rep='',index=0)

def write_table(data,name):
    output=name
    out=open(output,'w')
    for i in range(len(data)):
        for j in range(len(data[i])):
            out.write(data[i][j]+'\t')
        out.write('\n')
    out.close()  

def read_fa(faname,falist): # faname,falist = opts.f,d2["seq_name"].tolist()
    '''读取fasta为字典'''
    out={}
    name = seq = ''
    for ln in open(faname,'rt'):
        if ln[0]=='>':
           name=ln.strip()[1:]
           seq=''
        else:
           seq+=ln.strip()
        out[name]=seq
    # 取出
    data=[]
    for i in falist:
        data.append(['>'+i,out[i]])
    return data

def write_fasta(flist,outname):
    #写入表格
    out=open(outname,'w')
    for i in flist:
        out.write('\n'.join(i)+'\n')
    out.close()

def map_split(g1,g2): # g1,g2=opts.m,copy.deepcopy(table)
   ms1 = pd.read_csv(g1,sep  = "\t")
   ms2 = ListoDF(g2)
   ms3 = ms2.loc[:,[ms2.columns[0]]+ms1.iloc[:,0].tolist()+['taxonomy']]
   # 列求和
   Row_sum = ms3.iloc[:,range(1,ms3.shape[1]-1)].astype("int").apply(lambda x: x.sum(), axis=1)
   ms4 =  ms3.loc[Row_sum > 0,:]
   ms4.iloc[:,range(0,ms4.shape[1]-1)].to_csv("map.otu_table.xls",sep='\t',na_rep='',index=0)
   normlize_tb(ms4.iloc[:,range(0,ms4.shape[1]-1)]).to_csv("map.otu_table.percent.xls",sep='\t',na_rep='',index=0)

   return [ms4.columns.tolist()]+np.array(ms4).tolist()
     
def lefse_split(Data,gg):#Data,gg = copy.deepcopy(datax),opts.m
    DF = ListoDF(Data)
    DF.drop(['taxonomy'],axis=1,inplace=True)
#    for i in datal:
#        i.remove(i[taxwz])
#    DF = ListoDF(datal)
    DD = DataFrame()
    for k in liaa[:-1]:# k='phylum'
        tax_ = liaa.index(k)
        if any(DF[k] != ''):
            DF_1 = DF[~DF[k].isin([""])].iloc[:,list(range(taxwz+tax_+1))] # 去除包含空字符串的行 
            DF_1['taxon'] = DF_1.iloc[:,range(taxwz,(1+taxwz+tax_))].apply(lambda x:"|".join(x),axis = 1)
            DF_1.drop(liaa[:tax_+1],axis=1,inplace=True)
            DF_2 = DF_1.iloc[:,1:taxwz].apply(pd.to_numeric).groupby(DF_1['taxon']).sum()
            DD = pd.concat([DD,DF_2])
            
    if gg != "none":
        # read group 
        DD1 = pd.read_csv(gg,sep  = "\t",header=0).T
        DD1.columns = DD1.iloc[0,:]
        DD1.drop([r'#SampleID'],axis=0,inplace=True)
        # select data
        DD2 = pd.concat([DD1,DD.loc[:,DD1.columns]],axis=0)
        DD2.rename(index={'group':'class'},inplace=True)        
        DD2.to_csv("lefse.data.xls",sep='\t',na_rep='',float_format='%.1f',header=0)
    else:
        DD.to_csv("lefse.data.xls",sep='\t',na_rep='',float_format='%.1f')

def krona_data(kk,gg): # kk,gg = copy.deepcopy(data_),opts.m
    kd1 = ListoDF(kk)
    kd1['Col_sum'] = kd1.iloc[:,1:taxwz].apply(pd.to_numeric).apply(lambda x: x.sum(), axis=1)
    dd = list(range(taxwz+1,kd1.shape[1]-1))
    dd.insert(0,kd1.shape[1]-1)
    kd2 = kd1.iloc[:,dd]
    kd2.to_csv('Krona/krona.data.xls',sep='\t',na_rep='',index=0,header=0)
    
    if gg != "none":
        kd4 = pd.read_csv(gg,sep  = "\t",header=0) # 分组文件
        kdl = kd4.group.unique().tolist()
        for i in kdl: # i = kdl[0]
            kd3 = kd1.loc[:,[r"#OTU ID"] + kd4["#SampleID"][kd4.group.isin([i])].tolist() + liaa]
            coln = [i for i,x in enumerate(kd3.columns) if kd3.columns[i]=='kingdom']
            kd3['Col_sum'] = kd3.iloc[:,1:coln[0]].apply(pd.to_numeric).apply(lambda x: x.sum(), axis=1)
            kd3[kd3['Col_sum'] > 0].loc[:,["Col_sum"]+liaa].to_csv('Krona/krona.'+str(i)+".xls",sep='\t',na_rep='',index=0,header=0)

def unroot_tax(rr,tt,gg): # rr,tt,gg = copy.deepcopy(normlize_tb(data_)),opts.t,opts.r
    rr1 = ListoDF(rr)
    rtg = rr1.loc[rr1.iloc[:,1:taxwz].apply(lambda x:any(x>float(tt)),axis = 1),[r'#OTU ID']+liaa]
    rtg.to_csv("unrooted.otu_tree_anno."+tt+".xls",sep='\t',na_rep='',index=0)
    # 输出fasta
    if gg != "none":
        write_fasta(read_fa(gg,np.array(rtg.iloc[:,0]).tolist()),"map."+gg)  

def core_micro(cc,bb): # cc,bb = copy.deepcopy(data_),opts.m
    cb1 = ListoDF(cc)
    cbo = []
    cb1["core"] = cb1.iloc[:,1:taxwz].apply(pd.to_numeric).apply(lambda x:len(x[x>0])/len(x),axis=1)
    for i in [0.5,0.6,0.7,0.8,0.9,1]:#i = 0.5
        cb1.loc[cb1["core"] >= i,[r"#OTU ID","taxonomy","core"]].to_csv("Core_Microbime/core.microbiome."+"gt"+str(i)+".xls",sep='\t',na_rep='',index=0)
        cbo.append(cb1.loc[cb1["core"] >= i,["core"]].shape[0])
    
    cbd = pd.DataFrame({'rate' : np.array([0.5,0.6,0.7,0.8,0.9,1]),'All' : np.array(cbo)})
    cbd.index = cbd['rate']
    cbd.drop('rate',axis=1,inplace=True)
    
    if bb != "none":
        cb2 = pd.read_csv(bb,sep  = "\t",header=0) # 分组文件
        cbl = cb2.group.unique().tolist()
        for j in cbl:#j = cbl[0]
            cb3 = cb1.loc[:,[r"#OTU ID"] + cb2["#SampleID"][cb2.group.isin([j])].tolist() + ["taxonomy"]]
            cb3["core"] = cb3.iloc[:,1:cb3.shape[1]-1].apply(pd.to_numeric).apply(lambda x:len(x[x>0])/len(x),axis=1)
            cb4 = []
            for k in [0.5,0.6,0.7,0.8,0.9,1]:#k = 0.5
                cb3.loc[cb3["core"] >= k,[r"#OTU ID","taxonomy","core"]].to_csv("Core_Microbime/core.microbiome."+j+".gt"+str(k)+".xls",sep='\t',na_rep='',index=0)
                cb4.append(cb3.loc[cb3["core"] >= k,["core"]].shape[0])
            cbd[j] = cb4
 
    cbd.to_csv("Core_Microbime/core.microbiome.output.xls",sep='\t',na_rep='',index=0)
    
    # matplotlib 
    colors = [r"#B0C4DE",r"#CD3333",r"#483D8B",r"#458B00",r"#EEC900",r"#BFBC3B",r"#6E8B3D",r"#00688B",r"#C10077",r"#CAAA76",r"#EEEE00",r"#458B00",r"#8B4513",r"#008B8B",r"#6E8B3D",r"#8B7D6B",r"#7FFF00",r"#CDBA96",r"#ADFF2F"]
    plt.style.use('ggplot') #plt.style.available
    cbdp = cbd.plot(kind='barh',color=[colors[i] for i,j in enumerate(cbd.columns.tolist())]) #
    plt.grid(axis='x')
    plt.xlabel('rate',size=12)
    plt.ylabel('count',size=12)
    plt.xlim(0,max(cbd.iloc[0,:])+10)
    plt.title('Core_Microbime',size=20)
    plt.legend(fancybox=True,framealpha=1,shadow=True,borderpad=1)
    
    for rect in cbdp.patches:# rect = cbdp.patches[2]
        plt.text(rect.get_width()+3, rect.get_y() + rect.get_height() / 2, str(rect.get_width()),ha='center',va='center',fontsize=10)
        #plt.text(rect.get_x() + rect.get_width() / 2, rect.get_height() * 1.005, str(rect.get_height()), 
        #         ha='center', va='center', fontsize=10) # use bar
    plt.savefig('Core_Microbime/core.microbiome.barplot.png')

if __name__=='__main__':
    opts=option()
    # 分类级别 
    adict={'d__':'kingdom','p__':'phylum','c__':'class','o__':'order','f__':'family','g__':'genus','s__':'species'}
    liaa=['kingdom','phylum','class','order','family','genus','species','unclassified']
    
    # 初始OTU表格 
    table = []
    for i in open(opts.i,'rt'):
        table.append(i.strip('\n').split('\t'))
    # 如果存在group文件，则筛选OTU表格
    if opts.m != "none":
        table = map_split(opts.m,table)
    table[0].extend(liaa)
    # 补齐物种分类 
    data=[]
    for i in range(0,len(table)):     # 列表后补齐缺省值
        data.append(table[i]+[""]*(len(table[0])-len(table[i])))
    taxwz=data[0].index('taxonomy') # taxonomy的位置
    
    # 添加分类级别后的物种表格
    data_ = split_tax(data,"_")
    datax = split_tax(data,"")
    
    os.system('''mkdir Krona''')    
    # 分析用表格
    lefse_split(datax,opts.m)
    krona_data(data_,opts.m)
    otu_tax(data_,opts.x)
    if opts.r != "none":
        unroot_tax(normlize_tb(data_),opts.t,opts.r)
    # 门纲目科属种表格
    os.system('''mkdir Community''')
    taxon_split(data_,".xls")
    taxon_split(normlize_tb(data_),".percent.xls")
    # 核心微生物
    #os.system('''mkdir Core_Microbime''')
    #core_micro(data_,opts.m)
    
    write_table(data_,"split_"+opts.i)
    write_table(datax,"split_tax"+opts.i)
    normlize_tb(data_).to_csv("percent.split_"+opts.i,sep='\t',na_rep='',index=0)
    normlize_tb(datax).to_csv("percent.split_tax"+opts.i,sep='\t',na_rep='',index=0)
      
    print("this is after if __name__:%s"%__name__)
   