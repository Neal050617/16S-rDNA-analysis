# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 08:55:05 2018

@author: Colin Liu
"""

import os
import copy
import argparse
import re
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
    parser.add_argument('-f', help="是否要对物种表格进行矫正", default='T')
    args= parser.parse_args()
    return args

def tax_split(tax):
    data=tax.split(";")
    return data

# 取出键 
def handle_str(tax,m):# tax,m = ts[j],"_"
    if m == "_":
        if "__" in tax:
            taxx = tax[tax.index("__")+2:]
        if "(" and ")" in tax:
            taxx = tax[:tax.index('(')]
    else:
        taxx = copy.deepcopy(tax)
    return taxx

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
def split_tax(data,m):# data,m = copy.deepcopy(data),""
    datas = copy.deepcopy(data)
    for i in range(1,len(datas)):# i=1
        da=datas[i]
        ts=tax_split(da[taxwz])
        for j in range(len(ts)):# j=0            
            if "__" in ts[j]:
                k1=ts[j][:ts[j].index("__")+2]
                if k1 in adict:
                    da[datas[0].index(adict[k1])]=handle_str(ts[j],m)
                else:
                    da[-1] = handle_str(ts[j],m)
            elif 'unclassified' not in ts[j]:
                da[-1]==handle_str(ts[j],m)
    return datas

def otu_tax(tt,st): # tt,st = copy.deepcopy(data_),"genus"
    tt1 = ListoDF(tt)
    tt1.insert(0,r"OTU_tax",tt1.iloc[:,0]+r" ("+tt1[st]+")")
    tt1.drop("OTU_ID",axis=1,inplace=True)
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
            kd3 = kd1.loc[:,["OTU_ID"] + kd4["#SampleID"][kd4.group.isin([i])].tolist() + liaa]
            coln = [i for i,x in enumerate(kd3.columns) if kd3.columns[i]=='kingdom']
            kd3['Col_sum'] = kd3.iloc[:,1:coln[0]].apply(pd.to_numeric).apply(lambda x: x.sum(), axis=1)
            kd3[kd3['Col_sum'] > 0].loc[:,["Col_sum"]+liaa].to_csv('Krona/krona.'+str(i)+".xls",sep='\t',na_rep='',index=0,header=0)

def unroot_tax(rr,tt,gg): # rr,tt,gg = copy.deepcopy(normlize_tb(data_)),opts.t,opts.r
    rr1 = ListoDF(rr)
    rtg = rr1.loc[rr1.iloc[:,1:taxwz].apply(lambda x:any(x>float(tt)),axis = 1),["OTU_ID"]+liaa]
    rtg.to_csv("unrooted.otu_tree_anno."+tt+".xls",sep='\t',na_rep='',index=0)
    # 输出fasta
    if gg != "none":
        write_fasta(read_fa(gg,np.array(rtg.iloc[:,0]).tolist()),"map."+gg)  

def core_micro(cc,bb): # cc,bb = copy.deepcopy(data_),opts.m
    cb1 = ListoDF(cc)
    cbo = []
    cb1["core"] = cb1.iloc[:,1:taxwz].apply(pd.to_numeric).apply(lambda x:len(x[x>0])/len(x),axis=1)
    for i in [0.5,0.6,0.7,0.8,0.9,1]:#i = 0.5
        cb1.loc[cb1["core"] >= i,["OTU_ID","taxonomy","core"]].to_csv("Core_Microbime/core.microbiome."+"gt"+str(i)+".xls",sep='\t',na_rep='',index=0)
        cbo.append(cb1.loc[cb1["core"] >= i,["core"]].shape[0])
    
    cbd = pd.DataFrame({'rate' : np.array([0.5,0.6,0.7,0.8,0.9,1]),'All' : np.array(cbo)})
    cbd.index = cbd['rate']
    cbd.drop('rate',axis=1,inplace=True)
    
    if bb != "none":
        cb2 = pd.read_csv(bb,sep  = "\t",header=0) # 分组文件
        cbl = cb2.group.unique().tolist()
        for j in cbl:#j = cbl[0]
            cb3 = cb1.loc[:,["OTU_ID"] + cb2["#SampleID"][cb2.group.isin([j])].tolist() + ["taxonomy"]]
            cb3["core"] = cb3.iloc[:,1:cb3.shape[1]-1].apply(pd.to_numeric).apply(lambda x:len(x[x>0])/len(x),axis=1)
            cb4 = []
            for k in [0.5,0.6,0.7,0.8,0.9,1]:#k = 0.5
                cb3.loc[cb3["core"] >= k,["OTU_ID","taxonomy","core"]].to_csv("Core_Microbime/core.microbiome."+j+".gt"+str(k)+".xls",sep='\t',na_rep='',index=0)
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

def fix_tax(q):# q = copy.deepcopy(table)
    taxl = ['d__','p__','c__','o__','f__','g__','s__']
    Df1 = DataFrame(q[1:])
    Df1.columns = q[0]
    taxwz=Df1.columns.tolist().index('taxonomy') # taxonomy的位置
    
    # Unknown 替换
    Df1.iloc[:,taxwz].replace("__Unknown_\w*","__Unknown",regex=True,inplace = True)
    # taxonomy 拆分
    Df2 = [i.split(";") for i in Df1.iloc[:,taxwz]]
    bb = []
    for j in range(len(Df2)):
        bb.append(Df2[j]+[""]*(7-len(Df2[j])))
    Df3 = DataFrame(bb)
    Df3.columns = taxl
    regex = re.compile(r'(;s__[a-zA-Z0-9_-]*)',flags=re.IGNORECASE)
    Df3['tax'] = Df1.iloc[:,taxwz].replace(regex,"")
    Df3['OTU_ID'] = Df1.iloc[:,0]
    
    # 挑出需要修改的行 
    re1 = re.compile("uncultured|Incertae_Sedis|Unknown|Subgroup_|norank|Family_",flags=re.IGNORECASE) # 门到属包含指定字符的行
    bl1 = Df3['tax'].apply(lambda x:bool(re1.search(x)))
    re2 = re.compile(r'[dpcofgs]__',re.I)
    Df4_1 = Df3.loc[bl1,:].replace(re2,"").drop(['tax'], axis=1)
    #Df4_1.to_csv("Df4_1.txt",header=0,index=0)
    for i in range(1,Df4_1.shape[1]-2):# i = 4
        for k in range(Df4_1.shape[0]):# k = 278
            if bool(re.search(Df4_1.iloc[k,i]+"$",Df4_1.iloc[k,i-1],re.I)):
                Df4_1.iloc[k,i] = Df4_1.iloc[k,i-1]
            elif bool(re.search(re1,Df4_1.iloc[k,i])):
                Df4_1.iloc[k,i] = Df4_1.iloc[k,i-1] + "_" + Df4_1.iloc[k,i]
    
    for i in range(0,Df4_1.shape[1]-1):# i = 4
        for k in range(Df4_1.shape[0]):# k = 278
            if Df4_1.iloc[k,i] != "":
                Df4_1.iloc[k,i] = Df4_1.columns[i] + Df4_1.iloc[k,i]
    Df4_1['tax'] = Df4_1.iloc[:,:-1].apply(lambda row: ';'.join(row.values.astype(str)), axis=1)
    # 合并
    Df4_2 = Df3.loc[~bl1,:].replace(re2,"").drop(['tax'], axis=1)
    #Df4_2.to_csv("Df4_2.txt",header=0,index=0)
    for i in range(1,Df4_2.shape[1]-2):# i = 4
        for k in range(Df4_2.shape[0]):# k = 278
            if (Df4_2.iloc[k,i] == "") & ("_unclassified" not in Df4_2.iloc[k,i-1]):
                Df4_2.iloc[k,i] = Df4_2.iloc[k,i-1] + "_unclassified"
            elif (Df4_2.iloc[k,i] == "") & ("_unclassified" in Df4_2.iloc[k,i-1]):
                Df4_2.iloc[k,i] = Df4_2.iloc[k,i-1]
    
    for i in range(0,Df4_2.shape[1]-1):# i = 4
        for k in range(Df4_2.shape[0]):# k = 278
            if Df4_2.iloc[k,i] != "":
                Df4_2.iloc[k,i] = Df4_2.columns[i] + Df4_2.iloc[k,i]
    Df4_2['tax'] = Df4_2.iloc[:,:-1].apply(lambda row: ';'.join(row.values.astype(str)), axis=1)    
    
    regex2 = re.compile(r'(_[0-9]*$)')
    Df5 = pd.concat([Df4_1[Df4_2.columns.tolist()],Df4_2]).replace(regex2,"")
    Df5.index = Df5["OTU_ID"]
    Df6 = np.array(Df5.loc[Df1.iloc[:,0].tolist(),:]).tolist()
    Df7 = DataFrame([[Df6[pp][-1],";".join(Df6[pp][0:7])] for pp in range(len(Df6))])
    Df7.iloc[:,1] = DataFrame([[Df7.iloc[pp,1],";".join(Df3.loc[pp,'s__'])] for pp in range(len(Df6))]).iloc[:,0]
    Df7.replace(";+",";",regex=True,inplace = True)
    Df7.replace(";$","",regex=True,inplace = True)
    Df1.iloc[:,taxwz] = Df7.iloc[:,1]
    
    DF = copy.deepcopy(np.array(Df1).tolist())
    DF.insert(0,q[0])
    
    return [i[:(taxwz+1)] for i in DF]

if __name__=='__main__':
    opts=option()
    # 分类级别 
    adict={'d__':'kingdom','p__':'phylum','c__':'class','o__':'order','f__':'family','g__':'genus','s__':'species'}
    liaa=['kingdom','phylum','class','order','family','genus','species','unclassified']
    
    # 初始OTU表格 
    table = []
    for i in open(opts.i,'rt'):
        table.append(i.strip('\n').split('\t'))
        
    # 是否要对物种信息进行校正
    if opts.f == 'T':
        table = fix_tax(table)
    
    # 如果存在group文件，则筛选OTU表格
    if opts.m != "none":
        table = map_split(opts.m,table)
    ListoDF(table).to_csv("map."+opts.i,sep='\t',na_rep='',index=0)
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
    normlize_tb(datax).to_csv("percent.split_tax_"+opts.i,sep='\t',na_rep='',index=0)
      
    print("this is after if __name__:%s"%__name__)
   