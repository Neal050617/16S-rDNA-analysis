# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 11:57:41 2019

@author: Administrator
"""

import os
import argparse
import pandas as pd
import numpy as np
from pandas.core.frame import DataFrame
from scipy.stats import t

def option():
    parser= argparse.ArgumentParser(description="""
    整合alpha多样性文件
    需要补充一个置信区间                           
    """)
    parser.add_argument('-i', help="输入类别", default='rarefaction')
    parser.add_argument('-m', help="是否存在分组", default="none")
    parser.add_argument('-c', help="是否指定颜色color.txt", default="none")
    parser.add_argument('-b', help="均值中位数误差棒/False", default="False")
    args= parser.parse_args()
    return args

##### 提取mode #######
def get_name(name):
    a=name.split('.')
    return a[1]

def flat(nums): # 嵌套列表展开
    res = []
    for i in nums:
        if isinstance(i, list):
            res.extend(flat(i))
        else:
            res.append(i)
    return res

def dirs2list(l1): # l1 = opts.i    
    #列出目录下的所有文件和目录,提取出需要的文件名
    path = os.getcwd()
    List = os.listdir(path)
    lista=[]
    for i in List:
    #if '.rarefaction' in i:
        if '.' + l1  in i:
            lista.append(i)
    return lista

## 读取文件group  
def read_gp(name): # name = "map-group.txt"
    gp1 = pd.read_csv(name,sep  = "\t",header=0) # 分组文件
    allkey = gp1.group.unique().tolist()
    alldict = {gp1.iloc[i,0] : gp1.iloc[i,1] for i in range(gp1.shape[0])}
    return [allkey,alldict,gp1]

#提取文件的名字列，0.03列
def draw_col(name,lie1): # name,lie1 = "otus.fvc7.rarefaction","0.97"
    #print name
    data = pd.read_csv(name, sep  = "\t", header=0, index_col=False) # 分组文件
    lie = data.loc[:,lie1].tolist()
    return lie

#提取需要的数据tiqu_data(lista,out_xls,'numsampled',label)               
def draw_data(namelist,outname,lie1,lie2): # namelist,outname,lie1,lie2 = lista,out_xls,'numsampled',"0.97"
     #所有的文件名、输出文件名、第一行、后面行
    out1=open(outname,'w')
    data = []
    for i in namelist: # i= namelist[0]
        data.append([get_name(i)] + draw_col(i,lie2))    
    
    nm1 = [len(i) for i in data].index(max([len(i) for i in data]))
    data.insert(0,["sampleID"] + draw_col(namelist[nm1],lie1))
    data = [list(map(str,k)) for k in data]

    for i in data:
        for j in i:
           out1.write(j+'\t')
        out1.write('\n')
    out1.close()
    return data

# =============================================================================
# def COL(name): # name = opts.c
#     data = {}
#     with open(name,'rt') as df:
#         for cc in [d.strip().split('\t') for d in df]:
#             data[cc[0]] = cc[1]
#     return data
# =============================================================================

# =============================================================================
# def gp_col(gp,lista): # gp,lista = opts.m,lista
#      #生成样本文件名，编号，颜色，的list
#     if gp != 'none':
#         group=read_gp(gp)
#         dict_col={}
#         for i in group[0]:
#             dict_col[i]=color[group[0].index(i)-1]
#         yangben=[]
#         for i in lista: # i = lista[0]
#             colour=dict_col[group[1][i.split('.')[-2]]]
#             yangben.append([i.split('.')[-2],colour])
#     else:
#         yangben=[]
#         a=0
#         for i in lista:
#             yangben.append([get_name(i),color[a]])
#             a+=1
#     return yangben,group
# =============================================================================

def ggplot_alpha(Input,mode,gp,COL,sumary,NM): #  Input,mode,gp,COL,sumary = "rarefaction_plot_group_data.tsv",Data3,opts.m,opts.c,opts.b 
    #写R代码
    out=open(NM,'w')
    out.write('''
Input <- "%s"
mode <- "%s"
gp <- "%s"
COL <- "%s"
sumary <- "%s"
#library(data.table)
library(tidyverse)
library(stringr)

data1 <- read.table(Input,head= T ,sep="\t",comment.char = "",fileEncoding = "UTF-8")
if (str_detect(colnames(data1)[1],"ample")){
  colnames(data1)[1] <- c("SampleID")
}
# head(data1)
mycol <- c("#CD3333","#483D8B","#458B00","#EEC900","#BFBC3B","#6E8B3D","#00688B","#C10077","#CAAA76","#EEEE00","#458B00","#8B4513","#008B8B","#6E8B3D","#8B7D6B","#7FFF00","#CDBA96","#ADFF2F")
color=c("#FF0000","#00FF00","#0000FF","#8A2BE2","#FFFF00","#CD5555","#8B2323","#FF00FF","#98F5FF","#FFD700","#32CD32","#5CACEE","#6A5ACD","#7B68EE","#8B7500","#912CEE","#BA55D3","#BF3EFF","#CAFF70","#E066FF","#EE7621","#27408B","#C0FF3E","#FF83FA","#9A32CD","#FA8072","#FFC125","#B22222","#EEC900","#FF6EB4","#ADFF2F","#000000","#1874CD","#191970","#228B22","#458B00","#551A8B","#66CD00","#76EE00","#7D26CD","#8B0000","#8B636C","#8B6508","#8B8B00","#8DB6CD","#90EE90","#9370DB","#97FFFF","#9ACD32","#9AFF9A","#9B30FF","#9BCD9B","#9FB6CD","#A020F0","#A0522D","#A52A2A","#AB82FF","#B03060","#B0E2FF","#B23AEE","#B3EE3A","#B452CD","#B4CDCD","#B4EEB4","#B5B5B5","#B8860B","#B9D3EE","#BABABA","#BBFFFF","#BC8F8F","#BCD2EE","#BCEE68","#BDB76B","#C71585","#CD0000","#CD00CD","#CD6090","#CD6600","#CDCD00","#D02090","#D2691E","#DA70D6","#DAA520","#DBDBDB","#DEB887","#E0FFFF","#E9967A","#EE0000","#EE6AA7","#EE7600","#EE8262","#EEB422","#EEDC82","#EEEE00","#F08080","#F0E68C","#F4A460","#F5DEB3","#FF1493","#FF3030","#FF4500","#FF7F00","#FF82AB","#FFA500","#FFB90F","#FFD39B","#EEB4B4","#EEA2AD","#EE9A00","#EE799F","#EE5C42","#EE00EE","#EE30A7","#EE4000","#DEDEDE","#E0EEEE","#D2B48C","#D15FEE","#CDBE70","#CDB7B5","#CD9B1D","#CD69C9","#CD3333","#C6E2FF","#BFEFFF","#B0C4DE","#ADADAD","#A2CD5A","#A2B5CD","#9932CC","#8EE5EE","#8B8989","#8B7B8B","#8B5A00","#8B4726","#8B0A50","#87CEEB","#7FFF00")

if (COL != "none"){
  data2 <- read.table(COL,head= F ,sep="\t",comment.char = "",fileEncoding = "UTF-8")
  data1$group <- factor(data1$group,levels=unique(as.vector(data2[,1])))
  mycol <- as.vector(data2[,2])
}

if (mode == 'rarefaction'){
  xlab <- "Number of Sequences"
  ylab <- "Number of OTUs"
}else{
  xlab <- "Number of Sequences"
  ylab <- paste0(strsplit(mode,split='_')[2],"_Diversity")
}

mytheme1 <- theme_bw() + 
  theme(panel.grid.major=element_line(color="white"),
        panel.grid.minor=element_line(color="white"),
        plot.title = element_text(hjust = 0.5),
        legend.title=element_blank(),
        legend.background=element_rect(fill="white", colour="black"),
        legend.key=element_blank())

if (gp == "none"){ # 不包含分组
  p <- ggplot(data1, aes(x = num, y = count, group = SampleID, color = SampleID)) + 
    geom_line(size=0.8) +
    scale_colour_manual(values=color) +
    xlab(xlab)+ylab(ylab) +
    theme_bw() + 
    theme(panel.grid.major=element_line(color="white"),
          panel.grid.minor=element_line(color="white"),
          plot.title = element_text(hjust = 0.5),
          legend.title=element_blank(),
          legend.background=element_rect(fill="white", colour="white"),
          legend.key=element_blank())
  ggsave(paste0(mode,".pdf"), p, width = 8, height = 5)
}else if (gp != "none" & sumary == "False"){
  p <- ggplot(data1, aes(x = num, y = count, group = SampleID,color = group)) + 
    geom_line(size=0.8) +
    scale_colour_manual(values=mycol) +
    xlab(xlab)+ylab(ylab) +
    mytheme1
  ggsave(paste0(mode,"_groups.pdf"), p, width = 6, height = 5)
}else if (gp != "none" & sumary != "False"){
  # 按样品分组，按组上色
  p1 <- ggplot(data1, aes(x = num, y = mean, group = group, color = group)) + 
    geom_line(size=0.8) +
    scale_colour_manual(values=mycol,name="Group") +
    xlab(xlab)+ylab(ylab) +
    mytheme1
  ggsave(paste0(mode,"_mean.pdf"), p1, width = 6, height = 5)
  # 加误差棒
  p2 <- p1 + geom_errorbar(aes(ymin=mean-std, ymax=mean+std), width=.5)
  ggsave(paste0(mode,"_errorbar.pdf"), p2, width = 6, height = 5)
  
  # 按样品分组，按组上色
  p3 <- ggplot(data1, aes(x = num, y = median, group = group, color = group)) + 
    geom_line(size=0.8) +
    scale_colour_manual(values=mycol,name="Group") +
    xlab(xlab)+ylab(ylab) +
    mytheme1
  ggsave(paste0(mode,"_median.pdf"), p3, width = 6, height = 5)
}
'''%(Input,mode,gp,COL,sumary))
    out.close()
    
if __name__ == "__main__":
    opts=option()    
    lista = dirs2list(opts.i)
    Data1 = draw_data(lista,opts.i+'_data.tsv','numsampled',"0.97")      # 提取画图数据
      
    # 编辑长数据
    Data2 = []
    for i in range(1,len(Data1)): # i = 1
        for j in range(1,len(Data1[i])): # j = 1
            Data2.append([Data1[i][0]] + [Data1[0][j]] + [Data1[i][j]])
    Data3 = pd.DataFrame(Data2)
    Data3.columns = ["#SampleID","num","count"]
    Data3.to_csv(opts.i+'_plot_data.tsv',sep="\t",encoding='utf-8',index=False)

    if opts.m != "none":
        # 分组
        sample_group = read_gp(opts.m)[2]
        # 根据分组筛选
        Data4 = pd.merge(Data3,sample_group,how='right')
        Data4[['num','count']] = Data4[['num','count']].astype(float)
        if opts.b == 'False':
            Data4.to_csv(opts.i+'_plot_group_data.tsv',sep="\t",encoding='utf-8',index=False)
            ggplot_alpha(opts.i+'_plot_group_data.tsv',opts.i,opts.m,opts.c,opts.b,'plot_group_cmd.r')
            os.system('''Rscript plot_group_cmd.r''')
        else:
            # summary 
            Data4_mean = Data4[['group','num','count']].groupby(['group','num']).mean().rename(columns={'count':'mean'}).reset_index()
            Data4_std = Data4[['group','num','count']].groupby(['group','num']).std().rename(columns={'count':'std'}).reset_index()
            Data4_median = Data4[['group','num','count']].groupby(['group','num']).median().rename(columns={'count':'median'}).reset_index()
            Data5 = pd.concat([Data4_mean,Data4_std['std'],Data4_median['median']],axis=1).fillna(0)
            Data5.to_csv(opts.i+'_plot_group_summary_data.tsv',sep="\t",encoding='utf-8',index=False)
            ggplot_alpha(opts.i+'_plot_group_summary_data.tsv',opts.i,opts.m,opts.c,opts.b,'plot_group_summary_cmd.r')
            os.system('''Rscript plot_group_summary_cmd.r''')        
    else:
        ggplot_alpha(opts.i+'_plot_data.tsv',opts.i,opts.m,opts.c,opts.b,'plot_data_cmd.r')
        os.system('''Rscript plot_data_cmd.r''')
        
# t.ppf(1-(1-alpha)/2, N-1) * np.sqrt(1+1/N) * mean
