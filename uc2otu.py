# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 21:50:38 2018

@author: Administrator
"""

import argparse
import re
import pandas as pd

def option():
    parser= argparse.ArgumentParser(description="""
    说明：
    20190112-删除注释文件的最后一列；无根树注释文件
    """)
    parser.add_argument('-i', help="usearch uparse file", default='map.uc')
    parser.add_argument('-r', help="otu rep,不包含具体序列名;otu_reps.raw.fasta", default='otu_reps.raw.fasta')
    parser.add_argument('-t', help="taxonomy file;otu_reps_tax_assignments.txt", default='otu_reps_tax_assignments.txt')
    parser.add_argument('-f', help="fasta;meta.fasta", default='meta.fasta')
    parser.add_argument('-m', help="mode:otu:subsample|otu", default='otu')
    args= parser.parse_args()
    return args

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
    
def embedded_numbers(s): # s= files[0]
    re_digits = re.compile(r'(\d+)')
    pieces = re_digits.split(s)               # split number & string  
    pieces[1::2] = map(int, pieces[1::2])     # number to int  
    return pieces
    
def sort_sn(alist):
    #sort_strings_with_embedded_numbers
     return sorted(alist, key=embedded_numbers)
 
def write_fasta(flist,outname):
    #写入表格
    out=open(outname,'w')
    for i in flist:
        out.write('\n'.join(i)+'\n')
    out.close()
    
def write_data(flist,outname):
    #写入表格
    out=open(outname,'w')
    for i in flist:
        out.write(i[0]+"\t"+"\t".join(i[1])+'\n')
    out.close()
    
def analysis(data,rep,tax,fa,mode): # data,rep,tax,fa,mode = opts.i,opts.r,opts.t,opts.f,opts.m
    # read data
    d1 = pd.read_csv(data,sep='\t',header=None,comment='#')
    d2 = d1.loc[~d1.iloc[:,9].str.contains(u"\*"),[8,9]]
    d2["sample_ID"] = d2[8].apply(lambda x:x.split("_1000000")[0])
    d2['count1'] = int(1)
    d2.columns = ["seq_name","otu_name","sample_ID","count1"]
        
    # unique OTU & sample_ID
    ol = d2.loc[~d2.duplicated(['otu_name'],keep='first'),"otu_name"].tolist()
    sl = d2.loc[~d2.duplicated(["sample_ID"],keep='first'),"sample_ID"].tolist()
    
    # rep.fasta
    if rep != "none":
        write_fasta(read_fa(rep,sort_sn(ol)),mode+"_select_"+rep)
    
    # otu.fasta
    if fa != "none":
        write_fasta(read_fa(fa,d2["seq_name"].tolist()),mode+".fasta")
        
    # otu.groups
    d2.loc[:,["seq_name","sample_ID"]].to_csv(mode+'.groups',sep='\t',index=0,header=0)
        
    # otu_seqids.txt
    sn = []
    for i in ol: # i = ol[0]
        sn.append([i,d2.loc[d2["otu_name"].isin([i]),"seq_name"].tolist()])
    write_data(sn,mode+"_seqids.txt")    
    
    # reshape
    #d2.loc[:,["otu_name","sample_ID","count1"]].pivot(index="otu_name",columns="sample_ID",values="count1")
    d3 = d2.loc[:,["otu_name","sample_ID","count1"]].groupby(["otu_name","sample_ID"]).sum().reset_index()

    # groupby to dataframe
    dd = pd.DataFrame(0,columns=sl,index=ol)
    for i in ol: # i = ol[4]
        for j in sl: # j = sl[1]
            if ((d3["otu_name"]==i) & (d3["sample_ID"]==j)).any():
                dd.loc[i,j] = d3.loc[(d3["otu_name"]==i) & (d3["sample_ID"]==j),"count1"].tolist()[0]

    # subsample number
    if mode != "subsample":
        pd.DataFrame([min(dd.apply(lambda x: x.sum()))]).to_csv('subsample_num.txt',sep='\t',index=0,header=0) 
        
    # sort
    out = dd.loc[sort_sn(ol),sort_sn(sl)]
    out.insert(0,r"#OTU ID",out.index)
    out.to_csv(mode+'_table.xls',sep='\t',index=0,na_rep='')
    
    # shared
    shared = out.transpose()
    shared.drop(r"#OTU ID",axis=0,inplace=True)
    shared.insert(0,"numOtus",shared.apply(lambda x: len(x), axis=1))
    shared.insert(0,"Group",shared.index)        
    shared.insert(0,"label",0.97)
    shared.to_csv(mode+'.shared',sep='\t',index=0,na_rep='')
    
    # read tax 
    if tax != "none":
        t1 = pd.read_csv(tax,sep='\t',encoding='utf8',header=None)
        # 如果最后一列为数字，则删除最后一列
        if t1.iloc[:,-1].dtype ==  'float64':
            t1.drop(t1.shape[1]-1,axis=1,inplace=True)
        t1.columns = [r"#OTU ID","taxonomy"]
            
        t2 = pd.merge(out,t1,on=r"#OTU ID",how='left')
        t2.to_csv(mode+'_taxa_table.xls',sep='\t',index=0,na_rep='')
        
        if rep != "none":
            t2.loc[:,[r"#OTU ID","taxonomy"]].to_csv(mode+"_select_"+tax,sep='\t',index=0,header=None,na_rep='')
            
if __name__ == "__main__":
    opts=option()
    analysis(opts.i,opts.r,opts.t,opts.f,opts.m)
    print(('''生成：
          {}
          {}.fasta
          {}.groups
          {}_seqids.txt
          subsample_num.txt
          {}_table.xls
          {}.shared
          {}
          {}_taxa_table.xls
          ''').format(opts.m+"_select_"+opts.r,opts.m,opts.m,opts.m,opts.m,opts.m,opts.m+"_select_"+opts.t,opts.m))
