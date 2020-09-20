# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 08:59:58 2019

@author: Colin Liu
"""

import argparse
import pandas as pd
import re
import os,sys,copy
import shutil
import traceback

def option():
    parser= argparse.ArgumentParser(description="""
    说明：
    将多行fasta修改为一行序列名加上1000000
    """)
    parser.add_argument('-m', help="输入文件名", default='merge.list')
    parser.add_argument('-a', help="输入文件名", default='otu_reps_tax_assignments.txt')
    parser.add_argument('-s', help="usearch uparse file", default='otu_seqids.txt')
    parser.add_argument('-o', help="输出文件夹名", default='Rmseq')
    parser.add_argument('-r', help="删除的菌", default='rm.list')
    args= parser.parse_args()
    return args

def readtable(filename):    
    '''读取table文件 .append or .extend '''  
    data = []
    for ln in open(filename,'rt'):
        data.append(ln.strip().split('\t'))
    return data

def write_faq(fdict,outname):
    '''#写fq或fa'''
    out=open(outname,'w')
    for i in fdict.values():
        out.write('\n'.join(i)+'\n')
    out.close()
    
def write_fa(fdict,fa,outname):#fdict = data2
    '''#写fq或fa'''
    oo = []
    nn = 0
    if isinstance(fdict,dict):
        Values = fdict.values()
    elif isinstance(fdict,list):
        Values = copy.deepcopy(fdict)
    for i in Values:
        nn += 1
        oo.append(">"+fa+"_1000000"+str(nn)+';ee'+Ln(i[0])[1])
        oo.append(i[1])
        
    out=open(outname,'w')
    for j in oo:
        out.write(j+'\n')
    out.close()

def write_data(flist,outname):
    #写入表格
    out=open(outname,'w')
    for i in flist:
        out.write(i+'\n')
    out.close()

def Ln(lnn):#lnn = '>C70121:54:000000000-CL4MP:1:1101:13184:2054;ee=0.33;'
    ll = lnn.split(";ee")
    return ll
    
def Switch_name(faname):#faname = Merge_list[0][0] 
    Data = {}
    nn = 0
    with open("data/initial_filter_fa/"+faname+"_filtered.fasta",'rt') as df:
        for i in [h.strip() for h in df]:
            if i[0] == '>':
                nn+=1
                na = faname+"_1000000"+str(nn-1)+';ee'+Ln(i)[1]
                Data[na.split(";ee")[0]] = i[1:].split(";ee=")[0] 
    return Data
    
def read_seqids(ass,out2):#ass = args.s
    data = []
    with open(ass,'rt') as df:
        for i in [h.strip().split("\t") for h in df]:
            if i[0] in out2:
                data.append(i[1:])
    return data

def Split_LL(l1,out1):#l1 = out3
    data1 = pd.DataFrame([[i.split('_1000000')[0],i] for i in l1]).groupby([0])
    data1_1 = dict(list(data1))
    data2 = {}
    for key in data1_1:
        data2[key] = [out1[j].split(";ee=")[0] for j in data1_1[key][1].tolist()]        
    return data2

def Change_fq(fq,out5):# fq = list(out5)[0]
    with open('data/rawdata/'+fq+"_R1.fastq") as df:
        data1_1 = []
        for ln in df:
            data1_1.append(ln.strip())
        data1_2 = {}
        for i in range(int(len(data1_1)/4)):
            data1_2[data1_1[i*4][1:].split(' ')[0]] = data1_1[i*4:(i+1)*4]

    for key in out5[fq]:
        data1_2.pop(key)
        
    write_faq(data1_2,"Rmseq_"+fq+"_R1.fastq")
#############################################################
    with open('data/rawdata/'+fq+"_R2.fastq") as df:
        data2_1 = []
        for ln in df:
            data2_1.append(ln.strip())
        data2_2 = {}
        for i in range(int(len(data2_1)/4)):
            data2_2[data1_1[i*4][1:].split(' ')[0]] = data2_1[i*4:(i+1)*4]
            
    for key in out5[fq]:
        data2_2.pop(key)
                
    write_faq(data2_2,"Rmseq_"+fq+"_R2.fastq")

def Change_fa(fa,out5):# fa = list(out5)[0]
               
    data1 = []
    name = seq = ''
    with open("data/initial_filter_fa/"+fa+"_filtered.fasta",'rt') as df:
        for i in [h.strip() for h in df]:
            if i[0]=='>':
                data1.append([name,seq])
                name=i
                seq=''
            else:
                seq+=i
        data1.append([name,seq])
    del data1[0]     
    data2 = {}
    for i in range(int(len(data1))):#i = 0
         data2[data1[i][0][1:].split(';ee=')[0]] = data1[i]

    for key in out5[fa]:
        data2.pop(key)
        
    write_fa(data2,fa,"Rmseq_"+fa+"_filtered.fasta")

def fix_filter(faname):
    data = []
    name = seq = ''
    nn = 0
    with open("data/initial_filter_fa/"+faname+"_filtered.fasta",'rt') as df:
        for i in [h.strip() for h in df]:
            if i[0]=='>':
                data.append([name,seq])
                name=">"+faname+"_1000000"+str(nn)+';ee'+Ln(i)[1]
                seq=''
                nn+=1
            else:
                seq+=i
        data.append([name,seq])
    
    write_fa(data[1:],faname,"Remove_"+faname+"_filtered.fasta")

def move_file(srcfile,dstfile):
    if not os.path.isfile(srcfile):
        print ("%s not exist!"%(srcfile))

    try:
        fpath,fname=os.path.split(dstfile)    #分离文件名和路径
        if not os.path.exists(fpath):
            os.makedirs(fpath)                #创建路径
        shutil.move(srcfile,dstfile)
        print("move %s -> %s"%(srcfile,dstfile))
    except Exception as e:
        print('move_file ERROR: ',e)
        traceback.print_exc()
                           
def copy_file(srcfile,dstfile):
    if not os.path.isfile(srcfile):
        print("%s not exist!"%(srcfile))

    try:
        fpath,fname=os.path.split(dstfile)    #分离文件名和路径
        if not os.path.exists(fpath):
            os.makedirs(fpath)                #创建路径
        shutil.copy(srcfile,dstfile)
        print("copy %s -> %s"%(srcfile,dstfile))
    except Exception as e:
        print('copy_file ERROR: ',e)
        traceback.print_exc()

if __name__ == "__main__":
    opts=option()

    Merge_list = readtable(opts.m)    
    # i = Merge_list[0]
    
    path = os.getcwd()
    List = [i for i in os.listdir(path) if os.path.isdir(i)] #只列出文件夹
    try:
        pp = [s for s in List if 'analysis' in s]
        opts.a = os.path.join(path,pp[0],'process','otu_0.97',opts.a)
        opts.s = os.path.join(path,pp[0],'process','otu_0.97',opts.s)
    except Exception as e:
        print('analysis dictionary did not exists',e)
        traceback.print_exc()
        
    out1 = {}
    for i in Merge_list:
        out1.update(Switch_name(i[0]))
        
    out2 = []
    for j in sum(readtable(opts.r),[]):
        #out2.append([i[0] for i in readtable(opts.a) if j in i[1]])
        out2.append([i[0] for i in readtable(opts.a) if bool(re.search(j,i[1]))])
    out2 = sum(out2,[])
    #out2 = [i[0] for i in readtable(opts.a) if opts.r in i[1]]
    out3 = sum(read_seqids(opts.s,out2),[])
    out5 = Split_LL(out3,out1)
    
    #out4 = [i for i in out3 if re.search('A16_1000000|A12_1000000',i)]   
    #out5 = Split_LL(out4,out1) 
    for i in list(out5):
        Change_fa(i,out5)
        Change_fq(i,out5)
        Ld = os.getcwd()
        move_file(os.path.join(Ld,opts.o+"_"+i+"_R1.fastq"),os.path.join(Ld,opts.o,'rawdata',i+"_R1.fastq"))
        move_file(os.path.join(Ld,opts.o+"_"+i+"_R2.fastq"),os.path.join(Ld,opts.o,'rawdata',i+"_R2.fastq"))
        move_file(os.path.join(Ld,opts.o+"_"+i+"_filtered.fasta"),os.path.join(Ld,opts.o,'filter_fa',i+"_filtered.fasta"))

#    os.system('''
#              mkdir "%s" "%s"/filter_fa "%s"/rawdata
#              '''%(opts.o,opts.o,opts.o))
     
    out6 = set([i[0] for i in Merge_list]) - set(list(out5))
    if out6 == set():
        print("Nothing left")
    else:
        write_data(out6,'Left.list')

        for i in list(out6):#i = list(out6)[0]
            fix_filter(i)
            copy_file(os.path.join(Ld,'data','rawdata',i+"_R1.fastq"),os.path.join(Ld,opts.o,'rawdata',i+"_R1.fastq"))
            copy_file(os.path.join(Ld,'data','rawdata',i+"_R2.fastq"),os.path.join(Ld,opts.o,'rawdata',i+"_R2.fastq"))
            move_file(os.path.join(Ld,"Remove_"+i+"_filtered.fasta"),os.path.join(Ld,opts.o,'filter_fa',i+"_filtered.fasta"))

