# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 20:51:01 2019

@author: Administrator
"""

import argparse
import os
import shutil
import traceback

def option():
    parser= argparse.ArgumentParser(description="""
    说明：
    将多行fasta修改为一行序列名加上1000000
    """)
    parser.add_argument('-m', help="输入文件名", default='merge.list')
    args= parser.parse_args()
    return args


def readtable(filename):    
    '''读取table文件 .append or .extend '''  
    data = []
    for ln in open(filename,'rt'):
        data.append(ln.strip().split('\t'))
    return data


def Ln(lnn):#lnn = '>C70121:54:000000000-CL4MP:1:1101:13184:2054;ee=0.33;'
    ll = lnn.split(";ee")
    return ll

# def Switch_filter(faname):#faname = Merge_list[0][0]
#     Data = []
#     nn = 0
#     na = seq = ''
#     with open("data/filter_fa/"+faname+"_filtered.fasta",'rt') as df:
#         for i in [h.strip() for h in df]:
#             if i[0] == '>':
#                 Data.append(seq)
#                 seq = ''
#                 nn+=1
#                 na = faname+"_1000000"+str(nn)+';ee'+Ln(i)[1]
#                 Data.append('>'+na)
#             else:
#                 seq+=i    
#     return Data[1:]

def Switch_filter(faname):  #fasta文件名
    '''读取fasta'''
    data = []
    name = seq = ''
    nn = 0
    with open("data/filter_fa/"+faname+"_filtered.fasta",'rt') as df:
        for i in [h.strip() for h in df]:
            if i[0]=='>':
                data.append([name,seq])
                name=">"+faname+"_1000000"+str(nn)+';ee'+Ln(i)[1]
                seq=''
                nn+=1
            else:
                seq+=i
        data.append([name,seq])
    return sum(data[1:],[])

def write_fa(fa,name):
    out=open(name,'w')
    for j in fa:
        out.write(j+'\n')
    out.close()
    
def move_file(src_path, dst, file):
    try:
        # cmd = 'chmod -R +x ' + src_path
        # os.popen(cmd)
        f_src = os.path.join(src_path, file)
        f_dst = os.path.join(src_path, dst,file)
        shutil.move(f_src, f_dst)
    except Exception as e:
        print('move_file ERROR: ',e)
        traceback.print_exc()
                        
if __name__ == "__main__":#
    opts=option()
    
    Merge_list = readtable(opts.m)
    for i in [j[0] for j in Merge_list]:
        write_fa(Switch_filter(i),i+"_filtered.fasta")
    
    os.system('''
              mv data/filter_fa data/initial_filter_fa
              mkdir data/filter_fa
              mv *_filtered.fasta data/filter_fa              
    ''')
