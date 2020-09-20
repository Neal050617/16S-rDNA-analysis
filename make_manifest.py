# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 10:17:16 2019

@author: Administrator
"""

import os,re

re_digits = re.compile(r'(\d+)') 
def emb_numbers(s): # s = os.listdir(os.getcwd())[0]
    pieces=re_digits.split(s) 
    pieces[1::2]=map(int,pieces[1::2]) 
    return pieces 

def sortnumbers(alist):  
    return sorted(alist, key=emb_numbers)

def REFQ(line):
    # line = "C26_1_R1.fastq.gz"
    matchObj = re.match( r'(.*)_R([1-2]).f.*q.*', line, re.M|re.I)
    if matchObj:
        out1 = matchObj.group(1)
        if int(matchObj.group(2)) == 1:
            out2 = 'forward'
        else:
            out2 = 'reverse'
    else:
        print ("No match!!")
    return [out1,out2]

#i = os.listdir(os.getcwd())[0]
dirs = [[REFQ(i)[0],os.path.abspath(i),REFQ(i)[1]] for i in os.listdir(os.getcwd()) if ".fastq" in i or ".fq" in i]
sorted_list = sortnumbers([i[1] for i in dirs])
dict_dirs = {}
for i in dirs:
    dict_dirs[i[1]] = i

dirs2 = []
for k in sorted_list:
    dirs2.append(dict_dirs[k])

out=open("manifest.csv",'w')
out.writelines(r"sample-id,absolute-filepath,direction" + "\n")
for i in dirs2 :
       out.write((','.join(i)+'\n'))
out.close()
