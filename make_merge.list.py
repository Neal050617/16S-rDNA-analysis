# -*- coding: utf-8 -*-

import getopt,sys,os,random,re
opts, args = getopt.getopt(sys.argv[1:], "hio:", ["help", "imput", "output="])
def usage():
        print("""
        -o 输出文件名 'merge.list'         
        """)
output='merge.list'

for o, a in opts:
        if o in ('-i','--input'):
             imput = a
        elif o in ('-o','--output'):
             output = a
        elif h in ('-h','help'):
                usage()
                sys.exit()
   
def unique(old_list):
    newList = []
    for x in old_list:
        if x not in newList :
            newList.append(x)
    return newList

re_digits = re.compile(r'(\d+)') 
def emb_numbers(s): # s = os.listdir(os.getcwd())[0]
    pieces=re_digits.split(s) 
    pieces[1::2]=map(int,pieces[1::2]) 
    return pieces 

def sortnumbers(alist):  
    return sorted(alist, key=emb_numbers)

path = os.getcwd()
list = os.listdir(path)#列出目录下的所有文件和目录

list2=[]
for ii in list:
    if '.fastq' or '.fq' in ii and 'R' in ii:
       list2.append(ii)

#排序
list2=sorted(list2)
list3=[]
for i in range(int(len(list2)/2)):
#    replace=re.sub('_S(\d)*_L001_R(\d)_001.fastq',r'',list2[2*i])
    replace=re.sub(r'_R([1-2]).f.*q.*',r'',list2[2*i])
    list3.append([replace,list2[2*i],list2[2*i+1]])

sorted_list3 = sortnumbers([i[0] for i in list3])

dict_list3 = {}
for i in list3:
    dict_list3[i[0]] = i

list4 = []
for k in sorted_list3:
    list4.append(dict_list3[k])


out=open(output,'w')
for i in list4 :
       out.write(('\t'.join(i)+'\n'))
out.close()

print(u'ok------ok,please Check order,请检查一下文件，并加上样本名称')