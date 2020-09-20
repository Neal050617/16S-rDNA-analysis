#!/usr/bin/env python
# last updata 20150519 by yuguo
import argparse
import os

def option():
	par= argparse.ArgumentParser()
	par.add_argument("-i",metavar='[infiles]',required=True, help='Input abundance profile file')
	par.add_argument("-o",metavar='[ouput]',help="Output file name.",required=True)	
	par.add_argument("-g",metavar='[group file]',help="goup file with two columns ,first clo are samples,second indicates group label.",required=True)
	par.add_argument("-c",metavar='[c1-c2]',help="two group label to compare,sep with '-'. ",required=True)
	args = par.parse_args()	
	return args

def wilcox(infile,outfile,gfile,c1,c2):
	cmd='''
#library(qvalue)
source(\"/work/scripts/16s/qvalue.R\")
data <- read.table("'''+infile+'''",sep="\\t")
samp <- t(data[1,-1])
head=as.character(data[1,1])
data <-data[-1,]
rownames(data)<-data[,1]
data <-data[,-1]
colnames(data)<-samp

group=read.table("'''+gfile+'''",sep="\\t")

g1="'''+c1+'''"
g2="'''+c2+'''"
gsamp=group[which(group[,2] %in% c(g1,g2)),1]

gsamp1=group[which(group[,2] %in% g1),1]
gsamp2=group[which(group[,2] %in% g2),1]

da <- data
data <- apply(da,2,function(x) as.numeric(x)/sum(as.numeric(x)))
rownames(data)<-rownames(da)

data <-data[,which(samp %in% gsamp)]
samp <-samp[which(samp %in% gsamp)]

data <-data[apply(data,1,function(x)any(x>0)),]

set.seed(20190322)
out<-matrix(nrow=nrow(data),ncol=6)
for(i in 1:nrow(data)){
  d1=as.numeric(as.vector(unlist(data[i,which(samp %in% gsamp1)])))
  d2=as.numeric(as.vector(unlist(data[i,which(samp %in% gsamp2)])))
  wt <-wilcox.test(d1,d2,exact=F)
  me1 <-mean(d1)
  me2 <-mean(d2)
  sd1 <-sd(d1)
  sd2 <-sd(d2)
  out[i,]=c(rownames(data)[i],me1,sd1,me2,sd2,wt$p.value)
}
qv=qvalue(as.numeric(out[,6]),lambda=0.05)
out <-cbind(out,qv$qvalues)
colnames(out)=c(" ",paste("mean(",g1,")",sep=''),paste("sd(",g1,")",sep=''),paste("mean(",g2,")",sep=''),paste("sd(",g2,")",sep=''),"p-value","q-value")
out_order=out[order(out[,6]),]
write.table(out_order,"'''+outfile+'''",sep="\\t",col.names=T,row.names=F)
out_choosed=out_order[which(out_order[,6]<0.05 & out_order[,7]<0.01),]
data_choosed=data[out_choosed[,1],]
data_choosed=cbind(rownames(data_choosed),data_choosed)
colnames(data_choosed)[1]=head
write.table(data_choosed,paste("'''+outfile+'''",".filtered.xls",sep=""),sep="\t",col.names=T,row.names=F)
	'''
	#print(cmd)
	cmdfile=open("tmp.r",'w')
	cmdfile.write(cmd)
	cmdfile.close()
	os.system("Rscript tmp.r")

if __name__ == "__main__":
	### get options
	opts=option()
	wilcox(opts.i,opts.o,opts.g,opts.c.split('-')[0],opts.c.split('-')[1])
