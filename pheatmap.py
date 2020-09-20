# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 09:25:39 2018

@author: Colin Liu
"""

import getopt,sys,os

import argparse
parser= argparse.ArgumentParser(description="""
说明：
1. 画热图，pheatmap包
2. 距离计算方法可以取"correlation"(即peason), "euclidean", "maximum", "manhattan"
3. 聚类方法默认为 "complete",其他包括"ward", "single", "average", "mcquitty", "median", "centroid"
4. 输入文件和输出文件名必填
5. 修改方格的颜色有如下选择，输入名字如"wb"：
   rdb=c("navy", "white", "firebrick3"),
   rg=c("red","green"),
   wb=c("white","steelblue"),
   wg=c("white","lightgreen")
6. 添加color.txt文件修改颜色
""")

parser.add_argument('-i',help="输入文件",default="")
parser.add_argument('-tran',help="是否要转置，1要；0不要",default=0)
parser.add_argument('-cgc',help="是否要把行数据合并,导入文件",default="none")
parser.add_argument('-cd',help="分组文件名,map-grup.txt",default="none")
parser.add_argument('-rtop',help="挑选丰度前",default=50)
parser.add_argument('-cs',help="是否要标准化",default=1)
parser.add_argument('-rs',help="是否要标准化",default=0)
parser.add_argument('-o',help="输出文件名.如Heatmap.pdf",default="")
parser.add_argument('-drows',help="行距离计算方法",default="rank_correlation")
parser.add_argument('-dcols',help="列距离计算方法",default="bray")
parser.add_argument('-ccolor',help="定义每个小方格的颜色",default="")
parser.add_argument('-fsc',help="根据样本数量判断col字体大小",default=8)
parser.add_argument('-fsr',help="根据样本数量判断row字体大小",default=8)
parser.add_argument('-cw',help="每个小格子的row",default=8)
parser.add_argument('-ch',help="每个小格子的col",default=8)
parser.add_argument('-scale',help="是否标准化,row,column and none",default="none")
parser.add_argument('-cm',help="聚类方法",default="complete")
parser.add_argument('-clust_r',help="是否画行聚类树",default="T")
parser.add_argument('-clust_c',help="是否画列聚类树",default="T")
parser.add_argument('-cc',help="指定颜色：color.txt",default="none")

args= parser.parse_args()

cmd=open('pheatmapcmd.r','w')
cmd.write('''
          
################################### 定义变量 #####################################################
pm_i = "%s"
pm_tran = %s
pm_cgc = "%s"
pm_cd = "%s"
pm_rtop = %s
pm_cs = %s
pm_rs = %s
pm_o = "%s"
pm_drows = "%s"
pm_dcols = "%s"
pm_ccolor = "%s"
pm_fsc = %s
pm_fsr = %s
pm_cw = %s
pm_ch = %s
pm_scale = "%s"
pm_cm = "%s"
pm_clust_r = %s
pm_clust_c = %s
pm_cc = "%s"
################################## 读取OTU表格 ###################################################

otu <-read.table(file=pm_i,header=T,check.names=FALSE,sep="\t",fill=T)
rownames(otu) <- otu[,1]
rownames(otu) <-sapply(rownames(otu),function(x) gsub("_*{.+}","",x,perl = TRUE))
otu <-otu[,-1]
if(pm_tran==1){ otu <- t(otu) }#默认横轴为OTU，纵轴为物种，格式不符的可以转置一下下

################################ X轴名字旋转 #######################################################
## Edit body of pheatmap:::draw_colnames, customizing it to your liking
library("grid")
draw_colnames_45 <- function (coln, ...) {
  m = length(coln);  x = (1:m)/m - 1/2/m
  textGrob(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
           hjust = 1, rot = 90, gp = gpar(...))  ## Was 'hjust=0' and 'rot=270'
}
## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

############################## 是否根据分组合并列数据,顺便可以改一下名字 #######################

if(pm_cgc!="none"){
  group <- read.table(pm_cgc)
  glst <- lapply(1:length(unique(group[,2])),function(x)group[which(group[,2] %%in%% unique(group[,2])[x]),1])
  names(glst) <-unique(group[,2])
  tab <-sapply(1:length(glst),function(x) apply(otu[as.character(as.vector(glst[[x]]))],1,sum))
  otu <-tab[apply(tab,1,function(x)any(x>0)),]      
  colnames(otu) <-unique(group[,2])
}

################################ 分组map-group.txt ##################################################
mycol <-c("#CD3333","#483D8B","#458B00","#EEC900","#BFBC3B","#6E8B3D","#00688B","#C10077","#CAAA76","#EEEE00","#458B00","#8B4513","#008B8B","#6E8B3D","#8B7D6B","#7FFF00","#CDBA96","#ADFF2F")

if(pm_cd!="none"){
  de <-read.table(pm_cd,sep="\t",head= T,comment.char = "")
  #如果有color.txt文件
  if(pm_cc != "none"){
    sc <- read.table(pm_cc,sep="\t",comment.char = "",check.names = FALSE)
    sc <- sc[which(sc[,1] %%in%% unique(de$group)),]
    mycol <- as.vector(sc$V2)
    de$group <- factor(de$group,levels=as.vector(sc[,1]))
  }else{
    de$group <- factor(de$group,levels=as.vector(unique(de$group)))
  }
  
  dictpcol <- mycol[1:length(unique(de[,2]))]
  names(dictpcol) <- unique(de[,2])
  anncol = list(group=dictpcol)
  rownames(de) <- de[,1]
  mapping <- as.data.frame(de$group)
  rownames(mapping) <- rownames(de)
  colnames(mapping) <- colnames(de)[2]
  #	de <-de[order(de[,2]),]
  otu <-otu[,as.vector(de[,1])] # otu sort as de
}

otu <-otu[apply(otu,1,function(x)!all(x==0)),] 
otu <-otu[,apply(otu,2,function(x)!all(x==0))] 

al <- which(rownames(otu) %%in%% c("All"))
if(length(al)){ otu <-otu[-al,] }

############################## 是否挑选丰度排在前面的rtop个物种 #################################
if(pm_rtop >0){
  rsum <-sapply(1:nrow(otu),function(x) sum(otu[x,]))
  otu<-otu[order(rsum,decreasing=TRUE),]
  if(pm_rtop<=nrow(otu)){
    otu<-otu[1:pm_rtop,]}
}

############################# 是否标准化 #########################################

otu_hmp <-otu
if(pm_cs>0){
  otu_hmp <-apply(otu,2,function(x) x/sum(x)*pm_cs) 
}
if(pm_rs>0){
  otu_hmp <-t(apply(otu,1,function(x) x/sum(x)*pm_rs) )
}
colnames(otu_hmp) <-colnames(otu)
rownames(otu_hmp) <-rownames(otu)

ord_row <-c(nrow(otu_hmp):1)
ord_col <-c(1:ncol(otu_hmp))

otu_hmp <-otu_hmp[ord_row,ord_col] 
nr <- dim(otu_hmp)[1]
nc <- dim(otu_hmp)[2]
rowlab <-rownames(otu_hmp)
collab <-colnames(otu_hmp)
#print(otu_hmp)
outname <- paste(pm_o,".xls",sep="")
write.table(otu_hmp,outname,sep="\t",eol="\n",quote=FALSE)

#################################### 选择距离算法 #######################
dist_method <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "correlation" ,"rank_correlation" ,"bray")

if (pm_drows == "correlation"){
  drows = dist(otu_hmp,method = pm_drows)
} else if (pm_drows  == "rank_correlation"){
  library(vegan)
  dat.cor <- cor(t(otu_hmp),method="spearman")
  dat.cor.dist <-as.dist(1-dat.cor)
  write.table(file="dat.cor.xls",as.matrix(dat.cor),sep="	",quote=F)
  write.table(file="dat.cor.dist.xls",as.matrix(dat.cor.dist),sep="	",quote=F)
  drows <- dat.cor.dist
} else if (pm_drows  == "bray") {
  library(vegan)
  dat.t.bray <-vegdist(t(otu_hmp),method="bray")
  drows <- dat.cor.bray
} else {
 drows = dist(otu_hmp,method = pm_drows)
}

if (pm_dcols == "correlation"){
  dcols = dist(otu_hmp,method = pm_dcols)
} else if (pm_dcols  == "rank_correlation"){
  library(vegan)
  dat.cor <- cor(t(otu_hmp),method="spearman")
  dat.cor.dist <-as.dist(1-dat.cor)
  write.table(file="dat.cor.xls",as.matrix(dat.cor),sep="	",quote=F)
  write.table(file="dat.cor.dist.xls",as.matrix(dat.cor.dist),sep="	",quote=F)
  cols <- dat.cor.dist
} else if (pm_dcols  == "bray") {
  library(vegan)
  dat.t.bray <-vegdist(t(otu_hmp),method="bray")
  dcols <- dat.t.bray
} else {
  dcols = dist(otu_hmp,method = pm_dcols)
}

#################################### 选择聚类方法 ######################

Cluster_Method<-c( "ward", "single", "complete", "average", "mcquitty", "median", "centroid")

#################################### 方格颜色定义 ######################

if (pm_ccolor != ""){
  color_df <- list(rdb=c("navy", "white", "firebrick3"),
                   rg=c("red","green"),
                   wb=c("white","steelblue"),
                   wg=c("white","lightgreen")
  )
  cell_color = colorRampPalette(unlist(color_df[pm_ccolor]))(100)
}else{
  if (!require("RColorBrewer")) 
    install.packages("RColorBrewer")
  library(RColorBrewer)
  cell_color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
}

#################################### 大小定义 #############################
#根据样本数量判断col字体大小
if (pm_fsc==''){
  samplenum=ncol(otu_hmp)
  pm_fsc= 8
  if (samplenum>=20 ){pm_fsc='9'  }
  if (samplenum>=40 ){pm_fsc='8'  }
  if (samplenum>=45 ){pm_fsc='7.5'}
  if (samplenum>=50 ){pm_fsc='7'  }
  if (samplenum>=60 ){pm_fsc='6.5'}
  if (samplenum>=80 ){pm_fsc='6'  }
  if (samplenum>=100){pm_fsc='5'  }
  if (samplenum>=160){pm_fsc='9'  }
  if (samplenum>=180){pm_fsc='3'  }
  if (samplenum>=200){pm_fsc='2'  }
  if (samplenum>=300){pm_fsc='1'  }
}
pm_fsc<-as.numeric(pm_fsc)

#根据样本数量判断row字体大小
if (pm_fsr==''){
  samplenum=nrow(otu_hmp)
  pm_fsr= 8
  if (samplenum>=20 ){pm_fsr='9'  }
  if (samplenum>=40 ){pm_fsr='8'  }
  if (samplenum>=45 ){pm_fsr='7.5'}
  if (samplenum>=50 ){pm_fsr='7'  }
  if (samplenum>=60 ){pm_fsr='6.5'}
  if (samplenum>=80 ){pm_fsr='6'  }
  if (samplenum>=100){pm_fsr='5'  }
  if (samplenum>=160){pm_fsr='9'  }
  if (samplenum>=180){pm_fsr='3'  }
  if (samplenum>=200){pm_fsr='2'  }
  if (samplenum>=300){pm_fsr='1'  }
}
pm_fsr<-as.numeric(pm_fsr)

##判断每个方块的长宽
pm_cw=8
if (nrow(otu_hmp)<=200){cellhw=8}
if (nrow(otu_hmp)<=100){cellhw=10}
if (nrow(otu_hmp)<=50){cellhw=12}
if (nrow(otu_hmp)<=10){cellhw=20}


##判断每个方块的长宽
pm_ch=8
if (nrow(otu_hmp)<=200){cellhw=8}
if (nrow(otu_hmp)<=100){cellhw=10}
if (nrow(otu_hmp)<=50){cellhw=12}
if (ncol(otu_hmp)<=10){cellhw=20}

##################################### pheatmap ##########################
if (!require("pheatmap")) 
  install.packages("pheatmap")
library(pheatmap)

if (file.exists(pm_cd)){
  pheatmap(otu_hmp,#导入的数据
           color = cell_color,##颜色设置
           scale = pm_scale,
           cluster_row = pm_clust_r,#行是否聚类
           cluster_col = pm_clust_c,#列是否聚类
           clustering_distance_row = drows,
           clustering_distance_col = dcols,
           annotation = mapping,
           annotation_colors =anncol,
           main = "Heatmap",
           fontsize_row = pm_fsr,
           fontsize_col = pm_fsc,
           cellheight = pm_ch,
           cellwidth = pm_cw,
           filename = pm_o
  )
}else{
  pheatmap(otu_hmp,#导入的数据
           color = cell_color,##颜色设置
           scale = pm_scale,
           cluster_row = pm_clust_r,#行是否聚类
           cluster_col = pm_clust_c,#列是否聚类
           clustering_distance_row = drows,
           clustering_distance_col = dcols,
           main = "Heatmap",
           fontsize_row = pm_fsr,
           fontsize_col = pm_fsc,
           cellheight = pm_ch,
           cellwidth = pm_cw,
           filename = pm_o
  )
}
'''%(args.i,args.tran,args.cgc,args.cd,args.rtop,args.cs,args.rs,args.o,args.drows,args.dcols,args.ccolor,args.fsr,args.fsc,args.cw,args.ch,args.scale,args.cm,args.clust_r,args.clust_c,args.cc))

cmd.close()

os.system('Rscript pheatmapcmd.r')