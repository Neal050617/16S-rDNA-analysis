# -*- coding: utf-8 -*-
"""
Created on Sat Jun  9 14:48:34 2018

@author: Administrator
"""

import sys,os,time,shutil
import argparse
parser= argparse.ArgumentParser(description="""
说明：
PCA
NMDS
PCoA
color:group\tcolor
shape:sampleID\tgroup\tshaplegroup
20191107添加横竖误差棒选项0,1,2，默认0
""")

parser.add_argument('-i', help="输入文件", default='weighted_unifrac_dm.txt')
parser.add_argument('-md', help="NMDS,PCA,PCoA", default='PCoA')
parser.add_argument('-pc', help="1-2,轴", default='1-2')
parser.add_argument('-map', help="分组文件", default='none')
parser.add_argument('-map2',help="第二个分组文件、香农指数或者组内层次",default='none')
parser.add_argument('-gradt',help="渐变色设置red-white-blue,",default='none')
parser.add_argument('-col', help="指定颜色", default='none')
parser.add_argument('-lab', help="显示样本名", default='F')
parser.add_argument('-e', help="聚类圈", default='F')
parser.add_argument('-et', help="聚类圈的置信度", default=0.8)
parser.add_argument('-ep',help="聚类算法", default='t')
parser.add_argument('-base',help="碎石图", default='F')
parser.add_argument('-shp',help="形状", default='none')
parser.add_argument('-lp',help="legend位置", default='right')
parser.add_argument('-bx',help="加箱线图", default='F')
parser.add_argument('-sh',help="渐变",default='none')
parser.add_argument('-ct',help="误差十字架",default=0)

args= parser.parse_args()

cmd=open('cmd.r','w')
cmd.write(
'''
Input = "%s"
mode = "%s"
p_c = "%s"
map = "%s"
map2 = "%s"
gradt = "%s"
color = "%s"
label = "%s"
e = "%s"
et = %s
ep = "%s"
base = "%s"
shape = "%s"
legp = "%s"
bx = "%s"
shannon = "%s"
ct = %s

package_list <- c("vegan","reshape2","ggplot2","dplyr","tidyr","RColorBrewer","gridExtra","grid")

for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos="http://cran.r-project.org")
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

mycol <-c("#CD3333","#483D8B","#458B00","#EEC900","#BFBC3B","#6E8B3D","#00688B","#C10077","#CAAA76","#EEEE00","#458B00","#8B4513","#008B8B","#6E8B3D","#8B7D6B","#7FFF00","#CDBA96","#ADFF2F")
mypch <-c(16,15,17,18,7,8,9,10,11,12,13,14,21,22,23,24,25,3,4)

library("vegan")
library("reshape2")
library("ggplot2")
library("dplyr")
library("tidyr")
library("RColorBrewer")
library("gridExtra")
library("grid")

#确定分析的轴
pc_num =as.numeric(unlist(strsplit(p_c,"-")))
pc_x =pc_num[1]
pc_y =pc_num[2]
pccc<- paste0("pc",pc_x,"_",pc_y)

#读取文件_PCoA矩阵
da <-read.table(Input,sep="	",head=T,check.names=F,row.names = 1,stringsAsFactors = F)
#是否存在group
if (map!="none"){
  sd <- read.table(map,head=T,sep="	",comment.char = "",check.names = FALSE,stringsAsFactors = F)
  rownames(sd) <- as.character(sd[,1]);colnames(sd) <- c("sampleID","group")
  legend <- as.matrix(unique(sd$group))
  MM <- vector(length(legend),mode="list")#存储不同条件下的map
  FF <- vector(length(legend),mode="list")#用来存储因子
  for (i in 1:length(legend)){
    MM[[i]] <- sd[which(sd$group %%in%% legend[i,1]),]
    FF[[i]] <- unique(MM[[i]]$group)
  }
  if (mode == "PCA"){
    da <- da[,rownames(sd)]
  }else{
    da <- da[rownames(sd),rownames(sd)]
  }
}

if (shannon!="none"){
  Shan <- read.table(shannon,head=T,sep="	",comment.char = "",check.names = FALSE,stringsAsFactors = F)
  rownames(Shan) <- as.character(Shan[,1]);colnames(Shan) <- c("sampleID","shannon")
  if (mode == "PCA"){
    da <- da[,rownames(Shan)]
  }else{
    da <- da[rownames(Shan),rownames(Shan)]
  }
}

#如果存在第二个分组
if (map2!="none"){
  MAP <- read.table(map2,sep="	",stringsAsFactors = F)
  rownames(MAP) <- MAP[,1];colnames(MAP) <- c("sampleID","map2")
  MAP <- MAP[rownames(sd),]
  #字符或整数，每个分组8种levels以内
  for (j in 1:length(legend)){
    MM[[j]]$map2 <- MAP[which(MAP$sampleID %%in%% MM[[j]]$sampleID),2]
  }
  names(FF) <- as.vector(legend)
}


#设置颜色
if (color!="none"){
  sc <- read.table(color,sep="	",comment.char="",check.names=F,stringsAsFactors = F)
  sc <- sc[which(sc[,1] %%in%% legend),]
  rownames(sc) <- sc[,1]
  mycol <- as.vector(sc$V2)
}else{
  mycol <- mycol[1:length(legend)]
  if (gradt != "none"){mycol <- unlist(strsplit(gradt,"-"))}
}

#设置形状
if (shape!="none"){
  sp <- read.table(shape,sep="	",comment.char="",check.names=F,stringsAsFactors = F)
  sp <- sp[which(sp[,1] %%in%% sd[,1]),]
  sp[,2] <- factor(sp[,2],levels=unique(sp[,2]))
  mypch <- unique(sp[,3])
}

#合并数据
if (map!="none"){
  i=1
  MMM <- MM[[i]]
  while (i<length(MM)){
    MMM <- rbind(MMM,MM[[i+1]])
    i <- i+1
  }
}

#######################################################################################
#计算坐标
if(mode=="NMDS"){
  nmds <- metaMDS(da,k=2,trace=FALSE)
  pc12 <- as.data.frame(nmds$points)
  xlab=c("NMDS1")
  ylab=c("NMDS2")
  write.table(nmds$stress,paste0(mode,"_",Input,"_stress.xls"))
  write.table(nmds$points,paste0(mode,"_",Input,"_sites.xls"))
}else{
  if(mode=="PCA"){
    pca <- prcomp(da)
    pc12 <- as.data.frame(pca$rotation[,pc_num])
  }
  if(mode=="PCoA"){
    pca <- prcomp(da)
    pc12 <- as.data.frame(pca$x[,pc_num])
  }
  pc <-summary(pca)$importance[2,]*100
  xlab=paste("PC",pc_x,": ",round(pc[pc_x],2),"%%",sep="")
  ylab=paste("PC",pc_y," :  ",round(pc[pc_y],2),"%%",sep="")
  ##输出文件
  write.table(data.frame(sample_ID=rownames(pca$rotation),pca$rotation),file=paste0(mode,"_",Input,"_PC","_","rotation.xls"),quote=F,sep="	",row.names = F) #输出特征向量
  write.table(data.frame(sample_ID=rownames(predict(pca)),predict(pca)),file=paste0(mode,"_",Input,"_PC","_","sites.xls"),quote=F,sep="	",row.names = F) #输出新表
  pca.sum=summary(pca)
  write.table(data.frame(sample_ID=rownames(pca.sum$importance),pca.sum$importance),file=paste0(mode,"_",Input,"_PC","_","importance.xls"),quote=F,sep="	",row.names = F) #输出PC比重
}

#数据与分组合并
if (map!="none"){
  pc12$sampleID <- rownames(pc12)
  pc12 <- left_join(data.frame(pc12),MMM)
  rownames(pc12) <- pc12$sampleID
  colnames(pc12)[c(1,2)] <- c("PC1","PC2")
  #加上形状
  if (shape != "none"){
    colnames(sp) <- c("sampleID","group2","shape")
    pc12 <- left_join(pc12,sp)
  }
  #设置因子水平
  if (ncol(MM[[1]])>=2){
    if (color=="none"){
      pc12$group <- factor(pc12$group,levels=as.vector(legend))
    }else{
      pc12$group <- factor(pc12$group,levels=as.vector(sc[,1])) 
    }
  }
  
  if (ncol(MM[[1]])>=3){
    if (is.numeric(MMM[1,3])){
      pc12$map2 <- as.numeric(as.character(pc12$map2))
    }else{
      pc12$map2 <- factor(pc12$map2,ordered=T)
    }
  }
  
  if(shannon != "none"){
    pc12$sampleID <- rownames(pc12)
    pc12 <- left_join(data.frame(pc12),Shan)
    rownames(pc12) <- pc12$sampleID
    colnames(pc12)[c(1,2)] <- c("PC1","PC2")
    pc12$sampleID <- NULL
  }
  
  ##设置主题
  mytheme <- theme_bw() + 
    theme(panel.grid.major=element_line(color="white"),
          panel.grid.minor=element_line(color="white"),
          plot.title = element_text(hjust = 0.5),
          #legend.title=element_blank(),
          #legend.title=element_text("Group"),
          #legend.position=c(0,1),
          #legend.justification=c(0,1),
          legend.background=element_rect(fill="white", colour="black"),
          legend.key=element_blank(),
          legend.position=legp)
  
  centroids <- aggregate(cbind(PC1,PC2)~group,pc12,mean)
  f         <- function(z)sd(z)/sqrt(length(z)) # function to calculate std.err
  se        <- aggregate(cbind(se.x=PC1,se.y=PC2)~group,pc12,f)
  ct1 <- merge(pc12,aggregate(cbind(mean.x=PC1,mean.y=PC2)~group,pc12,mean),by="group")
  ct2 <- merge(centroids,se, by="group")# add std.err column to centroids
    
  if (ct == 1){
    p <- ggplot(ct1, aes(PC1,PC2,color=factor(group),shape=factor(group)))+
      geom_point(size=1)+
      geom_point(aes(x=mean.x,y=mean.y),size=3)+
      geom_segment(aes(x=mean.x, y=mean.y, xend=PC1, yend=PC2))+
      stat_ellipse(level=0.6)+ 
      scale_colour_manual(values=mycol,name="Group") +
      scale_shape_manual(values=mypch,name="Group")
    
  }else if (ct == 2){
    p <- ggplot(ct1, aes(PC1,PC2,color=factor(group),shape=factor(group)))+
      #geom_point(size=3)+
      geom_point(data=ct2, size=3)+
      geom_errorbar(data=ct2,aes(ymin=PC2-se.y,ymax=PC2+se.y),width=0.05)+
      geom_errorbarh(data=ct2,aes(xmin=PC1-se.x,xmax=PC1+se.x),height=0.05)+
      scale_colour_manual(values=mycol,name="Group") +
      scale_shape_manual(values=mypch,name="Group")
    
  }else if (ct == 0){
    if (map=="none" & shannon=="none"){
      p <- ggplot() +
        geom_point(data=pc12,mapping=aes(x=PC1,y=PC2,color='all'),size=3)+
        scale_colour_manual(values="blue",guide=FALSE)
    }else if (map=="none" & shannon!="none"){
      p <- ggplot() +
        geom_point(data=pc12,mapping=aes(x=PC1,y=PC2,color=shannon),size=3)+
        scale_colour_gradient2(low=mycol[1],mid=mycol[2],high=mycol[3],
                               midpoint=median(pc12$shannon),name="")
    }else if (map!="none" & shannon=="none"){
      if(ncol(MMM)==2){
        p<-ggplot()
        if (shape != "none"){
          p <- p+ 
            geom_point(data=pc12,aes(x=PC1,y=PC2,color=group,shape=group2),size=3)
        }else{
          p <- p + 
            geom_point(data=pc12,aes(x=PC1,y=PC2,color=group,shape=group),size=3)
        }
        p <- p + 
          scale_colour_manual(values=mycol,name="Group") +
          scale_shape_manual(values=mypch,name="Group")
        if (e=="T"){#是否加圈
          p <- p +
            stat_ellipse(data=pc12,mapping=aes(x=PC1,y=PC2,color=group),
                         level=et,type=ep)
        }
      }else{
        if(ncol(MMM)==3 & is.numeric(MMM[1,3])){
          p<-ggplot(pc12,aes(PC1,PC2,colour=map2)) +
            geom_point(size=3)+
            scale_colour_gradient2(low=mycol[1],mid=mycol[2],high=mycol[3],
                                   midpoint=median(pc12$map2),name="") +
            theme(legend.key.width =unit(8,'line'))
        }
        if(ncol(MMM)==3& is.character(MMM[1,3])){
          p<-ggplot() +
            geom_point(data=pc12,aes(x=PC1,y=PC2,color=group,
                                     shape=group,alpha=map2),size=3) +
            scale_colour_manual(values=mycol,name="Group") +
            scale_shape_manual(values=mypch,name="Group")
        }
        if (e=="T"){#是否加圈
          p <- p +
            stat_ellipse(data=pc12,mapping=aes(x=PC1,y=PC2,color=group),
                         level=et,type=ep,show.legend=TRUE)
        }
      }
    }
    
    #是否显示legend 
    if (label=="T"){
      p <- p+geom_text(data=pc12,mapping=aes(x=PC1, y=PC2),hjust=.5,vjust=-1,size = 3,
                       alpha=0.8,label=rownames(pc12))
    }
    
  }
  
  #坐标，内框线，标题，背景主题的统一设置  
  p <- p+labs(title=mode,x=xlab,y=ylab)+
    geom_vline(aes(xintercept=0),linetype="dotted")+
    geom_hline(aes(yintercept=0),linetype="dotted")+
    theme(legend.position='none')+
    mytheme
  
  if (bx=="T"){
    #导出legend的函数
    g_legend<-function(a.gplot){
      tmp <- ggplot_gtable(ggplot_build(a.gplot))
      Leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[Leg]]
      return(legend)
    }
    
    Legend <- g_legend(p)
    
    boxtheme <- theme_bw() + 
      theme(panel.grid.major=element_line(color="white"),
            panel.grid.minor=element_line(color="white"),
            plot.title = element_text(hjust = 0.5),
            legend.title=element_blank(),
            legend.background=element_rect(fill="white", colour="black"),
            legend.key=element_blank())
    
    p2 <- ggplot(data=pc12,aes(x=group,y=PC1)) + 
      geom_boxplot(aes(fill=group))+
      scale_fill_manual(values=mycol)+
      guides(fill=FALSE)+
      labs(x = "", y = "", title = "")+
      coord_flip()+
      scale_x_discrete(limits=rev(levels(pc12$group)))+
      boxtheme+
      theme(axis.text.x=element_blank(),
            axis.ticks.x = element_blank())
    
    p3 <- ggplot(data=pc12,aes(x=group,y=PC2)) + 
      geom_boxplot(aes(fill=group))+
      scale_fill_manual(values=mycol)+
      guides(fill=FALSE)+
      labs(x = "", y = "", title = "")+
      boxtheme+
      theme(axis.text.y=element_blank(),
            axis.ticks.y = element_blank())
    
    if (nlevels(pc12$group)/8 >0.5){
      hw <- 0.5
    }else{
      hw <- nlevels(pc12$group)/8
    }
    
    p<- grid.arrange((p+theme(legend.position = "none")),p3,p2,Legend,nrow=2,ncol=2,heights=c(1,hw),widths=c(1,hw))
  }
}

ggsave(paste0(mode,"_",Input,"_PC",p_c,".pdf"),p,width = 8.5,height = 8)

if (base=="T"){
  #接下来输出一些主要的记过表格以及柱状图和碎石图：
  pdf(file="pcaBarplot.pdf",width=15) #柱状图
  barplot(pca.sum$importance[2,1:10]*100,xlab="PC",ylab="percent",col="skyblue")
  dev.off()
  pdf(file="pcaPlot.pdf",width=15) #碎石图
  plot(pca.sum$importance[2,1:10]*100,type="o",col="red",xlab="PC",ylab="percent")
  dev.off()
}

'''%( args.i,args.md,args.pc,args.map,args.map2,args.gradt,args.col,args.lab,args.e,args.et,args.ep,args.base,args.shp,args.lp,args.bx,args.sh,args.ct))
cmd.close()
os.system('Rscript cmd.r')

