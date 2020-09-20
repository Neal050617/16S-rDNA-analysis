# -*- coding: utf-8 -*-
"""
Created on Thu May  3 09:53:42 2018

@author: Colin Liu
"""

import sys,os,time,shutil
import argparse

parser= argparse.ArgumentParser(description="""
说明：这个脚本是。。。。。 
1.gglpot2画柱形图
2.输入分组后可以按分组顺序排列或画平均值
3.输入两列的table可以用来改名，也可以用来挑选样本
4.输入分组计算组间差异显著性，并画柱形图或箱式图
5.箱式图添加log值与10倍的选项；
6.group添加组的分类
7.可以设置组的颜色
20190409-facet_grid(.~group,scale = "free",space = "free") 不会自动改变宽度
20190428-调整柱形图高度，ycor，使较低丰度的物种能够显示出;修正输出p-sign为NA的问题;修正单个差异物种，柱形图无法显示。
""")

parser.add_argument('-i',help="输入文件",default="")
parser.add_argument('-rp',help="rename.pick,修改样本名，也可以挑选样本,两列无表头",default="")
parser.add_argument('-di',help="按丰度排序，d为降序，i为升序",default="d")
parser.add_argument('-n',help="是否挑选每个样本丰度前n个的物种",default=-1)
parser.add_argument('-m',help="是否挑选丰度综合排在前m个的物种",default=-1)
parser.add_argument('-p',help="0.01或0.03或其他小数",default=-1)
parser.add_argument('-oth',help="是否保留others",default="T")
parser.add_argument('-map',help="map-group文件",default="")
parser.add_argument('-ave',help="1为组间冲击图，2为组内排序递增",default=0)
parser.add_argument('-bi',help="挑选的物种列表",default="none")
parser.add_argument('-pt',help="1为条形图，2为箱式图，p值检验",default=0)
parser.add_argument('-bb',help="1为条形图，2为气泡图，丰度图",default=1)
parser.add_argument('-col',help="colot.txt",default="none")
parser.add_argument('-Lg',help="箱式图log",default="T")

args= parser.parse_args()

cmd=open('gg_barplotcmd.r','w')
cmd.write('''
# 清理工作环境 clean enviroment object
#rm(list=ls())

bp_i	= "%s"
bp_rp	= "%s"
bp_di	= "%s"
bp_n	= %s
bp_m	= %s
bp_p	= %s
bp_oth	= "%s"
bp_map	= "%s"
bp_ave	= "%s"
bb_i = "%s"
pt = "%s"
bb = "%s"
paint = "%s"
Lg = "%s"


##检查依赖关系包的安装
package_list <- c("reshape2","ggplot2","dplyr")

for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos="http://cran.r-project.org")
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
###################### 载入包 #####################
require(reshape2)
require(ggplot2)
require(dplyr)

#显著性标记符号“*”
noteTrans <- function(dd,C1,C2,C3){
  P_list <- as.vector(dd)
  i <- 1
  note <- c()
  while (i <= length(P_list)){
    note[i] = ""
    if (is.na(P_list[i]) == TRUE){P_list[i] = 0}
    if (P_list[i] >0 & P_list[i] <= C1){note[i] = "*"}
    if (P_list[i] >0 & P_list[i] <= C2){note[i] = "**"}
    if (P_list[i] >0 &P_list[i] <= C3 ){note[i] = "***"}
    i = i + 1
  }
  return(note)
}
############################################# 读取数据 #############################################
tax_tb <- read.table(bp_i,sep="	",comment.char = "",fileEncoding = "UTF-8",stringsAsFactors = FALSE)
rownames(tax_tb) <- tax_tb[,1]
rownames(tax_tb) <- sapply(rownames(tax_tb),function(x) gsub("_*{.+}"," ",x,perl = TRUE)) 
colnames(tax_tb) <- as.character(unlist(tax_tb[1,]))
tax_tb <-tax_tb[,-1]

#tbl_df(tax_tb)
################################# 读取rename.pick文件,挑选样本,也可以修改样本名##########################
if(bp_rp!=""){
  pick = read.table(bp_rp)
  colnames(pick) <- c(rownames(tax_tb)[1],"re_pic")
  p_tax_tb <- tax_tb[,which(t(tax_tb)[,1] %%in%% pick[,1])]
  jp_tax_tb <- inner_join(pick,as.data.frame(t(p_tax_tb)))
  tax_tb <- as.data.frame(t(jp_tax_tb)[-1,])
  colnames(tax_tb) <- as.character(unlist(tax_tb[1,]))
  rm(p_tax_tb)
}

ntax_tb <- tax_tb[-1,]
########################################## 按丰度综合排序，默认降序 ######################################
ntax_tb <- sapply(ntax_tb, as.numeric)
rownames(ntax_tb) <- rownames(tax_tb)[-1]
Rowsum <-sapply(1:nrow(ntax_tb),function(x) sum(ntax_tb[x,]))
if(bp_di =="d"){
  ntax_tb<-ntax_tb[order(Rowsum,decreasing=TRUE),]
}else if(bp_di =="i"){
  ntax_tb<-ntax_tb[order(Rowsum,decreasing=FALSE),]
}

################ 百分化#######################
ntax_tb_pec <- ntax_tb
ntax_tb_pec <- sapply(1:ncol(ntax_tb),function(x) ntax_tb_pec[,x] <- ntax_tb[,x]/sum(ntax_tb[,x]))
colnames(ntax_tb_pec) <- colnames(ntax_tb)
##将总和为0的物种去掉
del <-unlist(sapply(1:ncol(ntax_tb_pec),function(x) if(sum(ntax_tb_pec[,x])==0) x))
if(!is.null(del)) {
  mydel <-colnames(ntax_tb_pec)[c(del)]
  ntax_tb_pec <-ntax_tb_pec[,-c(del)]
}
################ 判断是否有分组,并按照分组表格的顺序排列 ##############
if(bp_map != ""){
  mg <-read.table(bp_map,head=T,sep="	",comment.char = "",check.names = FALSE)
  mg[,2] <- factor(mg[,2],levels=as.vector(unique(mg[,2])))
  if(paint!="none"){
    sc <- read.table(paint,sep="\t",comment.char = "",check.names = FALSE)
    sc <- sc[which(sc[,1] %%in%% unique(mg[,2])),]
    mg[,2] <- factor(mg[,2],levels=as.vector(sc[,1]))
  }
  rownames(mg) <- mg[,1];colnames(mg)[1] <- c("sampleID")
  ntax_tb_pec <- ntax_tb_pec[,rownames(mg)]
  mg1 <- lapply(as.vector(unique(mg[,2])),function(x) ntax_tb_pec[,which(as.vector(mg[,2]) %%in%% x)])
  if (bp_ave == 1){
    avg <- matrix(NA,nrow=nrow(ntax_tb_pec),ncol = nlevels(mg[,2]))
    for (i in 1:nlevels(mg[,2])){
      for (j in 1:nrow(ntax_tb_pec)){
        avg[j,i] <- mean(mg1[[i]][j,])
      }
    }
    rownames(avg) <- rownames(ntax_tb_pec);colnames(avg) <- unique(mg[,2])
    ntax_tb_pec <- avg
  }else if (bp_ave == 2){
    mg2 <- lapply(1:nlevels(mg[,2]),function(i) mg1[[i]][,order(mg1[[i]][1,],decreasing = TRUE)])
    ntax_tb_pec <- as.matrix(as.data.frame(mg2))
    #colnames(ntax_tb_pec) <- gsub("^X","",colnames(ntax_tb_pec))
    colnames(ntax_tb_pec) <- colnames(ntax_tb_pec)
  }
}
####################################### 物种丰度挑选，有四种情况n'm'p‘s ###############################
if(bp_n>0){##是否挑选每个样本丰度前n个的物种进行画图
  max_names <- sapply(1:ncol(ntax_tb_pec),function(x) rownames(ntax_tb_pec)[order(ntax_tb_pec[,x],decreasing=T)[1:bp_n]])
  ntax_tb_pec <- ntax_tb_pec[which(rownames(ntax_tb_pec) %%in%% unique(as.vector(max_names))),]
}else if(bp_m>0){##是否挑选丰度综合排在前m个的物种进行画图
  ntax_tb_pec <- ntax_tb_pec[order(Rowsum,decreasing=TRUE),]
  ntax_tb_pec <- ntax_tb_pec[1:bp_m,] 
}else if(bp_p >0){##挑选丰度综合排在百分比为p的物种进行画图,p以下的为others，可以选择保留或删除
  if(bp_oth=="T"){
    minp <- sapply(1:nrow(ntax_tb_pec),function(y) all(ntax_tb_pec[y,]<=bp_p)) 	 
    if (table(minp)[2] ==1){
      ntax_tb_xp <- t(as.matrix(ntax_tb_pec[minp,]))
    }else{
      ntax_tb_xp <- ntax_tb_pec[minp,]
    }
    other <- sapply(1:ncol(ntax_tb_xp),function(y) sum(ntax_tb_xp[,y]))
    ntax_tb_pec <-rbind(ntax_tb_pec[!minp,],other)
    rownames(ntax_tb_pec)[nrow(ntax_tb_pec)] <-"Others"
  }else{
    minp <-sapply(1:nrow(ntax_tb_pec),function(y) all(ntax_tb_pec[y,]<=bp_p))          
    if (table(minp)[2] ==1){
      ntax_tb_xp <- t(as.matrix(ntax_tb_pec[minp,]))
    }else{
      ntax_tb_xp <- ntax_tb_pec[minp,]
    }
    other <- sapply(1:ncol(ntax_tb_xp),function(y) sum(ntax_tb_xp[,y]))
    ntax_tb_pec <-rbind(ntax_tb_pec[!minp,])
  }
}
#是否根据列表挑选物种
if(bb_i!="none"){
  std <- read.table(bb_i);rownames(std)<-std[,1]
  ntax_tb_pec <- ntax_tb_pec[rownames(std),]
}

############################################## 输出表格 ##############################
data <- ntax_tb_pec
data_temp<-cbind(row.names(data),data)
colnames(data_temp)<-c("Taxon",colnames(data))
write.table(data_temp,row.names=FALSE,paste0("percent.",bp_i),sep="	",eol="\n",quote=FALSE)

##文件名设置
ggbname <- unlist(strsplit(bp_i,split=".",fixed=TRUE))[1]
ggbn <- paste0(ggbname,"_barplot.pdf")
ggbu <- paste0(ggbname,"_bubble.pdf")
if(bp_map != ""){
  gbb <- paste0(as.vector(unique(mg$group)),collapse = "-")
  if(length(unique(mg$group))==2){gbbname <- paste0(ggbname,"_wilcox_",gbb)}
  if(length(unique(mg$group))>2){gbbname <- paste0(ggbname,"_kruskal_",gbb)}
  ggbn <- paste0(ggbname,"_",gbb,"_barplot.pdf")
  ggbu <- paste0(ggbname,"_",gbb,"_bubble.pdf")
}
########################################## 计算组间差异显著性 #################################
if(bp_map != ""&bp_ave==0){
  data2 <- left_join(data.frame(t(data),sampleID=colnames(data)),mg)
  rownames(data2) <- data2[,"sampleID"];data2[,"sampleID"]<-NULL
  #导出group内的分组名
  mgn<-as.vector(unique(mg[,2]))
  ###输出文件列名
  name1<-c()
  for (g in mgn){name1<-c(name1,paste(g,'Mean()',sep='-'),paste(g,'SE()',sep='-'))}
  name1<-c(name1,'P-value','sign')
  ###输出p值、mean、SE
  outdatal<-c()
  for (i in colnames(data2)[1:(ncol(data2)-1)]){
  set.seed(100)
    if (length(mgn)==2){
      d1 <- data2[which(data2$group %%in%% mgn[1]),i]
      d2 <- data2[which(data2$group %%in%% mgn[2]),i]
      pval<-wilcox.test(d1,d2,exact=F)$p.value
    }else{pval<-kruskal.test(data2[,i]~data2$group)$p.value}
    Fval<-''
    for (g in mgn){
      cdata<-data2[rownames(data2)[which(data2$group==g)],i]
      gmean<-mean(cdata);gse<-sd(cdata)/sqrt(length(cdata))
      outdatal<-c(outdatal,gmean,gse)
    }
    outdatal<-c(outdatal,pval,Fval)
  }
  outdatal <- as.numeric(outdatal)
  outdata<-matrix(outdatal,ncol = length(name1),byrow = T)     #转为矩阵
  rownames(outdata)<-colnames(data2)[-ncol(data2)];colnames(outdata)<-name1  #加行名列名
  outdata[,"sign"] <- noteTrans(outdata[,"P-value"],0.05,0.01,0.001)
  ## 输出数据
  write.table(outdata,paste0(gbbname,".xls"),col.names=NA,row.names=T,sep = '	')
  #print (c('output:',paste0(gbbname,".xls")))
  ################################### p值筛选###################################
  outdata<-format(outdata, scientific = FALSE)
  if(bb_i=="none"){#只显示差异显著的物种
    outdata <- outdata[which(outdata[,"P-value"]<=0.05),-ncol(outdata),drop=F]
    data2 <- data2[,c(rownames(outdata),"group"),drop=F]
  }
  
  SIGN <- noteTrans(outdata[,"P-value"],0.05,0.01,0.001)
  names(SIGN)<-rownames(outdata)
  
  bbcol <-c("#CD3333","#483D8B","#458B00","#EEC900","#BFBC3B","#6E8B3D","#00688B","#C10077","#CAAA76","#EEEE00","#458B00","#8B4513","#008B8B","#6E8B3D","#8B7D6B","#7FFF00","#CDBA96","#ADFF2F")
  if(paint!="none"){bbcol <- as.vector(sc[,2])}
  bbtheme <- theme_bw() + ##设置主题
    theme(plot.title=element_text(size=rel(1),hjust=0.5),
          plot.margin = unit(c(3, 3, 3, 3), "lines"),
          axis.title=element_text(size=rel(1)),
          axis.text.x=element_text(angle=45,hjust =1),
          legend.title=element_blank(),
          legend.text = element_text(size = 6),
          legend.key.size = unit(.4,'cm')
          #legend.position = "top"
    )
  ################################# 条形图 ##########################################
  if (pt==1){
    #重铸文件格式
    doutdata <- as.data.frame(outdata)
    if (nrow(doutdata)==1){
      outdata <- as.data.frame(matrix(apply(doutdata,2,as.numeric),nrow=1,dimnames=list(rownames(doutdata),colnames(doutdata))))
    } else {
      outdata <- as.data.frame(apply(apply(doutdata,2,as.character),2,as.numeric));dimnames(outdata) <- dimnames(doutdata)
    }
    sums<-rowSums(outdata[,1:(ncol(outdata)-2)]); outdata<-cbind(outdata,sums)
    outdata<-outdata[order(outdata[,'sums'],decreasing =T ),]
    SIGN <- SIGN[rownames(outdata)]
    
    pd<-c()
    for (i in rownames(outdata)){
      for (g in mgn){pd<-c(pd,i,outdata[i,paste(g,'Mean()',sep='-')],outdata[i,paste(g,'SE()',sep='-')],g)}
    }
    ##画图数据准备
    pdata<-as.data.frame(matrix(pd,ncol=4,byrow = T))
    colnames(pdata)<-c('x','mean','se','gp')
    pdata$x    <- factor(pdata$x,levels=unique(pdata$x))
    pdata$gp   <- factor(pdata$gp,levels=unique(pdata$gp))
    pdata$mean <- as.numeric(as.character(pdata$mean))
    pdata$se   <- as.numeric(as.character(pdata$se))
    ycor <- max(as.numeric(as.vector(rowSums(pdata[,2:3]))))+0.01
    
    bbar <- ggplot(data=pdata,aes(x=x,y=mean,fill=gp))+
      geom_errorbar(aes(ymax = mean + se, ymin = mean -  se),
                    position = position_dodge(0.9), width = 0.15)+
      geom_bar(stat = "identity", position = "dodge")+
      scale_fill_manual(values=bbcol)+
      xlab('')+
      #coord_trans(y = "log10")+
      #scale_y_log10() + annotation_logticks(sides = "l") +
      scale_y_continuous(paste("Abundance(mean±SE)"), expand = c(0, 0))+
      coord_cartesian(ylim=c(0,ycor))+
      annotate("text", x=names(SIGN), y=ycor-0.01, label=SIGN)+
      ggtitle(ggbname)+
      bbtheme
    ggsave(bbar,filename=paste0(gbbname,"_barplot.pdf"),width = (log2(nrow(outdata))+6),height = 6,limitsize = FALSE)
  }
  ########################################## 箱式图 ##############################################
  if(pt==2){
    data3 <- melt(data2)
    if(Lg=="T"){
      data3$value<-log2(data3$value)
    }else{
      data3$value<-100*data3$value
    }
    ycor <- max(as.numeric(as.vector(data3$value)))+0.1
    ##设置主题
    bx <- ggplot(data3,aes(x=variable,y=value,fill=group))+
      geom_boxplot()+
      scale_fill_manual(values=bbcol)+
      #coord_cartesian(ylim=c(0,ycor))+
      ylab("log2(Relative abundance(%%))") +
      scale_x_discrete(name="")+
      annotate("text", x=names(SIGN), y=ycor-0.01, label=SIGN)+
      ggtitle(ggbname)+
      bbtheme
    if(Lg!="T"){bx<-bx+coord_cartesian(ylim=c(0,ycor))+ylab("Relative abundance(%%)")}
    ggsave(bx,filename=paste0(gbbname,"_boxplot.pdf"),width = (log2(nrow(outdata))+6),height = 6,limitsize = FALSE)
  }
}  

######################################### 设置柱形图颜色 ######################################
#mycol <- c("#5B9BD5","#FF8C00","#B0C4DE","#FFD700","#4169E1","#32CD32","#4682B4","#8B4726","#696969","#B8860B","#6495ED","#006400","#87CEEB","#FF7F50","#DCDCDC","#F0E68C","#4682B4","#98FB98","#1E90FF","#D2691E","#778899","#DAA520")
#mycol <- rep(mycol,20)
#mycol<-rev(mycol[1:nrow(data)])

mycol <- c(35, 107, 51, 144, 36, 37, 7, 12, 84, 116, 121, 77, 56, 386, 373, 423, 435,471, 438, 512, 130, 52, 47, 6, 11, 43, 54, 367, 382, 422, 4, 8, 375, 124, 448, 419, 614, 401, 403, 613, 583, 652, 628, 633, 496, 638, 655, 132, 503, 24)
mycol <-colors()[rep(mycol,20)]
mycol<-rev(mycol[1:nrow(data)])

if(pt==0){
  ######################################### ggplot2 ###########################################
  mtax_tb <- melt(data)
  mtax_tb1 <- melt(data)
  if (bp_map != ""){
    colnames(mg)[1] <- c("Var2")
    mtax_tb <- left_join(mtax_tb,mg)
    mtax_tb$Var2 <- factor(as.character(mtax_tb$Var2), levels=unique(mtax_tb$Var2))
  }
  mtax_tb$Var1 <- factor(mtax_tb$Var1, levels=rev(unique(mtax_tb$Var1)))
  

  #if(bp_map != ""){ ##如果存在分组，facet
  #  for (i in 1:nrow(mtax_tb)){
  #    mtax_tb$group[i] <- mg[which(mg$`#sampleID` == mtax_tb$Var2[i]),2]
  #  }
  #}
  
  bartheme <- theme_bw() + ##设置主题
    theme(plot.title=element_text(size=rel(1),hjust=0.5),
          plot.margin = unit(c(1, 3, 1, 1), "lines"),
          axis.title=element_text(size=rel(1)),
          axis.text.x=element_text(size=rel(1),angle = 90, vjust = 0.5, hjust = 0.5,color = "black"),
          axis.ticks.x=element_blank(),
          axis.line.y=element_line(color="black"),
          panel.grid.major=element_line(color="white"),
          panel.grid.minor=element_line(color="white"),
          panel.border=element_rect(color="white"),
          legend.title=element_blank(),
          legend.text = element_text(size = 6),
          legend.key.size = unit(.4,'cm'),
          #legend.position = "bottom"
          legend.justification = c(1,1)
    )
  
  ##设置图片大小
  wdth=dim(data)[2]
  if (wdth<4){wdth=4}
  if (wdth>=4){wdth=5.2+wdth*.2}
  if (wdth>=200){wdth=45}
  hdth=dim(data)[1]
  if (hdth<10){hdth=5}
  if (hdth>=10&hdth<=30){hdth=6}
  if (hdth>30){hdth=floor(hdth/5)}
  
  if (bb==1){
    if (bp_ave==1){##面积图与柱形图合并
      mtax_tb$vv <- as.numeric(mtax_tb$Var2)
      mtax_tb2 <- rbind(transform(mtax_tb, vv = vv-.3),
                        transform(mtax_tb, vv = vv+.3))
      p <- ggplot(mtax_tb2,aes(x=vv,y=value,fill=Var1)) + 
        geom_area(alpha=.6)+
        geom_bar(data=mtax_tb, position="stack", stat="identity",width=0.6)+
        guides(fill=guide_legend(nrow=40,byrow=F))+
        scale_fill_manual(values=mycol) + ##设置颜色
        scale_y_continuous(labels=seq(0,100,25),expand = c(0, 0)) +
        scale_x_continuous(name ="",breaks=1:max(mtax_tb$vv),labels=levels(mtax_tb$Var2))+
        ylab("Relative abundance(%%)") +
        ggtitle(ggbname) +
        bartheme ##设置主题
      ggsave(p,filename = ggbn, dpi=300, height = 8,width=wdth,limitsize = FALSE)
    }else{
      p <- ggplot(mtax_tb,aes(x=Var2,y=value,fill=Var1)) + 
        geom_bar(stat="identity", width = 0.8) +
        scale_fill_manual(values=mycol) + ##设置颜色
        guides(fill = guide_legend(reverse=TRUE,nrow=40,byrow=F)) + ##将legend翻转显示
        scale_x_discrete(name="") +
        scale_y_continuous(labels=seq(0,100,25),expand = c(0, 0)) +
        ylab("Relative abundance(%%)") +
        ggtitle(ggbname) +
        bartheme ##背景主题等设置
      if(bp_ave==2){
        p<-p+facet_grid(.~group,scale = "free",space = "free")+
          theme(strip.text=element_text(colour = "black", face="bold"), 
                strip.background=element_rect(fill="white",colour = "white"))
      }
      ggsave(p,filename = ggbn, dpi=300, height = 8,width=wdth,limitsize = FALSE)
    }
  }
  #ggplot(mtax_tb, aes(1, fill=value))+geom_bar()+coord_polar(theta = "y")
  
  if (bb==2){#气泡图
    bubble <- ggplot(mtax_tb)+
      geom_hline(aes(x=Var2,y=Var1,yintercept = 1:nrow(mtax_tb)),linetype="dashed",alpha=.2)+
      geom_vline(aes(x=Var2,y=Var1,xintercept = 1:nrow(mtax_tb)),linetype="dashed",alpha=.2)+
      geom_point(aes(x=Var2,y=Var1,size=value,color=Var1),shape=16)+
      scale_size_area(max_size=10)+
      guides(fill = guide_legend(reverse=TRUE)) + ##将legend翻转显示
      labs(color="Taxonomy",size="Relative")+
      ylab("Taxonomy")+
      xlab("Sample")+
      theme_bw() +
      theme(
        panel.grid.major.x=element_line(linetype="dashed"),      #plot.margin=margin(5,5,5,5,unit="pt"),
        axis.text.x=element_text(size=rel(1),angle = 90, vjust = 0.5, hjust = 0.5,color = "black"),
        axis.text=element_text(size=15,hjust=0.5),
        legend.background=element_rect(fill="white", colour="black"),
        legend.key=element_blank(),
        legend.key.size=unit(.6,'cm')
      )
    ggsave(bubble,filename = paste0(ggbname,"_bubble.pdf"), dpi=300, height=(hdth+8) ,width=(wdth+10),limitsize = FALSE)
  }
}
'''%(args.i,args.rp,args.di,args.n,args.m,args.p,args.oth,args.map,args.ave,args.bi,args.pt,args.bb,args.col,args.Lg))

cmd.close()

os.system('Rscript gg_barplotcmd.r')
#os.system('rm gg_barplotcmd.r')
