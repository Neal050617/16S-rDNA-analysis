# Rscript plot_matrix.R -i bray_curtis_dm.txt -g map-group.txt
# 基于矩阵作图，盒状图和热图
rm(list=ls())
library(magrittr)
library(tidyverse)
library(vegan)
library(optparse)
library(reshape2)
library(pheatmap)
if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="bray_curtis_dm.txt",
                help="输入矩阵"),
    make_option(c("-g", "--group"), type="character", default="none",
                help="分组文件:map-group.txt"),
    make_option(c("-c", "--color"), type="character", default="none",
                help="指定颜色color.txt"),
    make_option(c("-t", "--tree"), type="logical", default=TRUE,
                help="热图是否画聚类树"),
    make_option(c("-j", "--jitter"), type="logical", default=TRUE,
                help="画小点点")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
  print(paste("The input file is ", opts$input,  sep = ""))
  print(paste("The group file is ", opts$group,  sep = ""))
  print(paste("The color file is ", opts$color,  sep = ""))
  print(paste("plot tree or not is ", opts$color,  sep = ""))
}

mycol <-c("#CD3333","#483D8B","#458B00","#EEC900","#BFBC3B","#6E8B3D","#00688B","#C10077","#CAAA76","#EEEE00","#458B00","#8B4513","#008B8B","#6E8B3D","#8B7D6B","#7FFF00","#CDBA96","#ADFF2F")

GGBOX <- function(data,jitter=opts$jitter,mode="box",name){# data = MPL;name = nm
  if (mode == "box"){
    p <- ggplot(data,aes(x=L1,y=value)) + 
      stat_boxplot(geom ='errorbar', width = 0.6) +
      geom_boxplot(aes(color=L1),notch = F)
    if (opts$jitter == TRUE){p <- p + geom_point(aes(color = L1),position="jitter")}
    
  }else if (mode == "dot"){
    p <- ggplot(data,aes(x=L1,y=value,fill=L1)) + 
      geom_dotplot(binaxis='y', stackdir='center',dotsize = 0.6) +
      stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                   geom="errorbar", color="black", width=0.2, show.legend = FALSE) +
      stat_summary(fun.y="mean", geom="point", color="black", show.legend = FALSE)
  }
  
  p <- p + scale_color_manual(values = mycol, guide = FALSE) + 
    scale_fill_manual(values = mycol) + 
    xlab("") +
    ylab("") +
    #ggtitle(colnames(data)[1]) +
    theme_bw() + ##设置主题
    theme(plot.title=element_text(size=rel(1),hjust=0.5),
          plot.margin = unit(c(1, 1, 1, 1), "lines"),
          axis.title=element_text(size=25,face="bold"),
          axis.text.x=element_text(angle=45,hjust =1),
          axis.text=element_text(size=18,face="bold"),
          legend.title = element_text(""),
          title= element_text(size=15,face= "bold", vjust=0.5, hjust=0.5))
  Width <- ifelse(log2((nlevels(data$L1))^(1/2))<=1,1,log2((nlevels(data$L1))^(1/2)))
  ggsave(paste0(name,".",mode,"plot.pdf"),p,width = 6*Width,height = 6)
}

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

# read file
Matrix <- read.table(opts$input,header = T,row.names =1,sep = "\t",stringsAsFactors = F)

#################################### 大小定义 #############################
#根据样本数量判断热图col字体大小
samplenum=ncol(Matrix)
fsc='10'
if (samplenum>=20){fsc='9'}
if (samplenum>=40){fsc='8'}
if (samplenum>=45){fsc='7.5'}
if (samplenum>=50){fsc='7'}
if (samplenum>=60){fsc='6.5'}
if (samplenum>=80){fsc='6'}
if (samplenum>=100){fsc='5'}
if (samplenum>=160){fsc='9'}
if (samplenum>=180){fsc='3'}
if (samplenum>=200){fsc='2'}
if (samplenum>=300){fsc='1'}
fsc <- as.numeric(fsc)

##判断热图每个方块的长宽
cellhw=6
if (nrow(Matrix)<=200){cellhw=8}
if (nrow(Matrix)<=100){cellhw=10}
if (nrow(Matrix)<=50){cellhw=12}
if (nrow(Matrix)<=10){cellhw=20}

nm <- paste0(strsplit(opts$input,".txt")[[1]])

if (opts$group != "none"){
  Map <- read.table(opts$group,header = T,sep = "\t",comment.char = "",stringsAsFactors = F)
  Matrix %<>% .[Map[,1],Map[,1]]
  
  # 转换成长数据
  mpl <- unique(Map$group)
  MPL <- vector(length(mpl),mode="list")
  for (i in 1:length(mpl)){#i <- 1
    mpll <- as.vector(Map[Map[,2] %in% mpl[i],1])
    MPL[[i]] <- as.vector((as.dist(Matrix[mpll,mpll])))
  }
  names(MPL) <- mpl
  MPL <- melt(MPL)
  
  #如果有color.txt文件
  if(opts$color != "none"){
    sc <- read.table(opts$color,sep="\t",comment.char = "",check.names = FALSE)
    sc <- sc[as.vector(sc[,1]) %in% unique(Map$group),]
    mycol <- as.vector(sc$V2)
    Map$group <- factor(Map$group,levels=as.vector(sc[,1]))
    MPL$L1 <- factor(MPL$L1,levels=as.vector(sc[,1]))
  }else{
    MPL$L1 <- factor(MPL$L1,levels=as.vector(unique(MPL$L1)))
    Map$group <- factor(Map$group,levels=unique(Map$group))
  }
  
  # 热图数据对应表
  mapping <- as.data.frame(Map$group)
  rownames(mapping) <- Map[,1]
  colnames(mapping) <- colnames(Map)[2]
  # 热图颜色
  anncol <- list(group=mycol[1:nlevels(Map[,2])])
  names(anncol$group) <- unique(Map[,2])
  # 文件名
  nm <- paste0(strsplit(opts$input,".txt")[[1]], ".", paste(mpl,collapse = "-"))
  # 盒状图
  GGBOX(MPL,name=nm,jitter=opts$jitter)
  if (opts$jitter){
    GGBOX(MPL,mode="dot",name=nm)
  }
  
  # 热图
  pheatmap(Matrix,#导入的数据
           scale = "none",
           cluster_row = opts$tree,#行是否聚类
           cluster_col = opts$tree,#列是否聚类
           main = "Heatmap",
           fontsize_row = fsc,
           fontsize_col = fsc,
           cellheight = cellhw,
           cellwidth = cellhw,
           filename = paste0(nm,".heatmap.pdf"),
           annotation_col = mapping,
           annotation_row = mapping,
           annotation_colors = anncol,
           annotation_names_row = F,
           annotation_names_col = F,
           clustering_method = "complete")
}else{
  pheatmap(Matrix,
           scale = "none",
           cluster_row = opts$tree,#行是否聚类
           cluster_col = opts$tree,#列是否聚类
           main = "Heatmap",
           fontsize_row = fsc,
           fontsize_col = fsc,
           cellheight = cellhw,
           cellwidth = cellhw,
           filename = paste0(nm,".heatmap.pdf"),
           clustering_method = "complete")
}
