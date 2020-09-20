# Rscript plor_alpha_box_dot.r -a alpha_rarefac.summary.xls -g map-group.txt
# 
library(optparse)
library(tidyverse)
library(magrittr)
library(reshape2)
library(ggplot2)
if (TRUE){
  option_list <- list(
    make_option(c("-a", "--alpha"), type="character", default="alpha_rarefac.summary.xls",
                help="输入的OTU表格"),
    make_option(c("-g", "--map"), type="character", default="map-group.txt",
                help="分组文件"),
    make_option(c("-c", "--color"), type="character", default="none",
                help="指定分组颜色:color.txt")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
}

mycol <-c("#CD3333","#483D8B","#458B00","#EEC900","#BFBC3B","#6E8B3D","#00688B","#C10077","#CAAA76","#EEEE00","#458B00","#8B4513","#008B8B","#6E8B3D","#8B7D6B","#7FFF00","#CDBA96","#ADFF2F")

a <- read.table(opts$alpha,header=T,sep="\t")
b <- read.table(opts$map,header=T,sep="\t",comment.char = "")
b[,2] <- factor(b[,2],levels=unique(b[,2]))
b$color <- mycol[as.numeric(b[,2])]

#如果有颜色定义的话，就可以直接定义了
if(opts$color != "none"){
  sc <- read.table(opts$color,sep="\t",comment.char = "",check.names = FALSE)
  sc <- sc[which(sc[,1] %in% unique(b[,2])),]
  mycol <- as.vector(sc[,2])
  b[,2] <- factor(b[,2],levels=as.vector(sc[,1]))
  b[,3] <- mycol[as.numeric(b[,2])]
}

colnames(b)[1] <- colnames(a)[1]
a <- left_join(a,b)

#添加表头
name1<-c()
for (g in as.vector(unique(b$group))){
  name1<-c(name1,paste(g,'Mean',sep='-'),paste(g,'SE',sep='-'))
}
name1<-c(name1,'P-value')

GGBOX <- function(data,jitter="TRUE",mode="box"){# data = a[,c(9,10,11)]
  if (mode == "box"){
    p <- ggplot(melt(data),aes(x=group,y=value)) + 
      stat_boxplot(geom ='errorbar', width = 0.6) +
      geom_boxplot(notch = F,aes(color = group))
    if (jitter == "TRUE"){p <- p + geom_point(aes(color = group),position="jitter")}
    
  }else if (mode == "dot"){
    p <- ggplot(melt(data),aes(x=group,y=value,fill=group)) + 
      geom_dotplot(binaxis='y', stackdir='center',dotsize = 0.6) +
      stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                   geom="errorbar", color="black", width=0.2, show.legend = FALSE) +
      stat_summary(fun.y="mean", geom="point", color="black", show.legend = FALSE)
    }

    p <- p + scale_color_manual(values = unique(data$color), guide = FALSE) + 
    scale_fill_manual(values = unique(data$color)) + 
    xlab("") +
    ylab("") +
    ggtitle(colnames(data)[1]) +
    theme_bw() + ##设置主题
    theme(plot.title=element_text(size=rel(1),hjust=0.5),
          plot.margin = unit(c(1, 1, 1, 1), "lines"),
          axis.title=element_text(size=25,face="bold"),
          axis.text.x=element_text(angle=45,hjust =1),
          axis.text=element_text(size=18,face="bold"),
          title= element_text(size=15,face= "bold", vjust=0.5, hjust=0.5))

  Width <- ifelse(log2((nlevels(data$group))^(1/2))<=1,1,log2((nlevels(data$group))^(1/2)))
  ggsave(paste0(colnames(data)[1],"-",mode,"plot.pdf"),p,width = 6*Width,height = 6)
}

for (i in c(4,5,6,7,9)){#c = c(4,5,6,7,9)[1]
  GGBOX(a[,c(as.numeric(i),10,11)])
  GGBOX(a[,c(as.numeric(i),10,11)],mode="dot")
}
