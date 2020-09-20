#Rscript core_microbiome.R -i otu.genus.xls -m map-group.txt

library(ggplot2)
library(tidyverse)
library(optparse)
library(magrittr)
library(plotrix)
library(RColorBrewer)
library(reshape2)

if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="otu.genus.xls",
                help="输入文件，OTU表格或种属表格"),
    make_option(c("-m", "--map"), type="character", default="none",
                help="分组文件：map-group.txt")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
}

#定义备选颜色
#ellipse_col <- c('#6181BD4E','#F348004E','#64A10E4E','#9300264E','#464E044E','#049a0b4E','#4E0C664E'#,'#D000004E','#FF6C004E','#FF00FF4E','#c7475b4E','#00F5FF4E','#BDA5004E','#A5CFED4E','#f0301c4E','#2B8BC34E'#,'#FDA1004E','#54adf54E','#CDD7E24E','#9295C14E')
#ellipse_col <- rep(ellipse_col,20)
colors <- c('#B0C4DE',"#CD3333","#483D8B","#458B00","#EEC900","#BFBC3B","#6E8B3D","#00688B","#C10077","#CAAA76","#EEEE00","#458B00","#8B4513","#008B8B","#6E8B3D","#8B7D6B","#7FFF00","#CDBA96","#ADFF2F")

#构建作图函数（参考自 https://www.cnblogs.com/xudongliang/p/7884667.html）
flower_plot <- function(sample, otu_num, core_otu, start=90, a=0.5, b=2, r=1, ellipse_col="", circle_col="white") {# sample <- names(dd);otu_num <- dd;core_otu <- 25
  #注：参数a和b用于设置花瓣椭圆的尺寸，ellipse_col用于设置花瓣椭圆的颜色；参数r用于设置中心圆圈尺寸，circle_col用于设置中心圆圈的颜色
  par( bty = 'n', ann = F, xaxt = 'n', yaxt = 'n', mar = c(1,1,1,1))
  plot(c(0,10),c(0,10),type='n')
  n   <- length(sample)
  deg <- 360 / n
  if (ellipse_col==""){
    mycol <- alpha(colorRampPalette(brewer.pal(11,'Spectral'))(length(sample)),alpha=0.5)
    ellipse_col <- mycol[rank(-otu_num,ties.method = "first")]
  }
  res <- lapply(1:n, function(t){
    draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180), 
                 y = 5 + sin((start + deg * (t - 1)) * pi / 180), 
                 col = ellipse_col[t],
                 border = ellipse_col[t],
                 a = a, b = b, angle = deg * (t - 1))
    text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
         y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
         otu_num[t])
    
    if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
      text(x = 5 + 4 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 4 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) - start,
           adj = 1,
           cex = 1
      )
    } else {
      text(x = 5 + 4 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 4 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) + start,
           adj = 0,
           cex = 1
      )
    }
  })
  draw.circle(x = 5, y = 5, r = r, col = circle_col, border = NA)
  text(x = 5, y = 5, paste('Core:', core_otu))
}

My_theme <- theme_bw() + ##设置主题
  theme(plot.title=element_text(size=rel(1),hjust=0.5),
        plot.margin = unit(c(1, 1, 1, 1), "lines"),
        axis.title=element_text(size=25,face="bold"),
        axis.text.x=element_text(angle=45,hjust =1),
        axis.text=element_text(size=18,face="bold"),
        legend.title = element_text(""),
        legend.text = element_text(size=18,face="bold"),
        title= element_text(size=15,face= "bold", vjust=0.5, hjust=0.5))

bar_group_plot <- function(data,Name){#data <- data3
  p <- ggplot(data,aes(x=Rate,y=Count,fill=Group,ymax=Count+10)) +
    geom_bar(stat = "identity", position = "dodge")+
    scale_fill_manual(values = colors) + 
    geom_text(aes(label=Count),position=position_dodge(width = 1),
              vjust = 0.5,hjust = -0.5,size = 3) +
    labs(fill = "Groups") +
    coord_flip() + My_theme

  ggsave(Name,p,width = 8,height = 6)
}

bar_plot <- function(data,Name){#data <- data3
  p <- ggplot(data,aes(x=Rate,y=Count,ymax=Count+10)) +
    geom_bar(stat = "identity", position = "dodge")+
    geom_text(aes(label=Count),position=position_dodge(width = 1),
              vjust = 0.5,hjust = -0.5,size = 3) +
    labs(fill = "Groups") +
    coord_flip() + My_theme
  
  ggsave(Name,p,width = 8,height = 6)
}

#读入做图文件，预处理 # opts$map <- c("map-group.txt")
data <- read.table(opts$input,head= T,sep="\t",comment.char = "",
                   row.names = 1,fileEncoding = "UTF-8")
if (opts$map != "none"){
  map <- read.table(opts$map,head= T,sep="\t",comment.char = "",fileEncoding = "UTF-8")
  data <- data[,map[,1]]
  data <- data[sapply(1:nrow(data),function(x) any(data[x,]>0)),]
}
# 每个OTU 样本所占比率
otu_coverage <- apply(data,1,function(x) length(x[x>0])/ncol(data))
# 导出表格
sapply(c(0.5,0.6,0.7,0.8,0.9,1),function(z) {#z <- 0.5
  oc <- sort(otu_coverage[rownames(data[otu_coverage>=z,])],decreasing = TRUE)
  write.table(file=paste0("core.microbiome.ALL.gt.",z,".xls"),
              data.frame(name=names(oc),rate=oc),
              sep="\t",quote = FALSE,row.names = FALSE, col.names = FALSE)
  z
})

# 汇总表格
rate=sapply(c(0.5,0.6,0.7,0.8,0.9,1),function(x) length(otu_coverage[otu_coverage>=x]))
names(rate)=c('0.5','0.6','0.7','0.8','0.9','1')
rate=as.matrix(rate);colnames(rate) <- c("All")
data2 <- apply(data,2,function(s) length(s[s>0]))

# barplot
data3 <- melt(rate);colnames(data3) <- c("Rate","Group","Count")
data3$Rate <- factor(data3$Rate,levels=unlist(data3$Rate))
bar_plot(data3,"core.microbiome.ALL.barplot.pdf")

# flower_plot
for (i in 1:nrow(rate)){# i <- 1
  png(paste0('core.microbiome.flower.all.gt',rownames(rate)[i],'.png'), width = 1500, height = 1500, res = 200, units = 'px')
  flower_plot(sample = names(data2), otu_num = data2, core_otu = as.integer(rate[i,1]))
  text(x=2,10,y=,paste0('core.microbiome.flower.all.gt',rownames(rate)[i]))
  dev.off()
}

if (opts$map != "none"){
  # 分组计算核心微生物
  l1 <- vector(nlevels(map[,2]),mode="list")
  l2 <- vector(nlevels(map[,2]),mode="list")
  l3 <- vector(nlevels(map[,2]),mode="list")
  unique_map <- as.vector(unique(map[,2]))
  for (i in 1:length(unique_map)){# i  <- 1
    mm1 <- data[,map[map[,2] %in% unique_map[i],1]]
    l1[[i]] <- mm1[sapply(1:nrow(mm1),function(y) any(mm1[y,]>0)),]
    l2[[i]] <- apply(l1[[i]],1,function(z) length(z[z>0])/ncol(l1[[i]]))
    l3[[i]] <- apply(l1[[i]],2,function(q) length(q[q>0]))
    
    # 导出表格
    sapply(c(0.5,0.6,0.7,0.8,0.9,1),function(z) {#z <- 0.5
      um <- sort(l2[[i]][l2[[i]]>=z],decreasing = TRUE)
      write.table(file=paste0("core.microbiome.",unique_map[i],".gt.",z,".xls"),
                  data.frame(name=names(um),rate=um),
                  sep="\t",quote = FALSE,row.names = FALSE, col.names = FALSE)
      z
    })
  }
  
  # 整合数据
  grate=lapply(l2, function(o)
    lapply(c(0.5,0.6,0.7,0.8,0.9,1), function(p)
      length(o[o>=p])))
  grate <- as.data.frame(matrix(unlist(grate),nrow=6))
  rownames(grate) <- c(0.5,0.6,0.7,0.8,0.9,1)
  colnames(grate) <- names(l3) <- unique_map
  
  bpd <- data.frame(Rate=as.matrix(rownames(rate)),rate,as.matrix(grate))
  write.table(file="core.microbiome.xls",bpd,
              sep="\t",quote = FALSE,row.names = FALSE, col.names = TRUE)
  
  # barplot
  library(reshape2)
  bpdata <- gather(bpd, colnames(bpd)[-1], key = "Group",value = "Count") %>% 
    mutate(Group=fct_inorder(Group)) %>% 
    mutate(Rate=fct_inorder(Rate))
  bar_group_plot(bpdata,"core.microbiome.groups.barplot.pdf")
  
  # flower plot
  for (j in 1:nrow(grate)){# j <- 1
    for (k in colnames(grate)){ # k <- colnames(grate)[1]
      png(paste0('core.microbiome.flower.',k,'.gt',rownames(grate)[j],'.png'), 
          width = 1500, height = 1500, res = 200, units = 'px')
      dd <- data2[map[map[,2] %in% k,1]]
      flower_plot(sample = names(dd), otu_num = dd, core_otu = as.integer(grate[j,k]))
      text(x=2,10,y=,paste0('core.microbiome.flower.',k,'.gt',rownames(grate)[j]))
      dev.off()
    }
  }
}else{
  write.table(file="core.microbiome.xls",
              data.frame(Rate=as.matrix(rownames(rate)),rate),
              sep="\t",quote = FALSE,row.names = FALSE, col.names = TRUE)
}


