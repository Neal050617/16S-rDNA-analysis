# 程序使用示例
#Rscript plot_tree.R -i bray_curtis_dm.txt -m map-group.txt -c color.txt

# 清理工作环境 clean enviroment object
rm(list=ls())

library(ggplot2)
library(ggtree)
library(colorspace)
library(tidyverse)
library(tidytree)
library(ape)
library(cowplot)
library(gridExtra)
library(grid)
library(optparse)
#library(dendextend)
library(readr)

if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="bray_curtis_dm.txt",
                help="输入矩阵"),
    make_option(c("-m", "--map"), type="character", default="none",
                help="分组信息"),
    make_option(c("-c", "--color"), type="character", default="none",
                help="颜色信息"),
    make_option(c("-x", "--xlim"), type="numeric", default=50,
                help="树的宽度"),    
    make_option(c("-t", "--tree_cluster"), type="character", default="average",
                help="ward.D, ward.D2, single, complete, average (= UPGMA), mcquitty (= WPGMA), median (= WPGMC) or centroid (= UPGMC)."),
    make_option(c("-e", "--height"), type="character", default="none",
                help="颜色信息"),
    make_option(c("-w", "--width"), type="character", default="none",
                help="颜色信息")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
}

#161124确定的颜色，用于分组画图
mycol <-c("#CD3333","#483D8B","#458B00","#EEC900","#BFBC3B","#6E8B3D","#00688B","#C10077","#CAAA76","#EEEE00","#458B00","#8B4513","#008B8B","#6E8B3D","#8B7D6B","#7FFF00","#CDBA96","#ADFF2F")

##dist文件画图
dist <- read.table(file=opts$input,head=T,row.names = 1,check.names=F)
hc <- hclust(as.dist(dist),method = opts$tree_cluster)
hc_tree <- as.phylo(hc)

#test1 <- as_tibble(hc_tree)
###判断PDF高度
if (opts$height == "none"){
  pdfh <- 6
  shul=nrow(dist)
  if (shul>20 )  {pdfh <- 7}
  if (shul>30 )  {pdfh <- 8}
  if (shul>40 )  {pdfh <- 9}
  if (shul>50 )  {pdfh <- 15}
  if (shul>70 )  {pdfh <- 17}
  if (shul>90 )  {pdfh <- 18}
  if (shul>100)  {pdfh <- 20}
  if (shul>120)  {pdfh <- 23}
  if (shul>140)  {pdfh <- 28}
  if (shul>150)  {pdfh <- 31}
  if (shul>200)  {pdfh <- 34}
}else{
  pdfh <- as.numeric(opts$height)
}

pdfw <- ifelse(opts$width == "none",6,as.numeric(opts$width))

print (c('pdfh:',pdfh))
print (c('pdfw:',pdfw))

pdfname <- strsplit(opts$input,".txt")[[1]]
p <- ggtree(hc_tree, branch.length="none") + 
  geom_tiplab() + # 标出样本名
  geom_treescale() + 
  theme_tree() +
  xlim(0, opts$xlim)

###根据group修改图形、颜色
if (opts$map != "none"){# opts$color <- c("color.txt")
  gp <- read.table(file=opts$map,head=T,comment.char = "",col.names = c("label","group"))
  dist <- dist[as.vector(gp[,1]),as.vector(gp[,1])]
  hc <- hclust(as.dist(dist),method = opts$tree_cluster)
  hc_tree <- as.phylo(hc)
  sc <- tibble(group=unique(as.vector(gp[,2])),color=mycol[1:nlevels(gp[,2])])
  # 如果存在颜色文件
  if(opts$color != "none"){
    sc <- read_tsv(opts$color,col_names = c("group","color")) %>% 
      filter(group %in% unique(gp$group))
  }
  # tree matrix 
  data1 <- as_tibble(hc_tree) %>% filter(node<=nrow(dist)) %>% 
    left_join(.,gp) %>% left_join(.,sc)
  
  data2 <- as_tibble(hc_tree) %>% filter(node>nrow(dist)) %>% 
    mutate(group = c("Others"),color = "#000000")
  
  data3 <- rbind(data1,data2)
  data3$group <- factor(data3$group,levels=c(unlist(sc[,1]),"Others"))
  mycol <- c(unlist(sc[,2]),"#000000");names(mycol) <- NULL
  
  pdfname <- paste0(pdfname, ".", paste(levels(gp$group),collapse = "-"))
  p <- ggtree(hc_tree, branch.length="none") + 
    geom_treescale() + 
    theme_tree() + 
    xlim(0, opts$xlim)
  
  p <- p %<+% data3 + 
    geom_tiplab(aes(color=group, label=label),offset=0.2) + 
    theme(legend.position = 'none')+
    scale_color_manual(values=mycol) + 
    aes(color=I(group)) + 
    #xlim(0, 0.06) +
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  
  # 制作legend
  data1$group <- factor(data1$group,levels=unlist(sc[,1]))
  COL <- unlist(sc[,2]);names(COL) <- NULL
  p2 <- ggplot(data1,aes(x=label,fill=group))+
    geom_bar()+
    scale_fill_manual(values=COL) +
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  
  #导出legend的函数
  Legend <- get_legend(p2)
  P<- grid.arrange(p,Legend,nrow=1,heights=c(10),widths=c(15,2))
  #cowplot::plot_grid(p, Legend, ncol=2)
  ggsave(P,filename = paste0(pdfname,".",opts$tree_cluster,".tree.pdf"), dpi=300, height = pdfh,width= pdfw,limitsize = FALSE)
}else{
  ggsave(p,filename = paste0(pdfname,".",opts$tree_cluster,".tree.pdf"), dpi=300, height = pdfh,width= pdfw,limitsize = FALSE)
}

