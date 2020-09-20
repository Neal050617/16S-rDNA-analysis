#Rscript unroot_tree.R -t unroot.otu_reps.raw_aligned.fasta.tre -a unrooted.otu_tree_anno.0.001.xls -s phylum
# 方法一 ggtree
library(ggplot2)
library(ggtree)
library(colorspace)
library(tidyverse)
library(ape)
library(cowplot)
library(gridExtra)
library(grid)
library(optparse)
library(RColorBrewer)
library(treeio)
if (TRUE){
  option_list <- list(
    make_option(c("-t", "--tree"), type="character", default="unroot.otu_reps.raw_aligned.fasta.tre",
                help="无根OTU聚类树"),
    make_option(c("-a", "--anno"), type="character", default="unrooted.otu_tree_anno.0.001.xls",
                help="注释信息"),
    make_option(c("-s", "--select_tax"), type="character", default="phylum",
                help="注释信息")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
}

#161124确定的颜色，用于分组画图
#mycol <- c("#CD3333","#483D8B","#458B00","#EEC900","#BFBC3B","#6E8B3D","#00688B","#C10077","#CAAA76","#EEEE00","#458B00","#8B4513","#008B8B","#6E8B3D","#8B7D6B","#7FFF00","#CDBA96","#ADFF2F")

# barplot(1:20,col = mycol)
# 读取树文件
tree <- read.newick(opts$t)

raw_tree <- readLines(opts$tree)
writeLines(gsub(");$","):1;",raw_tree),paste0("Rooted_",opts$tree))
tree <- read.tree(paste0("Rooted_",opts$tree))

# 读取树物种注释信息
tax <- read.table(opts$a,head= T ,sep="\t",comment.char = "",row.names = 1,fileEncoding = "UTF-8")
groupInfo <- split(row.names(tax), tax[,opts$s]) # OTU and phylum for group

mycol <- colorRampPalette(brewer.pal(9,'Set1'))(nlevels(tax[,opts$s]))

## 将分组信息添加到树中
tree <- groupOTU(tree, groupInfo)
#tree$node.label <- tree$node.label[tree$node.label != ""]
# 画树，按组上色
set.seed(2015-11-26)
p <- ggtree(tree, aes(color=group), layout='circular',branch.length='none',size=1) + 
  guides(fill=FALSE) +
  theme(legend.position = 'none')+
  #scale_fill_manual(values=mycol) +
  geom_tiplab(size=1.8, aes(angle=angle))
  #scale_color_manual(values=c("#000000",mycol))
  #ggtitle("unrooted Phylogenetic tree")

  if ("" %in% tree$node.label){
    p <- p + scale_color_manual(values=c("#000000",mycol)) 
  }else{
    p <- p + scale_color_manual(values=mycol) 
  }

# 制作legend
Bp <- as_tibble(tree) %>% filter(label %in% row.names(tax)) %>% filter(group != "")
Bp$group <- droplevels(Bp$group)
# as_tibble(tree)  %>% print(n = 52, width = Inf)
p2 <- ggplot(Bp,aes(x=label,fill=group))+
  geom_bar()+
  scale_fill_manual(values=mycol,name=opts$s)

#导出legend的函数
Legend <- get_legend(p2)
P<- grid.arrange(p,Legend,nrow=1,heights=c(10),widths=c(10,2))
ggsave(P,filename = paste0(opts$s,"-",opts$t,".pdf"), dpi=300, height = 15,width= 15,limitsize = FALSE)  

# p <- ggtree(tree, aes(color=group), ladderize=F,layout="daylight", size=1) +
#   guides(fill=FALSE) +
#   scale_fill_manual(values=color) +
#   ggtitle("unrooted Phylogenetic tree") +
#   # geom_tiplab(size=3) +
#   geom_tiplab2(aes(fill = group),
#                color = "white", # color for label font
#                geom = "label",  # labels not text
#                alpha = 0.6
#   )



################################################
## 方法二 ape
#library(ape)
#library(phytools)
#tree2 <- read.tree(file = "unroot.otu_reps.raw_aligned.fasta.tre")
#tax <- read.table("otu.unroot_tree_anno.xls",head= T ,sep="\t",comment.char = "",row.names = 1#,fileEncoding = "UTF-8")
#
### color
## OTU名
#tipcol <- rep('black', length(tree2$tip.label))
## make a vector of color we want:
#colorsList <-c("red", "darkolivegreen3", "blue", "orange", "yellow", "purple", "pink", "green")
## replace colours where grep gives "TK" as red, etc in a loop
#for(i in 1:nlevels(tax$phylum)){ # i <- 1
#  tipcol[tax$phylum %in% levels(tax$phylum)[i]] <- colorsList[i]
#}
## 线
#edgecol <- rep('black', nrow(tree$edge))
#num1 <- which(tree$edge[,2] <= length(tree$tip.label))
#num2 <- which(tree$edge[,2] > length(tree$tip.label))
#for (j in 1:length(num1)){ #j=1
#  edgecol[num1[j]] <- tipcol[j]
#}
#
## 画树
#plot(tree, #"u",
#     use.edge.length = FALSE,
#     tip.color=tipcol, 
#     edge.color = edgecol,
#     cex = 0.9)
#nodelabels()
#tiplabels()
