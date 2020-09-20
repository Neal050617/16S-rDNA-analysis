# Rscript Specaccum.R -i rarefac.otu_table.xls -m map-group.txt
library(tidyverse)
library(magrittr)
library(optparse)
library(vegan)
options("endocing"="UTF-8")
#install.packages("Matrix", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")#修改镜像源

if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="rarefac.otu_table.xls",
                help="丰度表格"),
    make_option(c("-m", "--map"), type="character", default="none",
                help="分组表格：map-group.txt"),
    make_option(c("-c", "--color"), type="character", default="none",
                help="颜色设置"),
    make_option(c("-s", "--shape"), type="character", default="none",
                help="形状设置"),
    make_option(c("-p", "--per"), type="integer", default=100,
                help="permutations"),
    make_option(c("-t", "--CI"), type="double", default=0.95,
                help="置信区间")
    )
  opts <- parse_args(OptionParser(option_list=option_list))
}

Specaccum <- function (comm, permutations = 100, CI = opts$CI) {
  x <- as.matrix(comm) %>% .[, colSums(.) > 0, drop = FALSE]
  n <- nrow(x) # OTU个数
  p <- ncol(x) # 样本数
  if (p == 1) {
    x <- t(x)
    n <- nrow(x) 
    p <- ncol(x) 
  }
  accumulator <- function(x, ind) {
    rowSums(apply(x[ind, ], 2, cumsum) > 0)
  }
  set.seed(20190517)
  permat <- vegan:::getPermuteMatrix(permutations, n)
  perm <- apply(permat, 1, accumulator, x = x)
  specaccum <- apply(perm, 1, mean)
  sdaccum <- apply(perm, 1, sd)
  seaccum <- sdaccum/sqrt(n)
  ciaccum <- seaccum * (qt(opts$CI/2 + .5, n-1))
  out <- list(sites =1:n, richness = specaccum, sd =  sdaccum, se = seaccum, ci = ciaccum) %>% bind_cols
  out2 <- perm %>% as_tibble %>% rowid_to_column %>% 
    reshape2::melt(id.vars = 'rowid') %>% 
    select(-variable) %>% as_tibble
  return(list(out,out2))
}

mytheme1 <- theme_bw() + 
  theme(panel.grid.major=element_line(color="white"),
        panel.grid.minor=element_line(color="white"),
        plot.title = element_text(hjust = 0.5),
        legend.title=element_blank(),
        legend.background=element_blank(),
        legend.key=element_blank(),
        axis.title=element_text(size=12,face="bold"),
        axis.text=element_text(size=8,face="bold"))

mapcol <-c("#CD3333","#483D8B","#458B00","#EEC900","#BFBC3B","#6E8B3D","#00688B","#C10077","#CAAA76","#EEEE00","#458B00","#8B4513","#008B8B","#6E8B3D","#8B7D6B","#7FFF00","#CDBA96","#ADFF2F")

mypch <-c(16,15,17,18,7,8,9,10,11,12,13,14,21,22,23,24,25,3,4)

data1 <- read_tsv(opts$input)

if (opts$map != "none"){
  map <- read_tsv(opts$map) %>% rename(SampleID = `#SampleID`)
  
  # color.txt
  if (opts$color != "none"){# opts$color = c("color.txt")
    sc <- read.table(opts$color,sep="\t",comment.char = "",check.names = FALSE)
    sc <- sc[which(as.vector(sc[,1]) %in% unique(unlist(map$group))),]
    map$group <- factor(map$group,levels = as.vector(sc[,1]))
  } else{
    map$group <- factor(map$group,levels = unique(unlist(map$group)))
    sc <- cbind(levels(map$group),mapcol[1:nlevels(map$group)]) %>% 
      as.data.frame(stringsAsFactors = FALSE)
  }
  sc$V2 <- factor(sc$V2,levels = as.vector(sc[,2]))
  
  #设置形状
  if (opts$shape!="none"){
    sp <- read.table(opts$shape,sep="\t",comment.char="",check.names=F,stringsAsFactors = F)
    sp <- sp[which(as.vector(sp[,1]) %in% unique(unlist(map$group))),]
    sp[,2] <- factor(sp[,2],levels=unique(sp[,2]))
    mypch <- unique(sp[,3])
  }
  
  DD <- data1 %>% .[,-1] %>% t %>% 
    as.data.frame %>% rownames_to_column %>% as_tibble %>% 
    right_join(map,by = c("rowname" = "SampleID")) %>%
    mutate(group = fct_inorder(group)) %>% select(-rowname) %>%
    group_by(group) %>% nest
  
  colnames(DD$data[[1]]) <- data1$`OTU ID`
  
  DD1 <- DD %>% mutate(data = map(data,function(x) 
    Specaccum(comm = x,permutations = opts$per)[[1]])) %>% unnest
  DD2 <- DD %>% mutate(data = map(data,function(x) 
    Specaccum(comm = x,permutations = opts$per)[[2]])) %>% unnest

  p1 <- ggplot() + 
    geom_ribbon(DD1,mapping = aes(x = sites,ymin = richness-ci, ymax=richness+ci, 
                                  group = group, fill = group), alpha = 0.2,show.legend = FALSE) +
    geom_line(DD1,mapping = aes(x = sites, y = richness, 
                                group = group, color = group),size=0.6,show.legend = FALSE) +
    geom_point(DD1,mapping = aes(x = sites, y = richness, 
                                group = group, color = group, shape = group),size=1) +
    scale_shape_manual(values=mypch,name="Group") +
    scale_colour_manual(values=as.vector(sc$V2),name="Group") +
    scale_fill_manual(values=as.vector(sc$V2),name="Group") +
    geom_errorbar(DD1,mapping = aes(x = sites, ymin=richness-ci, ymax=richness+ci, 
                                    color = group), width=.2, size = 0.1,show.legend = FALSE) +
    xlab("Number of Samples")+ylab(str_c("Number of OTUs(Mean ","\U00B1"," 95%CI",")")) +
    mytheme1 +
    theme(legend.position=c(1,0), legend.justification=c(1,0))
  
  ggsave("specaccum.groups.pdf", p1, width = 6, height = 5)
}else{
  comm <- data1 %>% .[,-1] %>% t %>% 
    as.data.frame %>% rownames_to_column %>% 
    as_tibble %>% select(-rowname) %>%
    Specaccum(.,permutations = opts$per)

  p2 <- ggplot() + 
    geom_line(comm[[1]],mapping = aes(x = sites, y = richness), size=0.6) +
    stat_boxplot(comm[[2]],mapping = aes(x = rowid, y = value, group = rowid),
                 geom ='errorbar', width = 0.6) +
    geom_boxplot(comm[[2]],mapping = aes(x = rowid, y = value, group = rowid), 
                 width = .5, color = "#2177C7", fill = "#F9703C") +
    xlab("Number of Samples")+ylab("Number of OTUs") +
    mytheme1
  
  ggsave("specaccum.all.boxplot.pdf", p2, width = 6, height = 5)
}
