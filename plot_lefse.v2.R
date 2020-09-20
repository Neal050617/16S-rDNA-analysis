#Rscript plot_lefse.R -i lefse.data.xls -r LDA.xls -l 2 -m otu
#Rscript plot_lefse.R -i KEGG.lefse.data.tsv -r LDA.xls -l 2 -m kegg
# 方法一 ggtree
# 20190510 添加kegg画图选项；添加lm辅助修改点大小
library(ggplot2)
library(ggtree)
library(colorspace)
library(tidyverse)
library(ape)
library(optparse)
library(microbiomeViz)
library(treeio)
library(grid)
library(gridExtra)
library(cowplot)
library(magrittr)

if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="lefse.data.xls",
                help="准备好的|分隔文件;lefse.data.xls;KEGG.lefse.data.tsv"),
    make_option(c("-r", "--res"), type="character", default="LDA.xls",
                help="pyhton-lefse结果"),
    make_option(c("-g", "--group"), type="character", default="map-group.txt",
                help="分组文件"),
    make_option(c("-l", "--lda"), type="numeric", default=2,
                help="lda阈值"),
    make_option(c("-c", "--color"), type="character", default="none",
                help="color.txt"),
    make_option(c("-m", "--mode"), type="character", default="otu",
                help="otu,kegg,ko")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
}

# microbiomeViz需要 R 3.5 以上，依赖包安装
#library(devtools)
#setRepositories(ind=1:2)
#devtools::install_github("lch14forever/microbiomeViz")

## functions
# point size
LM <- function(oo,y=c(0.1,2)){#oo <- data2 %>% colSums
  oo <- log10(oo)
  x = c(min(oo),max(oo))
  y = y
  o <- lm(y ~ x)
  predict(o,data.frame(x = oo))
}

ParseMetaphlanTSV <- function(taxtab,node.size.from=0.1, node.size.to=5, mode = 'otu'){ 
  # taxtab <- dat ; mode <- 'kegg'
  names(taxtab) <- c('tax','rel_abun')
  #taxtab %<>% dplyr::slice(grep('unclassified', .[,1], invert=TRUE)) # remove unclassified taxa
  if (mode == "otu"){
    tax_chars <- c('k', 'p', 'c', 'o', 'f', 'g', 's', 't')
    sp <- '__'
  }else if (mode == "kegg"){
    tax_chars <- c('L0','L1','L2','L3')
    sp <- '_'
  }
  tax_split <- strsplit(taxtab$tax, "\\|")    ## split into different taxonomy levels
  # tip
  child <- vapply(tax_split, tail, n=1, '')
  tax_class <- do.call(rbind, str_split(child, sp, n = 2))[,1]
  # parent
  parent <- vapply(tax_split, function(x) ifelse(length(x)>1, x[length(x)-1], 'root'), '')
  
  isTip <- !child %in% parent # 标出非共有的
  index <- c()
  index[isTip] <- 1:sum(isTip)
  index[!isTip] <- (sum(isTip)+1):length(isTip)
  ## tips comes first
  mapping <- data.frame(node=index, row.names=child, isTip, taxaAbun=taxtab$rel_abun)
  edges <- cbind(mapping[parent,]$node, mapping$node)
  edges <- edges[!is.na(edges[,1]),]
  
  mapping$nodeSize <- LM(mapping$taxaAbun,y = c(node.size.from,node.size.to))
  mapping$nodeClass <- factor(tax_class, levels = rev(tax_chars))
  mapping <- mapping[order(mapping$node),]
  
  node.label = rownames(mapping)[!mapping$isTip]
  phylo <- structure(list(edge = edges,
                          node.label = node.label,
                          tip.label = rownames(mapping[mapping$isTip,]),
                          edge.length=rep(1, nrow(edges)),
                          Nnode = length(node.label)),class = "phylo")
  
  d <- mapping %>% dplyr::select(-isTip)
  treedata(phylo = phylo, data = as_tibble(d)) # treeio
}

Skeleton <- function (tree, size = 2, layout = "circular", shape = 21, fill = "white", color = "black"){#tree <- tr
  ggtree(tree, size = size, layout = layout) + 
    geom_point(aes(size = I(nodeSize)),shape = shape, fill = fill, color = color)
}

geom_cladelabel3 <- function (node, label, offset = 0, offset.text = 0, extend = 0, 
                              align = FALSE, barsize = 0.5, fontsize = 3.88, angle = 0, 
                              geom = "text", hjust = 0, color = NULL, fill = NA, family = "sans", 
                              parse = FALSE, ...) {
  mapping <- NULL
  data <- NULL
  position <- "identity"
  show.legend <- NA
  na.rm <- TRUE
  inherit.aes <- FALSE
  if (!is.null(color)) {
    if (length(color) > 2) {
      stop("color should be of length 1 or 2")
    }
    if (length(color) == 0) {
      color = NULL
    }
    else if (length(color) == 1) {
      barcolor <- color
      labelcolor <- color
    }
    else {
      barcolor <- color[1]
      labelcolor <- color[2]
    }
  }
  
  if (is.null(color)) {
    if (geom == "text") {
      layer_text = ggtree:::stat_cladeText(node = node, label = label, 
                                           offset = offset + offset.text, align = align, 
                                           size = fontsize, angle = angle, family = family, 
                                           mapping = mapping, data = data, geom = geom, 
                                           hjust = hjust, position = position, show.legend = show.legend, 
                                           inherit.aes = inherit.aes, na.rm = na.rm, parse = parse, 
                                           ...)
    }
    else {
      layer_text = ggtree:::stat_cladeText(node = node, label = label, 
                                           offset = offset + offset.text, align = align, 
                                           size = fontsize, angle = angle, fill = fill, 
                                           family = family, mapping = mapping, data = data, 
                                           geom = geom, hjust = hjust, position = position, 
                                           show.legend = show.legend, inherit.aes = inherit.aes, 
                                           na.rm = na.rm, parse = parse, ...)
    }
  }
  else {
    if (geom == "text") {
      layer_text = ggtree:::stat_cladeText(node = node, label = label, 
                                           offset = offset + offset.text, align = align, 
                                           size = fontsize, angle = angle, color = labelcolor, 
                                           family = family, mapping = mapping, data = data, 
                                           geom = geom, hjust = hjust, position = position, 
                                           show.legend = show.legend, inherit.aes = inherit.aes, 
                                           na.rm = na.rm, parse = parse, ...)
    }
    else {
      layer_text = ggtree:::stat_cladeText(node = node, label = label, 
                                           offset = offset + offset.text, align = align, 
                                           size = fontsize, angle = angle, color = labelcolor, 
                                           fill = fill, family = family, mapping = mapping, 
                                           data = data, geom = geom, hjust = hjust, position = position, 
                                           show.legend = show.legend, inherit.aes = inherit.aes, 
                                           na.rm = na.rm, parse = parse, ...)
    }
  }
  list(layer_text)
}

Clade.anno <- function (gtree, anno.data, alpha = 0.3, anno.depth = 6) {
  # gtree <- p;anno.data <- lefse_lists
  short.labs <- expand.grid(0:9,letters) %>% 
    mutate(labs=paste0(Var2,Var1)) %>% .$labs %>% c(letters,.)
  
  get_offset <- function(x) {
    (x * 0.2 + 0.2)^2
  }
  
  get_angle <- function(node) {#node <- n
    data <- gtree$data
    sp <- offspring(data, node)$node
    sp2 <- c(sp, node)
    sp.df <- data[match(sp2, data$node), ]
    mean(range(sp.df$angle))
  }
  
  colnames(anno.data) <- c("node","color")
  anno.data <- arrange(anno.data, node)
  hilight.color <- anno.data$color
  node_list <- anno.data$node
  node_ids <- (gtree$data %>% dplyr::filter(label %in% node_list) %>% 
                 arrange(label))$node
  anno <- rep("white", nrow(gtree$data))
  # clour the clade
  for (i in 1:length(node_ids)) {# i <- 1
    n <- node_ids[i]
    color <- as.vector(hilight.color)[i]
    anno[n] <- color
    mapping <- gtree$data %>% dplyr::filter(node == n)
    nodeClass <- as.numeric(mapping$nodeClass)
    offset <- get_offset(nodeClass)
    gtree <- gtree + geom_hilight(node = n, fill = color, 
                                  alpha = alpha, extend = offset)
  }
  gtree$layers <- rev(gtree$layers) # point layer up
  gtree <- gtree + geom_point2(aes(size = I(nodeSize)), fill = anno, shape = 21)
  
  short.labs.anno <- NULL
  ni <- 0
  for (i in 1:length(node_ids)) {# i <- 1
    n <- node_ids[i]
    mapping <- gtree$data %>% dplyr::filter(node == n)
    nodeClass <- as.numeric(mapping$nodeClass)
    if (nodeClass <= anno.depth) {
      ni <- ni + 1
      lab <- short.labs[ni]
      if (is.null(short.labs.anno)) {
        short.labs.anno = data.frame(lab = lab, annot = mapping$label, stringsAsFactors = F)
      }else {
        short.labs.anno = rbind(short.labs.anno, c(lab, mapping$label))
      }
    }else {
      lab <- mapping$label
    }
    offset <- get_offset(nodeClass) - 0.4
    angle <- get_angle(n) + 90
    gtree <- gtree + geom_cladelabel3(node = n, label = lab, 
                                      angle = angle, fontsize = 1.5 + sqrt(nodeClass), 
                                      offset = offset, barsize = 0, hjust = 0.5)
  }
  
  if (is.null(short.labs.anno)) {
    return(gtree)
  }
  
  anno_shapes = as.numeric(fct_inorder(short.labs.anno$lab)) + 96
  short.labs.anno.cbin <- short.labs.anno %>% 
    mutate(Cbin = str_c(lab,annot,sep = ": "),
           num = 0.2) %>% 
    left_join(anno.data,c("annot" = "node"))
  
  p2 <- ggplot(short.labs.anno.cbin)+
    geom_bar(aes(x=fct_inorder(Cbin),fill=fct_inorder(Cbin)),color="black") +
    scale_fill_manual(values=as.vector(short.labs.anno.cbin$color)) +
    guides(fill = guide_legend(nrow=ifelse(opts$mode == 'otu',30,50),byrow=F,title="Annotation:")) +
    theme(legend.title = element_text(size=15,face="bold"),
          legend.text = element_text(size=12),
          legend.position = c(0, 1),
          legend.justification = c(0,1),
          legend.key.size = unit(.6, "cm"))
  
  #导出legend的函数
  Legend <- get_legend(p2)
  grid.arrange(gtree,Legend,nrow=1,heights=c(10),widths=c(10,10))
  #plot_grid(gtree,Legend,rel_heights = c(1, 1))
  
  #gtree + 
  #  geom_point(data = short.labs.anno.cbin, 
  #             mapping = aes(x = 0, y = 0, color = fct_inorder(Cbin)), size = 0, stroke = 0) + 
  #  guides(color = guide_legend(override.aes = list(size = 3,color = "white"),
  #                              nrow=40,byrow=F)) + 
  #  theme(legend.position = c(1.2, 0.5), legend.title = element_blank())
}

# 读取python-lefse生成的数据
Map <- read_tsv(opts$group) %>% 
  rename("SampleID" = colnames(.)[1])

lefse <- read_tsv(opts$res, 
                  col_types=cols(`P_value` = col_double()),
                  skip_empty_rows = TRUE) %>% 
  drop_na() %>% dplyr::filter(LDA_value >= opts$lda) %>%
  filter(`P_value` <= 0.05) %>%
  mutate(Biomaker_names = sapply(1:nrow(.),function(x){
    if (opts$mode == 'otu'){
      str_replace_all(.[x,1],c(".p__" = "|p__",".c__" = "|c__",
                               ".o__" = "|o__",".f__" = "|f__",
                               ".g__" = "|g__",".s__" = "|s__"))
      
    }else{
      str_replace_all(.[x,1],c("___" = "_")) %>% 
        str_replace_all(.,c("__" = "_")) %>%
        str_replace_all(.,c(".L2_" = "|L2_",".L3_" = "|L3_"))
    }
  })) %>% 
  mutate(Biomaker_names = sapply(1:nrow(.),function(x) 
    tail(str_split(.[x,1],"\\|")[[1]], n=1, ''))
  ) %>%
  mutate(Groups = factor(Groups,levels = unique(Map$group)[unique(Map$group) %in% unique(Groups)]))

Map <- Map %>% filter(group %in% levels(lefse$Groups)) %>% 
  mutate(group = factor(group,levels=levels(lefse$Groups))) %>% arrange(group)
mycol <-c("#CD3333","#483D8B","#458B00","#EEC900","#BFBC3B","#6E8B3D","#00688B","#C10077","#CAAA76","#EEEE00","#458B00","#8B4513","#008B8B","#6E8B3D","#8B7D6B","#7FFF00","#CDBA96","#ADFF2F")

# color.txt
if(opts$color != "none"){# opts$color = c("color.txt")
  sc <- read.table(opts$color,sep="\t",comment.char = "",check.names = FALSE)
  sc <- sc[which(as.vector(sc[,1]) %in% levels(Map$group)),]
  lefse$Groups <- factor(lefse$Groups,levels = as.vector(sc[,1]))
}else{
  lefse$Groups <- factor(lefse$Groups,levels = levels(Map$group))
  sc <- cbind(levels(lefse$Groups),mycol[1:nlevels(lefse$Groups)]) %>% as.data.frame(stringsAsFactors = FALSE)
}
sc <- sc %>% as_tibble %>% 
  mutate(V2 = factor(sc$V2,levels = as.vector(sc[,2]))) %>% 
  filter(V1 %in% lefse$Groups) %>%
  mutate(V1 = as.vector(V1))

lefse_lists <- lefse %>% full_join(sc,c("Groups"="V1")) %>% .[,c("Biomaker_names","V2")]
######################
## LDA barplot data
lefse_bar <- lefse 
if (length(unique(sc$V1)) == 2){
  Sign <- as.numeric(lefse_bar$Groups) %>% sapply(.,function(x) ifelse(x>1,x-3,x))
  lefse_bar$LDA_value <- lefse_bar$LDA_value * Sign
}

lefse_bar %<>% arrange(Groups,desc(`LDA_value`)) %>% 
  mutate(hjust = ifelse(.$LDA_value>0, 1.3,-1.3)) %>% 
  full_join(sc,c("Groups"="V1")) %>% 
  mutate(Groups = factor(.$Groups,levels(lefse$Groups))) %>%
  mutate(Biomaker_names = fct_rev(fct_inorder(Biomaker_names)))

ymin <- ifelse(min(lefse_bar$LDA_value)>=0,0,floor(min(lefse_bar$LDA_value)))
ymax <- ifelse(max(lefse_bar$LDA_value)>=0, ceiling(max(lefse_bar$LDA_value)), 
               floor(max(lefse_bar$LDA_value)))
## plot
q <- ggplot() +
  geom_hline(yintercept = ymin:ymax, color="grey88", linetype="longdash", size=0.5) +
  geom_bar(data = lefse_bar,mapping = aes(Biomaker_names, `LDA_value`, fill = Groups),stat = "identity", color = "black")

if (length(unique(sc$V1)) == 2){
  q <- q + 
    geom_text(data = filter(lefse_bar,LDA_value > 0),
              mapping = aes(x = Biomaker_names,y = 0, label = Biomaker_names),
              colour = "black", hjust = 1.05, vjust = 0.2, size = 3) +
    geom_text(data = filter(lefse_bar,LDA_value < 0),
              mapping = aes(x = Biomaker_names,y = 0, label = Biomaker_names),
              colour = "black", hjust = -0.05, vjust = 0.2, size = 3) +
    scale_x_discrete(breaks = NULL) +
    theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
}else{
  q <- q + scale_y_continuous(labels = seq(ymin,ymax,1),breaks=seq(ymin,ymax,1)) +
    theme(axis.text.y = element_text(colour = "gray0",size = 15))
}

q <- q + scale_fill_manual(values=as.vector(sc$V2)) +
  coord_flip() + xlab("") +
  ylab("LDA SCORE (log 10)") +
  theme_bw() + 
  theme(panel.background = element_rect(fill = "white"),
        plot.margin=unit(c(0.5,0.5,0.5,0.5), "lines"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(colour = "gray0",size = 12),
        axis.title = element_text(colour = "gray0",size = 15, face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size=18,face="bold"),
        legend.text = element_text(size=15)
  )

Width <- ifelse((ymax+abs(ymin)) < 10,10,ymax+abs(ymin))
Height <- nrow(lefse_bar)^(1/2) + log2(nrow(lefse_bar))
name <- str_c(str_c(sc$V1,collapse = "-"),"_lda-gt-",opts$lda,"_lefse_")

ggsave(str_c(name,"lda.pdf"),q,width = Width,height = Height,limitsize = FALSE)

###################################### cladogram #####################################
if (opts$mode != 'ko'){
  # 加载lefse数据
  df <- read_tsv(opts$input,col_names = FALSE) %>% 
    .[c(TRUE,(.[-1,-1] %>% apply(.,2,as.numeric) %>% rowMeans) > 0),]
  
  if (opts$mode == 'kegg'){
    df %<>% mutate(X1 = sapply(1:nrow(.),function(x){
      str_replace_all(.[x,1],"\\(|\\)|\\'|\\-","_") %>% 
        str_replace_all(.,c("___" = "_")) %>% 
        str_replace_all(.,c("__" = "_"))
    })) %>% mutate(X1 = paste0("L0_KEGG|",X1)) %>% 
      #add_row(.,(c("L0_KEGG",replicate(ncol(df),"1000000")),.before = 2)
      rbind(c("L0_KEGG",replicate(ncol(df)-1,"1000000"))) %>% 
      .[c(1,nrow(.),2:(nrow(.)-1)),]
    df[1,1] <- "class"
  }
  ## 计算均值用于呈现结点大小
  MEAN <- df %>% .[-1,-1] %>% apply(.,2,as.numeric) %>% rowMeans
  dat <- data.frame(V1=df[-1,1], V2=MEAN, stringsAsFactors = FALSE)
  
  # 用物种和丰度生成树骨架
  tr <- ParseMetaphlanTSV(dat,node.size.from = 0.1,node.size.to = 2 , mode = opts$mode)
  p <- Skeleton(tr, size=0.5)
  # Clade
  p <- Clade.anno(p,lefse_lists, alpha = 0.3, anno.depth = ifelse(opts$mode == 'otu',6,2))
  
  ggsave(str_c(name,"cladogram.pdf"),p,width = Height * 2,height = Height,limitsize = FALSE)

}
