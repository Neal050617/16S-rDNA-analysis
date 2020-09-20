# 示例：
# Rscript /work/users/chaoliu/scripts/make_community_boxplot.r -i phylum.percents.xls -m map-group.txt -c color.txt -l T
library(ggplot2)
library(optparse)
library(tidyverse)
library(patchwork)
library(cowplot)

if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="rarefac.otu_genus.xls",
                help="blast结果"),
    make_option(c("-m", "--map"), type="character", default="map-group.txt",
                help="分组文件"),
    make_option(c("-f", "--unif"), type="logical", default=F,
                help="要不要归一化"),
    make_option(c("-s", "--select"), type="character", default="none",
                help="筛选表格:select.tsv"),
    make_option(c("-c", "--color"), type="character", default="none",
                help="颜色"),
    make_option(c("-l", "--log"), type="double", default=2,
                help="box图取log值:0表示no;log2、log10"),
    make_option(c("-p", "--pt"), type="logical", default=T,
                help="要不要单独成图")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
}
#######
bbtheme <- function(){
  theme_bw() + ##设置主题
    theme(plot.title=element_text(size=rel(1),hjust=0.5),
          plot.margin = unit(c(3, 3, 3, 3), "lines"),
          axis.title=element_text(size=rel(1)),
          axis.text.x=element_text(angle=45,hjust =1),
          legend.title=element_blank(),
          legend.text = element_text(size = 6),
          legend.key.size = unit(.4,'cm'))
}

GGbox <- function(data,mycol,Log,nm,pt=FALSE){# data <- Data3$data[[i]] 
  
  data$value <- 100*data$value
  #ycor <- max(as.numeric(as.vector(data$value)))+1
  
  bx <- ggplot(data,aes(x=group,y=value,fill=group))+
    geom_boxplot()+
    geom_point()+ 
    geom_line(aes(group=Index)) +
    scale_fill_manual(values=mycol)+
    #coord_cartesian(ylim=c(0,ycor))+
    ylab("Relative abundance(%)") +
    scale_x_discrete(name="")+
    #annotate("text", x=names(SIGN), y=ycor-0.05, label=SIGN)+
    ggtitle(nm) +
    bbtheme()
  
  if(Log!=0){
    library(scales)
    if(Log==2){
      bx <- bx + 
        scale_y_continuous(trans = log2_trans(),labels = scientific)
      #bx <- bx + scale_y_continuous(trans = 'log2')
    }else if(Log==10){
      bx <- bx + 
        scale_y_continuous(trans = log10_trans(),labels = scientific) +
        annotation_logticks(sides="l")
    }
  }
  
  if (pt == TRUE){
    ggsave(bx,filename=paste0(nm,"_boxplot.pdf"),width = 6,height = 6,limitsize = FALSE)
  }
  
  return(bx)
}

#########
mycol <-c("#CD3333","#483D8B","#458B00","#EEC900","#BFBC3B","#6E8B3D","#00688B","#C10077","#CAAA76","#EEEE00","#458B00","#8B4513","#008B8B","#6E8B3D","#8B7D6B","#7FFF00","#CDBA96","#ADFF2F")

Gp <- read_tsv(opts$map) %>% rename(SampleID = colnames(.)[1]) %>%
  mutate_if(is.character,~fct_inorder(.))

Data <- read_tsv(opts$input) %>% rename(ID = colnames(.)[1]) %>%
  gather(SampleID,value,-ID) %>%
  mutate_if(is.character,~fct_inorder(.)) %>%
  spread(ID,value) %>%
  inner_join(.,Gp[,1])

# Uniform
if (opts$unif){
  Data[,2:ncol(Data)] <- sapply(1:nrow(Data),function(x) # x <- 1
    as.numeric(Data[x,2:ncol(Data)])/sum(as.numeric(Data[x,2:ncol(Data)])))
}

################################# 读取挑选物种 #########
if (opts$select != "none") { # opts$select <- "select.txt"
  pick <- read_tsv(opts$select,col_names = FALSE)
  Data <- Data[,c("SampleID",intersect(colnames(Data),pick$X1))]
}

# color.txt
if (opts$color != "none"){# opts$color = "color.txt"
  sc <- read_tsv(opts$color,col_names = FALSE)
  sc <- sc[sc$X1 %in% levels(Gp$group),]
  Gp$group <- factor(Gp$group,levels = as.vector(sc[,1]))
} else{
  sc <- cbind(levels(Gp$group),mycol[1:nlevels(Gp$group)]) %>% 
    as.data.frame(stringsAsFactors = FALSE)
}
sc$V2 <- factor(sc$V2,levels = as.vector(sc[,2]))

Data2 <- Data %>% inner_join(Gp) %>% gather(ID,value,-group,-SampleID) %>% 
  #mutate(paired = as.numeric(SampleID)) %>%
  mutate(ID = fct_inorder(ID)) %>%
  group_by(group,ID) %>% nest %>% 
  mutate(data = map(data,function(z){# z <- Data2$data[[1]]
    z %>% mutate(Index = 1:nrow(.))
  })) %>% unnest

Data3 <- Data2 %>% group_by(ID) %>% nest

PLOT <- lapply(1:nrow(Data3),function(x){# x <- 1
  GGbox(Data3$data[[x]],mycol,opts$log,as.character(Data3$ID[x]),pt=opts$pt)
})

if (!opts$pt){
  # 把图合并哇
  Q <- PLOT[[1]]
  if (nrow(Data3) >= 2){
    Legend <- get_legend(Q)
    for (i in 2:nrow(Data3)){
      Q <- Q + theme(legend.position = "none") + 
        PLOT[[i]] + theme(legend.position = "none")
    }
    Q <- Q + Legend + plot_layout()
  }
  
  ggsave("plot.pdf", Q, dpi = 600, device = cairo_pdf, 
         width = log(nrow(Data3),2)*8, height = log(nrow(Data3),2)*7)
  ## ggsave("plot.pdf", Q, dpi = 600, device = cairo_pdf, 
  ##        width = 45, height = 45)
}


## bx <- ggplot(Data2,aes(x=interaction(group,ID),y=value)) +
##   geom_boxplot(aes(fill=group)) +
##   geom_line(aes(group = interaction(Index,ID))) +
##   geom_point() + #aes(fill=group),size=2,shape=21
##   theme(legend.position = "none")
## if(opts$log!=0){
##   library(scales)
##   if(opts$log==2){
##     bx <- bx + 
##       scale_y_continuous(trans = log2_trans(),labels = scientific)
##     #bx <- bx + scale_y_continuous(trans = 'log2')
##   }else if(opts$log==10){
##     bx <- bx + 
##       scale_y_continuous(trans = log10_trans(),labels = scientific) +
##       annotation_logticks(sides="l")
##   }
## }
## 
## bx <- bx + facet_grid(SampleID~ID)
## 
