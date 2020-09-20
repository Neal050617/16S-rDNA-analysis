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
                help="指定分组颜色:color.txt"),
    make_option(c("-t", "--test"), type="character", default="rank",
                help="差异物种计算方法rank、t、aov"),
    make_option(c("-r", "--paired"), type="logical", default=FALSE,
                help="是否成对检验"),
    make_option(c("-m", "--mode"), type="character", default="box",
                help="箱式图box-点状密度图dot"),
    make_option(c("-b", "--base"), type="logical", default=FALSE,
                help="是否单独画各图")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
}

# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
Base_Statistics <- function(TT,ts=FALSE,ss="rank"){# TT <- mp3
  tt <- nlevels(TT$group)
  TT %<>% as.data.frame
  set.seed(20190306)
  if (ss == "rank" & tt == 2)
    TEST <- wilcox.test(value ~ group,data = TT,paired=ts,exact = TRUE)
  if (ss == "rank" & tt > 2)
    TEST <- kruskal.test(value ~ group,data = TT)
  if (ss == "t" & tt == 2)
    TEST <- t.test(value ~ group,as.data.frame(TT),paired=ts,
                   var.equal = FALSE,alternative = "two.sided")
  if (ss == "aov")
    TEST <- aov(value ~ group,TT) %>% broom::glance(.)
  return(TEST$p.value)
}

pp <- c(0.25, 0.5, 0.75)
p_names <- map_chr(pp, ~paste0("IQR",.*100, "%"))
p_funs <- map(pp, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
  set_names(nm = p_names)

Minus <- function(x,n){# 四舍五入
  d1 <- 10^(-n)
  ifelse(x >= d1,as.vector(round(x,n)),paste0("<",d1))
}
# 整理基础统计表格
Pvalue <- function(f,ts=FALSE,ss="rank"){ # f <- data2$data[[1]]
  gp1 <- levels(f$group)
  mpn <- length(gp1)
  mp2 <- vector(mpn-1,mode = "list")
  for (i in 2:mpn){# i <- 2
    gp2 <- combn(gp1,i)
    mp2[[i-1]] <- vector(ncol(gp2),mode="list")
    for (j in 1:ncol(gp2)){# j <- 1
      # 数据筛选阶段
      mp3 <- f %>% filter(group %in% gp2[,j]) %>% # 筛选出group文件
        mutate(group = factor(group,levels=gp1[gp1 %in% unique(.$group)]))
      gp <- levels(mp3$group)
      # 计算p-value;合并表格
      mp2[[i-1]][[j]] <- mp3 %>% Base_Statistics(.,ts,ss) %>% 
        enframe(value = "p_value") %>% 
        mutate(Sig_mark = ifelse(p_value > 0.05,"",
                                 ifelse(p_value > 0.01,"*",
                                        ifelse(p_value > 0.001,"**","***")))) %>%
        mutate(name = str_c(gp,collapse = "-"))
      if (mpn > 2)
        mp2[[i-1]][[j]]$Sig_mark %<>% 
        if_else(mp2[[i-1]][[j]]$name == str_c(gp1,collapse = "-"),
                str_replace_all(.,"*","#"),.)
    }
  }
  mp2 <- mp2 %>% melt(.) %>% as_tibble
  return(mp2)
}
# 整理基础统计表格
stat_merge <- function(s){ # s <- stat_test$data[[1]]
  S <- str_c(round(s$Mean,4),"(",round(s$SE,4),")")
  names(S) <- as.vector(s$group)
  as_tibble(t(as.matrix(S)))
}

Ptable <- function(t,value){ # t <- data2$data[[1]]
  outt <- t %>% select(.,name,value) %>% select(-name) %>% t %>% 
    as_tibble(.name_repair = "minimal")
  colnames(outt) <- t$name
  return(outt)
}

Sig_table <- function(t){ # t <- Mp2$data[[1]]
  outt <- t %>% select(.,name,Sig_mark) %>% select(-name) %>% t %>% 
    as_tibble(.name_repair = "minimal")
  colnames(outt) <- t$name
  return(outt)
}

Sig_mark <- function(t,maxy,map=mp1,value="p-value",v=1){# value <- "p-value";v <- 0.05;map<-mp1
  #t <- data2$data[[1]];maxy <- max(data1[,as.vector(data2$Alpha)[[4]]])
  tt <- t %>% filter(.[[value]] <= v)
  if (nrow(tt) != 0){
    tt <- tt %>% 
      mutate(L1 = as.character(L1)) %>%
      mutate(L1 = fct_inorder(L1)) %>% 
      mutate(L1 = as.numeric(L1)) %>% 
      mutate(x_lab = sapply(1:nrow(.),function(m){ # m <- 1
        mm <- strsplit(.$name[m],split = "-")[[1]]
        M <- factor(mm,levels = levels(map$group)) %>% as.numeric %>% as.character
        str_c(M[1],M[length(M)],sep="-")})) %>% # m <- 1
      mutate(x_from = sapply(1:nrow(.),function(m) 
        as.numeric(str_split(.$x_lab[m],"-")[[1]][1]))) %>%
      mutate(x_to = sapply(1:nrow(.),function(m) 
        as.numeric(str_split(.$x_lab[m],"-")[[1]][2]))) %>%
      mutate(x_from_g = levels(map$group)[x_from]) %>%
      mutate(x_to_g = levels(map$group)[x_to]) %>%
      mutate(x_med = (x_from + x_to)/2) %>%
      mutate(y_lab = seq(from = ifelse(maxy>100,maxy+5,ifelse(maxy>1,maxy+0.5,maxy+0.05)), 
                         by = ifelse(maxy>100,20,ifelse(maxy>1,0.5,0.05)),length.out = nrow(.))) %>%
      mutate(y_lab2 = seq(from = ifelse(maxy>100,maxy+15,ifelse(maxy>1,maxy+0.75,maxy+0.075)),
                          by = ifelse(maxy>100,20,ifelse(maxy>1,0.5,0.05)),length.out = nrow(.) ))
  }
  return(tt)
}

GGBOX <- function(data,jitter=TRUE,mode="box",Log=FALSE,
                  ylab,anno,anno2,tp="base",pt=FALSE){# y <- 2
  # data <- box$data[[y]];ylab <- as.vector(box$Alpha[y]);anno <- ddata[[9]]$data[[y]]
  # jitter <- TRUE;mode <- "box";Log <- FALSE;anno2 <- ddata[[4]][y,-1];
  if (mode == "box"){
    q <- ggplot() + 
      stat_boxplot(data,mapping=aes(x=group,y=value,group=group,colour=group),
                   geom ='errorbar', width = 0.6) +
      geom_boxplot(data,mapping=aes(x=group,y=value,group=group,colour=group),notch = F)
    if (jitter == TRUE){
      q <- q + geom_point(data,mapping=aes(x=group,y=value,colour=group),position="jitter")
    }
  }else if (mode == "dot"){
    q <- ggplot() + 
      geom_dotplot(data,mapping=aes(x=group,y=value,group=group,colour=group),
                   binaxis='y', stackdir='center',dotsize = 0.6) +
      stat_summary(data,mapping=aes(x=group,y=value,group=group,colour=group),
                   fun.data="mean_sdl", fun.args = list(mult=1), 
                   geom="errorbar", color="black", width=0.2, show.legend = FALSE) +
      stat_summary(data,mapping=aes(x=group,y=value,group=group,colour=group),
                   fun.y="mean", geom="point", color="black", show.legend = FALSE)
  }
  if (Log == TRUE)
    q <- q + scale_y_log10() + annotation_logticks(sides = "l")
  #q <- annotate("segment", x = 1, xend = 2, y = 15, yend = 15) +
  if (!is.null(anno)){
    q <- Anno("segment",anno,q)
    q <- Anno("text",a=anno,b=anno2,q)
    q <- q + ylim(0,max(anno$y_lab2) + ifelse(max(data$value)>100,1,
                                              ifelse(max(data$value)>1,0.5,0.2)))
    q <- q + geom_point(data=anno,mapping=aes(x=x_med,y=y_lab),colour = "grey40",
                        shape=17, show.legend = FALSE)
  }
  
  q <- q + scale_color_manual(values = as.vector(sc$V2), guide = FALSE) + 
    scale_fill_manual(values = as.vector(sc$V2), guide = FALSE) + 
    xlab("") +
    ylab(ylab) +
    #ggtitle(colnames(data)[1]) +
    theme_bw() + ##设置主题
    theme(plot.title=element_text(size=rel(1),hjust=0.5),
          plot.margin = unit(c(1, 1, 1, 1), "lines"),
          legend.text = element_text(size=10,face="bold"),
          axis.title=element_text(size=25,face="bold"),
          axis.text.x=element_text(angle=45,hjust =1),
          axis.text=element_text(size=15,face="bold"),
          title= element_text(size=15,face= "bold", vjust=0.5, hjust=0.5))
  
  Width <- ifelse(log2((nlevels(data$group))^(1/2))<=1,1,log2((nlevels(data$group))^(1/2)))
  if (pt == TRUE){
    ggsave(paste0(tp,".",ylab,"-",mode,"plot.pdf"),q,width = 6*Width,height = 6)
  }
  return(q)
}

Anno <- function(Geom,a,plot,b){# a<-ddata[[9]]$data[[2]];b<-ddata[[10]][2,-1]
  #Geom <- "text";a <- anno;b <- anno2;plot <- q;
  if (Geom == "segment"){
    for (z in 1:nrow(a)) {#z <- 1
      plot <- plot + annotate(Geom, x = a$x_from[z], xend = a$x_to[z], 
                              y = a$y_lab[z], yend = a$y_lab[z])
    }
  }else if (Geom == "text"){
    for (z in 1:nrow(a)) {#z <- 1
      plot <- plot + annotate(Geom, x = a$x_med[z], y = a$y_lab2[z],
                              label = b[1,a$name][1,z], size = 4)
    }
  }
  return(plot)
}

# test <- seq(from = maxy + 1,by = ifelse(maxy>1,20,1),length.out = 5)
################################################## 读取文件 #################################################################
mycol <-c("#CD3333","#483D8B","#458B00","#EEC900","#BFBC3B","#6E8B3D","#00688B","#C10077","#CAAA76","#EEEE00","#458B00","#8B4513","#008B8B","#6E8B3D","#8B7D6B","#7FFF00","#CDBA96","#ADFF2F")

mp1 <- read_tsv(opts$map) %>% rename("SampleID" = colnames(.)[1])
#如果有颜色定义的话，就可以直接定义了
if (opts$color != "none"){# opts$color = c("color.txt")
  sc <- read_tsv(opts$color,col_names = FALSE)
  sc <- sc[which(unlist(sc[,1]) %in% unique(unlist(mp1$group))),]
  mp1$group <- factor(mp1$group,levels = unlist(sc[,1]))
} else{
  sc <- cbind(unique(mp1$group),mycol[1:length(unique(mp1$group))]) %>% 
    as.data.frame(stringsAsFactors = FALSE) %>% as_tibble
  mp1$group <- factor(mp1$group,levels = unique(mp1$group))
}
sc <- sc %>% rename_all(vars(c("V1","V2"))) %>%
  map_at(vars(c("V1","V2")),function(x) fct_inorder(x)) %>% as_tibble

mpn <- nlevels(mp1$group)
gp1 <- as.vector(unique(mp1$group))
gpn <- unlist(lapply(1:length(gp1),function(g){# g = 1
  unlist(lapply(1:3,function(p){paste0(gp1[g],"-",c("median","mean","se")[p])}))
}))

# 宽数据
data1 <- read_tsv(opts$alpha)  %>%
  rename("SampleID" = colnames(.)[1]) %>%
  right_join(.,mp1) %>%
  select(-c(group,reads,label,coverage))

# 长数据
data1_2 <- data1 %>% left_join(mp1) %>% 
  gather(Alpha,value,-group,-SampleID) %>%
  mutate(Alpha = fct_inorder(Alpha))

# 计算基础统计量 中位数-均值-标准误 data1_3
###############################################################
data1_3 <- data1_2 %>% group_by(Alpha,group) %>% 
  summarize_at(vars(value),funs(mean, se = sd(.)/sqrt(n()),!!!p_funs)) %>% 
  ungroup() %>%
  mutate(IQR = sapply(1:nrow(.),function(x){
    str_c(Minus(.$`IQR50%`[[x]],2),"(",
          Minus(.$`IQR25%`[[x]],2),",",
          Minus(.$`IQR75%`[[x]],2),")")
  })) %>% select_at(vars(-contains("%"))) %>%
  mutate(Alpha = fct_inorder(Alpha),
         group = factor(group,levels = levels(mp1$group))) %>% 
  arrange(Alpha,group) %>% 
  group_by(Alpha) %>% 
  nest %>% 
  mutate(data = map(data,function(x) x %>% select(-group) %>% unlist(.) %>% 
                      .[order(as.numeric(gsub("([A-Za-z]+)([0-9]+)", "\\2", names(.))),
                              gsub("([A-Za-z]+)([0-9]+)", "\\1", names(.)))] %>% 
                      as.data.frame %>% t %>% as_tibble)) %>%# 保留四位小数
  unnest %>%
  rename_all(vars(c("Alpha",gpn)))

#######################################################################################
# 计算p-value、q-value
data2 <- data1_2 %>% 
  group_by(Alpha) %>% 
  nest %>% 
  mutate(Alpha = fct_inorder(Alpha),
         data = map(data,function(x) Pvalue(x,ts=opts$paired,ss=opts$test))) %>%
  mutate(data = map(data,function(x) # x <- data2$data[[1]]
    x %>% mutate(`q-value`= p.adjust(.$value,method = "fdr")) %>%
      select(-variable) %>% rename("p-value" = "value")))

ddata <- list()
ddata[[1]] <- data2 %>% 
  mutate(data = map(data,function(x) Ptable(x,"p-value"))) %>% 
  unnest %>% mutate_if(is.numeric,~Minus(.x,4))
ddata[[2]] <- data2 %>% 
  mutate(data = map(data,Sig_table)) %>% 
  unnest
ddata[[3]] <- ddata[[1]] %>% filter_at(vars(contains("-")), any_vars(. <= 0.05)) %>% 
  arrange_at(vars(contains("-"))) %>% 
  arrange(Alpha) %>%
  mutate(Alpha = factor(Alpha,levels=as.vector(Alpha))) %>%
  mutate_if(is.numeric,~Minus(.x,4))
ddata[[4]] <- ddata[[2]] %>% filter(Alpha %in% ddata[[3]]$Alpha) %>%
  mutate(Alpha = factor(Alpha,levels = levels(ddata[[3]]$Alpha))) %>%
  arrange(Alpha)
ddata[[5]] <- data2 %>% 
  mutate(data = map(data,function(x) Ptable(x,"q-value"))) %>% 
  unnest 
ddata[[6]] <- sapply(1:nrow(ddata[[1]]),function(x){# x <- 2
  sapply(2:ncol(ddata[[1]]),function(y){
    str_c("P=",Minus(ddata[[1]][x,y],4),";Q=",Minus(ddata[[5]][x,y],4))
  })
}) %>% t %>% as_tibble %>% rownames_to_column(.) %>% 
  rename_all(vars(colnames(ddata[[1]]))) %>% 
  mutate(Alpha = ddata[[1]]$Alpha)

# 是否显著
ddata[[7]] <- data2 %>% 
  mutate(Significant = map_chr(.$data,function(x){
    ifelse(any(x$`p-value`<=0.05),"TRUE","FALSE")})) %>%
  filter(Significant == "TRUE") %>% .$Alpha %>% as.vector

# Sig_mark position
ddata[[8]] <- data2 %>% 
  arrange(Alpha) %>% 
  mutate(data = lapply(1:nrow(.),function(x){#x <- 1
  Sig_mark(.$data[[x]],
           max(data1[,as.vector(.$Alpha)[[x]]]),map=mp1,value="p-value",v=1)
}))

# Sig_mark position
ddata[[9]] <- data2 %>% 
  arrange(Alpha) %>% 
  mutate(data = lapply(1:nrow(.),function(x){#x <- 1
    Sig_mark(.$data[x][[1]],
             max(data1[,as.vector(.$Alpha)[[x]]]),map=mp1,value="p-value",v=0.05)
  })) %>% filter(Alpha %in% ddata[[7]])
################################# 画图 ##################################
# boxplot_data 
PLOT <- function(Dd=data1_2,AN1=ddata[[8]]$data,
                 AN2=ddata[[6]],md="box",NM="",pt=FALSE){
  # Dd <- data1_2;AN2 <- ddata[[4]];AN1 <- ddata[[9]]
  box <- Dd %>% filter(Alpha %in% as.vector(AN2$Alpha)) %>% 
    mutate(Alpha = factor(Alpha,levels=levels(AN2$Alpha))) %>%
    group_by(Alpha) %>% nest %>% 
    mutate(Plot = lapply(1:nrow(.),function(y){# y <- 1
      GGBOX(data = .$data[[y]],ylab = as.vector(.$Alpha[y]),
            anno = AN1[[y]],anno2 = AN2[y,-1],tp=NM,pt=pt,mode=md)#
    }))
  # 把图合并哇
  Q <- box$Plot[[1]]
  for (i in 2:nrow(box)){
    Q <- Q + box$Plot[[i]]
  }
  Q <- Q + plot_layout()
  ggsave(str_c(NM,".",md,".plot.pdf"),Q,dpi = 600,
         width = nrow(box)*4,height = nrow(box)*3, device = cairo_pdf)
}

PLOT(AN1=ddata[[8]]$data,AN2=ddata[[6]],md=opts$mode,NM="P-Q",pt=opts$base)
PLOT(AN1=ddata[[8]]$data,AN2=ddata[[1]],md=opts$mode,NM="pvalue",pt=opts$base)
PLOT(AN1=ddata[[9]]$data,AN2=ddata[[4]],md=opts$mode,NM="p-sign",pt=opts$base)

cbind(data1_3,ddata[[1]] %>% select(-Alpha) %>% 
                      rename_all(vars(map_chr(colnames(ddata[[2]])[-1],~ str_c(.x,".Pvalue")))),
          ddata[[2]] %>% select(-Alpha) %>% 
            rename_all(vars(map_chr(colnames(ddata[[2]])[-1],~ str_c(.x,".Sig")))),
      ddata[[5]] %>% select(-Alpha) %>% 
        rename_all(vars(map_chr(colnames(ddata[[2]])[-1],~ str_c(.x,".Qvalue"))))) %>% 
  write_tsv("alpha_rarefac.test.xls")

