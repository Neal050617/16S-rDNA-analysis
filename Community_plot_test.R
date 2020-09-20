#Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i phylum.xls -m map-group.txt -p 0 -t rank -c color.txt
#Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i phylum.xls -m map-group.txt -p 0 -t rank -c color.txt
#Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i class.xls -m map-group.txt -p 0.01 -t rank -c color.txt
#Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i genus.xls -m map-group.txt -p 0.01 -t rank -c color.txt

#Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i phylum.xls -m map-group.txt -p 0 -t rank -w 0.6,0.3,0.15,0.5  -c color.txt
#Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i class.xls -m map-group.txt -p 0.01 -t rank -w 0.6,0.3,0.15,0.5 -c color.txt
#Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i genus.xls -m map-group.txt -p 0.01 -t rank -w 0.4,0.3,0.4,0.5 -c color.txt

#Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i phylum.xls -m map-group.txt -p 0.01 -a 1
#Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i phylum.xls -m map-group.txt -p 0.01 -a 1 -b 2

library(tidyverse)
library(magrittr)
library(optparse)
library(gridExtra)
library(grid)
library(reshape2)
options("endocing"="UTF-8")
options(scipen = 200)

if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="phylum.xls",
                help="丰度表格"),
    make_option(c("-s", "--select"), type="character", default="none",
                help="筛选表格:select.tsv"),
    make_option(c("-j", "--rename"), type="character", default="none",
                help="rename.list"),
    make_option(c("-m", "--map"), type="character", default="none",
                help="分组文件"),
    make_option(c("-n", "--nn"), type="integer", default=0,
                help="n=1,挑选每个样本丰度前p个的物种进行画图;
                n=2,挑选综合丰度排在前p个的物种进行画图；
                n=0,挑选丰度综合排在百分比为p的物种进行画图,p以下的为others"),
    make_option(c("-f", "--unif"), type="logical", default=T,
                help="要不要归一化"),
    make_option(c("-p", "--per"), type="double", default=0,
                help="筛选阈值"),
    make_option(c("-o", "--other"), type="logical", default=T,
                help="per筛选后是否保留others"),
    make_option(c("-a", "--average"), type="integer", default=0,
                help="0:不处理；1：组内求均值，阴影连接；2：组内排序，分层画图"),
    make_option(c("-v", "--CI"), type="double", default=0.95,
                help="置信区间"),
    make_option(c("-l", "--log"), type="double", default=2,
                help="box图取log值:0表示no;log2、log10"),
    make_option(c("-c", "--color"), type="character", default="none",
                help="指定颜色：color.txt"),
    make_option(c("-e", "--shape"), type="character", default="none",
                help="形状设置"),
    make_option(c("-b", "--bb"), type="integer", default=1,
                help="条形图1；气泡图2"),
    make_option(c("-t", "--test"), type="character", default="none",
                help="差异物种计算方法none、rank、t、aov"),
    make_option(c("-x", "--complex"), type="double", default=10,
                help="子母图阈值"),
    make_option(c("-w", "--cowplot"), type="character", default="none",
                help="子图尺寸修改,逗号分隔:x,y,width,heigth;0.5,0.5,0.4,0.4"),
    make_option(c("-z", "--size"), type="double", default="0",
                help="指定图形的宽度"),
    make_option(c("-y", "--sort"), type="logical", default=F,
                help="是否要排序"),
    make_option(c("-r", "--paired"), type="logical", default=F,
                help="是否成对检验"),
    make_option(c("-u", "--row"), type="double", default=30,
                help="legend一列的个数")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
}

####################################### 函数 ###########################################
Base_Statistics <- function(TT,ts=FALSE,ss=opts$test){# 
  tt <- length(unique(as.vector(TT$group)))
  TT <- TT %>% .[,-1]
  colnames(TT)[2] <- "data"
  set.seed(20190515)
  if (ss == "rank" & tt == 2)
    TEST <- wilcox.test(data ~ group,data = TT,paired=ts,exact = TRUE)
  if (ss == "rank" & tt > 2)
    TEST <- kruskal.test(data ~ group,data = TT)
  if (ss == "t" & tt == 2)
    TEST <- t.test(data ~ group,as.data.frame(TT),paired=ts,
                   var.equal = FALSE,alternative = "two.sided")
  if (ss == "aov")
    TEST <- aov(data ~ group,TT) %>% broom::glance(.)
  return(TEST$p.value)
}
#broom::glance(TEST)
#broom::tidy(TEST)
#broom::augment(model,data)

# 计算均值、标准差、置信区间
MSC <- function(T){# T <- tax_tbs$data[[1]]
  TT <- T %>% reshape2::melt(.,variable.name = "SampleID") %>% 
    as_tibble %>% left_join(Map) %>% select(-SampleID)
  TTT <- split(as.vector(unlist(TT$value)),TT$group)
  lapply(seq_along(TTT[[1]]),function(x){
    lapply(seq_along(TTT[[2]]),function(y){
      TTT[[1]][x] - TTT[[2]][y]
    })
  }) %>% unlist %>% as.numeric %>% as_tibble %>% 
    summarise_each(funs(count = n(), mean, 
                        se = sd(.)/sqrt(n()), 
                        ci = se * (qt(opts$CI/2 + .5, count-1))))  
}

bartheme <- function(){
  return(theme_bw() + ##设置主题
           theme(plot.title=element_text(size=rel(1),hjust=0.5),
                 plot.margin = unit(c(1, 3, 1, 1), "lines"),
                 axis.title=element_text(size=rel(1)),
                 axis.text.x=element_text(size=rel(1),angle = 90, vjust = 0.5, 
                                          hjust = 0.5,color = "black"),
                 axis.ticks.x=element_blank(),
                 axis.line.y=element_line(color="black"),
                 panel.grid.major=element_line(color="white"),
                 panel.grid.minor=element_line(color="white"),
                 panel.border=element_rect(color="white"),
                 legend.title=element_blank(),
                 legend.text = element_text(size = 6),
                 legend.key.size = unit(.4,'cm'),
                 legend.spacing.x = unit(0.1, 'cm'),
                 #legend.position = "bottom"
                 legend.justification = c(1,1)
           ))
}

bbtheme <- function(){
  return(theme_bw() + ##设置主题
           theme(plot.title=element_text(size=rel(1),hjust=0.5),
                 plot.margin = unit(c(1, 1, 1, 1), "lines"),
                 panel.grid.major=element_line(color="white"),
                 panel.grid.minor=element_line(color="white"),
                 axis.title=element_text(size=rel(1)),
                 axis.text.x=element_text(angle=30,hjust =1),
                 legend.title=element_blank(),
                 legend.text = element_text(size = 6),
                 legend.key.size = unit(.4,'cm'),
                 legend.spacing.x = unit(0.1, 'cm')
           ))
}

pp <- c(0.25, 0.5, 0.75)
p_names <- map_chr(pp, ~paste0("IQR",.*100, "%"))
p_funs <- map(pp, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
  set_names(nm = p_names)

Minus <- function(x,n){
  d1 <- 10^(-n)
  ifelse(x >= d1,as.vector(round(x,n)),paste0("<",d1))
}

################################## 读取数据 ############################################
tax_tb <- read_tsv(opts$input) %>% 
  rename('ID' = colnames(.)[1]) %>% 
  mutate(SUM = sapply(1:nrow(.),function(x) sum(.[x,2:ncol(.)]))) %>% 
  filter(SUM > 0) %>% arrange(desc(SUM)) %>% select(-SUM)

################################# 读取rename.pick文件,挑选样本,也可以修改样本名#########
if (opts$rename != "none") {
  pick <- read_tsv(opts$rename,col_names = FALSE) %>% 
    mutate(X1 = fct_inorder(X1)) %>%
    rename_all(~c("SampleID","Change_nm"))
  tax_tb <- tax_tb %>% gather(SampleID,value,-ID) %>% 
    inner_join(.,pick) %>%
    mutate_if(is.character,~ fct_inorder(.x)) %>%
    select(-SampleID) %>% 
    spread(Change_nm,value) %>%
    mutate(ID = as.character(ID))
}

# Uniform
if (opts$unif){
  tax_tb[,2:ncol(tax_tb)] <- sapply(2:ncol(tax_tb),function(x) tax_tb[,x]/sum(tax_tb[,x]))
}

################################# 指定物种 #########
if (opts$select != "none") {
  picks <- read_tsv(opts$select,col_names = FALSE) %>% mutate(X1 = fct_inorder(X1))
  tax_tb <- tax_tb %>% .[unlist(.$ID) %in% picks$X1,] %>% mutate(ID = factor(ID,levels = levels(picks$X1)))
}

################################# 读取group文件 ########################################
if (opts$map != "none"){
  Map <- read_tsv(opts$map) %>% 
    rename("SampleID" = colnames(.)[1],"group" = colnames(.)[2]) %>%
    mutate(group = fct_inorder(group))
  
  if (opts$rename != "none") {
    Map <- Map %>% inner_join(pick,.) %>% select(-SampleID) %>% rename("SampleID" = colnames(.)[1])
  }
  
  # 分组设计文件
  gp <- as.vector(unique(unlist(Map[,"group"])))
  
  if (opts$test != "none"){
    mapcol <- c("#CD3333","#483D8B","#458B00","#EEC900","#BFBC3B","#6E8B3D","#00688B","#C10077","#CAAA76","#EEEE00","#458B00","#8B4513","#008B8B","#6E8B3D","#8B7D6B","#7FFF00","#CDBA96","#ADFF2F")
    mypch <-c(21,22,24,23,25,3,4,16,15,17,18,7,8,9,10,11,12,13,14)
    
    # color.txt
    if (opts$color != "none"){# opts$color = c("color.txt")
      sc <- read.table(opts$color,sep="\t",comment.char = "",check.names = FALSE)
      sc <- sc[which(as.vector(sc[,1]) %in% unique(unlist(Map$group))),]
      Map$group <- factor(Map$group,levels = as.vector(sc[,1]))
    } else{
      sc <- cbind(levels(Map$group),mapcol[1:nlevels(Map$group)]) %>% 
        as.data.frame(stringsAsFactors = FALSE)
    }
    sc$V2 <- factor(sc$V2,levels = as.vector(sc[,2]))
    
    #设置形状
    if (opts$shape != "none"){
      sp <- read.table(opts$shape,sep="\t",comment.char="",check.names=F)
      sp <- sp[which(as.vector(sp[,1]) %in% unique(unlist(Map$group))),]
      sp[,2] <- factor(sp[,2],levels=unique(sp[,2]))
      mypch <- unique(sp[,3])
    }
  }
  
  tax_tb <- tax_tb %>% .[,c('ID',as.vector(Map$SampleID))] %>% 
    mutate(SUM = sapply(1:nrow(.),function(x) sum(.[x,2:ncol(.)]))) %>% 
    filter(SUM > 0) %>% arrange(desc(SUM)) %>% select(-SUM)
}

##文件名设置
ggbname <- unlist(strsplit(opts$input,split=".",fixed=TRUE))[1]
if(opts$map != "none"){
  if(opts$average == 1 & opts$test == "none")
    ggbname <- str_c(ggbname,".average.")
  if(opts$average == 2 & opts$test == "none")
    ggbname <- str_c(ggbname,".facet.")
  if(opts$average == 0 & opts$test != "none")
    if (opts$test == "rank"){
      gbbname <- str_c(ggbname,".",ifelse(length(gp)==2,
                                          ifelse(opts$paired == TRUE,
                                                 "Wilcoxon_rank_sum_paired",
                                                 "Wilcoxon_rank_sum_unpaired"),
                                          "Kruskal-Wallis_rank_sum"))
    }else{
      gbbname <- str_c(ggbname,".",opts$test,".")
    }
}
################################ 物种丰度挑选，有四种情况n'm'p‘s #######################
if (opts$nn == 1){
  max_names <- sapply(2:ncol(tax_tb),function(x) #x <- 2
    tax_tb$ID[order(as.numeric(unlist(tax_tb[,x])),
                    decreasing=T)[1:ifelse(nrow(tax_tb)>=opts$per,
                                           opts$per,
                                           nrow(tax_tb))]])
  if (opts$other){
    Ot <- tax_tb[!(tax_tb$ID %in% unique(as.vector(max_names))),]
    other <- sapply(2:ncol(Ot),function(y)sum(Ot[,y]))
    tax_tb <- tax_tb[which(tax_tb$ID %in% unique(as.vector(max_names))),] %>%
      add_row(ID = "Others")
    tax_tb[nrow(tax_tb),2:ncol(tax_tb)] <- other
  }else{
    tax_tb <- tax_tb[which(tax_tb$ID %in% unique(as.vector(max_names))),]
  }
}else if (opts$nn == 2){
  if (opts$other){
    Ot <- tax_tb[ifelse(nrow(tax_tb)>=opts$per,opts$per+1,nrow(tax_tb)):nrow(tax_tb),]
    other <- sapply(2:ncol(Ot),function(y)sum(Ot[,y]))
    tax_tb <- tax_tb[1:ifelse(nrow(tax_tb)>=opts$per,opts$per,nrow(tax_tb)),] %>%
      add_row(ID = "Others")
    tax_tb[nrow(tax_tb),2:ncol(tax_tb)] <- other
  }else{
    tax_tb <- tax_tb[1:ifelse(nrow(tax_tb)>=opts$per,opts$per,nrow(tax_tb)),]
  }
}else if (opts$nn == 0){
  minp <- sapply(1:nrow(tax_tb),function(y) #y <- 1
    all(tax_tb[y,-1]<=opts$per)|tax_tb$ID[y]=="Others")
  if (any(minp)){
    other <- sapply(2:ncol(tax_tb),function(y)sum(tax_tb[minp,y]))
    if(opts$other){
      tax_tb <- tax_tb[!minp,] %>% add_row(ID = "Others")
      tax_tb[nrow(tax_tb),2:ncol(tax_tb)] <- other
    }else{
      tax_tb <- tax_tb[!minp,]
    }
  }
}

write_tsv(tax_tb,str_c("percent.",opts$input))

################################ 是否计算均值 ###################################
if (opts$map != "none" & opts$average != 0){
  mg1 <- lapply(as.vector(unique(Map$group)),function(x) 
    tax_tb[,c(TRUE,as.vector(Map$group) %in% x)])
  
  if (opts$average == 1){# 阴影柱形图
    avg <- matrix(NA,nrow=nrow(tax_tb),ncol = nlevels(Map$group))
    for (i in 1:nlevels(Map$group)){#i <- 1
      for (j in 1:nrow(tax_tb)){#j <- 1
        avg[j,i] <- mean(as.numeric(unlist(mg1[[i]][j,-1])))
      }
    }
    rownames(avg) <- unlist(tax_tb$ID);colnames(avg) <- unique(Map$group)
    tax_tb <- avg %>% as.data.frame %>% 
      rownames_to_column %>% rename('ID' = 'rowname')
    write_tsv(tax_tb,str_c("average.percent.",opts$input))
  }
  if (opts$average == 2){#i <- 1 # 分层柱形图
    tax_tb <- lapply(1:nlevels(Map$group),function(i) 
      mg1[[i]][,c(1,order(as.numeric(unlist(mg1[[i]][1,-1])),
                          decreasing = TRUE)+1)]) %>% bind_cols %>%
      select(-matches("ID[0-9]+"))
    write_tsv(tax_tb,str_c("sort.percent.",opts$input))
  }
}

################################ 差异物种 ######################################
if (opts$map != "none" & opts$average == 0 & opts$test != "none"){
  
  # 带有分组的丰度表格
  tax_tb_gp <- tax_tb %>% .[,-1] %>% t %>% 
    as.data.frame %>% rownames_to_column %>% as_tibble %>%
    rename(SampleID = rowname) %>%
    left_join(Map) %>% 
    mutate(group = factor(group,levels = levels(Map$group)))
  colnames(tax_tb_gp) <- c("SampleID",as.vector(tax_tb$ID),'group')
  
  # 列名
  gpn <- unlist(lapply(1:length(gp),function(g){# g = 1
    unlist(lapply(1:3,function(p){paste0(gp[g],"-",c("median","mean","se")[p])}))
  }))
  
  # 计算均值和标准误# x <- colnames(tax_tb_gp)[2]
  Stats <- (lapply(colnames(tax_tb_gp)[2:(ncol(tax_tb_gp)-1)],function(x){#x <- colnames(tax_tb_gp)[2]
    stats <- tax_tb_gp[,c(x,"group")] %>% group_by(group) %>% 
      summarise_each(funs(mean, se = sd(.)/sqrt(n()),!!!p_funs),x) %>% 
      mutate(IQR = sapply(seq_along(.$group),function(x){
        str_c(Minus(.$`IQR50%`[x],4),"(",
              Minus(.$`IQR25%`[x],4),",",
              Minus(.$`IQR75%`[x],4),")")
      })) %>% select_at(vars(-contains("%"))) %>%
      dplyr::select(-group) %>% unlist(.) %>% 
      .[order(as.numeric(gsub("([A-Za-z]+)([0-9]+)", "\\2", names(.))),
              gsub("([A-Za-z]+)([0-9]+)", "\\1", names(.)))]
    names(stats) <- gpn
    return(stats)
  }))
  names(Stats) <- colnames(tax_tb_gp)[2:(ncol(tax_tb_gp)-1)]
  
  # 计算p-value、q-value ;合并表格
  Stats_p <- as.data.frame(t(data.frame(Stats))) %>% 
    mutate(`p-value` = (reshape2::melt(tax_tb_gp,
                                       id=c("SampleID","group"),
                                       variable.name='ID') %>% 
                          group_by(ID) %>% nest %>% 
                          mutate(data = map(data,function(x) 
                            Base_Statistics(TT=x,ts=opts$paired,ss=opts$test))) %>%
                          unnest(data) %>% rename(`p-value` = data) %>% .$`p-value`)) %>%
    mutate(`z-score` = qnorm(`p-value`/2)) %>%
    mutate(Sig_mark = ifelse(`p-value`>0.05,"",
                             ifelse(`p-value`>0.01,"*",
                                    ifelse(`p-value`>0.001,"**","***")))) %>% 
    mutate(`q-value`= p.adjust(as.numeric(unlist(.$`p-value`)),
                               n = nrow(.),method = "fdr")) %>% 
    rownames_to_column %>% as_tibble %>% rename('ID' = 'rowname') %>% 
    mutate(ID = colnames(tax_tb_gp)[2:(ncol(tax_tb_gp)-1)])
  
  ## 输出数据
  write_tsv(Stats_p,str_c(gbbname,".xls"))
  
  ##################################### 条形图 #####################################
  # 挑出显著性标识
  if (opts$select == "none"){
    SIGN <- Stats_p %>%
      filter(`p-value` <= 0.05) %>% 
      select("ID","Sig_mark")
  }else{
    SIGN <- Stats_p %>%
      select("ID","Sig_mark")
  }
  
  # 挑出柱形图绘图数据
  plot_tb <- tax_tb_gp %>% select(-SampleID) %>% group_by(group) %>% 
    nest %>% mutate(data = lapply(seq_along(data),function(y) {#y <- 1
      lapply(colnames(data[[y]]),function(z){#z <- rownames(Stats_p)[1]
        data[[y]][z] %>% summarise_each(funs(mean, se= sd(.)/sqrt(n())))
      })%>% bind_rows %>% mutate(ID = Stats_p$ID)
    })) %>% 
    unnest(data) %>% 
    right_join((Stats_p %>% filter(`p-value` <= 0.05) %>% .[,"ID"])) %>%
    mutate(ID = fct_inorder(ID)) %>%
    mutate(mean = mean*100,se = se*100)
  
  # 排序排序排序排序
  if (isTRUE(opts$sort)){
    # 判断递增还是递减
    z1 <- str_c(1:nlevels(plot_tb$group),collapse = "-")
    z2 <- str_c(nlevels(plot_tb$group):1,collapse = "-")
    zz <- plot_tb %>% group_by(ID) %>%
      nest %>%
      #mutate(Rank = sapply(seq_along(ID),function(x){# x <- 2
      #  oz <- zz$data[[x]]
      #  str_c(rank(oz$mean),collapse = "-")
      #})) %>% 
      mutate(Rank = map_chr(data,function(x){str_c(rank(x$mean),collapse = "-")})) %>%
      filter(Rank == z1|Rank == z2) %>%
      #mutate(Max = sapply(seq_along(ID),function(y){
      #  max(.$data[[y]]$mean)
      #})) %>% 
      mutate(Max = map_dbl(data,function(y){max(y$mean)})) %>%
      arrange(Max) %>%
      # mutate(Mean = sapply(seq_along(ID),function(y){
      #   mean(.$data[[y]]$mean)
      # })) %>% arrange(Mean) %>%
      mutate(Rank = factor(Rank,levels = c(z1,z2))) %>%
      arrange(Rank) %>%
      mutate(nRank = as.numeric(Rank) %>% 
               map_dbl(~ ifelse(.x == 1,-1,
                                ifelse(.x == 2,1,0)))) %>% 
      # mutate(Mean = Mean*nRank) %>%
      mutate(Max = Max*nRank) %>%
      group_by(Rank) %>% nest %>% 
      # mutate(data = data %>% map(function(x) arrange(x,desc(Mean)))) %>%
      mutate(data = data %>% map(function(x) arrange(x,desc(Max)))) %>%
      unnest(data) %>%
      .$ID %>% as.character
    
    plot_tb <- plot_tb %>% 
      filter(ID %in% zz) %>% 
      mutate(ID = factor(ID,levels=zz)) %>% 
      arrange(ID)
    
    SIGN <- SIGN %>% right_join(.,zz %>% enframe %>% select(-name) %>% rename(ID = value))
  }
  
  # 计算图形边界
  ycor1 <- max(as.numeric(as.vector(rowSums(plot_tb[,2:3]))))
  Width <- ifelse(opts$size==0,(log2(nrow(SIGN))*(1+length(gp)/10)+3),opts$size)
  Height <- 6
  
  bbar1 <- ggplot(data=plot_tb,aes(x=ID,y=mean,fill=group))+
    geom_errorbar(aes(ymax = mean + se, ymin = mean -  se),
                  position = position_dodge(0.9), width = 0.15, size = 0.2)+
    geom_bar(stat = "identity", position = "dodge")+
    scale_fill_manual(values=as.vector(unlist(sc$V2)))+
    xlab('')+ylab(str_c("Abundance(mean","\U00B1","SE)%"))+
    scale_y_continuous(expand = c(0, 0))+
    ggtitle(ggbname)+
    bbtheme()+
    theme(plot.margin = unit(c(1, 1, 1, 4), "lines"))
  
  ggsave((bbar1 + coord_cartesian(ylim=c(0,ycor1*0.03+ycor1))+
            annotate("text", x=SIGN$ID, y=ycor1*0.005+ycor1, label=SIGN$Sig_mark))
         ,filename=paste0(gbbname,"barplot.pdf"),
         width = Width,height = Height,limitsize = FALSE)
  
  ################################## 箱式图 ###########################################
  ## plot_box <- tax_tb_gp %>% select(-SampleID) %>% 
  ##   .[,c(SIGN$ID,"group")] %>%
  ##   melt(.) %>% as_tibble %>% 
  ##   mutate(variable = fct_inorder(variable),
  ##          value = sapply(1:nrow(.),function(x) 
  ##            ifelse(opts$log != 0,log(.$value[x],base = opts$log),.$value[x]*100)))
  ## plot_bx <- plot_tb %>% ungroup() %>%
  ##   mutate(mean = sapply(1:nrow(.),function(x) 
  ##     ifelse(opts$log != 0,log(.$mean[x]/100,base = opts$log),.$mean[x]))
  ##   ) %>%
  ##   select(ID,mean,group)
  ## ###### log转化2
  ## plot_box2 <- tax_tb_gp %>% select(c("SampleID",SIGN$ID,"group")) %>% 
  ##   mutate_if(is.numeric,function(x) map_dbl(x+10^(-8),function(y)
  ##     log(y,base = opts$log)))
  ## plot_box22 <- plot_box2 %>% gather(taxa,value,-group,-SampleID) %>%
  ##   mutate(taxa = fct_inorder(taxa))
  ## plot_bx2 <- plot_box2 %>% select(-SampleID) %>% group_by(group) %>% nest %>% 
  ##   mutate(data = map(data,function(x) x %>% summarise_each(mean))) %>%
  ##   unnest(data) %>% gather(taxa,value,-group) %>%
  ##   inner_join(.,Stats_p %>% filter(`p-value` <= 0.05) %>% .[,"ID"],by=c("taxa"="ID"))
  
  ###### 不做log转化
  plot_box <- tax_tb_gp %>% select(c("SampleID",SIGN$ID,"group")) %>%
    gather(taxa,value,-group,-SampleID) %>%
    mutate(taxa = fct_inorder(taxa),
           value = value*100)
  plot_bx <- plot_tb %>% ungroup() %>%
    select(ID,mean,group) %>% filter(mean > 0)
  ycor4 <- max(as.numeric(as.vector(plot_box$value)))+0.1
  
  bx <- ggplot()+
    geom_boxplot(plot_box,mapping = aes(x=taxa,y=value,color=group),
                 outlier.colour = NULL,outlier.size = 0.6,show.legend = F) +
    geom_point(plot_bx, position = position_dodge(width = 0.75),color="grey40",
               mapping = aes(x=ID,y=mean,group=group,fill=group,shape=group),size=1.5) +
    scale_fill_manual(values=as.vector(unlist(sc$V2)))+
    scale_color_manual(values=as.vector(unlist(sc$V2)))+
    scale_shape_manual(values=mypch)+
    #coord_cartesian(ylim=c(0,ycor4))+
    scale_x_discrete(name="")+
    annotate("text", x=SIGN$ID, y=ycor4-0.004, label=SIGN$Sig_mark)+
    ggtitle(ggbname)+
    bbtheme()+
    theme(plot.margin = unit(c(1, 1, 1, 4), "lines"))+
    ylab("Relative abundance(%)")

  if(opts$log!=0){
    library(scales)
    if(opts$log==2){
      bx <- bx + 
        scale_y_continuous(trans = log2_trans(),labels = scientific)
      #bx <- bx + scale_y_continuous(trans = 'log2')
    }else if(opts$log==10){
      bx <- bx + 
        scale_y_continuous(trans = log10_trans(),labels = scientific) +
        annotation_logticks(sides="l")
    }
  }

  ggsave(bx,filename=paste0(gbbname,"boxplot.pdf"),
         width = Width,height = Height,limitsize = FALSE)
  
  ################################# 子母图数据筛选 ####################################
  plot_tf <- plot_tb %>% group_by(ID) %>% nest %>% 
    mutate(data = lapply(seq_along(data),function(l) mean(data[[l]]$mean) )) %>%
    unnest(data) %>% mutate(complex = (max(.$data)/data >= opts$complex)) %>% 
    select(-data) %>% left_join(plot_tb,.) %>% 
    filter(complex) %>% ungroup()
  
  if (nrow(plot_tf) != 0){
    SIGN2 <- SIGN %>% .[which(.$ID %in% as.vector(unique(plot_tf$ID))),]
    ycor2 <- max(as.numeric(as.vector(rowSums(plot_tf[,2:3]))))
    
    bbar2 <- ggplot(data=plot_tf,aes(x=ID,y=mean,fill=group))+
      geom_errorbar(aes(ymax = mean + se, ymin = mean -  se),
                    position = position_dodge(0.9), width = 0.15, size = 0.2)+
      geom_bar(stat = "identity", position = "dodge", show.legend = F)+
      scale_fill_manual(values=as.vector(unlist(sc$V2)))+
      xlab('') + ylab('') + ggtitle('') +
      scale_y_continuous(expand = c(0, 0))+
      coord_cartesian(ylim=c(0,ycor2*0.06+ycor2))+
      annotate("text", x=SIGN2$ID, y=ycor2*0.01+ycor2, label=SIGN2$Sig_mark, size=2)+
      theme_bw() + ##设置主题
      theme(plot.margin = unit(c(0, 0, 0, 0), "lines"),
            axis.title=element_text(size=0.6),
            axis.text.x=element_text(angle=30,hjust =1,size = 6),
            axis.text.y=element_text(hjust =1,size = 6)
      )
    
    #pdf(paste0(gbbname,"barplot_complex.pdf"),width = Width,height = Height)
    #sub <- viewport(x = ifelse(1 - 0.45*nrow(SIGN2)/nrow(SIGN)>0.7,0.7,
    #                           1 - 0.45*nrow(SIGN2)/nrow(SIGN)), 
    #                y = 0.7, 
    #                width = nrow(SIGN2)/nrow(SIGN), height = 0.4)
    #bbar1
    #print(bbar2,vp = sub)
    #dev.off()
    if (opts$cowplot == "none"){
      x1 <- max((1-nrow(SIGN2)/nrow(SIGN))*0.8,nrow(SIGN2)/nrow(SIGN))
      y1 <- max((1-max(plot_tf$mean)/ycor1)*0.6,max(plot_tf$mean)/ycor1)
      x2 <- (1-x1)*0.5#weight
      y2 <- (1-y1)*0.7#height
      cc <- c(x1,y1,x2,y2)
    }else{#opts$cowplot <- '0.5,0.5,0.3,0.4'
      cc <- as.numeric(str_split(opts$cowplot,",")[[1]])
    }
    
    #cc <- as.numeric(c("0.2","0.5","0.1","0.4"))#x1,y1,x_,y_
    library(cowplot)
    sub <- ggdraw() +
      draw_plot((bbar1 + coord_cartesian(ylim=c(0,ycor1*0.1+ycor1)) + 
                   annotate("text", x=SIGN$ID, y=ycor1*0.05+ycor1, label=SIGN$Sig_mark)),
                0,0,1,1) +
      draw_plot(bbar2,cc[1],cc[2],cc[3],cc[4])
    ggsave(str_c(gbbname,"barplot_complex.pdf"),sub,width = Width,height = Height)
  }
  
  ###################################### STAMP 分析 ####################################
  yMin <- function(x){
    ifelse(x>0,x*0.9,x*1.1)
  }
  
  if (length(gp) == 2 & nrow(plot_tb) > 0){
    tax_tbs <- tax_tb %>% filter(ID %in% SIGN$ID) %>% 
      mutate(ID = factor(ID,levels = levels(plot_tb$ID))) %>% 
      arrange(ID) %>% group_by(ID) %>% nest %>% 
      mutate(data = map(data,MSC)) %>% unnest(data) %>% 
      mutate(mean = mean*100,se= se*100,ci=ci*100) %>% 
      rownames_to_column %>% 
      rename(x = rowname) %>%
      mutate(x = as.numeric(x)) %>%
      mutate(group = ifelse(mean > 0 ,gp[1],gp[2])) %>%
      mutate(group = factor(group,levels = levels(Map$group))) %>%
      ungroup()
    
    plot_tb <- plot_tb %>% 
      left_join((tax_tbs %>% select(ID) %>% rownames_to_column %>% 
                   rename(x = rowname))) %>% mutate(x = as.numeric(x))
    
    shading_data <- tax_tbs %>% select(x,ID) %>% 
      .[which(.$x %%2 == 0),] %>% mutate(min=0,max=max(x)+1) %>% 
      mutate(x_min=x-0.5,x_max=x+0.5,shade_color="grey44")
    
    bbar3 <- ggplot()+
      geom_rect(data = shading_data, 
                mapping = aes(xmin = x_min, xmax = x_max,ymin = -Inf, ymax = Inf),
                fill = shading_data$shade_color, alpha = 0.25) +
      geom_bar(data = plot_tb,
               aes(x=x,y=mean,fill=group),#show.legend = FALSE,
               stat = "identity", position = "dodge", color="black") +
      scale_x_reverse(labels=as.vector(tax_tbs$ID),breaks = tax_tbs$x, 
                      name = NULL,expand = c(0,0))+
      scale_y_continuous(expand = c(0,0),limits = c(0, max(plot_tb$mean)*1.1)) +
      scale_fill_manual(values=as.vector(unlist(sc$V2))) +
      coord_flip() +
      ylab('Proportion(%)')+
      theme_bw() +
      labs(title = "STAMP") +
      bbtheme() +
      theme(plot.title=element_text(size=rel(1), lineheight=.9, color = "white"),
            legend.position = "left",
            panel.grid = element_blank(),
            panel.border = element_blank())
    
    ycor3 <- tax_tbs %>% mutate(ymax = mean+ci,ymin = mean-ci) %>% select(ymax,ymin)
    
    bbar4 <- ggplot()+
      geom_hline(yintercept = 0, color="grey40", linetype="longdash", size=0.5) +
      geom_rect(data = shading_data, mapping = aes(xmin = x_min, xmax = x_max,
                                                   ymin = -Inf, ymax = Inf),
                fill = shading_data$shade_color, alpha = 0.25) +
      geom_errorbar(data = tax_tbs,aes(x = x,ymax = mean + ci, ymin = mean - ci),
                    position = position_dodge(0.9), width = 0.15, size = 0.2) +
      geom_point(data = tax_tbs,aes(x=x,y=mean,color=group),
                 show.legend = FALSE,size = 2) +
      scale_x_reverse(labels = SIGN$Sig_mark, breaks = tax_tbs$x, 
                      position = 'left',expand = expand_scale(add = c(0,0.4))) +
      scale_y_continuous(limits = c(min(yMin(ycor3$ymin)), max(ycor3$ymax)*1.1),
                         expand = c(0,0)) +
      scale_color_manual(values=as.vector(unlist(sc$V2))) + 
      ylab('Difference between proportion(%)') +
      xlab('p-value') +
      coord_flip() +
      theme_bw() +
      bbtheme() +
      labs(title = "95% confidence intervals") +
      theme(plot.title=element_text(size=rel(1), lineheight=.9),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            axis.ticks.y = element_blank())
    
    BBar <- grid.arrange(bbar3,bbar4,nrow=1,widths = c(1.2,1))
    ggsave(paste0(gbbname,"_barplot_stamp.pdf"),BBar,
           width = 8+(max(sapply(as.vector(plot_tb$ID),nchar)))^(1/3),
           height = (log2(nrow(SIGN))*(1+length(gp)/10)+2))
  }
}

############################### 物种分布条形图 ################################
if (opts$test == "none"){
  ##设置图片大小
  if (opts$size == 0){
    wdth=dim(tax_tb)[2]-1
    if (wdth<4){wdth=4}
    if (wdth>=4){wdth=5.2+wdth*.2}
    if (wdth>=200){wdth=45}
  }else{
    wdth=opts$size
  }
  
  hdth=dim(tax_tb)[1]
  if (hdth<10){hdth=5}
  if (hdth>=10&hdth<=30){hdth=6}
  if (hdth>30){hdth=floor(hdth/5)}
  
  mycol <- c(35, 107, 51, 144, 36, 37, 7, 12, 84, 116, 121, 77, 56, 
             386, 373, 423, 435,471, 438, 512, 130, 52, 47, 6, 11, 43, 
             54, 367, 382, 422, 4, 8, 375, 124, 448, 419, 614, 401, 403, 
             613, 583, 652, 628, 633, 496, 638, 655, 132, 503, 24)
  mycol <- colors()[rep(mycol,20)]
  mycol <- rev(mycol[1:nrow(tax_tb)])
  
  tax_tb_bar <- tax_tb %>% melt(.) %>% as_tibble %>% 
    mutate(ID = factor(ID,levels=rev(unique(ID))))
  
  if (opts$map != "none" & opts$average == 2){
    tax_tb_bar <- tax_tb_bar %>% left_join(.,Map,by=c("variable"="SampleID")) %>% 
      mutate(variable = fct_inorder(variable))
  }else if(opts$map != "none" & opts$average == 1){
    tax_tb_bar <- tax_tb_bar %>% mutate(variable = fct_inorder(variable))
  }else if(opts$map != "none" & opts$average == 0){
    tax_tb_bar <- tax_tb_bar %>% left_join(.,Map,by=c("variable"="SampleID")) %>% 
      mutate(variable = fct_inorder(variable))
  }else if(opts$map == "none" & opts$average == 0){
    tax_tb_bar <- tax_tb_bar %>% 
      mutate(variable = factor(variable, 
                               levels = c(str_sort(colnames(tax_tb)[-1], 
                                                   numeric = TRUE))))
  }
  
  if (opts$bb==1){
    if (opts$average==1){##面积图与柱形图合并
      tax_tb_bar$vv <- as.numeric(tax_tb_bar$variable)
      tax_tb_bar2 <- rbind(transform(tax_tb_bar, vv = vv-.3),
                           transform(tax_tb_bar, vv = vv+.3)) %>% as_tibble
      B <- ggplot(tax_tb_bar2,aes(x=vv,y=value,fill=ID)) + 
        geom_area(alpha=.6)+
        geom_bar(data=tax_tb_bar,aes(x=vv,y=value,fill=ID),color = "grey88",
                 position="stack", stat="identity", width=0.6, size=0.1)+
        guides(fill=guide_legend(nrow=opts$row,byrow=F,reverse=TRUE))+
        scale_fill_manual(values=mycol) + ##设置颜色
        scale_y_continuous(labels=seq(0,100,25),expand = c(0, 0)) +
        scale_x_continuous(name ="",breaks=1:max(tax_tb_bar$vv),
                           labels=levels(tax_tb_bar$variable))+
        ylab("Relative abundance(%)") +
        ggtitle(ggbname) +
        bartheme() ##设置主题
      ggsave(B,filename = str_c(ggbname,".barplot.pdf"), dpi=300, 
             height = 6,width= wdth,limitsize = FALSE)
    }else{
      B <- ggplot(tax_tb_bar,aes(x=variable,y=value,fill=ID)) + 
        geom_bar(stat="identity", width = 0.8) +
        scale_fill_manual(values=mycol) + ##设置颜色
        guides(fill = guide_legend(reverse=TRUE,nrow=opts$row,byrow=F)) + ##将legend翻转显示
        scale_x_discrete(name="") +
        scale_y_continuous(labels=seq(0,100,25),expand = c(0, 0)) +
        ylab("Relative abundance(%%)") +
        ggtitle(ggbname) +
        bartheme() ##背景主题等设置
      if(opts$map != "none"){
        B <- B + facet_grid(.~group,scale = "free",space = "free")+
          theme(strip.text=element_text(colour = "black", face="bold"), 
                strip.background=element_rect(fill="white",colour = "white"))
      }
      ggsave(B,filename = str_c(ggbname,".barplot.pdf"), dpi=300, 
             height = 6,width=wdth,limitsize = FALSE)
    }
  }
  
  if (opts$bb==2){#气泡图
    bubble <- ggplot(tax_tb_bar)+
      geom_hline(aes(x=variable,y=ID,yintercept = 1:nrow(tax_tb_bar)),
                 linetype="dashed",alpha=.2)+
      geom_vline(aes(x=variable,y=ID,xintercept = 1:nrow(tax_tb_bar)),
                 linetype="dashed",alpha=.2)+
      geom_point(aes(x=variable,y=ID,size=value,color=ID),shape=16)+
      scale_size_area(max_size=10)+
      guides(fill = guide_legend(reverse=TRUE)) + ##将legend翻转显示
      labs(color="Taxonomy",size="Relative")+
      ylab("Taxonomy")+
      xlab("Sample")+
      theme_bw() +
      theme(panel.grid.major.x=element_line(linetype="dashed"),
            #plot.margin=margin(5,5,5,5,unit="pt"),
            axis.text.x=element_text(size=rel(1),angle = 90, 
                                     vjust = 0.5, hjust = 0.5,color = "black"),
            axis.text=element_text(size=15,hjust=0.5),
            legend.background=element_rect(fill="white", colour="black"),
            legend.key=element_blank(),
            legend.key.size=unit(.6,'cm')
      )
    ggsave(bubble,filename = str_c(ggbname,".bubble.pdf"), dpi=300, 
           height=(hdth+8) ,width=(wdth+10),limitsize = FALSE)
  }
}
