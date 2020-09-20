#丰度筛选
#差异筛选
#核心微生物筛选
############################################# 载入包 ###################################
library(optparse)
library(corrplot) #载入corrplot这个包
library(Hmisc) #载入Hmisc这个包
library(gplots)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(reshape2)
library(PerformanceAnalytics)
library(tidyverse)
library(WGCNA)
library(igraph)
library(gridExtra)
library(grid)
options("endocing"="UTF-8")

if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="rarefac.otu_genus.xls",
                help="丰度表格;rarefac.otu_genus.xls"),
    make_option(c("-s", "--select"), type="character", default="none",
                help="筛选表格:select.txt挑选物种"),
    make_option(c("-m", "--map"), type="character", default="map-group.txt",
                help="分组文件:map-group.txt"),
    make_option(c("-e", "--env"), type="character", default="env.txt",
                help="生理数据:env.txt"),
    make_option(c("-t", "--test"), type="character", default="rarefac.Wilcoxon_rank_sum_unpaired.xls",
                help="p值文件,wilcox.otu.HC-DN.xls"),
    make_option(c("-p", "--per"), type="double", default=0.01,
                help="丰度筛选"),
    make_option(c("-o", "--order"), type="character", default="AOE",
                help="original,AOE,FPC,hclust,aplhabet"),
    make_option(c("-r", "--rho"), type="double", default=0,
                help="rho threshhold,0.3;0.5;"),
    make_option(c("-v", "--value"), type="double", default=0.2,
                help="核心微生物筛选"),
    make_option(c("-f", "--tranfm"), type="logical", default=FALSE,
                help="画图时转置"),
    make_option(c("-w", "--width"), type="double", default=30,
                help="宽"),
    make_option(c("-x", "--height"), type="double", default=30,
                help="高"),
    make_option(c("-a", "--mar"), type="character", default="2-20-20-2",
                help="bottom, left, top, right"),
    make_option(c("-u", "--feature"), type="character", default="feature_importance_scores-all.txt",
                help="feature_importance_scores.txt"),
    make_option(c("-y", "--Meande"), type="double", default=0.0001,
                help="feature_importance_scores 阈值"),
    make_option(c("-j", "--pvalue"), type="double", default=0.05,
                help="pvalue"),
    make_option(c("-k", "--qvalue"), type="double", default=0,
                help="pvalue")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
}

###################################### 函数 ######################################
Tidy_env <- function(x){# x <- env
  colSums(x*0+1,na.rm = T) >= nrow(x)*0.8 %>% as.vector
}

Fix_na <- function(x){#x <- env$data[[1]][,3]
  if (!any(is.na(x))){
    x
  }else{
    ss <- x %>% unlist %>% as.numeric %>%
      mean(.,na.rm = TRUE)
    x %>% replace(., is.na(.),ss) 
  }
}

#先定义一个标记*号的函数
noteTrans <- function(matrix,C1,C2,C3){
  P_list <- as.vector(matrix)
  i <- 1
  note <- c()
  while (i <= length(P_list)){
    note[i] = ""
    if (is.na(P_list[i]) == TRUE){P_list[i] = 0}
    if (P_list[i] >0 & P_list[i] <= C1){note[i] = "*"}
    if (P_list[i] >0 & P_list[i] <= C2){note[i] = "**"}
    if (P_list[i] >0 &P_list[i] <= C3 ){note[i] = "***"}
    i = i + 1	
  }
  noteM <- matrix(note,byrow=F,nrow=dim(matrix)[1])
  return(noteM)	
}

# point size
LM <- function(oo,y=c(1,3)){#oo <- data2 %>% colSums
  oo <- log10(oo)
  x = c(min(oo),max(oo))
  y = y
  o <- lm(y ~ x)
  predict(o,data.frame(x = oo))
}

# Genus_threshold
Genus_threshold <- function(nc,otu_coverage,genus){
  genus_threshold <- enframe(sapply(1:length(nc),function(y) 
    sum(as.numeric(otu_coverage >= nc[y])))) %>% 
    mutate(remian = length(otu_coverage) - value,
           coverage = nc) %>%
    mutate(per = sapply(1:length(nc),function(o){
      sum(unlist(genus[otu_coverage >= .$coverage[o],2:(ncol(genus))])) / 
        sum(unlist(genus[,2:(ncol(genus))]))
    })) %>% select(-name)
  
  p1 <- ggplot(genus_threshold,aes(coverage,per,color = "blue")) + 
    geom_vline(xintercept = opts$value, linetype="dotted") + 
    geom_point() + geom_line() + guides(color = FALSE)
  p2 <- ggplot(genus_threshold,aes(coverage,value,color = "red")) + 
    geom_vline(xintercept = opts$value, linetype="dotted") + 
    geom_point() + geom_line() + guides(color = FALSE)
  pp <- cowplot::plot_grid(p1, p2, labels = c('A', 'B'))
  ggsave("threshold.png",pp,width = 12,height = 8)
}

# 计算相关性系数 ： 不包含环境因子
COR0 <- function(genus,nm,WW=12,HH=12,TF=FALSE){#nm <- str_c(unique(Map$group),collapse = "-")
  Data <- genus %>% 
    mutate(SampleID = fct_inorder(SampleID)) %>% 
    arrange(SampleID) %>% .[,-1] %>% as.matrix %>% apply(.,2,as.numeric)
  set.seed(20190506)
  Cor <- Hmisc::rcorr(Data,type = "spearman")
  Cor$P <- Cor$P %>% replace_na(0)
  write_tsv(data.frame("rho"=rownames(Cor$r),Cor$r),str_c(nm,".rho.xls"))
  write_tsv(data.frame("p-value"=rownames(Cor$P),Cor$P),str_c(nm,".p-value.xls"))
  
  mynote <- noteTrans(Cor$P,0.05,0.01,0.001)
  if (TF == TRUE){
    mynote <- t(mynote)
    Cor$r <- t(Cor$r)
    Cor$P <- t(Cor$P)
  }
  Name = paste(nm,".corrplot.pdf",sep="")
  pdf(Name,w=WW,h=HH)
  corrplot(Cor$r, type="upper", order= "hclust", p.mat = Cor$P,method = 'circle',
           insig = "blank",# 不显著的cor值设置为“不显示”
           sig.level = 0.05,  # P值显著性设置
           tl.cex= 1,  #字体大小
           tl.srt= 45#字体角度
  )  
  dev.off()
  
}

# 计算相关性系数 
COR1 <- function(DATA,genus,hcolor,nm,ord,Rho,WW=15,HH=10,Mar="2-10-10-2",TF=FALSE){
  #DATA <- env;hcolor <- "RdYlGn";nm <- "test1";ord <- opts$order;Rho <- 0;
  #WW<-15;HH<-10;Mar<-"2-10-10-2";TF<-opts$tranfm
  data1 <- DATA[,-1] %>% mutate_if(is.character,function(x) 
    factor(x) %>% as.numeric) %>% as.matrix
  data2 <- genus %>% 
    mutate(SampleID = factor(SampleID,levels = levels(DATA$SampleID))) %>% 
    arrange(SampleID) %>% .[,-1] %>% as.matrix %>% apply(.,2,as.numeric)
  set.seed(20190506)
  Cor <- Hmisc::rcorr(data2,data1,type = "spearman")
  GCor <- Hmisc::rcorr(data2,type = "spearman")
  ####################################### rho out ###############################
  hnum <- dim(data2)[2]
  Cor$r <- Cor$r[1:hnum,(hnum+1):dim(Cor$r)[1]]
  Cor$P <- Cor$P[1:hnum,(hnum+1):dim(Cor$P)[1]]
  write_tsv(data.frame("rho"=rownames(Cor$r),Cor$r),str_c(nm,".cross.rho.xls"))
  write_tsv(data.frame("p-value"=rownames(Cor$P),Cor$P),str_c(nm,".cross.p-value.xls"))
  
  write_tsv(data.frame("rho"=rownames(GCor$r),GCor$r),str_c(nm,".rho.xls"))
  write_tsv(data.frame("p-value"=rownames(GCor$P),GCor$P),str_c(nm,".p-value.xls"))
  # data plot preparation
  n1 <- ncol(data1)
  n2 <- ncol(data2)
  grp <- colnames(data1)
  subx <- seq(-(n1-1),0,by=1)# 分组的X坐标
  suby <- sort(seq(from = 1,by=n2/n1,length.out=n1),decreasing = T) # 分组的Y坐标
  
  Df <- data.frame(
    grp = rep(grp, each = n2), # 分组名称，每个重复n次
    subx = rep(subx, each = n2), # 组X坐标，每个重复n次
    suby = rep(suby, each = n2), # 组Y坐标，每个重复n次
    x = rep(0:(n2 - 1) - 0.5, n1), # 变量连接点X坐标
    y = rep(n2:1, n1), # 变量连接点Y坐标
    Genus = colnames(data2)
  ) %>% as_tibble
  
  df_segment <- cbind(Df,
                      "rho"=melt(Cor$r,varnames = c("genus","grp"),
                                 value.name = "rho")$rho,
                      "p"=melt(Cor$P,varnames = c("genus","grp"),
                               value.name = "p")$p) %>% as_tibble
  #barplot(1:3,col =c("#FF8E8A","#00CCCC","#999999"))
  df_segment <- df_segment %>% 
    mutate(
      lcol = ifelse(p > 0.05, 0, 1), 
      lwd = ifelse(rho >0,"#FF8E8A",ifelse(rho<0,"#00CCCC","#999999")),
      rank = sapply(1:nrow(.),function(i){
        I <- seq(1,0,-0.1)
        max(I[abs(.$rho)[i] >= I])*10})
    ) %>%
    mutate(grp = fct_inorder(grp))
  
  point_size <- LM(data2 %>% colSums,y=c(1,2))
  text_posi <- df_segment %>% distinct(grp,subx,suby)
  ################################### corrplot作图 ################################
  cor_name = paste(nm,".corrplot.pdf",sep="")
  #insig = c('pch', 'p-value', 'blank', 'n', 'label_sig')
  #method = c('circle', 'square','ellipse','number','shade','color','pie')
  pdf(cor_name,w=WW,h=HH)
  set.seed(20190506)
  CC <- corrplot(GCor$r, type="upper", p.mat = GCor$p,method = 'circle',order= ord,
                 addrect = 1,rect.col = "blue",rect.lwd = 2,tl.pos = "n",
                 mar = as.numeric(str_split(Mar,"-")[[1]]),#c(2,10,10,2),
                 col = rev(colorRampPalette(brewer.pal(9,"RdBu"))(200)),
                 insig = 'blank',# 不显著的cor值设置为“不显示”
                 sig.level = 0.05,  # P值显著性设置
                 tl.cex= 1,  #字体大小
                 tl.srt= 45#字体角度
  )
  #dev.off()
  if (opts$order != "original"){
    df_segment <- bind_cols(df_segment[,c(4,5)],
                            full_join(df_segment[,c(1,2,3)] %>% 
                                        distinct(grp,.keep_all = T),
                                      df_segment[,c(1,6,7,8,9,10,11)] %>% 
                                        mutate(Genus = factor(Genus,levels=rownames(CC))) %>% 
                                        arrange(Genus)))
  }
  
  write_tsv(df_segment,str_c(nm,".corrplot.data.xls"))
  
  df_segment2 <- df_segment %>% filter(lcol > 0) %>% 
    filter(p <= 0.05 & p >0.01) %>% filter(abs(rho)>=Rho)
  df_segment3 <- df_segment %>% filter(lcol > 0) %>% 
    filter(p <= 0.01) %>% filter(abs(rho)>=Rho)
  
  segments(df_segment2$subx, df_segment2$suby, df_segment2$x, df_segment2$y, 
           lty = 'dotdash', lwd = df_segment2$rank, col = df_segment2$lwd, xpd = TRUE)
  segments(df_segment3$subx, df_segment3$suby, df_segment3$x, df_segment3$y, 
           lty = 'solid', lwd = df_segment3$rank, col = df_segment3$lwd, xpd = TRUE)
  points(df_segment$subx, df_segment$suby, pch = 24, col = 'black', 
         bg = 'black', cex = 1.5, xpd  = TRUE)
  points(df_segment$x,df_segment$y,pch=21,col="black",bg="black",
         cex=point_size,xpd=TRUE)
  text(text_posi$subx - 0.5, text_posi$suby - 0.5, labels = text_posi$grp, 
       adj = c(0.8, 0.5),cex = 1.2, xpd = TRUE, font = 2)
  
  xx <- min(text_posi$subx)
  yy <- head(tail(text_posi$suby, n=2), n=1)
  
  labels01 <- c('<= 0.01','0.01 < x <= 0.05')
  labels02 <- c('positive','negative')
  labels_x <- rep(xx-2, 2)
  labels_y <- c(yy, yy-0.5)
  
  text(labels_x[1]-0.5, yy+1, "Spearman's rho", adj = c(0, 0.5), cex = 1.2, font = 2, xpd = TRUE)
  text(labels_x, labels_y+0.5, labels02, adj = c(0, 0.5), cex = 1.2, xpd = TRUE)
  points(labels_x - 0.5, labels_y+0.5, pch = 20, 
         col = c("#FF8E8A","#00CCCC"),
         cex = 3, xpd = TRUE)
  
  #lines_x <- c(-16, -16)
  #lines_y <- c(2.5,1)
  #text(-18, 4, "p-value", adj = c(0, 0.5), cex = 1.2, font = 2, xpd = TRUE)
  lines_x <- rep(xx-2, 2)
  lines_y <- c(yy-1, yy-1.5)
  text(labels_x[1]-0.5, yy-0.5, "p-value", adj = c(0, 0.5), cex = 1.2, font = 2, xpd = TRUE)
  
  text(lines_x, lines_y, labels01, adj = c(0, 0.5), cex = 1.2, xpd = TRUE)
  segments(lines_x[1] - 0.5, lines_y[1], lines_x[1], lines_y[1], lwd = 2.5, 
           lty = 'solid', col = '#B3B3B3', xpd = TRUE)
  segments(lines_x[2] - 0.5, lines_y[2], lines_x[2], lines_y[2], lwd = 2.5, 
           lty = 'dotdash', col = '#B3B3B3', xpd = TRUE)
  
  text(seq_along(1:n2),n2+0.8,colnames(CC), adj = c(0, 0.5), cex = 0.8, 
       font = 2, xpd = TRUE, srt = 30)
  
  dev.off()
  
  ###################################### color #################################
  if (hcolor != ""){#"RdYlBu"
    mycol = colorRampPalette(rev(brewer.pal(n = 7, name = hcolor)))(25)
  }else{
    mycol <- colorpanel(25,"#87CEEB","white",high="#B03060")
  }
  ###################################### heatmap.2 画图 ########################
  mynote <- noteTrans(Cor$P,0.05,0.01,0.001)
  Name = paste(nm,".heatmap_env.pdf",sep="")
  pdf(Name,w=WW,h=HH)
  if (TF == TRUE){
    mynote <- t(mynote)
    Cor$r <- t(Cor$r)
    Cor$P <- t(Cor$P)
  }
  heatmap.2(Cor$r, trace="none",density="none",col=mycol,
            margins=as.numeric(str_split(Mar,"-")[[1]])[2:3], #图主体的位置，主要是避免名字太长的row或者column
            cexRow=1,# row字体的大小
            cexCol=1, # Column字体的大小
            srtCol=30,#Column字体的角度
            srtRow=45,
            cellnote=mynote, Colv = FALSE,
            notecex=3,  #注释信息的大小
            notecol="black", #注释信息的颜色
            dendrogram = 'row' #画树否 "both","row","column","none"
  )
  dev.off()
}

# 计算相关性系数 
COR2 <- function(DATA,genus,hcolor,nm,TF=FALSE){
  #DATA <- env;hcolor <- "RdYlGn";nm <- "genus_env"
  OUT <- inner_join(DATA,genus) %>% select(-SampleID) %>% as.matrix
  set.seed(20190506)
  Cor <- Hmisc::rcorr(OUT,type = "spearman")
  write_tsv(data.frame("rho"=rownames(Cor$r),Cor$r),str_c(nm,".rho.xls"))
  write_tsv(data.frame("p-value"=rownames(Cor$P),Cor$P),str_c(nm,".p-value.xls"))
  ###################################### color #################################
  if (hcolor != ""){#"RdYlBu"
    mycol = colorRampPalette(rev(brewer.pal(n = 7, name = hcolor)))(25)
  }else{
    mycol <- colorpanel(25,"#87CEEB","white",high="#B03060")
  }
  ####################################### 只展示物种和环境因子 #################
  hnum <- dim(DATA)[2]-1
  RR <- Cor$r[(hnum+1):dim(Cor$r)[1],1:hnum]
  PP <- Cor$P[(hnum+1):dim(Cor$P)[1],1:hnum]
  ### heatmap.2 画图 
  mynote <- noteTrans(PP,0.05,0.01,0.001)
  Name = paste(nm,".heatmap_env_part.pdf",sep="")
  pdf(Name,w=10,h=15)
  heatmap.2(RR, trace="none",density="none",col=mycol,
            margins=as.numeric(str_split(Mar,"-")[[1]])[2:3], #图主体的位置，主要是避免名字太长的row或者column
            cexRow=0.7,# row字体的大小
            cexCol=1, # Column字体的大小
            srtCol=30,#Column字体的角度
            srtRow=45,
            cellnote=mynote, 
            notecex=2,  #注释信息的大小
            notecol="black", #注释信息的颜色
            dendrogram = 'none' #画树否 "both","row","column","none"
  )
  dev.off()
  
  ###################################### heatmap.2 画图 ########################
  mynote <- noteTrans(Cor$P,0.05,0.01,0.001)
  if (TF == TRUE){
    mynote <- t(mynote)
    Cor$r <- t(Cor$r)
    Cor$p <- t(Cor$p)
  }
  Name = paste(nm,".heatmap_env.pdf",sep="")
  pdf(Name,w=20,h=20)
  heatmap.2(Cor$r, trace="none",density="none",col=mycol,
            margins=as.numeric(str_split(Mar,"-")[[1]])[2:3], #图主体的位置，主要是避免名字太长的row或者column
            cexRow=0.7,# row字体的大小
            cexCol=0.7, # Column字体的大小
            srtCol=30,#Column字体的角度
            srtRow=45,
            cellnote=mynote, Rowv = F, Colv = F,
            notecex=1,  #注释信息的大小
            notecol="black", #注释信息的颜色
            dendrogram = 'none' #画树否 "both","row","column","none"
  )
  dev.off()
}

#################################### read in ########################
# 环境因子整理
Map <- read_tsv(opts$map) %>% rename(SampleID = colnames(.)[1])

genus <- read_tsv(opts$input) %>% 
  rename(SampleID = colnames(.)[1]) %>% 
  .[,c("SampleID",Map$SampleID)]

if (file.exists(opts$env)) {
  env <- read_tsv(opts$env) %>% 
    rename(SampleID = colnames(.)[1]) %>% 
    inner_join(Map[,1]) %>%
    mutate_if(is.character,funs(fct_inorder(.))) %>%
    mutate_at(c(1),funs(as.character)) %>%
    mutate_if(is.factor,funs(as.numeric)) %>%
    #mutate_if(is.numeric,funs(Tidy_env)) %>%
    .[,c(TRUE,Tidy_env(.[,2:(ncol(.)-1)]),TRUE)] %>%
    inner_join(Map) %>% mutate(group = fct_inorder(group)) %>%
    #group_by(group) %>%
    #nest %>% 
    #mutate(data = map_dbl(data,function(y) 
    #  y %>% mutate_at(colnames(y)[2:ncol(y)],Fix_na))) %>%
    #unnest %>%
    select(-group) %>% 
    mutate(SampleID = fct_inorder(SampleID))
  
  Map <- Map %>% inner_join(env[,1]) %>% mutate_if(is.character,funs(fct_inorder(.)))
  genus <- genus[,c("SampleID",as.character(env$SampleID))]
}

# Uniform
genus[,2:ncol(genus)] <- sapply(2:ncol(genus),function(x) genus[,x]/sum(genus[,x]))

# 丰度筛选
if(opts$per != 0){
  genus <- genus %>% 
    mutate(SELECT = sapply(1:nrow(.),function(x){
      any(.[x,2:ncol(.)]>=opts$per)})) %>% 
    .[.$SELECT,] %>% as_tibble %>% select(-SELECT)
}

###################################### core microbiom #####################
otu_coverage <- apply(genus[,2:ncol(genus)],1,function(x) 
  length(x[x>0])/(ncol(genus)-1))
# threshold
nc <- sort(unique(otu_coverage))
Genus_threshold(nc,otu_coverage,genus)
# output data
genus <- genus[otu_coverage >= opts$value,]
#data <- data[otu_coverage >= 0.2,]

##################################### feature—selection #######################
if (opts$feature != "none"){
  Meande <- read_tsv(opts$feature) %>% 
    filter(Mean_decrease_in_accuracy >= opts$Meande) %>% .$Feature_id
  
  genus_md <- sapply(genus$SampleID,function(x) strsplit(x," ", fixed=TRUE)[[1]][1])
  genus <- genus[genus_md %in% Meande,]
}

#################################### analysis #################################
if(opts$test == "none" & opts$select == "none"){
  # 转置
  genus <- genus %>% gather(tax,value,-SampleID) %>% 
    mutate(SampleID = fct_inorder(SampleID),
           tax = fct_inorder(tax)) %>% 
    spread(SampleID, value) %>%
    rename(SampleID = colnames(.)[1])
  
  COR0(genus,str_c(unique(Map$group),collapse = "-"),opts$width,opts$height,opts$tranfm)
}

# select genus
if (opts$select != "none" & opts$test == "none"){
  selectg <- read_tsv(opts$select,col_names = FALSE) %>% 
    rename(SampleID = X1)
  
  genus <- genus %>% inner_join(selectg[,1]) %>% 
    gather(tax,value,-SampleID) %>% 
    mutate(SampleID = fct_inorder(SampleID),
           tax = fct_inorder(tax)) %>% 
    spread(SampleID, value) %>%
    rename(SampleID = colnames(.)[1])
  
  COR1(env,genus,"RdYlGn",
       str_c(str_c(unique(Map$group),collapse = "-"),".select."),
       opts$order,opts$rho,opts$width,opts$height,opts$mar,opts$tranfm)
  COR2(env,genus,"RdYlGn",str_c(str_c(unique(Map$group),collapse = "-"),".select."),opts$tranfm)
}

# wilcox genus
if (opts$test != "none"){
  test <- read_tsv(opts$test) %>% 
    filter(`p-value` <= opts$pvalue) %>% 
    rename(SampleID = colnames(.)[1]) #%>% filter(SampleID %in% colnames(genus)[-1])
  
  if (opts$qvalue != 0){
    test <- test %>% filter(`q-value` <= opts$qvalue)
  }
  
  if (opts$select == "none"){
    #nm <- unlist(inner_join(genus[,1],test[,1]))
    genus <- genus %>% inner_join(test[,1]) %>% 
      gather(tax,value,-SampleID) %>% 
      mutate(SampleID = fct_inorder(SampleID),
             tax = fct_inorder(tax)) %>% 
      spread(SampleID, value) %>%
      rename(SampleID = colnames(.)[1])
    COR1(env,genus,"RdYlGn",
         str_c(str_c(unique(Map$group),collapse = "-"),".wilcox."),
         opts$order,opts$rho,opts$width,opts$height,opts$mar,opts$tranfm)
  }
}

# 合并
if (opts$select != "none" & opts$test != "none"){
  st <- as.numeric(table(selectg$SampleID %in% test$SampleID))
  if (st[1] * st[2] != 0){
    genus <- genus %>% right_join(union(test[,1],selectg[,1])) %>% .[,-1] %>% t %>% 
      as.data.frame %>% rownames_to_column %>% as_tibble
    colnames(genus) <- c("SampleID",unlist(union(test[,1],selectg[,1])))
    
    COR1(env,genus,"RdYlGn",
         str_c(str_c(unique(Map$group),collapse = "-"),".union."),
         opts$order,opts$rho,opts$width,opts$height,opts$mar,opts$tranfm)
  }
}

write_tsv(opts %>% as_tibble,
          str_c("Parameter",
                str_replace_all(as.character(date())," ","_") %>% str_replace_all(":","_"),
                ".xls"),
          col_names = TRUE)

