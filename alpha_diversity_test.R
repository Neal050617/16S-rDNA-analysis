# 支持多组排列组合
# 支持成对或不成对 选择
# 支持显著性标记
# Rscript alpha_diversity_test.R -a alpha_rarefac.summary.xls -g map-group.txt -t FALSE

library(optparse)
library(tidyverse)
library(magrittr)
library(reshape2)
library(gtools)

if (TRUE){
  option_list <- list(
    make_option(c("-a", "--alpha"), type="character", default="alpha_rarefac.summary.xls",
                help="输入的OTU表格"),
    make_option(c("-g", "--map"), type="character", default="map-group.txt",
                help="分组文件"),
    make_option(c("-t", "--test"), type="character", default="FALSE",
                help="成对检验TRUE；非成对检验（曼惠特尼检验）FALSE")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
}

Base_Statistics <- function(TT,ts){# TT <- mp4[,c(4,7)]
  tt <- length(unique(as.vector(TT[,2])))
  colnames(TT)[1] <- c("data")
  TT %<>% as.data.frame
  set.seed(20190306)
  if (tt == 2){
    if (ts == "FALSE"){
      TEST <- wilcox.test(data ~ group,TT,paired=FALSE,exact = TRUE)
    }else if (ts == "TRUE"){
      TEST <- wilcox.test(data ~ group,TT,paired=TRUE,exact = TRUE)
    }
  }else if (tt > 2){
    TEST <- kruskal.test(data ~ group,data = TT)
  }
  return(TEST$p.value)
}

pp <- c(0.25, 0.5, 0.75)
p_names <- map_chr(pp, ~paste0("IQR",.*100, "%"))
p_funs <- map(pp, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
  set_names(nm = p_names)

Minus <- function(x,n){
  d1 <- 10^(-n)
  ifelse(x >= d1,as.vector(round(x,n)),paste0("<",d1))
}

################################################## 读取文件 #################################################################

mp1 <- read.table(opts$map,header=T,sep="\t",comment.char = "") %>%
  rename("SampleID" = colnames(.)[1])
mpn <- nlevels(mp1[,2])
gp1 <- as.vector(unique(mp1$group))

data <- read.table(opts$alpha,header=T,sep="\t")  %>%
  rename("SampleID" = colnames(.)[1]) %>%
  right_join(.,mp1) %>%
  select(-group)

mp2 <- vector(mpn-1,mode="list")
for (i in 2:mpn){# i=2
  gp2 <- combn(gp1,i)
  mp2[[i-1]] <- vector(ncol(gp2),mode="list")
  for (j in 1:ncol(gp2)){# j=2
    mp3 <- mp1[mp1$group %in% gp2[,j],] # 筛选出group文件
    mp4 <- data[data[,1] %in% as.vector(mp3[,1]),c(1,4,5,6,7,9)] %>% 
      left_join(.,mp3) %>% 
      mutate(group = factor(group,levels=unique(as.vector(group)))) %>%
      arrange(group)

    gp <- as.vector(unique(unlist(mp4[,c("group")])))
    gpn <- unlist(lapply(1:length(gp),function(g){# g = 1
      unlist(lapply(1:3,function(p){paste0(gp[g],"-",c("median","mean","se")[p])}))
    }))
    # 计算均值和标准误
    mp6 <- (lapply(colnames(mp4)[2:6],function(x){# x = colnames(mp4)[2]
      mp5 <- mp4[,c("SampleID",x,"group")] %>% group_by(group) %>% 
        summarise_each(funs(mean, se = sd(.)/sqrt(n()),!!!p_funs),x) %>% 
      mutate(IQR = sapply(seq_along(.$group),function(x){
        str_c(Minus(.$`IQR50%`[x],4),"(",
              Minus(.$`IQR25%`[x],4),",",
                        Minus(.$`IQR75%`[x],4),")")
      })) %>% select_at(vars(-contains("%"))) %>%
        dplyr::select(-group) %>% unlist(.) %>% 
        .[order(as.numeric(gsub("([A-Za-z]+)([0-9]+)", "\\2", names(.))),
                              gsub("([A-Za-z]+)([0-9]+)", "\\1", names(.)))]
      names(mp5) <- gpn
      return(mp5)
    }))
    names(mp6) <- colnames(mp4)[2:6]
    # 计算p-value、q-value ;合并表格
    mp7 <- as.data.frame(t(data.frame(mp6))) %>% 
      mutate(`p-value` = unlist(lapply(2:6,function(m){
        Base_Statistics(mp4[,c(m,7)],opts$test)}))) %>%
      mutate(`z-score` = qnorm(`p-value`/2)) %>%
      mutate(Sig_mark = ifelse(`p-value`>0.05,"",
                               ifelse(`p-value`>0.01,"*",
                                      ifelse(`p-value`>0.001,"**","***"))))
    # %>% mutate(`q-value`= p.adjust(as.numeric(.[,5]),lambda=0.05)$qvalues)
    rownames(mp7) <- colnames(mp4)[2:6]
    mp2[[i-1]][[j]] <- mp7
  }
}

# 写入文件
for (ll in 1:length(mp2)){# ll <- 1
  nm <- paste0("map",as.character(ll+1),".alpha_rarefac.",
           ifelse(ll==1,"Mann-Whitney","Kruskal-Wallis_rank_sum_test"),".xls")
  if (file.exists(nm)) # 删除文件
    file.remove(nm)
  for (lll in 1:length(mp2[[ll]])){# lll <- 1
    write.table(file=nm,
                cbind(Alpha_diversity=rownames(mp2[[ll]][lll][[1]]),
                              mp2[[ll]][lll][[1]]),
                sep="\t",quote = FALSE,append = TRUE,row.names = F)
  }
}

#if(!require(coin)){install.packages("coin")}
#set.seed(169)
#n = 16
#N = n + n
#A = runif(n,1,10)
#B = A + 5 + rnorm(n,0,4)
#Group = factor(c(rep("A", length(A)), rep("B", length(B))))
#Y = c(A,B)
#boxplot(Y ~ Group)
#
#library(coin)
#MWa = wilcox.test(Y ~ Group, exact=FALSE)
#MWa
#Za = qnorm(MWa$p.value/2)
#Za
#
#qnorm(MWa$p.value/2)
#
#ra = abs(Za)/sqrt(N)
#names(ra) = "ra"
#ra
