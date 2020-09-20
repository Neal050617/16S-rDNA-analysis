library(tidyverse)
library(optparse)
library(magrittr)
library(pROC)
library(broom)

if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="Part1-rarefac.otu_genus.xls",
                help="otu表格而已，带有分类学信息的表格最好别用"),
    make_option(c("-l", "--lda"), type="character", default="LDA.xls",
                help="otu表格而已，带有分类学信息的表格最好别用"),
    make_option(c("-m", "--map"), type="character", default="Part.list",
                help="分组数据"),
    make_option(c("-d", "--mode"), type="numeric", default="1",
                help="测试集1还是验证集2"),
    make_option(c("-a", "--num"), type="numeric", default=2,
                help="lda阈值")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
}

LDA <- read_tsv(opts$lda) %>% 
  mutate(P_value = as.numeric(P_value)) %>% 
  filter(LDA_value>=opts$num,P_value <= 0.05) %>%
  .[,1] %>% write_tsv("lda.select.xls")

data <- read_tsv(opts$input,col_names = F) %>% 
  right_join((LDA[,1] %>% 
               add_row(.,.before = 1,Biomaker_names="class")),
             by=c("X1"="Biomaker_names"))

Map <- read_tsv(opts$map,col_names = FALSE) %>% 
  filter(X2 == opts$mode) %>% select(-X2) %>%
  rename("SampleID"=colnames(.)[1],"group"=colnames(.)[2]) %>% 
  mutate(group=fct_inorder(group))

#分组
out <- vector(nlevels(Map$group),mode="list")
for (i in 1:length(out)){# i=1
  out[[i]] <- data[-1,c(TRUE,which(data[1,] %in% levels(Map$group)[i]))]
  out[[i]]$Rowsum <- sapply(1:nrow(out[[i]]),function(x) # x <- 1
    sum(as.numeric(out[[i]][x,2:ncol(out[[i]])])))
}

#分类
v <- c()
for (j in 1:nrow(LDA)){# j=1
  if (out[[1]]$Rowsum[j] >= out[[2]]$Rowsum[j]){
    v <- c(v,levels(Map$group)[1])
  }else{
    v <- c(v,levels(Map$group)[2])
  }
}
LDA %<>% mutate(group = v)

# MAI指数
MAI <- c()
for (k in 1:length(out)){# k = 1
  MAIsum <- sapply(2:ncol(out[[k]]),function(x)  # x = 2
    sum(as.numeric(unlist(out[[k]][LDA$group %in% levels(Map$group)[1],x])))-
      sum(as.numeric(unlist(out[[k]][LDA$group %in% levels(Map$group)[2],x]))))
  MAI <- c(MAI,MAIsum[-length(MAIsum)])
}
MAI <- MAI %>% as.data.frame(.)
colnames(MAI) <- c("MAI")
rownames(MAI) <- colnames(LDA$Biomaker_names)
MAI$group <- Map$group

# train ROC
pdf(str_c(ifelse(opts$mode==1,"train","test"),".ROC.pdf"))
roc <- roc(MAI$group,MAI$MAI,plot=T,col="black",ci=F,auc.polygon=F,print.thres=F,print.auc=F,percent=F,xlab="Specificity(%)",ylab="Sensitivity(%)")
# 秩和检验
set.seed(20190618)
U <- wilcox.test(MAI[MAI$group %in% levels(Map$group)[1],1],
                 MAI[MAI$group %in% levels(Map$group)[2],1])
U %>% tidy %>% write.csv(file = str_c(ifelse(opts$mode==1,"train","test"),
                                      ".wilcox.ROC.txt"))
# 置信区间
sens.ci <- ci.se(roc)
plot(sens.ci, type="s", col="lightblue",border="white")
# 标注
legend("bottomright",bty="n",paste0("AUC: ",round(ci(roc)[2],4),"\n","95% CI: ",round(ci(roc)[1],4),"-",round(ci(roc)[3],4),"\n","p-value: ",round(U$p.value,4),"\n"))
tt <- paste0(round(coords(roc,"best")[1],4),"\n","(",round(coords(roc,"best")[2],4),",",round(coords(roc,"best")[3],4),")")
points(coords(roc,"best")[2],coords(roc,"best")[3],pch=16,col="red",cex=1.5,font=1)
text(coords(roc,"best")[2]-0.1,coords(roc,"best")[3]-0.1,tt,cex=1.5,pos=4,col="black")
dev.off()
