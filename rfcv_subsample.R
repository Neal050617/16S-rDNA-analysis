library(tidyverse)
library(magrittr)
library(randomForest)
library(pROC)
library(optparse)
library(broom)
library(patchwork)

if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="rarefac.otu_genus.xls",
                help="otu表格而已，带有分类学信息的表格最好别用"),
    make_option(c("-w", "--wilcox"), type="character", default="wilcox.otu.AD-CN.xls",
                help="wilcox检验"),
    make_option(c("-m", "--map"), type="character", default="map-group.txt",
                help="otu表格而已，带有分类学信息的表格最好别用"),
    make_option(c("-n", "--cvnumber"), type="numeric", default=5,
                help="n倍的交叉验证"),
    make_option(c("-t", "--cvtime"), type="numeric", default=5,
                help="重复几遍n倍交叉验证"),
    make_option(c("-p", "--per"), type="numeric", default=0.005,
                help="OTU丰度筛选"),
    make_option(c("-k", "--markernum"), type="numeric", default=0,
                help="挑选的标志物个数"),
    make_option(c("-f", "--feature"), type="character", default="none",
                help="feature_importance_scores.txt"),
    make_option(c("-c", "--color"), type="character", default="none",
                help="设置颜色"),
    make_option(c("-s", "--seed"), type="character", default="1234",
                help="设置种子,默认1234"),
    make_option(c("-x", "--Meande"), type="double", default=0.0001,
                help="feature_importance_scores 阈值"),
    make_option(c("-y", "--pvalue"), type="double", default=0.05,
                help="pvalue"),
    make_option(c("-z", "--qvalue"), type="double", default=0,
                help="qvalue"),
    make_option(c("-v", "--valid"), type="character", default="valid.rarefac.otu_table.xls",
                help="验证集数据"),
    make_option(c("-g", "--gp"), type="character", default="none",
                help="验证集数据分组map2.txt")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
}

i <- opts$input
w <- opts$wilcox
m <- opts$map
cvn <- opts$cvnumber
cvt <- opts$cvtime
p <- opts$per
f <- opts$feature
c <- opts$color
marker.num <- opts$markernum

options(scipen = 200)

## 设置随机数种子，保证结果可重复
set.seed(opts$seed) # set.seed(2345) # set.seed(3456)

# 修改后的rfcv函数
my.rfcv <- function (trainx, trainy, cv.fold = 5, scale = "log", step = 0.5, 
                     mtry = function(p) max(1, floor(sqrt(p))), recursive = FALSE, ...) 
{
  classRF <- is.factor(trainy)
  n <- nrow(trainx)
  p <- ncol(trainx)
  if (scale == "log") {
    k <- floor(log(p, base = 1/step))
    n.var <- round(p * step^(0:(k - 1)))
    same <- diff(n.var) == 0
    if (any(same)) 
      n.var <- n.var[-which(same)]
    if (!1 %in% n.var) 
      n.var <- c(n.var, 1)
  }
  else {
    n.var <- seq(from = p, to = 1, by = step)
  }
  k <- length(n.var)
  cv.pred <- vector(k, mode = "list")
  for (i in 1:k) cv.pred[[i]] <- trainy
  if (classRF) {
    f <- trainy
  }
  else {
    f <- factor(rep(1:5, length = length(trainy))[order(order(trainy))])
  }
  nlvl <- table(f)
  idx <- numeric(n)
  for (i in 1:length(nlvl)) {
    idx[which(f == levels(f)[i])] <- sample(rep(1:cv.fold, length = nlvl[i]))
  }
  
  res = list() 
  
  for (i in 1:cv.fold) {
    all.rf <- randomForest(trainx[idx != i, , drop = FALSE], 
                           trainy[idx != i], trainx[idx == i, , drop = FALSE], 
                           trainy[idx == i], mtry = mtry(p), importance = TRUE, ...)
    cv.pred[[1]][idx == i] <- all.rf$test$predicted
    #impvar <- (1:p)[order(all.rf$importance[, 1], decreasing = TRUE)]
    impvar <- (1:p)[order(randomForest::importance(all.rf, type = 1), decreasing = TRUE)]
    res[[i]] <- impvar
    for (j in 2:k) {
      imp.idx <- impvar[1:n.var[j]]
      sub.rf <- randomForest(trainx[idx != i, imp.idx, 
                                    drop = FALSE], trainy[idx != i], trainx[idx == i, imp.idx, drop = FALSE],
                             trainy[idx == i], mtry = mtry(n.var[j]), importance = recursive, ...)
      cv.pred[[j]][idx == i] <- sub.rf$test$predicted
      if (recursive) {
        #impvar <- (1:length(imp.idx))[order(sub.rf$importance[, 1], decreasing = TRUE)]
        impvar <- (1:length(imp.idx))[order(randomForest::importance(sub.rf, type = 1), decreasing = TRUE)]
      }
      NULL
    }
    NULL
  }
  if (classRF) {
    error.cv <- sapply(cv.pred, function(x) mean(trainy != x))
  }
  else {
    error.cv <- sapply(cv.pred, function(x) mean((trainy - x)^2))
  }
  names(error.cv) <- names(cv.pred) <- n.var
  list(n.var = n.var, error.cv = error.cv, predicted = cv.pred, res = res)
}

Minus <- function(x,n){
  d1 <- 10^(-n)
  ifelse(x >= d1,round(x,4),paste0("< ",d1))
}

data1 <- read.table(m,head= T ,sep="\t",comment.char = "",fileEncoding = "UTF-8")
colnames(data1) <- c("SampleID","group")

data2 <- read.table(i,head= T ,sep="\t",comment.char = "",row.names = 1,fileEncoding = "UTF-8")
data2_1 <- data2[,data1$SampleID]
data2_1 <- sapply(1:ncol(data2_1),function(x) data2_1[,x] <- data2_1[,x]/sum(data2_1[,x]))
rownames(data2_1) <- rownames(data2)
colnames(data2_1) <- data1$SampleID
d2p <- sapply(1:nrow(data2_1),function(y) any(data2_1[y,]>=p))	 
data2_2 <- data2_1[d2p,]

data3 <- read.table(w,head= T ,sep="\t",comment.char = "",row.names = 1,fileEncoding = "UTF-8")
d3p <- rownames(data3)[data3$p.value <= opts$pvalue]
data2_3 <- data2_2[rownames(data2_2) %in% d3p,]
if (opts$qvalue != 0){
  d3q <- rownames(data3)[data3$q.value <= opts$qvalue]
  data2_3 <- data2_3[rownames(data2_3) %in% d3q,]
}
Data <- cbind(t(data2_3),data1)

if (f != "none"){
  data4 <- read.table(f,head= T ,sep="\t",comment.char = "",row.names = 1,fileEncoding = "UTF-8")
  data4 <- data4[data4$Mean_decrease_in_accuracy >= opts$Meande,]
  data2_4 <- sapply(rownames(data2_3),function(x) strsplit(x," ", fixed=TRUE)[[1]][1])
  data2_5 <- data2_3[data2_4 %in% rownames(data4),]
  nrow(data2_5)
  Data <- cbind(t(data2_5),data1)
}

Data <- Data[,-(ncol(Data)-1)]

#如果有颜色定义的话，就可以直接定义了
if (c != "none"){
  sc <- read.table(c,sep="\t",comment.char = "",check.names = FALSE)
  sc <- sc[which(sc[,1] %in% unique(Data$group)),]
  mycol <- as.vector(sc[,2])
  Data$group <- factor(Data$group,levels=as.vector(sc[,1]))
}else{
  Data$group <- factor(Data$group,levels=unique(as.vector(Data$group)))
}
Data <- dplyr::arrange(Data,group)

## 设置随机数种子，保证结果可重复
#set.seed(1234) # set.seed(2345) # set.seed(3456)
## 随机选取2/3预测，1/3验证
##随机1-2有放回抽样150次，概率分别为0.67和0.33，用于获取训练集
ind=sample(2,nrow(Data),replace = TRUE,prob = c(0.67,0.33))
data1 %>% mutate(Split = ind) %>% write_tsv("split.map-group.txt")
# 2/3寻找最优
result <- replicate(cvt, my.rfcv(Data[ind==1,-ncol(Data)], Data[ind==1,"group"],cv.fold = cvn,step=0.5), simplify=FALSE)
error.cv <- sapply(result, "[[", "error.cv")

error.cv.rm <- rowMeans(error.cv)
id <- error.cv.rm < min(error.cv.rm) + sd(error.cv.rm)
error.cv[id, ]
if (marker.num == 0) { 
  marker.num <- min(as.numeric(names(error.cv.rm)[id]))
}

prefix <- paste(unique(as.vector(data1$group)),collapse = "-")
pdf.dir1=paste0(prefix, "_vars.pdf")
pdf.dir2=paste0(prefix, "_train_boxplot.pdf") 
pdf.dir3=paste0(prefix, "_train_roc.pdf")
pdf.dir4=paste0(prefix, "_test_boxplot.pdf") 
pdf.dir5=paste0(prefix, "_test_roc.pdf")

pdf(pdf.dir1) 
matplot(result[[1]]$n.var, error.cv, type = "l", log = "x", col = rep(1, cvt),
        #main = paste("select", marker.num, "Vars"),
        xlab = "Number of vars",  
        ylab = "CV Error", lty = 1) 
lines(result[[1]]$n.var, error.cv.rm, lwd = 2)
text(x = marker.num,y = min(error.cv),labels = marker.num,adj = c(1, 0.5),cex = 1.2, xpd = TRUE, font = 2, col = "pink")
abline(v = marker.num, col = "pink", lwd = 2)
dev.off()

# pick marker by corossvalidation 
marker.t <- table(unlist(lapply(result, function(x) { 
  lapply(x$res, "[", 1:marker.num) 
}))) 
marker.t <- sort(marker.t, d = T) 
names(marker.t) <- colnames(Data)[as.numeric(names(marker.t))] 
marker.dir <- paste0(prefix, "_marker.txt") 
write.table(marker.t, marker.dir, col.names = F, sep = "\t", quote = F) 
marker.p <- names(marker.t)[1:marker.num] 

# 2/3作预测建模
set.seed(0) # set.seed(10) # set.seed(100)
train.rf=randomForest(Data[ind==1,marker.p],Data[ind==1,"group"],ntree=1000,proximity=TRUE,importance=TRUE)
train.p <- predict(train.rf, type = "prob")
pdf(pdf.dir2)
if (c != "none"){
  boxplot(train.p[, 2] ~ Data[ind==1,"group"], col = mycol, main = "Probability", 
          names = levels(Data[,"group"]), xlab = "Groups", ylab = "POD")
}else{
  boxplot(train.p[, 2] ~ Data[ind==1,"group"], col = 3:2, main = "Probability", 
          names = levels(Data[,"group"]), xlab = "Groups", ylab = "POD") 
}
dev.off() 
pr.dir <- paste0(prefix, "_train_probability.txt") 
write.table(train.p[, 2], pr.dir, sep = "\t", quote = F, col.names = F)

# varImPlot::MeanDecreaseAccuracy(平均减少度）,即没有这个Feature，分类准确度下降的程度，相当于常用的分类贡献度的概念。
#pdf(prefix, "_train_varImPlot.pdf")
#varImpPlot(train.rf)
#dev.off()
bbtheme <- function(){
  return(theme_bw() + ##设置主题
           theme(plot.title=element_text(size=rel(1),hjust=0.5),
                 plot.margin = unit(c(1, 1, 1, 1), "lines"),
                 panel.grid.major=element_line(color="white"),
                 panel.grid.minor=element_line(color="white"),
                 axis.title=element_text(size=rel(1))#,
                 #axis.text.x=element_text(angle=30,hjust =1),
                 #legend.title=element_blank(),
                 #legend.text = element_text(size = 6),
                 #legend.key.size = unit(.4,'cm'),
                 #legend.spacing.x = unit(0.1, 'cm')
           ))
}

WHat <- varImpPlot(train.rf) %>% as.data.frame %>% 
  rownames_to_column() %>% as_tibble %>%
  rename(OTU_genus = colnames(.)[1]) %>%
  gather(value,ID,-OTU_genus) %>%
  group_by(value) %>% nest %>% 
  mutate(data = map(data,function(x) 
    x %>% arrange(ID) %>%
      mutate(OTU_genus = fct_inorder(OTU_genus))))

GGplot <- function(x,yb){
  ggplot(x,aes(x=OTU_genus,y=ID)) +
    geom_point(size=3,color="#FF4500") +
    #geom_segment(aes(x=OTU_genus,xend=OTU_genus,y=0,yend=MeanDecreaseAccuracy))+
    #bbtheme() + #geom_line() +
    theme_bw() +
    coord_flip() + xlab('') + ylab(yb) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(colour="grey60", linetype="dashed"),
          axis.text.y = element_text(hjust = 0))
}

Pt <- map(seq_along(WHat),function(x) GGplot(WHat$data[[x]],WHat$value[[x]]))

ggsave(str_c(prefix,"_varImpPlot.pdf"),Pt[[1]] + Pt[[2]],
       width = 12,
       height = log2(nrow(WHat$data[[1]]))*1.2+3
)

data2_3 %>% as.data.frame %>% rownames_to_column() %>% as_tibble %>% 
  rename(`OTU ID` = colnames(.)[1]) %>%
  filter(`OTU ID` %in% marker.p) %>%
  write_tsv(str_c("marker.p.",opts$input))

# train ROC
pdf(pdf.dir3) 
roc <- roc(Data[ind==1,"group"],train.p[, 2],plot=T,col="black",ci=F,auc.polygon=F,print.thres=F,print.auc=F,percent=F,xlab="Specificity(%)",ylab="Sensitivity(%)")
# 秩和检验
U <- wilcox.test(train.p[Data[ind==1,"group"] %in% levels(Data[ind==1,"group"])[1],2],
                 train.p[Data[ind==1,"group"] %in% levels(Data[ind==1,"group"])[2],2])
U %>% tidy %>% write.csv(file = "wilcox.train.ROC.txt")
# 置信区间
sens.ci <- ci.se(roc)
plot(sens.ci, type="s", col="lightblue",border="white")
# 标注
legend("bottomright",bty="n",paste0("AUC: ",round(ci(roc)[2],4),"\n","95% CI: ",round(ci(roc)[1],4),"-",round(ci(roc)[3],4),"\n","p-value: ",Minus(U$p.value,4),"\n"," "))
tt <- paste0(round(coords(roc,"best")[1],4),"\n","(",round(coords(roc,"best")[2],4),",",round(coords(roc,"best")[3],4),")")
points(coords(roc,"best")[2],coords(roc,"best")[3],pch=16,col="red",cex=1.5,font=1)
text(coords(roc,"best")[2]-0.1,coords(roc,"best")[3]-0.1,tt,cex=1.5,pos=4,col="black")
dev.off()

#stats <- wilcox.test(pred[obs == 1], pred[obs == 0], alternative = "great")

## Using different symbols for the classes:
#MDSplot(train.rf, Data[ind==1,"group"])
#plot(margin(train.rf))

# ## 分类结果主坐轴分析和可视化 Do PCoA/MDS on 1-proximity:
# train.mds = cmdscale(1-train.rf$proximity,eig=TRUE)
# ## 设置显示样本点，而不是变量点
# op=par(pty="s")
# ## 散点图展示每个Feature与PCoA1/2轴间的相互分布
# pairs(cbind(Data[ind==1,marker.p],train.mds$points),cex=0.6,gap=0,col=c("red","green")[as.numeric(Data[ind# ==1,"group"])],main="Predictors and MDS of Proximity Based on RandomForest")
# par(op)
# print(train.mds$GOF)

if (opts$gp == "none"){
  # test predict 
  # 1/3作预测模型
  test.p <- predict(train.rf, Data[ind==2,-ncol(Data)], type = "prob") 
  pr.dir <- paste0(prefix, "_test_probability.txt") 
  write.table(test.p[,2], pr.dir, sep = "\t", quote = F, col.names = F)
  #test.result=cbind(Data[ind==2,"group"],test.p[,2]) 
  #write.table(test.result, pr.dir, sep = "\t", quote = F, col.names = F)
  
  pdf(pdf.dir4) 
  if (c != "none"){
    boxplot(test.p[, 2] ~ Data[ind==2,"group"], col = mycol, main = "Probability", 
            names = levels(Data[,"group"]), xlab = "Groups", ylab = "POD")
  }else{
    boxplot(test.p[, 2] ~ Data[ind==2,"group"], col = 3:2, main = "Probability", 
            names = levels(Data[,"group"]), xlab = "Groups", ylab = "POD") 
  }
  dev.off() 
  
  # predict plot 
  #p.col <- ifelse(is.na(Data[ind==2,-ncol(Data)]), 4, as.numeric(Data[ind==2,-ncol(Data)]) + 1)
  #plot(rank(test.p[, 2]), test.p[, 2], col = 4, pch = 16, xlab = "", ylab = "Probability", main = "Testset") 
  #txt <- levels(Data[ind==1,"group"])
  #if (length(levels(Data[ind==1,"group"])) > 2) { 
  #  txt <- c(txt, "the rest") 
  #}
  #legend("bottomright", txt, col = 2:4, pch = 16) 
  #abline(h = 0.5) 
  
  # test ROC 
  pdf(pdf.dir5)
  roc <- roc(Data[ind==2,"group"],test.p[, 2],plot=T,col="black",ci=F,auc.polygon=F,print.thres=F,print.auc=F,percent=F,xlab="Specificity(%)",ylab="Sensitivity(%)")
  # 秩和检验
  U <- wilcox.test(test.p[Data[ind==2,"group"] %in% levels(Data[ind==2,"group"])[1],2],
                   test.p[Data[ind==2,"group"] %in% levels(Data[ind==2,"group"])[2],2])
  U %>% tidy %>% write.csv(file = "wilcox.test.ROC.txt")
  # 置信区间
  sens.ci <- ci.se(roc)
  plot(sens.ci, type="s", col="lightblue",border="white")
  # 标注
  legend("bottomright",bty="n",paste0("AUC: ",round(ci(roc)[2],4),"\n","95% CI: ",round(ci(roc)[1],4),"-",round(ci(roc)[3],4),"\n","p-value: ",Minus(U$p.value,4),"\n"," "))
  tt <- paste0(round(coords(roc,"best")[1],4),"\n","(",round(coords(roc,"best")[2],4),",",round(coords(roc,"best")[3],4),")")
  points(coords(roc,"best")[2],coords(roc,"best")[3],pch=16,col="red",cex=1.5,font=1)
  text(coords(roc,"best")[2]-0.1,coords(roc,"best")[3]-0.1,tt,cex=1.5,pos=4,col="black")
  dev.off()
  
}else{
  ## 验证
  testgp <- read_tsv(opts$gp) %>% 
    rename(SampleID = colnames(.)[1],group = colnames(.)[2]) %>%
    mutate(group = factor(.$group,levels=levels(data1$group)))
  
  testo <- marker.p %>% enframe %>% 
    mutate(backup = value) %>%
    separate(backup,into = c("SampleID","genus"),sep = " ") %>%
    .[,c("SampleID","genus")]
  
  test <- read_tsv(opts$valid) %>% 
    rename("SampleID" = colnames(.)[1]) %>%
    .[,c("SampleID",as.vector(testgp$SampleID))] %>% 
    mutate_at(vars(2:ncol(.)),function(x) x/sum(x)) %>%
    right_join(testo) %>% 
    replace(., is.na(.), 0) %>% 
    mutate(SampleID = str_c(SampleID,genus,sep = " ")) %>%
    mutate(SampleID = factor(SampleID, 
                             levels=unique(SampleID))) %>%
    select(-genus) %>%
    gather(var, value, -SampleID) %>% 
    mutate(var = factor(var,levels=unique(var))) %>%
    spread(SampleID, value) %>% 
    rename(SampleID = var) %>%
    right_join(.,testgp)
  
  test$group <- factor(test$group,levels = levels(Data$group))
  write_tsv(test %>% select(-group),"test_output.xls")
  cc <- test$SampleID
  Test <- test %>% as.data.frame %>% .[,-1]
  rownames(Test) <- cc
  
  # test predict 
  # 1/3作预测模型
  test.p <- predict(train.rf, Test[,-ncol(Test)], type = "prob")
  pr.dir <- paste0(prefix, "_test_probability.txt") 
  write.table(test.p[,2], pr.dir, sep = "\t", quote = F, col.names = F)
  
  pdf(pdf.dir4) 
  if (c != "none"){
    boxplot(test.p[, 2] ~ Test[,"group"], col = mycol, main = "Probability", 
            names = levels(Data[,"group"]), xlab = "Groups", ylab = "POD")
  }else{
    boxplot(test.p[, 2] ~ Test[,"group"], col = 3:2, main = "Probability", 
            names = levels(Data[,"group"]), xlab = "Groups", ylab = "POD") 
  }
  dev.off() 
  
  # test ROC 
  pdf(pdf.dir5)
  roc <- roc(Test[,"group"],test.p[, 2],plot=T,col="black",ci=F,auc.polygon=F,
             print.thres=F,print.auc=F,percent=F,xlab="Specificity(%)",ylab="Sensitivity(%)")
  # 秩和检验
  U <- wilcox.test(test.p[Test[,"group"] %in% levels(Test[,"group"])[1],2],
                   test.p[Test[,"group"] %in% levels(Test[,"group"])[2],2])
  U %>% tidy %>% write.csv(file = "wilcox.test.ROC.txt")
  # 置信区间
  sens.ci <- ci.se(roc)
  plot(sens.ci, type="s", col="lightblue",border="white")
  # 标注
  legend("bottomright",bty="n",paste0("AUC: ",round(ci(roc)[2],4),"\n","95% CI: ",round(ci(roc)[1],4),"-",round(ci(roc)[3],4),"\n","p-value: ",Minus(U$p.value,4),"\n"," "))
  tt <- paste0(round(coords(roc,"best")[1],4),"\n","(",round(coords(roc,"best")[2],4),",",round(coords(roc,"best")[3],4),")")
  points(coords(roc,"best")[2],coords(roc,"best")[3],pch=16,col="red",cex=1.5,font=1)
  text(coords(roc,"best")[2]-0.1,coords(roc,"best")[3]-0.1,tt,cex=1.5,pos=4,col="black")
  dev.off()
  
}

##################################################################################################

write_tsv(opts %>% as_tibble,
          str_c("Parameter",
                str_replace_all(as.character(date())," ","_") %>% str_replace_all(":","_"),
                ".xls"),
          col_names = TRUE)
