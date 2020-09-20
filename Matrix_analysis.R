# Rscript Matrix_analysis.R -i bray_curtis_dm.txt -m F -g map-group.txt
# 基于矩阵的检验方法，不但可以输出检验显著性结果（p值），还有程度结果(R，R2，A)，可以用来判断分组贡献度大小
# adonis(permanova),anosim,mrpp
rm(list=ls())
library(magrittr)
library(tidyverse)
library(vegan)
library(optparse)
if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="bray_curtis_dm.txt",
                help="输入矩阵"),
    make_option(c("-m", "--mode"), type="logical", default=F,
                help="排列组合分析T,或常规分析F"),
    make_option(c("-g", "--group"), type="character", default="map-group.txt",
                help="分组文件"),
    make_option(c("-c", "--color"), type="character", default="",
                help="指定颜色color.txt"),
    make_option(c("-p", "--plot"), type="logical", default=F,
                help="是否画图")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
  print(paste("The input file is ", opts$input,  sep = ""))
  print(paste("combn or not mode is ", opts$mode,  sep = ""))
  print(paste("The group file is ", opts$group,  sep = ""))
  print(paste("The color file is ", opts$color,  sep = ""))
  print(paste("plot or not is ", opts$plot,  sep = ""))
}

mycol <-c('#B0C4DE',"#CD3333","#483D8B","#458B00","#EEC900","#BFBC3B","#6E8B3D","#00688B","#C10077","#CAAA76","#EEEE00","#458B00","#8B4513","#008B8B","#6E8B3D","#8B7D6B","#7FFF00","#CDBA96","#ADFF2F")
# read file
Map <- read.table(opts$group,header = T,sep = "\t",comment.char = "")
Matrix <- read.table(opts$input,header = T,row.names =1,sep = "\t") %>% .[as.vector(Map[,1]),as.vector(Map[,1])]

mpn <- nlevels(Map$group)
gp1 <- as.vector(unique(Map$group))
ADONIS_out <- c();ANOSIM_out <- c();MRPP_out <-c()

# 判断分析模式，排列组合还是指定最大分组
if (opts$mode){
  mpn1 <- 2;mpn2 <- mpn
}else {
  mpn1 <- mpn2 <- mpn
}

for (i in mpn1:mpn2){#i=9
  gp2 <- combn(gp1,i)
  for (j in 1:ncol(gp2)){#j=1
    mp1 <- Map[as.vector(Map$group) %in% gp2[,j],]
    Matrix1 <- Matrix[as.vector(mp1[,1]),as.vector(mp1[,1])]
    #如果有color.txt文件
    if(opts$color != ""){
      sc <- read.table(opts$color,sep="\t",comment.char = "",check.names = FALSE)
      sc <- sc[which(sc[,1] %in% unique(mp1$group)),]
      mycol <- c('#B0C4DE',as.vector(sc$V2))
      mp1[,2] <- factor(mp1$group,levels=as.vector(sc[,1]))
    }else{
      mp1[,2] <- factor(mp1$group,levels=as.vector(unique(mp1$group)))
    }
    Name1 <- paste(as.vector(unique(mp1$group)),collapse = "-")
    Name2 <- paste0(strsplit(opts$input,".txt")[[1]],".",Name1)
    set.seed(20190314)
    ADONIS <- adonis(Matrix1 ~ group, permutations = 9999,data = mp1)
    out_adonis <- as.data.frame(ADONIS$aov.tab)[1,]
    rownames(out_adonis) <- Name1
    ADONIS_out <- rbind(ADONIS_out,out_adonis)
    #p.adjusted <- p.adjust(ADONIS$aov.tab[1,6],method="fdr")
    
    set.seed(20190314)
    ANOSIM <- anosim(Matrix1,mp1$group, permutations = 9999)
    out_anosim <- t(as.matrix(c("test_statistic_name" = "R" ,"sample_size" = nrow(Matrix1),
                                "number_of_groups" = length(unique(mp1[,2])),
                                "R-value"	= ANOSIM$statistic,
                                "p-value"	= ANOSIM$signif,
                                "number_of_permutations" = ANOSIM$permutations))) 
    rownames(out_anosim) <- Name1
    ANOSIM_out <- rbind(ANOSIM_out,out_anosim)
    
    if (opts$plot == T){
      #画图
      ANOSIM$class.vec <- factor(ANOSIM$class.vec,levels=c("Between",as.vector(levels(mp1$group))))
      pdf(paste(Name2,'.anosim.pdf',sep = ''))
      plot(ANOSIM,col=mycol,ylab="Rank Dissimilarity",main="ANOSIM")
      dev.off()
    }
    
    set.seed(20190314)
    MRPP <- mrpp(Matrix1,mp1$group, permutations = 9999)
    out_mrpp <- t(as.matrix(c("Observed_delta" = MRPP$delta,
                              "Expected_delta" = MRPP$E.delta,
                              "A" = MRPP$A,#A = 1 -delta/E(delta)
                              "Significance_of_delta" = MRPP$Pvalue,
                              "Number_of_permutations" = MRPP$permutations)))
    # E(delta) is the expected delta assessed as the average of permutations
    # A > 0，表示组间差异大于组内差异
    # A < 0，表示组内差异大于组间差异
    rownames(out_mrpp) <- Name1
    MRPP_out <- rbind(MRPP_out,out_mrpp)
    if (opts$plot == T){
      # 画图
      md <- meandist(as.dist(Matrix1), mp1$group)
      pdf(paste(Name2,'.mrpp.pdf',sep = ''))
      plot(md)
      dev.off()
    }
  }
}

# 显著性标记
ROWNAME <- rownames(MRPP_out)

MRPP_out %<>% as.data.frame(.)
MRPP_out$`Significance_of_delta` <- as.numeric(as.vector(MRPP_out$`Significance_of_delta`))
MRPP_out %<>% mutate(Sig_mark = ifelse(`Significance_of_delta`>0.05,"",
                                                           ifelse(`Significance_of_delta`>0.01,"*",
                                                                  ifelse(`Significance_of_delta`>0.001,"**","***"))))

ANOSIM_out %<>% as.data.frame(.)
ANOSIM_out$`p-value` <- as.numeric(as.vector(ANOSIM_out$`p-value`))
ANOSIM_out %<>% mutate(Sig_mark = ifelse(`p-value`>0.05,"",
                                         ifelse(`p-value`>0.01,"*",
                                                ifelse(`p-value`>0.001,"**","***"))))

ADONIS_out$`Pr(>F)` <- as.numeric(as.vector(ADONIS_out$`Pr(>F)`))
ADONIS_out %<>% as.data.frame(.) %>% mutate(Sig_mark = ifelse(`Pr(>F)`>0.05,"",
                                                             ifelse(`Pr(>F)`>0.01,"*",
                                                                    ifelse(`Pr(>F)`>0.001,"**","***"))))

write.table(file=paste0(Name2,'.mrpp.xls'),cbind(Groups = ROWNAME,MRPP_out),
            sep="\t",quote = FALSE,row.names = F,col.names = T)
write.table(file=paste0(Name2,'.anosim.xls'),cbind(Groups = ROWNAME,ANOSIM_out),
            sep="\t",quote = FALSE,row.names = F,col.names = T)
write.table(file=paste0(Name2,'.adonis.xls'),cbind(Groups = ROWNAME,ADONIS_out),
            sep="\t",quote = FALSE,row.names = F,col.names = T)

#library(formattable)
#
#P_formatter <- 
#  formatter("span", 
#            style = x ~ style(
#              font.weight = "bold", 
#              color = ifelse(x <0.05, customGreen, "black")))
#
#formattable(ADONIS_out,align =c("l","c","c","c","c","r"),
#            `Indicator Name` = formatter("span", 
#  style = ~ style(color = "grey",font.weight = "bold"),
#  `R2` = color_bar(customRed),
#  `PR(>F)` = P_formatter))
