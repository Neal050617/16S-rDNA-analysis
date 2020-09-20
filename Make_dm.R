# Rscript make_dm.R -i rarefac.otu_table.xls -m hellinger -p otu_reps_aligned.fasta.tre
#Useage: 代替qiime脚本生成(WU)unifric_dm、bray_curtis_dm、jaccard_dm
#
rm(list=ls())
#
package_list <- c("picante","GUniFrac","ade4","ape","vegan","reshape2","ggplot2","dplyr","optparse","magrittr","exactRankTests","nlme")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p,repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
#
if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="rarefac.otu_table.xls",
                help="otu表格而已，带有分类学信息的表格最好别用"),
    make_option(c("-m", "--method"), type="character", default="hellinger",
                help="c('none','total','max','normalize','standardize','chi.square','hellinger','log')"),
    make_option(c("-p", "--tree"), type="character", default="otu_reps_aligned.fasta.tre",
                help="树文件为构建系统发育树所输出的文件，格式为.tre 、.nwk等有根进化树，OTU个数必须与otu表格中otu的个数一致")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
  print(paste("The input file is ", opts$input,  sep = ""))
  print(paste("The input tree is ", opts$tree,  sep = ""))
}

# read file
otu_tab <- read.table(opts$input,header = T,row.names =1,sep = "\t")
if(all(substr(colnames(otu_tab),1,1)=="X")){
  colnames(otu_tab) <- substr(colnames(otu_tab),2,nchar(colnames(otu_tab)))
}

set.seed(20190313)
# 矩阵转化
if (opts$method!="none"){
  otu_tabt <- t(decostand(otu_tab,opts$method))
}else{
  otu_tabt <- t(otu_tab)
}

## otu表格抽平，因为本项目的OTU表格默认为抽平至1000000，所以第一步跳过；
## Otu_tab_rff <- Rarefy(t(otu_tab))$otu.tab.rff   
#
#计算bray-curtis 与 jaccard 距离矩阵
bc <- as.matrix(vegdist(otu_tabt,method='bray')) #bray-curtis
jc <- as.matrix(vegdist(otu_tabt,method='jaccard',binary=TRUE)) #jaccard
ec <- as.matrix(vegdist(otu_tabt,method='euclidean')) #euclidean

# 文件夹名
Name <- paste0("Beta_diversity",ifelse(opts$method=="none","",paste0("_",opts$method)))
if (!file.exists(Name)){dir.create(Name)}

if (opts$tree != "none"){
  raw_tree <- readLines(opts$tree)
  writeLines(gsub(");$","):1;",raw_tree),paste0("Rooted_",opts$tree))
  rtree <- read.tree(paste0("Rooted_",opts$tree))
  
  ## 计算unifric距离矩阵
  unifracs <- GUniFrac(otu_tabt,rtree,alpha=c(0, 0.5, 1))$unifracs
  #
  dw <- unifracs[, , "d_1"]		# Weighted UniFrac
  du <- unifracs[, , "d_UW"]		# Unweighted UniFrac	
  #dv <- unifracs[, , "d_VAW"]		# Variance adjusted weighted UniFrac
  #alpha	：Parameter controlling weight on abundant lineages
  #d0 <- unifracs[, , "d_0"]     	# GUniFrac with alpha 0  
  #d5 <- unifracs[, , "d_0.5"]   	# GUniFrac with alpha 0.5 
  write.table(cbind(rownames(dw),dw),paste0(Name,"/weighted_unifrac_dm.txt"),sep = "\t",row.names = F, quote =FALSE)
  write.table(cbind(rownames(du),du),paste0(Name,"/unweighted_unifrac_dm.txt"),sep = "\t",row.names = F, quote =FALSE)
}
write.table(cbind(rownames(bc),bc),paste0(Name,"/bray_curtis_dm.txt"),sep = "\t",row.names = F, quote =FALSE)
write.table(cbind(rownames(jc),jc),paste0(Name,"/jaccard-binary_dm.txt"),sep = "\t",row.names = F, quote =FALSE)
write.table(cbind(rownames(ec),ec),paste0(Name,"/euclidean_dm.txt"),sep = "\t",row.names = F, quote =FALSE)
#
print(paste0("There will be output *dm in ",Name))
