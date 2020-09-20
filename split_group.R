# Rscript split_group.R -g map-group.txt -m none -n 0
# 排列组合，得到所有可能的分组结果
library(optparse)
library(magrittr)
library(tidyverse)

if (TRUE){
  option_list <- list(
    make_option(c("-g", "--map"), type="character", default="map-group.txt",
                help="分组文件"),
    make_option(c("-m", "--mode"), type="character", default="none",
                help="列出所有组合:none；venn要求小于等于5组"),
    make_option(c("-n", "--num"), type="numeric", default=0,
                help="0为自动检测，划分分组；或指定")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
  print(paste("The input file is ", opts$map,  sep = ""))
  print(paste("The mode is ", opts$mode,  sep = ""))
  print(paste("The number is ", opts$num,  sep = ""))
}
# 读取分组文件
mp1 <- read.table(opts$map,head= T ,sep="\t",comment.char = "",fileEncoding = "UTF-8")
rownames(mp1) <- mp1[,1]

# 判断分组个数
if (opts$num == 0 & opts$mode == "venn"){
  mpn <- ifelse(nlevels(mp1$group) >= 5,5,nlevels(mp1$group))
}else if (opts$num == 0 & opts$mode != "venn"){
  mpn <- nlevels(mp1$group)
}else if (opts$num != 0 & opts$mode == "venn"){
  mpn <- ifelse(nlevels(mp1$group) >= 5,5,nlevels(mp1$group))
}else if (opts$num != 0 & opts$mode != "venn"){
  mpn <- opts$num
}

gp1 <- as.vector(unique(mp1$group))
for (i in 2:mpn){#i=2
  gp2 <- combn(gp1,i)
  for (j in 1:ncol(gp2)){#j=1
    mp1 %<>% mutate(pp = unlist(lapply(1:nrow(mp1),function(k){# k = 1
      ifelse(as.vector(mp1[k,2]) %in% gp2[,j],as.vector(mp1[k,2]),"")
    })))
    colnames(mp1)[ncol(mp1)] <- paste(gp2[,j],collapse="-")
    
    mp3 <- mp1[mp1$group %in% gp2[,j],c(1,2)]
    nm <- paste0("map",length(gp2[,j]),".",paste(gp2[,j],collapse="-"),".txt")
    write("#sample\tgroup",file=nm)
    write.table(file=nm,mp3,sep="\t",quote = FALSE,row.names = FALSE, col.names = FALSE,append = TRUE)
  }
}
colnames(mp1)[1] <- "#SampleID"
write.table(file="mapping.txt",mp1,sep="\t",quote = FALSE,row.names = FALSE)
print("OK!")
