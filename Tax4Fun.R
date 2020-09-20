# require(devtools)
# install_version("qiimer", version = "0.9.4",repos = "http://cran.us.r-project.org")
# install_version("biom", version = "0.3.12",repos = "http://cran.us.r-project.org")
# install.packages("D:/Database/tax4fun_database/Tax4Fun_0.3.1.tar.gz", repos = NULL, type = "source")

library(Tax4Fun)
library(tidyverse)
library(stringr)
# 导入OTUs表
QIIMESingleData <- importQIIMEData("otu_table_tax.txt")

# 可选步骤
# 发现OTU表为按Taxonomy合并的结果
otu_table <- QIIMESingleData$otuTable

# 根据Tax4Fun提供的SILVA123最新数据库进行预测，要求此数据的压缩包拉于此工作目录，这个命令得出来的是KO号的各种酶的基因丰度
Tax4FunOutput1 <- Tax4Fun(QIIMESingleData, "/work/users/chaoliu/database/Tax4Fun/SILVA123_Tax4Fun", 
                         fctProfiling = TRUE, refProfile = "UProC",
                         shortReadMode = TRUE, normCopyNo = TRUE) 

# 提取KO表，生成6508个KO相关的通路
KO_table1 <- t(Tax4FunOutput1$Tax4FunProfile)
# 所有样品标准化为1, 没有原始数据，可以用anova和Limma，无法使用edgeR和DESeq2
# 输出KO表，表头写个制表符，用于对齐表头
write.table("ID\t", file="KO_table.txt",append = FALSE, quote = FALSE, sep="\t",eol = "", na = "NA", dec = ".", row.names = F,col.names = F)
write.table(KO_table1, file="KO_table.txt",append = T, quote = FALSE, sep="\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE)

KO_table1 %>% as.data.frame %>% rownames_to_column %>% as_tibble %>%
  mutate(rowname = sapply(1:nrow(.),function(x){# x <- test$rowname[993]
    .$rowname[x] %>% 
      ifelse(str_detect(.,","),str_replace_all(.,",",""),.) %>%
      ifelse(str_detect(.," \\/ "),str_replace_all(.,"\\/","_"),.) %>% 
      ifelse(str_detect(.,"\\+"),str_replace_all(.,"\\+","_"),.) %>%
      ifelse(str_detect(.,"\\/"),str_replace_all(.,"\\/","_"),.) %>%
      ifelse(str_detect(.,"\\'"),str_replace_all(.,"\\'","_"),.) %>%
      ifelse(str_detect(.,"-"),str_replace_all(.,"-","_"),.) %>%
      #str_remove_all("^K[0-9].*; \\[EC.*\\]") %>% 
      #str_remove_all("^K[0-9].*; ") %>% 
      str_remove_all(" \\[EC.*\\]") %>% 
      str_replace_all(.,";","") %>% 
      ifelse(str_detect(.," "),str_replace_all(.," ","_"),.) %>%
      str_replace_all(.,"\\(","_") %>% str_replace_all(.,"\\)","_") %>%
      str_replace_all(.,c("___" = "_")) %>% 
      str_replace_all(.,c("__" = "_")) %>%
      str_replace_all(.,"_$","")
  })) %>% write_tsv("Tax4Fun_Pathways_ko.tsv")

# /work/users/chaoliu/database/Tax4Fun/SILVA123_Tax4Fun
# 上步结果中KO名称，只关注代谢通路的变化，调整fctProfiling=FALSE
Tax4FunOutput2 <- Tax4Fun(QIIMESingleData, "/work/users/chaoliu/database/Tax4Fun/SILVA123_Tax4Fun", 
                          fctProfiling = FALSE, refProfile = "UProC",
                          shortReadMode = TRUE, normCopyNo = TRUE) 

# 提取KO表，生成275个代谢相关的通路
KO_table2 <- t(Tax4FunOutput2$Tax4FunProfile)
# 表头写个制表符，用于对齐表头
write.table("ID\t", file="KO_table_fct.txt",append = FALSE, quote = FALSE, sep="\t",eol = "", na = "NA", dec = ".", row.names = F,col.names = F)

write.table(KO_table2, file="KO_table_fct.txt",append = T, quote = FALSE, sep="\t",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE)

KO_table2 %>% as.data.frame %>% rownames_to_column %>% as_tibble %>%
  mutate(rowname = sapply(1:nrow(.),function(x){# x <- test$rowname[993]
    .$rowname[x] %>% 
      ifelse(str_detect(.,","),str_replace_all(.,",",""),.) %>%
      ifelse(str_detect(.," \\/ "),str_replace_all(.,"\\/","_"),.) %>% 
      ifelse(str_detect(.,"\\+"),str_replace_all(.,"\\+","_"),.) %>%
      ifelse(str_detect(.,"\\/"),str_replace_all(.,"\\/","_"),.) %>%
      ifelse(str_detect(.,"\\'"),str_replace_all(.,"\\'","_"),.) %>%
      ifelse(str_detect(.,"-"),str_replace_all(.,"-","_"),.) %>%
      str_remove_all("^ko[0-9].*; ") %>% 
      str_replace_all(.,";","") %>% 
      ifelse(str_detect(.," "),str_replace_all(.," ","_"),.) %>%
      str_replace_all(.,"\\(","_") %>% str_replace_all(.,"\\)","_") %>%
      str_replace_all(.,c("___" = "_")) %>% 
      str_replace_all(.,c("__" = "_")) %>%
      str_replace_all(.,"_$","")
  })) %>% write_tsv("Tax4Fun_Pathways_Metabolites.tsv")

