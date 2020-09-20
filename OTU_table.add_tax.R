library(tidyverse)
library(optparse)

if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="rarefac.otu_taxa_table.xls",
                help="otu表格带有分类学信息"),
    make_option(c("-m", "--map"), type="character", default="map-group.txt",
                help="分组文件")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
}

data1 <- read_tsv(opts$map) %>% 
  rename("SampleID" = colnames(.)[1]) %>%
  mutate(group = fct_inorder(group)) %>% arrange(group)

data2_1 <- read_tsv(opts$input) %>%
  rename("OTUID" = colnames(.)[1]) %>%
  .[,c("OTUID",data1$SampleID,"taxonomy")] %>%
  mutate(SUM = sapply(1:nrow(.),function(x) sum(.[x,2:(ncol(.)-1)]))) %>% 
  filter(SUM > 0) %>% select(-SUM)
# 均一化
data2_2 <- data2_1
data2_2[,2:ncol(data2_2)] <- sapply(2:(ncol(data2_2)-1),function(x) data2_2[,x]/sum(data2_2[,x]))

# adict={'d__':'kingdom','p__':'phylum','c__':'class','o__':'order','f__':'family','g__':'genus','s__':'species'}
liaa <- c('kingdom','phylum','class','order','family','genus','species')

Tax <- data2_1$taxonomy %>% str_split(.,";")
for (i in seq_along(Tax)){# i <- 1
  if (length(Tax[[i]]) < 7){
    Tax[[i]] <- c(Tax[[i]],rep("",7-length(Tax[[i]])))
  }
}

Tax <- Tax %>% bind_cols(.) %>% t %>% as_tibble %>% rename_all(vars(liaa))
Out <- list()
for (j in 1:nrow(Tax)){
  if (any(as.logical(Tax[j,] != "") & !str_detect(Tax[j,],"__uncultured"))){
    Out[[j]] <- tail(unlist(Tax[j,as.logical(Tax[j,] != "") & 
                                  !str_detect(Tax[j,],"__uncultured") ]),1)
  }else{
    Out[[j]] <- ""
  }
}

data3_1 <- data2_1 %>% mutate(taxonomy = as.vector(unlist(Out))) %>% 
  mutate(OTUID =paste0(OTUID," ",taxonomy)) %>% select(-taxonomy) %>%
  write_tsv("rarefac.otu_taxa.xls")

data3_2 <- data2_2 %>% mutate(taxonomy = as.vector(unlist(Out))) %>% 
  mutate(OTUID =paste0(OTUID," ",taxonomy)) %>% select(-taxonomy) %>%
  write_tsv("rarefac.otu_taxa.percent.xls")
