library(tidyverse)
library(optparse)

if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="rarefac.otu_genus.xls",
                help="otu表格而已，带有分类学信息的表格最好别用"),
    make_option(c("-n", "--num"), type="numeric", default=0.67,
                help="拆分数据"),
    make_option(c("-m", "--map"), type="character", default="none",
                help="分组数据")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
}

data <- read_tsv(opts$input) %>% rename("SampleID"=colnames(.)[1]) %>% 
  rename_all(funs(str_replace(., "-", "_"))) %>% 
  mutate(SampleID = sapply(seq_along(.$SampleID),function(x) 
    str_replace_all(.$SampleID[x]," ","_") %>% 
      str_replace_all(.,"\\(","") %>%
      str_replace_all(.,"\\)","")))

if (opts$map != "none"){
  Map <- read_tsv(opts$map) %>% rename("SampleID"=colnames(.)[1])
  
  MAP <- Map %>% t %>% as_tibble %>% rownames_to_column %>%
    mutate(rowname = c("SampleID","class"))
  colnames(MAP) <- MAP[1,]
  MAP <- MAP %>% .[-1,]
  
  data <- data %>% .[,c("SampleID",Map$SampleID)] %>% 
    gather(var, value, -SampleID) %>% 
    mutate(SampleID = factor(SampleID, 
                             levels=str_sort(unique(SampleID), numeric = TRUE))) %>% 
    spread(SampleID, value) %>%
    rename(SampleID = var) %>%
    left_join(Map) %>%
    .[,c(TRUE,(.[,c(-1,-ncol(.))] %>% colSums) > 0,TRUE)] %>% 
    .[,-ncol(.)] %>% 
    gather(var, value, -SampleID) %>% 
    mutate(SampleID = factor(SampleID, 
                             levels=str_sort(unique(SampleID), numeric = TRUE))) %>% 
    spread(SampleID, value) %>% rename(SampleID = var) %>%
    rbind(MAP,.)
}else{
  MAP <- colnames(data) %>% matrix(.,nrow=1) %>% as_tibble
  colnames(MAP) <- MAP[1,]
  data <- data %>% rbind(MAP,.)
}

set.seed(20190618)
ind=sample(2,(ncol(data)-1),replace = TRUE,prob = c(opts$num,1-opts$num))

data1 <- data %>% .[,c(TRUE,ind==1)] %>% 
  write_tsv(str_c("Part1-",opts$input),col_names = FALSE)
data2 <- data %>% .[,c(TRUE,ind==2)] %>% 
  write_tsv(str_c("Part2-",opts$input),col_names = FALSE)
data3 <- colnames(data)[-1] %>% enframe %>% 
  cbind(ind) %>% .[,-1]

if (opts$map != "none"){
  data3 <- data3 %>% rename(SampleID = value) %>% 
    left_join(Map) %>% arrange(ind)
}

write_tsv(data3,"Part.list",col_names = FALSE)



