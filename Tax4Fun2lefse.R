library(tidyverse)
library(optparse)

if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="Tax4Fun_Pathways_ko.tsv",
                help="输入的OTU表格"),
    make_option(c("-g", "--map"), type="character", default="map-group.txt",
                help="分组文件")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
}

Gp <- read_tsv(opts$map) %>% rename("SampleID" = colnames(.)[1])

data1 <- read_tsv(opts$input) %>% rename("ID" = colnames(.)[1]) %>% 
  gather(SampleID,value,-ID) %>% right_join(.,Gp) %>% 
  rename("class" = "group") %>% spread(ID,value,convert = TRUE) %>% 
  t %>% as.data.frame %>% rownames_to_column %>% as_tibble %>% .[-1,]

data2 <- data1 %>% .[-1,] %>% mutate_if(is.factor,function(x) as.character(x) %>% as.numeric(.)) %>%
  mutate(SUM = sapply(1:nrow(.),function(x) sum(.[x,2:(ncol(.)-1)]))) %>% 
  filter(SUM > 0) %>% select(-SUM)

colnames(data2) <- as.character(data1[1,] %>% t)

NM <- str_split(opts$input,"_")[[1]][3] %>% str_split(.,"\\.") %>% .[[1]] %>% .[1]

write_tsv(data2,str_c(NM,".lefse.data.tsv"))


