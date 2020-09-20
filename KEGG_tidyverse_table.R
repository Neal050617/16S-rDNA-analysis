library(tidyverse)
library(magrittr)
library(optparse)

if (TRUE){
  option_list <- list(
    make_option(c("-o", "--KO"), type="character", default="metagenome_predictions.tsv",
                help="ko表格"),
    make_option(c("-a", "--ko_anno"), type="character", default="ko00001-anno.xls",
                help="ko00001-anno.xls"),
    make_option(c("-g", "--kegg"), type="character", default="predicted_metagenome_l3.tsv",
                help="kegg表格"),
    make_option(c("-m", "--map"), type="character", default="map-group.txt",
                help="分组文件")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
}

###############################################  function  ########################################################################
.community <- function(dd,nm=''){# dd <- kegg
  lineage=c("KEGG_Pathways_Level_1","KEGG_Pathways_Level_2","KEGG_Pathways_Level_3")
  sapply(1:length(lineage),function(l){ # l <- 1
    Name <- colnames(dd)[ncol(dd)-3+l]
    ddd <- dd[,c(ncol(dd)-3+l,2:(ncol(dd)-4))] %>% rename("KEGG" = Name) %>% 
      group_by(KEGG) %>% 
      nest %>% mutate(data = map(data,colSums)) %>% 
      mutate(data = map(data,as.matrix)) %>%
      mutate(data = map(data,t)) %>% 
      mutate(data = map(data,as_tibble)) %>% 
      unnest(data)
    write_tsv(ddd,str_c(lineage[l],nm,".tsv"))
  })
  return(nm)
}

.lefse_data <- function(ll,nm='',map){#ll <- kegg
  lineage=c("KEGG_Pathways_Level_1","KEGG_Pathways_Level_2","KEGG_Pathways_Level_3")
  ddd <- lapply(1:length(lineage),function(l){ # l <- 1
    ll %>% mutate(KEGG_Pathways = sapply(1:nrow(.),function(s){#s <- 1
      ifelse(l>1,
             str_c(unlist(.[s,(ncol(ll)-2):(ncol(ll)-3+l)]),collapse = "|"),
             as.vector(unlist(.[s,ncol(ll)-2])))})) %>% 
      .[,c(2:(ncol(ll)-3))] %>% 
      group_by(KEGG_Pathways) %>% 
      nest %>% mutate(data = map(data,colSums)) %>% 
      mutate(data = map(data,as.matrix)) %>%
      mutate(data = map(data,t)) %>% 
      mutate(data = map(data,as_tibble)) %>% 
      unnest(data)
  }) %>% bind_rows
  colnames(ddd) <- c("class",as.vector(map$group))
    
  write_tsv(ddd,str_c("KEGG.lefse.data",nm,".tsv"))
}

map <- read_tsv(opts$map) %>% 
  mutate(`#SampleID` = fct_inorder(`#SampleID`),group = fct_inorder(group))

ka <- read_tsv(opts$ko_anno,col_names = FALSE) %>% rename(`#OTU ID`=X1,EC=X2)

ko <- read_tsv(opts$KO,skip = 1) %>% 
  rename_all(funs(str_replace_all(., "-", "_"))) %>% 
  .[,c("#OTU ID",as.vector(map$`#SampleID`),"KEGG_Pathways")] %>%
  mutate(SUM = sapply(1:nrow(.),function(s){
    sum(.[s,2:(ncol(.)-1)])
  })) %>% filter(SUM > 0) %>% select(-SUM) %>% 
  left_join(ka) %>% filter(!is.na(EC)) %>%
  .[,c(ncol(.),2:(ncol(.)-2))] %>% 
  group_by(EC) %>% nest %>% 
  mutate(data = map(data,colSums)) %>% 
  mutate(data = map(data,as.matrix)) %>%
  mutate(data = map(data,t)) %>% 
  mutate(data = map(data,as_tibble)) %>% 
  unnest(data)

ko_norm <- ko %>% .[,2:ncol(ko)] %>% ungroup() %>%
  lapply(.,function(n){n/sum(n)}) %>% 
  cbind(ko[,1],.) %>%
  as_tibble

ko_lefse <- ko 
colnames(ko_lefse) <- c("class",as.vector(map$group))

write_tsv(ko,"KEGG_Pathways_Level_ko.tsv")
write_tsv(ko_norm,"KEGG_Pathways_Level_ko.percentage.tsv")
write_tsv(ko_lefse,"ko.lefse.data.tsv")

#######################################################################
kegg <- read_tsv(opts$kegg,skip = 1) %>% 
  rename_all(funs(str_replace_all(., "-", "_"))) %>%
  mutate(`#OTU ID` = sapply(seq_along(.$`#OTU ID`),function(x) 
    str_replace_all(.$`#OTU ID`[x]," ","_"))) %>%
  mutate(KEGG_Pathways = sapply(seq_along(.$KEGG_Pathways),function(y) 
    str_replace_all(.$KEGG_Pathways[y],"; ",";") %>% 
      str_replace_all(.," ","_")
    )) %>% .[,c("#OTU ID",as.vector(map$`#SampleID`),"KEGG_Pathways")] %>%
  mutate(kegg_pathway = KEGG_Pathways) %>%
  separate(kegg_pathway,sep = ";",into = c("L1","L2","L3")) %>% 
  arrange(L1,L2,L3) %>% 
  mutate(SUM = sapply(1:nrow(.),function(z) sum(.[z,2:(ncol(.)-4)]))) %>%
  filter(SUM > 0) %>% select(-SUM)

kegg_rank <- kegg %>% 
  mutate(L1 = sapply(1:nrow(.),function(r) str_c("L1_",.$L1[r]))) %>%
  mutate(L2 = sapply(1:nrow(.),function(r) str_c("L2_",.$L2[r]))) %>%
  mutate(L3 = sapply(1:nrow(.),function(r) str_c("L3_",.$L3[r]))) %>%
  mutate(KEGG_Pathways = paste0(L1,L2,L3,sep = ";"))

kegg_norm <- kegg[,2:(ncol(kegg)-4)] %>% 
  sapply(.,function(n){n/sum(n)}) %>% 
  cbind(kegg[,1],.,kegg[,(ncol(kegg)-3):ncol(kegg)]) %>%
  as_tibble

kegg_rank_norm <- kegg_rank[,2:(ncol(kegg_rank)-4)] %>% 
  sapply(.,function(n){n/sum(n)}) %>% 
  cbind(kegg_rank[,1],.,kegg_rank[,(ncol(kegg_rank)-3):ncol(kegg_rank)]) %>%
  as_tibble

.community(kegg)
.community(kegg_norm,".percentage")
.lefse_data(kegg_rank,nm='',map)

write_tsv(kegg,"KEGG_Pathways.report.xls")
write_tsv(kegg_norm,"KEGG_Pathways.report.percentage.xls")
write_tsv(kegg_rank,"KEGG_Pathways.rank.xls")
write_tsv(kegg_rank_norm,"KEGG_Pathways.rank.percentage.xls")



#mutate(KEGG_Pathways = 
#         (.[,ncol(.)] %>% unlist %>% 
#            sapply(.,function(k){
#              k %>% str_replace_all("\"","") %>% 
#                str_replace_all(", ",",") %>% 
#                str_replace_all(" ","_") %>% 
#                str_replace_all("\\[\\[|\\]\\]","") %>%
#                str_split("\\],\\[") %>%
#                .[[1]] %>% .[length(.)] %>% str_split(",") %>% 
#                .[[1]] %>% .[length(.)]
#            }))) %>% mutate(KO = paste(`#OTU ID`,KEGG_Pathways,sep="__"))#