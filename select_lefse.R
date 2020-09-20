#! /usr/lib64/R/bin/Rscript --vanilla

library(tidyverse)
library(optparse)
library(stringr)

if (TRUE){
  option_list <- list(
    make_option(c("-s", "--select"), type="character", default="g__", help="筛选"),
    make_option(c("-i", "--input"), type="character", default="data.res", help="输入")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
}

data1 <- read_tsv(opts$input,col_names = FALSE) %>% 
  #mutate_at(colnames(.)[2:ncol(.)],~as.character)
  #mutate_if(is.numeric,~replace_na(""))
  #mutate_all(~replace_na("")) %>%
  .[sapply(.[,1],function(x) str_detect(x,opts$select)),] %>%
  write_tsv(.,str_c(opts$select,"data.res"),col_names = FALSE,na = "")
