library(optparse)
library(tidyverse)
library(magrittr)
library(reshape2)
library(optparse)
library(formattable)
library(tidyverse)
library(htmltools)
library(webshot)

if (TRUE){
  option_list <- list(
    make_option(c("-a", "--aa"), type="character", default="seqkit.R1.out",
                help="原始序列统计"),
    make_option(c("-b", "--bb"), type="character", default="seqstat.xls",
                help="质控统计表"),
    make_option(c("-c", "--cc"), type="character", default="otu_table.xls",
                help="OTU表格"),
    make_option(c("-d", "--dd"), type="character", default="otu_taxa_table.xls",
                help="OTU注释"),
    make_option(c("-o", "--out"), type="character", default="16S_rRNA_sequencing.xls",
                help="指定图形的宽度")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
}

data1 <- read_tsv(opts$aa) %>% .[,c(1,4)] %>% 
  rename("SampleID"=file,"Raw_Tags"=num_seqs) %>%
  mutate(SampleID = str_replace_all(SampleID,c("_R1.fastq"="")))

data2 <- read_tsv(opts$bb) %>% .[,c(1,2)] %>% 
  rename("SampleID"=X1,"Clean_Tags"=Seq_num)

dd <- getwd()
d1 <- list.dirs() #查看当前目录的子目录
d2 <- d1[str_detect(d1,"analysis.*process.*otu_0.97$")]

data3 <- read_tsv(str_c(d2,opts$cc,sep = "/")) %>%
  rename("OTU_ID" = colnames(.)[1]) %>%
  gather(var, value, -OTU_ID) %>% 
  mutate(OTU_ID = factor(OTU_ID, 
                         levels=str_sort(unique(OTU_ID), numeric = TRUE))) %>%
  spread(OTU_ID, value) %>% 
  rename(SampleID = var) %>%
  mutate(Final_Tags = rowSums(.[,2:ncol(.)])) %>%
  .[,c(1,ncol(.))]

data4 <- read_tsv(str_c(d2,opts$dd,sep = "/")) %>%
  rename("OTU_ID" = colnames(.)[1]) %>%
  select(-taxonomy) %>%
  gather(var, value, -OTU_ID) %>% 
  mutate(OTU_ID = factor(OTU_ID, 
                         levels=str_sort(unique(OTU_ID), numeric = TRUE))) %>%
  spread(OTU_ID, value) %>% 
  rename(SampleID = var) %>%
  mutate_at(vars(2:ncol(.)),function(x){
    sapply(x,function(y)ifelse(y>0,1,0))
  }) %>%
  mutate(OTU = rowSums(.[,2:ncol(.)])) %>%
  .[,c(1,ncol(.))]


Data <- left_join(left_join(left_join(data4,data3),data2),data1) %>%
  .[,c("SampleID","Raw_Tags","Clean_Tags","Final_Tags","OTU")] %>%
  write_tsv(opts$out)

## library(gridExtra)
## pdf("mypdf.pdf", height=6, width=4)
## grid.table(read.table("16S_rRNA_sequencing.xls"))
## dev.off()

## webshot::install_phantomjs()
export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}

FT <- formattable(Data,align = rep("c", NCOL(data1) ))
export_formattable(FT,"16S_rRNA_sequencing.png")




