# 程序使用示例
#

# 清理工作环境 clean enviroment object
rm(list=ls())

##检查依赖关系包的安装
package_list <- c("optparse","tidyverse","PerformanceAnalytics","corrplot","dplyr","digest","Hmisc",
                  "psych","gplots","RColorBrewer","igraph","ggplot2","reshape2","ggraph","tidygraph")

for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos="http://cran.r-project.org")
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
############################################# 载入包 #######################################
if (TRUE){
  option_list <- list(
    make_option(c("-p", "--per"), type="double", default=0.05,
                help="丰度筛选"),
    make_option(c("-r", "--rho"), type="double", default=0,
                help="rho threshhold,0.3;0.5;"),
    make_option(c("-v", "--value"), type="double", default=0.2,
                help="核心微生物筛选"),
    make_option(c("-j", "--pvalue"), type="double", default=0.05,
                help="pvalue"),
    make_option(c("-k", "--qvalue"), type="double", default=0.05,
                help="pvalue"),
    make_option(c("-u", "--feature"), type="character", default="feature_importance_scores-all.txt",
                help="feature_importance_scores.txt"),
    make_option(c("-y", "--Meande"), type="double", default=0.0001,
                help="feature_importance_scores 阈值")
    )
  opts <- parse_args(OptionParser(option_list=option_list))
}

##################################### function #############################################
## Genus_threshold
Genus_threshold <- function(nc,otu_coverage,genus){
  genus_threshold <- enframe(sapply(1:length(nc),function(y) 
    sum(as.numeric(otu_coverage >= nc[y])))) %>% 
    mutate(remian = length(otu_coverage) - value,
           coverage = nc) %>%
    mutate(per = sapply(1:length(nc),function(o){
      sum(unlist(genus[otu_coverage >= .$coverage[o],2:(ncol(genus))])) / 
        sum(unlist(genus[,2:(ncol(genus))]))
    })) %>% select(-name)
  
  p1 <- ggplot(genus_threshold,aes(coverage,per,color = "blue")) + 
    geom_vline(xintercept = opts$value, linetype="dotted") + 
    geom_point() + geom_line() + guides(color = FALSE)
  p2 <- ggplot(genus_threshold,aes(coverage,value,color = "red")) + 
    geom_vline(xintercept = opts$value, linetype="dotted") + 
    geom_point() + geom_line() + guides(color = FALSE)
  pp <- cowplot::plot_grid(p1, p2, labels = c('A', 'B'))
  ggsave("threshold.png",pp,width = 12,height = 8)
}

add0 <- function(a){# a <- test$rowname
  Max <- max(nchar(a))
  a <- as.numeric(a)
  sapply(a,function(aa) case_when(nchar(aa) < Max ~ paste0(rep("0",Max-nchar("1")),as.character(aa)),
                                  TRUE ~ as.character(aa)))
}

# point size
LM <- function(oo,y=c(1,3)){#oo <- data2 %>% colSums
  oo <- log10(oo)
  x = c(min(oo),max(oo))
  y = y
  o <- lm(y ~ x)
  predict(o,data.frame(x = oo))
}
##################################### Lipid ######################################################
data1 <- read_tsv("fecal.txt") %>% 
  rename_all(~c("F1","F2","F3","F4","GROUPS","F5","F6")) %>%
  #filter(GROUPS == "train") %>% 
  select(-GROUPS) %>% 
  filter(!is.na(F5)) %>%
  filter(F4 %in% c("CP","CPR")) %>%
  group_by(F4) %>% nest
  
data2 <- read_tsv("shetai.txt") %>% 
  rename_all(~c("T1","T2","T3","T4","GROUPS","T5","T6")) %>%
  #filter(GROUPS == "train") %>% 
  select(-GROUPS) %>% 
  filter(!is.na(T5)) %>%
  filter(T4 %in% c("CP","CPR")) %>%
  group_by(T4) %>% nest

data12 <- map2(data1$data,data2$data,function(a,b){# a <- data1$data[[1]];b <- data2$data[[1]]
  inner_join(a,b,by=c("F2" = "T2"))
}) %>% map(.,function(x) # x <- data12[[1]]
  x %>% mutate(group = map_chr(F3,function(x) str_replace_all(x,"[0-9]*",""))) %>%
  as.data.frame() %>% rownames_to_column() %>% 
    mutate(rowname = add0(rowname)) %>%
  mutate(SampleID = paste0(group,"_",rowname)) %>% select(-rowname)) %>%
  bind_rows() %>%
  select(F5,F6,T6,SampleID,group) %>%
  rename_all(~c("ID","Fecal","Tongue","SampleID","group"))

data3 <- read_tsv("Lipid/CP_CPR-diff-heatmap.txt") %>%
  rename("ID" = colnames(.)[1]) %>% 
  mutate(ID = map_chr(ID,function(x){#x <- data1$SampleID[[1]]
    x %>% str_split("-") %>% .[[1]] %>% .[2] %>% str_split("_") %>% .[[1]] %>% .[1]
  }))
#%>% gather(ID,value,-SampleID) %>% spread(SampleID,value)

data4 <- read_tsv("Lipid/Volcano_CP_CPR.txt") %>% 
  filter(VIP == "VIP >=1" & log_t.test_p.value_BHcorrect >= 2) %>% 
  filter(logFC >= 1 | logFC <= -1) %>%
  arrange(desc(VIP_num),desc(abs(logFC)))

SS <- data4 %>% 
  inner_join(.,colnames(data3) %>% as.matrix %>% as.data.frame %>% 
               as_tibble %>% mutate(V1 = as.character(V1)) %>% rename(LipidIon = colnames(.)[1])) %>%
  .$LipidIon %>% .[1:10]

data34 <- data3 %>% select(c(ID,SS)) %>% inner_join(data12,.) %>% select(-ID) %>% as_tibble

data34 %>%  select(-Fecal,-Tongue,-group) %>% write_tsv("Out/env.txt")

data34 %>% select(Fecal,Tongue,SampleID,group) %>% write_tsv("Out/Rename_Fecal_Tongue_Lipid.xls")

##################################### Fecal ###################################################
fata1 <- read_tsv("CP_CPR-Fecal/map-group.txt") %>% 
  right_join(.,data34[,c(1,3)],by=c(`#SampleID`="Fecal")) %>%
  select(SampleID,group) %>% rename(`#SampleID` = colnames(.)[1]) %>%
  write_tsv("CP_CPR-Fecal/Rename-map-group.txt")

fata2 <- read_tsv("CP_CPR-Fecal/rarefac.otu_genus.xls") %>%
  select(c(`OTU ID`,data34$Fecal)) %>%
  gather(Fecal,value,-`OTU ID`) %>%
  mutate_if(is.character,~fct_inorder(.x)) %>% 
  left_join(data34[,c(1,3)],.) %>% select(-Fecal) %>%
  spread(SampleID,value) %>% write_tsv("CP_CPR-Fecal/Rename-rarefac.otu_genus.xls") %>%
  rename("OTUID" = "OTU ID")

############## fata2 ################
# Uniform
fata2[,2:ncol(fata2)] <- sapply(2:ncol(fata2),function(x) fata2[,x]/sum(fata2[,x]))

# 丰度筛选
if(opts$per != 0){
  fata2 <- fata2 %>% 
    mutate(SELECT = sapply(1:nrow(.),function(x){
      any(.[x,2:ncol(.)]>=opts$per)})) %>% 
    .[.$SELECT,] %>% as_tibble %>% select(-SELECT)
}
# core microbiom
otu_coverage <- apply(fata2[,2:ncol(fata2)],1,function(x) 
  length(x[x>0])/(ncol(fata2)-1))
nc <- sort(unique(otu_coverage))
Genus_threshold(nc,otu_coverage,fata2)
fata2 <- fata2[otu_coverage >= opts$value,]
############## fata3-4 ######################
fata4 <- read_tsv("CP_CPR-Fecal/wilcox.otu.CP-CPR.xls") %>%
  filter(`p-value` <= 0.05)  %>% 
  rename(OTUID = colnames(.)[1])

if (opts$qvalue != 0){
  fata4 <- fata4 %>% filter(`q-value` <= opts$qvalue)
}
fata2 <- fata2 %>% inner_join(fata4[,1])

# randomforest
if (opts$feature != "none"){
  Meande <- read_tsv(str_c("CP_CPR-Fecal/",opts$feature)) %>% 
    filter(Mean_decrease_in_accuracy >= opts$Meande) %>% .$Feature_id
  
  fata2_md <- sapply(fata2$OTUID,function(x) strsplit(x," ", fixed=TRUE)[[1]][1])
  fata2 <- fata2[fata2_md %in% Meande,]
}

fata3 <- read_tsv("CP_CPR-Fecal/split_rarefac.otu_taxa_table.xls") %>%
  rename(OTUID = "OTU_ID") %>%
  select(OTUID,taxonomy,kingdom,phylum,class,order,family,genus) %>%
  left_join(fata2 %>% mutate(OTUID = sapply(OTUID,function(x) strsplit(x," ", fixed=TRUE)[[1]][1])),.) %>%
  write_tsv("CP_CPR-Fecal/Rename-split_rarefac.otu_taxa_table.xls")

fata2 <- left_join(fata2 %>% rename(OTUo = OTUID) %>% mutate(OTUID = map_chr(OTUo,function(x) 
  str_split(x," ") %>% .[[1]] %>% .[[1]])),
          fata3 %>% mutate(OTU_genus = paste0(OTUID," (",genus,")")) %>% select(OTUID,OTU_genus)) %>%
  mutate(OTUo = OTU_genus) %>% select(-OTU_genus,-OTUID) %>%
  rename(OTUID = OTUo) %>% 
  gather(tax,value,-OTUID) %>% 
  mutate(OTUID = fct_inorder(OTUID),
         tax = fct_inorder(tax)) %>% 
  spread(OTUID, value) %>%
  rename(SampleID = colnames(.)[1])
##################################### Tongue ##########################################################

tata1 <- read_tsv("CP_CPR-Tongue/map-group.txt") %>% 
  right_join(.,data34[,c(2,3)],by=c(`#SampleID`="Tongue")) %>%
  select(SampleID,group) %>% rename(`#SampleID` = colnames(.)[1]) %>%
  write_tsv("CP_CPR-Tongue/Rename-map-group.txt")

tata2 <- read_tsv("CP_CPR-Tongue/rarefac.otu_genus.xls") %>%
  select(c(`OTU ID`,data34$Tongue)) %>%
  gather(Tongue,value,-`OTU ID`) %>%
  mutate_if(is.character,~fct_inorder(.x)) %>% 
  left_join(data34[,c(2,3)],.) %>% select(-Tongue) %>%
  spread(SampleID,value) %>% write_tsv("CP_CPR-Tongue/Rename-rarefac.otu_genus.xls") %>%
  rename("OTUID" = "OTU ID")

############## tata2 ################
# Uniform
tata2[,2:ncol(tata2)] <- sapply(2:ncol(tata2),function(x) tata2[,x]/sum(tata2[,x]))

# 丰度筛选
if(opts$per != 0){
  tata2 <- tata2 %>% 
    mutate(SELECT = sapply(1:nrow(.),function(x){
      any(.[x,2:ncol(.)]>=opts$per)})) %>% 
    .[.$SELECT,] %>% as_tibble %>% select(-SELECT)
}
# core microbiom
otu_coverage <- apply(tata2[,2:ncol(tata2)],1,function(x) 
  length(x[x>0])/(ncol(tata2)-1))
nc <- sort(unique(otu_coverage))
Genus_threshold(nc,otu_coverage,tata2)
tata2 <- tata2[otu_coverage >= opts$value,]

# randomforest
if (opts$feature != "none"){
  Meande <- read_tsv(str_c("CP_CPR-Tongue/",opts$feature)) %>% 
    filter(Mean_decrease_in_accuracy >= opts$Meande) %>% .$Feature_id
  
  tata2_md <- sapply(tata2$OTUID,function(x) strsplit(x," ", fixed=TRUE)[[1]][1])
  tata2 <- tata2[tata2_md %in% Meande,]
}
############## tata3-4 ######################
tata4 <- read_tsv("CP_CPR-Tongue/wilcox.otu.CP-CPR.xls") %>%
  filter(`p-value` <= 0.05)  %>% 
  rename(OTUID = colnames(.)[1])

if (opts$qvalue != 0){
  tata4 <- tata4 %>% filter(`q-value` <= opts$qvalue)
}
tata2 %>% inner_join(tata4[,1])

# randomforest
if (opts$feature != "none"){
  Meande <- read_tsv(str_c("CP_CPR-Tongue/",opts$feature)) %>% 
    filter(Mean_decrease_in_accuracy >= opts$Meande) %>% .$Feature_id
  
  tata2_md <- sapply(tata2$OTUID,function(x) strsplit(x," ", fixed=TRUE)[[1]][1])
  tata2 <- tata2[tata2_md %in% Meande,]
}

tata3 <- read_tsv("CP_CPR-Tongue/split_rarefac.otu_taxa_table.xls") %>%
  rename(OTUID = "OTU_ID") %>%
  select(OTUID,taxonomy,kingdom,phylum,class,order,family,genus) %>%
  left_join(tata2 %>% mutate(OTUID = sapply(OTUID,function(x) strsplit(x," ", fixed=TRUE)[[1]][1])),.) %>%
  write_tsv("CP_CPR-Tongue/Rename-split_rarefac.otu_taxa_table.xls",na = "")


tata2 <- left_join(tata2 %>% rename(OTUo = OTUID) %>% mutate(OTUID = map_chr(OTUo,function(x) 
  str_split(x," ") %>% .[[1]] %>% .[[1]])),
  tata3 %>% mutate(OTU_genus = paste0(OTUID," (",genus,")")) %>% select(OTUID,OTU_genus)) %>%
  mutate(OTUo = OTU_genus) %>% select(-OTU_genus,-OTUID) %>%
  rename(OTUID = OTUo) %>% 
  gather(tax,value,-OTUID) %>% 
  mutate(OTUID = fct_inorder(OTUID),
         tax = fct_inorder(tax)) %>% 
  spread(OTUID, value) %>%
  rename(SampleID = colnames(.)[1])
############################# 合并 ####################

fout1 <- fata2 %>% rename_at(2:ncol(.),function(x) str_c("F_",x))
fout2 <- fout1 %>% gather(OTUID,value,-SampleID) %>% mutate(OTUID = fct_inorder(OTUID)) %>% 
  group_by(OTUID) %>% summarise(mean(value)) %>% rename(value = colnames(.)[2]) %>%
  separate(OTUID,sep="[^[0-9a-zA-Z]]",c("Class","OTUID")) %>% mutate(value = value*100)
fout3 <- fata3 %>% select_at(vars(-contains("_"))) %>% 
  select(-taxonomy) %>% add_column(Class = "F",.before = 1)

tout1 <- tata2 %>% rename_at(2:ncol(.),function(x) str_c("T_",x))
tout2 <- tout1 %>% gather(OTUID,value,-SampleID) %>% mutate(OTUID = fct_inorder(OTUID)) %>% 
  group_by(OTUID) %>% summarise(mean(value)) %>% rename(value = colnames(.)[2]) %>%
  separate(OTUID,sep="[^[0-9a-zA-Z]]",c("Class","OTUID")) %>% mutate(value = value*100)
tout3 <- tata3 %>% select_at(vars(-contains("_"))) %>% select(-taxonomy) %>% 
  add_column(Class = "T",.before = 1)

lout <- data34 %>% select(-Fecal,-Tongue,-group) %>% 
  gather(OTUID,value,-SampleID) %>% 
  mutate(OTUID = fct_inorder(OTUID)) %>% 
  group_by(OTUID) %>% summarise(mean(value)) %>% rename(value = colnames(.)[2]) %>%
  add_column(Class = "L",.before = 1) %>% mutate(phylum = "zLipid",value_scale = 8)

OUT0 <- data34 %>% select(-Fecal,-Tongue,-group) %>% cbind(.,fout1[,-1],tout1[,-1])
rownames(OUT0) <- OUT0$SampleID
OUT <- OUT0[,-1] %>% as.matrix
################################## spearman ############################################
Cor <- Hmisc::rcorr(OUT,type = "spearman")
write_tsv(data.frame("rho"=rownames(Cor$r),Cor$r),"rho.xls")
write_tsv(data.frame("p-value"=rownames(Cor$P),Cor$P),"p-value.xls")

rho <- Cor$r
pvalue <- Cor$P

ig1 <- rho %>% as.data.frame %>% rownames_to_column %>% melt %>% as_tibble %>% rename_all(~c("from","to","rho"))
ig2 <- pvalue %>% as.data.frame %>% rownames_to_column %>% melt %>% as_tibble %>% rename_all(~c("from","to","p"))
ig_v <- left_join(ig1,ig2) %>% 
  filter(!(is.na(p) & rho == 1)) %>% filter(p != 0) %>%
  mutate(rho_pn = map_chr(rho,~ case_when(.x > 0 ~ "positive",TRUE ~ "negative"))) %>%
  mutate(rho_abs = abs(rho)) %>%
  mutate(rank = map_chr(rho_abs,~ case_when(.x >= 0.8 ~ "0.8",.x >= 0.5 ~ "0.5",.x >= 0.3 ~ "0.3",TRUE ~ "0"))) %>%
  mutate(`p_0.05` = map_dbl(p,~ case_when(.x <= 0.05 ~ 1,TRUE ~ 0))) %>%
  mutate(`p_0.01` = map_dbl(p,~ case_when(.x <= 0.01 ~ 2,.x <= 0.05 ~ 1,TRUE ~ 0))) %>% 
  mutate(p_lg = map_dbl(p,~ log(.x,10) * (-1))) %>%
  mutate(to = as.character(to)) %>%
  mutate(CB = sapply(1:nrow(.),function(x) str_c(sort(as.character(.[x,1:2])),collapse = "-"))) %>%
  distinct(CB,.keep_all = TRUE) %>% select(-CB)

ig_n <- rbind(left_join(fout2,fout3),left_join(tout2,tout3)) %>%
  select(Class,OTUID,value,phylum) %>% 
  left_join(.,c(colnames(fout1)[-1],colnames(tout1)[-1]) %>% as_tibble %>% rename(ID = value) %>%
              mutate(OTU_genus = map_chr(ID,~ str_split(.x,"_") %>% .[[1]] %>% .[2])) %>% 
              mutate(OTUID = map_chr(OTU_genus,~ str_split(.x," ") %>% .[[1]] %>% .[1]))) %>%
  mutate(value_scale = LM(value)) %>% 
  select(Class,ID,value,phylum,value_scale) %>% rename("OTUID" = "ID") %>%
  bind_rows(.,lout) %>% mutate_if(is.character,~ fct_inorder(.x)) %>% 
  mutate(phylum = factor(phylum,levels=unique(phylum)[unique(phylum) != ""])) %>%
  select(OTUID,Class,phylum,value,value_scale)

##############################################################################
mycol <- c("#FF8E8A","#00CCCC","#F9F25C","#13678A","#00B631","#E55FAD","#E6271B","#999999","#000000")
col=c("#FF8E8A","#999999")
mypch <-c(21,22,23,24,25,3,4,16,15,17,18,7,8,9,10,11,12,13,14)
#
graph <- graph_from_data_frame(d=ig_v %>% filter(p <= 0.05,rho_abs>=0.3),vertices = ig_n)
## 
# Not specifying the layout - defaults to "auto"
extrafont::loadfonts()
p <- ggraph(graph, layout = 'circle') + 
  geom_edge_fan(aes(colour = rho_pn,edge_width = rho_abs,edge_alpha = p_lg),
                check_overlap= TRUE) + 
  geom_node_point(aes(fill = phylum, size = value, shape = Class)) + 
  #geom_node_point(aes(fill = phylum, size = value_scale, shape = Class)) + 
  geom_node_text(aes(label=name),size=3,check_overlap=TRUE,repel=TRUE) +
  scale_edge_width(range = c(0.1,2)) +
  scale_fill_manual(values=c(mycol[1:(nlevels(ig_n$phylum)-1)],"#FFFFFF")) +
  #scale_colour_manual(values=c(mycol[1:(nlevels(ig_n$phylum)-1)],"#FFFFFF")) +
  scale_shape_manual(values=mypch) +
  scale_size_continuous(range = c(1,8))+
  theme_graph() +
  guides(shape = guide_legend(override.aes = list(size = 3)),
         #colour = guide_legend(override.aes = list(size = 3)),
         fill = guide_legend(override.aes = list(size = 3)))

# save the modified plot to an object
ggsave("ggraph2.pdf",p,width = 11.55,height = 10)

################################ igraph ######################################
# 颜色设置;可以按照富集的分组或物种等级
mycol <- c("#FF8E8A","#00CCCC","#F9F25C","#13678A","#00B631","#E55FAD","#E6271B","#999999","#000000")
col=c("#FF8E8A","#999999")

net <- graph_from_data_frame(d=ig_v, vertices=ig_n, directed=T) # colnames(ig_n)
#choose edges levels to show
net <- delete_edges(net, E(net)[rank<2])

# 颜色设置;可以按照富集的分组或物种等级
ig_col <- c("#FF8E8A","#00CCCC","#F9F25C","#13678A","#00B631","#E55FAD","#E6271B","#999999","#000000")
col=c("#FF8E8A","#999999")
V(net)$color <- ig_col[factor(V(net)$phylum,levels=unique(V(net)$phylum))]

# 节点大小设置
V(net)$size <- V(net)$value_scale

#设置线宽
E(net)$width <- E(net)$rank
#设置连接线的颜色
edge.col <- col[E(net)$rho_pn]

#出多张图
#layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 
#layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]
#par(mfrow=c(3,3), mar=c(1,1,1,1))
#for (layout in layouts) {
#  print(layout)
#  l <- do.call(layout, list(net)) 
#plot(net,edge.arrow.size=0,edge.color=edge.col,layout=l)
#}

#draw the picture
par(mfrow=c(1,3), mar=c(1,1,1,1))
pdf("igraph.pdf")
plot(net,edge.arrow.size=0,edge.color=edge.col)
plot(net,edge.arrow.size=0,edge.color=edge.col,edge.lty=E(net)$p_005)
plot(net,edge.arrow.size=0,edge.color=edge.col,edge.lty=E(net)$p_001)
##加legend
legend(x=-1, y=-1, unique(V(net)$phylum), pch=21,
       col="#777777", pt.bg=mycol, pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()




















