## Part1. ASV
```shell
mkdir preparation
/usr/local/R-3.6.0/bin/R -e '.libPaths("/usr/local/R-3.6.0/lib64/R/library");
library(readr);library(tidyverse);
otu <- read_tsv("merged/otu_table.xls") 
n.otu <- cbind(ID=data.frame(qq="ASV",ww=c(1:nrow(otu))) %>% mutate(ee=paste(qq,ww,sep="")) %>% .[,3],
               otu %>% rename_at(1,~"name"))
n.otu[is.na(n.otu)] <- 0
n.otu[,-2] %>% rename_at(1,~"ASV ID") %>% 
  write_tsv("preparation/otu_table.xls")

Tax <- read_tsv("merged/taxonomy.tsv")
Rename <- n.otu[,1:2] %>% write_tsv("ASV.rename")

Rename %>% left_join(Tax,by=c("name"="Feature ID")) %>% 
  select(-"name") %>% rename_at(1,~"ASV ID") %>% 
  write_tsv("preparation/rep-seqs-taxonomy.tsv")

data2 <- sapply(3:ncol(n.otu),function(x) sum(as.numeric(unlist(n.otu[,x]))))
data3 <- data.frame("ID" = colnames(n.otu)[-c(1:2)], "count" = data2) %>% mutate(rank = rank(count))
write_tsv(data3 %>% arrange(count), "preparation/p-sampling-sort.xls")
'

biom convert -i preparation/otu_table.xls -o 04.feature-table.biom --to-hdf5 --table-type="OTU table"
qiime tools import --input-path 04.feature-table.biom --type 'FeatureTable[Frequency]' --output-path 04.otu_table.qza --input-format BIOMV210Format

# 更改代表序列名称：
cd preparation
ln -s ../merged/dna-sequences.fasta ./
less dna-sequences.fasta | paste - - | sed '1i ASVID\tseq' > rep.fa
/usr/local/R-3.6.0/bin/R -e '.libPaths("/usr/local/R-3.6.0/lib64/R/library");pacman::p_load(tidyverse,stringr,magrittr);Rename <- read_tsv("../ASV.rename");rep <- read_tsv("rep.fa") %>% mutate(ASVID = map_chr(ASVID,~str_replace_all(.x,">",""))) %>% inner_join(Rename,.,by=c("name"="ASVID")) %>% select(-name) %>% rename_at(1,~"ASVID") %>% mutate(ASVID = map_chr(ASVID,~str_c(">",.x)));write_tsv(rep,"rep.xls")'
less rep.xls | sed '1d' | sed "s/\r//g" | tr "\t" "\n" > rep-seqs.fasta
cd ..

# 将代表序列转换成qza格式：
qiime tools import --type 'FeatureData[Sequence]' --input-path preparation/rep-seqs.fasta --output-path 04.rep-seqs.qza

# 代表序列统计
qiime feature-table tabulate-seqs --i-data 04.rep-seqs.qza --o-visualization 04.rep-seqs.qzv

# 系统发育树能够服务于后续多样性分析
qiime phylogeny align-to-tree-mafft-fasttree \
      --i-sequences 04.rep-seqs.qza \
      --o-alignment 04.rep-seqs_aligned.qza \
      --o-masked-alignment 04.rep_seqs_masked.qza \
      --p-n-threads 20 \
      --o-tree preparation/unrooted-tree.qza \
      --o-rooted-tree preparation/rooted-tree.qza

qiime tools export --input-path preparation/rooted-tree.qza --output-path preparation # rooted-tree.nwk 
mv preparation/tree.nwk preparation/rooted-tree.nwk 
qiime tools export --input-path preparation/unrooted-tree.qza --output-path preparation # unrooted-tree.nwk 
mv preparation/tree.nwk preparation/unrooted-tree.nwk 

#########################################################################################
NOW=`env LANG=en_US.UTF-8 date +%a_%b_%d_%Y-%m-%d_%H_%M_%S`
PWD=`pwd`
PWD2=$PWD"/analysis_"$NOW
mkdir ${PWD2}
mkdir ${PWD2}/results
mkdir ${PWD2}/results/OTU_Taxa 
mkdir ${PWD2}/results/Estimators
mkdir ${PWD2}/results/Rarefactions
mkdir ${PWD2}/results/Community
mkdir ${PWD2}/results/Shannon_rarefac
mkdir ${PWD2}/results/Rank_abundance
mkdir ${PWD2}/results/Specaccum
mkdir ${PWD2}/results/Heatmap
mkdir ${PWD2}/results/Beta_diversity
mkdir ${PWD2}/results/Pca
mkdir ${PWD2}/results/Pcoa
mkdir ${PWD2}/results/Nmds
mkdir ${PWD2}/results/Hcluster_tree
mkdir ${PWD2}/results/Hclust_bar
mkdir ${PWD2}/results/Picrust2
mkdir ${PWD2}/results/Tax4Fun2

mkdir ${PWD2}/process ${PWD2}/process/ASV

cd ${PWD2}/process/ASV
\cp $PWD2/../preparation/rep-seqs.fasta ./otu_reps.fasta
\cp $PWD2/../preparation/otu_table.xls ./
\cp $PWD2/../preparation/rep-seqs-taxonomy.tsv ./otu_reps_tax_assignments.txt
\cp $PWD2/../preparation/rooted-tree.nwk ./otu_reps.raw_aligned.fasta.tre
\cp $PWD2/../preparation/unrooted-tree.nwk ./

## fix taxonomy
sed -i 's/ //g' otu_reps_tax_assignments.txt
sed -i 's/"//g' otu_table.xls
awk 'NR==FNR{ a[$1] }NR>FNR{ if($1 in a) print $0}' otu_table.xls otu_reps_tax_assignments.txt > otu_select_fix_cluster_tax_assignments.txt

/usr/local/R-3.6.0/bin/R -e '.libPaths("/usr/local/R-3.6.0/lib64/R/library"); library(tidyverse); taxa <- "otu_select_fix_cluster_tax_assignments.txt" %>% read_tsv(col_names=F) %>% rename_all(~c("ASV ID", "taxonomy" ,"score")) ; table <- "otu_table.xls" %>% read_tsv() %>% left_join(taxa) %>% mutate(taxonomy = map_chr(taxonomy,~replace_na(.x,""))) %>% select(-score) %>% write_tsv("otu_taxa_table.xls") '

source /root/anaconda3/etc/profile.d/conda.sh
conda activate qiime1

biom convert -i otu_table.xls  -o otu_table.biom --table-type "OTU table" --to-hdf5
less otu_table.xls | grep --color 'ASV ID' | sed 's/ASV ID//' | xargs -n1 | sort -n | awk '{print $1"\t"$1}' > sample_order

#subsample
less otu_table.xls | awk 'BEGIN{OFS="\t"}{for(i=2;i<=NF;i++)$i=(a[i]+=$i)}END{$1="";print}' | python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" | sort -n | awk 'NR==1{print}' > subsample.num

#判断
/usr/local/R-3.6.0/bin/R -e '.libPaths("/usr/local/R-3.6.0/lib64/R/library");library(data.table);library(tidyverse)
subsample.num <- read_tsv("subsample.num",col_names = F) %>% as.numeric
if(subsample.num > 10000){
  print("the lost is more 10000")
}else{
  otu_table <- fread("otu_table.xls")
  #超过10000的要抽
  otu_table1 <- otu_table %>% select_if(function(x) !is.numeric(x)|(is.numeric(x)&&sum(x)>10000)) %>% fwrite("otu_table1.xls",sep = "\t")
  #低于10000的保留
  otu_table2 <- otu_table %>% select_if(function(x) !is.numeric(x)|(is.numeric(x)&&sum(x)<=10000)) %>% fwrite("otu_table2.xls",sep = "\t")
  print("嗷嗷嗷 最少的比10000小耶")
}
'

# 将otu_table1.xls的样本单独进行处理
if [ -f otu_table1.xls ];then
# 按照超过10000部分的最小值去抽
less otu_table1.xls | awk 'BEGIN{OFS="\t"}{for(i=2;i<=NF;i++)$i=(a[i]+=$i)}END{$1="";print}' | python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" | sort -n | awk 'NR==1{print}' > subsample2.num
biom convert -i otu_table1.xls -o otu_table1.biom --table-type "OTU table"  --to-hdf5
single_rarefaction.py -i otu_table1.biom -o rarefac.otu_table1.biom -d `less subsample2.num` 
# 将第一个表转化并将两表合并
biom convert -i rarefac.otu_table1.biom -o rarefac.otu_table1.txt  --table-type "OTU table" --to-tsv
# rarefac.otu_table1.xls
/usr/local/R-3.6.0/bin/R -e '.libPaths("/usr/local/R-3.6.0/lib64/R/library"); library(tidyverse);library(data.table); bb <- fread("rarefac.otu_table1.txt") %>% rename_at(1,~"OTU ID"); nn <- bb[order(bb$`OTU ID`),]; arrange(nn,nchar(nn$`OTU ID`)) %>% rename_at(1,~"ASV ID") %>% fwrite("rarefac.otu_table1.xls",sep = "\t",scipen = 999)'

# rarefac.otu_table.xls 去掉了合并后行和为0的行
/usr/local/R-3.6.0/bin/R -e '.libPaths("/usr/local/R-3.6.0/lib64/R/library"); library(tidyverse);library(data.table);
rarefac.otu_table1 <- fread("rarefac.otu_table1.xls");rarefac.otu_table2 <- fread("otu_table2.xls")
rarefac.otu_table <- full_join(rarefac.otu_table1,rarefac.otu_table2);rarefac.otu_table[is.na(rarefac.otu_table)] <- 0
rarefac.otu_table %>% mutate(sum=apply(rarefac.otu_table[,-1],1,sum)) %>% filter(sum != 0) %>% select(-sum) %>%
  fwrite("rarefac.otu_table.xls",sep = "\t")
'
else
single_rarefaction.py -i otu_table.biom -o rarefac.otu_table.biom -d `less subsample.num`
# rarefac.otu_table.xls
biom convert -i rarefac.otu_table.biom -o rarefac.otu_table.txt  --table-type "OTU table" --to-tsv
/usr/local/R-3.6.0/bin/R -e '.libPaths("/usr/local/R-3.6.0/lib64/R/library"); library(tidyverse);library(data.table); bb <- fread("rarefac.otu_table.txt") %>% rename_at(1,~"OTU ID"); nn <- bb[order(bb$`OTU ID`),]; arrange(nn,nchar(nn$`OTU ID`)) %>% rename_at(1,~"ASV ID") %>% fwrite("rarefac.otu_table.xls",sep = "\t",scipen = 999)'
fi

# rarefac.otu_table.biom
biom convert -i rarefac.otu_table.xls -o rarefac.otu_table.biom --table-type="OTU table" --to-json

## rarefac.otu_taxa_table.xls
/usr/local/R-3.6.0/bin/R -e '.libPaths("/usr/local/R-3.6.0/lib64/R/library");
library(tidyverse);library(dplyr);library(data.table);
ASV_tax <- fread("otu_reps_tax_assignments.txt") %>% rename_all(~c("ASV ID","taxonomy","V3"))
rarefac.otu_table <- fread("rarefac.otu_table.xls") %>% left_join(ASV_tax[,-3]) %>% rename_at(1,~"ASV ID") %>% fwrite("rarefac.otu_taxa_table.xls",sep = "\t")
'

## rarefac.otu_genus.xls
biom convert -i rarefac.otu_taxa_table.xls -o rarefac.otu_taxa_table.biom --table-type "OTU table" --process-obs-metadata taxonomy --to-hdf5
vip_taxa_table_format.pl rarefac.otu_taxa_table.xls > rarefac.otu_genus.xls

summarize_taxa.py -i rarefac.otu_taxa_table.biom -o tax_summary_a -L 1,2,3,4,5,6,7,8 -a
summarize_taxa.py -i rarefac.otu_taxa_table.biom -o tax_summary_r -L 1,2,3,4,5,6,7,8
for ((i=1;i<=8;i+=1)){
       less tax_summary_a/rarefac.otu_taxa_table_L$i.txt|sed 's/Other/Unclassified/g'|grep -v 'Constructed from biom file'|sed 's/#OTU ID/Taxon/' >tax_summary_a/rarefac.otu_taxa_table_L$i.txt.1
       mv tax_summary_a/rarefac.otu_taxa_table_L$i.txt.1 tax_summary_a/rarefac.otu_taxa_table_L$i.txt
       less tax_summary_r/rarefac.otu_taxa_table_L$i.txt|sed 's/Other/Unclassified/g' |grep -v 'Constructed from biom file'|sed 's/#OTU ID/Taxon/' >tax_summary_r/rarefac.otu_taxa_table_L$i.txt.1
       mv tax_summary_r/rarefac.otu_taxa_table_L$i.txt.1 tax_summary_r/rarefac.otu_taxa_table_L$i.txt
}

levels=("L1" "L2" "L3" "L4" "L5" "L6" "L7" )
names=("kindom" "phylum" "class" "order" "family" "genus" "species")
for i in ${!levels[@]}; do
    level=${levels[$i]}
    name=${names[$i]}
sum_tax.pl -i tax_summary_a/rarefac.otu_taxa_table_${level}.txt -o tax_summary_a/${name}.xls.tmp
OTU_table_sortBySam.pl  tax_summary_a/${name}.xls.tmp sample_order tax_summary_a/${name}.xls
sum_tax.pl -i tax_summary_r/rarefac.otu_taxa_table_${level}.txt -o tax_summary_r/${name}.percents.xls.tmp
OTU_table_sortBySam.pl  tax_summary_r/${name}.percents.xls.tmp sample_order tax_summary_r/${name}.percents.xls
rm tax_summary_a/${name}.xls.tmp
rm tax_summary_r/${name}.percents.xls.tmp
done

#less otu_reps.fasta |awk -F '_' '{print $1}' > otu_reps.raw.fasta
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Make_dm.R -i rarefac.otu_table.xls -p otu_reps.raw_aligned.fasta.tre -m none

#### Alpha Diversity & Rarefaction curves  ###############################################
sum=$(awk '{sum += $2} END {print sum}' rarefac.otu_table.xls)
if [ "$sum" -eq 1000000 ]; then
  python /work/users/chaoliu/scripts/otu2shared.py -i otu_table.xls -o otus.shared
else
  python /work/users/chaoliu/scripts/otu2shared.py -i rarefac.otu_table.xls -o otus.shared
fi
mkdir alpha_rarefac
cp otus.shared alpha_rarefac/
cd alpha_rarefac
mothur "#summary.single(shared=otus.shared,calc=ace-chao-shannon-simpson-coverage,groupmode=f)"
mothur "#rarefaction.single(shared=otus.shared,calc=sobs-chao-shannon-simpson,groupmode=f,freq=100,processors=30)"

python /work/users/chaoliu/scripts/plot_alpha_diversity.py -i rarefaction
python /work/users/chaoliu/scripts/plot_alpha_diversity.py -i r_shannon

shannon-ace-table.pl -d . -o estimators.html
cd ../

#### picrust2 ###################################################
source /root/anaconda3/etc/profile.d/conda.sh
conda activate picrust2

picrust2_pipeline.py -s otu_reps.fasta -i rarefac.otu_table.biom -o Picrust2 --processes 30 --in_traits COG,EC,KO,PFAM,TIGRFAM

cd Picrust2
find . -name "*.gz" | xargs gunzip

add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv -m EC  -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv
add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv -m KO  -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv            
add_descriptions.py -i COG_metagenome_out/pred_metagenome_unstrat.tsv -m COG  -o COG_metagenome_out/pred_metagenome_unstrat_descrip.tsv           
add_descriptions.py -i PFAM_metagenome_out/pred_metagenome_unstrat.tsv -m PFAM  -o PFAM_metagenome_out/pred_metagenome_unstrat_descrip.tsv             
add_descriptions.py -i TIGRFAM_metagenome_out/pred_metagenome_unstrat.tsv -m TIGRFAM  -o TIGRFAM_metagenome_out/pred_metagenome_unstrat_descrip.tsv             
add_descriptions.py -i pathways_out/path_abun_unstrat.tsv -m METACYC -o pathways_out/path_abun_unstrat_descrip.tsv

# 生成kegg pathway 丰度表 # https://github.com/picrust/picrust2.git
pathway_pipeline.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv \
    -o KEGG_pathways_out --no_regroup \
    --map /work/softwares/picrust2-2.3.0-b/picrust2/default_files/pathway_mapfiles/KEGG_pathways_to_KO.tsv
# 添加功能描述
add_descriptions.py -i KEGG_pathways_out/path_abun_unstrat.tsv.gz \
    --custom_map_table /work/softwares/picrust2-2.3.0-b/picrust2/default_files/description_mapfiles/KEGG_pathways_info.tsv.gz \
    -o KEGG_pathways_out/path_abun_unstrat_descrip.tsv

# cp KEGG_modules_info.tsv.gz KEGG_modules_info_377.tsv.gz
# gunzip KEGG_modules_info_377.tsv.gz
# gzip KEGG_modules_info_377.tsv
# 生成kegg modules 丰度表 # https://github.com/picrust/picrust2.git
pathway_pipeline.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv \
    -o KEGG_modules_out --no_regroup \
    --map /work/softwares/picrust2-2.3.0-b/picrust2/default_files/pathway_mapfiles/KEGG_modules_to_KO.tsv
# 添加功能描述
add_descriptions.py -i KEGG_modules_out/path_abun_unstrat.tsv.gz \
    --custom_map_table /work/softwares/picrust2-2.3.0-b/picrust2/default_files/description_mapfiles/KEGG_modules_info_377.tsv.gz \
    -o KEGG_modules_out/path_abun_unstrat_descrip.tsv

cd ..
##### Tax4Fun2 ########################
source /root/anaconda3/etc/profile.d/conda.sh
conda deactivate
/usr/local/R-3.6.0/bin/R -e '.libPaths("/usr/local/R-3.6.0/lib64/R/library");pacman::p_load(tidyverse,stringr,magrittr); library(Tax4Fun2);runRefBlast(path_to_otus = "otu_reps.fasta", path_to_reference_data = "/work/users/chaoliu/database/Tax4Fun2_ReferenceData_v2", path_to_temp_folder = "Tax4Fun2", database_mode = "Ref99NR", use_force = T, num_threads = 30);makeFunctionalPrediction(path_to_otu_table = "otu_table.xls", path_to_reference_data = "/work/users/chaoliu/database/Tax4Fun2_ReferenceData_v2", path_to_temp_folder = "Tax4Fun2", database_mode = "Ref99NR", normalize_by_copy_number = TRUE, min_identity_to_reference = 0.97, normalize_pathways = FALSE)'

######## collect results  ##########
## now in analysis_Tue_Dec__1_18_22_16_2020/process/
cp tax_summary_a/*.xls          ../../results/Community/
cp tax_summary_r/*.percents.xls ../../results/Community/
cp -r otu_table.biom otu_table.xls otu_reps.fasta  otu_taxa_table.xls  rarefac.otu_table.biom rarefac.otu_table.xls rarefac.otu_taxa_table.biom rarefac.otu_taxa_table.xls tax_summary_a tax_summary_r rarefac.otu_genus.xls ../../results/OTU_Taxa/

cd alpha_rarefac/
cp otus.*.summary estimators.html ../../../results/Estimators
rename summary summary.xls ../../../results/Estimators/*.summary
cp -r otus.*.rarefaction rarefaction.*.pdf ../../../results/Rarefactions/
rename rarefaction rarefaction.xls ../../../results/Rarefactions/*.rarefaction
cp -r otus.*.r_shannon r_shannon.*.pdf ../../../results/Shannon_rarefac/
rename r_shannon r_shannon.xls ../../../results/Shannon_rarefac/*.r_shannon
cd ../../

cp ASV/Beta_diversity/*.txt ../results/Beta_diversity/
cp -r ASV/Picrust2/*_out ../results/Picrust2/
cp -r ASV/Tax4Fun2/*_prediction.txt  ../results/Tax4Fun2/
cd ..
```

## Part2.Masslin3
```R
library(optparse)
library(tidyverse)
#library(Maaslin2)
#library(Maaslin2, lib.loc = "/work/users/chaoliu/R/x86_64-pc-linux-gnu-library/3.6/")
library(reshape2)
library(pheatmap)
library(cowplot)
library(data.table)

# 从Bioconductor安装
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("maaslin3")
library(maaslin3)

if (TRUE) {
  option_list <- list(
    make_option(c("-i", "--input"),            type = "character", default = "inner_otu.genus.xls", help = "丰度表格;rarefac.otu_genus.xls"),
    make_option(c("-s", "--select"),           type = "character", default = "none",                help = "筛选表格:select.txt挑选物种"),
    make_option(c("-m", "--map"),              type = "character", default = "../map-group.txt",       help = "分组文件:map-group.txt"),
    make_option(c("-e", "--env"),              type = "character", default = "../age-gender-area-depth.txt",      help = "生理数据:env.txt"),
    make_option(c("-f", "--unif"),             type = "logical",   default = F,                     help="要不要归一化"),
    make_option(c("-p", "--per"),              type = "double",    default = 0,                     help="丰度筛选"),
    make_option(c("-v", "--value"),            type = "double",    default = 0,                   help="核心微生物筛选"),
    make_option(c("-j", "--pvalue"),           type = "double",    default = 0.05,                  help = "pvalue"),
    
    #maaslin3内部参数
    make_option(c("-o", "--output"),           type = "character", default = "./",                  help = "输出文件夹名"),
    make_option(c("--formula"),                type = "character", default = "~ age + gender + region + age:gender + age:region",                  help = "输出文件夹名"),
    #make_option(c("--fixed_effects"),          type = "character", default = c("age", "gender", "region"),   help = "固定效应"),
    make_option(c("--reference"),              type = "character", default = c("gender,F","region,central"),            help = "对照组 组名+对照组"),
    make_option(c("--min_abundance"),          type = "double",    default = 0,                  help = "校正后的q值小于"),
    make_option(c("--min_prevalence"),         type = "double",    default = 0.0001,                  help = "校正后的q值小于"),
    make_option(c("--zero_threshold"),         type = "double",    default = 0,                  help = "校正后的q值小于"),
    make_option(c("--min_variance"),           type = "double",    default = 0.1,                  help = "校正后的q值小于"),
    make_option(c("--max_significance"),       type = "double",    default = 0.1,                  help = "校正后的q值小于"),
    make_option(c("--normalization"),          type = "character", default = "CLR",                 help = "标准化方法TSS（总求和标准化）, CLR（中心对数比）,NONE（不进行标准化）"),
    make_option(c("--transform"),              type = "character", default = "NONE",                help = "LOG（base 2）,PLOG,NONE"),
    make_option(c("--standardize"),            type = "logical",   default = F,                     help = "要不要标准化"),
    make_option(c("--warn_prevalence"),        type = "logical",   default = F,                     help = "仅丰度分析就是FALSE；真实流行率效应和丰度效应就是TRUE"),
    make_option(c("--evaluate_only"),          type = "character", default = NULL,           help = "仅评估丰度（abundance）模型还是流行率（prevalence）模型 默认是NULL"),
    make_option(c("--plot_summary_plot"),      type = "logical",   default = T,                     help = "显著关联汇总图summary_plot.pdf"),
    make_option(c("--summary_plot_first_n"),   type = "double",    default = 25,                    help = "展示前几个特征"),
    make_option(c("--coef_plot_vars"),         type = "character", default = NULL,                  help = "NULL或者c('gender M','gender F')"),
    make_option(c("--heatmap_vars"),           type = "character", default = NULL,                  help = "展不展示热图NULL或者c('gender M','gender F')"),
    make_option(c("--plot_associations"),      type = "logical",   default = T,                     help = "散点图"),
    make_option(c("--max_pngs"),               type = "double",    default = 30,                    help = "展示前几个特征的png"),
    make_option(c("--cores"),                  type = "double",    default = 1,                     help = "使用 4 个核心并行计算"),
    make_option(c("--summary_plot_balanced"),  type = "logical",   default = F,                     help = "平均分配coef_plot_vars")
    
  )
  opts <- parse_args(OptionParser(option_list = option_list))
}

############################ Read in #################
Map <- read_tsv(opts$map) %>%
  rename_at(c(1), ~"SampleID") %>%
  mutate(group = fct_inorder(group))

# 行为样本，列为菌群
#SampleID,Bacteroides,Prevotella,Lactobacillus,Roseburia,Clostridium,Streptococcus
#S1,1500,200,50,300,0,100
#S2,3000,50,150,0,400,80
#S3,1200,300,30,250,100,200
#S4,800,150,200,100,50,0

genus <- read_tsv(opts$input) %>% 
  rename_at(1,~"SampleID") %>% 
  .[,c("SampleID",Map$SampleID)]


#genus_longer <- pivot_longer(genus, 
#             cols = colnames(genus)[-1],  # 明确指定需要转换的列
#             names_to = "Sample",       # 新列名，用于存储原始列名
#             values_to = "Value") 

#genus_wider <- pivot_wider(genus_longer, 
#            id_cols = Sample,       # 保留的列
#            names_from = SampleID,       # 从"Sample"列获取新列名
#            values_from = Value)    

#genus <- genus_wider %>% rename_at(1,~"SampleID")
# Uniform
if (opts$unif){
  genus[,2:ncol(genus)] <- sapply(2:ncol(genus),function(x) genus[,x]/sum(genus[,x]))
}

# 丰度筛选
if(opts$per != 0){
  genus <- genus %>% 
    mutate(SELECT = sapply(1:nrow(.),function(x){
      any(.[x,2:ncol(.)]>=opts$per)})) %>% 
    .[.$SELECT,] %>% as_tibble %>% dplyr::select(-SELECT)
}


# 核心微生物筛选
if(opts$value != 0){
  otu_coverage <- apply(genus[,2:ncol(genus)],1,function(x) 
    length(x[x>0])/(ncol(genus)-1))
  genus <- genus[otu_coverage >= opts$value,]
}


## # 调整ASV名
## if(opts$input == "../final.CLR.new.otu_taxa_table.xls"){
##   ASV_name <- fread("../../5.分段年龄性别地区/ASV/inner_otu.genus.xls") %>% .[,1]
##   rename <- separate(ASV_name, "OTU_tax", into = c("ASV", "Taxonomy"), sep = " ", remove = FALSE) %>% 
##     select(-Taxonomy) %>% rename_at(2,~"SampleID")
##   genus0 <- genus
##   genus <- dplyr::left_join(rename,genus) %>% select(-SampleID) %>% rename_at(1,~"SampleID")
## }


# 物中挑选
if(opts$select != "none"){
  ss <- read_tsv(opts$select,col_names = F) %>% rename_at(1,~"SampleID")
  genus <- inner_join(genus,ss)
}

Env <- read_tsv(opts$env) %>%
  rename_at(c(1), ~"SampleID") %>%
  dplyr::inner_join(., Map) %>% select(-group) %>% rename_at(4,~"region")


################################################################################
# 在Maaslin3中，输入数据的格式要求是：
# 行 = 样本
# 列 = 特征（物种）

data1 <- as.matrix(genus[,-1] %>% t)
colnames(data1) <- genus$SampleID

data2 <- as.data.frame(Env[, -1])
rownames(data2) <- Env$SampleID

#data2 <- Map
#rownames(data2) <- Map$SampleID
#data1 <- data1[,colSums(data1) != 0]

####################### maaslin3 #####################
set.seed(20190731)

fit_data3 <- maaslin3::maaslin3(
  input_data = data1,                   # 微生物丰度数据
  input_metadata = data2,               # 环境因子数据
  output = opts$output,      # 输出路径
  #fixed_effects = opts$fixed_effects,              # 固定效应 自变量：年龄分组（分类变量）
  formula =  opts$formula ,
  #formula =  ~ age + gender + region + age:gender + age:region ,
  reference = opts$reference,                     #对照组 不设置就默认第一个
  min_abundance = opts$min_abundance   ,                 # 经过normalization参数指定的方法处理后的值
  min_prevalence = opts$min_prevalence   ,                 # 经过normalization参数指定的方法处理后的值
  zero_threshold = opts$zero_threshold   ,                 # 经过normalization参数指定的方法处理后的值
  min_variance = opts$min_variance   ,                 # 经过normalization参数指定的方法处理后的值
  max_significance = opts$max_significance  ,               # 校正后的q值小于等于opts$pvalue 默认0.1
  #max_significance = 1  ,               # 校正后的q值小于等于opts$pvalue 默认0.1
  normalization = opts$normalization,                # 标准化方法（推荐 CLR 处理组成性数据）
  transform = opts$transform,                   # 数据已转换则无需再变换
  standardize = opts$standardize,                 # 不进行标准化
  warn_prevalence = opts$warn_prevalence,              # 仅丰度分析就是FALSE；真实流行率效应和丰度效应就是TRUE
  evaluate_only = opts$evaluate_only,          # 仅评估丰度（“abundance”）模型还是流行率（“prevalence”）模型
  plot_summary_plot = opts$plot_summary_plot ,               # 显著关联汇总图summary_plot.pdf
  summary_plot_first_n = opts$summary_plot_first_n,            # 展示前几个特征
  #coef_plot_vars = opts$coef_plot_vars,       # 每次手动调一下 记得中间有空格
  #heatmap_vars = opts$heatmap_vars,          # 热图
  plot_associations = opts$plot_associations           ,     # 画散点图
  max_pngs = opts$max_pngs                   ,    #排名前max_pngs的会用png展示出来
  cores = opts$cores,                           # 使用 24 个核心并行计算
  summary_plot_balanced = opts$summary_plot_balanced         
)


#原始 p 值（pval_individual）未经过多重检验校正，可能高估显著性，因此优先以校正后的qval_individual作为判断标准。
#常用0.25
#严格0.05
#宽松0.5
#然后再结合coef

#summary_plot_balanced
#展示coef_plot_vars中每个变量的前 N 个特征，其中 N 等于：ceiling (summary_plot_first_n /length (coef_plot_vars))
#当分析涉及多个变量（如同时纳入age、genderF、genderM等）时，默认情况下summary_plot_first_n（如 25）会展示所有变量中排名前 25 的特征，可能导致某一变量占据多数展示位置，其他变量被忽略。
#而summary_plot_balanced = TRUE会让每个变量 “公平分配” 展示名额
#例如，若summary_plot_first_n = 25且coef_plot_vars包含 3 个变量，则每个变量展示ceiling(25/3) = 9个特征（总展示 27 个，略多于 25），确保每个变量都有足够的特征被呈现。
#使用前提：
#必须同时设置coef_plot_vars（指定要展示的变量），否则会报错。例如：

#coef_plot_vars = c("age", "genderF", "genderM"),  # 3个变量
#summary_plot_first_n = 25,
#summary_plot_balanced = TRUE  # 每个变量展示9个特征



#适用场景：
#适合多变量分析（如同时研究年龄、性别等），希望在汇总图中均衡展示每个变量的重要关联特征，避免某一变量的特征 “垄断” 图表，更全面地呈现不同变量的关联模式。




fit_data3$fit_data_abundance$results %>%
  as_tibble() %>% # group_by(feature) %>% nest %>%
  #filter(pval_individual < opts$pvalue) %>%
  filter(qval_individual < opts$max_significance) %>%
  pull(feature) %>%
  unique() %>%
  as_tibble() %>%
  write_tsv(str_c(opts$output, "maaslin3.select.list"), col_names = FALSE)

##########################################################################################
# 找到最长参数的长度
max_len <- max(sapply(opts, length))

# 每个参数填充至相同长度
opts_padded <- lapply(opts, function(x) {
  max_len <- max(sapply(opts, length))  # 按最长参数长度对齐
  if (length(x) < max_len) c(x, rep(NA, max_len - length(x))) else x
})

write_tsv(opts_padded %>% as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column() %>% as_tibble() %>% mutate(across(everything(), ~replace_na(., ""))),
          str_c(
            "Parameter",
            str_replace_all(as.character(date()), " ", "_") %>% str_replace_all(":", "_"),
            ".xls"
          ),
          col_names = FALSE
)

```
## Part3. ANCOM-bc2
```R
library(optparse)
library(tidyverse)
#library(Maaslin2)
#library(Maaslin2, lib.loc = "/work/users/chaoliu/R/x86_64-pc-linux-gnu-library/3.6/")
library(reshape2)
library(pheatmap)
library(cowplot)
library(data.table)

# 从Bioconductor安装
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("maaslin3")
library(maaslin3)

if (TRUE) {
  option_list <- list(
    make_option(c("-i", "--input"),            type = "character", default = "inner_otu.genus.xls", help = "丰度表格;rarefac.otu_genus.xls"),
    make_option(c("-s", "--select"),           type = "character", default = "none",                help = "筛选表格:select.txt挑选物种"),
    make_option(c("-m", "--map"),              type = "character", default = "../map-group.txt",       help = "分组文件:map-group.txt"),
    make_option(c("-e", "--env"),              type = "character", default = "../age-gender-area-depth.txt",      help = "生理数据:env.txt"),
    make_option(c("-f", "--unif"),             type = "logical",   default = F,                     help="要不要归一化"),
    make_option(c("-p", "--per"),              type = "double",    default = 0,                     help="丰度筛选"),
    make_option(c("-v", "--value"),            type = "double",    default = 0,                   help="核心微生物筛选"),
    make_option(c("-j", "--pvalue"),           type = "double",    default = 0.05,                  help = "pvalue"),
    
    #maaslin3内部参数
    make_option(c("-o", "--output"),           type = "character", default = "./",                  help = "输出文件夹名"),
    make_option(c("--formula"),                type = "character", default = "~ age + gender + region + age:gender + age:region",                  help = "输出文件夹名"),
    #make_option(c("--fixed_effects"),          type = "character", default = c("age", "gender", "region"),   help = "固定效应"),
    make_option(c("--reference"),              type = "character", default = c("gender,F","region,central"),            help = "对照组 组名+对照组"),
    make_option(c("--min_abundance"),          type = "double",    default = 0,                  help = "校正后的q值小于"),
    make_option(c("--min_prevalence"),         type = "double",    default = 0.0001,                  help = "校正后的q值小于"),
    make_option(c("--zero_threshold"),         type = "double",    default = 0,                  help = "校正后的q值小于"),
    make_option(c("--min_variance"),           type = "double",    default = 0.1,                  help = "校正后的q值小于"),
    make_option(c("--max_significance"),       type = "double",    default = 0.1,                  help = "校正后的q值小于"),
    make_option(c("--normalization"),          type = "character", default = "CLR",                 help = "标准化方法TSS（总求和标准化）, CLR（中心对数比）,NONE（不进行标准化）"),
    make_option(c("--transform"),              type = "character", default = "NONE",                help = "LOG（base 2）,PLOG,NONE"),
    make_option(c("--standardize"),            type = "logical",   default = F,                     help = "要不要标准化"),
    make_option(c("--warn_prevalence"),        type = "logical",   default = F,                     help = "仅丰度分析就是FALSE；真实流行率效应和丰度效应就是TRUE"),
    make_option(c("--evaluate_only"),          type = "character", default = NULL,           help = "仅评估丰度（abundance）模型还是流行率（prevalence）模型 默认是NULL"),
    make_option(c("--plot_summary_plot"),      type = "logical",   default = T,                     help = "显著关联汇总图summary_plot.pdf"),
    make_option(c("--summary_plot_first_n"),   type = "double",    default = 25,                    help = "展示前几个特征"),
    make_option(c("--coef_plot_vars"),         type = "character", default = NULL,                  help = "NULL或者c('gender M','gender F')"),
    make_option(c("--heatmap_vars"),           type = "character", default = NULL,                  help = "展不展示热图NULL或者c('gender M','gender F')"),
    make_option(c("--plot_associations"),      type = "logical",   default = T,                     help = "散点图"),
    make_option(c("--max_pngs"),               type = "double",    default = 30,                    help = "展示前几个特征的png"),
    make_option(c("--cores"),                  type = "double",    default = 1,                     help = "使用 4 个核心并行计算"),
    make_option(c("--summary_plot_balanced"),  type = "logical",   default = F,                     help = "平均分配coef_plot_vars")
    
  )
  opts <- parse_args(OptionParser(option_list = option_list))
}

############################ Read in #################
Map <- read_tsv(opts$map) %>%
  rename_at(c(1), ~"SampleID") %>%
  mutate(group = fct_inorder(group))

# 行为样本，列为菌群
#SampleID,Bacteroides,Prevotella,Lactobacillus,Roseburia,Clostridium,Streptococcus
#S1,1500,200,50,300,0,100
#S2,3000,50,150,0,400,80
#S3,1200,300,30,250,100,200
#S4,800,150,200,100,50,0

genus <- read_tsv(opts$input) %>% 
  rename_at(1,~"SampleID") %>% 
  .[,c("SampleID",Map$SampleID)]


#genus_longer <- pivot_longer(genus, 
#             cols = colnames(genus)[-1],  # 明确指定需要转换的列
#             names_to = "Sample",       # 新列名，用于存储原始列名
#             values_to = "Value") 

#genus_wider <- pivot_wider(genus_longer, 
#            id_cols = Sample,       # 保留的列
#            names_from = SampleID,       # 从"Sample"列获取新列名
#            values_from = Value)    

#genus <- genus_wider %>% rename_at(1,~"SampleID")
# Uniform
if (opts$unif){
  genus[,2:ncol(genus)] <- sapply(2:ncol(genus),function(x) genus[,x]/sum(genus[,x]))
}

# 丰度筛选
if(opts$per != 0){
  genus <- genus %>% 
    mutate(SELECT = sapply(1:nrow(.),function(x){
      any(.[x,2:ncol(.)]>=opts$per)})) %>% 
    .[.$SELECT,] %>% as_tibble %>% dplyr::select(-SELECT)
}


# 核心微生物筛选
if(opts$value != 0){
  otu_coverage <- apply(genus[,2:ncol(genus)],1,function(x) 
    length(x[x>0])/(ncol(genus)-1))
  genus <- genus[otu_coverage >= opts$value,]
}


## # 调整ASV名
## if(opts$input == "../final.CLR.new.otu_taxa_table.xls"){
##   ASV_name <- fread("../../5.分段年龄性别地区/ASV/inner_otu.genus.xls") %>% .[,1]
##   rename <- separate(ASV_name, "OTU_tax", into = c("ASV", "Taxonomy"), sep = " ", remove = FALSE) %>% 
##     select(-Taxonomy) %>% rename_at(2,~"SampleID")
##   genus0 <- genus
##   genus <- dplyr::left_join(rename,genus) %>% select(-SampleID) %>% rename_at(1,~"SampleID")
## }


# 物中挑选
if(opts$select != "none"){
  ss <- read_tsv(opts$select,col_names = F) %>% rename_at(1,~"SampleID")
  genus <- inner_join(genus,ss)
}

Env <- read_tsv(opts$env) %>%
  rename_at(c(1), ~"SampleID") %>%
  dplyr::inner_join(., Map) %>% select(-group) %>% rename_at(4,~"region")


################################################################################
# 在Maaslin3中，输入数据的格式要求是：
# 行 = 样本
# 列 = 特征（物种）

data1 <- as.matrix(genus[,-1] %>% t)
colnames(data1) <- genus$SampleID

data2 <- as.data.frame(Env[, -1])
rownames(data2) <- Env$SampleID

#data2 <- Map
#rownames(data2) <- Map$SampleID
#data1 <- data1[,colSums(data1) != 0]

####################### maaslin3 #####################
set.seed(20190731)

fit_data3 <- maaslin3::maaslin3(
  input_data = data1,                   # 微生物丰度数据
  input_metadata = data2,               # 环境因子数据
  output = opts$output,      # 输出路径
  #fixed_effects = opts$fixed_effects,              # 固定效应 自变量：年龄分组（分类变量）
  formula =  opts$formula ,
  #formula =  ~ age + gender + region + age:gender + age:region ,
  reference = opts$reference,                     #对照组 不设置就默认第一个
  min_abundance = opts$min_abundance   ,                 # 经过normalization参数指定的方法处理后的值
  min_prevalence = opts$min_prevalence   ,                 # 经过normalization参数指定的方法处理后的值
  zero_threshold = opts$zero_threshold   ,                 # 经过normalization参数指定的方法处理后的值
  min_variance = opts$min_variance   ,                 # 经过normalization参数指定的方法处理后的值
  max_significance = opts$max_significance  ,               # 校正后的q值小于等于opts$pvalue 默认0.1
  #max_significance = 1  ,               # 校正后的q值小于等于opts$pvalue 默认0.1
  normalization = opts$normalization,                # 标准化方法（推荐 CLR 处理组成性数据）
  transform = opts$transform,                   # 数据已转换则无需再变换
  standardize = opts$standardize,                 # 不进行标准化
  warn_prevalence = opts$warn_prevalence,              # 仅丰度分析就是FALSE；真实流行率效应和丰度效应就是TRUE
  evaluate_only = opts$evaluate_only,          # 仅评估丰度（“abundance”）模型还是流行率（“prevalence”）模型
  plot_summary_plot = opts$plot_summary_plot ,               # 显著关联汇总图summary_plot.pdf
  summary_plot_first_n = opts$summary_plot_first_n,            # 展示前几个特征
  #coef_plot_vars = opts$coef_plot_vars,       # 每次手动调一下 记得中间有空格
  #heatmap_vars = opts$heatmap_vars,          # 热图
  plot_associations = opts$plot_associations           ,     # 画散点图
  max_pngs = opts$max_pngs                   ,    #排名前max_pngs的会用png展示出来
  cores = opts$cores,                           # 使用 24 个核心并行计算
  summary_plot_balanced = opts$summary_plot_balanced         
)


#原始 p 值（pval_individual）未经过多重检验校正，可能高估显著性，因此优先以校正后的qval_individual作为判断标准。
#常用0.25
#严格0.05
#宽松0.5
#然后再结合coef

#summary_plot_balanced
#展示coef_plot_vars中每个变量的前 N 个特征，其中 N 等于：ceiling (summary_plot_first_n /length (coef_plot_vars))
#当分析涉及多个变量（如同时纳入age、genderF、genderM等）时，默认情况下summary_plot_first_n（如 25）会展示所有变量中排名前 25 的特征，可能导致某一变量占据多数展示位置，其他变量被忽略。
#而summary_plot_balanced = TRUE会让每个变量 “公平分配” 展示名额
#例如，若summary_plot_first_n = 25且coef_plot_vars包含 3 个变量，则每个变量展示ceiling(25/3) = 9个特征（总展示 27 个，略多于 25），确保每个变量都有足够的特征被呈现。
#使用前提：
#必须同时设置coef_plot_vars（指定要展示的变量），否则会报错。例如：

#coef_plot_vars = c("age", "genderF", "genderM"),  # 3个变量
#summary_plot_first_n = 25,
#summary_plot_balanced = TRUE  # 每个变量展示9个特征



#适用场景：
#适合多变量分析（如同时研究年龄、性别等），希望在汇总图中均衡展示每个变量的重要关联特征，避免某一变量的特征 “垄断” 图表，更全面地呈现不同变量的关联模式。




fit_data3$fit_data_abundance$results %>%
  as_tibble() %>% # group_by(feature) %>% nest %>%
  #filter(pval_individual < opts$pvalue) %>%
  filter(qval_individual < opts$max_significance) %>%
  pull(feature) %>%
  unique() %>%
  as_tibble() %>%
  write_tsv(str_c(opts$output, "maaslin3.select.list"), col_names = FALSE)

##########################################################################################
# 找到最长参数的长度
max_len <- max(sapply(opts, length))

# 每个参数填充至相同长度
opts_padded <- lapply(opts, function(x) {
  max_len <- max(sapply(opts, length))  # 按最长参数长度对齐
  if (length(x) < max_len) c(x, rep(NA, max_len - length(x))) else x
})

write_tsv(opts_padded %>% as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column() %>% as_tibble() %>% mutate(across(everything(), ~replace_na(., ""))),
          str_c(
            "Parameter",
            str_replace_all(as.character(date()), " ", "_") %>% str_replace_all(":", "_"),
            ".xls"
          ),
          col_names = FALSE
)

```
## Part4.Group_analysis
```shell
#!/bin/bash
#set -euxo pipefail
##########  argparser  ###########
## -c color.txt # 指定分组颜色
## -p TRUE      # 是否成对分析
## -b TRUE      # 是否排列组合
# 20220614 v3.0.1.6 修改了lefse脚本，可以指定颜色了
# v3.0.1.7 matrix_analysis合并adonis-anosim-mrpp-manova分析
# v3.0.1.8 plot_community_boxplot.R修改单物种箱式图脚本

until [ $# -eq 0 ]
do
  name=${1:1}; shift;
  if [[ -z "$1" || $1 == -* ]] ; then eval "export $name=true"; else eval "export $name=$1"; shift; fi
done
########## set default  ###########
if [ ${#c} == 0 ]; then Colour=none;  else Colour=`pwd`"/"$c; fi
if [ ${#p} == 0 ]; then Paired=FALSE; else Paired=$p;         fi
if [ ${#b} == 0 ]; then Combn=TRUE;   else Combn=$p;          fi
if [ ${#k} == 0 ]; then Paint=none;  else Paint=$k;          fi

echo "Define colors?"  $Colour
echo "Paired test?"    $Paired
echo "Combn annlysis?" $Combn
echo "Paint Bucket?"   $Paint

######################################################################################################################
a=`cat map-group.txt | wc -l`
if [ $a -ge 150 ]
then
    export Heatmap_w="30"
    export Heatmap_h="10"
    export Heatmap_keyh="14"
    export Heatmap_marble="4-0-0-20"
    export Hclust_bar_h="30"
    export Hcluster_tree_h="30"
elif [ $a -ge 100 ]
then
    export Heatmap_w="20"
    export Heatmap_h="10"
    export Heatmap_keyh="14"
    export Heatmap_marble="4-0-0-14"
    export Hclust_bar_h="20"
    export Hcluster_tree_h="20"
elif [ $a -ge 75 ]
then
    export Heatmap_w="16"
    export Heatmap_h="10"
    export Heatmap_keyh="13"
    export Heatmap_marble="4-0-0-14"
    export Hclust_bar_h="15"
    export Hcluster_tree_h="15"
elif [ $a -ge 50 ]
then
    export Heatmap_w="12"
    export Heatmap_h="10"
    export Heatmap_keyh="12"
    export Heatmap_marble="4-0-0-14"
    export Hclust_bar_h="10"
    export Hcluster_tree_h="10"
elif [ $a -ge 25 ]
then
    export Heatmap_w="10"
    export Heatmap_h="10"
    export Heatmap_keyh="10"
    export Heatmap_marble="4-0-0-10"
    export Hclust_bar_h="8"
    export Hcluster_tree_h="8"
elif [ $a -gt 0 ]
then
    export Heatmap_w="8"
    export Heatmap_h="9"
    export Heatmap_keyh="9"
    export Heatmap_marble="4-0-0-8"
    export Hclust_bar_h="5"
    export Hcluster_tree_h="7"
fi

# need order; map-group.txt groups.txt File
# run directory: results
cat map-group.txt | grep -v "SampleID.*group" | awk '{print $1"\t"$1}' >order
grep -v "SampleID.*group" map-group.txt | sed '1i#SampleID\tgroup' > map-group.txt.tmp
rm map-group.txt
mv map-group.txt.tmp map-group.txt
#cat map-group.txt |sed '1d'|awk -F"\t" '{if($2==last){group=group","$1}else{print group;group=$2":"$1};last=$2}END{print group}'|sed '1d' >map.txt
NUM=`awk -F "\t" '{if(ARGIND==1) {val[$1]}else{if($1 in val)  delete val[$1]}}END{for(i in val) print i}' order ../process/ASV/sample_order | wc -l`
if [ $NUM -gt 0 ]; then echo "command failed"; exit 1; fi

mkdir Split_groups
cd Split_groups
ln -s ../map-group.txt .
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/split_group.R -g map-group.txt -m none -n 0
cd ..

mkdir ASV_analysis
cd ASV_analysis
#------ 原代码 ------
#cp ../../process/ASV/rarefac.otu_taxa_table.xls .
#sed -i 's/ASV ID/OTU_ID/g' rarefac.otu_taxa_table.xls
# 然后用R把这个文件修改一下
#------ ------ ------

#------ 修改后 ------
cp ../../process/ASV/rarefac.otu_taxa_table.xls .
sed -i 's/ASV ID/OTU_ID/g' rarefac.otu_taxa_table.xls
ln -s /Mobio/users/yanlv/Project/wuzhongwen/20250115005/analysis_Thu_Jan_23_2025-01-23_15_38_26/4.group/4.7.Tax_test/3/otu_reps_tax_assignments.txt ./
R -e 'library(data.table)
rarefac.otu_taxa_table <- fread("rarefac.otu_taxa_table.xls")
fix_tax <- read_tsv("otu_reps_tax_assignments.txt") %>% .[,-3] %>% rename_all(~c("OTU_ID","taxonomy2"))
left_join(rarefac.otu_taxa_table,fix_tax) %>% mutate(taxonomy=taxonomy2) %>% select(-taxonomy2) %>%
  fwrite("fix.rarefac.otu_taxa_table.xls",sep = "\t")
'
#对比看一下 看下第一行最后面有没有制表符
#head -n 10 rarefac.otu_taxa_table.xls > 01.xls
#head -n 10 fix.rarefac.otu_taxa_table.xls > 02.xls
#rm 01.xls 02.xls
rm otu_reps_tax_assignments.txt rarefac.otu_taxa_table.xls
mv fix.rarefac.otu_taxa_table.xls rarefac.otu_taxa_table.xls
#------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------

cp ../../process/ASV/otu_reps.fasta .
cp ../map-group.txt .
less otu_reps.fasta |awk -F '_' '{print $1}' > otu_reps.raw.fasta

python /work/users/chaoliu/scripts/tax_split.v2.0.2.py -i rarefac.otu_taxa_table.xls -m map-group.txt -r otu_reps.raw.fasta -t 0 -f F

#mkdir ASV_tree
#cd ASV_tree
#ln -s ../map.otu_reps.raw.fasta ./ASV_reps.raw.fasta
#ln -s ../unrooted.otu_tree_anno.0.001.xls ./unrooted.ASV_tree_anno.0.001.xls
#muscle -in ASV_reps.raw.fasta -out ASV_reps.raw_aligned.fasta
#FastTree -nt ASV_reps.raw_aligned.fasta > unroot.ASV_reps.raw_aligned.fasta.tre
#/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/unroot_tree.R -t unroot.ASV_reps.raw_aligned.fasta.tre -a unrooted.ASV_tree_anno.0.001.xls -s phylum
#rm Rplots.pdf
#cd ..

cd Krona
for i in *.xls; do ktImportText $i -o $i.html;done
cd ..

mkdir Core_Microbiome
cd Core_Microbiome
ln -s ../../OTU_Taxa/rarefac.otu_genus.xls
ln -s ../map-group.txt .
sed 's/\t$//' rarefac.otu_genus.xls > fix.rarefac.ASV_genus.xls
# 此时此刻 刚刚生成的fix.rarefac.ASV_genus.xls的第一行最后面还是制表符 手动删一下
sed -i 's/ //g' fix.rarefac.ASV_genus.xls
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/core_microbiome.R -i fix.rarefac.ASV_genus.xls -m map-group.txt -c $Colour

cd ../../
mv ASV_analysis/Core_Microbiome  ASV_analysis/Krona ASV_analysis/ASV_tree .
#rm -rf OTU_analysis

rm -rf OTU_Taxa
mv ASV_analysis ASV_Taxa
cd ASV_Taxa
rm rarefac.otu_taxa_table.xls otu_reps.raw.fasta otu_reps.fasta
mv map.otu_reps.raw.fasta ASV_reps.raw.fasta
mv map.otu_table.xls rarefac.ASV_table.xls
mv map.otu_table.percent.xls rarefac.ASV_table.percent.xls
mv map.rarefac.otu_taxa_table.xls rarefac.ASV_taxa_table.xls
mv otu.genus.xls rarefac.ASV_genus.xls

/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -r $Paired -o F -i rarefac.ASV_genus.xls -n 2 -p 20  -m map-group.txt  -c $Colour --combn $Combn
cp ../results/OTU_Taxa/otu_table.xls ./
cp ../results/OTU_Taxa/otu_taxa_table.xls ./

cd ..

#Specaccum
cd Specaccum
ln -s ../ASV_Taxa/rarefac.ASV_table.xls .
ln -s ../map-group.txt .
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Specaccum.R -i rarefac.ASV_table.xls -m map-group.txt -c $Colour -d ASV
rm *.xls
cd ../

#Venn
mkdir Venn
cd Venn
ln -s ../ASV_Taxa/rarefac.ASV_genus.xls .
ln -s ../map-group.txt
Rscript /work/users/chaoliu/scripts/venn_upset.R -i rarefac.ASV_genus.xls -g map-group.txt -m all -c $Colour
cd ../

mkdir Corrplot
cd Corrplot
ln -s ../ASV_Taxa/rarefac.ASV_genus.xls .
ln -s ../map-group.txt .
Rscript /work/users/chaoliu/scripts/Spearman_corrplot.R -i rarefac.ASV_genus.xls -m map-group.txt -o ASV_top50
cd ..

#Rank_abundance
cd Rank_abundance
ln -s ../ASV_Taxa/rarefac.ASV_table.xls
perl /work/users/chaoliu/scripts/rank_abundance.pl -i rarefac.ASV_table.xls -gd ../map-group.txt  -o rankabundance.group.pdf -w 6 -h 5 -color $Colour -mode ASV
rm rarefac.ASV_table.xls
cd ../

#Alpha_rarefac
mkdir Alpha_rarefac
cp ../process/ASV/alpha_rarefac/* Alpha_rarefac -r
cd Alpha_rarefac/
rm map-group.txt
#这个之前出过图 把之前出的图和表都删了
ln -s ../map-group.txt ./
rename 's/\.rarefaction/\.rarefaction_asv/g' otus.*.rarefaction
#python /work/users/chaoliu/scripts/plot_alpha_diversity.py -i rarefaction_asv
#python /work/users/chaoliu/scripts/plot_alpha_diversity.py -i r_shannon
python /work/users/chaoliu/scripts/plot_alpha_diversity.py -i rarefaction_asv -m map-group.txt -c $Colour
python /work/users/chaoliu/scripts/plot_alpha_diversity.py -i r_shannon -m map-group.txt -c $Colour
rm ../Rarefactions/*.pdf
rm ../Shannon_rarefac/*.pdf
mv rarefaction_* ../Rarefactions
mv r_shannon_* ../Shannon_rarefac

ln -s ../map-group.txt .
rename 's/\.rarefaction_asv/\.rarefaction/g' otus.*.rarefaction_asv
alpha_rarefac.pl map-group.txt > alpha_rarefac.summary.xls
sed -i 's/observed_otus/observed_asvs/'  alpha_rarefac.summary.xls

#Wilcox-alpha
#/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/alpha_diversity_test.R -a alpha_rarefac.summary.xls -g map-group.txt -t $Paired
/usr/local/R-4.3.1/bin/Rscript /work/users/yuren/Scripts/alpha_diversity_test.pro.R -a alpha_rarefac.summary.xls -g map-group.txt -t unpara -r $Paired --pv 0.05 --pj 1 --combn $Combn

#/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Alpha_diversity_box_plot.R -a alpha_rarefac.summary.xls -g map-group.txt -c $Colour
/usr/local/R-4.3.1/bin/Rscript /work/users/yuren/Scripts/plot_index_boxplot.pro.R -i alpha_rarefac.summary.xls -m map-group.txt -c $Colour -d alpha --test Alpha_diversity.combn.xls -f FALSE --bb cloud -l 0 --cores 1
rm *.r mothur* otus.* *.txt
cd ../

#Hcluster_tree 
cd Hcluster_tree 
ln -s ../Beta_diversity/bray_curtis_dm.txt ../Beta_diversity/unweighted_unifrac_dm.txt ../Beta_diversity/weighted_unifrac_dm.txt ../map-group.txt ./ 
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plot_tree.R -i bray_curtis_dm.txt -m map-group.txt -t average -c $Colour
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plot_tree.R -i unweighted_unifrac_dm.txt -m map-group.txt -t average -c $Colour
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plot_tree.R -i weighted_unifrac_dm.txt -m map-group.txt -t average -c $Colour
rm *.txt *.r Rplots.pdf
cd ../

#Beta_diversity
cd Beta_diversity
ln -s ../map-group.txt
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Matrix_analysis.R -i bray_curtis_dm.txt        -g map-group.txt  -c $Colour  -m $Combn
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Matrix_analysis.R -i unweighted_unifrac_dm.txt -g map-group.txt  -c $Colour  -m $Combn
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Matrix_analysis.R -i weighted_unifrac_dm.txt   -g map-group.txt  -c $Colour  -m $Combn
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Matrix_analysis.R -i jaccard-binary_dm.txt     -g map-group.txt  -c $Colour  -m $Combn
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Matrix_analysis.R -i euclidean_dm.txt          -g map-group.txt  -c $Colour  -m $Combn

/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plot_matrix.R -j FALSE -i bray_curtis_dm.txt        -g map-group.txt -c $Colour
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plot_matrix.R -j FALSE -i unweighted_unifrac_dm.txt -g map-group.txt -c $Colour
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plot_matrix.R -j FALSE -i weighted_unifrac_dm.txt   -g map-group.txt -c $Colour
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plot_matrix.R -j FALSE -i jaccard-binary_dm.txt     -g map-group.txt -c $Colour
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plot_matrix.R -j FALSE -i euclidean_dm.txt          -g map-group.txt -c $Colour

mkdir Anosim Adonis MRPP MANOVA
mv *anosim* Anosim
mv *adonis* Adonis
mv *mrpp* MRPP
mv *MANOVA* MANOVA
cd ..

mkdir Beta_diversity_clr
cd Beta_diversity_clr
ln -s ../map-group.txt ./
ln -s ../../process/ASV/ASV_absolute_analysis/aitchison_dm.txt ./
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Matrix_analysis.R -i aitchison_dm.txt          -g map-group.txt  -c $Colour  -m $Combn
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plot_matrix.R -j FALSE -i aitchison_dm.txt          -g map-group.txt -c $Colour

mkdir Anosim Adonis MRPP MANOVA
mv *anosim* Anosim
mv *adonis* Adonis
mv *mrpp* MRPP
mv *MANOVA* MANOVA

# 创建不同可视化目录
for dir in noname 3D name ellipse box crossbar central; do
  mkdir $dir
done

# noname - 基础PCoA图
for pc in "1-2" "1-3" "2-3"; do
  python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i aitchison_dm.txt -md PCoA -pc $pc -map map-group.txt -col $Colour
done
mv *.pdf noname/

# 3D - 三维PCoA图
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plotly_PCoA-3D.R -i aitchison_dm.txt -c $Colour
mv *.pdf 3D/

# name - 带样本名的PCoA图
for pc in "1-2" "1-3" "2-3"; do
  python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i aitchison_dm.txt -md PCoA -pc $pc -map map-group.txt -col $Colour
done
mv *.pdf name/

# ellipse - 带置信椭圆的PCoA图
for pc in "1-2" "1-3" "2-3"; do
  python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i aitchison_dm.txt -md PCoA -pc $pc -map map-group.txt -col $Colour
done
mv *.pdf ellipse/

# box - 带箱线图和置信椭圆的PCoA图
for pc in "1-2" "1-3" "2-3"; do
  python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i aitchison_dm.txt -md PCoA -pc $pc -map map-group.txt -col $Colour
done
mv *.pdf box/

# crossbar - 带十字线的PCoA图
for pc in "1-2" "1-3" "2-3"; do
  python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i aitchison_dm.txt -ct 2 -md PCoA -pc $pc -map map-group.txt -col $Colour
done
mv *.pdf crossbar/

# central - 带中心点的PCoA图
for pc in "1-2" "1-3" "2-3"; do
  python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i aitchison_dm.txt -ct 1 -md PCoA -pc $pc -map map-group.txt -col $Colour
done
mv *.pdf central/

# NMDS
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i aitchison_dm.txt   -md NMDS -pc 1-2 -map map-group.txt -col $Colour

rm cmd.r 
cd ../

#Community
rm -rf Community/
mv ASV_Taxa/Community .
cd Community
ln -s ../order ./sample
sed -i 's/ //g' *.xls
if [[ -f "../paint.txt" ]];then ln -s ../paint.txt ;fi

mkdir Community_barplot
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i phylum.xls  -n 2 -p 50 --paint $Paint
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i class.xls   -n 2 -p 50 --paint $Paint
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i order.xls   -n 2 -p 50 --paint $Paint
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i family.xls  -n 2 -p 50 --paint $Paint
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i genus.xls   -n 2 -p 50 --paint $Paint
mv *.pdf percent.*.xls ./Community_barplot/

ln -s ../order ./sample
cp ../map-group.txt ./
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -a 2 -i phylum.xls -n 2 -p 50 -m map-group.txt --paint $Paint
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -a 2 -i class.xls  -n 2 -p 50 -m map-group.txt --paint $Paint
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -a 2 -i order.xls  -n 2 -p 50 -m map-group.txt --paint $Paint
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -a 2 -i family.xls -n 2 -p 50 -m map-group.txt --paint $Paint
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -a 2 -i genus.xls  -n 2 -p 50 -m map-group.txt --paint $Paint
mkdir Community_barplot_groups/
mv *.pdf sort.percent.*.xls ./Community_barplot_groups/

mkdir Community_average
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -a 1 -i phylum.xls -n 2 -p 50 -m map-group.txt --paint $Paint
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -a 1 -i class.xls  -n 2 -p 50 -m map-group.txt --paint $Paint
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -a 1 -i order.xls  -n 2 -p 50 -m map-group.txt --paint $Paint
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -a 1 -i family.xls -n 2 -p 50 -m map-group.txt --paint $Paint
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -a 1 -i genus.xls  -n 2 -p 50 -m map-group.txt --paint $Paint -z 10
mv *.pdf average.percent.*.xls ./Community_average

mkdir Community_test
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -r $Paired -o F -i phylum.xls -p 0     -m map-group.txt -c $Colour --combn $Combn
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -r $Paired -o F -i class.xls  -p 0.005 -m map-group.txt -c $Colour --combn $Combn
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -r $Paired -o F -i order.xls  -p 0.005 -m map-group.txt -c $Colour --combn $Combn
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -r $Paired -o F -i family.xls -p 0.005 -m map-group.txt -c $Colour --combn $Combn
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -r $Paired -o F -i genus.xls  -p 0.005 -m map-group.txt -c $Colour --combn $Combn
mv *.pdf *rank_sum*.xls Community_test
rm *.r percent*.xls Rplots.pdf

mkdir Heatmap_tax
plot-heatmap.pl -i Community_barplot_groups/sort.percent.genus.xls  -o heatmap.genus.pdf  -rtop 50 -ct 0 -slas 2 -rlc 0.7 -clc 0.7  -w $Heatmap_w -h $Heatmap_keyh -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt -color $Colour
plot-heatmap.pl -i Community_barplot_groups/sort.percent.family.xls -o heatmap.family.pdf -rtop 50 -ct 0 -slas 2 -rlc 0.7 -clc 0.7  -w $Heatmap_w -h $Heatmap_keyh -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt -color $Colour
plot-heatmap.pl -i Community_barplot_groups/sort.percent.order.xls  -o heatmap.order.pdf  -rtop 50 -ct 0 -slas 2 -rlc 0.7 -clc 0.7  -w $Heatmap_w -h $Heatmap_keyh -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt -color $Colour
plot-heatmap.pl -i Community_barplot_groups/sort.percent.class.xls  -o heatmap.class.pdf  -rtop 50 -ct 0 -slas 2 -rlc 0.7 -clc 0.7  -w $Heatmap_w -h $Heatmap_keyh -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt -color $Colour
plot-heatmap.pl -i Community_barplot_groups/sort.percent.phylum.xls -o heatmap.phylum.pdf -rtop 50 -ct 0 -slas 2 -rlc 0.7 -clc 0.7  -w $Heatmap_w -h $Heatmap_keyh -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt -color $Colour
#python /work/users/chaoliu/scripts/pheatmap.py -i Community_barplot_groups/sort.percent.genus.xls  -o pheatmap.genus.pdf -cd map-group.txt -clust_c F -cs 0 -rs 1 -scale none -cc $Colour
rm *.pdf.xls cmd.r dat.cor.dist.xls dat.cor.xls
mv heatmap*.pdf Heatmap_tax

#Community_boxplot
mkdir Community_boxplot
cd Community_boxplot
mkdir -p phylum genus family class order species
cd ..
if [ $(awk 'NR>1 {print $2}' map-group.txt | sort | uniq | wc -l) -gt 2 ]
then
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plot_community_boxplot.R -i Community_barplot/percent.species.xls -m map-group.txt -l 10 -p T --bb bx -c $Colour -o "Community_boxplot/species" --test "Community_test/species.Kruskal-Wallis_rank_sum.combn.xls"
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plot_community_boxplot.R -i Community_barplot/percent.genus.xls   -m map-group.txt -l 10 -p T --bb bx -c $Colour -o "Community_boxplot/genus"   --test "Community_test/genus.Kruskal-Wallis_rank_sum.combn.xls"  
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plot_community_boxplot.R -i Community_barplot/percent.family.xls  -m map-group.txt -l 10 -p T --bb bx -c $Colour -o "Community_boxplot/family"  --test "Community_test/family.Kruskal-Wallis_rank_sum.combn.xls" 
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plot_community_boxplot.R -i Community_barplot/percent.order.xls   -m map-group.txt -l 10 -p T --bb bx -c $Colour -o "Community_boxplot/order"   --test "Community_test/order.Kruskal-Wallis_rank_sum.combn.xls"  
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plot_community_boxplot.R -i Community_barplot/percent.class.xls   -m map-group.txt -l 10 -p T --bb bx -c $Colour -o "Community_boxplot/class"   --test "Community_test/class.Kruskal-Wallis_rank_sum.combn.xls"  
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plot_community_boxplot.R -i Community_barplot/percent.phylum.xls  -m map-group.txt -l 10 -p T --bb bx -c $Colour -o "Community_boxplot/phylum"  --test "Community_test/phylum.Kruskal-Wallis_rank_sum.combn.xls" 
fi

if [ $(awk 'NR>1 {print $2}' map-group.txt | sort | uniq | wc -l) == 2 ]
then
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plot_community_boxplot.R -i Community_barplot/percent.species.xls -m map-group.txt -l 10 -p T --bb bx -c $Colour -o "Community_boxplot/species" --test "Community_test/species.Wilcoxon_rank_sum_unpaired.combn.xls"
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plot_community_boxplot.R -i Community_barplot/percent.genus.xls   -m map-group.txt -l 10 -p T --bb bx -c $Colour -o "Community_boxplot/genus"   --test "Community_test/genus.Wilcoxon_rank_sum_unpaired.combn.xls"
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plot_community_boxplot.R -i Community_barplot/percent.family.xls  -m map-group.txt -l 10 -p T --bb bx -c $Colour -o "Community_boxplot/family"  --test "Community_test/family.Wilcoxon_rank_sum_unpaired.combn.xls"
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plot_community_boxplot.R -i Community_barplot/percent.order.xls   -m map-group.txt -l 10 -p T --bb bx -c $Colour -o "Community_boxplot/order"   --test "Community_test/order.Wilcoxon_rank_sum_unpaired.combn.xls"
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plot_community_boxplot.R -i Community_barplot/percent.class.xls   -m map-group.txt -l 10 -p T --bb bx -c $Colour -o "Community_boxplot/class"   --test "Community_test/class.Wilcoxon_rank_sum_unpaired.combn.xls"
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plot_community_boxplot.R -i Community_barplot/percent.phylum.xls  -m map-group.txt -l 10 -p T --bb bx -c $Colour -o "Community_boxplot/phylum"  --test "Community_test/phylum.Wilcoxon_rank_sum_unpaired.combn.xls"
fi

cd ../

############################################################# Beta ####################################################################
# PCoA
mkdir PCoA
rm -rf Pcoa
cd PCoA
rm *
ln -s ../map-group.txt ./
ln -s ../Beta_diversity/bray_curtis_dm.txt ./
ln -s ../Beta_diversity/weighted_unifrac_dm.txt ./
ln -s ../Beta_diversity/unweighted_unifrac_dm.txt ./
ln -s ../Beta_diversity/jaccard-binary_dm.txt ./
#noname
mkdir noname
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt        -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt        -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt        -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt   -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt   -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt   -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i jaccard-binary_dm.txt -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i jaccard-binary_dm.txt -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i jaccard-binary_dm.txt -md PCoA -pc 2-3 -map map-group.txt -col $Colour

mv *.pdf noname/

# 3D
mkdir 3D
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plotly_PCoA-3D.R -i bray_curtis_dm.txt        -c $Colour
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plotly_PCoA-3D.R -i weighted_unifrac_dm.txt   -c $Colour
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plotly_PCoA-3D.R -i unweighted_unifrac_dm.txt -c $Colour
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plotly_PCoA-3D.R -i jaccard-binary_dm.txt -c $Colour
mv *.pdf 3D

#name
mkdir name
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i bray_curtis_dm.txt        -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i bray_curtis_dm.txt        -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i bray_curtis_dm.txt        -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i weighted_unifrac_dm.txt   -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i weighted_unifrac_dm.txt   -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i weighted_unifrac_dm.txt   -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i unweighted_unifrac_dm.txt -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i unweighted_unifrac_dm.txt -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i unweighted_unifrac_dm.txt -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i jaccard-binary_dm.txt     -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i jaccard-binary_dm.txt     -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i jaccard-binary_dm.txt     -md PCoA -pc 2-3 -map map-group.txt -col $Colour
mv *.pdf name/

#ellipse
mkdir ellipse
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i bray_curtis_dm.txt        -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i bray_curtis_dm.txt        -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i bray_curtis_dm.txt        -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i weighted_unifrac_dm.txt   -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i weighted_unifrac_dm.txt   -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i weighted_unifrac_dm.txt   -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i unweighted_unifrac_dm.txt -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i unweighted_unifrac_dm.txt -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i unweighted_unifrac_dm.txt -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i jaccard-binary_dm.txt     -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i jaccard-binary_dm.txt     -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i jaccard-binary_dm.txt     -md PCoA -pc 2-3 -map map-group.txt -col $Colour
mv *.pdf ellipse/

#box
mkdir box
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i bray_curtis_dm.txt        -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i bray_curtis_dm.txt        -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i bray_curtis_dm.txt        -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i weighted_unifrac_dm.txt   -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i weighted_unifrac_dm.txt   -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i weighted_unifrac_dm.txt   -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i unweighted_unifrac_dm.txt -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i unweighted_unifrac_dm.txt -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i unweighted_unifrac_dm.txt -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i jaccard-binary_dm.txt     -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i jaccard-binary_dm.txt     -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i jaccard-binary_dm.txt     -md PCoA -pc 2-3 -map map-group.txt -col $Colour
mv *.pdf box/

# crossbar
mkdir crossbar
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt         -ct 2 -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt         -ct 2 -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt         -ct 2 -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt    -ct 2 -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt    -ct 2 -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt    -ct 2 -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt  -ct 2 -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt  -ct 2 -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt  -ct 2 -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i jaccard-binary_dm.txt      -ct 2 -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i jaccard-binary_dm.txt      -ct 2 -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i jaccard-binary_dm.txt      -ct 2 -md PCoA -pc 2-3 -map map-group.txt -col $Colour
mv *.pdf crossbar/

# central
mkdir central
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt         -ct 1 -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt         -ct 1 -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt         -ct 1 -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt    -ct 1 -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt    -ct 1 -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt    -ct 1 -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt  -ct 1 -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt  -ct 1 -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt  -ct 1 -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i jaccard-binary_dm.txt      -ct 1 -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i jaccard-binary_dm.txt      -ct 1 -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i jaccard-binary_dm.txt      -ct 1 -md PCoA -pc 2-3 -map map-group.txt -col $Colour
mv *.pdf central/

rm cmd.r 
cd ../
############################################################ CAP ##########################################################
mkdir CAP
cd CAP
ln -s ../map-group.txt ./
ln -s ../Beta_diversity/bray_curtis_dm.txt ./
ln -s ../Beta_diversity/weighted_unifrac_dm.txt ./
ln -s ../Beta_diversity/unweighted_unifrac_dm.txt ./
ln -s ../Beta_diversity/jaccard-binary_dm.txt ./
#noname
mkdir noname
python /work/users/chaoliu/scripts/gg_CAP.py -i bray_curtis_dm.txt        -md CAP -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_CAP.py -i weighted_unifrac_dm.txt   -md CAP -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_CAP.py -i unweighted_unifrac_dm.txt -md CAP -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_CAP.py -i jaccard-binary_dm.txt     -md CAP -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf noname/

#name
mkdir name
python /work/users/chaoliu/scripts/gg_CAP.py -lab T -i bray_curtis_dm.txt        -md CAP -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_CAP.py -lab T -i weighted_unifrac_dm.txt   -md CAP -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_CAP.py -lab T -i unweighted_unifrac_dm.txt -md CAP -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_CAP.py -lab T -i jaccard-binary_dm.txt     -md CAP -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf name/

#ellipse
mkdir ellipse
python /work/users/chaoliu/scripts/gg_CAP.py -e T -i bray_curtis_dm.txt        -md CAP -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_CAP.py -e T -i weighted_unifrac_dm.txt   -md CAP -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_CAP.py -e T -i unweighted_unifrac_dm.txt -md CAP -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_CAP.py -e T -i jaccard-binary_dm.txt     -md CAP -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf ellipse/

#box
mkdir box
python /work/users/chaoliu/scripts/gg_CAP.py -bx T -e T -i bray_curtis_dm.txt        -md CAP -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_CAP.py -bx T -e T -i weighted_unifrac_dm.txt   -md CAP -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_CAP.py -bx T -e T -i unweighted_unifrac_dm.txt -md CAP -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_CAP.py -bx T -e T -i jaccard-binary_dm.txt     -md CAP -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf box/

# crossbar
mkdir crossbar
python /work/users/chaoliu/scripts/gg_CAP.py -i bray_curtis_dm.txt         -ct 2 -md CAP -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_CAP.py -i weighted_unifrac_dm.txt    -ct 2 -md CAP -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_CAP.py -i unweighted_unifrac_dm.txt  -ct 2 -md CAP -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_CAP.py -i jaccard-binary_dm.txt      -ct 2 -md CAP -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf crossbar/

# central
mkdir central
python /work/users/chaoliu/scripts/gg_CAP.py -i bray_curtis_dm.txt         -ct 1 -md CAP -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_CAP.py -i weighted_unifrac_dm.txt    -ct 1 -md CAP -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_CAP.py -i unweighted_unifrac_dm.txt  -ct 1 -md CAP -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_CAP.py -i jaccard-binary_dm.txt      -ct 1 -md CAP -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf central/

rm cmd.r 
cd ../

############################################################ PCA ##########################################################
mkdir PCA
rm Pca -rf
cd PCA
ln -s ../map-group.txt ./
ln -s ../ASV_Taxa/rarefac.ASV_table.xls ./

#noname
mkdir noname
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i rarefac.ASV_table.xls -md PCA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i rarefac.ASV_table.xls -md PCA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i rarefac.ASV_table.xls -md PCA -pc 2-3 -map map-group.txt -col $Colour
mv PCA*.pdf noname/

/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/plotly_PCoA-3D.R -i rarefac.ASV_table.xls -d PCA

#name
mkdir name
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i rarefac.ASV_table.xls -md PCA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i rarefac.ASV_table.xls -md PCA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i rarefac.ASV_table.xls -md PCA -pc 2-3 -map map-group.txt -col $Colour
mv PCA*.pdf name/

#ellipse
mkdir ellipse
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i rarefac.ASV_table.xls -md PCA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i rarefac.ASV_table.xls -md PCA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i rarefac.ASV_table.xls -md PCA -pc 2-3 -map map-group.txt -col $Colour
mv PCA*.pdf ellipse/

#box
mkdir box
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i rarefac.ASV_table.xls -md PCA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i rarefac.ASV_table.xls -md PCA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i rarefac.ASV_table.xls -md PCA -pc 2-3 -map map-group.txt -col $Colour
mv PCA*.pdf box/

rm cmd.r 
cd ../

############################################################ NMDS ##########################################################
mkdir NMDS
rm Nmds -rf
cd NMDS
ln -s ../map-group.txt ./
ln -s ../Beta_diversity/bray_curtis_dm.txt ./
ln -s ../Beta_diversity/weighted_unifrac_dm.txt ./
ln -s ../Beta_diversity/unweighted_unifrac_dm.txt ./
ln -s ../Beta_diversity/jaccard-binary_dm.txt ./

#noname
mkdir noname
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt        -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt   -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i jaccard-binary_dm.txt -md NMDS -pc 1-2 -map map-group.txt -col $Colour

mv *.pdf noname/

#name
mkdir name
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i bray_curtis_dm.txt        -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i weighted_unifrac_dm.txt   -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i unweighted_unifrac_dm.txt -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i jaccard-binary_dm.txt -md NMDS -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf name/

#ellipse
mkdir ellipse
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i bray_curtis_dm.txt        -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i weighted_unifrac_dm.txt   -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i unweighted_unifrac_dm.txt -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i jaccard-binary_dm.txt -md NMDS -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf ellipse/

#box
mkdir box
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i bray_curtis_dm.txt        -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i weighted_unifrac_dm.txt   -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i unweighted_unifrac_dm.txt -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i jaccard-binary_dm.txt -md NMDS -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf box/

# crossbar
mkdir crossbar
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt         -ct 2 -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt    -ct 2 -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt  -ct 2 -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i jaccard-binary_dm.txt      -ct 2 -md NMDS -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf crossbar/

# central
mkdir central
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt         -ct 1 -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt    -ct 1 -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt  -ct 1 -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i jaccard-binary_dm.txt      -ct 1 -md NMDS -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf central/

rm cmd.r 
cd ../

##################################################### plsda  #############
mkdir PLS-DA
cd PLS-DA
ln -s ../map-group.txt ./
ln -s ../ASV_Taxa/rarefac.ASV_table.xls ./

perl /work/users/chaoliu/scripts/plsda.pl -i rarefac.ASV_table.xls -m map-group.txt -o ./ -g group -l F  -c $Colour
cd ..
#### treebar ###################################################
cd Hclust_bar/
plot-treebar.pl -otu ../ASV_Taxa/rarefac.ASV_table.xls -tax ../Community/phylum.xls -o treebar_ASV_table_phylum.pdf -h $Hclust_bar_h -lcex 1.3
cd ..

######################################### Lefse #################################
source /root/anaconda3/etc/profile.d/conda.sh
conda activate qiime1
export R_HOME_DIR=/root/anaconda3/envs/qiime1/lib/R/library
export R_LIBS=$R_HOME_DIR/lib64/R/library

mkdir Lefse
cd Lefse/
cp ../ASV_Taxa/rarefac.ASV_taxa_table.xls ./
sed -i 's/;s__.*//g' rarefac.ASV_taxa_table.xls
python2 /work/scripts/16s/tax_split.py -i rarefac.ASV_taxa_table.xls

sh /work/users/chaoliu/scripts/S_S/001.Community.OTU.ASV.lefse.sh  -c $Colour
cd ../
####################################################### Random_Forest
mkdir Random_Forest
cd Random_Forest/
ln -s ../ASV_Taxa/rarefac.ASV_table.xls ./
ln -s ../ASV_Taxa/rarefac.ASV_genus.xls ./
ln -s ../map-group.txt ./
ln -s ../order ./
biom convert -i rarefac.ASV_table.xls -o rarefac.ASV_table.biom --table-type "OTU table" --to-hdf5

random_forest4key_out_select_ASV.pl rarefac.ASV_table.biom rarefac.ASV_genus.xls  map-group.txt 0.001
sortSample4otu_table.pl OTU-extract.all.xls order > OTU-extract.all.sort.xls
cd ../

####################################################### Heatmap ####################################################################
cd Heatmap
rm ./*
ln -s ../ASV_Taxa/rarefac.ASV_genus.xls ./
ln -s ../map-group.txt ./
plot-heatmap.pl -i rarefac.ASV_genus.xls -o heatmap.ASV.top50.pdf -rtop 50 -ct 0 -slas 2 -rlc 0.7 -clc 0.7 -w $Heatmap_w -h $Heatmap_h -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt -color $Colour
ln -s ../Random_Forest/OTU-extract.all.sort.xls ./
bar_pie.pl -i OTU-extract.all.sort.xls -pie F -p 0.01
sed -i '/Others/d' ALL.new.OTU-extract.all.sort.xls
mv ALL.new.OTU-extract.all.sort.xls ALL.new.OTU-extract.all.sort.percent001.xls
bar_pie.pl -i OTU-extract.all.sort.xls -pie F -p 0.03
sed -i '/Others/d' ALL.new.OTU-extract.all.sort.xls
mv ALL.new.OTU-extract.all.sort.xls ALL.new.OTU-extract.all.sort.percent003.xls
plot-heatmap.pl -i OTU-extract.all.sort.xls                    -o heatmap.keyASV.pdf            -ct 0 -slas 2 -rlc 0.7 -clc 0.7 -w $Heatmap_w -h $Heatmap_keyh -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt -color $Colour
plot-heatmap.pl -i ALL.new.OTU-extract.all.sort.percent001.xls -o heatmap.keyASV.percent001.pdf -ct 0 -slas 2 -rlc 0.7 -clc 0.7 -w $Heatmap_w -h $Heatmap_h    -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt -color $Colour
plot-heatmap.pl -i ALL.new.OTU-extract.all.sort.percent003.xls -o heatmap.keyASV.percent003.pdf -ct 0 -slas 2 -rlc 0.7 -clc 0.7 -w $Heatmap_w -h $Heatmap_h    -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt -color $Colour
rm ALL*.xls dat* bar* cmd.r ./*ct*
cd ..

mkdir Heatmap_PRETTY
cd Heatmap_PRETTY
ln -s ../map-group.txt ./
ln -s ../ASV_Taxa/rarefac.ASV_genus.xls ./
python /work/users/chaoliu/scripts/pheatmap.py -i rarefac.ASV_genus.xls -cd map-group.txt -o pheatmap.ASV.top50.pdf -rtop 50 -clust_c F -cs 0 -rs 1 -scale none -cc $Colour
ln -s ../Random_Forest/OTU-extract.all.sort.xls ./
bar_pie.pl -i OTU-extract.all.sort.xls -pie F -p 0.01
sed -i '/Others/d' ALL.new.OTU-extract.all.sort.xls
mv ALL.new.OTU-extract.all.sort.xls ALL.new.OTU-extract.all.sort.percent001.xls
bar_pie.pl -i OTU-extract.all.sort.xls -pie F -p 0.03
sed -i '/Others/d' ALL.new.OTU-extract.all.sort.xls
mv ALL.new.OTU-extract.all.sort.xls ALL.new.OTU-extract.all.sort.percent003.xls
python /work/users/chaoliu/scripts/pheatmap.py -i OTU-extract.all.sort.xls                    -o pheatmap.keyASV.pdf            -cd map-group.txt -clust_c F -cs 0 -rs 1 -scale none -cc $Colour
python /work/users/chaoliu/scripts/pheatmap.py -i ALL.new.OTU-extract.all.sort.percent001.xls -o pheatmap.keyASV.percent001.pdf -cd map-group.txt -clust_c F -cs 0 -rs 1 -scale none -cc $Colour
python /work/users/chaoliu/scripts/pheatmap.py -i ALL.new.OTU-extract.all.sort.percent003.xls -o pheatmap.keyASV.percent003.pdf -cd map-group.txt -clust_c F -cs 0 -rs 1 -scale none -cc $Colour
rm ALL*.xls dat* bar* cmd.r ./*ct*
cd ../

###################################################### Tax4Fun2 ###############################################################################################
if [ -d Tax4Fun2 ];then
cd Tax4Fun2/
sh /work/users/chaoliu/scripts/S_S/Tax4Fun2.lefse.sh -c $Colour
cd ..
fi
##################################################### Picrust2 ################################################################################################################
if [ -d Picrust2 ];then
cd Picrust2/
sh /work/users/chaoliu/scripts/S_S/Picrust2.lefse.sh -c $Colour
cd ..
fi

######################################################################### change_color ####################################################################
#if [ "$Colour" != "none" ];then
#sh /work/users/chaoliu/scripts/S_S/002.functional_prediction.OTU_ASV.change_color.sh -c $Colour
#fi
##################################################################################################################
conda deactivate

export R_HOME_DIR=/usr/local/R-3.6.0
export R_LIBS=$R_HOME_DIR/lib64/R/library

if [ -d Tax4Fun2 ];then
cd  Tax4Fun2/
mkdir Tax4Fun_pathway_summary
mv Tax4Fun.pathways*_out.tsv Tax4Fun.KO_out.tsv Tax4Fun_pathway_summary
cd Tax4Fun_pathway_summary
rename 's/Tax4Fun.pathways/Tax4Fun_pathways/g' *.tsv
rename 's/Tax4Fun.KO/Tax4Fun_KO/g' *.tsv
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i Tax4Fun_pathwaysLevel1_out.tsv -n 2 -p 50 -c $Colour --combn $Combn
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i Tax4Fun_pathwaysLevel2_out.tsv -n 2 -p 50 -c $Colour --combn $Combn
#/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i Tax4Fun_pathwaysLevel3_out.tsv -n 2 -p 50 -c $Colour --combn $Combn
#/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i Tax4Fun_KO_out.tsv             -n 2 -p 50 -c $Colour --combn $Combn

/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -o F -i Tax4Fun_pathwaysLevel1_out.tsv -n 0 -p 0.005 -m ../../map-group.txt -c $Colour --combn $Combn
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -o F -i Tax4Fun_pathwaysLevel2_out.tsv -n 0 -p 0.005 -m ../../map-group.txt -c $Colour --combn $Combn
#/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -o F -i Tax4Fun_pathwaysLevel3_out.tsv -n 0 -p 0.005 -m ../../map-group.txt -c $Colour --combn $Combn
#/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -o F -i Tax4Fun_KO_out.tsv             -n 0 -p 0.005 -m ../../map-group.txt -c $Colour --combn $Combn
rm Rplots.pdf
cd ../../
fi
##################################################################################################################################################################################
if [ -d Picrust2 ];then
cd Picrust2
cd KEGG_pathways_out
mkdir KEGG_pathway_summary
mv KEGG_pathwaysL* KEGG_pathway_summary
cd KEGG_pathway_summary
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i KEGG_pathwaysLevel1_out.tsv -n 2 -p 50 -c $Colour --combn $Combn
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i KEGG_pathwaysLevel2_out.tsv -n 2 -p 50 -c $Colour --combn $Combn
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i KEGG_pathwaysLevel3_out.tsv -n 2 -p 50 -c $Colour --combn $Combn

/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -o F -i KEGG_pathwaysLevel1_out.tsv -n 0 -p 0.005  -m ../../map-group.txt -c $Colour --combn $Combn
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -o F -i KEGG_pathwaysLevel2_out.tsv -n 0 -p 0.005  -m ../../map-group.txt -c $Colour --combn $Combn
/usr/local/R-3.6.0/bin/Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -o F -i KEGG_pathwaysLevel3_out.tsv -n 0 -p 0.005  -m ../../map-group.txt -c $Colour --combn $Combn
rm Rplots.pdf
cd ../../../
fi

####################################################################################################################################################################################
rm NMDS/*.txt PCoA/*.txt PCA/*.txt  Community/*ALL* Community/*sh Community/sample*
rm -r Estimators
rm Rarefactions/*.rarefaction Shannon_rarefac/*.r_shannon
rm Heatmap/ALL.new.* Heatmap/map-group.txt 
rm Heatmap/OTU-extract.all.sort.xls Heatmap/dat.cor.* Heatmap/cmd.r Heatmap/bar.ALL.OTU-extract.all.sort.xls.pdf Heatmap/map-group.txt Heatmap/OTU-extract.all.sort.xls Heatmap/rarefac.otu_genus.xls
rm Heatmap/-ct Heatmap/-ct.xls 
rm Hclust_bar/cmd.r OTU_Taxa/otu_seqids.txt Venn/rarefac.otu_genus.xls Community/percent.*.xls

#find ./ -name "*.biom" |xargs rm -rf
find ./ -type d -name "tax_summary_a"| xargs rm -rf
find ./ -type d -name "tax_summary_r"| xargs rm -rf
find . -name "Rplots.pdf" | xargs rm -rf
find . -name "*cmd.r"| xargs rm -rf
#find . -name "otu_seqids.txt"| xargs rm -rf
# 删除失效软连接
for file in `find . -type l`
do
    if [ ! -e $file ]
    then
        echo "rm $file"
        rm -f $file
    fi
done
find . -name "*" -type f -size 0c | xargs -n 1 rm -f

```
## Part5.liner_regression
```R
library(ggplot2)
library(MASS)
library(splines)
require("tidyverse")
library(fs)
library(optparse)
library(data.table)
library(openxlsx, lib.loc = "/usr/local/R-4.0.5/lib64/R/library")
rm(list = ls())
########################### 01.set parameters ###########################
if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"),   type="character",  default="genus.xls", help="input"),
    make_option(c("-s", "--select"),  type="character",  default="none",      help="select.list"),
    make_option(c("-e", "--env"),     type="character",  default="env.txt",   help="env"),
    make_option(c("--log"),           type="character",  default="2",         help="0 for no transfer"),
    make_option(c("--map"),           type="character",  default="map.txt",   help="group"),
    make_option(c("-v", "--value"),   type="double",     default=0.75,        help="core_microbiome"),
    make_option(c("-p", "--pvalue"),  type="double",     default=0.05,        help="pvalue"),
    make_option(c("--bra"),           type="character",  default="5",         help="mode")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
}

########################## 02.rean-in ############################
genus <- fread(opts$input)
if(opts$log=="1"){
  #install.packages("zCompositions", version = "1.3.4")
  library(zCompositions)
  #installed.packages("devtools")
  #installed.packages("rJava")
  library(devtools)
  library(rJava)
  #devtools::install_github("OHDSI/Achilles")
  library(Achilles)
  library(DatabaseConnector)
  #devtools::install_github("ggloor/CoDaSeq/CoDaSeq")
  library(CoDaSeq)
  genus.CZM <- cmultRepl(genus[,-1],method = "CZM")
  genus2 <- cbind(genus[,1],codaSeq.clr(genus.CZM)) %>% as.data.frame %>% rename_at(1,~"genus") 
  #genus2[genus==0]<-0
  genus <- genus2
}else if(opts$log=="2"){
  genus0<-genus[,-1]+1
  genus <- cbind(genus[,1],codaSeq.clr(genus0)) %>% as.data.frame %>% rename_at(1,~"genus") 
  save_data <- function(df,wd, name) {
    addWorksheet(wb, name)
    writeDataTable(wb, name, df)
  }
  
  wb <- openxlsx::createWorkbook() 
  save_data(data.frame(cbind(genus[,1],genus0)),wd,"genus1")
  save_data(data.frame(genus) ,wd,"genus2")
  saveWorkbook(wb,"genus.plus.xls", overwrite = T) 
}

if(opts$select!="none"){
  ss <- read_tsv("select.list",col_names = F) %>% rename_at(1,~"genus")
  genus <- genus %>% inner_join(ss) 
}

DATA <- genus %>% #mutate_if(is.numeric,~map_dbl(.x,function(y) y/ sum(.x))) %>%
  rename_at(1,~"ID") %>% gather(SampleID,value,-ID) %>% mutate_if(is.character,~fct_inorder(.x)) %>%
  spread(ID,value) %>% mutate_if(is.factor,~as.character(.x))

Env <- read_tsv(opts$env) %>% rename_at(1,~"SampleID")

DATA <- DATA %>% inner_join(Env[,1]) %>% as_tibble
Env <- Env %>% inner_join(DATA[,1],.)
colnames(Env) <- str_replace(colnames(Env),"/","_")
if(opts$map != "none"){
  map <- read_tsv(opts$map) %>% rename_all(~c("SampleID","gender"))
}

Minus <- function(x,n){
  d1 <- 10^(-n)
  ifelse(x >= d1,round(x,4),paste0("< ",d1))
}

mytheme <- function(){
  theme_bw() + 
    theme(plot.title=element_text(size=rel(1),hjust=0.5),
          plot.margin = unit(c(3, 3, 3, 3), "lines"),
          axis.title=element_text(size=rel(1)),
          legend.title=element_blank(),
          legend.text = element_text(size = 6),
          legend.key.size = unit(.4,'cm')
    )
}

## 1.Analysis ä¸€ä¸ªè�Œä¸‰ä¸ªå¹´é¾„æ®µåœ¨å�Œä¸€å¼ å›¾ä¸Šçš„çº¿æ€§å›¾ ###########
if(opts$bra == "1"){
  Analysis <- function(huge){ # x <- colnames(Env)[2];y <- colnames(DATA)[2];huge <- bind_cols(DATA[,y],Env[,x])
    NM <- colnames(huge)
    huge <- huge %>% rename_all(~c("PA","DA")) %>% mutate(group=map_chr(DA,~if_else(.x<=15,"A0_15",if_else(.x<=70,"A16_70","A71_120")))) %>% 
      mutate(group=fct_inorder(group)) %>% group_by(group) %>% nest() %>% mutate(Data2=pmap(list(group,data),function(z,u){ #z<-huge$group[[1]];u<-huge$data[[1]]
        Data_lm <- lm(PA~DA,data=u )
        Data_sum <- summary(Data_lm)
        
        r1 <- paste("=",as.character(round(Data_sum$r.squared,4)))
        p1 <- paste("; p=",as.character(Minus(Data_sum$coefficients[2,4],4)))
        s1 <- str_c("; y = ",as.character(round(Data_sum$coefficients[2,1],4)),"x")
        t1 <- ifelse(round(Data_sum$coefficients[1,1],4) < 0,
                     as.character(round(Data_sum$coefficients[1,1],4)),
                     str_c("+",as.character(round(Data_sum$coefficients[1,1],4))))
        new_T <- seq(min(u$DA), max(u$DA), 0.01)
        pred_wt <- data.frame(predict(Data_lm, newdata = data.frame(DA = new_T),
                                      interval = "confidence"), new_T = new_T)
        p <- ggplot() +
          geom_point(data = u, mapping = aes(x = DA, y = PA),color="#B7C59D",size=3) +
          geom_line(data = pred_wt, mapping = aes(x = new_T, y = fit), 
                    color = "blue", size = 1) +
          geom_ribbon(data = pred_wt, mapping = aes(x = new_T, ymin = lwr, ymax = upr), 
                      fill = "grey", alpha = 0.5) +
          labs(title = bquote(R^"2" ~ .(r1) ~ .(p1) ~ .(s1) ~ .(t1))) + #
          xlab(NM[2]) + ylab(NM[1]) +
          mytheme()
        return(list(tibble(R2=as.character(round(Data_sum$r.squared,4)),
                           pvalue=as.character(round(Data_sum$coefficients[2,4],4)),
                           slope=as.character(round(Data_sum$coefficients[2,1],4)),
                           Intercept=t1),p))
      }))
    
    n1 <- as.character(huge$group[[1]])
    d1 <- huge$Data2[[1]][[1]]
    q1 <- huge$Data2[[1]][[2]]
    
    for (i in 2:nrow(huge)){# i <- 2
      as.character(huge$group[[i]])
      d1 <- bind_rows(d1,huge$Data2[[i]][[1]])
      q1 <- q1 %>% aplot::insert_right(huge$Data2[[i]][[2]])
    }
    
    ggsave(str_c(str_c(NM,collapse = "--"),".output.all.pdf"),q1,width = 15,height = 5)   
    write_tsv(d1 %>% add_column(.before = 1,Groups = huge$group),str_c(str_c(NM,collapse = "--"),".output.all.tsv"))    
    
  }
  for (x in colnames(Env)[-1]){# x <- colnames(Env)[2]
    for (y in colnames(DATA)[-1]){# y <- colnames(DATA)[2]
      Data <- bind_cols(DATA[,y],Env[,x])
      if(Data[,-2] %>% as.matrix %>% apply(.,2,sum)!=0){
        Analysis(Data)
        NM <- colnames(Data)
        Data <- Data %>% rename_all(~c("PA","DA"))
        Data_lm <- lm(PA~DA,data=Data)
        Data_sum <- summary(Data_lm)
      }
    }
  }
}

## 2.Analysis2ä¸€ä¸ªè�Œotuä¸°åº¦å’Œå¹´é¾„çš„æ•£ç‚¹å›¾å’Œæ¯�ä¸€å²�å¹´é¾„ä¸‹å¹³å�‡ä¸°åº¦çš„æŠ˜çº¿ï¼ˆæ¯�ä¸€å²�ä¸‹å¦‚æžœå°‘äºŽ3ä¸ªæ ·æœ¬å°±ä¸�çº³å…¥ï¼‰ ###########
## ï¼ˆ0-3å²�å››èˆ�äº”å…¥ï¼Œ0.5ä»¥ä¸‹ç®—å�š0å²�ï¼Œå�ªç”¨å�‡å€¼è¿žçº¿ï¼Œåˆ†ç”·å¥³ï¼Œä¸€ä¸ªå›¾ä¿©çº¿ï¼‰
## æ¯”å¦‚3å²�çš„äººæœ‰6ä¸ªï¼Œå…¶ä¸­2å¥³4ç”·ï¼Œä¸¤ä¸ªéƒ½å±•ç¤º
if(opts$bra == "2"){
  Analysis2 <- function(huge,ym,bee){ # x <- colnames(Env)[2];y <- colnames(DATA)[2];huge <- bind_cols(DATA[,y],Env[,x]);ym <- map ;bee <- DATA[,1]
    NM <- colnames(huge)
    huge <- bind_cols(huge,bee) %>% rename_all(~c("PA","DA","SampleID")) %>% 
      left_join(map)
    jingtian <- huge %>% mutate(group=DA) %>% group_by(group) %>% nest() %>% mutate(Data2=map_chr(data,function(u){#u<-jingtian$data[[1]]
      if_else(nrow(u)<3,"nrow<3","nrow>=3") })) %>% #filter(Data2=="nrow>=3") %>% 
      dplyr::select(-Data2) %>% mutate(M=map_chr(data,function(gg){#gg<-jingtian$data[[3]]
        gg %>% filter(gender=="M") %>% select(PA) %>% apply(.,2,mean) })) %>% mutate(F=map_chr(data,function(mm){
          mm %>% filter(gender=="F") %>% select(PA) %>% apply(.,2,mean)
        })) %>% select(-data) %>% mutate(M=as.numeric(M),F=as.numeric(F)) %>%
      gather(key="gender",value="otu",-group) %>% mutate(gender=fct_inorder(gender)) %>% 
      subset(otu != "NaN")
    
    p <- ggplot() +
      geom_point(data = huge, mapping = aes(x = DA, y = PA ,color=gender,group=gender),alpha=0.3,size=1) +
      geom_line(data = jingtian, mapping = aes(x = group, y = otu, color = gender), size = 0.8) +
      xlab(NM[2]) + ylab(NM[1]) +
      scale_x_continuous(breaks = seq(0,100,by=5)) +
      mytheme()
    
    ggsave(str_c(str_c(NM,collapse = "--"),".pdf"),p,width = 12,height = 5)
  }
  for (x in colnames(Env)[-1]){# x <- colnames(Env)[2]
    for (y in colnames(DATA)[-1]){# y <- colnames(DATA)[2]
      Data <- bind_cols(DATA[,y],Env[,x]) 
      if(Data[,-2] %>% as.matrix %>% apply(.,2,sum)!=0){
        Analysis2(Data,map,DATA[,1])
      }
    }
  }
}

## 3.Analysis3å±žæ°´å¹³å�šä¸€ä¸‹ä¸‰ä¸ªå¹´é¾„æ®µè�Œå±žä¸Žå¹´é¾„çš„çº¿æ€§ç›¸å…³å›¾ï¼Œç»˜åˆ¶çº¿æ€§ç›¸å…³å›¾çš„å±žè¦�æ±‚ï¼š ###################
##ï¼ˆ1ï¼‰åœ¨ä»»æ„�ä¸€ç»„ä¸­æœ‰è‡³å°‘75%çš„æ ·æœ¬ä¸­èƒ½æ£€æµ‹åˆ°ï¼›
##ï¼ˆ2ï¼‰è‡³å°‘åœ¨ä¸€ä¸ªç»„ä¸­ä¸Žå¹´é¾„æœ‰æ˜¾è‘—ç›¸å…³æ€§ï¼›å…­æ�¡çº¿ éƒ½ä¸�æ˜¾è‘— p>0.05 å°±ä¸�è¦�å›¾
##ï¼ˆ3ï¼‰å‡ºå›¾æ—¶è¦�æ±‚æ ·æœ¬ç‚¹æŒ‰ç…§ç”·å¥³æ€§åˆ«åˆ†åˆ«ç”¨ä¸�å�Œé¢œè‰²å±•ç¤ºï¼Œä¸”ç”·å¥³åˆ†åˆ«å�šä¸¤æ�¡æ‹Ÿå�ˆç›´çº¿ï¼›
##ï¼ˆ4ï¼‰æŠŠæ‰€æœ‰ç”»å›¾çš„å±žåœ¨ä¸‰ç»„ä¸¤ä¸ªæ€§åˆ«ä¸­çš„çº¿æ€§ç›¸å…³æ•°å€¼ï¼ˆR2ã€�på€¼ã€�æ–œçŽ‡ï¼‰éƒ½æ•´ç�†åˆ°ä¸€ä¸ªæ•°æ�®è¡¨ä¸­ï¼›
if(opts$bra == "3"){
  Analysis3 <- function(huge,ym,hee){ # x <- colnames(Env)[2];y <- colnames(DATA)[3];huge <- bind_cols(DATA[,y],Env[,x]);ym <- map ;hee <- DATA[,1]
    NM <- colnames(huge)
    huge2 <- bind_cols(huge,hee) %>% rename_all(~c("PA","DA","SampleID")) %>% mutate(group=map_chr(DA,~if_else(.x<=15,"A0_15",if_else(.x<=70,"A16_70","A71_100")))) %>% 
      left_join(ym) %>% mutate(group=fct_inorder(group),gender=fct_inorder(gender)) %>% group_by(group) %>% nest() %>% 
      mutate(value=map_chr(data,function(ww){ #ww<-huge$data[[1]]
        colSums(ww != 0)[1]/nrow(ww) }))
    
    if(max(huge2$value)>=opts$value){
      huge3 <- huge2 %>% dplyr::select(-value) %>% mutate(pvalue=map_dbl(data,function(kk){ #kk<-huge2$data[[1]] 
        gg_lm <- lm(PA~DA,data=kk %>% filter(gender=="M"))
        gg_sum <- summary(gg_lm)
        mm_lm <- lm(PA~DA,data=kk %>% filter(gender=="F"))
        mm_sum <- summary(mm_lm)
        min(gg_sum$coefficients[2,4],mm_sum$coefficients[2,4])
      }))
      
      if(min(huge3$pvalue)<=opts$pvalue){
        huge4 <- huge3 %>% mutate(Data2=pmap(list(group,data),function(jj,kk){ #jj<-huge3$group[[1]]; kk<-huge3$data[[1]]
          group<-as.character(jj)
          gg_lm <- lm(PA~DA,data=kk %>% filter(gender=="M"))
          gg_sum <- summary(gg_lm)
          mm_lm <- lm(PA~DA,data=kk %>% filter(gender=="F"))
          mm_sum <- summary(mm_lm)
          min(gg_sum$coefficients[2,4],mm_sum$coefficients[2,4])
          
          
          r1 <- round(gg_sum$r.squared,4)
          p1 <- Minus(gg_sum$coefficients[2,4],4) %>% as.character
          s1 <- round(gg_sum$coefficients[2,1],4)
          t1 <- ifelse(round(gg_sum$coefficients[1,1],4) < 0,
                       as.character(round(gg_sum$coefficients[1,1],4)),
                       str_c("+",as.character(round(gg_sum$coefficients[1,1],4))))
          r2 <- round(mm_sum$r.squared,4)
          p2 <- Minus(mm_sum$coefficients[2,4],4) %>% as.character
          s2 <- round(mm_sum$coefficients[2,1],4)
          t2 <- ifelse(round(mm_sum$coefficients[1,1],4) < 0,
                       as.character(round(mm_sum$coefficients[1,1],4)),
                       str_c("+",as.character(round(mm_sum$coefficients[1,1],4))))
          
          new_T <- seq(min(kk$DA), max(kk$DA), 0.0001)
          gg_pred_wt <- data.frame(predict(gg_lm, newdata = data.frame(DA = new_T),
                                           interval = "confidence"), new_T = new_T)
          mm_pred_wt <- data.frame(predict(mm_lm, newdata = data.frame(DA = new_T),
                                           interval = "confidence"), new_T = new_T)
          p <- ggplot() +
            geom_point(data = kk, mapping = aes(x = DA, y = PA ,color=gender),alpha=0.6,size=2) + #alpha=0.6,
            geom_line(data = gg_pred_wt, mapping = aes(x = new_T, y = fit), color = "#f8766d" ,size = 1) +
            geom_line(data = mm_pred_wt, mapping = aes(x = new_T, y = fit), color = "#00bfc4" ,size = 1) +
            #geom_ribbon(data = pred_wt, mapping = aes(x = new_T, ymin = lwr, ymax = upr), 
            #            fill = "grey", alpha = 0.5) +
            #labs(title = bquote(R^"2" ~ .(r1) ~ .(p1) ~ .(s1) ~ .(t1))) + #
            xlab(NM[2]) + ylab(NM[1]) +
            mytheme()
          return(list(tibble(genus=c(NM[1],NM[1]),group=c(group,group),R2=c(r1,r2),pvalue=c(p1,p2),slope=c(s1,s2),Intercept=c(t1,t2),gender=c("M","F")),p))
        }))
        
        d1 <- huge4$Data2[[1]][[1]]
        q1 <- huge4$Data2[[1]][[2]]
        
        for (i in 2:nrow(huge4)){# i <- 2
          d1 <- bind_rows(d1,huge4$Data2[[i]][[1]])
          q1 <- q1 %>% aplot::insert_right(huge4$Data2[[i]][[2]])
        }
        ggsave(str_c(str_c(NM,collapse = "--"),".output.all.pdf"),q1,width = 15,height = 5)   
        write_tsv(d1 %>% add_column(.before = 1,Groups = huge$group),str_c(str_c(NM,collapse = "--"),".output.all.tsv"))    
        
      }
      
    }  
    
  }
  for (x in colnames(Env)[-1]){# x <- colnames(Env)[2]
    for (y in colnames(DATA)[-1]){# y <- colnames(DATA)[3]
      Data <- bind_cols(DATA[,y],Env[,x]) 
      if(Data[,-2] %>% as.matrix %>% apply(.,2,sum)!=0){
        Analysis3(Data,map,DATA[,1])
      }
    }
  }
}

## 4.ASV #############################
if(opts$bra == "4"){
  Analysis4 <- function(huge,ym,hee){ # x <- colnames(Env)[2];y <- colnames(DATA)[3];huge <- bind_cols(DATA[,y],Env[,x]);ym <- map ;hee <- DATA[,1]
    NM <- colnames(huge)
    huge2 <- bind_cols(huge,hee) %>% rename_all(~c("PA","DA","SampleID")) %>% left_join(ym) %>% mutate(gender=fct_inorder(gender)) #%>% group_by(gender) %>% nest() 
    gg <- huge2 %>% filter(gender=="M") ; mm <- huge2 %>% filter(gender=="F")
    Data_lm1 <- lm(PA~poly(DA,2),data=gg)
    Data_sum1 <- summary(Data_lm1)
    
    new_T1 <- seq(min(gg$DA), max(gg$DA), 0.01)
    pred_wt1 <- data.frame(predict(Data_lm1, newdata = data.frame(DA = new_T1),
                                   interval = "confidence"), new_T = new_T1)
    
    r.squared1 <- round(Data_sum1$r.squared,4)
    adj.r.squared1 <- round(Data_sum1$adj.r.squared,4)
    p_x_1 <- Minus(Data_sum1$coefficients[2,4],4)
    p_x2_1 <- Minus(Data_sum1$coefficients[3,4],4)
    c1 <- ifelse(round(Data_sum1$coefficients[1,1],4) < 0,
                 as.character(round(Data_sum1$coefficients[1,1],4)),
                 str_c("+",as.character(round(Data_sum1$coefficients[1,1],4))))
    b1 <- ifelse(round(Data_sum1$coefficients[2,1],4) < 0,
                 as.character(round(Data_sum1$coefficients[2,1],4)),
                 str_c("+",as.character(round(Data_sum1$coefficients[2,1],4))))
    a1 <- ifelse(round(Data_sum1$coefficients[3,1],4) < 0,
                 as.character(round(Data_sum1$coefficients[3,1],4)),
                 str_c("+",as.character(round(Data_sum1$coefficients[3,1],4))))
    
    Data_lm2 <- lm(PA~poly(DA,2),data=mm)
    Data_sum2 <- summary(Data_lm2)
    
    new_T2 <- seq(min(mm$DA), max(mm$DA), 0.01)
    pred_wt2 <- data.frame(predict(Data_lm2, newdata = data.frame(DA = new_T2),
                                   interval = "confidence"), new_T = new_T2)
    
    r.squared2 <- round(Data_sum2$r.squared,4)
    adj.r.squared2 <- round(Data_sum2$adj.r.squared,4)
    p_x_2 <- Minus(Data_sum2$coefficients[2,4],4)
    p_x2_2 <- Minus(Data_sum2$coefficients[3,4],4)
    c2 <- ifelse(round(Data_sum2$coefficients[1,1],4) < 0,
                 as.character(round(Data_sum2$coefficients[1,1],4)),
                 str_c("+",as.character(round(Data_sum2$coefficients[1,1],4))))
    b2 <- ifelse(round(Data_sum2$coefficients[2,1],4) < 0,
                 as.character(round(Data_sum2$coefficients[2,1],4)),
                 str_c("+",as.character(round(Data_sum2$coefficients[2,1],4))))
    a2 <- ifelse(round(Data_sum2$coefficients[3,1],4) < 0,
                 as.character(round(Data_sum2$coefficients[3,1],4)),
                 str_c("+",as.character(round(Data_sum2$coefficients[3,1],4))))
    
    p <- ggplot() +
      geom_point(data = huge2, mapping = aes(x = DA, y = PA,color=gender),alpha=0.3,size=1) + #alpha=0.3,
      geom_line(data = pred_wt1, mapping = aes(x = new_T, y = fit), color = "#f8766d", size = 0.8) +
      geom_line(data = pred_wt2, mapping = aes(x = new_T, y = fit), color = "#00bfc4", size = 0.8) +
      xlab(NM[2]) + ylab(NM[1]) +
      mytheme()
    
    gold <- tibble(genus=c(NM[1],NM[1]),a=c(a1,a2),b=c(b1,b2),Intercept=c(c1,c2),
                   R2=c(r.squared1,r.squared2),adj_R2=c(adj.r.squared1,adj.r.squared2),
                   p1=c(p_x_1,p_x_2),p2=c(p_x2_1,p_x2_2),
                   gender=c("M","F"))
    ggsave(str_c(str_c(NM,collapse = "--"),".output.all.pdf"),p,width = 6,height = 5)   
    write_tsv(gold,str_c(str_c(NM,collapse = "--"),".output.all.tsv"))    
    
  }
  for (x in colnames(Env)[-1]){# x <- colnames(Env)[2]
    for (y in colnames(DATA)[-1]){# y <- colnames(DATA)[3]
      Data <- bind_cols(DATA[,y],Env[,x]) 
      if(Data[,-2] %>% as.matrix %>% apply(.,2,sum)!=0){
        Analysis4(Data,map,DATA[,1])
      }
    }
  }
}


## 5.Analysis ###################
if(opts$bra == "5"){
  Dooo <- fread(opts$input)
  Doo1 <- Dooo %>% #mutate_if(is.numeric,~map_dbl(.x,function(y) y/ sum(.x))) %>%
    rename_at(1,~"ID") %>% gather(SampleID,value,-ID) %>% mutate_if(is.character,~fct_inorder(.x)) %>%
    spread(ID,value) %>% mutate_if(is.factor,~as.character(.x)) %>% 
    inner_join(Env[,1]) %>% as_tibble
  
  Analysis3 <- function(huge,ym,hee){ # x <- colnames(Env)[2];y <- colnames(DATA)[3];huge <- bind_cols(DATA[,y],Env[,x],Doo1[,y] %>% rename_at(1,~"QAZ"));ym <- map ;hee <- DATA[,1]
    NM <- colnames(huge)[-3]
    huge2 <- bind_cols(huge,hee) %>% rename_all(~c("PA","DA","QAZ","SampleID")) %>% mutate(group=map_chr(DA,~if_else(.x<=15,"A0_15",if_else(.x<=70,"A16_70","A71_100")))) %>% 
      left_join(ym) %>% mutate(group=fct_inorder(group),gender=fct_inorder(gender)) %>% group_by(group) %>% nest() %>% 
      mutate(value=map_chr(data,function(ww){ #ww<-huge2$data[[1]]
        colSums(ww != 0)[3]/nrow(ww) }))
    
    if(max(huge2$value)>=opts$value){
      huge3 <- huge2 %>% dplyr::select(-value) %>% mutate(pvalue=map_dbl(data,function(kk){ #kk<-huge2$data[[1]] 
        gg_lm <- lm(PA~DA,data=kk %>% filter(gender=="M"))
        gg_sum <- summary(gg_lm)
        mm_lm <- lm(PA~DA,data=kk %>% filter(gender=="F"))
        mm_sum <- summary(mm_lm)
        min(gg_sum$coefficients[2,4],mm_sum$coefficients[2,4])
      }))
      
      if(min(huge3$pvalue)<=opts$pvalue){
        huge4 <- huge3 %>% mutate(Data2=pmap(list(group,data),function(jj,kk){ #jj<-huge3$group[[1]]; kk<-huge3$data[[1]]
          group<-as.character(jj)
          gg_lm <- lm(PA~DA,data=kk %>% filter(gender=="M"))
          gg_sum <- summary(gg_lm)
          mm_lm <- lm(PA~DA,data=kk %>% filter(gender=="F"))
          mm_sum <- summary(mm_lm)
          min(gg_sum$coefficients[2,4],mm_sum$coefficients[2,4])
          
          
          r1 <- round(gg_sum$r.squared,4)
          p1 <- Minus(gg_sum$coefficients[2,4],4) %>% as.character
          s1 <- round(gg_sum$coefficients[2,1],4)
          t1 <- ifelse(round(gg_sum$coefficients[1,1],4) < 0,
                       as.character(round(gg_sum$coefficients[1,1],4)),
                       str_c("+",as.character(round(gg_sum$coefficients[1,1],4))))
          r2 <- round(mm_sum$r.squared,4)
          p2 <- Minus(mm_sum$coefficients[2,4],4) %>% as.character
          s2 <- round(mm_sum$coefficients[2,1],4)
          t2 <- ifelse(round(mm_sum$coefficients[1,1],4) < 0,
                       as.character(round(mm_sum$coefficients[1,1],4)),
                       str_c("+",as.character(round(mm_sum$coefficients[1,1],4))))
          
          new_T <- seq(min(kk$DA), max(kk$DA), 0.0001)
          gg_pred_wt <- data.frame(predict(gg_lm, newdata = data.frame(DA = new_T),
                                           interval = "confidence"), new_T = new_T)
          mm_pred_wt <- data.frame(predict(mm_lm, newdata = data.frame(DA = new_T),
                                           interval = "confidence"), new_T = new_T)
          p <- ggplot() +
            geom_point(data = kk, mapping = aes(x = DA, y = PA ,color=gender),alpha=0.6,size=2) + #alpha=0.6,
            geom_line(data = gg_pred_wt, mapping = aes(x = new_T, y = fit), color = "#f8766d" ,size = 1) +
            geom_line(data = mm_pred_wt, mapping = aes(x = new_T, y = fit), color = "#00bfc4" ,size = 1) +
            xlab(NM[2]) + ylab(NM[1]) +
            mytheme()
          return(list(tibble(genus=c(NM[1],NM[1]),group=c(group,group),R2=c(r1,r2),pvalue=c(p1,p2),slope=c(s1,s2),Intercept=c(t1,t2),gender=c("M","F")),p))
        }))
        
        d1 <- huge4$Data2[[1]][[1]]
        q1 <- huge4$Data2[[1]][[2]]
        
        for (i in 2:nrow(huge4)){# i <- 2
          d1 <- bind_rows(d1,huge4$Data2[[i]][[1]])
          q1 <- q1 %>% aplot::insert_right(huge4$Data2[[i]][[2]])
        }
        ggsave(str_c(str_c(NM,collapse = "--"),".output.all.pdf"),q1,width = 15,height = 5)   
        write_tsv(d1 %>% add_column(.before = 1,Groups = huge$group),str_c(str_c(NM,collapse = "--"),".output.all.tsv"))    
        
      }
      
    }  
    
  }
  for (x in colnames(Env)[-1]){# x <- colnames(Env)[2]
    for (y in colnames(DATA)[-1]){# y <- colnames(DATA)[3]
      Data <- bind_cols(DATA[,y],Env[,x],Doo1[,y] %>% rename_at(1,~"QAZ")) 
      if(Data[,"QAZ"] %>% as.matrix %>% apply(.,2,sum)!=0){
        Analysis3(Data,map,DATA[,1])
      }
    }
  }
}

write_tsv(opts %>% as.data.frame %>% t %>% as.data.frame %>% rownames_to_column() %>% as_tibble,
          str_c("Parameter_spearman",
                str_replace_all(as.character(date())," ","_") %>% str_replace_all(":","_"),
                ".xls"),
          col_names = FALSE)
```

## Part6.XGB-sharp
```markdown
# XGBoost-SHAP 微生物组年龄预测分析脚本说明

## 简介

本脚本 `xgb-shap251216.py` 是一个完整的机器学习分析工具，专门用于基于微生物组数据预测人类年龄。它结合了 XGBoost 回归算法和 SHAP (SHapley Additive exPlanations) 解释框架，不仅能够建立准确的预测模型，还能深入解释哪些微生物特征对年龄预测最为重要。

## 功能特点

### 1. 数据处理与预处理
- 自动加载微生物组数据（CLR转换后的特征）和样本元数据（年龄信息）
- 数据质量控制和一致性验证
- 自动创建必要的输出目录结构

### 2. 模型训练与优化
- 使用嵌套交叉验证（Nested Cross-Validation）确保模型泛化能力
- 集成递归特征消除（RFE）进行特征选择
- 通过随机搜索（RandomizedSearchCV）优化超参数
- 支持模型保存和加载功能

### 3. 模型评估
- 计算多种评估指标：R²、MSE、RMSE、MAE
- 生成实际值vs预测值散点图
- 残差分析图
- 训练集和测试集性能对比

### 4. SHAP可解释性分析
- 全面的SHAP值计算和可视化
- 特征重要性分析（XGBoost内置和SHAP全局重要性）
- 个体预测解释（Waterfall图、Force图、Decision图）
- 特征依赖关系分析（Dependence图、散点图）
- 年龄组特异性分析（Age Group Aggregation Analysis）

### 5. 高级可视化功能
- SHAP摘要图（条形图、蜜蜂图）
- SHAP热力图
- 特征交互分析（依赖图对）
- 年龄组间差异分析

## 输入文件

### 必需文件
1. `input/clr_variance_filtered.xls` - 微生物组特征数据（已进行CLR转换）
   - 行代表样本，列代表微生物分类单元（ASV）
   - 第一列为OTU_tax标识符

2. `input/age.txt` - 样本元数据
   - 包含样本ID和对应年龄信息
   - 列名：#SampleID, age

## 输出文件

脚本运行后会在 `output/` 目录下生成以下内容：

### 模型相关
- `models/final_pipeline.pkl` - 训练好的完整模型管道
- `results/final_features.txt` - 最终选择的特征列表
- `results/final_model_params.json` - 最优超参数配置

### 评估结果
- `results/cross_validation_summary.csv` - 交叉验证结果汇总
- `results/sample_predictions.csv` - 每个样本的实际值和预测值

### 可视化图表
- `plots/` - 各类图表文件
  - 实际值vs预测值散点图
  - 残差分析图
  - XGBoost特征重要性图
  - SHAP各类解释图（PDF格式）
- `shap_values/shap_values.npy` - 计算的SHAP值

### SHAP分析专用
- `csv_paths/` - SHAP分析数据
  - SHAP全局重要性表
  - SHAP局部值表
  - 年龄组聚合分析结果

## 使用方法

### 基本运行
```bash
python xgb-shap251216.py
```
