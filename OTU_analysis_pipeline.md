# Part.1 OTU cluster

```shell

######################################## 01.quality control #######################################################################

seqkit stats *_R1.fastq -a -j 20 -T | sed 's/_R1.fastq//g' > seqkit.R1.out.xls
seqkit stats *_R2.fastq -a -j 20 -T | sed 's/_R2.fastq//g' > seqkit.R2.out.xls

mkdir qct
fastqc *fastq -t 10 -o qct
multiqc qct/*

less merge.list | awk '{print"cutadapt --cut -5 -o "$2".cut "$2"\ncutadapt --cut -50 -o "$3".cut "$3""}' > cut.sh
sh cut.sh

mv *.fastq rawdata/
rename 's/.cut//g' *

python ~/scripts/run_usearch_mergefilter.py -i merge.list -m 16 -r TRUE -l 400

mkdir cutadapt merge filter_fa filter_fq
mv *_merged.fastq merge
mv *_filtered.fasta filter_fa
mv *_filtered.fastq filter_fq
mv *fastq cutadapt
mv merge filter_fa filter_fq rawdata cutadapt data

cd data/filter_fa
ln -s ../../merge.list .
seqstat2excel.pl merge.list > ../../seqstat.xls 
cd ../../

cat seqstat.xls | sed '1d' | awk '{if($2>=100)print$1}' > analysis.list
cat seqstat.xls | sed '1d' | awk '{if($2<100)print$0}' > not_analysis.list

cat analysis.list | awk '{print"ln -s data/filter_fa/"$1"_filtered.fasta ./"}'>ln.sh
sh ln.sh
cat *_filtered.fasta >trimed.fasta
rm *_filtered.fasta

mothur "#pcr.seqs(fasta=trimed.fasta,oligos=~/scripts/oligo_v3v4.txt,keepprimer=TRUE,pdiffs=2,rdiffs=2,processors=10)"
sed -i "s/\t\tfpdiffs=.*rpdiffs=.*//g" trimed.pcr.fasta

############################################## 02.Analysis ######################################################################
mkdir Analysis_File_Directory
mkdir Analysis_File_Directory/process
cat trimed.pcr.fasta|awk -F ";" '{print $1}'|awk '{if(/>/){if(NR!=1){print SEQ "\n"$0;SEQ="";}else{print $0}}else{SEQ=SEQ""$0}}END{print SEQ}' >Analysis_File_Directory/process/meta.fasta
cd Analysis_File_Directory/process/
~/mothur/vsearch -derep_prefix meta.fasta -output meta_derepprefix.fasta -sizeout 
~/mothur/vsearch -sortbysize meta_derepprefix.fasta -output meta_derepprefix_sorted.fasta -minsize 2
mkdir otu_0.97
uparse -cluster_otus meta_derepprefix_sorted.fasta -otus otu_0.97/cluster.fasta -otu_radius_pct 3
cd otu_0.97

source /root/anaconda3/etc/profile.d/conda.sh
conda activate qiime2
python ~/scripts/fasta.multi2one.row.py -f cluster.fasta -o rep-seqs.fasta
qiime tools import --type 'FeatureData[Sequence]' --input-path rep-seqs.fasta --output-path rep-seqs.qza
qiime feature-classifier classify-sklearn \
                --i-classifier ~/database/qiime2_database/silva-138-99-nb-classifier.qza  \
                --i-reads rep-seqs.qza \
                --o-classification cluster_tax_assignments.qza \
                --p-n-jobs 20 \
                --verbose
qiime tools export --input-path cluster_tax_assignments.qza --output-path assign_taxonomy
sed '1d' assign_taxonomy/taxonomy.tsv | sed 's/ //g' > assign_taxonomy/cluster_tax_assignments.txt
source /root/anaconda3/etc/profile.d/conda.sh
conda deactivate

cat assign_taxonomy/cluster_tax_assignments.txt|awk '{print "OTU"NR"\t"$2"\t"$3;print "OTU"NR"\t"$1>"otu_reps.rename"}'>otu_reps_tax_assignments.txt
#map reads to otus
awk 'ARGIND==1{A[$2]=$1}ARGIND==2{if(/>/){s=0;gsub(/>/,"",$1);if(A[$1]){s=1;print ">"A[$1]}}else{if(s==1){print}}}'  otu_reps.rename  cluster.fasta >otu_reps.fasta
uparse -usearch_global ../meta.fasta -db otu_reps.fasta -strand plus -id 0.97 -uc map.uc
awk 'ARGIND==1{A[$2]=$1}ARGIND==2{if(/>/){s=0;gsub(/>/,"",$1);if(A[$1]){s=1;print ">"A[$1]"	"$1}}else{if(s==1){print}}}'  otu_reps.rename  cluster.fasta |awk -F ";size" '{print $1}' >otu_reps.fasta
sed 's/\t/_/' -i otu_reps.fasta
cat otu_reps.fasta |grep --color '>'|sed 's/>//' |sed 's/^OTU//'|sort -n|sed 's/^/OTU/' >otu_reps.accnos
choose_seqs.pl -f otu_reps.fasta -l otu_reps.accnos -o otu_reps.sort.fasta 
mv otu_reps.sort.fasta  otu_reps.fasta
uc2otuseqids.pl -i map.uc -o otu_seqids.txt 
cat otu_seqids.txt|awk '{split($0,line,"\t");new=line[1];for(i=2;i<NF+1;i++){match(line[i],/_[^_]+$/);smp=substr(line[i],1,RSTART-1);id=substr(line[i],RSTART+1,RLENGTH);nsmp=smp;gsub(/_/,".",nsmp);new=new"\t"nsmp"_"id;print nsmp"\t"smp;print line[i]"\t"smp >"otus.groups"};print new >"otu_seqids.tmp";}'|sort|uniq >name.check
make_otu_table.py -i otu_seqids.tmp  -o otu_table.biom
cat name.check|awk '{gsub(/\./,"\\.",$1);print "sed '\''s/\""$1"\"/\""$2"\"/g'\''      otu_table.biom >otu_table.biom.tmp\nmv otu_table.biom.tmp otu_table.biom";}'>name.check.sh
sed -i 's/"//g' name.check.sh
sh name.check.sh
biom convert -i otu_table.biom -o otu_table.txt  --table-type "OTU table" --to-tsv
cat otu_table.txt|sed -n '2p'|sed 's/#//' >otu_table.xls.tmp
cat otu_table.txt|sed -n '3,$p'|sort -V |sed 's/\.0//g'>>otu_table.xls.tmp
less otu_table.xls.tmp |grep --color 'OTU ID'|sed 's/OTU ID//'|xargs -n1|sort -n|awk '{print $1"\t"$1}'>sample_order
OTU_table_sortBySam.pl otu_table.xls.tmp sample_order otu_table.xls

## fix taxonomy
sed -i 's/ //g' otu_reps_tax_assignments.txt
awk 'NR==FNR{ a[$1] }NR>FNR{ if($1 in a) print $0}' otu_table.xls otu_reps_tax_assignments.txt > otu_select_fix_cluster_tax_assignments.txt

make_otu_table.py -i otu_seqids.tmp -t otu_select_fix_cluster_tax_assignments.txt -o otu_taxa_table.biom 
cat name.check|awk '{gsub(/\./,"\\.",$1);print "sed '\''s/\""$1"\"/\""$2"\"/g'\''      otu_taxa_table.biom >otu_taxa_table.biom.tmp\nmv otu_taxa_table.biom.tmp otu_taxa_table.biom";}' >name.check.sh
sed -i 's/"//g' name.check.sh
sh name.check.sh 
biom convert -i otu_taxa_table.biom -o otu_taxa_table.txt --header-key taxonomy  --table-type "OTU table" --to-tsv
cat otu_taxa_table.txt|sed -n '2p'|sed 's/#//' >otu_taxa_table.xls.tmp
cat otu_taxa_table.txt|sed -n '3,$p'|sort -V|sed 's/\.0//g' >>otu_taxa_table.xls.tmp
cp sample_order sample_tax_order
echo -e "taxonomy\ttaxonomy">>sample_tax_order
OTU_table_sortBySam.pl otu_taxa_table.xls.tmp sample_tax_order otu_taxa_table.xls

/usr/local/R-3.6.0/bin/R -e '.libPaths("/usr/local/R-3.6.0/lib64/R/library"); library(tidyverse);  da <- read_tsv("otu_taxa_table.xls") %>% rename("ASVID" = colnames(.)[1]) %>% select(-taxonomy) %>% gather(ID,value,-ASVID) %>% group_by(ID) %>% summarise_at(vars(c("value")),sum) %>% pull(value) %>% min ; oo <- c(10000,100,8000,125,5000,200,4000,250,2000,500,1000,1000,500,2000,400,2500,200,5000,100,10000) %>% matrix(ncol=2,byrow=T) %>% as_tibble %>% filter(V1 <= da) %>% .[1,] %>% write_tsv("rarefactions_even",col_names=F);'
####################################### 03.rarefactions_even_depth ###################################################
multiple_rarefactions_even_depth.py -i otu_table.biom -o rarefied_otu_tables -d `less rarefactions_even | awk '{print$1}'` -n `less rarefactions_even | awk '{print$2}'`
cd rarefied_otu_tables/
ls rarefaction_`less ../rarefactions_even | awk '{print$1}'`* | xargs -n `less ../rarefactions_even | awk '{print$2}'` | sed 's/ /,/g' > rare.biom.list
merge_otu_tables.py -i `cat rare.biom.list` -o ../rarefaction.merged_otu_table.biom 
cd ../
biom convert -i rarefaction.merged_otu_table.biom -o rarefaction.merged_otu_table.txt  --table-type "OTU table" --to-tsv
cat rarefaction.merged_otu_table.txt|sed -n '2p'|sed 's/#//' >rarefaction.merged_otu_table.xls.tmp
cat rarefaction.merged_otu_table.txt|sed -n '3,$p'|sort -V|sed 's/\.0//g' >>rarefaction.merged_otu_table.xls.tmp
OTU_table_sortBySam.pl rarefaction.merged_otu_table.xls.tmp sample_order rarefaction.merged_otu_table.xls

less rarefaction.merged_otu_table.xls|grep -v 'OTU ID' > rarefaction.merged_otu_table.tmp
less otu_taxa_table.xls|grep 'OTU ID' > rarefaction.merged_otu_taxa_table.xls
paste rarefaction.merged_otu_table.tmp otu_select_fix_cluster_tax_assignments.txt|awk '{$(NF-2)="";$(NF)="";print}'|awk -vOFS="\t" '{$1=$1}1' >> rarefaction.merged_otu_taxa_table.xls
biom convert -i rarefaction.merged_otu_taxa_table.xls -o rarefaction.merged_otu_taxa_table.biom --table-type "OTU table" --process-obs-metadata taxonomy --to-hdf5
rename 's/rarefaction.merged_otu/rarefac.otu/' rarefaction.merged_otu_*
vip_taxa_table_format.pl rarefac.otu_taxa_table.xls >rarefac.otu_genus.xls

################################### 04.Beta Diversity & Rarefaction curves  ###############################################

less otu_reps.fasta |awk -F '_' '{print $1}' >otu_reps.raw.fasta
align_seqs.py -i otu_reps.raw.fasta -m muscle -o .
make_phylogeny.py -i otu_reps.raw_aligned.fasta -o otu_reps_aligned.fasta.tre
cp /work/scripts/16s/qiime_parameters_97.txt ./

echo -e "#SampleID\tgroup">map.txt 
less rarefac.otu_table.xls|grep --color 'OTU ID'|sed 's/OTU ID//'|xargs -n1|sort -n|awk '{print $1"\t"$1}' >>map.txt
beta_diversity_through_plots.py -i rarefac.otu_table.biom -m map.txt -o beta_diversity -p qiime_parameters_97.txt -t otu_reps_aligned.fasta.tre

################################### 05.Alpha Diversity & Rarefaction curves  ###############################################
format_list.pl -i  otu_seqids.txt -l 0.97>otus.list
mothur "#make.shared(list=otus.list,group=otus.groups)"
rm otus.*.rabund

mkdir alpha_rarefac
cp otus.shared alpha_rarefac/
cd alpha_rarefac
mothur "#summary.single(shared=otus.shared,calc=ace-chao-shannon-simpson-coverage,groupmode=f)"
mothur "#rarefaction.single(shared=otus.shared,calc=sobs-chao-shannon-simpson,groupmode=f,freq=100,processors=30)"
python ~/scripts/plot_alpha_diversity.py -i rarefaction
python ~/scripts/plot_alpha_diversity.py -i r_shannon
shannon-ace-table.pl -d . -o estimators.html
cd ../

############################################### 06.picrust2 ##################################################################
source /root/anaconda3/etc/profile.d/conda.sh
conda activate picrust2
picrust2_pipeline.py -s otu_reps.raw.fasta -i otu_table.biom -o Picrust2 --processes 30 --in_traits COG,EC,KO,PFAM,TIGRFAM

cd Picrust2
find . -name "*.gz" | xargs gunzip

add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv      -m EC      -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv
add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv      -m KO      -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv
add_descriptions.py -i COG_metagenome_out/pred_metagenome_unstrat.tsv     -m COG     -o COG_metagenome_out/pred_metagenome_unstrat_descrip.tsv
add_descriptions.py -i PFAM_metagenome_out/pred_metagenome_unstrat.tsv    -m PFAM    -o PFAM_metagenome_out/pred_metagenome_unstrat_descrip.tsv
add_descriptions.py -i TIGRFAM_metagenome_out/pred_metagenome_unstrat.tsv -m TIGRFAM -o TIGRFAM_metagenome_out/pred_metagenome_unstrat_descrip.tsv
add_descriptions.py -i pathways_out/path_abun_unstrat.tsv                 -m METACYC -o pathways_out/path_abun_unstrat_descrip.tsv

pathway_pipeline.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv \
    -o KEGG_pathways_out --no_regroup \
    --map /work/softwares/picrust2-2.3.0-b/picrust2/default_files/pathway_mapfiles/KEGG_pathways_to_KO.tsv
# æ·»åŠ åŠŸèƒ½æè¿°
add_descriptions.py -i KEGG_pathways_out/path_abun_unstrat.tsv.gz \
    --custom_map_table /work/softwares/picrust2-2.3.0-b/picrust2/default_files/description_mapfiles/KEGG_pathways_info.tsv.gz \
    -o KEGG_pathways_out/path_abun_unstrat_descrip.tsv

# ç”Ÿæˆkegg modules ä¸°åº¦è¡¨ # https://github.com/picrust/picrust2.git
pathway_pipeline.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv \
    -o KEGG_modules_out --no_regroup \
    --map /work/softwares/picrust2-2.3.0-b/picrust2/default_files/pathway_mapfiles/KEGG_modules_to_KO.tsv
# æ·»åŠ åŠŸèƒ½æè¿°
add_descriptions.py -i KEGG_modules_out/path_abun_unstrat.tsv.gz \
    --custom_map_table /work/softwares/picrust2-2.3.0-b/picrust2/default_files/description_mapfiles/KEGG_modules_info_377.tsv.gz \
    -o KEGG_modules_out/path_abun_unstrat_descrip.tsv

cd ..

######## collect results  ##########
## now in Analysis_File_Directory/process/
cd ../../../
mkdir Analysis_File_Directory/results
mkdir Analysis_File_Directory/results/OTU_Taxa 
mkdir Analysis_File_Directory/results/Estimators
mkdir Analysis_File_Directory/results/Rarefactions
mkdir Analysis_File_Directory/results/Community
mkdir Analysis_File_Directory/results/Shannon_rarefac
mkdir Analysis_File_Directory/results/Beta_diversity
mkdir Analysis_File_Directory/results/Picrust2

cd Analysis_File_Directory/process/
cd  otu_0.97
cp  tax_summary_a/*.xls tax_summary_a/*.pdf  ../../results/Community/
cp tax_summary_r/*.percents.xls ../../results/Community/
rm  tax_summary_a/*.pdf tax_summary_a/ALL.new.*.xls 
cp -r otu_table.biom otu_taxa_table.biom otu_table.xls otu_reps.fasta otu_seqids.txt  otu_taxa_table.xls  rarefac.otu_table.biom rarefac.otu_table.xls rarefac.otu_taxa_table.biom rarefac.otu_taxa_table.xls tax_summary_a tax_summary_r rarefac.otu_genus.xls ../../results/OTU_Taxa/

cd alpha_rarefac/
cp otus.*.summary estimators.html ../../../results/Estimators
rename summary summary.xls ../../../results/Estimators/*.summary
cp -r otus.*.rarefaction rarefaction.*.pdf ../../../results/Rarefactions/
rename rarefaction rarefaction.xls ../../../results/Rarefactions/*.rarefaction
cp -r otus.*.r_shannon r_shannon.*.pdf ../../../results/Shannon_rarefac/
rename r_shannon r_shannon.xls ../../../results/Shannon_rarefac/*.r_shannon
cd ../../

cp otu_0.97/beta_diversity/*.txt ../results/Beta_diversity/
cp otu_0.97/heatmap.otu.top50.pdf otu_0.97/rarefac.otu_genus.xls ../results/Heatmap/
cp -r otu_0.97/Picrust2/*_out  ../results/Picrust2/

cd ../../

```

# Part.2 group analysis

```shell

########## set default  ###########
export Colour=none
export Paired=FALSE
export Combn=TRUE
export Paint=none
export Threshold=0.97
export Heatmap_w="30"
export Heatmap_h="10"
export Heatmap_keyh="14"
export Heatmap_marble="4-0-0-20"
export Hclust_bar_h="30"
export Hcluster_tree_h="30"

cat map-group.txt | grep -v "SampleID.*group" | awk '{print $1"\t"$1}' >order
grep -v "SampleID.*group" map-group.txt | sed '1i#SampleID\tgroup' > map-group.txt.tmp
rm map-group.txt
mv map-group.txt.tmp map-group.txt
NUM=`awk -F "\t" '{if(ARGIND==1) {val[$1]}else{if($1 in val)  delete val[$1]}}END{for(i in val) print i}' order ../process/otu_$Threshold/sample_order | wc -l`
if [ $NUM -gt 0 ]; then echo "command failed"; exit 1; fi

mkdir Split_groups
cd Split_groups
ln -s ../map-group.txt .
/usr/local/R-3.6.0/bin/Rscript ~/scripts/split_group.R -g map-group.txt -m none -n 0
cd ..

mkdir OTU_analysis
cd OTU_analysis
ln -s $PP1/process/otu_$Threshold/rarefac.otu_taxa_table.xls ./rarefac.otu_taxa_table1.xls
ln -s $PP2/process/rarefac.otu_taxa_table.xls                ./rarefac.otu_taxa_table2.xls
/opt/R-3.6.3/bin/R -e '
library(tidyverse)
library(data.table)
otu1 <- fread("rarefac.otu_taxa_table1.xls")
otu2 <- fread("rarefac.otu_taxa_table2.xls")
otu <- full_join(otu1,otu2)
otu[is.na(otu)]<-0
otu <- cbind(otu %>% select(-taxonomy),otu[,"taxonomy"])
write_tsv(otu,"rarefac.otu_taxa_table.xls")
'

sed -i 's/OTU ID/OTU_ID/g' rarefac.otu_taxa_table.xls
cp $PP1/process/otu_$Threshold/otu_reps.fasta .
cp ../map-group.txt .
less otu_reps.fasta |awk -F '_' '{print $1}' > otu_reps.raw.fasta
python ~/scripts/tax_split.v2.0.2.py -i rarefac.otu_taxa_table.xls -m map-group.txt -r otu_reps.raw.fasta -t 0.001

mkdir OTU_tree
cd OTU_tree
ln -s ../map.otu_reps.raw.fasta ./otu_reps.raw.fasta
ln -s ../unrooted.otu_tree_anno.0.001.xls .
muscle -in otu_reps.raw.fasta -out otu_reps.raw_aligned.fasta
FastTree -nt otu_reps.raw_aligned.fasta > unroot.otu_reps.raw_aligned.fasta.tre
/usr/local/R-3.6.0/bin/Rscript ~/scripts/unroot_tree.R -t unroot.otu_reps.raw_aligned.fasta.tre -a unrooted.otu_tree_anno.0.001.xls -s phylum
rm Rplots.pdf
cd ..

cd Krona
for i in *.xls; do ktImportText $i -o $i.html;done
cd ..

mkdir Core_Microbiome
cd Core_Microbiome
ln -s ../otu.genus.xls ./rarefac.otu_genus.xls
ln -s ../map-group.txt .
sed 's/\t$//' rarefac.otu_genus.xls > fix.rarefac.otu_genus.xls
sed -i 's/ //g' fix.rarefac.otu_genus.xls
/usr/local/R-3.6.0/bin/Rscript ~/scripts/core_microbiome.R -i fix.rarefac.otu_genus.xls -m map-group.txt

cd ../../
mv OTU_analysis/Core_Microbiome  OTU_analysis/Krona OTU_analysis/OTU_tree .

rm -rf OTU_Taxa
mv OTU_analysis OTU_Taxa
cd OTU_Taxa
rm rarefac.otu_taxa_table.xls otu_reps.raw.fasta otu_reps.fasta
mv map.otu_reps.raw.fasta otu_reps.raw.fasta
mv map.otu_table.xls rarefac.otu_table.xls
mv map.otu_table.percent.xls rarefac.otu_table.percent.xls
mv map.rarefac.otu_taxa_table.xls rarefac.otu_taxa_table.xls
mv otu.genus.xls rarefac.otu_genus.xls

/usr/local/R-3.6.0/bin/Rscript ~/scripts/Community_plot_test.R -t rank -r $Paired -o F -i rarefac.otu_genus.xls  -n 2 -p 30  -m map-group.txt  -c $Colour --combn $Combn

cd ..

#Specaccum
mkdir Specaccum
cd Specaccum
ln -s ../OTU_Taxa/rarefac.otu_table.xls .
ln -s ../map-group.txt .
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Specaccum.R -i rarefac.otu_table.xls -m map-group.txt -c $Colour
rm *.xls
cd ../

#Venn
mkdir Venn
cd Venn
ln -s ../OTU_Taxa/rarefac.otu_genus.xls .
ln -s ../map-group.txt
/usr/local/R-4.0.5/bin/Rscript ~/scripts/venn_upset.R -i rarefac.otu_genus.xls -g map-group.txt -m all -c $Colour
cd ../

mkdir Corrplot
cd Corrplot
ln -s ../OTU_Taxa/rarefac.otu_genus.xls .
ln -s ../map-group.txt .
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Spearman_corrplot.R -i rarefac.otu_genus.xls -m map-group.txt -o otu_top50
cd ..

#Rank_abundance
mkdir Rank_abundance
cd Rank_abundance
ln -s ../OTU_Taxa/rarefac.otu_table.xls
perl ~/scripts/rank_abundance.pl -i rarefac.otu_table.xls -gd ../map-group.txt  -o rankabundance.group.pdf -w 6 -h 5 -color $Colour
rm rarefac.otu_table.xls
cd ../

#Alpha_rarefac
mkdir Alpha_rarefac
cp $PP1/process/otu_0.97hhh/otu_0.97/subsample/alpha_subsample/*.r_chao       Alpha_rarefac -rl
cp $PP1/process/otu_0.97hhh/otu_0.97/subsample/alpha_subsample/*.r_shannon    Alpha_rarefac -rl
cp $PP1/process/otu_0.97hhh/otu_0.97/subsample/alpha_subsample/*.r_simpson    Alpha_rarefac -rl
cp $PP1/process/otu_0.97hhh/otu_0.97/subsample/alpha_subsample/*.rabund       Alpha_rarefac -rl
cp $PP1/process/otu_0.97hhh/otu_0.97/subsample/alpha_subsample/*.rarefaction  Alpha_rarefac -rl
cp $PP1/process/otu_0.97hhh/otu_0.97/subsample/alpha_subsample/*.summary      Alpha_rarefac -rl
cp $PP2/process/subsample/alpha_subsample/*.r_chao       Alpha_rarefac -rl
cp $PP2/process/subsample/alpha_subsample/*.r_shannon    Alpha_rarefac -rl
cp $PP2/process/subsample/alpha_subsample/*.r_simpson    Alpha_rarefac -rl
cp $PP2/process/subsample/alpha_subsample/*.rabund       Alpha_rarefac -rl
cp $PP2/process/subsample/alpha_subsample/*.rarefaction  Alpha_rarefac -rl
cp $PP2/process/subsample/alpha_subsample/*.summary      Alpha_rarefac -rl

cd Alpha_rarefac/
ln -s ../map-group.txt ./
python ~/scripts/plot_alpha_diversity.py -i rarefaction
python ~/scripts/plot_alpha_diversity.py -i r_shannon
python ~/scripts/plot_alpha_diversity.py -i rarefaction -m map-group.txt -c $Colour
python ~/scripts/plot_alpha_diversity.py -i r_shannon -m map-group.txt -c $Colour
rm ../Rarefactions/*.pdf ../Rarefactions/otus.*.rarefaction
rm ../Shannon_rarefac/*.pdf ../Shannon_rarefac/otus.*.r_shannon
mkdir ../Rarefactions ../Shannon_rarefac
mv rarefaction_* ../Rarefactions
mv r_shannon_* ../Shannon_rarefac

ln -s ../map-group.txt .
alpha_rarefac.pl map-group.txt > alpha_rarefac.summary.xls
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Alpha_diversity_box_plot.R -a alpha_rarefac.summary.xls -g map-group.txt -c $Colour

#Wilcox-alpha
/usr/local/R-3.6.0/bin/Rscript ~/scripts/alpha_diversity_test.R -a alpha_rarefac.summary.xls -g map-group.txt -t $Paired
rm *.r mothur* otus.* *.txt
cd ../

cd $PP2/process
cp $PP1/process/otu_0.97/map.txt ./
less rarefac.otu_table.xls|grep --color 'OTU ID'|sed 's/OTU ID//'|xargs -n1|sort -n|awk '{print $1"\t"$1}' >>map.txt
ln -s $PP1/process/otu_0.97/rarefac.otu_table.xls rarefac.otu_table1.xls
python ~/scripts/merge_table.py -i1 rarefac.otu_table1.xls -i2 rarefac.otu_table.xls -d1 none -d2 none -o new.rarefac.otu_table.xls
biom convert -i new.rarefac.otu_table.xls -o new.rarefac.otu_table.biom --table-type="OTU table" --to-hdf5

ln -s $PP1/process/otu_0.97/otu_reps.raw_aligned.fasta.tre ./
cp /work/scripts/16s/qiime_parameters_97.txt ./qiime_parameters_970.txt
beta_diversity_through_plots.py -i new.rarefac.otu_table.biom -m map.txt -o beta_diversity -p qiime_parameters_970.txt -t otu_reps.raw_aligned.fasta.tre
cd -

conda deactivate

rm Beta_diversity/ -r
mkdir Beta_diversity
cp $PP2/process/beta_diversity/*.txt Beta_diversity/

#Hcluster_tree 
mkdir Hcluster_tree
cd Hcluster_tree 
ln -s ../Beta_diversity/bray_curtis_dm.txt ../Beta_diversity/unweighted_unifrac_dm.txt ../Beta_diversity/weighted_unifrac_dm.txt ../map-group.txt ./ 
/usr/local/R-3.6.0/bin/Rscript ~/scripts/plot_tree.R -i bray_curtis_dm.txt -m map-group.txt -t average -c $Colour
/usr/local/R-3.6.0/bin/Rscript ~/scripts/plot_tree.R -i unweighted_unifrac_dm.txt -m map-group.txt -t average -c $Colour
/usr/local/R-3.6.0/bin/Rscript ~/scripts/plot_tree.R -i weighted_unifrac_dm.txt -m map-group.txt -t average -c $Colour
rm *.txt *.r Rplots.pdf
cd ../

#Beta_diversity
cd Beta_diversity
ln -s ../map-group.txt
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Matrix_analysis.R -i bray_curtis_dm.txt        -g map-group.txt -c $Colour  -m $Combn
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Matrix_analysis.R -i unweighted_unifrac_dm.txt -g map-group.txt -c $Colour  -m $Combn
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Matrix_analysis.R -i weighted_unifrac_dm.txt   -g map-group.txt -c $Colour  -m $Combn

/usr/local/R-3.6.0/bin/Rscript ~/scripts/plot_matrix.R -j FALSE -i bray_curtis_dm.txt        -g map-group.txt -c $Colour
/usr/local/R-3.6.0/bin/Rscript ~/scripts/plot_matrix.R -j FALSE -i unweighted_unifrac_dm.txt -g map-group.txt -c $Colour
/usr/local/R-3.6.0/bin/Rscript ~/scripts/plot_matrix.R -j FALSE -i weighted_unifrac_dm.txt   -g map-group.txt -c $Colour

mkdir Anosim Adonis MRPP MANOVA
mv *anosim* Anosim
mv *adonis* Adonis
mv *mrpp*   MRPP
mv *MANOVA* MANOVA
cd ..

#Community
rm -rf Community/
mv OTU_Taxa/Community .
cd Community
ln -s ../order ./sample
sed -i 's/ //g' *.xls
if [[ -f "../paint.txt" ]];then ln -s ../paint.txt ;fi

mkdir Community_barplot
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Community_plot_test.R -i phylum.xls  -n 2 -p 50 --paint $Paint
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Community_plot_test.R -i class.xls   -n 2 -p 50 --paint $Paint
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Community_plot_test.R -i order.xls   -n 2 -p 50 --paint $Paint
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Community_plot_test.R -i family.xls  -n 2 -p 50 --paint $Paint
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Community_plot_test.R -i genus.xls   -n 2 -p 50 --paint $Paint
mv *.pdf percent.*.xls ./Community_barplot/

ln -s ../order ./sample
cp ../map-group.txt ./
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Community_plot_test.R -a 2 -i phylum.xls -n 2 -p 50 -m map-group.txt --paint $Paint
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Community_plot_test.R -a 2 -i class.xls  -n 2 -p 50 -m map-group.txt --paint $Paint
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Community_plot_test.R -a 2 -i order.xls  -n 2 -p 50 -m map-group.txt --paint $Paint
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Community_plot_test.R -a 2 -i family.xls -n 2 -p 50 -m map-group.txt --paint $Paint
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Community_plot_test.R -a 2 -i genus.xls  -n 2 -p 50 -m map-group.txt --paint $Paint
mkdir Community_barplot_groups/
mv *.pdf sort.percent.* ./Community_barplot_groups/

mkdir Community_average
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Community_plot_test.R -a 1 -i phylum.xls -n 2 -p 50 -m map-group.txt --paint $Paint
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Community_plot_test.R -a 1 -i class.xls  -n 2 -p 50 -m map-group.txt --paint $Paint
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Community_plot_test.R -a 1 -i order.xls  -n 2 -p 50 -m map-group.txt --paint $Paint
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Community_plot_test.R -a 1 -i family.xls -n 2 -p 50 -m map-group.txt --paint $Paint
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Community_plot_test.R -a 1 -i genus.xls  -n 2 -p 50 -m map-group.txt --paint $Paint
mv *.pdf average.percent.* ./Community_average


mkdir Community_bubble
/usr/local/R-3.6.0/bin/Rscript ~/Scripts/Community_plot_test.R -b 2 -i phylum.xls -p 0   
/usr/local/R-3.6.0/bin/Rscript ~/Scripts/Community_plot_test.R -b 2 -i class.xls  -p 0.01
/usr/local/R-3.6.0/bin/Rscript ~/Scripts/Community_plot_test.R -b 2 -i order.xls  -p 0.01
/usr/local/R-3.6.0/bin/Rscript ~/Scripts/Community_plot_test.R -b 2 -i family.xls -p 0.01
/usr/local/R-3.6.0/bin/Rscript ~/Scripts/Community_plot_test.R -b 2 -i genus.xls  -p 0.01
mv *.pdf ./Community_bubble

mkdir Community_test
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Community_plot_test.R -t rank -r $Paired -o F -i phylum.xls -p 0.01 -m map-group.txt -c $Colour --combn $Combn
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Community_plot_test.R -t rank -r $Paired -o F -i class.xls  -p 0.01 -m map-group.txt -c $Colour --combn $Combn
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Community_plot_test.R -t rank -r $Paired -o F -i order.xls  -n 2 -p 30 -m map-group.txt -c $Colour --combn $Combn
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Community_plot_test.R -t rank -r $Paired -o F -i family.xls -n 2 -p 30 -m map-group.txt -c $Colour --combn $Combn
/usr/local/R-3.6.0/bin/Rscript ~/scripts/Community_plot_test.R -t rank -r $Paired -o F -i genus.xls  -n 2 -p 30 -m map-group.txt -c $Colour --combn $Combn
mv *.pdf *rank_sum*.xls Community_test
rm *.r percent*.xls

mkdir Heatmap_tax
plot-heatmap.pl -i Community_barplot_groups/sort.percent.genus.xls  -o heatmap.genus.pdf  -rt 0 -ct 0 -slas 2 -rlc 0.7 -clc 0.7 -w $Heatmap_w -h $Heatmap_keyh -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt -color $Colour
plot-heatmap.pl -i Community_barplot_groups/sort.percent.family.xls -o heatmap.family.pdf -rt 0 -ct 0 -slas 2 -rlc 0.7 -clc 0.7 -w $Heatmap_w -h $Heatmap_keyh -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt -color $Colour
plot-heatmap.pl -i Community_barplot_groups/sort.percent.order.xls  -o heatmap.order.pdf  -rt 0 -ct 0 -slas 2 -rlc 0.7 -clc 0.7 -w $Heatmap_w -h $Heatmap_keyh -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt -color $Colour
plot-heatmap.pl -i Community_barplot_groups/sort.percent.class.xls  -o heatmap.class.pdf  -rt 0 -ct 0 -slas 2 -rlc 0.7 -clc 0.7 -w $Heatmap_w -h $Heatmap_keyh -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt -color $Colour
plot-heatmap.pl -i Community_barplot_groups/sort.percent.phylum.xls -o heatmap.phylum.pdf -rt 0 -ct 0 -slas 2 -rlc 0.7 -clc 0.7 -w $Heatmap_w -h $Heatmap_keyh -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt -color $Colour

rm *.pdf.xls
mv heatmap*.pdf Heatmap_tax

#Community_boxplot
mkdir Community_boxplot
cd Community_boxplot
ln -s ../../map-group.txt ./
cp ../phylum.percent.xls .
cp ../genus.percent.xls  .
cp ../family.percent.xls .
cp ../order.percent.xls  .
cp ../class.percent.xls  .
mkdir phylum genus family class order
/usr/local/R-3.6.0/bin/Rscript ~/Projects/yangzhifu/chenyue/20220720002/analysis_Wed_Jul_20_12_32_59_2022/20220908001/æŒ‰GP7.FS-FM-FHH-FFMTåˆ†ç»„ç»“æžœ/Community/Community_boxplot/make_community_boxplot.r -i phylum.percent.xls -m map-group.txt -l 10 -p T -c $Colour
mv *.pdf phylum
/usr/local/R-3.6.0/bin/Rscript ~/Projects/yangzhifu/chenyue/20220720002/analysis_Wed_Jul_20_12_32_59_2022/20220908001/æŒ‰GP7.FS-FM-FHH-FFMTåˆ†ç»„ç»“æžœ/Community/Community_boxplot/make_community_boxplot.r -i genus.percent.xls  -m map-group.txt -l 10 -p T -c $Colour
mv *.pdf genus
/usr/local/R-3.6.0/bin/Rscript ~/Projects/yangzhifu/chenyue/20220720002/analysis_Wed_Jul_20_12_32_59_2022/20220908001/æŒ‰GP7.FS-FM-FHH-FFMTåˆ†ç»„ç»“æžœ/Community/Community_boxplot/make_community_boxplot.r -i order.percent.xls  -m map-group.txt -l 10 -p T -c $Colour
mv *.pdf order
/usr/local/R-3.6.0/bin/Rscript ~/Projects/yangzhifu/chenyue/20220720002/analysis_Wed_Jul_20_12_32_59_2022/20220908001/æŒ‰GP7.FS-FM-FHH-FFMTåˆ†ç»„ç»“æžœ/Community/Community_boxplot/make_community_boxplot.r -i family.percent.xls -m map-group.txt -l 10 -p T -c $Colour
mv *.pdf family
/usr/local/R-3.6.0/bin/Rscript ~/Projects/yangzhifu/chenyue/20220720002/analysis_Wed_Jul_20_12_32_59_2022/20220908001/æŒ‰GP7.FS-FM-FHH-FFMTåˆ†ç»„ç»“æžœ/Community/Community_boxplot/make_community_boxplot.r -i class.percent.xls  -m map-group.txt -l 10 -p T -c $Colour
mv *.pdf  class
rm *
cd ../../

# PCoA
mkdir PCoA
rm -rf Pcoa
cd PCoA
rm *
ln -s ../map-group.txt ./
ln -s ../Beta_diversity/bray_curtis_dm.txt ./
ln -s ../Beta_diversity/weighted_unifrac_dm.txt ./
ln -s ../Beta_diversity/unweighted_unifrac_dm.txt ./

#noname
mkdir noname
python ~/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt        -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt        -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt        -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt   -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt   -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt   -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt -md PCoA -pc 2-3 -map map-group.txt -col $Colour
mv *.pdf noname/

# 3D
mkdir 3D
/usr/local/R-3.6.0/bin/Rscript ~/scripts/plotly_PCoA-3D.R -i bray_curtis_dm.txt        -c $Colour
/usr/local/R-3.6.0/bin/Rscript ~/scripts/plotly_PCoA-3D.R -i weighted_unifrac_dm.txt   -c $Colour
/usr/local/R-3.6.0/bin/Rscript ~/scripts/plotly_PCoA-3D.R -i unweighted_unifrac_dm.txt -c $Colour
mv *.pdf 3D

#name
mkdir name
python ~/scripts/gg_PCoA_NMDS_PCA.py -lab T -i bray_curtis_dm.txt        -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -lab T -i bray_curtis_dm.txt        -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -lab T -i bray_curtis_dm.txt        -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -lab T -i weighted_unifrac_dm.txt   -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -lab T -i weighted_unifrac_dm.txt   -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -lab T -i weighted_unifrac_dm.txt   -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -lab T -i unweighted_unifrac_dm.txt -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -lab T -i unweighted_unifrac_dm.txt -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -lab T -i unweighted_unifrac_dm.txt -md PCoA -pc 2-3 -map map-group.txt -col $Colour
mv *.pdf name/

#ellipse
mkdir ellipse
python ~/scripts/gg_PCoA_NMDS_PCA.py -e T -i bray_curtis_dm.txt        -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -e T -i bray_curtis_dm.txt        -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -e T -i bray_curtis_dm.txt        -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -e T -i weighted_unifrac_dm.txt   -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -e T -i weighted_unifrac_dm.txt   -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -e T -i weighted_unifrac_dm.txt   -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -e T -i unweighted_unifrac_dm.txt -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -e T -i unweighted_unifrac_dm.txt -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -e T -i unweighted_unifrac_dm.txt -md PCoA -pc 2-3 -map map-group.txt -col $Colour
mv *.pdf ellipse/

#box
mkdir box
python ~/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i bray_curtis_dm.txt        -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i bray_curtis_dm.txt        -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i bray_curtis_dm.txt        -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i weighted_unifrac_dm.txt   -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i weighted_unifrac_dm.txt   -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i weighted_unifrac_dm.txt   -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i unweighted_unifrac_dm.txt -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i unweighted_unifrac_dm.txt -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i unweighted_unifrac_dm.txt -md PCoA -pc 2-3 -map map-group.txt -col $Colour
mv *.pdf box/

# crossbar
mkdir crossbar
python ~/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt         -ct 2 -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt         -ct 2 -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt         -ct 2 -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt    -ct 2 -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt    -ct 2 -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt    -ct 2 -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt  -ct 2 -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt  -ct 2 -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt  -ct 2 -md PCoA -pc 2-3 -map map-group.txt -col $Colour
mv *.pdf crossbar/

# central
mkdir central
python ~/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt         -ct 1 -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt         -ct 1 -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt         -ct 1 -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt    -ct 1 -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt    -ct 1 -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt    -ct 1 -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt  -ct 1 -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt  -ct 1 -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt  -ct 1 -md PCoA -pc 2-3 -map map-group.txt -col $Colour
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
#noname
mkdir noname
python ~/scripts/gg_CAP.py -i bray_curtis_dm.txt        -md CAP -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_CAP.py -i weighted_unifrac_dm.txt   -md CAP -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_CAP.py -i unweighted_unifrac_dm.txt -md CAP -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf noname/

#name
mkdir name
python ~/scripts/gg_CAP.py -lab T -i bray_curtis_dm.txt        -md CAP -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_CAP.py -lab T -i weighted_unifrac_dm.txt   -md CAP -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_CAP.py -lab T -i unweighted_unifrac_dm.txt -md CAP -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf name/

#ellipse
mkdir ellipse
python ~/scripts/gg_CAP.py -e T -i bray_curtis_dm.txt        -md CAP -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_CAP.py -e T -i weighted_unifrac_dm.txt   -md CAP -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_CAP.py -e T -i unweighted_unifrac_dm.txt -md CAP -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf ellipse/

#box
mkdir box
python ~/scripts/gg_CAP.py -bx T -e T -i bray_curtis_dm.txt        -md CAP -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_CAP.py -bx T -e T -i weighted_unifrac_dm.txt   -md CAP -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_CAP.py -bx T -e T -i unweighted_unifrac_dm.txt -md CAP -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf box/

# crossbar
mkdir crossbar
python ~/scripts/gg_CAP.py -i bray_curtis_dm.txt         -ct 2 -md CAP -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_CAP.py -i weighted_unifrac_dm.txt    -ct 2 -md CAP -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_CAP.py -i unweighted_unifrac_dm.txt  -ct 2 -md CAP -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf crossbar/

# central
mkdir central
python ~/scripts/gg_CAP.py -i bray_curtis_dm.txt         -ct 1 -md CAP -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_CAP.py -i weighted_unifrac_dm.txt    -ct 1 -md CAP -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_CAP.py -i unweighted_unifrac_dm.txt  -ct 1 -md CAP -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf central/

rm cmd.r 
cd ../

############################################################ PCA ##########################################################
mkdir PCA
rm Pca -rf
cd PCA
ln -s ../map-group.txt ./
ln -s ../OTU_Taxa/rarefac.otu_table.xls ./

#noname
mkdir noname
python ~/scripts/gg_PCoA_NMDS_PCA.py -i rarefac.otu_table.xls -md PCA -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i rarefac.otu_table.xls -md PCA -pc 1-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i rarefac.otu_table.xls -md PCA -pc 2-3 -map map-group.txt -col $Colour
mv PCA*.pdf noname/

/usr/local/R-3.6.0/bin/Rscript ~/scripts/plotly_PCoA-3D.R -i rarefac.otu_table.xls -d PCA

#name
mkdir name
python ~/scripts/gg_PCoA_NMDS_PCA.py -lab T -i rarefac.otu_table.xls -md PCA -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -lab T -i rarefac.otu_table.xls -md PCA -pc 1-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -lab T -i rarefac.otu_table.xls -md PCA -pc 2-3 -map map-group.txt -col $Colour
mv PCA*.pdf name/

#ellipse
mkdir ellipse
python ~/scripts/gg_PCoA_NMDS_PCA.py -e T -i rarefac.otu_table.xls -md PCA -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -e T -i rarefac.otu_table.xls -md PCA -pc 1-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -e T -i rarefac.otu_table.xls -md PCA -pc 2-3 -map map-group.txt -col $Colour
mv PCA*.pdf ellipse/

#box
mkdir box
python ~/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i rarefac.otu_table.xls -md PCA -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i rarefac.otu_table.xls -md PCA -pc 1-3 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i rarefac.otu_table.xls -md PCA -pc 2-3 -map map-group.txt -col $Colour
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

#noname
mkdir noname
python ~/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt        -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt   -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt -md NMDS -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf noname/

#name
mkdir name
python ~/scripts/gg_PCoA_NMDS_PCA.py -lab T -i bray_curtis_dm.txt        -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -lab T -i weighted_unifrac_dm.txt   -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -lab T -i unweighted_unifrac_dm.txt -md NMDS -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf name/

#ellipse
mkdir ellipse
python ~/scripts/gg_PCoA_NMDS_PCA.py -e T -i bray_curtis_dm.txt        -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -e T -i weighted_unifrac_dm.txt   -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -e T -i unweighted_unifrac_dm.txt -md NMDS -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf ellipse/

#box
mkdir box
python ~/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i bray_curtis_dm.txt        -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i weighted_unifrac_dm.txt   -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python ~/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i unweighted_unifrac_dm.txt -md NMDS -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf box/

rm cmd.r 
cd ../

##################################################### plsda  #############
mkdir PLS-DA
cd PLS-DA
ln -s ../map-group.txt ./
ln -s ../OTU_Taxa/rarefac.otu_table.xls ./

perl ~/scripts/plsda.pl -i rarefac.otu_table.xls -m map-group.txt -o ./ -g group -l F -c $Colour
cd ..
#### treebar ###################################################
mkdir Hclust_bar
cd Hclust_bar/
plot-treebar.pl -otu ../OTU_Taxa/rarefac.otu_table.xls -tax ../Community/phylum.xls -o treebar_otu_table_phylum.pdf -h $Hclust_bar_h -lcex 1.3
cd ..

######################################### Lefse #################################
source /root/anaconda3/etc/profile.d/conda.sh
conda activate qiime1
export R_HOME_DIR=/root/anaconda3/envs/qiime1/lib/R/library
export R_LIBS=$R_HOME_DIR/lib64/R/library

mkdir Lefse
cd Lefse/
cp ../OTU_Taxa/rarefac.otu_taxa_table.xls ./
sed -i 's/;s__.*//g' rarefac.otu_taxa_table.xls
python2 /work/scripts/16s/tax_split.py -i rarefac.otu_taxa_table.xls

sh ~/scripts/S_S/001.Community.OTU.ASV.lefse.sh  -c $Colour
cd ../
####################################################### Random_Forest
mkdir Random_Forest
cd Random_Forest/
ln -s ../OTU_Taxa/rarefac.otu_table.xls ./
ln -s ../OTU_Taxa/rarefac.otu_genus.xls ./
ln -s ../map-group.txt ./
ln -s ../order ./
biom convert -i rarefac.otu_table.xls -o rarefac.otu_table.biom --table-type "OTU table" --to-hdf5

random_forest4key_out_select.pl rarefac.otu_table.biom rarefac.otu_genus.xls  map-group.txt 0.001
sortSample4otu_table.pl OTU-extract.all.xls order > OTU-extract.all.sort.xls
cd ../

####################################################### Heatmap ####################################################################
cd Heatmap
rm ./*
ln -s ../OTU_Taxa/rarefac.otu_genus.xls ./
ln -s ../map-group.txt ./
plot-heatmap.pl -i rarefac.otu_genus.xls -o heatmap.otu.top50.pdf -rtop 50 -ct 0 -slas 2 -rlc 0.7 -clc 0.7 -w $Heatmap_w -h $Heatmap_h -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt -color $Colour
ln -s ../Random_Forest/OTU-extract.all.sort.xls ./
bar_pie.pl -i OTU-extract.all.sort.xls -pie F -p 0.01
sed -i '/Others/d' ALL.new.OTU-extract.all.sort.xls
mv ALL.new.OTU-extract.all.sort.xls ALL.new.OTU-extract.all.sort.percent001.xls
bar_pie.pl -i OTU-extract.all.sort.xls -pie F -p 0.03
sed -i '/Others/d' ALL.new.OTU-extract.all.sort.xls
mv ALL.new.OTU-extract.all.sort.xls ALL.new.OTU-extract.all.sort.percent003.xls
plot-heatmap.pl -i OTU-extract.all.sort.xls                    -o heatmap.keyotu.pdf            -ct 0 -slas 2 -rlc 0.7 -clc 0.7 -w $Heatmap_w -h $Heatmap_keyh -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt -color $Colour
plot-heatmap.pl -i ALL.new.OTU-extract.all.sort.percent001.xls -o heatmap.keyotu.percent001.pdf -ct 0 -slas 2 -rlc 0.7 -clc 0.7 -w $Heatmap_w -h $Heatmap_h    -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt -color $Colour
plot-heatmap.pl -i ALL.new.OTU-extract.all.sort.percent003.xls -o heatmap.keyotu.percent003.pdf -ct 0 -slas 2 -rlc 0.7 -clc 0.7 -w $Heatmap_w -h $Heatmap_h    -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt -color $Colour
rm ALL*.xls dat* bar* cmd.r ./*ct*
cd ..

mkdir Heatmap_PRETTY
cd Heatmap_PRETTY
ln -s ../map-group.txt ./
ln -s ../OTU_Taxa/rarefac.otu_genus.xls ./
python ~/scripts/pheatmap.py -i rarefac.otu_genus.xls -cd map-group.txt -o pheatmap.otu.top50.pdf -rtop 50 -clust_c F -cs 0 -rs 1 -scale none -cc $Colour
ln -s ../Random_Forest/OTU-extract.all.sort.xls ./
bar_pie.pl -i OTU-extract.all.sort.xls -pie F -p 0.01
sed -i '/Others/d' ALL.new.OTU-extract.all.sort.xls
mv ALL.new.OTU-extract.all.sort.xls ALL.new.OTU-extract.all.sort.percent001.xls
bar_pie.pl -i OTU-extract.all.sort.xls -pie F -p 0.03
sed -i '/Others/d' ALL.new.OTU-extract.all.sort.xls
mv ALL.new.OTU-extract.all.sort.xls ALL.new.OTU-extract.all.sort.percent003.xls
python ~/scripts/pheatmap.py -i OTU-extract.all.sort.xls                    -o pheatmap.keyotu.pdf            -cd map-group.txt -clust_c F -cs 0 -rs 1 -scale none -cc $Colour
python ~/scripts/pheatmap.py -i ALL.new.OTU-extract.all.sort.percent001.xls -o pheatmap.keyotu.percent001.pdf -cd map-group.txt -clust_c F -cs 0 -rs 1 -scale none -cc $Colour
python ~/scripts/pheatmap.py -i ALL.new.OTU-extract.all.sort.percent003.xls -o pheatmap.keyotu.percent003.pdf -cd map-group.txt -clust_c F -cs 0 -rs 1 -scale none -cc $Colour
rm ALL*.xls dat* bar* cmd.r ./*ct*
cd ../
##################################################### Picrust2 ################################################################################################################
if [ -d Picrust2 ];then
cd Picrust2/
sh ~/scripts/S_S/Picrust2.lefse.sh -c $Colour
cd ../
fi

conda deactivate
export R_HOME_DIR=/usr/local/R-3.6.0
export R_LIBS=$R_HOME_DIR/lib64/R/library

if [ -d Picrust2 ];then
if [ -d Picrust2/KEGG_pathways_out ];then
cd Picrust2/KEGG_pathways_out
mkdir KEGG_pathway_summary
mv KEGG_pathwaysL* KEGG_pathway_summary
cd KEGG_pathway_summary
/usr/local/R-3.6.0/bin/Rscript ~/Scripts/Community_plot_test.R -i KEGG_pathwaysLevel1_out.tsv -n 2 -p 49 -c $Colour --combn $Combn
/usr/local/R-3.6.0/bin/Rscript ~/Scripts/Community_plot_test.R -i KEGG_pathwaysLevel2_out.tsv -n 2 -p 49 -c $Colour --combn $Combn
/usr/local/R-3.6.0/bin/Rscript ~/Scripts/Community_plot_test.R -i KEGG_pathwaysLevel3_out.tsv -n 2 -p 49 -c $Colour --combn $Combn

/usr/local/R-3.6.0/bin/Rscript ~/Scripts/Community_plot_test.R -t rank -o F -i KEGG_pathwaysLevel1_out.tsv -n 0 -p 0.005  -m ../../map-group.txt -c $Colour --combn $Combn
/usr/local/R-3.6.0/bin/Rscript ~/Scripts/Community_plot_test.R -t rank -o F -i KEGG_pathwaysLevel2_out.tsv -n 0 -p 0.005  -m ../../map-group.txt -c $Colour --combn $Combn
/usr/local/R-3.6.0/bin/Rscript ~/Scripts/Community_plot_test.R -t rank -o F -i KEGG_pathwaysLevel3_out.tsv -n 0 -p 0.005  -m ../../map-group.txt -c $Colour --combn $Combn
rm Rplots.pdf
cd ../../../
fi
fi

```

# Part.3 linear regression

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

## 1.Analysis ä¸€ä¸ªèŒä¸‰ä¸ªå¹´é¾„æ®µåœ¨åŒä¸€å¼ å›¾ä¸Šçš„çº¿æ€§å›¾ ###########
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

## 2.Analysis2ä¸€ä¸ªèŒotuä¸°åº¦å’Œå¹´é¾„çš„æ•£ç‚¹å›¾å’Œæ¯ä¸€å²å¹´é¾„ä¸‹å¹³å‡ä¸°åº¦çš„æŠ˜çº¿ï¼ˆæ¯ä¸€å²ä¸‹å¦‚æžœå°‘äºŽ3ä¸ªæ ·æœ¬å°±ä¸çº³å…¥ï¼‰ ###########
## ï¼ˆ0-3å²å››èˆäº”å…¥ï¼Œ0.5ä»¥ä¸‹ç®—åš0å²ï¼Œåªç”¨å‡å€¼è¿žçº¿ï¼Œåˆ†ç”·å¥³ï¼Œä¸€ä¸ªå›¾ä¿©çº¿ï¼‰
## æ¯”å¦‚3å²çš„äººæœ‰6ä¸ªï¼Œå…¶ä¸­2å¥³4ç”·ï¼Œä¸¤ä¸ªéƒ½å±•ç¤º
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

## 3.Analysis3å±žæ°´å¹³åšä¸€ä¸‹ä¸‰ä¸ªå¹´é¾„æ®µèŒå±žä¸Žå¹´é¾„çš„çº¿æ€§ç›¸å…³å›¾ï¼Œç»˜åˆ¶çº¿æ€§ç›¸å…³å›¾çš„å±žè¦æ±‚ï¼š ###################
##ï¼ˆ1ï¼‰åœ¨ä»»æ„ä¸€ç»„ä¸­æœ‰è‡³å°‘75%çš„æ ·æœ¬ä¸­èƒ½æ£€æµ‹åˆ°ï¼›
##ï¼ˆ2ï¼‰è‡³å°‘åœ¨ä¸€ä¸ªç»„ä¸­ä¸Žå¹´é¾„æœ‰æ˜¾è‘—ç›¸å…³æ€§ï¼›å…­æ¡çº¿ éƒ½ä¸æ˜¾è‘— p>0.05 å°±ä¸è¦å›¾
##ï¼ˆ3ï¼‰å‡ºå›¾æ—¶è¦æ±‚æ ·æœ¬ç‚¹æŒ‰ç…§ç”·å¥³æ€§åˆ«åˆ†åˆ«ç”¨ä¸åŒé¢œè‰²å±•ç¤ºï¼Œä¸”ç”·å¥³åˆ†åˆ«åšä¸¤æ¡æ‹Ÿåˆç›´çº¿ï¼›
##ï¼ˆ4ï¼‰æŠŠæ‰€æœ‰ç”»å›¾çš„å±žåœ¨ä¸‰ç»„ä¸¤ä¸ªæ€§åˆ«ä¸­çš„çº¿æ€§ç›¸å…³æ•°å€¼ï¼ˆR2ã€på€¼ã€æ–œçŽ‡ï¼‰éƒ½æ•´ç†åˆ°ä¸€ä¸ªæ•°æ®è¡¨ä¸­ï¼›
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

## 4.OTUä¸°åº¦éƒ½å’Œå¹´é¾„åšä¸€ä¸‹ã€Šå¤šé¡¹å¼å›žå½’ã€‹ #############################
#ï¼ˆå…¨å¹´é¾„ä¸€èµ·ï¼Œä¸åˆ†å¹´é¾„æ®µï¼Œç”·å¥³åˆ†åˆ«åšå›žå½’ï¼Œç”»åœ¨ä¸€ä¸ªå›¾ä¸Šï¼‰
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


## 5.Analysis5 æ˜¯å‡çº§ç‰ˆçš„Analysis3 ###################
##ï¼ˆ1ï¼‰åœ¨ä»»æ„ä¸€ç»„ä¸­æœ‰è‡³å°‘75%çš„æ ·æœ¬ä¸­èƒ½æ£€æµ‹åˆ°ï¼›
##ï¼ˆ2ï¼‰è‡³å°‘åœ¨ä¸€ä¸ªç»„ä¸­ä¸Žå¹´é¾„æœ‰æ˜¾è‘—ç›¸å…³æ€§ï¼›å…­æ¡çº¿ éƒ½ä¸æ˜¾è‘— p>0.05 å°±ä¸è¦å›¾
##ï¼ˆ3ï¼‰å‡ºå›¾æ—¶è¦æ±‚æ ·æœ¬ç‚¹æŒ‰ç…§ç”·å¥³æ€§åˆ«åˆ†åˆ«ç”¨ä¸åŒé¢œè‰²å±•ç¤ºï¼Œä¸”ç”·å¥³åˆ†åˆ«åšä¸¤æ¡æ‹Ÿåˆç›´çº¿ï¼›
##ï¼ˆ4ï¼‰æŠŠæ‰€æœ‰ç”»å›¾çš„å±žåœ¨ä¸‰ç»„ä¸¤ä¸ªæ€§åˆ«ä¸­çš„çº¿æ€§ç›¸å…³æ•°å€¼ï¼ˆR2ã€på€¼ã€æ–œçŽ‡ï¼‰éƒ½æ•´ç†åˆ°ä¸€ä¸ªæ•°æ®è¡¨ä¸­ï¼›
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

# Part.4 random\_forest

```R
library(tidyverse)
library(magrittr)
library(randomForest)
library(pROC)
library(optparse)
library(broom)
library(patchwork)
library(MASS)
#library(randomForest, lib.loc = "/usr/local/R-4.0.5/lib64/R/library")
#library(pROC, lib.loc = "/usr/local/R-4.0.5/lib64/R/library")
######################### 01.parameters ##############################
rm(list=ls())
if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"),     type="character", default="genus.xls",help="input"),
    make_option(c("-w", "--wilcox"),    type="character", default="none", help="wilcox test"),
    make_option(c(      "--select"),    type="character", default="select.list", help="select.list"),
    make_option(c("-m", "--map"),       type="character", default="map-group.txt",help="group"),
    make_option(c("-n", "--cvnumber"),  type="numeric",   default=5,help="x fold"),
    make_option(c("-t", "--cvtime"),    type="numeric",   default=5, help="x times"),
    make_option(c("-p", "--per"),       type="numeric",   default=0,help="OTU filter"),
    make_option(c("-k", "--markernum"), type="numeric",   default=0,help="pick"),
    make_option(c("-f", "--feature"),   type="character", default="none", help="feature_importance_scores.txt"),
    make_option(c("-c", "--color"),     type="character", default="none",help="color"),
    make_option(c("-s", "--seed"),      type="character", default="1234", help="seed 1234"),
    make_option(c("-x", "--Meande"),    type="double",    default=0.001, help="feature_importance_scores threshold"),
    make_option(c("-y", "--pvalue"),    type="double",    default=0.05, help="pvalue"),
    make_option(c("-z", "--qvalue"),    type="double",    default=1,  help="qvalue"),
    make_option(c(      "--value"),     type="numeric",   default=0, help="core_microbiome"),
    make_option(c(      "--part"),      type="character", default="2/3", help="split 2/3 or 3/4 or 0.8"),
    make_option(c(      "--nnew"),      type="logical",   default=T, help="split.map-group.txt"),
    make_option(c("-v", "--valid"),     type="character", default="none", help="valid.rarefac.otu_table.xls"),
    make_option(c("-g", "--gp"),        type="character", default="none", help="map2.txt"),
    make_option(c("-r", "--paired"),    type="logical",   default=FALSE, help="paireded testï¼Ÿï¼Ÿï¼Ÿ")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
}

i <- opts$input
w <- opts$wilcox
m <- opts$map
cvn <- opts$cvnumber
cvt <- opts$cvtime
p <- opts$per
f <- opts$feature
c <- opts$color
marker.num <- opts$markernum

options(scipen = 200)
set.seed(opts$seed) # set.seed(2345) # set.seed(3456)

######################### 02.function ##############################
my.rfcv <- function (trainx, trainy, cv.fold = 5, scale = "log", step = 0.5, 
                     mtry = function(p) max(1, floor(sqrt(p))), recursive = FALSE, ...) 
{
  classRF <- is.factor(trainy)
  n <- nrow(trainx)
  p <- ncol(trainx)
  if (scale == "log") {
    k <- floor(log(p, base = 1/step))
    n.var <- round(p * step^(0:(k - 1)))
    same <- diff(n.var) == 0
    if (any(same)) 
      n.var <- n.var[-which(same)]
    if (!1 %in% n.var) 
      n.var <- c(n.var, 1)
    print(n.var)
  }
  else {
    n.var <- seq(from = p, to = 1, by = step)
  }
  k <- length(n.var)
  cv.pred <- vector(k, mode = "list")
  for (i in 1:k) cv.pred[[i]] <- trainy
  if (classRF) {
    f <- trainy
  }
  else {
    f <- factor(rep(1:5, length = length(trainy))[order(order(trainy))])
  }
  nlvl <- table(f)
  idx <- numeric(n)
  for (i in 1:length(nlvl)) {
    idx[which(f == levels(f)[i])] <- sample(rep(1:cv.fold, length = nlvl[i]))
  }
  
  res = list() 
  
  for (i in 1:cv.fold) {
    all.rf <- randomForest(trainx[idx != i, , drop = FALSE], 
                           trainy[idx != i], trainx[idx == i, , drop = FALSE], 
                           trainy[idx == i], mtry = mtry(p), importance = TRUE, ...)
    cv.pred[[1]][idx == i] <- all.rf$test$predicted
    #impvar <- (1:p)[order(all.rf$importance[, 1], decreasing = TRUE)]
    impvar <- (1:p)[order(randomForest::importance(all.rf, type = 1), decreasing = TRUE)]
    res[[i]] <- impvar
    for (j in 2:k) {
      imp.idx <- impvar[1:n.var[j]]
      sub.rf <- randomForest(trainx[idx != i, imp.idx, 
                                    drop = FALSE], trainy[idx != i], trainx[idx == i, imp.idx, drop = FALSE],
                             trainy[idx == i], mtry = mtry(n.var[j]), importance = recursive, ...)
      cv.pred[[j]][idx == i] <- sub.rf$test$predicted
      if (recursive) {
        #impvar <- (1:length(imp.idx))[order(sub.rf$importance[, 1], decreasing = TRUE)]
        impvar <- (1:length(imp.idx))[order(randomForest::importance(sub.rf, type = 1), decreasing = TRUE)]
      }
      NULL
    }
    NULL
  }
  if (classRF) {
    error.cv <- sapply(cv.pred, function(x) mean(trainy != x))
  }
  else {
    error.cv <- sapply(cv.pred, function(x) mean((trainy - x)^2))
  }
  names(error.cv) <- names(cv.pred) <- n.var
  list(n.var = n.var, error.cv = error.cv, predicted = cv.pred, res = res)
}

Minus <- function(x,n){
  d1 <- 10^(-n)
  ifelse(x >= d1,round(x,4),paste0("< ",d1))
}

Genus_threshold <- function(nc,otu_coverage,data2_1){
  genus_threshold <- enframe(sapply(1:length(nc),function(y) 
    sum(as.numeric(otu_coverage >= nc[y])))) %>% 
    mutate(remian = length(otu_coverage) - value,
           coverage = nc) %>%
    mutate(per = sapply(1:length(nc),function(o){
      sum(unlist(data2_1[otu_coverage >= .$coverage[o],1:(ncol(data2_1))])) / 
        sum(unlist(data2_1[,1:(ncol(data2_1))]))
    })) %>% dplyr::select(-name)
}

if (file.exists("projectimage.RData")){

  load(file = "projectimage.RData")
  
}else{
  data1 <- read.table(m,head= T ,sep="\t",comment.char = "",fileEncoding = "UTF-8")
  colnames(data1) <- c("SampleID","group")
  
  data2 <- read.table(i,head= T ,sep="\t",comment.char = "",row.names = 1,fileEncoding = "UTF-8")
  
  data2_1 <- data2[,as.character(data1$SampleID)]
  data2_1 <- sapply(1:ncol(data2_1),function(x) data2_1[,x] <- data2_1[,x]/sum(data2_1[,x]))
  rownames(data2_1) <- rownames(data2)
  colnames(data2_1) <- data1$SampleID
  
  otu_coverage <- apply(data2_1[,1:ncol(data2_1)],1,function(x) #20221212
    length(x[x>0])/ncol(data2_1))
  nc <- sort(unique(otu_coverage))
  
  Genus_threshold(nc,otu_coverage,data2_1)
  data2_1 <- data2_1[otu_coverage >= opts$value,] #20221212
  
  d2p <- sapply(1:nrow(data2_1),function(y) any(data2_1[y,]>=p))
  data2_2 <- data2_1[d2p,]
  
  if(opts$wilcox != "none"){
    data3 <- read.table(w,head= T ,sep="\t",comment.char = "",row.names = 1,fileEncoding = "UTF-8")
    d3p <- rownames(data3)[data3$p.value <= opts$pvalue]
    data2_3 <- data2_2[rownames(data2_2) %in% d3p,]
    
    if (opts$qvalue != 0){
      d3q <- rownames(data3)[data3$q.value <= opts$qvalue]
      data2_3 <- data2_3[rownames(data2_3) %in% d3q,]
    }
    Data <- cbind(t(data2_3),data1)
  }
  
  if (f != "none"){
    data4 <- read.table(f,head= T ,sep="\t",comment.char = "",row.names = 1,fileEncoding = "UTF-8")
    data4 <- data4[data4$Mean_decrease_in_accuracy >= opts$Meande,]
    data2_4 <- sapply(rownames(data2_3),function(x) strsplit(x," ", fixed=TRUE)[[1]][1])
    data2_5 <- data2_3[data2_4 %in% rownames(data4),]
    nrow(data2_5)
    Data <- cbind(t(data2_5),data1)
  }
  
  if(opts$select != "none"){ #20221221
    data5 <- read_tsv(opts$select,col_names = F)
    data2_6 <- data2_2[as.character(data5 %>% t),]
    Data <- cbind(t(data2_6),data1)
  }
  
  Data <- Data[,-(ncol(Data)-1)]
  
  if (c != "none"){
    sc <- read.table(c,sep="\t",comment.char = "",check.names = FALSE)
    sc <- sc[which(sc[,1] %in% unique(Data$group)),]
    mycol <- as.vector(sc[,2])
    Data$group <- factor(Data$group,levels=as.vector(sc[,1]))
  }else{
    Data$group <- factor(Data$group,levels=unique(as.vector(Data$group)))
  }
  Data <- dplyr::arrange(Data,group)
  
  ######################### 04.marker ##############################
  if(opts$nnew){ #20221212
    if (!opts$paired){
      if(grepl("/",opts$part)){ #20230109
        aa <- strsplit(opts$part,'/')[[1]][1] %>% as.numeric ; bb <- strsplit(opts$part,'/')[[1]][2] %>% as.numeric 
      }else{
        aa <- strsplit(MASS::fractions(opts$part %>% as.numeric) %>% as.character,'/')[[1]][1] %>% as.numeric 
        bb <- strsplit(MASS::fractions(opts$part %>% as.numeric) %>% as.character,'/')[[1]][2] %>% as.numeric 
      }
      #ind <- table(sample(2,nrow(Data),replace = TRUE,prob = c(0.67,0.33)))
      SPLIT <- function(DD){# DD <- Data
        set.seed(opts$seed)
        ind1 <- sample(1:nrow(DD),round(nrow(DD)*aa/bb),replace=FALSE,prob=1:nrow(DD))
        ind2 <- setdiff(1:nrow(DD),ind1)
        bind_rows(ind1 %>% as_tibble() %>% mutate(Split=1),ind2 %>% as_tibble() %>% mutate(Split=2)) %>% arrange(value) %>% pull(Split)
      }
      zz <- data1 %>% as_tibble() %>% group_by(group) %>% nest
      ind <- lapply(1:nrow(zz),function(x) SPLIT(zz$data[[x]]) %>% as_tibble()) %>% unlist() %>% unname()
    }else{
      ind <- rep(sample(2,nrow(Data)/2,replace = TRUE,prob = c(0.67,0.33)),2)
    }
    if(opts$part!="1/1"){data1 %>% mutate(Split = ind) %>% write_tsv("split.map-group.txt")}
  }else{
    ind <-read_tsv("split.map-group.txt",col_select = 3) %>% t %>% as.vector
  }
  
  # 2/3 best marker
  result <- replicate(cvt, my.rfcv(Data[ind==1,-ncol(Data),drop=FALSE], Data[ind==1,"group"],cv.fold = cvn,step=0.5), simplify=FALSE)
  error.cv <- sapply(result, "[[", "error.cv")
  
  error.cv.rm <- rowMeans(error.cv)
  id <- error.cv.rm < min(error.cv.rm) + sd(error.cv.rm)
  error.cv[id, ]
  if (marker.num == 0) { 
    marker.num <- min(as.numeric(names(error.cv.rm)[id]))
  }
  
  prefix <- paste(unique(as.vector(data1$group)),collapse = "-")
  pdf.dir1=paste0(prefix, "_vars.pdf")
  pdf.dir2=paste0(prefix, "_train_boxplot.pdf") 
  pdf.dir3=paste0(prefix, "_train_roc.pdf")
  if(!(opts$part=="1/1" & opts$valid=="none")){
    pdf.dir4=paste0(prefix, "_test_boxplot.pdf") 
    pdf.dir5=paste0(prefix, "_test_roc.pdf")
  }
  
  pdf(pdf.dir1) 
  matplot(result[[1]]$n.var, error.cv, type = "l", log = "x", col = rep(1, cvt),
          xlab = "Number of vars",  
          ylab = "CV Error", lty = 1) 
  lines(result[[1]]$n.var, error.cv.rm, lwd = 2)
  text(x = marker.num,y = min(error.cv),labels = marker.num,adj = c(1, 0.5),cex = 1.2, xpd = TRUE, font = 2, col = "pink")
  abline(v = marker.num, col = "pink", lwd = 2)
  dev.off()
  
  # pick marker by corossvalidation 
  marker.t <- table(unlist(lapply(result, function(x) { 
    lapply(x$res, "[", 1:marker.num) 
  }))) 
  marker.t <- sort(marker.t, d = T) 
  names(marker.t) <- colnames(Data)[as.numeric(names(marker.t))] 
  marker.dir <- paste0(prefix, "_marker.txt") 
  write.table(marker.t, marker.dir, col.names = F, sep = "\t", quote = F) 
  marker.p <- names(marker.t)[1:marker.num] 
  
  # 2/3 for train
  set.seed(0) # set.seed(10) # set.seed(100)
  train.rf=randomForest(Data[ind==1,marker.p,drop=FALSE],Data[ind==1,"group"],ntree=1000,proximity=TRUE,importance=TRUE)
  train.p <- predict(train.rf, type = "prob")
  pdf(pdf.dir2)
  if (c != "none"){
    boxplot(train.p[, 2] ~ Data[ind==1,"group"], col = mycol, main = "Probability", 
            names = levels(Data[,"group"]), xlab = "Groups", ylab = "POD")
  }else{
    boxplot(train.p[, 2] ~ Data[ind==1,"group"], col = 3:2, main = "Probability", 
            names = levels(Data[,"group"]), xlab = "Groups", ylab = "POD") 
  }
  dev.off() 
  pr.dir <- paste0(prefix, "_train_probability.txt") 
  write.table(train.p[, 2], pr.dir, sep = "\t", quote = F, col.names = F)
  
  # varImPlot::MeanDecreaseAccuracy(å¹³å‡å‡å°‘åº¦ï¼‰,å³æ²¡æœ‰è¿™ä¸ªFeatureï¼Œåˆ†ç±»å‡†ç¡®åº¦ä¸‹é™çš„ç¨‹åº¦ï¼Œç›¸å½“äºŽå¸¸ç”¨çš„åˆ†ç±»è´¡çŒ®åº¦çš„æ¦‚å¿µã€‚
  #pdf(prefix, "_train_varImPlot.pdf")
  #varImpPlot(train.rf)
  #dev.off()
  
  if (length(marker.p) > 1){
    bbtheme <- function(){
      return(theme_bw() + 
               theme(plot.title=element_text(size=rel(1),hjust=0.5),
                     plot.margin = unit(c(1, 1, 1, 1), "lines"),
                     panel.grid.major=element_line(color="white"),
                     panel.grid.minor=element_line(color="white"),
                     axis.title=element_text(size=rel(1))#,
                     #axis.text.x=element_text(angle=30,hjust =1),
                     #legend.title=element_blank(),
                     #legend.text = element_text(size = 6),
                     #legend.key.size = unit(.4,'cm'),
                     #legend.spacing.x = unit(0.1, 'cm')
               ))
    }
    
    GGplot <- function(x,yb){
      ggplot(x,aes(x=OTU_genus,y=ID)) +
        geom_point(size=3,color="#FF4500") +
        #geom_segment(aes(x=OTU_genus,xend=OTU_genus,y=0,yend=MeanDecreaseAccuracy))+
        #bbtheme() + #geom_line() +
        theme_bw() +
        coord_flip() + xlab('') + ylab(yb) +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.y = element_line(colour="grey60", linetype="dashed"),
              axis.text.y = element_text(hjust = 0))
    }
    
    WHat <- varImpPlot(train.rf) %>% as.data.frame %>% 
      rownames_to_column() %>% as_tibble %>%
      rename_at(1,~"OTU_genus") %>% #rename(OTU_genus = colnames(.)[1]) %>%
      gather(value,ID,-OTU_genus) %>%
      group_by(value) %>% nest %>% 
      mutate(data = map(data,function(x) 
        x %>% arrange(ID) %>%
          mutate(OTU_genus = fct_inorder(OTU_genus))))
    
    Pt <- map(seq_along(WHat),function(x) GGplot(WHat$data[[x]],WHat$value[[x]]))
    
    ggsave(str_c(prefix,"_varImpPlot.pdf"),Pt[[1]] + Pt[[2]],
           width = 12,
           height = log2(nrow(WHat$data[[1]]))*1.2+3
    )
  }
  
  Data %>% dplyr::select(-group) %>% rownames_to_column() %>% rename_at(1,~"SampleID") %>% #20230110
    gather(OTUID,value,-SampleID) %>%
    mutate_if(is.character,~fct_inorder(.x)) %>%
    spread(SampleID,value) %>%
    rename_at(1,~"OTU ID") %>% #rename(`OTU ID` = colnames(.)[1]) %>%
    filter(`OTU ID` %in% marker.p) %>%
    write_tsv(str_c("marker.p.",opts$input))
  
  # train ROC
  pdf(pdf.dir3) 
  roc <- roc(Data[ind==1,"group"],train.p[, 2],plot=T,col="black",ci=F,auc.polygon=F,print.thres=F,print.auc=F,percent=F,xlab="Specificity(%)",ylab="Sensitivity(%)")
  
  U <- wilcox.test(train.p[Data[ind==1,"group"] %in% levels(Data[ind==1,"group"])[1],2],
                   train.p[Data[ind==1,"group"] %in% levels(Data[ind==1,"group"])[2],2])
  U %>% tidy %>% write.csv(file = "wilcox.train.ROC.txt")
  
  sens.ci <- ci.se(roc)
  plot(sens.ci, type="s", col="lightblue",border="white")

  legend("bottomright",bty="n",paste0("AUC: ",round(ci(roc)[2],4),"\n","95% CI: ",round(ci(roc)[1],4),"-",round(ci(roc)[3],4),"\n","p-value: ",Minus(U$p.value,4),"\n"," "))
  tt <- paste0(round(coords(roc,"best")[1,1],4),"\n","(",round(coords(roc,"best")[1,2],4),",",round(coords(roc,"best")[1,3],4),")")
  points(coords(roc,"best")[1,2],coords(roc,"best")[1,3],pch=16,col="red",cex=1.5,font=1)
  text(coords(roc,"best")[1,2]-0.1,coords(roc,"best")[1,3]-0.1,tt,cex=1.5,pos=4,col="black")
  dev.off()
  aoteman_train_AUC<-round(ci(roc)[2],4)
  aoteman_train_pvalue<-Minus(U$p.value,4)



  save.image(file = "projectimage.RData")
}

if (opts$part!="1/1" & opts$gp == "none"){
  # test predict 
  # 1/3
  test.p <- predict(train.rf, Data[ind==2,-ncol(Data)], type = "prob") 
  pr.dir <- paste0(prefix, "_test_probability.txt") 
  write.table(test.p[,2], pr.dir, sep = "\t", quote = F, col.names = F)
  #test.result=cbind(Data[ind==2,"group"],test.p[,2]) 
  #write.table(test.result, pr.dir, sep = "\t", quote = F, col.names = F)
  
  pdf(pdf.dir4) 
  if (c != "none"){
    boxplot(test.p[, 2] ~ Data[ind==2,"group"], col = mycol, main = "Probability", 
            names = levels(Data[,"group"]), xlab = "Groups", ylab = "POD")
  }else{
    boxplot(test.p[, 2] ~ Data[ind==2,"group"], col = 3:2, main = "Probability", 
            names = levels(Data[,"group"]), xlab = "Groups", ylab = "POD") 
  }
  dev.off() 
  
  # test ROC 
  pdf(pdf.dir5)
  roc <- roc(Data[ind==2,"group"],test.p[, 2],plot=T,col="black",ci=F,auc.polygon=F,print.thres=F,print.auc=F,percent=F,xlab="Specificity(%)",ylab="Sensitivity(%)")
  
  U <- wilcox.test(test.p[Data[ind==2,"group"] %in% levels(Data[ind==2,"group"])[1],2],
                   test.p[Data[ind==2,"group"] %in% levels(Data[ind==2,"group"])[2],2])
  U %>% tidy %>% write.csv(file = "wilcox.test.ROC.txt")
  
  sens.ci <- ci.se(roc)
  plot(sens.ci, type="s", col="lightblue",border="white")
  legend("bottomright",bty="n",paste0("AUC: ",round(ci(roc)[2],4),"\n","95% CI: ",round(ci(roc)[1],4),"-",round(ci(roc)[3],4),"\n","p-value: ",Minus(U$p.value,4),"\n"," "))
  tt <- paste0(round(coords(roc,"best")[1,1],4),"\n","(",round(coords(roc,"best")[1,2],4),",",round(coords(roc,"best")[1,3],4),")")
  points(coords(roc,"best")[1,2],coords(roc,"best")[1,3],pch=16,col="red",cex=1.5,font=1)
  text(coords(roc,"best")[1,2]-0.1,coords(roc,"best")[1,3]-0.1,tt,cex=1.5,pos=4,col="black")
  dev.off()
  
  Data[ind==2,] %>% dplyr::select(-group) %>% rownames_to_column() %>% rename_at(1,~"SampleID") %>%
    gather(OTUID,value,-SampleID) %>%
    mutate_if(is.character,~fct_inorder(.x)) %>%
    spread(SampleID,value) %>%
    rename_at(1,~"OTU ID") %>% #rename(`OTU ID` = colnames(.)[1]) %>%
    filter(`OTU ID` %in% marker.p) %>%
    write_tsv(str_c("marker.p.test.",opts$input))
  
}else if(opts$gp != "none"){
  testgp <- read_tsv(opts$gp) %>% 
    rename_all(~c("SampleID","group")) %>% #rename(SampleID = colnames(.)[1],group = colnames(.)[2]) %>%
    mutate(group = factor(group,levels=unique(data1$group)))
  
  testo <- marker.p %>% enframe %>% 
    mutate(backup = value) %>%
    separate(backup,into = c("SampleID","genus"),sep = " ") %>%
    .[,c("SampleID","genus")]
  
  test <- read_tsv(opts$valid) %>% 
    rename_at(1,~"SampleID") %>% #rename("SampleID" = colnames(.)[1]) %>%
    .[,c("SampleID",as.vector(testgp$SampleID))] %>% 
    mutate_at(vars(2:ncol(.)),function(x) x/sum(x)) %>%
    mutate(SampleID = map_chr(SampleID,~str_split(.x," ")[[1]][1])) %>%
    right_join(testo) %>% 
    replace(., is.na(.), 0) %>% 
    mutate(SampleID = str_c(SampleID,genus,sep = " ")) %>%
    mutate(SampleID = factor(SampleID, 
                             levels=unique(SampleID))) %>%
    dplyr::select(-genus) %>%
    gather(var, value, -SampleID) %>% 
    mutate(var = factor(var,levels=unique(var))) %>%
    spread(SampleID, value) %>% 
    rename_at(1,~"SampleID") %>% #rename(SampleID = var) %>%
    right_join(.,testgp)
  
  test$group <- factor(test$group,levels = levels(Data$group))
  write_tsv(test %>% dplyr::select(-group),"test_output.xls")
  cc <- test$SampleID
  Test <- test %>% as.data.frame %>% .[,-1]
  rownames(Test) <- cc
  
  # test predict 
  # 1/3
  if(ncol(Test)!=2){ 
    test.p <- predict(train.rf, Test[,-ncol(Test)], type = "prob")
    pr.dir <- paste0(prefix, "_test_probability.txt") 
    write.table(test.p[,2], pr.dir, sep = "\t", quote = F, col.names = F)
  }else{
    TTest <- data.frame(Test[,-ncol(Test)])
    rownames(TTest) <- rownames(Test)
    colnames(TTest) <- colnames(Test)[1]
    test.p <- predict(train.rf, TTest , type = "prob")
    pr.dir <- paste0(prefix, "_test_probability.txt") 
    write.table(test.p[,2], pr.dir, sep = "\t", quote = F, col.names = F)
  }
  
  pdf(pdf.dir4) 
  if (c != "none"){
    boxplot(test.p[, 2] ~ Test[,"group"], col = mycol, main = "Probability", 
            names = levels(Data[,"group"]), xlab = "Groups", ylab = "POD")
  }else{
    boxplot(test.p[, 2] ~ Test[,"group"], col = 3:2, main = "Probability", 
            names = levels(Data[,"group"]), xlab = "Groups", ylab = "POD") 
  }
  dev.off() 
  
  # test ROC 
  pdf(pdf.dir5)
  roc <- roc(Test[,"group"],test.p[, 2],plot=T,col="black",ci=F,auc.polygon=F,
             print.thres=F,print.auc=F,percent=F,xlab="Specificity(%)",ylab="Sensitivity(%)")
  
  U <- wilcox.test(test.p[Test[,"group"] %in% levels(Test[,"group"])[1],2],
                   test.p[Test[,"group"] %in% levels(Test[,"group"])[2],2])
  U %>% tidy %>% write.csv(file = "wilcox.test.ROC.txt")
  
  sens.ci <- ci.se(roc)
  plot(sens.ci, type="s", col="lightblue",border="white")
  # æ ‡æ³¨
  legend("bottomright",bty="n",paste0("AUC: ",round(ci(roc)[2],4),"\n","95% CI: ",round(ci(roc)[1],4),"-",round(ci(roc)[3],4),"\n","p-value: ",Minus(U$p.value,4),"\n"," "))
  tt <- paste0(round(coords(roc,"best")[1,1],4),"\n","(",round(coords(roc,"best")[1,2],4),",",round(coords(roc,"best")[1,3],4),")")
  points(coords(roc,"best")[1,2],coords(roc,"best")[1,3],pch=16,col="red",cex=1.5,font=1)
  text(coords(roc,"best")[1,2]-0.1,coords(roc,"best")[1,3]-0.1,tt,cex=1.5,pos=4,col="black")
  
  dev.off()
  
  test %>% dplyr::select(-group) %>% gather(OTUID,value,-SampleID) %>%
    mutate_if(is.character,~fct_inorder(.x)) %>%
    spread(SampleID,value) %>%
    rename_at(1,~"OTU ID") %>% #rename(`OTU ID` = colnames(.)[1]) %>%
    filter(`OTU ID` %in% marker.p) %>%
    write_tsv(str_c("marker.p.",opts$valid))
}

if(!(opts$part=="1/1" & opts$valid=="none")){
  aoteman_test_AUC<-round(ci(roc)[2],4)
  aoteman_test_pvalue<-Minus(U$p.value,4)
}else{
  aoteman_test_AUC<-NA
  aoteman_test_pvalue<-NA
}
   
##################################################################################################
write_tsv(opts %>% as.matrix %>% as.data.frame %>% rownames_to_column() %>% 
            as_tibble %>% mutate(V1 = as.character(V1)),
          str_c("Parameter",
                str_replace_all(as.character(date())," ","_") %>% str_replace_all(":","_"),
                ".xls"),
          col_names = FALSE)

if(!file.exists("Parameter.xls")){
  newopts<-opts %>% as.matrix %>% as.data.frame %>% rownames_to_column() %>% 
    as_tibble %>% mutate(V1 = as.character(V1)) %>% rbind(c("train_AUC",aoteman_train_AUC),c("train_pvalue",aoteman_train_pvalue),c("test_AUC",aoteman_test_AUC),c("test_pvalue",aoteman_test_pvalue),c("marker.num",marker.num),c("Rscript","rfcv_subsample.R"))
  
  write_tsv(newopts %>% as.matrix %>% as.data.frame %>% 
              as_tibble %>% mutate(V1 = as.character(V1)),
            str_c("Parameter", ".xls"),
            col_names = FALSE)
}else{
  oldopts<-read_tsv("Parameter.xls",col_names = F)
  
  newopts<-opts %>% as.matrix %>% as.data.frame %>% rownames_to_column() %>% 
    as_tibble %>% mutate(V1 = as.character(V1)) %>% rbind(c("train_AUC",aoteman_train_AUC),c("train_pvalue",aoteman_train_pvalue),c("test_AUC",aoteman_test_AUC),c("test_pvalue",aoteman_test_pvalue),c("marker.num",marker.num),c("Rscript","rfcv_subsample.R"))
  
  nnewopts<-cbind(oldopts,newopts[,2])
  write_tsv(nnewopts %>% as.matrix %>% as.data.frame %>% 
              as_tibble %>% mutate(V1 = as.character(V1)),
            str_c("Parameter", ".xls"),
            col_names = FALSE)
}

```

