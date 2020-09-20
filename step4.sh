a=`cat map-group.txt | wc -l`
if [ $a -ge 150 ]
then
    export Community_b_lcex="0.5"
    export Community_b_ncol="6"
    export Community_b_w="30"
    export Community_b_bo="5"
    export Heatmap_w="30"
    export Heatmap_h="10"
    export Heatmap_keyh="14"
    export Heatmap_marble="4-0-0-20"
    export Hclust_bar_h="30"
    export Hcluster_tree_h="30"
elif [ $a -ge 100 ]
then
    export Community_b_lcex="0.5"
    export Community_b_ncol="6"
    export Community_b_w="20"
    export Community_b_bo="5"
    export Heatmap_w="20"
    export Heatmap_h="10"
    export Heatmap_keyh="14"
    export Heatmap_marble="4-0-0-14"
    export Hclust_bar_h="20"
    export Hcluster_tree_h="20"
elif [ $a -ge 75 ]
then
    export Community_b_lcex="0.5"
    export Community_b_ncol="6"
    export Community_b_w="16"
    export Community_b_bo="5"
    export Heatmap_w="16"
    export Heatmap_h="10"
    export Heatmap_keyh="13"
    export Heatmap_marble="4-0-0-14"
    export Hclust_bar_h="15"
    export Hcluster_tree_h="15"
elif [ $a -ge 50 ]
then
    export Community_b_lcex="0.5"
    export Community_b_ncol="5"
    export Community_b_w="12"
    export Community_b_bo="5"
    export Heatmap_w="12"
    export Heatmap_h="10"
    export Heatmap_keyh="12"
    export Heatmap_marble="4-0-0-14"
    export Hclust_bar_h="10"
    export Hcluster_tree_h="10"
elif [ $a -ge 25 ]
then
    export Community_b_lcex="0.5"
    export Community_b_ncol="5"
    export Community_b_w="10"
    export Community_b_bo="5"
    export Heatmap_w="10"
    export Heatmap_h="10"
    export Heatmap_keyh="10"
    export Heatmap_marble="4-0-0-10"
    export Hclust_bar_h="8"
    export Hcluster_tree_h="8"
elif [ $a -gt 0 ]
then
    export Community_b_lcex="0.5"
    export Community_b_ncol="5"
    export Community_b_w="8"
    export Community_b_bo="4"
    export Heatmap_w="8"
    export Heatmap_h="9"
    export Heatmap_keyh="9"
    export Heatmap_marble="4-0-0-8"
    export Hclust_bar_h="5"
    export Hcluster_tree_h="7"
fi

# need order; map-group.txt groups.txt File;
# run directory: results
cat map-group.txt|awk '{print $1"\t"$1}' >order
cp order groups.txt
cp /work/scripts/16s/less_than_100/head.txt ./
cat head.txt map-group.txt |awk -F" " '{print $1"\t"$2}' >map-group.txt2
mv map-group.txt2 map-group.txt
cat map-group.txt |sed '1d'|awk -F"\t" '{if($2==last){group=group","$1}else{print group;group=$2":"$1};last=$2}END{print group}'|sed '1d' >map.txt

######################################### Lefse #################################
source /root/miniconda3/bin/activate
source activate qiime1

mkdir Lefse
cd Lefse/
mkdir Lefse_genus Lefse_species
cd Lefse_species
cp ../../OTU_Taxa/rarefac.otu_taxa_table.xls ./
sed -i 's/ //g' rarefac.otu_taxa_table.xls
python /work/scripts/16s/tax_split.py -i rarefac.otu_taxa_table.xls
#准备数据
python /work/scripts/16s/make_lefse_data.py -i rarefac.otu_taxa_split_index_table.xls -g ../../map-group.txt
sed -i '2d' lefse.data.xls
#处理数据分析
python2 /work/scripts/16s/format_input.py lefse.data.xls data.in -c 1 -o 1000000
python2 /work/scripts/16s/run_lefse.py data.in data.res -l 2
#画图
python2 /work/scripts/16s/plot_res.py data.res lefse.pdf  --format pdf
python2 /work/scripts/16s/plot_cladogram.py data.res lefse.cladogram.pdf --format pdf
cat data.res|awk 'BEGIN{printf "Biomaker_names\tLogarithm value\tGroups\tLDA_value\tP_value\n"}{print $0}'>LDA.xls
#显示种
python2 /work/scripts/16s/plot_cladogram.py data.res lefse.cladogram_species.pdf --format pdf --labeled_start_lev 1 --labeled_stop_lev 7 --abrv_start_lev 1 --abrv_stop_lev 7

cd ../Lefse_genus
cp ../../OTU_Taxa/rarefac.otu_taxa_table.xls ./
sed -i 's/ //g' rarefac.otu_taxa_table.xls
sed -i 's/;s__.*//g' rarefac.otu_taxa_table.xls
python2 /work/scripts/16s/tax_split.py -i rarefac.otu_taxa_table.xls
#准备数据
python2 /work/scripts/16s/make_lefse_data.py -i rarefac.otu_taxa_split_index_table.xls -g ../../map-group.txt
sed -i '2d' lefse.data.xls
#处理数据分析
python2 /work/scripts/16s/format_input.py lefse.data.xls data.in -c 1 -o 1000000
python2 /work/scripts/16s/run_lefse.py data.in data.res -l 2
#画图
python2 /work/scripts/16s/plot_res.py data.res lefse.pdf  --format pdf
python2 /work/scripts/16s/plot_cladogram.py data.res lefse.cladogram.pdf --format pdf
cat data.res|awk 'BEGIN{printf "Biomaker_names\tLogarithm value\tGroups\tLDA_value\tP_value\n"}{print $0}'>LDA.xls
cd ../../

#Random_Forest
mkdir Random_Forest
cd Random_Forest/
ln -s ../OTU_Taxa/rarefac.otu_table.biom ./
ln -s ../OTU_Taxa/rarefac.otu_genus.xls ./
ln -s ../map-group.txt ./
ln -s ../order ./
random_forest4key_out_select.pl  rarefac.otu_table.biom rarefac.otu_genus.xls  map-group.txt 0.001
sortSample4otu_table.pl OTU-extract.all.xls order  >OTU-extract.all.sort.xls
cd ../

source deactivate

#Alpha_rarefac
mkdir Alpha_rarefac
cp ../process/otu_0.97/alpha_rarefac/* Alpha_rarefac -r
cd Alpha_rarefac/
#ln -s ../Shannon_rarefac/*.r_shannon ./
#ln -s ../Rarefactions/*.rarefaction ./
#ln -s ../Estimators/*.summary ./
#cp -r ../map.txt ./
#plot-rarefaction.pl -d . -a l -l 0.97 -g map.txt
#plot-rarefaction.pl -d . -a l -l 0.97 -g map.txt -m shannon
cp ../map-group.txt ./
cp ../color.txt ./
python /work/scripts/16s/plot_alpha_point.py -m rarefaction -g map-group.txt -gc color.txt
python /work/scripts/16s/plot_alpha_point.py -m r_shannon -g map-group.txt -gc color.txt
mv rarefaction.pdf ../Rarefactions 
mv r_shannon.pdf ../Shannon_rarefac

cp ../map-group.txt ./
cp ../color.txt ./
alpha_rarefac.pl map-group.txt  >alpha_rarefac.summary.xls
plot-alpha_rarefac.pl alpha_rarefac.summary.xls map-group.txt -col color.txt

#Wilcox-alpha
cat alpha_rarefac.summary.xls|awk '{print $1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$9}' > alpha_rarefac.summary.sub.xls
transposition.pl alpha_rarefac.summary.sub.xls > alpha.xls
ln -s ../Random_Forest/map.*.txt ./
ls map.*-*.txt|awk -F "." '{print "wilcox_nounif.py -i alpha.xls -c "$2" -g map."$2".txt -o wilcox.alpha."$2".xls"}' >wilcox.sh
sh wilcox.sh 
rm *.filtered.xls
sed -i 's/"//g' *.xls
cat wilcox.alpha.*.xls|sed '/mean/{x;p;x;}' > wilcox.alpha.xls
rm wilcox.alpha.*.xls  map-group.txt *sh *.r mothur* alpha_rarefac.summary.sub.xls map.*.txt  otus.* alpha.xls
cd ../

#Community
cd Community/
rm *
ln -s ../order ./sample
sub_Sample4otu_table.pl ../../results/Community/phylum.xls sample > phylum.xls
sub_Sample4otu_table.pl ../../results/Community/class.xls sample > class.xls
sub_Sample4otu_table.pl ../../results/Community/order.xls sample > order.xls
sub_Sample4otu_table.pl ../../results/Community/family.xls sample > family.xls
sub_Sample4otu_table.pl ../../results/Community/genus.xls sample  > genus.xls

sub_Sample4otu_table.pl ../../results/Community/phylum.percents.xls sample > phylum.percents.xls
sub_Sample4otu_table.pl ../../results/Community/class.percents.xls sample > class.percents.xls
sub_Sample4otu_table.pl ../../results/Community/order.percents.xls sample > order.percents.xls
sub_Sample4otu_table.pl ../../results/Community/family.percents.xls sample > family.percents.xls
sub_Sample4otu_table.pl ../../results/Community/genus.percents.xls sample  > genus.percents.xls

sed -i 's/ //g' *.xls

mkdir Community_barplot
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -i phylum.xls  -p 0
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -i class.xls  -p 0.01
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -i order.xls  -p 0.01
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -i family.xls  -p 0.01
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -i genus.xls  -p 0.01
mv *.pdf ./Community_barplot/

ln -s ../order ./sample
cp ../map-group.txt ./
cp ../color.txt ./
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -col color.txt -ave 2 -i phylum.xls  -p 0 -map map-group.txt
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -col color.txt -ave 2 -i class.xls  -p 0.01  -map map-group.txt
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -col color.txt -ave 2 -i order.xls  -p 0.01  -map map-group.txt
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -col color.txt -ave 2 -i family.xls  -p 0.01  -map map-group.txt
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -col color.txt -ave 2 -i genus.xls  -p 0.01  -map map-group.txt
mkdir Community_barplot_groups/
mv *.pdf ./Community_barplot_groups/

mkdir Community_average
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -ave 1 -i phylum.xls  -p 0 -map map-group.txt
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -ave 1 -i class.xls  -p 0.01  -map map-group.txt
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -ave 1 -i order.xls  -p 0.01  -map map-group.txt
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -ave 1 -i family.xls  -p 0.01  -map map-group.txt
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -ave 1 -i genus.xls  -p 0.01  -map map-group.txt
mv *.pdf ./Community_average

mkdir Community_bubble
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -bb 2 -i phylum.xls  -p 0 
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -bb 2 -i class.xls  -p 0.01  
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -bb 2 -i order.xls  -p 0.01  
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -bb 2 -i family.xls  -p 0.01 
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -bb 2 -i genus.xls  -p 0.01  
mv *.pdf ./Community_bubble

mkdir Community_bubble_average
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -ave 1 -bb 2 -i phylum.xls  -p 0 -map map-group.txt
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -ave 1 -bb 2 -i class.xls  -p 0.01  -map map-group.txt
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -ave 1 -bb 2 -i order.xls  -p 0.01  -map map-group.txt
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -ave 1 -bb 2 -i family.xls  -p 0.01  -map map-group.txt
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py -ave 1 -bb 2 -i genus.xls  -p 0.01  -map map-group.txt
mv *.pdf ./Community_bubble_average

mkdir Community_test_barplot
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py  -col color.txt -oth F -pt 1 -i phylum.xls  -p 0 -map map-group.txt
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py  -col color.txt -oth F -pt 1 -i class.xls  -p 0.01  -map map-group.txt
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py  -col color.txt -oth F -pt 1 -i order.xls  -p 0.01  -map map-group.txt
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py  -col color.txt -oth F -pt 1 -i family.xls  -p 0.01  -map map-group.txt
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py  -col color.txt -oth F -pt 1 -i genus.xls  -p 0.01  -map map-group.txt
mv *.pdf ./Community_test_barplot

mkdir Community_test_boxplot
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py  -col color.txt -oth F -pt 2 -i phylum.xls  -p 0 -map map-group.txt
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py  -col color.txt -oth F -pt 2 -i class.xls  -p 0.01  -map map-group.txt
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py  -col color.txt -oth F -pt 2 -i order.xls  -p 0.01  -map map-group.txt
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py  -col color.txt -oth F -pt 2 -i family.xls  -p 0.01  -map map-group.txt
python /work/users/chaoliu/scripts/ggplot2_bar_box_bubble.py  -col color.txt -oth F -pt 2 -i genus.xls  -p 0.01  -map map-group.txt
mv *.pdf ./Community_test_boxplot

#mkdir Lefse
#mkdir Lefse/otu Lefse/genus Lefse/phylum Lefse/class Lefse/order Lefse/family
#sub_Sample4otu_table.pl ../../process/otu_0.97/rarefac.otu_genus.xls ../groups.txt  > rarefac.otu_genus.xls
#sed -i 's/(/_/' rarefac.otu_genus.xls
#sed -i 's/)//' rarefac.otu_genus.xls
#plot-lda.pl -i rarefac.otu_genus.xls -o Lefse/otu -m ../map-group.txt -g group
#plot-lda.pl -i genus.xls -o Lefse/genus -m ../map-group.txt -g group
#plot-lda.pl -i phylum.xls -o Lefse/phylum -m ../map-group.txt -g group
#plot-lda.pl -i class.xls -o Lefse/class -m ../map-group.txt -g group
#plot-lda.pl -i order.xls -o Lefse/order -m ../map-group.txt -g group
#plot-lda.pl -i family.xls -o Lefse/family -m ../map-group.txt -g group

#Community_boxplot
mkdir Community_boxplot
cd Community_boxplot
ln -s ../../map-group.txt ./
ln -s ../../color.txt
sed -i 's/\#//g' map-group.txt
transposition.pl ../phylum.percents.xls >phylum.percents.xls
transposition.pl ../genus.percents.xls >genus.percents.xls
transposition.pl ../family.percents.xls >family.percents.xls
transposition.pl ../order.percents.xls >order.percents.xls
transposition.pl ../class.percents.xls >class.percents.xls 
mkdir phylum genus family class order
Rscript /work/users/chaoliu/scripts/make_community_boxplot.r -i phylum.percents.xls -m map-group.txt -c color.txt -l T
mv *.pdf phylum
Rscript /work/users/chaoliu/scripts/make_community_boxplot.r -i genus.percents.xls -m map-group.txt -c color.txt -l T
mv *.pdf genus
Rscript /work/users/chaoliu/scripts/make_community_boxplot.r -i order.percents.xls -m map-group.txt -c color.txt -l T
mv *.pdf order
Rscript /work/users/chaoliu/scripts/make_community_boxplot.r -i family.percents.xls -m map-group.txt -c color.txt -l T
mv *.pdf family
Rscript /work/users/chaoliu/scripts/make_community_boxplot.r -i class.percents.xls -m map-group.txt -c color.txt -l T
mv *.pdf  class
rm *
cd ../../

mkdir PCoA
rm -rf Pcoa
cd PCoA
rm *
ln -s ../map-group.txt ./
ln -s ../Beta_diversity/bray_curtis_dm.txt ./
ln -s ../Beta_diversity/weighted_unifrac_dm.txt ./
ln -s ../Beta_diversity/unweighted_unifrac_dm.txt ./
ln -s ../color.txt ./
#noname
mkdir noname
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -i bray_curtis_dm.txt -md PCoA -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -i bray_curtis_dm.txt -md PCoA -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -i bray_curtis_dm.txt -md PCoA -pc 2-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -i weighted_unifrac_dm.txt -md PCoA -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -i weighted_unifrac_dm.txt -md PCoA -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -i weighted_unifrac_dm.txt -md PCoA -pc 2-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -i unweighted_unifrac_dm.txt -md PCoA -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -i unweighted_unifrac_dm.txt -md PCoA -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -i unweighted_unifrac_dm.txt -md PCoA -pc 2-3 -map map-group.txt
mv *.pdf *.xls noname/

#name
mkdir name
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -lab T -i bray_curtis_dm.txt -md PCoA -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -lab T -i bray_curtis_dm.txt -md PCoA -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -lab T -i bray_curtis_dm.txt -md PCoA -pc 2-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -lab T -i weighted_unifrac_dm.txt -md PCoA -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -lab T -i weighted_unifrac_dm.txt -md PCoA -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -lab T -i weighted_unifrac_dm.txt -md PCoA -pc 2-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -lab T -i unweighted_unifrac_dm.txt -md PCoA -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -lab T -i unweighted_unifrac_dm.txt -md PCoA -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -lab T -i unweighted_unifrac_dm.txt -md PCoA -pc 2-3 -map map-group.txt
mv *.pdf *.xls name/

#ellipse
mkdir ellipse
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -e T -i bray_curtis_dm.txt -md PCoA -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -e T -i bray_curtis_dm.txt -md PCoA -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -e T -i bray_curtis_dm.txt -md PCoA -pc 2-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -e T -i weighted_unifrac_dm.txt -md PCoA -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -e T -i weighted_unifrac_dm.txt -md PCoA -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -e T -i weighted_unifrac_dm.txt -md PCoA -pc 2-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -e T -i unweighted_unifrac_dm.txt -md PCoA -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -e T -i unweighted_unifrac_dm.txt -md PCoA -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -e T -i unweighted_unifrac_dm.txt -md PCoA -pc 2-3 -map map-group.txt
mv *.pdf *.xls ellipse/

#box
mkdir box
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -bx T -e T -i bray_curtis_dm.txt -md PCoA -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -bx T -e T -i bray_curtis_dm.txt -md PCoA -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -bx T -e T -i bray_curtis_dm.txt -md PCoA -pc 2-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -bx T -e T -i weighted_unifrac_dm.txt -md PCoA -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -bx T -e T -i weighted_unifrac_dm.txt -md PCoA -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -bx T -e T -i weighted_unifrac_dm.txt -md PCoA -pc 2-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -bx T -e T -i unweighted_unifrac_dm.txt -md PCoA -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -bx T -e T -i unweighted_unifrac_dm.txt -md PCoA -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -bx T -e T -i unweighted_unifrac_dm.txt -md PCoA -pc 2-3 -map map-group.txt
mv *.pdf *.xls box/

cd ../
############################################################ PCA ##########################################################
mkdir PCA
cd PCA
rm -rf pca
rm *
ln -s ../map-group.txt ./
ln -s ../../process/otu_0.97/rarefac.otu_table.xls ./
ln -s ../color.txt ./

#noname
mkdir noname
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -i rarefac.otu_table.xls -md PCA -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -i rarefac.otu_table.xls -md PCA -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -i rarefac.otu_table.xls -md PCA -pc 2-3 -map map-group.txt
mv PCA* noname/

#name
mkdir name
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -lab T -i rarefac.otu_table.xls -md PCA -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -lab T -i rarefac.otu_table.xls -md PCA -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -lab T -i rarefac.otu_table.xls -md PCA -pc 2-3 -map map-group.txt
mv PCA* name/

#ellipse
mkdir ellipse
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -e T -i rarefac.otu_table.xls -md PCA -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -e T -i rarefac.otu_table.xls -md PCA -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -e T -i rarefac.otu_table.xls -md PCA -pc 2-3 -map map-group.txt
mv PCA* ellipse/

#box
mkdir box
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -bx T -e T -i rarefac.otu_table.xls -md PCA -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -bx T -e T -i rarefac.otu_table.xls -md PCA -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -bx T -e T -i rarefac.otu_table.xls -md PCA -pc 2-3 -map map-group.txt
mv PCA*  box/

cd ../

############################################################ NMDS ##########################################################
mkdir NMDS
rm -rf nmds
cd NMDS
rm *
ln -s ../color.txt ./
ln -s ../map-group.txt ./
ln -s ../Beta_diversity/bray_curtis_dm.txt ./
ln -s ../Beta_diversity/weighted_unifrac_dm.txt ./
ln -s ../Beta_diversity/unweighted_unifrac_dm.txt ./

#noname
mkdir noname
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -i bray_curtis_dm.txt -md NMDS -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -i bray_curtis_dm.txt -md NMDS -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -i bray_curtis_dm.txt -md NMDS -pc 2-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -i weighted_unifrac_dm.txt -md NMDS -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -i weighted_unifrac_dm.txt -md NMDS -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -i weighted_unifrac_dm.txt -md NMDS -pc 2-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -i unweighted_unifrac_dm.txt -md NMDS -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -i unweighted_unifrac_dm.txt -md NMDS -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -i unweighted_unifrac_dm.txt -md NMDS -pc 2-3 -map map-group.txt
mv *.pdf *.xls noname/

#name
mkdir name
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -lab T -i bray_curtis_dm.txt -md NMDS -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -lab T -i bray_curtis_dm.txt -md NMDS -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -lab T -i bray_curtis_dm.txt -md NMDS -pc 2-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -lab T -i weighted_unifrac_dm.txt -md NMDS -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -lab T -i weighted_unifrac_dm.txt -md NMDS -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -lab T -i weighted_unifrac_dm.txt -md NMDS -pc 2-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -lab T -i unweighted_unifrac_dm.txt -md NMDS -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -lab T -i unweighted_unifrac_dm.txt -md NMDS -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -lab T -i unweighted_unifrac_dm.txt -md NMDS -pc 2-3 -map map-group.txt
mv *.pdf *.xls name/

#ellipse
mkdir ellipse
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -e T -i bray_curtis_dm.txt -md NMDS -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -e T -i bray_curtis_dm.txt -md NMDS -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -e T -i bray_curtis_dm.txt -md NMDS -pc 2-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -e T -i weighted_unifrac_dm.txt -md NMDS -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -e T -i weighted_unifrac_dm.txt -md NMDS -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -e T -i weighted_unifrac_dm.txt -md NMDS -pc 2-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -e T -i unweighted_unifrac_dm.txt -md NMDS -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -e T -i unweighted_unifrac_dm.txt -md NMDS -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -e T -i unweighted_unifrac_dm.txt -md NMDS -pc 2-3 -map map-group.txt
mv *.pdf *.xls ellipse/

#box
mkdir box
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -bx T -e T -i bray_curtis_dm.txt -md NMDS -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -bx T -e T -i bray_curtis_dm.txt -md NMDS -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -bx T -e T -i bray_curtis_dm.txt -md NMDS -pc 2-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -bx T -e T -i weighted_unifrac_dm.txt -md NMDS -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -bx T -e T -i weighted_unifrac_dm.txt -md NMDS -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -bx T -e T -i weighted_unifrac_dm.txt -md NMDS -pc 2-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -bx T -e T -i unweighted_unifrac_dm.txt -md NMDS -pc 1-2 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -bx T -e T -i unweighted_unifrac_dm.txt -md NMDS -pc 1-3 -map map-group.txt
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -col color.txt -bx T -e T -i unweighted_unifrac_dm.txt -md NMDS -pc 2-3 -map map-group.txt
mv *.pdf *.xls box/

cd ../


#Heatmap
cd Heatmap
ln -s ../map-group.txt ./
cp ../color.txt ./
plot-heatmap.pl -color color.txt -i rarefac.otu_genus.xls -o heatmap.otu.top50.pdf -rtop 50 -ct 0 -slas 2 -rlc 0.7 -clc 0.7 -w $Heatmap_w -h $Heatmap_h -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt

ln -s ../Random_Forest/OTU-extract.all.sort.xls ./
bar_pie.pl -i OTU-extract.all.sort.xls -pie F -p 0.01
sed -i '/Others/d' ALL.new.OTU-extract.all.sort.xls
mv ALL.new.OTU-extract.all.sort.xls ALL.new.OTU-extract.all.sort.percent001.xls
bar_pie.pl -i OTU-extract.all.sort.xls -pie F -p 0.03
sed -i '/Others/d' ALL.new.OTU-extract.all.sort.xls
mv ALL.new.OTU-extract.all.sort.xls ALL.new.OTU-extract.all.sort.percent003.xls
plot-heatmap.pl -color color.txt -i OTU-extract.all.sort.xls -o heatmap.keyotu.pdf  -ct 0 -slas 2 -rlc 0.7 -clc 0.7 -w $Heatmap_w -h $Heatmap_keyh -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt
plot-heatmap.pl -color color.txt -i ALL.new.OTU-extract.all.sort.percent001.xls -o heatmap.keyotu.percent001.pdf  -ct 0 -slas 2 -rlc 0.7 -clc 0.7 -w $Heatmap_w -h $Heatmap_h -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt
plot-heatmap.pl -color color.txt -i ALL.new.OTU-extract.all.sort.percent003.xls -o    -ct 0 -slas 2 -rlc 0.7 -clc 0.7 -w $Heatmap_w -h $Heatmap_h -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt


#mkdir pheatmap
#ln -s ../map-group.txt ./
#cp ../../process/otu_0.97/rarefac.otu_genus.xls ./
#python /work/users/chaoliu/python_test/20180108_热图/pheatmap.py -i rarefac.otu_genus.xls -cd map-group.txt -o pheatmap.otu.top50.pdf -rtop 50 -clust_c F
#
#ln -s ../Random_Forest/OTU-extract.all.sort.xls ./
#bar_pie.pl -i OTU-extract.all.sort.xls -pie F -p 0.01
#sed -i '/Others/d' ALL.new.OTU-extract.all.sort.xls
#mv ALL.new.OTU-extract.all.sort.xls ALL.new.OTU-extract.all.sort.percent001.xls
#bar_pie.pl -i OTU-extract.all.sort.xls -pie F -p 0.03
#sed -i '/Others/d' ALL.new.OTU-extract.all.sort.xls
#mv ALL.new.OTU-extract.all.sort.xls ALL.new.OTU-extract.all.sort.percent003.xls
#
#python /work/users/chaoliu/python_test/20180108_热图/pheatmap.py -i OTU-extract.all.sort.xls -o pheatmap.keyotu.pdf -cd map-group.txt
#python /work/users/chaoliu/python_test/20180108_热图/pheatmap.py -i ALL.new.OTU-extract.all.sort.percent001.xls -o pheatmap.keyotu.percent001.pdf -cd map-group.txt
#python /work/users/chaoliu/python_test/20180108_热图/pheatmap.py -i ALL.new.OTU-extract.all.sort.percent003.xls -o pheatmap.keyotu.percent001.pdf -cd map-group.txt
#rm ALL.new.* map-group.txt OTU-extract.all.sort.xls dat.cor.* cmd.r bar.ALL.OTU-extract.all.sort.xls.pdf map-group.txt OTU-extract.all.sort.xls rarefac.otu_genus.xls
cd ../

#### treebar ###################################################
cd Hclust_bar/
plot-treebar.pl -otu ../../process/otu_0.97/rarefac.otu_table.xls -tax ../../process/otu_0.97/tax_summary_a/phylum.xls -o treebar_otu_table_phylum.pdf -h $Hclust_bar_h -lcex 1.3
cd ..


# Wilcox
mkdir Wilcox
cd Wilcox/
ln -s ../Community/phylum.xls ./
ln -s ../Community/class.xls ./
ln -s ../Community/order.xls ./
ln -s ../Community/family.xls ./
ln -s ../Community/genus.xls ./
ln -s ../Random_Forest/map.* ./

sub_Sample4otu_table.pl ../../process/otu_0.97/rarefac.otu_genus.xls ../groups.txt > rarefac.otu_genus.xls
ls map.*|awk -F "." '{print "wilcox.py -i phylum.xls -c "$2" -g map."$2".txt -o wilcox.phylum."$2".xls"}' >wilcox.sh
ls map.*|awk -F "." '{print "wilcox.py -i class.xls -c "$2" -g map."$2".txt -o wilcox.class."$2".xls"}' >>wilcox.sh
ls map.*|awk -F "." '{print "wilcox.py -i order.xls -c "$2" -g map."$2".txt -o wilcox.order."$2".xls"}' >>wilcox.sh
ls map.*|awk -F "." '{print "wilcox.py -i family.xls -c "$2" -g map."$2".txt -o wilcox.family."$2".xls"}' >>wilcox.sh
ls map.*|awk -F "." '{print "wilcox.py -i genus.xls -c "$2" -g map."$2".txt -o wilcox.genus."$2".xls"}' >>wilcox.sh
ls map.*|awk -F "." '{print "wilcox.py -i rarefac.otu_genus.xls -c "$2" -g map."$2".txt -o wilcox.otu."$2".xls"}' >>wilcox.sh
sh wilcox.sh
sed -i 's/"//g' *.xls
rm *.filtered.xls map.*.txt genus.xls phylum.xls class.xls order.xls family.xls wilcox.sh rarefac.otu_genus.xls *.r
cd ../


# Adonis
mkdir Adonis
cd Adonis
ln -s ../Random_Forest/*.samples.list ../Random_Forest/map.*.txt ./
ln -s ../Beta_diversity/bray_curtis_dm.txt ../Beta_diversity/unweighted_unifrac_dm.txt ../Beta_diversity/weighted_unifrac_dm.txt ./
for i in *.list; do filter_distance_matrix.py -i bray_curtis_dm.txt -o $i.bray_curtis_dm.txt --sample_id_fp $i; done
for i in *.list; do filter_distance_matrix.py -i unweighted_unifrac_dm.txt -o $i.unweighted_unifrac_dm.txt --sample_id_fp $i; done
for i in *.list; do filter_distance_matrix.py -i weighted_unifrac_dm.txt -o $i.weighted_unifrac_dm.txt --sample_id_fp $i; done

ls map.*.txt|awk -F "." '{print "compare_categories.py --method adonis -i "$2".samples.list.unweighted_unifrac_dm.txt -m map."$2".txt -c group -o adonis_"$2"_unweighted_unifrac"}' > adonis.sh
ls map.*.txt|awk -F "." '{print "compare_categories.py --method adonis -i "$2".samples.list.weighted_unifrac_dm.txt -m map."$2".txt -c group -o adonis_"$2"_weighted_unifrac"}' >> adonis.sh
ls map.*.txt|awk -F "." '{print "compare_categories.py --method adonis -i "$2".samples.list.bray_curtis_dm.txt -m map."$2".txt -c group -o adonis_"$2"_bray_curtis"}' >> adonis.sh

ls map.*.txt |sed 's/map.//'|sed 's/.txt//'|awk '{print "sed -i \"s/qiime\\.data\\$map\\[\\[opts\\$category\\]\\]/"$1"/\" adonis_"$1"_bray_curtis/adonis_results.txt"}'>>adonis.sh
ls map.*.txt |sed 's/map.//'|sed 's/.txt//'|awk '{print "sed -i \"s/qiime\\.data\\$map\\[\\[opts\\$category\\]\\]/"$1"/\" adonis_"$1"_unweighted_unifrac/adonis_results.txt"}'>>adonis.sh
ls map.*.txt |sed 's/map.//'|sed 's/.txt//'|awk '{print "sed -i \"s/qiime\\.data\\$map\\[\\[opts\\$category\\]\\]/"$1"/\" adonis_"$1"_weighted_unifrac/adonis_results.txt"}'>>adonis.sh
sh adonis.sh

cp /work/scripts/16s/adonis.head adonis.weighted_unifrac.xls
cp /work/scripts/16s/adonis.head adonis.bray_curtis.xls
cp /work/scripts/16s/adonis.head adonis.unweighted_unifrac.xls
cat adonis_*_weighted_unifrac/*|sed -n '/Df/,+1p' |grep -v Df|sed 's/\s\+/\t/g' >>adonis.weighted_unifrac.xls
#cat /work/scripts/16s/adonis.tail >>adonis.weighted_unifrac.xls
cat adonis_*_bray_curtis/*|sed -n '/Df/,+1p' |grep -v Df|sed 's/\s\+/\t/g' >>adonis.bray_curtis.xls
#cat /work/scripts/16s/adonis.tail >>adonis.bray_curtis.xls
cat adonis_*_unweighted_unifrac/*|sed -n '/Df/,+1p' |grep -v Df|sed 's/\s\+/\t/g' >>adonis.unweighted_unifrac.xls
#cat /work/scripts/16s/adonis.tail >>adonis.unweighted_unifrac.xls

rm -r adonis_*_weighted_unifrac adonis_*_bray_curtis adonis_*_unweighted_unifrac *.list  adonis.sh
cd ../

# Distance_boxplot
mkdir Distance_boxplot
cd Distance_boxplot
cat ../groups.txt|cut -f 1 > sample.list
filter_distance_matrix.py -i ../Beta_diversity/bray_curtis_dm.txt -o bray_curtis_dm.txt --sample_id_fp sample.list
filter_distance_matrix.py -i ../Beta_diversity/unweighted_unifrac_dm.txt -o unweighted_unifrac_dm.txt --sample_id_fp sample.list
filter_distance_matrix.py -i ../Beta_diversity/weighted_unifrac_dm.txt -o weighted_unifrac_dm.txt --sample_id_fp sample.list
make_distance_boxplots.py -d bray_curtis_dm.txt -m ../map-group.txt -f group -o bray_curtis -n 999 --suppress_all_between --suppress_all_within --save_raw_data --suppress_individual_between 
make_distance_boxplots.py -d unweighted_unifrac_dm.txt -m ../map-group.txt -f group -o unweighted_unifrac -n 999 --suppress_all_between --suppress_all_within --save_raw_data --suppress_individual_between
make_distance_boxplots.py -d weighted_unifrac_dm.txt -m ../map-group.txt -f group -o weighted_unifrac -n 999 --suppress_all_between --suppress_all_within --save_raw_data --suppress_individual_between
rm */*.pdf
cd bray_curtis
cp ../../map-group.txt ./
transposition.pl group_Distances.txt >group_Distances.tran.xls
distance_boxplot.pl -i group_Distances.tran.xls -w 4
cd ../unweighted_unifrac
cp ../../map-group.txt ./
transposition.pl group_Distances.txt >group_Distances.tran.xls
distance_boxplot.pl -i group_Distances.tran.xls -w 4
cd ../weighted_unifrac
cp ../../map-group.txt ./
transposition.pl group_Distances.txt >group_Distances.tran.xls
distance_boxplot.pl -i group_Distances.tran.xls -w 4
cd ../../

# Anosim
mkdir Anosim
cd Anosim
ln -s ../Random_Forest/*.samples.list ../Random_Forest/map.*.txt ./
ln -s ../Beta_diversity/bray_curtis_dm.txt ../Beta_diversity/unweighted_unifrac_dm.txt ../Beta_diversity/weighted_unifrac_dm.txt ./
for i in *.list; do filter_distance_matrix.py -i bray_curtis_dm.txt -o $i.bray_curtis_dm.txt --sample_id_fp $i; done
for i in *.list; do filter_distance_matrix.py -i unweighted_unifrac_dm.txt -o $i.unweighted_unifrac_dm.txt --sample_id_fp $i; done
for i in *.list; do filter_distance_matrix.py -i weighted_unifrac_dm.txt -o $i.weighted_unifrac_dm.txt --sample_id_fp $i; done
ls map.*.txt|awk -F "." '{print "compare_categories.py --method anosim -i "$2".samples.list.unweighted_unifrac_dm.txt -m map."$2".txt -c group -o anosim_"$2"_unweighted_unifrac"}' > anosim.sh
ls map.*.txt|awk -F "." '{print "compare_categories.py --method anosim -i "$2".samples.list.weighted_unifrac_dm.txt -m map."$2".txt -c group -o anosim_"$2"_weighted_unifrac"}' >> anosim.sh
ls map.*.txt|awk -F "." '{print "compare_categories.py --method anosim -i "$2".samples.list.bray_curtis_dm.txt -m map."$2".txt -c group -o anosim_"$2"_bray_curtis"}' >> anosim.sh

ls map.*.txt |sed 's/map.//'|sed 's/.txt//'|awk '{print "sed -i \"s/ANOSIM/"$1"/\" anosim_"$1"_bray_curtis/anosim_results.txt"}'>>anosim.sh
ls map.*.txt |sed 's/map.//'|sed 's/.txt//'|awk '{print "sed -i \"s/ANOSIM/"$1"/\" anosim_"$1"_unweighted_unifrac/anosim_results.txt"}'>>anosim.sh
ls map.*.txt |sed 's/map.//'|sed 's/.txt//'|awk '{print "sed -i \"s/ANOSIM/"$1"/\" anosim_"$1"_weighted_unifrac/anosim_results.txt"}'>>anosim.sh

ls map.*.txt |sed 's/map.//'|sed 's/.txt//'|awk '{print "transposition.pl anosim_"$1"_bray_curtis/anosim_results.txt > anosim_"$1"_bray_curtis/anosim_results.txt1"}' >>anosim.sh
ls map.*.txt |sed 's/map.//'|sed 's/.txt//'|awk '{print "transposition.pl anosim_"$1"_weighted_unifrac/anosim_results.txt > anosim_"$1"_weighted_unifrac/anosim_results.txt1"}'>>anosim.sh
ls map.*.txt |sed 's/map.//'|sed 's/.txt//'|awk '{print "transposition.pl anosim_"$1"_unweighted_unifrac/anosim_results.txt > anosim_"$1"_unweighted_unifrac/anosim_results.txt1"}'>>anosim.sh

#ls map.*.txt|awk -F "." '{print "Rscript /work/users/chaoliu/R_test/20180522_anosim/anosim.r -f "$2".samples.list.unweighted_unifrac_dm.txt -g map."$2".txt "}' >> anosim.sh
#ls map.*.txt|awk -F "." '{print "Rscript /work/users/chaoliu/R_test/20180522_anosim/anosim.r -f "$2".samples.list.weighted_unifrac_dm.txt -g map."$2".txt "}' >> anosim.sh
#ls map.*.txt|awk -F "." '{print "Rscript /work/users/chaoliu/R_test/20180522_anosim/anosim.r -f "$2".samples.list.bray_curtis_dm.txt -g map."$2".txt "}' >> anosim.sh

sh anosim.sh

cp /work/scripts/16s/anosim.head ./
cat anosim.head anosim_*_weighted_unifrac/anosim_results.txt1|grep -v 'method name' >anosim.weighted_unifrac.xls
cat anosim.head anosim_*_bray_curtis/anosim_results.txt1|grep -v 'method name' >anosim.bray_curtis.xls
cat anosim.head anosim_*_unweighted_unifrac/anosim_results.txt1|grep -v 'method name' >anosim.unweighted_unifrac.xls
rm -r anosim_*_weighted_unifrac anosim_*_bray_curtis anosim_*_unweighted_unifrac *.list *.txt anosim.sh anosim.head
cd ../

#Venn
mkdir Venn
cd Venn
sub_Sample4otu_table.pl ../../process/otu_0.97/rarefac.otu_genus.xls ../groups.txt > rarefac.otu_genus.xls
ln -s ../Random_Forest/map.* ./ 
ln -s ../map-group.txt
for i in *.txt; do venn_abundance.pl -t rarefac.otu_genus.xls -g $i ; done
rename sets.rarefac. "" *
rename venn.rarefac. venn. *
rename .txt. . *
cd ../

#Rank_abundance
cd Rank_abundance
ln -s ../OTU_Taxa/rarefac.otu_table.xls
rank_abundance.pl -i rarefac.otu_table.xls -gd ../map-group.txt  -o rankabundance.group.pdf -w 6 -h 5
rm rarefac.otu_table.xls
cd ../


#Hcluster_tree 
cd Hcluster_tree 
ln -s ../Beta_diversity/bray_curtis_dm.txt ../Beta_diversity/unweighted_unifrac_dm.txt ../Beta_diversity/weighted_unifrac_dm.txt ./
plot-hcluster_tree.pl -i bray_curtis_dm.txt -o .
plot-hcluster_tree.pl -i unweighted_unifrac_dm.txt -o .
plot-hcluster_tree.pl -i weighted_unifrac_dm.txt -o .
plot-tree.pl -i hcluster_tree_bray_curtis_dm.txt_average.tre -o hcluster_tree_bray_curtis_dm.txt_average.pdf -d ../map-group.txt -w 8 -h $Hcluster_tree_h
plot-tree.pl -i hcluster_tree_unweighted_unifrac_dm.txt_average.tre -o hcluster_tree_unweighted_unifrac_dm.txt_average.pdf -d ../map-group.txt -w 8 -h $Hcluster_tree_h
plot-tree.pl -i hcluster_tree_weighted_unifrac_dm.txt_average.tre -o hcluster_tree_weighted_unifrac_dm.txt_average.pdf -d ../map-group.txt -w 8 -h $Hcluster_tree_h
rm *.txt *.r
cd ../

#Specaccum
cd Specaccum

cd ../

#Cytoscape
#mkdir Cytoscape
#cd Cytoscape
#cat ../Heatmap/ALL.new.OTU-extract.all.sort.percent001.xls |cut -f 1|awk -F" " '{print $1}'|sed 1d >Key.otu.list
#sub_OTU4otu_table.pl ../OTU_Taxa/rarefac.otu_table.xls Key.otu.list > Key_otu_table.xls
#sub_Sample4otu_table.pl Key_otu_table.xls ../groups.txt > Key_otu_table.sub.xls
#transposition.pl Key_otu_table.sub.xls >Key_otu_table.trans.xls
#cor.spearman.pl -tax Key_otu_table.trans.xls -bc Key_otu_table.trans.xls > Key_otu.Spearman.cor.xls
#sed -i '2d' Key_otu.Spearman.cor.xls
#cat Key_otu.Spearman.cor.xls|awk '{$2="";print}'|sed 's/ /\t/g'|sed 's/\t\t/\t/' > Key_otu.Spearman.cor.new.xls
#mv Key_otu.Spearman.cor.new.xls Key_otu.Spearman.cor.xls
#matrix2Related_R_edge4CytoScape.pl Key_otu.Spearman.cor.xls >Key_otu.R.edge.xls
#ln -s ../map-group.txt ./
#richest4group.pl Key_otu_table.xls map-group.txt
#echo "" >Clinical.list
#taxa_table2Taxon_info.pl ../OTU_Taxa/rarefac.otu_taxa_table.xls >Taxon.info.xls
#toTaxon4CytoScape.pl Taxon.info.xls Key_otu.CAG.list Clinical.list >Key_otu.abundance.taxon.xls
#rm Taxon.info.xls Key.otu.list Key_otu_table.sub.xls Key_otu_table.trans.xls tax_bc_cor.txt  tax_bc_qva.txt  tmp.relation.xls Taxon.info.xls Clinical.list cmd.r map-group.txt
#cd ../


#Picrust
#unset PYTHONPATH
#source /work/users/chaoliu/miniconda3/bin/activate
#source activate python2
#
#mkdir  Picrust
#cd Picrust
#cp ../OTU_Taxa/rarefac.otu_taxa_table.xls ./
#sed -i 's/OTU ID/\#OTU ID/g' rarefac.otu_taxa_table.xls
#ln -s ../../process/otu_0.97/otu_reps.raw.fasta 
#picrust.v1.pl -otu rarefac.otu_taxa_table.xls -rep otu_reps.raw.fasta -o .
#sh predict.sh
#
#cd ko
#cp ../../map-group.txt ./
#sed -i "s/ /_/g" `grep " " -rl ./*.xls`
#run_lefse4fungene.v1.pl predictions_ko.xls map-group.txt /work/users/qiangli/Meta-Pipeline/mgpipev1/KEGG/ko00001-anno.xls Lefse_ko00001-anno
#plot-lda.v1.pl -i predictions_ko.L1.xls -o L1 -m map-group.txt -g group
#plot-lda.v1.pl -i predictions_ko.L2.xls -o L2 -m map-group.txt -g group
#plot-lda.v1.pl -i predictions_ko.L3.xls -o L3 -m map-group.txt -g group
#cd ../
#rm OTU2greengee.txt otu_reps.raw.fasta new_otu_table.biom  predict.sh rarefac.otu_taxa_table.xls
#cd ../

#cd cog 
#mkdir bar
#cog_profile2annotation.pl predictions_cog.xls > cog.category.function.xls
#cog_bar_eachSam.pl catalog_profile.xls bar/bar
#rm bar/bar.Catalog_Function.cog.pdf Rplots.pdf

#source deactivate
#source deactivate
cd ../
rm Nmds/* Pcoa/* Pca/*  Community/*ALL* Community/*sh Community/sample* Adonis/*.txt
rm -r Estimators Rarefactions/*xls Shannon_rarefac/*xls Adonis/*.txt
rm Heatmap/ALL.new.* Heatmap/map-group.txt Heatmap/OTU-extract.all.sort.xls Heatmap/dat.cor.* Heatmap/cmd.r Heatmap/bar.ALL.OTU-extract.all.sort.xls.pdf Heatmap/map-group.txt Heatmap/OTU-extract.all.sort.xls Heatmap/rarefac.otu_genus.xls

rm Heatmap/-ct Heatmap/-ct.xls Hclust_bar/cmd.r OTU_Taxa/otu_seqids.txt Venn/*.txt Venn/rarefac.otu_genus.xls Community/percent.*.xls
find ./ -name "*.biom" |xargs rm -rf

find ./ -type d -name "tax_summary_a" -exec rm -rf {} \;
find ./ -type d -name "tax_summary_r" -exec rm -rf {} \;

