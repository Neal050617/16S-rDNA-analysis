Color=${1:-none}

if [[ $Color != none ]]
then 
    Colour=`pwd`"/"$Color
else
    Colour=$Color
fi

echo $Colour

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
cat map-group.txt|awk '{print $1"\t"$1}' >order
cp order groups.txt
cp /work/scripts/16s/less_than_100/head.txt ./
cat head.txt map-group.txt |awk -F" " '{print $1"\t"$2}' >map-group.txt2
mv map-group.txt2 map-group.txt
cat map-group.txt |sed '1d'|awk -F"\t" '{if($2==last){group=group","$1}else{print group;group=$2":"$1};last=$2}END{print group}'|sed '1d' >map.txt

mkdir Split_groups
cd Split_groups
ln -s ../map-group.txt .
Rscript /work/users/chaoliu/scripts/split_group.R -g map-group.txt -m none -n 0
cd ..

mkdir OTU_analysis
cd OTU_analysis
cp ../../process/otu_0.97/rarefac.otu_taxa_table.xls .
sed -i 's/OTU ID/OTU_ID/g' rarefac.otu_taxa_table.xls
cp ../../process/otu_0.97/otu_reps.fasta .
cp ../map-group.txt .
less otu_reps.fasta |awk -F '_' '{print $1}' > otu_reps.raw.fasta
python /work/users/chaoliu/scripts/tax_split.v2.0.2.py -i rarefac.otu_taxa_table.xls -m map-group.txt -r otu_reps.raw.fasta -t 0.001

mkdir OTU_tree
cd OTU_tree
ln -s ../map.otu_reps.raw.fasta ./otu_reps.raw.fasta
ln -s ../unrooted.otu_tree_anno.0.001.xls .
muscle -in otu_reps.raw.fasta -out otu_reps.raw_aligned.fasta
FastTree -nt otu_reps.raw_aligned.fasta > unroot.otu_reps.raw_aligned.fasta.tre
Rscript /work/users/chaoliu/scripts/unroot_tree.R -t unroot.otu_reps.raw_aligned.fasta.tre -a unrooted.otu_tree_anno.0.001.xls -s phylum
rm Rplots.pdf
cd ..

cd Krona
for i in *.xls; do ktImportText $i -o $i.html;done
cd ..

mkdir Core_Microbiome
cd Core_Microbiome
ln -s ../../OTU_Taxa/rarefac.otu_genus.xls
ln -s ../map-group.txt .
sed 's/\t$//' rarefac.otu_genus.xls > fix.rarefac.otu_genus.xls
sed -i 's/ //g' fix.rarefac.otu_genus.xls
Rscript /work/users/chaoliu/scripts/core_microbiome.R -i fix.rarefac.otu_genus.xls -m map-group.txt

cd ../../
mv OTU_analysis/Core_Microbiome  OTU_analysis/Krona OTU_analysis/OTU_tree .
#rm -rf OTU_analysis

rm -rf OTU_Taxa
mv OTU_analysis OTU_Taxa
cd OTU_Taxa
rm rarefac.otu_taxa_table.xls otu_reps.raw.fasta otu_reps.fasta
mv map.otu_reps.raw.fasta otu_reps.raw.fasta
mv map.otu_table.xls rarefac.otu_table.xls
mv map.otu_table.percent.xls rarefac.otu_table.percent.xls
mv map.rarefac.otu_taxa_table.xls rarefac.otu_taxa_table.xls
mv map.otu_reps.raw.fasta otu_reps.raw.fasta
mv otu.genus.xls rarefac.otu_genus.xls

Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -o F -i rarefac.otu_genus.xls  -p 0.01  -m map-group.txt

cd ..

#Specaccum
cd Specaccum
ln -s ../OTU_Taxa/rarefac.otu_table.xls .
ln -s ../map-group.txt .
Rscript /work/users/chaoliu/scripts/Specaccum.R -i rarefac.otu_table.xls -m map-group.txt -c $Colour
rm *.xls
cd ../

#Venn
mkdir Venn
cd Venn
ln -s ../OTU_Taxa/rarefac.otu_genus.xls .
ln -s ../map-group.txt
Rscript /work/users/chaoliu/scripts/venn_upset.R -i rarefac.otu_genus.xls -g map-group.txt -m all -c $Colour
cd ../

#Rank_abundance
cd Rank_abundance
ln -s ../OTU_Taxa/rarefac.otu_table.xls
rank_abundance.pl -i rarefac.otu_table.xls -gd ../map-group.txt  -o rankabundance.group.pdf -w 6 -h 5
rm rarefac.otu_table.xls
cd ../

#Alpha_rarefac
mkdir Alpha_rarefac
cp ../process/otu_0.97/alpha_rarefac/* Alpha_rarefac -r
cd Alpha_rarefac/
ln -s ../map-group.txt ./
python /work/users/chaoliu/scripts/plot_alpha_diversity.py -i rarefaction
python /work/users/chaoliu/scripts/plot_alpha_diversity.py -i r_shannon
python /work/users/chaoliu/scripts/plot_alpha_diversity.py -i rarefaction -m map-group.txt -c $Colour
python /work/users/chaoliu/scripts/plot_alpha_diversity.py -i r_shannon -m map-group.txt -c $Colour
rm ../Rarefactions/*.pdf
rm ../Shannon_rarefac/*.pdf
mv rarefaction_* ../Rarefactions
mv r_shannon_* ../Shannon_rarefac

ln -s ../map-group.txt .
alpha_rarefac.pl map-group.txt > alpha_rarefac.summary.xls
Rscript /work/users/chaoliu/scripts/Alpha_diversity_box_plot.R -a alpha_rarefac.summary.xls -g map-group.txt -c $Colour

#Wilcox-alpha
Rscript /work/users/chaoliu/scripts/alpha_diversity_test.R -a alpha_rarefac.summary.xls -g map-group.txt -t FALSE
rm *.r mothur* otus.* *.txt
cd ../

#Hcluster_tree 
cd Hcluster_tree 
ln -s ../Beta_diversity/bray_curtis_dm.txt ../Beta_diversity/unweighted_unifrac_dm.txt ../Beta_diversity/weighted_unifrac_dm.txt ../map-group.txt ./ 
Rscript /work/users/chaoliu/scripts/plot_tree.R -i bray_curtis_dm.txt -m map-group.txt -t average -c $Colour
Rscript /work/users/chaoliu/scripts/plot_tree.R -i unweighted_unifrac_dm.txt -m map-group.txt -t average -c $Colour
Rscript /work/users/chaoliu/scripts/plot_tree.R -i weighted_unifrac_dm.txt -m map-group.txt -t average -c $Colour
rm *.txt *.r Rplots.pdf
cd ../

#Beta_diversity
cd Beta_diversity
ln -s ../map-group.txt
Rscript /work/users/chaoliu/scripts/Matrix_analysis.R -i bray_curtis_dm.txt -m T -g map-group.txt -p T
Rscript /work/users/chaoliu/scripts/Matrix_analysis.R -i unweighted_unifrac_dm.txt -m T -g map-group.txt -p T
Rscript /work/users/chaoliu/scripts/Matrix_analysis.R -i weighted_unifrac_dm.txt -m T -g map-group.txt -p T

Rscript /work/users/chaoliu/scripts/plot_matrix.R -j FALSE -i bray_curtis_dm.txt -g map-group.txt -c $Colour
Rscript /work/users/chaoliu/scripts/plot_matrix.R -j FALSE -i unweighted_unifrac_dm.txt -g map-group.txt -c $Colour
Rscript /work/users/chaoliu/scripts/plot_matrix.R -j FALSE -i weighted_unifrac_dm.txt -g map-group.txt -c $Colour

mkdir Anosim Adonis MRPP
mv *anosim* Anosim
mv *adonis* Adonis
mv *mrpp* MRPP

cd ..

#Community
rm -rf Community/
mv OTU_Taxa/Community .
cd Community
ln -s ../order ./sample
sed -i 's/ //g' *.xls

mkdir Community_barplot
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i phylum.xls -p 0
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i class.xls  -p 0.01
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i order.xls  -p 0.01
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i family.xls -p 0.01
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i genus.xls  -p 0.01
mv *.pdf ./Community_barplot/

ln -s ../order ./sample
cp ../map-group.txt ./
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -a 2 -i phylum.xls -p 0    -m map-group.txt -c $Colour
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -a 2 -i class.xls  -p 0.01 -m map-group.txt -c $Colour
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -a 2 -i order.xls  -p 0.01 -m map-group.txt -c $Colour
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -a 2 -i family.xls -p 0.01 -m map-group.txt -c $Colour
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -a 2 -i genus.xls  -p 0.01 -m map-group.txt -c $Colour
mkdir Community_barplot_groups/
mv *.pdf ./Community_barplot_groups/

mkdir Community_average
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -a 1 -i phylum.xls -p 0     -m map-group.txt -c $Colour
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -a 1 -i class.xls  -p 0.01  -m map-group.txt -c $Colour
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -a 1 -i order.xls  -p 0.01  -m map-group.txt -c $Colour
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -a 1 -i family.xls -p 0.01  -m map-group.txt -c $Colour
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -a 1 -i genus.xls  -p 0.01  -m map-group.txt -c $Colour -z 10
mv *.pdf ./Community_average

mkdir Community_bubble
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -b 2 -i phylum.xls -p 0    
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -b 2 -i class.xls  -p 0.01 
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -b 2 -i order.xls  -p 0.01 
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -b 2 -i family.xls -p 0.01 
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -b 2 -i genus.xls  -p 0.01 
mv *.pdf ./Community_bubble

mkdir Community_test_barplot Community_test_boxplot Community_test_stamp Community_test_complex
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -o F -i phylum.xls -p 0     -m map-group.txt -c $Colour
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -o F -i class.xls  -p 0.01  -m map-group.txt -c $Colour
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -o F -i order.xls  -p 0.01  -m map-group.txt -c $Colour
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -o F -i family.xls -p 0.01  -m map-group.txt -c $Colour
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -o F -i genus.xls  -p 0.01  -m map-group.txt -c $Colour
mv *barplot.pdf Community_test_barplot
mv *boxplot.pdf Community_test_boxplot
mv *stamp.pdf Community_test_stamp
mv *complex.pdf Community_test_complex
rm *.r percent*.xls

mkdir Heatmap_tax
plot-heatmap.pl -i sort.percent.genus.xls  -o heatmap.genus.pdf  -rt 0 -ct 0 -slas 2 -rlc 0.7 -clc 0.7 -w $Heatmap_w -h $Heatmap_keyh -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt -c $Colour
plot-heatmap.pl -i sort.percent.family.xls -o heatmap.family.pdf -rt 0 -ct 0 -slas 2 -rlc 0.7 -clc 0.7 -w $Heatmap_w -h $Heatmap_keyh -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt -c $Colour
plot-heatmap.pl -i sort.percent.order.xls  -o heatmap.order.pdf  -rt 0 -ct 0 -slas 2 -rlc 0.7 -clc 0.7 -w $Heatmap_w -h $Heatmap_keyh -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt -c $Colour
plot-heatmap.pl -i sort.percent.class.xls  -o heatmap.class.pdf  -rt 0 -ct 0 -slas 2 -rlc 0.7 -clc 0.7 -w $Heatmap_w -h $Heatmap_keyh -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt -c $Colour
plot-heatmap.pl -i sort.percent.phylum.xls -o heatmap.phylum.pdf -rt 0 -ct 0 -slas 2 -rlc 0.7 -clc 0.7 -w $Heatmap_w -h $Heatmap_keyh -lh 1:0.2:7:1 -marble $Heatmap_marble -cd map-group.txt -c $Colour

#python /work/users/chaoliu/scripts/pheatmap.py -i sort.percent.genus.xls -o pheatmap.genus.pdf -cd map-group.txt -clust_c F -cs 0 -rs 1 -scale none
rm *.pdf.xls
mv heatmap*.pdf Heatmap_tax

#Community_boxplot
mkdir Community_boxplot
cd Community_boxplot
ln -s ../../map-group.txt ./
sed -i 's/\#//g' map-group.txt
transposition.pl ../phylum.percent.xls >phylum.percents.xls
transposition.pl ../genus.percent.xls >genus.percents.xls
transposition.pl ../family.percent.xls >family.percents.xls
transposition.pl ../order.percent.xls >order.percents.xls
transposition.pl ../class.percent.xls >class.percents.xls
mkdir phylum genus family class order
Rscript /work/users/chaoliu/scripts/make_community_boxplot.r -i phylum.percents.xls -m map-group.txt -l T -c $Colour
mv *.pdf phylum
Rscript /work/users/chaoliu/scripts/make_community_boxplot.r -i genus.percents.xls  -m map-group.txt -l T -c $Colour
mv *.pdf genus
Rscript /work/users/chaoliu/scripts/make_community_boxplot.r -i order.percents.xls  -m map-group.txt -l T -c $Colour
mv *.pdf order
Rscript /work/users/chaoliu/scripts/make_community_boxplot.r -i family.percents.xls -m map-group.txt -l T -c $Colour
mv *.pdf family
Rscript /work/users/chaoliu/scripts/make_community_boxplot.r -i class.percents.xls  -m map-group.txt -l T -c $Colour
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
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt        -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt        -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt        -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt   -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt   -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt   -md PCoA -pc 2-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt -md PCoA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt -md PCoA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt -md PCoA -pc 2-3 -map map-group.txt -col $Colour
mv *.pdf noname/

Rscript /work/users/chaoliu/scripts/plotly_PCoA-3D.R -i bray_curtis_dm.txt        -c $Colour
Rscript /work/users/chaoliu/scripts/plotly_PCoA-3D.R -i weighted_unifrac_dm.txt   -c $Colour
Rscript /work/users/chaoliu/scripts/plotly_PCoA-3D.R -i unweighted_unifrac_dm.txt -c $Colour

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
mv *.pdf box/

rm cmd.r 
cd ../

############################################################ PCA ##########################################################
mkdir PCA
rm Pca -rf
cd PCA
ln -s ../map-group.txt ./
ln -s ../../process/otu_0.97/rarefac.otu_table.xls ./

#noname
mkdir noname
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i rarefac.otu_table.xls -md PCA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i rarefac.otu_table.xls -md PCA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i rarefac.otu_table.xls -md PCA -pc 2-3 -map map-group.txt -col $Colour
mv PCA*.pdf noname/

Rscript /work/users/chaoliu/scripts/plotly_PCoA-3D.R -i rarefac.otu_table.xls -d PCA

#name
mkdir name
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i rarefac.otu_table.xls -md PCA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i rarefac.otu_table.xls -md PCA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i rarefac.otu_table.xls -md PCA -pc 2-3 -map map-group.txt -col $Colour
mv PCA*.pdf name/

#ellipse
mkdir ellipse
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i rarefac.otu_table.xls -md PCA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i rarefac.otu_table.xls -md PCA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i rarefac.otu_table.xls -md PCA -pc 2-3 -map map-group.txt -col $Colour
mv PCA*.pdf ellipse/

#box
mkdir box
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i rarefac.otu_table.xls -md PCA -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i rarefac.otu_table.xls -md PCA -pc 1-3 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i rarefac.otu_table.xls -md PCA -pc 2-3 -map map-group.txt -col $Colour
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
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i bray_curtis_dm.txt        -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i weighted_unifrac_dm.txt   -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -i unweighted_unifrac_dm.txt -md NMDS -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf noname/

#name
mkdir name
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i bray_curtis_dm.txt        -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i weighted_unifrac_dm.txt   -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -lab T -i unweighted_unifrac_dm.txt -md NMDS -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf name/

#ellipse
mkdir ellipse
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i bray_curtis_dm.txt        -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i weighted_unifrac_dm.txt   -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -e T -i unweighted_unifrac_dm.txt -md NMDS -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf ellipse/

#box
mkdir box
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i bray_curtis_dm.txt        -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i weighted_unifrac_dm.txt   -md NMDS -pc 1-2 -map map-group.txt -col $Colour
python /work/users/chaoliu/scripts/gg_PCoA_NMDS_PCA.py -bx T -e T -i unweighted_unifrac_dm.txt -md NMDS -pc 1-2 -map map-group.txt -col $Colour
mv *.pdf box/

rm cmd.r 
cd ../

#### treebar ###################################################
cd Hclust_bar/
plot-treebar.pl -otu ../OTU_Taxa/rarefac.otu_table.xls -tax ../Community/phylum.xls -o treebar_otu_table_phylum.pdf -h $Hclust_bar_h -lcex 1.3
cd ..

# Picrust
cd Picrust
cd ko
ln  -s /work/users/chaoliu/scripts/ko00001-anno.xls .
ln  -s ../../map-group.txt .
Rscript /work/users/chaoliu/scripts/KEGG_tidyverse_table.R
mkdir KEGG_table KEGG_lefse ko_lefse
mv KEGG_Pathways_Level_* KEGG_table
cd ../../

cd Tax4Fun/
ln -s ../map-group.txt
Rscript  /work/users/chaoliu/scripts/Tax4Fun2lefse.R -i Tax4Fun_Pathways_ko.tsv
Rscript  /work/users/chaoliu/scripts/Tax4Fun2lefse.R -i Tax4Fun_Pathways_Metabolites.tsv
cd ../
######################################### Lefse #################################
#source /root/anaconda3/bin/activate
# 激活环境
source activate
# 退出环境
conda deactivate

conda activate qiime1

mkdir Lefse
cd Lefse/

cp ../OTU_Taxa/rarefac.otu_taxa_table.xls ./
sed -i 's/;s__.*//g' rarefac.otu_taxa_table.xls
#python /work/users/chaoliu/scripts/fix_otutable_tax.py -i rarefac.otu_taxa_table.xls
python2 /work/scripts/16s/tax_split.py -i rarefac.otu_taxa_table.xls
#准备数据
python2 /work/scripts/16s/make_lefse_data.py -i rarefac.otu_taxa_split_index_table.xls -g ../map-group.txt
sed -i '2d' lefse.data.xls

#处理数据分析
python2 /work/scripts/16s/format_input.py lefse.data.xls data.in -c 1 -o 1000000
python2 /work/scripts/16s/run_lefse.py data.in data.res -l 2
#画图
cat data.res|awk 'BEGIN{printf "Biomaker_names\tLogarithm value\tGroups\tLDA_value\tP_value\n"}{print $0}'>LDA.xls
python2 /work/scripts/16s/plot_res.py data.res lefse.pdf  --format pdf
python2 /work/scripts/16s/plot_cladogram.py data.res lefse.cladogram.pdf --format pdf
####################################################### lefes_genus
mkdir Lefse_all
mv ./* Lefse_all
ln -s Lefse_all/data.res .
Rscript /work/users/chaoliu/scripts/select_lefse.R -s g__
python2 /work/scripts/16s/plot_res.py g__data.res lefse.pdf  --format pdf
ln -s Lefse_all/lefse.data.xls .
cat g__data.res|awk 'BEGIN{printf "Biomaker_names\tLogarithm value\tGroups\tLDA_value\tP_value\n"}{print $0}'>LDA.xls
#ln  -s ../color.txt .
ln -s ../map-group.txt .

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

# Wilcox
mkdir Wilcox
cd Wilcox/
ln -s ../Community/phylum.xls ./
ln -s ../Community/class.xls ./
ln -s ../Community/order.xls ./
ln -s ../Community/family.xls ./
ln -s ../Community/genus.xls ./
ln -s ../Random_Forest/map.* ./
ln -s ../OTU_Taxa/rarefac.otu_genus.xls ./

ls map.*|awk -F "." '{print "python /work/users/chaoliu/scripts/wilcox.py -i phylum.xls -c "$2" -g map."$2".txt -o wilcox.phylum."$2".xls"}' >wilcox.sh
ls map.*|awk -F "." '{print "python /work/users/chaoliu/scripts/wilcox.py -i class.xls -c "$2" -g map."$2".txt -o wilcox.class."$2".xls"}' >>wilcox.sh
ls map.*|awk -F "." '{print "python /work/users/chaoliu/scripts/wilcox.py -i order.xls -c "$2" -g map."$2".txt -o wilcox.order."$2".xls"}' >>wilcox.sh
ls map.*|awk -F "." '{print "python /work/users/chaoliu/scripts/wilcox.py -i family.xls -c "$2" -g map."$2".txt -o wilcox.family."$2".xls"}' >>wilcox.sh
ls map.*|awk -F "." '{print "python /work/users/chaoliu/scripts/wilcox.py -i genus.xls -c "$2" -g map."$2".txt -o wilcox.genus."$2".xls"}' >>wilcox.sh
ls map.*|awk -F "." '{print "python /work/users/chaoliu/scripts/wilcox.py -i rarefac.otu_genus.xls -c "$2" -g map."$2".txt -o wilcox.otu."$2".xls"}' >>wilcox.sh
sh wilcox.sh
sed -i 's/"//g' *.xls
rm *.filtered.xls map.*.txt genus.xls phylum.xls class.xls order.xls family.xls wilcox.sh rarefac.otu_genus.xls *.r
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
python /work/users/chaoliu/scripts/pheatmap.py -i rarefac.otu_genus.xls -cd map-group.txt -o pheatmap.otu.top50.pdf -rtop 50 -clust_c F -cs 0 -rs 1 -scale none -cc $Colour
ln -s ../Random_Forest/OTU-extract.all.sort.xls ./
bar_pie.pl -i OTU-extract.all.sort.xls -pie F -p 0.01
sed -i '/Others/d' ALL.new.OTU-extract.all.sort.xls
mv ALL.new.OTU-extract.all.sort.xls ALL.new.OTU-extract.all.sort.percent001.xls
bar_pie.pl -i OTU-extract.all.sort.xls -pie F -p 0.03
sed -i '/Others/d' ALL.new.OTU-extract.all.sort.xls
mv ALL.new.OTU-extract.all.sort.xls ALL.new.OTU-extract.all.sort.percent003.xls
python /work/users/chaoliu/scripts/pheatmap.py -i OTU-extract.all.sort.xls                    -o pheatmap.keyotu.pdf            -cd map-group.txt -clust_c F -cs 0 -rs 1 -scale none -cc $Colour
python /work/users/chaoliu/scripts/pheatmap.py -i ALL.new.OTU-extract.all.sort.percent001.xls -o pheatmap.keyotu.percent001.pdf -cd map-group.txt -clust_c F -cs 0 -rs 1 -scale none -cc $Colour
python /work/users/chaoliu/scripts/pheatmap.py -i ALL.new.OTU-extract.all.sort.percent003.xls -o pheatmap.keyotu.percent003.pdf -cd map-group.txt -clust_c F -cs 0 -rs 1 -scale none -cc $Colour
rm ALL*.xls dat* bar* cmd.r ./*ct*
cd ../

###################################################### Tax4Fun ###############################################################################################
cd Tax4Fun/
mkdir Tax4Fun_table Metabolites_lefse ko_lefse

#处理数据分析
python2 /work/scripts/16s/format_input.py ko.lefse.data.tsv ko.lefse.data.in -c 1 -o 1000000
python2 /work/scripts/16s/run_lefse.py ko.lefse.data.in ko.lefse.data.res -l 2
#画图
python2 /work/scripts/16s/plot_res.py ko.lefse.data.res ko.lefse.data.pdf  --format pdf
cat ko.lefse.data.res | awk 'BEGIN{printf "Biomaker_names\tLogarithm value\tGroups\tLDA_value\tP_value\n"}{print $0}'>ko_LDA.xls

mv ko.lefse.data* ko_LDA.xls ko_lefse/

#处理数据分析
python2 /work/scripts/16s/format_input.py Metabolites.lefse.data.tsv Metabolites.lefse.data.in -c 1 -o 1000000
python2 /work/scripts/16s/run_lefse.py Metabolites.lefse.data.in Metabolites.lefse.data.res -l 2
#画图
python2 /work/scripts/16s/plot_res.py Metabolites.lefse.data.res Metabolites.lefse.data.pdf  --format pdf
cat Metabolites.lefse.data.res | awk 'BEGIN{printf "Biomaker_names\tLogarithm value\tGroups\tLDA_value\tP_value\n"}{print $0}'>Metabolites.lefse_LDA.xls

mv Metabolites.lefse.data* Metabolites.lefse_LDA.xls Metabolites_lefse/

cd ../
##################################################### Picrust ##############################################################
cd Picrust/ko/
#处理数据分析
python2 /work/scripts/16s/format_input.py KEGG.lefse.data.tsv KEGG.lefse.data.in -c 1 -o 1000000
python2 /work/scripts/16s/run_lefse.py KEGG.lefse.data.in KEGG.lefse.data.res -l 2
#画图
python2 /work/scripts/16s/plot_res.py KEGG.lefse.data.res KEGG.lefse.data.pdf  --format pdf
python2 /work/scripts/16s/plot_cladogram.py KEGG.lefse.data.res KEGG.lefse.data.cladogram.pdf --format pdf --labeled_start_lev 1 --labeled_stop_lev 3 --abrv_start_lev 2 --abrv_stop_lev 3
cat KEGG.lefse.data.res | awk 'BEGIN{printf "Biomaker_names\tLogarithm value\tGroups\tLDA_value\tP_value\n"}{print $0}'>LDA.xls
mv KEGG.lefse.data* *.pdf LDA.xls KEGG_lefse/

#处理数据分析
python2 /work/scripts/16s/format_input.py ko.lefse.data.tsv ko.lefse.data.in -c 1 -o 1000000
python2 /work/scripts/16s/run_lefse.py ko.lefse.data.in ko.lefse.data.res -l 2
#画图
python2 /work/scripts/16s/plot_res.py ko.lefse.data.res ko.lefse.data.pdf  --format pdf
cat ko.lefse.data.res | awk 'BEGIN{printf "Biomaker_names\tLogarithm value\tGroups\tLDA_value\tP_value\n"}{print $0}'>ko_LDA.xls

mv ko.lefse.data* *.pdf ko_LDA.xls ko_lefse/
################################# L3 
ln -s KEGG_lefse/KEGG.lefse.data.res .
Rscript /work/users/chaoliu/scripts/select_lefse.R -s L3_ -i KEGG.lefse.data.res
python2 /work/scripts/16s/plot_res.py L3_data.res lefse.pdf  --format pdf
ln -s KEGG_lefse/KEGG.lefse.data.tsv .
cat L3_data.res|awk 'BEGIN{printf "Biomaker_names\tLogarithm value\tGroups\tLDA_value\tP_value\n"}{print $0}'>LDA.xls

##################################################################################################################
conda deactivate

cd KEGG_table
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i KEGG_Pathways_Level_1.tsv   -p 0    -c $Colour
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i KEGG_Pathways_Level_2.tsv   -p 0.01 -c $Colour
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i KEGG_Pathways_Level_3.tsv   -p 0.01 -c $Colour
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i KEGG_Pathways_Level_ko.tsv  -p 0.01 -c $Colour

ln -s ../../../map-group.txt
mkdir Community_test_barplot Community_test_boxplot Community_test_stamp Community_test_complex
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -o F -i KEGG_Pathways_Level_1.tsv  -p 0     -m map-group.txt -c $Colour
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -o F -i KEGG_Pathways_Level_2.tsv  -p 0.01  -m map-group.txt -c $Colour
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -o F -i KEGG_Pathways_Level_3.tsv  -p 0.01  -m map-group.txt -c $Colour
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -o F -i KEGG_Pathways_Level_ko.tsv -p 0.01  -m map-group.txt -c $Colour
mv *barplot.pdf Community_test_barplot
mv *boxplot.pdf Community_test_boxplot
mv *stamp.pdf Community_test_stamp
mv *complex.pdf Community_test_complex

cd ../../../

cd  Tax4Fun/Tax4Fun_table/
mv ../Tax4Fun_Pathways*.tsv .
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i Tax4Fun_Pathways_Metabolites.tsv   -p 0.01 -c $Colour
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -i Tax4Fun_Pathways_ko.tsv            -p 0.01 -c $Colour

ln -s ../../map-group.txt
mkdir Community_test_barplot Community_test_boxplot Community_test_stamp Community_test_complex
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -o F -i Tax4Fun_Pathways_Metabolites.tsv -p 0.01 -m map-group.txt -c $Colour
Rscript /work/users/chaoliu/scripts/Community_plot_test.R -t rank -o F -i Tax4Fun_Pathways_ko.tsv          -p 0.01 -m map-group.txt -c $Colour
mv *barplot.pdf Community_test_barplot
mv *boxplot.pdf Community_test_boxplot
mv *stamp.pdf Community_test_stamp
mv *complex.pdf Community_test_complex

cd ../../

cd Lefse/Lefse_all/
Rscript /work/users/chaoliu/scripts/plot_lefse.v2.R -i lefse.data.xls -r LDA.xls -l 2 -g ../map-group.txt  -c $Colour
cd ../
Rscript /work/users/chaoliu/scripts/plot_lefse.v2.R -i lefse.data.xls -r LDA.xls -l 2 -g map-group.txt  -c $Colour
mkdir Lefse_genus
ls | grep -v Lefse_ | xargs -t -I '{}' mv {} Lefse_genus
#mv !(Lefse_*) Lefse_genus
cd ../Picrust/ko/KEGG_lefse
Rscript /work/users/chaoliu/scripts/plot_lefse.v2.R -i KEGG.lefse.data.tsv -r LDA.xls -l 2 -m kegg -g ../map-group.txt -c $Colour
cd ../
Rscript /work/users/chaoliu/scripts/plot_lefse.v2.R -i KEGG.lefse.data.tsv -r LDA.xls -l 2 -g map-group.txt -m kegg -c $Colour
cd ko_lefse/
Rscript /work/users/chaoliu/scripts/plot_lefse.v2.R -i ko.lefse.data.tsv -r ko_LDA.xls -l 2 -m ko -g ../map-group.txt -c $Colour
cd ../../../Tax4Fun/Metabolites_lefse/
Rscript /work/users/chaoliu/scripts/plot_lefse.v2.R -i Metabolites.lefse.data.tsv -r Metabolites.lefse_LDA.xls -l 2 -m ko  -g ../map-group.txt -c $Colour
cd ../ko_lefse
Rscript /work/users/chaoliu/scripts/plot_lefse.v2.R -i ko.lefse.data.tsv -r ko_LDA.xls -l 2 -m ko  -g ../map-group.txt -c $Colour
cd ../../

rm NMDS/*.txt PCoA/*.txt PCA/*.txt  Community/*ALL* Community/*sh Community/sample*
rm -r Estimators
rm Rarefactions/*.rarefaction Shannon_rarefac/*.r_shannon
rm Heatmap/ALL.new.* Heatmap/map-group.txt 
rm Heatmap/OTU-extract.all.sort.xls Heatmap/dat.cor.* Heatmap/cmd.r Heatmap/bar.ALL.OTU-extract.all.sort.xls.pdf Heatmap/map-group.txt Heatmap/OTU-extract.all.sort.xls Heatmap/rarefac.otu_genus.xls
rm Heatmap/-ct Heatmap/-ct.xls 
rm Hclust_bar/cmd.r OTU_Taxa/otu_seqids.txt Venn/*.txt Venn/rarefac.otu_genus.xls Community/percent.*.xls

#find ./ -name "*.biom" |xargs rm -rf
find ./ -type d -name "tax_summary_a"| xargs rm -rf
find ./ -type d -name "tax_summary_r"| xargs rm -rf
find . -name "Rplots.pdf" | xargs rm -rf
find . -name "*cmd.r"| xargs rm -rf
find . -name "otu_seqids.txt"| xargs rm -rf
# 删除失效软连接
for file in `find . -type l`
do
    if [ ! -e $file ]
    then
        echo "rm $file"
        rm -f $file
    fi
done
