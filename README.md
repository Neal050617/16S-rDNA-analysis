# 16S-rDNA-analysis
This repository was designed to save the pipeline and scripts of analysis of 16 SrRNA gene datasets.

# R 
R version 3.6.0 (2019-04-26)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so

locale:
 [1] LC_CTYPE=zh_CN.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=zh_CN.UTF-8        LC_COLLATE=zh_CN.UTF-8    
 [5] LC_MONETARY=zh_CN.UTF-8    LC_MESSAGES=zh_CN.UTF-8   
 [7] LC_PAPER=zh_CN.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=zh_CN.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] compiler_3.6.0

> .packages(all.available=T)
  [1] "dynamicTreeCut"       "fastcluster"          "foreach"             
  [4] "iterators"            "abind"                "acepack"             
  [7] "ade4"                 "animation"            "AnnotationDbi"       
 [10] "ape"                  "askpass"              "assertthat"          
 [13] "backports"            "base"                 "base64enc"           
 [16] "bayesm"               "BH"                   "Biobase"             
 [19] "BiocFileCache"        "BiocGenerics"         "BiocManager"         
 [22] "BiocVersion"          "biom"                 "biomaRt"             
 [25] "bit"                  "bit64"                "bitops"              
 [28] "blob"                 "boot"                 "brew"                
 [31] "broom"                "Cairo"                "callr"               
 [34] "car"                  "carData"              "caTools"             
 [37] "cellranger"           "checkmate"            "class"               
 [40] "cli"                  "clipr"                "clisymbols"          
 [43] "cluster"              "clusterGeneration"    "clusterSim"          
 [46] "coda"                 "codetools"            "colorspace"          
 [49] "combinat"             "commonmark"           "compiler"            
 [52] "compositions"         "corrplot"             "covr"                
 [55] "cowplot"              "crayon"               "crosstalk"           
 [58] "curl"                 "data.table"           "datasets"            
 [61] "DBI"                  "dbplyr"               "dendextend"          
 [64] "DEoptimR"             "desc"                 "devtools"            
 [67] "digest"               "dlstats"              "doParallel"          
 [70] "dplyr"                "DT"                   "e1071"               
 [73] "ellipse"              "ellipsis"             "evaluate"            
 [76] "exactRankTests"       "expm"                 "factoextra"          
 [79] "FactoMineR"           "fansi"                "farver"              
 [82] "fastmap"              "fastmatch"            "flashClust"          
 [85] "forcats"              "forecast"             "foreign"             
 [88] "formatR"              "formattable"          "Formula"             
 [91] "fracdiff"             "fs"                   "furrr"               
 [94] "futile.logger"        "futile.options"       "future"              
 [97] "gdata"                "generics"             "getopt"              
[100] "gganimate"            "ggforce"              "ggimage"             
[103] "gglayer"              "ggplot2"              "ggplotify"           
[106] "ggpubr"               "ggraph"               "ggrepel"             
[109] "ggsci"                "ggsignif"             "ggtern"              
[112] "ggtree"               "gh"                   "git2r"               
[115] "githubinstall"        "globals"              "glue"                
[118] "GO.db"                "gplots"               "graphics"            
[121] "graphlayouts"         "grDevices"            "greybox"             
[124] "grid"                 "gridExtra"            "gridGraphics"        
[127] "gtable"               "gtools"               "GUniFrac"            
[130] "haven"                "hexbin"               "highr"               
[133] "Hmisc"                "hms"                  "htmlTable"           
[136] "htmltools"            "htmlwidgets"          "httpuv"              
[139] "httr"                 "igraph"               "impute"              
[142] "ini"                  "IRanges"              "jpeg"                
[145] "jsonlite"             "KernSmooth"           "knitr"               
[148] "labeling"             "lambda.r"             "lamW"                
[151] "later"                "latex2exp"            "lattice"             
[154] "latticeExtra"         "lavaan"               "lazyeval"            
[157] "leaps"                "lifecycle"            "listenv"             
[160] "lme4"                 "lmtest"               "lubridate"           
[163] "magick"               "magrittr"             "manipulateWidget"    
[166] "maps"                 "maptools"             "markdown"            
[169] "MASS"                 "Matrix"               "MatrixModels"        
[172] "matrixStats"          "meme"                 "memoise"             
[175] "methods"              "mgcv"                 "microbiomeViz"       
[178] "mime"                 "miniUI"               "minqa"               
[181] "mnormt"               "mockery"              "modelr"              
[184] "munsell"              "nlme"                 "nloptr"              
[187] "nnet"                 "numDeriv"             "openssl"             
[190] "openxlsx"             "optparse"             "parallel"            
[193] "patchwork"            "pbivnorm"             "pbkrtest"            
[196] "PerformanceAnalytics" "permute"              "phangorn"            
[199] "pheatmap"             "phytools"             "picante"             
[202] "pillar"               "pkgbuild"             "pkgconfig"           
[205] "pkgload"              "plogr"                "plotly"              
[208] "plotrix"              "plyr"                 "png"                 
[211] "polyclip"             "polynom"              "pracma"              
[214] "praise"               "preprocessCore"       "prettyunits"         
[217] "pROC"                 "processx"             "progress"            
[220] "promises"             "proto"                "ps"                  
[223] "psych"                "purrr"                "qiimer"              
[226] "qrcode"               "quadprog"             "quantmod"            
[229] "quantreg"             "R.cache"              "R.methodsS3"         
[232] "R.oo"                 "R.utils"              "R2HTML"              
[235] "R6"                   "randomForest"         "rappdirs"            
[238] "rcmdcheck"            "RColorBrewer"         "Rcpp"                
[241] "RcppArmadillo"        "RcppEigen"            "RcppParallel"        
[244] "readr"                "readxl"               "rematch"             
[247] "remotes"              "reprex"               "reshape2"            
[250] "rex"                  "rgl"                  "rio"                 
[253] "RISmed"               "rJava"                "RJSONIO"             
[256] "rlang"                "rmarkdown"            "robustbase"          
[259] "roxygen2"             "rpart"                "rprojroot"           
[262] "RSQLite"              "rstudioapi"           "rsvg"                
[265] "rvcheck"              "rversions"            "rvest"               
[268] "S4Vectors"            "scales"               "scatterplot3d"       
[271] "scholar"              "selectr"              "sessioninfo"         
[274] "shiny"                "showtext"             "showtextdb"          
[277] "smooth"               "sourcetools"          "sp"                  
[280] "SparseM"              "spatial"              "splines"             
[283] "statmod"              "stats"                "stats4"              
[286] "stringi"              "stringr"              "survival"            
[289] "sys"                  "sysfonts"             "Tax4Fun"             
[292] "tcltk"                "tensorA"              "testthat"            
[295] "tibble"               "tidygraph"            "tidyr"               
[298] "tidyselect"           "tidytree"             "tidyverse"           
[301] "timeDate"             "tinytex"              "tools"               
[304] "treeio"               "tseries"              "TTR"                 
[307] "tweenr"               "UpSetR"               "urca"                
[310] "usethis"              "utf8"                 "utils"               
[313] "vctrs"                "vegan"                "VennDiagram"         
[316] "venneuler"            "viridis"              "viridisLite"         
[319] "webshot"              "WGCNA"                "whisker"             
[322] "withr"                "xfun"                 "xlsx"                
[325] "xlsxjars"             "XML"                  "xml2"                
[328] "xopen"                "xtable"               "xts"                 
[331] "yaml"                 "yyplot"               "zeallot"             
[334] "zip"                  "zoo" 



