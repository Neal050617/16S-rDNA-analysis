package_list <- c("vegan","reshape2","ggplot2","tidyverse","RColorBrewer","gridExtra",
                  "grid","magrittr","optparse","patchwork","plotly","scatterplot3d","rgl", "car")

for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos="http://cran.r-project.org")
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

if (TRUE){
  option_list <- list(
    make_option(c("-t", "--type"), type="logical", default=T,
                help="用现成的还是重新计算"),
    make_option(c("-i", "--input"), type="character", default="bray_curtis_dm.txt",
                help="丰度表格,或是矩阵"),
    make_option(c("-g", "--Gp"), type="character", default="map-group.txt",
                help="分组文件"),
    make_option(c("-d", "--mode"), type="character", default="PCoA",
                help="分析模式"),
    make_option(c("-c", "--color"), type="character", default="none",
                help="指定颜色：color.txt"),
    make_option(c("-s", "--shape"), type="character", default="none",
                help="指定形状设置"),
    make_option(c("-r", "--range"), type="logical", default=T,
                help="长宽高一样")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
}

options("endocing"="UTF-8")
options(scipen = 200)

mycol <- c("#CD3333","#483D8B","#458B00","#EEC900","#BFBC3B","#6E8B3D","#00688B","#C10077","#CAAA76","#EEEE00","#458B00","#8B4513","#008B8B","#6E8B3D","#8B7D6B","#7FFF00","#CDBA96","#ADFF2F")
mypch <- c(16,15,17,18,7,8,9,10,11,12,13,14,21,22,23,24,25,3,4)

#分组文件
if (opts$Gp != "none"){
  Gp <- read_tsv(opts$Gp) %>% 
    rename("SampleID" = colnames(.)[1],"group" = colnames(.)[2]) %>%
    mutate(group = fct_inorder(group),SampleID = fct_inorder(SampleID))
  # 分组设计文件
  gp <- as.vector(unique(unlist(Gp[,"group"])))
}else{
  Gp <- read_tsv(opts$input,n_max=1) %>% colnames(.) %>% .[-1] %>% 
    as_tibble %>% rename("SampleID" = colnames(.)[1]) %>%
    mutate(SampleID = fct_inorder(SampleID))
}

# 颜色和形状
if (opts$Gp != "none"){
  if (opts$color != "none"){# opts$color = c("color.txt")
    sc <- read_tsv(opts$color,col_names = FALSE) %>%
      rename_all(~c("group","Color")) %>%
      filter(group %in% levels(Gp$group)) %>%
      mutate(group = fct_inorder(group))
    
    Gp$group <- factor(Gp$group,levels = as.vector(sc$group))
  } else {
    sc <- cbind(levels(Gp$group),mycol[1:nlevels(Gp$group)]) %>% 
      as.data.frame(stringsAsFactors = FALSE) %>%
      rename_all(~c("group","Color")) %>%
      mutate(group = factor(group,levels=levels(Gp$group)))
  }
  Gp <- Gp %>% inner_join(sc)
  
  #设置形状
  if (opts$shape != "none"){
    sp <- read_tsv(opts$shape,col_names = FALSE) %>% 
      rename_all(~c("SampleID","group","shape")) %>% 
      mutate(SampleID = factor(SampleID,levels=levels(Gp$SampleID)),
             group = factor(group,levels=levels(Gp$group)))
  } else{
    sp <- Gp %>% inner_join(cbind(levels(Gp$group),mypch[1:nlevels(Gp$group)]) %>% 
                              as.data.frame(stringsAsFactors = FALSE) %>% rename_all(~c("group","shape"))) %>%
      mutate(SampleID = factor(SampleID,levels=levels(Gp$SampleID)),
             group = factor(group,levels=levels(Gp$group)),
             shape = as.numeric(shape))
  }
  Gp <- Gp %>% inner_join(sp)
}

if (!opts$type){
  #读取文件_PCoA矩阵
  da <- read_tsv(opts$input) %>% rename("ID" = colnames(.)[1]) %>%
    select(c("ID",Gp$SampleID)) %>% filter(ID %in% colnames(.)) %>%
    mutate(ID = factor(ID,levels = colnames(.)[-1])) %>% arrange(ID) %>%
    as.data.frame(.) %>% .[,-1]
  
  rownames(da) <- colnames(da)
  
  if(opts$mode=="PCA"){
    pca <- prcomp(da)
    pc12 <- as.data.frame(pca$rotation[,c(1,2,3)])
  }
  if(opts$mode=="PCoA"){
    pca <- prcomp(da)
    pc12 <- as.data.frame(pca$x[,c(1,2,3)])
  }
  pc <- summary(pca)$importance[2,]*100
  pca.sum <- summary(pca)
  xlab <- paste("PC1_",round(pc[1],2),"%",sep="")
  ylab <- paste("PC2_",round(pc[2],2),"%",sep="")
  zlab <- paste("PC3_",round(pc[3],2),"%",sep="")
  ##输出文件
  write.table(data.frame(sample_ID=rownames(pca$rotation),pca$rotation),quote=F,sep="	",row.names = F,
              file=paste0(opts$mode,"_",opts$input,"_PC","_","rotation.xls"))
  write.table(data.frame(sample_ID=rownames(predict(pca)),predict(pca)),quote=F,sep="	",row.names = F,
              file=paste0(opts$mode,"_",opts$input,"_PC","_","sites.xls"))
  write.table(data.frame(sample_ID=rownames(pca.sum$importance),pca.sum$importance),quote=F,sep="	",row.names = F,
              file=paste0(opts$mode,"_",opts$input,"_PC","_","importance.xls")) 
  
  # 聚类点
  if (opts$Gp != "none"){
    # 跟分组合并
    pc12 <- pc12 %>% rownames_to_column() %>% rename("SampleID" = colnames(.)[1]) %>% 
      inner_join(.,Gp) %>% as_tibble
    
    centroids <- aggregate(cbind(PC1,PC2,PC3)~group,pc12,mean)
    f         <- function(z)sd(z)/sqrt(length(z)) # function to calculate std.err
    se        <- aggregate(cbind(se.x=PC1,se.y=PC2,se.z=PC3)~group,pc12,f)
    ct1 <- merge(pc12,aggregate(cbind(mean.x=PC1,mean.y=PC2,mean.z=PC3)~group,pc12,mean),by="group")
    ct2 <- merge(centroids,se, by="group")# add std.err column to centroids
  }
  
}else{
  pc <- read_tsv(str_c(opts$mode,"_",opts$input,"_PC","_","importance.xls")) %>%
    .[2,2:4] %>% as.character %>% as.numeric
  xlab <- paste("PC1_",round(pc[1]*100,2),"%",sep="")
  ylab <- paste("PC2_",round(pc[2]*100,2),"%",sep="")
  zlab <- paste("PC3_",round(pc[3]*100,2),"%",sep="")
  if(opts$mode=="PCA"){
    pc12 <- read_tsv(str_c(opts$mode,"_",opts$input,"_PC","_","rotation.xls")) %>% 
      .[,1:4] %>% rename("SampleID" = colnames(.)[1])
  }
  if(opts$mode=="PCoA"){
    pc12 <- read_tsv(str_c(opts$mode,"_",opts$input,"_PC","_","sites.xls")) %>% 
      .[,1:4] %>% rename("SampleID" = colnames(.)[1])
  }
  # 聚类点
  if (opts$Gp != "none"){
    # 跟分组合并
    pc12 <- pc12 %>% inner_join(.,Gp) %>% as_tibble
    
    centroids <- aggregate(cbind(PC1,PC2,PC3)~group,pc12,mean)
    f         <- function(z)sd(z)/sqrt(length(z)) # function to calculate std.err
    se        <- aggregate(cbind(se.x=PC1,se.y=PC2,se.z=PC3)~group,pc12,f)
    ct1 <- merge(pc12,aggregate(cbind(mean.x=PC1,mean.y=PC2,mean.z=PC3)~group,pc12,mean),by="group")
    ct2 <- merge(centroids,se, by="group")# add std.err column to centroids
  }
}

############################################### 画图函数#####################################
# 3D
GGplot8 <- function(pc12,Range){
  # custom grid style
  axx <- list(
    #gridcolor='rgb(255, 255, 255)',
    #zerolinecolor='rgb(255, 255, 255)',
    #showbackground=TRUE,
    #backgroundcolor='rgb(230, 230,230)',
    autorange = F, 
    aspectmode = 'manual'
  )
  p <- plot_ly(pc12, x = ~PC1, y = ~PC2, z = ~PC3, color = ~group, 
               colors = unique(pc12$Color), size = 10) %>%
    add_markers() %>%
    layout(scene = list(xaxis = append(list(title = xlab,range = c(-Range[1],Range[1])),axx),
                        yaxis = append(list(title = ylab,range = c(-Range[2],Range[2])),axx),
                        zaxis = append(list(title = zlab,range = c(-Range[3],Range[3])),axx),
                        aspectratio = list(x = 1, y = 1, z = 1)))
}

addgrids3d <- function(x, y=NULL, z=NULL, grid = TRUE,
                       col.grid = "grey", lty.grid = par("lty"),
                       lab = par("lab"), lab.z = mean(lab[1:2]),
                       scale.y = 1, angle = 40,
                       xlim=NULL, ylim=NULL, zlim=NULL, Rangexyz=NULL){

  if(inherits(x, c("matrix", "data.frame"))){
    x <- as.data.frame(x)
    y <- unlist(x[,2])
    z <- unlist(x[,3])
    x <- unlist(x[,1])
  }
  
  p.lab <- par("lab")
  
  angle <- (angle%%360)/90
  yz.f <- scale.y * abs(if (angle < 1) angle else if (angle >3) angle - 4 else 2 - angle)
  yx.f <- scale.y * (if (angle < 2) 1 - angle else angle - 3)
  
  # x axis range
  if (is.null(Rangexyz)){
    x.range <- range(x[is.finite(x)], xlim)
  }else{
    x.range <- c(-Rangexyz[1],Rangexyz[1])
  }
  x.prty <- pretty(x.range, n = lab[1], min.n = max(1, min(0.5 *lab[1], p.lab[1])))
  x.scal <- round(diff(x.prty[1:2]), digits = 12)
  x <- x/x.scal
  x.range <- range(x.prty)/x.scal
  x.max <- ceiling(x.range[2])
  x.min <- floor(x.range[1])
  if (!is.null(xlim)) {
    x.max <- max(x.max, ceiling(xlim[2]/x.scal))
    x.min <- min(x.min, floor(xlim[1]/x.scal))
  }
  x.range <- range(x.min, x.max)
  
  # y axis range
  if (is.null(Rangexyz)){
    y.range <- range(y[is.finite(y)], ylim)
  }else{
    y.range <- c(-Rangexyz[2],Rangexyz[2])
  }
  y.prty <- pretty(y.range, n = lab[2], min.n = max(1, min(0.5 *lab[2], p.lab[2])))
  y.scal <- round(diff(y.prty[1:2]), digits = 12)
  y.add <- min(y.prty)
  y <- (y - y.add)/y.scal
  y.max <- (max(y.prty) - y.add)/y.scal
  if (!is.null(ylim))
    y.max <- max(y.max, ceiling((ylim[2] - y.add)/y.scal))
  
  # Z axis range
  if (is.null(Rangexyz)){
    z.range <- range(z[is.finite(z)], zlim)
  }else{
    z.range <- c(-Rangexyz[3],Rangexyz[3])
  }
  z.prty <- pretty(z.range, n = lab.z, min.n = max(1, min(0.5 *lab.z, p.lab[2])))
  z.scal <- round(diff(z.prty[1:2]), digits = 12)
  z <- z/z.scal
  z.range <- range(z.prty)/z.scal
  z.max <- ceiling(z.range[2])
  z.min <- floor(z.range[1])
  if (!is.null(zlim)) {
    z.max <- max(z.max, ceiling(zlim[2]/z.scal))
    z.min <- min(z.min, floor(zlim[1]/z.scal))
  }
  z.range <- range(z.min, z.max)
  
  # Add grid
  if ("xy" %in% grid || grid == TRUE) {
    i <- x.min:x.max
    segments(i, z.min, i + (yx.f * y.max), yz.f * y.max + 
               z.min, col = col.grid, lty = lty.grid)
    i <- 0:y.max
    segments(x.min + (i * yx.f), i * yz.f + z.min, x.max + 
               (i * yx.f), i * yz.f + z.min, col = col.grid, lty = lty.grid)
  }
  
  if ("xz" %in% grid) {
    i <- x.min:x.max
    segments(i + (yx.f * y.max), yz.f * y.max + z.min, 
             i + (yx.f * y.max), yz.f * y.max + z.max, 
             col = col.grid, lty = lty.grid)
    temp <- yx.f * y.max
    temp1 <- yz.f * y.max
    i <- z.min:z.max
    segments(x.min + temp,temp1 + i, 
             x.max + temp,temp1 + i , col = col.grid, lty = lty.grid)
    
  }
  
  if ("yz" %in% grid) {
    i <- 0:y.max
    segments(x.min + (i * yx.f), i * yz.f + z.min,  
             x.min + (i * yx.f) ,i * yz.f + z.max,  
             col = col.grid, lty = lty.grid)
    temp <- yx.f * y.max
    temp1 <- yz.f * y.max
    i <- z.min:z.max
    segments(x.min + temp,temp1 + i, 
             x.min, i , col = col.grid, lty = lty.grid)
  }
}
###############################################################################
# 基础拼图
Range <- pc12 %>% .[,2:4] %>% map(summary) %>% enframe %>% 
  mutate(value = map(value,function(x) x %>% as.matrix %>% t %>% as_tibble)) %>%
  unnest(value) %>% select(name,Min.,Max.) %>% mutate(abs_min = abs(Min.),abs_max = abs(Max.)) %>%
  mutate(Mmax = sapply(1:nrow(.),function(x) max(.[x,-1]))) %>% .$Mmax

Range <- Range + Range*0.01

############## rgl ######################
pdf(str_c(opts$mode,"_",opts$input,"_3D.pdf"),width=6,height = 6)
if (opts$range){
  s3d <- scatterplot3d(pc12[,2:4], pch = pc12$shape, color = pc12$Color,#type = "h",
                       xlab=xlab,ylab=ylab,zlab=zlab,xlim=c(-Range[1],Range[1]),
                       ylim=c(-Range[2],Range[2]),zlim=c(-Range[3],Range[3]), box=F)
  addgrids3d(pc12[,2:4], grid = c("xy", "xz", "yz"),Rangexyz = Range)
  s3d$points3d(pc12[,2:4],pch = pc12$shape, col = pc12$Color)
}else{

  s3d <- scatterplot3d(pc12[,2:4], pch = pc12$shape, color = pc12$Color,#type = "h",
                       xlab=xlab,ylab=ylab,zlab=zlab, box=F)

  addgrids3d(pc12[,2:4], grid = c("xy", "xz", "yz"))
  s3d$points3d(pc12[,2:4],pch = pc12$shape, col = pc12$Color)

}
legend("bottom", legend = levels(sp$group),
       col = unique(sp$Color), pch = unique(sp$shape), 
       inset = -0.24, xpd = TRUE, horiz = TRUE)
dev.off()
###################################
## scatter3d(x = pc12$PC1, y = pc12$PC2, z = pc12$PC3, groups = pc12$group,
##           surface=FALSE, grid = FALSE, ellipsoid = TRUE,xlab=xlab,ylab=ylab,zlab=zlab,
##           surface.col = unique(sp$Color),axis.col = c("black", "black", "black"))
## rgl.postscript(str_c(opts$mode,"_",opts$input,"_3D-typeII.pdf"),fmt="pdf")
## rgl.snapshot(str_c(opts$mode,"_",opts$input,"_3D-typeII.png"),filename = "plot.png")

p <- GGplot8(pc12,Range)
Sys.setenv("plotly_username"="Colin_Liu_92")
Sys.setenv("plotly_api_key"="f3Yw9vnrm7GTTtGmhBZf")
#plotly_IMAGE(p, width = 1000, height = 1000, format = "png", scale = 1,
#             out_file = str_c(opts$mode,"_",opts$input,"_3D.png"))
htmlwidgets::saveWidget(p, file = str_c(opts$mode,"_",opts$input,"_3D.html"))
