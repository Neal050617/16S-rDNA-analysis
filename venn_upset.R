# Rscript venn_upset.R -i rarefac.otu_table.xls -g map-group.txt -m all
# 排列组合，得到所有可能的分组结果
# 20190313-修改图像长宽;排序规则
library(optparse)
library(tidyverse)
library(gridExtra)
library(grid)
library(ggforce)
library(VennDiagram)
library(venneuler)#rjava更新
library(magrittr)

if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="rarefac.otu_genus.xls",
                help="输入的OTU表格"),
    make_option(c("-g", "--map"), type="character", default="map-group.txt",
                help="分组文件"),
    make_option(c("-m", "--mode"), type="character", default="all",
                help="venn <= 5组; upset; all"),
    make_option(c("-c", "--color"), type="character", default="none",
                help="指定分组颜色:color.txt")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
  print(paste("The input file is ", opts$input,  sep = ""))
  print(paste("The group file is ", opts$map,  sep = ""))
  print(paste("The mode is ", opts$mode,  sep = ""))
  print(paste("The color file is ", opts$color,  sep = ""))
}

######## 函数 ############################################################################
# 制作真值表
make.truth.table <- function(mp2){
  mp3 <- vector(length(mp2),mode="list")
  for (i in 1:length(mp2)){mp3[[i]] <- c(TRUE,FALSE)}
  vend <- expand.grid(mp3) %>% .[-nrow(.),]
  for (i in 2:nrow(vend)){# i = 2
    vv <- as.vector(unlist(vend[i,]))
    as.matrix(colnames(vend))[vv,]
  }
  names(vend) <- names(mp2)
  return(vend)
}
# 得到venn图需要的数据
get.venn.table <- function (x,keep.order = F) { # x = mp2
  stopifnot(typeof(x) == "list")
  emptyInds <- unlist(lapply(x, is.null))
  if (any(emptyInds)) {
    warning("removing NULL elements in list.")
    x <- x[!emptyInds]
  }
  out <- make.truth.table(x) # 真值表
  # 列好组合名
  setNames <- apply(out, 1, function(categories) { # categories <- as.matrix(out[2,])
    include <- paste(names(x)[categories], collapse = "∩")
    if (all(categories)) {
      return(include)
    }
    include <- paste0("(", include, ")")
    exclude <- paste0("(", paste(names(x)[!categories], collapse = "∪"), ")")
    paste(include, exclude, sep = "&")
  })
  # 找不同~啦啦啦
  setValues <- apply(out, 1, function(categories) {
    include <- Reduce(intersect, x[categories])
    exclude <- Reduce(union, x[!categories])
    setdiff(include, exclude)
  })
  # upsetname
  #setupset <- apply(out, 1, function(categories) {paste(names(x)[categories], collapse = "&")})
  setEle <- sapply(setValues,function(x) paste0(x,collapse = ","))
  setNum <- unlist(lapply(setValues, length))
  
  out <- cbind(out, setNames, setNum, setEle)
  colnames(out)[(ncol(out) - 2):(ncol(out))] <- c("set","count","values")
  out$set <- as.character(out$set)
  Encoding(out$set) <- "UTF-8"
  
  # 生成表格
  #out2 <- as.vector(setNum)
  #names(out2) <- as.character(setupset)
  #Encoding(names(out2)) <- "UTF-8"
  #return(list(out,out2))
  return(out)
}
# 主图柱形图
make_main_bar  <- function (Main_bar_data, show_num = "Yes"){# Main_bar_data = out2
  bottom_margin <- -0.65
  ymax <- max(Main_bar_data$freq) + (max(Main_bar_data$freq)) * 0.1
  
  Main_bar_plot <- (ggplot(data = Main_bar_data, aes_string(x = "x",y = "freq")) + 
                      ylim(0, ymax) + 
                      geom_bar(stat = "identity", width = 0.6,fill = Main_bar_data$color) + 
                      scale_x_continuous(limits = c(0,(nrow(Main_bar_data) + 1)), 
                                         expand = c(0, 0), breaks = NULL) + 
                      xlab(NULL) + 
                      ylab("Intersection Size") + 
                      labs(title = NULL) + 
                      geom_vline(xintercept = 0,color = "gray0") + 
                      geom_hline(yintercept = 0, color = "gray0")) +
    theme(panel.background = element_rect(fill = "white"),
          plot.margin = unit(c(0.5, 0.5, bottom_margin, 0.5), "lines"),
          panel.border = element_blank(), 
          axis.title.y = element_text(vjust = -0.8,size = 8.3), 
          axis.text.y = element_text(vjust = 0.3,size = 7 ))
  
  if ((show_num == "yes") || (show_num == "Yes")) {
    Main_bar_plot <- (Main_bar_plot + geom_text(aes_string(label = "freq"), 
                                                size = 3 , vjust = -1, 
                                                angle = 0, colour = Main_bar_data$color))
  }
  
  Main_bar_plot <- ggplotGrob(Main_bar_plot)
  return(Main_bar_plot)
}
# 矩阵点线图
make_matrix_plot <- function(Mat_data, shading_data, labels){
  # Mat_data = out3;shading_data = out4;labels = colnames(out2_1)[1:ncol(out)]
  Matrix_plot <- (ggplot() + 
                    geom_rect(data = shading_data, aes_string(xmin = "min", xmax = "max",
                                                              ymin = "y_min", ymax = "y_max"),
                              fill = shading_data$shade_color, alpha = 0.25) +
                    geom_line(data= Mat_data, aes_string(group = "Intersection", x="x", y="y",
                                                         colour = "line_col"), size = 0.7) +
                    geom_point(data= Mat_data, aes_string(x= "x", y= "y"), colour = Mat_data$color,
                               size= 5, alpha = Mat_data$alpha, shape=16) +
                    scale_color_identity() +
                    xlab(NULL) + ylab("   ") + 
                    scale_y_continuous(breaks = c(1:length(unique(Mat_data$y))),
                                       limits = c(0.5,(length(unique(Mat_data$y)) +0.5)),
                                       labels = labels, expand = c(0,0)) + 
                    scale_x_continuous(limits = c(0,Mat_data$x[nrow(Mat_data)]+1), expand = c(0,0)) +
                    theme(panel.background = element_rect(fill = "white"),
                          plot.margin=unit(c(-0.2,0.5,0.5,0.5), "lines"),
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.text.y = element_text(colour = "gray0",
                                                     size = 7, hjust = 0.4),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())
  )
  
  Matrix_plot <- ggplotGrob(Matrix_plot)
  return(Matrix_plot)
}
# 左下独立组内OTU数
make_size_plot <- function(Set_size_data, sbar_color){# Set_size_data = out2_3;sbar_color = mycol
  Size_plot <- (ggplot(data = Set_size_data, aes_string(x = "x",y = "y")) + 
                  geom_bar(stat = "identity", colour = sbar_color, width = 0.4, 
                           fill = sbar_color, position = "identity") + 
                  scale_x_continuous(limits = c(0.5, (nrow(Set_size_data) + 0.5)), 
                                     breaks = c(0, max(Set_size_data)), expand = c(0,0)) + 
                  theme(panel.background = element_rect(fill = "white"),
                        plot.margin = unit(c(-0.11, -1.3, 0.5, 0.5), "lines"),
                        axis.title.x = element_text(size = 8.3),
                        axis.text.x = element_text(size = 7, angle = 0, vjust = 1, hjust = 0.5),
                        axis.line = element_line(colour = "gray0"), axis.line.y = element_blank(), 
                        axis.line.x = element_line(colour = "gray0", size = 0.3),
                        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
                  xlab(NULL) + ylab("Set Size") + coord_flip() + 
                  scale_y_continuous(trans = "reverse"))
  
  Size_plot <- ggplotGrob(Size_plot)
  return(Size_plot)
}
# ggplot2-venn
make_venn_plot <- function(data,matrix.col=mycol){# data = out1
  gvd1 <- venneuler(data)
  gvd2 <- data.frame(gvd1$centers,diameters = gvd1$diameters,labels = gvd1$labels,stringsAsFactors = FALSE)
  gvd2$labels <- factor(gvd2$labels,levels=unique(gvd2$labels))
  GV <- ggplot(gvd2) +
    geom_circle(aes_(x0 = ~x, y0 = ~y,r = ~diameters/2, fill = ~labels ,color = "white"),alpha = 0.5) +
    #geom_label_repel(aes_(x = ~x, y = ~y, label = ~labels)) +
    #guides(fill=FALSE) + 
    scale_color_manual(values = "white", guide = FALSE) +
    scale_fill_manual(values = matrix.col) +
    xlab(NULL) + ylab(NULL) + labs(title = NULL) + coord_fixed() +
    theme(panel.background = element_rect(fill = "white"),
          axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.title = element_blank(),legend.spacing.x = unit(0.2, 'cm')
    )
  
  make_venn_plot <- ggplotGrob(GV)
}
# 拼图
get.upset.table <- function (x,matrix.col=mycol,sort_by="combn",hratios=0.7) {#x=mp2;matrix.col=mycol[1:length(gp1)];0.7
  stopifnot(typeof(x) == "list")
  emptyInds <- unlist(lapply(x, is.null))
  if (any(emptyInds)) {
    warning("removing NULL elements in list.")
    x <- x[!emptyInds]
  }
  out <- make.truth.table(x) # 真值表
  # 找不同~啦啦啦
  setValues <- apply(out, 1, function(categories) {
    include <- Reduce(intersect, x[categories])
    exclude <- Reduce(union, x[!categories])
    setdiff(include, exclude)
  })
  setupset <- apply(out, 1, function(categories) {paste(names(x)[categories], collapse = "&")})
  setNum <- unlist(lapply(setValues, length))
  # 生成"&"分隔的统计表格
  out1 <- as.vector(setNum)
  names(out1) <- as.character(setupset)
  Encoding(names(out1)) <- "UTF-8"
  # upset-matrix右下方图的数据
  out2 <- as.data.frame(apply(out,2,as.numeric)) %>% 
    mutate(freq = setNum) %>% 
    filter(freq > 0) %>% 
    mutate(x = seq(nrow(.),1,by = -1),color = "gray23")
  # 排序  
  if (sort_by=="freq"){
    out2 %<>% arrange(desc(freq))
  }else{
    print("sort by combn")
  }
  # 修改行名
  rownames(out2) <- out2$x
  out2$x <- seq(1, nrow(out2))
  
  out2_1 <- out2
  # 补全组名
  colnames(out2_1)[1:ncol(out)] <- unlist(lapply(1:ncol(out),function(p) {
    ifelse(nchar(colnames(out2)[p])<7,
           paste0(paste0(replicate(7-nchar(colnames(out2)[p]), " "),collapse = ""),colnames(out2)[p]),
           colnames(out2)[p]) # p=1
  }))
  out2_2 <- t(out2_1[,1:ncol(out)])
  # each groups' OTU counts
  out2_3 <- sapply(1:ncol(out),function(g) out2[,g]*out2[,ncol(out)+1]) %>% 
    colSums(.) %>% data.frame(y=.,x=seq(ncol(out)))
  # matrix data
  out3 <- data.frame(expand.grid(y = seq(nrow(out2_2)), x = seq(ncol(out2_2))), 
                     value = as.vector(out2_2)) %>% 
    mutate(color=ifelse(value>0,matrix.col[y],"gray83")) %>%
    mutate(alpha=ifelse(value>0,1,0.5)) %>%
    mutate(Intersection=ifelse(value>0,paste(x, "yes", sep = ""),
                               paste(rownames(.), "No", sep = ""))) %>%
    mutate(line_col='black')
  
  # matrix shade
  out4 <- unique(out3$y) %>% .[which(.%%2 != 0)] %>% 
    data.frame(y=.) %>% mutate(min=0,max=max(out3$x)+1) %>% 
    mutate(y_min=y-0.5,y_max=y+0.5,shade_color="gray88")
  
  Main_bar_plot <- make_main_bar(out2)
  Matrix_plot <- make_matrix_plot(out3,out4,colnames(out2_1)[1:ncol(out)])
  Size_plot <- make_size_plot(out2_3,matrix.col)
  Venn_plot <- make_venn_plot(out1,matrix.col)
  
  Main_bar_plot$widths <- Matrix_plot$widths
  Matrix_plot$heights <- Size_plot$heights
  size_plot_height <- ((hratios + 0.01) * 100) # 根据比例设置大图的高度0.7/0.3
  
  Height <- ifelse(100*log2((max(out3$x))^(1/3))<=100,100,100*log2((max(out3$x))^(1/3)))
  
  pdf(paste0(paste(colnames(out),collapse = "-"),"_sort_by_",sort_by,"_upset.pdf"), height = 8, width = 8 * Height/100)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(100, Height)))
  
  vp = viewport(layout.pos.row = 1:100, layout.pos.col = 21:Height)
  pushViewport(vp)
  grid.draw(arrangeGrob(Main_bar_plot, Matrix_plot, heights = c(hratios,1-hratios)))
  popViewport()
  
  vp = viewport(layout.pos.row = size_plot_height:100, layout.pos.col = 1:20)
  pushViewport(vp)
  grid.draw(arrangeGrob(Size_plot))
  popViewport()
  
  vp = viewport(layout.pos.row = 1:50,layout.pos.col = (Height-50):Height)
  pushViewport(vp)
  grid.draw(arrangeGrob(Venn_plot))
  popViewport()
  
  dev.off()
  #return(out2)
}

######################################################################################################

mycol <-c("#CD3333","#483D8B","#458B00","#EEC900","#BFBC3B","#6E8B3D","#00688B","#C10077","#CAAA76","#EEEE00","#458B00","#8B4513","#008B8B","#6E8B3D","#8B7D6B","#7FFF00","#CDBA96","#ADFF2F")

# 读取分组文件
mp1 <- read.table(opts$map,head= T ,sep="\t",comment.char = "",fileEncoding = "UTF-8",stringsAsFactors = F)
colnames(mp1)[1] <- "SampleID"
rownames(mp1) <- mp1[,1]
mp1$group <- fct_inorder(mp1$group)

# 判断是否存在颜色文件
if (opts$color != "none"){
  cp1 <- read.table(opts$color,head= F ,sep="\t",comment.char = "",fileEncoding = "UTF-8")
  rownames(cp1) <- cp1[,1]
  cp2 <- as.vector(cp1[which(cp1[,1] %in% mp1[,2]),1])
  mycol <- as.vector(cp1[which(cp1[,1] %in% mp1[,2]),2])
  mp1 <- mp1[unlist(lapply(1:length(cp2),function(x) which(mp1[,2] %in% cp2[x]))),]
}

# 读取文件
df1 <- read.table(opts$i,head= T ,sep="\t",comment.char = "",row.names = 1)
df1 <- df1[,rownames(mp1)] %>% mutate(SUM = rowSums(.),names = rownames(.)) %>% 
  filter(SUM>0) %>% select(-SUM)
rownames(df1) <- df1$names
df1 %<>% select(-names)
# 转化为二进制矩阵
df2 <- ifelse(df1>0, 1, 0)

# 根据group文件排序；真值表转化为NA或行名
df2 <- df2[,rownames(mp1)]
df3 <- sapply(1:ncol(df2),function(x) ifelse(df2[,x]>0,rownames(df2),NA)) 
colnames(df3) <- colnames(df2)

# venn图的输入数据
gp1 <- as.vector(unique(mp1$group))
mp2 <- vector(nlevels(mp1$group),mode="list")
for (i in 1:length(gp1)){# i = 1
  mp2[[i]] <- unique(as.vector(na.omit(as.vector(df3[,mp1[mp1$group %in% gp1[i],1]]))))
}
names(mp2) <- gp1
f <- get.venn.table(mp2)
write.table(f,paste0(paste(gp1,collapse = "-"),"_venn_diagramm.txt"),row.names=F,col.names=T,quote=F, sep="\t")

if ((opts$mode == "venn" | opts$mode == "all") & length(mp2) <= 5){
  venn.diagram(mp2,
               filename = paste0(paste(gp1,collapse = "-"),"_venn_diagramm.tiff"),
               category.names = gp1,  #名字
               col = "black",
               fill = mycol[1:length(gp1)],  #颜色
               cat.col = mycol[1:length(gp1)],
               margin=0.05, main.cex = 2,  sub.cex = 1.3,
               cat.fontface = 2,
               alpha = 0.50,    # 透明度
               height = 480 , width = 480 , 
               resolution = 300, #清晰度
               lwd = 2,lty = 'blank',
               cex = 0.5,
               #fontface = "bold",fontfamily = "sans",
               cat.cex = 0.6,
               cat.default.pos = "outer",scaled = FALSE,inverted = FALSE,reverse = FALSE
  )
  file.remove(dir(path=".",pattern='.log$')) # 删除log文件
  
}

if (opts$mode == "upset" | opts$mode == "all"){
  get.upset.table(mp2,mycol[1:length(gp1)],sort_by="freq",0.7)
  get.upset.table(mp2,mycol[1:length(gp1)],sort_by="combn",0.7)
}

# library(UpSetR)
# upset(fromExpression(u),order.by="freq",nsets=length(gp1),nintersects=NA,keep.order=T,
#       sets = gp1,
#       number.angles = 0, point.size = 6, 
#       line.size = 1,text.scale = 2,
#       sets.bar.color = mycol[1:length(gp1)]
#       #matrix.color = mycol[1:length(gp1)]
#       )


