# Rscript plor_alpha_box_dot.r -a alpha_rarefac.summary.xls -g map-group.txt
# 
library(optparse)
library(tidyverse)
library(magrittr)
library(reshape2)
library(ggplot2)
if (TRUE){
  option_list <- list(
    make_option(c("-a", "--alpha"), type="character", default="alpha_rarefac.summary.xls",
                help="输入的OTU表格"),
    make_option(c("-g", "--map"), type="character", default="map-group.txt",
                help="分组文件"),
    make_option(c("-c", "--color"), type="character", default="none",
                help="指定分组颜色:color.txt"),
    make_option(c("-v", "--CI"), type="double", default=0.95,
                help="置信区间")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
}

# Check if required packages are installed ----
packages <- c("cowplot", "readr", "ggplot2", "dplyr", "lavaan", "smooth", "Hmisc")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))
}

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(
                ymin = min(y),
                ymax = max(y),
                xmin = x,
                xmax = x + width / 2
              )
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data,
                              xminv = x,
                              xmaxv = x + violinwidth * (xmax - x)
            )
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(
              plyr::arrange(transform(data, x = xminv), y),
              plyr::arrange(transform(data, x = xmaxv), -y)
            )
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1, ])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(
            weight = 1, colour = "grey20", fill = "white", size = 0.5,
            alpha = NA, linetype = "solid"
          ),
          
          required_aes = c("x", "y")
  )
##################################################################################################
mycol <-c("#CD3333","#483D8B","#458B00","#EEC900","#BFBC3B","#6E8B3D","#00688B","#C10077","#CAAA76","#EEEE00","#458B00","#8B4513","#008B8B","#6E8B3D","#8B7D6B","#7FFF00","#CDBA96","#ADFF2F")

a <- read.table(opts$alpha,header=T,sep="\t")
b <- read.table(opts$map,header=T,sep="\t",comment.char = "")
b[,2] <- factor(b[,2],levels=unique(b[,2]))
b$color <- mycol[as.numeric(b[,2])]

#如果有颜色定义的话，就可以直接定义了
if(opts$color != "none"){
  sc <- read.table(opts$color,sep="\t",comment.char = "",check.names = FALSE)
  sc <- sc[which(sc[,1] %in% unique(b[,2])),]
  mycol <- as.vector(sc[,2])
  b[,2] <- factor(b[,2],levels=as.vector(sc[,1]))
  b[,3] <- mycol[as.numeric(b[,2])]
}

colnames(b)[1] <- colnames(a)[1]
a <- inner_join(a,b) %>% mutate(group = fct_inorder(group))

#添加表头
name1<-c()
for (g in as.vector(unique(b$group))){
  name1<-c(name1,paste(g,'Mean',sep='-'),paste(g,'SE',sep='-'))
}
name1 <- c(name1,'P-value')

GGBOX <- function(data,jitter="TRUE",Mode="box",CI){
  # data <- a[,c(9,10,11)];CI <- opts$CI;Mode <- "cloud"
  d <- melt(data) %>% group_by(group) %>% 
    summarize_at(vars(value),list(mean = ~mean(.), count = ~ n(), sd = ~sd(.), se = ~ sd/sqrt(count),
                                  ci = ~ se * (qt(CI/2 + .5, count-1))))
  if (Mode == "box"){
    p <- ggplot(melt(data),aes(x=group,y=value)) + 
      stat_boxplot(geom ='errorbar', width = 0.6)
    
    if (jitter == "TRUE"){
      p <- p + geom_boxplot(notch = F,aes(color = group),outlier.size = 0,outlier.shape = NA) +
        geom_point(aes(color = group),position = position_jitterdodge())
    }else{
      p <- p + geom_boxplot(notch = F,aes(color = group))
    }
    
    #if (jitter == "TRUE"){p <- p + geom_point(aes(color = group),position="jitter")}
  }else if (Mode == "dot"){
    p <- ggplot(melt(data),aes(x=group,y=value,fill=group)) + 
      geom_dotplot(binaxis='y', stackdir='center',dotsize = 0.6) +
      stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                   geom="errorbar", color="black", width=0.2, show.legend = FALSE) +
      stat_summary(fun.y="mean", geom="point", color="black", show.legend = FALSE)
  }else if (Mode == "violin"){
    p <- ggplot(melt(data),aes(x=group,y=value,fill=group)) + 
      geom_violin(alpha=.5,show.legend = FALSE) +
      geom_boxplot(width=.1) +
      scale_fill_manual(values=unique(data$color))+
      geom_jitter(show.legend = FALSE)
  }else if (Mode == "cloud"){
    p <- ggplot(melt(data),aes(x=group,y=value,fill=group)) +
      geom_flat_violin(position = position_nudge(x=.1,y = 0), 
                       adjust = 1.5, trim = TRUE, alpha = .3, colour = NA) +
      geom_point(aes(x = as.numeric(group)-.15, y = value, colour = group),
                 position = position_jitter(width = .05), size = 1, shape = 16,show.legend = FALSE) +
      geom_boxplot(aes(x = group, y = value, fill = group),#outlier.shape = NA,
                   alpha = .5, width = .1, colour = "black",show.legend = FALSE) +
      geom_point(data = d, shape = 16,size = 3,show.legend = FALSE, 
                 aes(x = as.numeric(group)+.1, y = mean, group = group, color = group)) +
      geom_errorbar(data = d, aes(x = as.numeric(group)+.1, y = mean, 
                                  group = group, colour = group, ymin = mean-se, 
                                  ymax = mean+se), width = .05,show.legend = FALSE) +
      coord_flip()
  }
  
  p <- p + scale_color_manual(values = unique(data$color), guide = FALSE) + 
    scale_fill_manual(values = unique(data$color)) + 
    xlab("") +
    ylab("") +
    ggtitle(colnames(data)[1]) +
    theme_bw() + ##设置主题
    theme(plot.title=element_text(size=rel(1),hjust=0.5),
          plot.margin = unit(c(1, 1, 1, 1), "lines"),
          axis.title=element_text(size=25,face="bold"),
          axis.text.x=element_text(hjust =1),#angle=45,
          axis.text.y=element_text(angle=45),
          axis.text=element_text(size=18,face="bold"),
          title= element_text(size=15,face= "bold", vjust=0.5, hjust=0.5))
  
  if (Mode != "cloud"){
    Width <- ifelse(log2((nlevels(data$group))^(1/2))<=1,1,log2((nlevels(data$group))^(1/2)))*6
    Height <- 6
  }else{
    Height <- log2(nlevels(data$group))*4
    Width <- 8
  }
  
  ggsave(paste0(colnames(data)[1],"-",Mode,"plot.pdf"),p,width = Width,height = Height)
}

for (i in c(4,5,6,7,9)){#i <- c(4,5,6,7,9)[1]
  #GGBOX(a[,c(as.numeric(i),10,11)])
  GGBOX(a[,c(as.numeric(i),10,11)],Mode="cloud",CI=opts$CI)
}
