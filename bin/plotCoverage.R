## Check if necessary packages are installed
packages <- c('ggplot2','egg','patchwork','dplyr', "argparse","rvg","officer","ggpubr")
install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)
invisible(suppressPackageStartupMessages(lapply(packages, library, character.only = TRUE)))
parser <- ArgumentParser(description='Produce Coverage Plot',add_help=T)
parser$add_argument('-c', "--csv",help="Provide CSV file containing the path of mosdepth and feature file", default=NULL,nargs=1)
parser$add_argument('-w', "--wtPrefix",help="Prefix used to represent Wild Type samples", default="Trt",nargs=1)
parser$add_argument('-t', "--trtPrefix",help="Prefix used to represent Treatment samples", default="Ctrl",nargs=1)
parser$add_argument('-x', "--xCoverage",help="Coverage threshold for plots", default=20,nargs=1)
parser$add_argument('-p', "--pname",help="Name of the presentation", default="output",nargs=1)

args <- parser$parse_args()
file <- read.csv(args$csv)
ctrlpfx <- args$wtPrefix
trtpfx <- args$trtPrefix
coverageTh <- as.numeric(args$xCoverage)
## Make full paths and remove file name column
file[,1] <- paste0(file[,1],"/",file[,2])
file[,3] <- paste0(file[,3],"/",file[,4])
file <- file[,c(1,3,4,5)]


## We wish to split the dataframe w.r.t feature files since there will be control and treatment per feature

split_df <- file %>% group_by_at(3) %>% group_split()
#names(split_df) <- gsub(".feature.tsv","",unique(file$File2))

plotCov <- function(ctrl,trt,feature,geneID){
  ctrl <- read.csv(ctrl,header = F, sep = '\t')
  trt <- read.csv(trt,header = F, sep = '\t')
  feature <- read.table(feature,sep = "\t")
  feature$gene <- geneID
  feature[,1] <- trimws(feature[,1])
  plt.ctrl <- ggplot(ctrl) +
    geom_rect(aes(xmin = V2, xmax = V3, ymin = 0, ymax = V4), col="blue", fill="blue")+ 
    coord_cartesian(ylim = c(0, coverageTh),xlim = c(min(ctrl$V2),max(ctrl$V3)), expand = FALSE)+theme_minimal()+
    theme(axis.ticks.x=element_blank(),axis.text.x=element_blank())+ylab(ctrlpfx)
  plt.trt <- ggplot(trt) +
    geom_rect(aes(xmin = V2, xmax = V3, ymin = 0, ymax = V4), col="maroon", fill="maroon")+ 
    coord_cartesian(ylim = c(0, coverageTh),xlim = c(min(ctrl$V2),max(ctrl$V3)), expand = FALSE)+theme_minimal()+
    theme(axis.ticks.x=element_blank(),axis.text.x=element_blank())+ylab(trtpfx)

  ## Took inspiration from: https://zhiganglu.com/post/ggplot-visualising-gene-features/
  plt.feature <- ggplot(data = feature, aes(gene, max(V3),width = 1))+geom_bar(stat="identity", fill="white", width = 0.05)+coord_flip(expand = FALSE)+labs(x="Gene", y="Position")+ylim(0,max(feature$V3))+geom_segment(data=feature, aes(x=gene, xend=gene, y=V2, yend=V3, color=V1), size=5)+theme_classic()+labs(color="Gene")+theme(aspect.ratio = 0.1/1,axis.ticks.y = element_blank(),axis.line.y = element_blank(),axis.line = element_line(size = 1),legend.position = 'bottom')+labs(colour="Features")
  #grid::grid.draw(egg::ggarrange(plots=list(plt.ctrl, plt.trt,plt.feature)))
 egg::ggarrange(plots=list(plt.ctrl, plt.trt,plt.feature))
  # plot_grid(plt.ctrl, plt.trt, plt.feature,
  #           ncol = 1, align = "v")
}


pptx <- read_pptx()
for (i in seq(length(split_df))) {
  pptx<- add_slide(pptx,layout = "Title and Content", master = "Office Theme")
}

plots <- lapply(1:length(split_df), function(x){
  controlSample <- as.character(subset(split_df[[x]], split_df[[x]][,4]==ctrlpfx)[,1])
  trtSample <- as.character(subset(split_df[[x]], split_df[[x]][,4]==trtpfx)[,1])
  featurefile <- as.character(unique(split_df[[x]][,2]))
  temp <- plotCov(ctrl=controlSample,trt=trtSample,feature=featurefile,geneID=unique(gsub(".feature.tsv","",split_df[[x]][,3] %>% pull())))
  #captured_plot <- recordPlot()
  pptx %>% 
    #add_slide(layout = "Title and Content", master = "Office Theme") %>% 
    on_slide(index = x) %>% 
  ph_with(value = ggpubr::as_ggplot(temp), location = ph_location_type(type = "body")) %>%
    # ph_with(value = temp, location = ph_location_type(type = "body")) %>%
    print(target = paste0(args$pname,".pptx"))
})

#saveRDS(plots,"plots.rds")