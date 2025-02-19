# load required packages
library(CATALYST)
library(cowplot)
library(flowCore)
library(ggplot2)
library(dittoSeq)
library(dplyr)
library(reshape2)
library(SingleCellExperiment)
source('fun.R')
# Load Data
sce <- readRDS('sce.rds')

### Fig1B
markernames <- c("PanCK","CK8","FAP","aSMA","CD31",
                 "CD45","CD19","CD27","CD28",
                 "CD3","CD4","CD8","CD45RO","CD25","TCRgd","CCR7","CXCR3","CCR6","CCR4","CD103",
                 "CD56",
                 "CD11b","CD11c","CD14","CD68","CD16",
                 "CD66b","CRTH2","FceRI","CD123","HLA_DR",
                 "Ki67","PD_L1","PD_1","acetylated","lactyllysine","p_tyrosine","ubiquitin")


mycolor1 <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", 
              "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", 
              "#33A02C", "#A6761D","#B2DF8A", "#55A1B1", "#8DD3C7",  
              "#E6AB02", "#7570B3", "#BEAED4","#6b473c", "#999999", 
              "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000", 
              "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00",
              "#abc88b","#9b95c9","#ac6767","#dbce8f","#666666")

mycolor2 <- c("#DC0000FF","#6950a1","#4DBBD5FF","#00A087FF",
              "#7E6148FF","#F39B7FFF","#91D1C2FF","#8491B4FF")

###################################################
row_title_gp.col <- setNames(mycolor2[1:length(levels(sce$celltype))],
                             levels(sce$celltype))[c("T",
                                                     "NK","B",
                                                     "Myeloid",'Other',
                                                     "Endothelial",
                                                     "Fibroblast")]

p <- plotExprComplexHeatmap.v2(object = sce,
                               features = markernames,
                               metacluster = 'cluster_id',
                               het_col =cmap(75),
                               cluster_rows = T,
                               ann_bar_form ='Celltypes',
                               gap_var= 'celltype',
                               row_title_gp.col=row_title_gp.col,
                               ann_bar_col = mycolor1)
dev.off()
pdf('Fig1B.pdf',he=8,wi=15)
print(p)
dev.off()

### Fig1C
p <- dittoDimPlot(sce,
                  labels.highlight = F,
                  var = "celltype", 
                  reduction.use = "tSNE", 
                  labels.size = 5,
                  do.label = TRUE) + 
  scale_color_manual("Celltype",values = mycolor1) +
  ggtitle("Celltype on tSNE") + 
  tidydr::theme_dr() + 
  theme(panel.grid =element_blank(), 
        axis.title = element_text(hjust = 0.08,size = 18),
        plot.title = element_text(size = 18,hjust = 0.5),
        axis.line =element_line(size = 0.5),
        legend.position ='none')


ggsave(filename = "Fig1C.pdf",plot = p,width = 6,height = 6)

### Fig1D
p <- dittoDimPlot(sce,
                  labels.highlight = F,
                  var = "Celltypes", 
                  reduction.use = "tSNE", 
                  labels.size = 5,
                  do.label = TRUE) + 
  scale_color_manual("Celltype",values = mycolor1) +
  ggtitle("Celltype on tSNE") + 
  tidydr::theme_dr() + 
  theme(panel.grid =element_blank(), 
        axis.title = element_text(hjust = 0.08,size = 18),
        plot.title = element_text(size = 18,hjust = 0.5),
        axis.line =element_line(size = 0.5),
        legend.position ='none')


ggsave(filename = "Fig1D.pdf",plot = p,width = 6,height = 6)

### Fig1E
library(viridisLite)

plot_list <- multi_dittoDimPlot(sce.down, 
                                var = c("acetylated","lactyllysine","p_tyrosine","ubiquitin"),
                                assay = "exprs",
                                reduction.use = "tSNE",
                                list.out = TRUE,
                                legend.show =T,
                                split.by ='Group')  

plot_list <- lapply(plot_list, function(x) x + 
                      scale_color_gradientn(colors = inferno(10)) +
                      theme(strip.background = element_blank(),
                            plot.title = element_text(size = 18),
                            strip.text = element_text(size = 16)))

plot_list <- lapply(plot_list, function(x) x + theme(plot.title = element_text(hjust = 0.5)))

p<-  patchwork::wrap_plots(plotlist = plot_list) + plot_layout(guides = 'collect')

ggsave(filename = "Fig1E.pdf",plot = p,width = 14,height = 8)

### Fig1F
p<- Dotplot_compare.V2(object=sce,
                       assays='data',
                       vars='celltype',
                       features = c("Acetylation","Lactylation","Phosphorylation","Ubiquitination"),
                       summarise_fun=c('mean'),
                       method=c('t.test'),
                       group.by=c('sample_id','Group'),
                       colours_bar=colorRampPalette(c("#426ab3","#FFFFB3","#f15a22"))(256))           

ggsave('Fig1F.pdf',plot = p,width = 6,height = 4.5)

### Fig1G
p<- Dotplot_compare.V2(object=sce,
                       assays='data',
                       vars='Celltypes',
                       features = c("Acetylation","Lactylation","Phosphorylation","Ubiquitination"),
                       summarise_fun=c('mean'),
                       method=c('t.test'),
                       group.by=c('sample_id','Group'),
                       colours_bar=colorRampPalette(c("#426ab3","#FFFFB3","#f15a22"))(256))+
  coord_flip()

ggsave('Fig1G.pdf',plot = p,width =13,height = 4.5)


