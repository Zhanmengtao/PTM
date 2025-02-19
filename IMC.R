## IMC 
library(imcRtools)
library(cytomapper)
library(stringr) 
library(dittoSeq)
library(CATALYST)
library(tidyverse)
library(ggrepel)
library(EBImage)
library(RColorBrewer)
library(cowplot)
library(scater)
library(tidyverse)
library(ggrepel)
library(dittoSeq)
library(viridis)
library(reshape2)
library(ggpubr)
source("/data/IMC/fun.R")

### load data
spe <- readRDS('spe.rds')
masks<- readRDS('masks.rds')

### normalize
combined_data_raw <- cbind(t(as.data.frame(counts(spe))),as.character(spe$sample_id))
combined_data_raw <- as.data.frame(combined_data_raw)
colnames(combined_data_raw)[ncol(combined_data_raw)] <- "patient_id"
combined_data_raw[,-ncol(combined_data_raw)] <- apply(combined_data_raw[,-ncol(combined_data_raw)],2,as.numeric)
### log2 data
assay(spe,'data')<- t(log2(combined_data_raw[colnames(spe),1:37]+1))

### cytofAsinh
for(i in unique(combined_data_raw$patient_id)){
  data <- combined_data_raw[combined_data_raw$patient_id==i,]
  data[,c(1:(ncol(data)-1))] <- apply(data[,c(1:(ncol(data)-1))],2,cytofAsinh)
  max_value=as.numeric(quantile(unlist(data[,markernames]),0.998))
  min_value=as.numeric(quantile(unlist(data[,markernames]),0.003))
  for (j in 1:length(markernames)){
    data[,markernames[j]]=(data[,markernames[j]]-min_value)/(max_value-min_value)
  }
  combined_data_raw[combined_data_raw$patient_id==i,] <- data
}

# combined_data_raw[,c(1:(ncol(combined_data_raw)-1))] <- apply(combined_data_raw[,c(1:(ncol(combined_data_raw)-1))],2,cytofAsinh)
# max_value=as.numeric(quantile(unlist(combined_data_raw[,markernames]),0.998))
# min_value=as.numeric(quantile(unlist(combined_data_raw[,markernames]),0.002))
# for (j in 1:length(markernames)){
#   combined_data_raw[,markernames[j]]=(combined_data_raw[,markernames[j]]-min_value)/(max_value-min_value)
# }

###############################
combined_data_raw[,markernames][combined_data_raw[,markernames]>1] = 1
combined_data_raw[,markernames][combined_data_raw[,markernames]<0] = 0
combined_data <- t(combined_data_raw[,-ncol(combined_data_raw)])
rm(combined_data_raw)
### save Exprs
assay(spe, "Exprs") <- combined_data


### Fig3A
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
set.seed(22)

type_markers <- c("Vimentin",'CD45',
                  "CD3","CD4","CD8","Foxp3",
                  "CD33","CD11b","CD68","CD163","CD11c","HLA-DR","GZMB",
                  "Pan-keratin","a-SMA")


state_markers <- c("LAG-3","CD24","CD28","CTLA-4","CD47","CD80","CD86","CD206",
                   "MERTK","PD-1","PD-L1","Siglec-10","SirPa","TIM3","THBS1")

###
mycolor <- c("#843900", "#FB8072", 
             "#7BAFDE","#1e90ff",
             "#BF0A3D","#EF8ECC","#F4800C","#33A02C",
             "#882E72",
             "#B2DF8A", "#55A1B1", "#8DD3C7",
             "#A6761D","#d9d6c3","#999999")


#### Heatmap body color 
col_exprs <- colorRamp2(c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8), 
                        viridis(option = "B",n=9))


p <- plotExprComplexHeatmap(object = spe,
                            type_markers = type_markers,
                            state_markers = state_markers,
                            het_col =col_exprs,
                            cluster_rows = T,
                            ann_bar_form ='celltype',
                            ann_bar_col =  mycolor)
dev.off()
pdf('Fig3A.pdf',he=8,wi=9)
print(p)
dev.off()

### Fig3B
somemarker_count <- as.data.frame(t(assay(spe, "Exprs")))
somemarker_count$celltype <- colData(spe)[rownames(somemarker_count),]$celltype
somemarker_count$sample_id <- colData(spe)[rownames(somemarker_count),]$sample_id
somemarker_count$Group <- colData(spe)[rownames(somemarker_count),]$Group
somemarker_count <- somemarker_count[somemarker_count$celltype%in%levels(somemarker_count$celltype)[9:14],]
somemarker_count$Celltype <- ifelse(somemarker_count$`p-HSP27`>=0.3,'p-HSP27+','p-HSP27-')

celltype_list <- list(Celltype=setNames(c('#F39B7FFF','#8491B4FF'),c('p-HSP27+','p-HSP27-')),celltype=setNames(mycolor[9:14],levels(spe$celltype)[9:14]))

cur_spe <- spe[,spe$celltype%in%levels(spe$celltype)[9:14]]
cur_spe$Celltype <- somemarker_count[colnames(cur_spe),]$Celltype

pdf(paste0('Fig3B-1.pdf'),width = 12, height = 12)
plotCells( mask = masks,
           exprs_values = 'Exprs',
           object = cur_spe, 
           cell_id = "ObjectNumber", 
           img_id = "Sample",
           colour_by = 'Celltype',
           outline_by = "celltype",
           colour = celltype_list,
           margin = 25,
           thick = TRUE)
dev.off()

#### 
somemarker_count <- as.data.frame(t(assay(spe, "Exprs")))
somemarker_count$celltype <- colData(spe)[rownames(somemarker_count),]$celltype
somemarker_count$sample_id <- colData(spe)[rownames(somemarker_count),]$sample_id
somemarker_count$Group <- colData(spe)[rownames(somemarker_count),]$Group
somemarker_count <- somemarker_count[somemarker_count$celltype%in%levels(somemarker_count$celltype)[9:14],]
somemarker_count$Celltype <- ifelse(somemarker_count$`p-Stat1`>=0.3,'p-HSP27+','p-HSP27-')

celltype_list <- list(Celltype=setNames(c("#BF0A3D","#55A1B1"),c('p-Stat1+','p-Stat1-')),celltype=setNames(mycolor[9:14],levels(spe$celltype)[9:14]))

cur_spe$Celltype <- somemarker_count[colnames(cur_spe),]$Celltype

pdf(paste0('Fig3B-2.pdf'),width = 12, height = 12)
plotCells( mask = masks,
           exprs_values = 'Exprs',
           object = cur_spe, 
           cell_id = "ObjectNumber", 
           img_id = "Sample",
           colour_by = 'Celltype',
           outline_by = "celltype",
           colour = celltype_list,
           margin = 25,
           thick = TRUE)
dev.off()

####
somemarker_count <- as.data.frame(t(assay(spe, "Exprs")))
somemarker_count$celltype <- colData(spe)[rownames(somemarker_count),]$celltype
somemarker_count$sample_id <- colData(spe)[rownames(somemarker_count),]$sample_id
somemarker_count$Group <- colData(spe)[rownames(somemarker_count),]$Group
somemarker_count <- somemarker_count[somemarker_count$celltype%in%levels(somemarker_count$celltype)[9:14],]
somemarker_count$Celltype <- ifelse(somemarker_count$`p-TRIM28`>=0.3,'p-TRIM28+','p-TRIM28-')

celltype_list <- list(Celltype=setNames(c("#843900","#33A02C"),c('p-TRIM28+','p-TRIM28-')),celltype=setNames(mycolor[9:14],levels(spe$celltype)[9:14]))

cur_spe$Celltype <- somemarker_count[colnames(cur_spe),]$Celltype

pdf(paste0('Fig3B-3.pdf'),width = 12, height = 12)
plotCells( mask = masks,
           exprs_values = 'Exprs',
           object = cur_spe, 
           cell_id = "ObjectNumber", 
           img_id = "Sample",
           colour_by = 'Celltype',
           outline_by = "celltype",
           colour = celltype_list,
           margin = 25,
           thick = TRUE)
dev.off()

### Fig3D
somemarker_count <- as.data.frame(t(assay(spe, "Exprs")))
somemarker_count$celltype <- colData(spe)[rownames(somemarker_count),]$celltype
somemarker_count$sample_id <- colData(spe)[rownames(somemarker_count),]$sample_id
somemarker_count$Group <- colData(spe)[rownames(somemarker_count),]$Group
somemarker_count <- somemarker_count[somemarker_count$celltype%in%levels(somemarker_count$celltype)[9:14],]

dat1 <- somemarker_count %>% dplyr::group_by(sample_id,Group) %>%
  dplyr::summarize(
    cellcount=n()) %>%
  select_at(vars(one_of('cellcount','sample_id','Group')))%>% unique() %>%
  data.frame()

somemarker_count <- merge(somemarker_count,dat1[,1:2],by='sample_id')

dat1 <- somemarker_count %>% dplyr::group_by(sample_id,celltype,Group,cellcount) %>%
  dplyr::summarize(pTRIM28_Freq=sum(`p-TRIM28`>0.3)/cellcount*100,
                   pStat1_Freq=sum(`p-Stat1`>0.3)/cellcount*100,
                   pHSP27_Freq=sum(`p-HSP27`>0.3)/cellcount*100,) %>%
  select_at(vars(one_of('cellcount','sample_id','Group','pTRIM28_Freq',
                        'pStat1_Freq','pHSP27_Freq','celltype')))%>% unique() %>%
  data.frame()

p_value <- compare_means(pTRIM28_Freq ~ Group,data=dat1,group.by = "celltype")

p<- ggerrorplot(dat1, x = "celltype", y = "pTRIM28_Freq",
                color = "Group", palette = "Paired",
                error.plot = "pointrange",
                add = 'jitter',
                add.params = list(alpha=0.3,size=1,width = 0.5),
                desc_stat = "mean_sd",
                position = position_dodge(0.8)) + ylab('Positive ratio')+
  ggtitle('pTRIM28')+ 
  stat_compare_means(aes(group = Group),
                     method  = "wilcox.test",
                     size = 5,
                     label ="p.signif",
                     vjust=0.5)+
  coord_flip() +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.ticks = element_line(),
        plot.title = element_text(hjust = 0.5,size = 16),
        axis.text = element_text(color = "black",size = 14))  +
  ggsci::scale_color_npg()+xlab('')

ggsave(filename = 'Fig3D/pTRIM28.celltype.pdf',plot = p,width = 5,height = 4)


p_value <- compare_means(pStat1_Freq ~ Group,data=dat1,group.by = "celltype")

p<- ggerrorplot(dat1, x = "celltype", y = "pStat1_Freq",
                color = "Group", palette = "Paired",
                error.plot = "pointrange",
                add = 'jitter',
                add.params = list(alpha=0.3,size=1,width = 0.5),
                desc_stat = "mean_sd",
                position = position_dodge(0.8)) + ylab('Positive ratio')+
  ggtitle('pSTAT1')+ 
  stat_compare_means(aes(group = Group),
                     method  = "wilcox.test",
                     size = 5,
                     label ="p.signif",
                     vjust=0.5)+
  coord_flip() +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.ticks = element_line(),
        plot.title = element_text(hjust = 0.5,size = 16),
        axis.text = element_text(color = "black",size = 14))  +
  ggsci::scale_color_npg()+xlab('')

ggsave(filename = 'Fig3D/pStat1.celltype.pdf',plot = p,width = 5,height = 4)


p_value <- compare_means(pHSP27_Freq ~ Group,data=dat1,group.by = "celltype")

p<- ggerrorplot(dat1, x = "celltype", y = "pHSP27_Freq",
                color = "Group", palette = "Paired",
                error.plot = "pointrange",
                add = 'jitter',
                add.params = list(alpha=0.3,size=1,width = 0.5),
                desc_stat = "mean_sd",
                position = position_dodge(0.8)) + ylab('Positive ratio')+
  ggtitle('pHSP27')+ 
  stat_compare_means(aes(group = Group),
                     method  = "wilcox.test",
                     size = 5,
                     label ="p.signif",
                     vjust=0.5)+
  coord_flip() +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.ticks = element_line(),
        plot.title = element_text(hjust = 0.5,size = 16),
        axis.text = element_text(color = "black",size = 14))  +
  ggsci::scale_color_npg()+xlab('')

ggsave(filename = 'Fig3D/pHSP27.celltype.pdf',plot = p,width = 5,height = 4)

### Fig3C
somemarker_count <- as.data.frame(t(assay(spe, "data")))
somemarker_count$celltype <- colData(spe)[rownames(somemarker_count),]$celltype
somemarker_count$sample_id <- colData(spe)[rownames(somemarker_count),]$sample_id
somemarker_count$Group <- colData(spe)[rownames(somemarker_count),]$Group
somemarker_count <- somemarker_count[somemarker_count$celltype%in%levels(somemarker_count$celltype)[9:14],]


dat1 <- somemarker_count %>% dplyr::group_by(sample_id,Group) %>%
  summarise_if(is.numeric,'median',na.rm=TRUE) %>%
  select_at(vars(one_of('p-HSP27','p-Stat1','p-TRIM28','sample_id','Group')))%>% unique() %>%
  data.frame()

colnames(dat1) <- gsub('\\.','-',colnames(dat1))
p1 <- ggplot(dat1,
             aes(x = Group,
                 y = `p-HSP27`,
                 color= Group)) + 
  geom_bar(stat="summary",
           fun='mean',
           fill = 'white',
           position = position_dodge(.8)) +
  geom_jitter(alpha=0.2,size=1)+
  stat_summary(color='black',fun.data = "mean_cl_normal",
               geom = "errorbar",
               width = .4) +
  stat_summary(color='black',size=1,fun = "mean", geom = "point") +
  scale_color_manual("Group",values = c("#E64B35FF","#4DBBD5FF")) + 
  ggtitle("Myeloid cells
(p-HSP27)") +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        axis.ticks = element_line(),legend.position = 'none',
        plot.title = element_text(hjust = 0.5,size = 12),
        axis.text = element_text(color = "black")) +
  labs(x = "" , y = "Expression level") +
  stat_compare_means(label ="p.signif",
                     comparisons = list(c("CT","NT")),
                     vjust=0.5) + ylim(c(0,max(dat1$`p-HSP27`)+0.2))  

p2 <- ggplot(dat1,
             aes(x = Group,
                 y = `p-Stat1`,
                 color= Group)) + 
  geom_bar(stat="summary",
           fun='mean',
           fill = 'white',
           position = position_dodge(.8)) +
  geom_jitter(alpha=0.2,size=1)+
  stat_summary(color='black',fun.data = "mean_cl_normal",
               geom = "errorbar",
               width = .4) +
  stat_summary(color='black',size=1,fun = "mean", geom = "point") +
  scale_color_manual("Group",values = c("#E64B35FF","#4DBBD5FF")) + 
  ggtitle("Myeloid cells
(p-Stat1)") +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        axis.ticks = element_line(),legend.position = 'none',
        plot.title = element_text(hjust = 0.5,size = 12),
        axis.text = element_text(color = "black")) +
  labs(x = "" , y = "Expression level") +
  stat_compare_means(label ="p.signif",
                     comparisons = list(c("CT","NT")),
                     vjust=0.5) + ylim(c(0,max(dat1$`p-Stat1`)+0.2))  

p3 <- ggplot(dat1,
             aes(x = Group,
                 y = `p-TRIM28`,
                 color= Group)) + 
  geom_bar(stat="summary",
           fun='mean',
           fill = 'white',
           position = position_dodge(.8)) +
  geom_jitter(alpha=0.2,size=1)+
  stat_summary(color='black',fun.data = "mean_cl_normal",
               geom = "errorbar",
               width = .4) +
  stat_summary(color='black',size=1,fun = "mean", geom = "point") +
  scale_color_manual("Group",values = c("#E64B35FF","#4DBBD5FF")) + 
  ggtitle("Myeloid cells
(p-TRIM28)") +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        axis.ticks = element_line(),legend.position = 'none',
        plot.title = element_text(hjust = 0.5,size = 12),
        axis.text = element_text(color = "black")) +
  labs(x = "" , y = "Expression level") +
  stat_compare_means(label ="p.signif",
                     comparisons = list(c("CT","NT")),
                     vjust=0.5) + ylim(c(0,max(dat1$`p-TRIM28`)+0.2))  


p<- grid_arrange_shared_legend_plotlist(plotlist=list(p1,p2,p3),ncol = 3,position='right')

ggsave("Fig3C.pdf", plot = p, width =5.5, height =2.5,limitsize = FALSE)

### Fig4 A
spe_Mac <- readRDS('spe_Mac.rds')

somemarkers <- c("MHCI",'MHCII','HLA-DR',"CD80","CD86",
                 "LAG-3","CD24","CD28","CTLA-4","CD47",
                 "MERTK","PD-1","PD-L1","Siglec-10","SirPa","TIM3","THBS1",
                 "p-TRIM28","p-Stat1","p-HSP27")

col_exprs <- colorRamp2(c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8),colorRampPalette(c('blue','white','red'))(12))
mycolor2 <- setNames(ggsci::pal_igv()(8),levels(spe_Mac$celltype))
p <- plotExprComplexHeatmap.V1(object = spe_Mac,
                            marker = somemarkers,
                            het_col =col_exprs,
                            cluster_rows = T,
                            ann_bar_form ='celltype',
                            ann_bar_col =  mycolor2)
dev.off()
pdf('Fig4A-1.pdf',he=8,wi=9)
print(p)
dev.off()

### Fig4 B
p<- dittoDimPlot(spe_Mac, 
                 labels.highlight = F,
                 labels.size = 4,
                 var = "cluster_id", 
                 do.label =T, 
                 reduction.use = "tSNE", 
                 size = 0.5) +
  ggtitle("Mac M2 on tSNE") + tidydr::theme_dr() + 
  scale_color_manual('Celltype',values = ggsci::pal_igv()(8)) + 
  theme(plot.title = element_text(size = 14,hjust = 0.5),
        legend.position = 'none',
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border =element_blank())

ggsave(filename = "Fig4A-2.pdf",plot = p,width = 4,height = 4)

p<- dittoDimPlot(spe_Mac, 
                 labels.highlight = F,labels.size = 4,
                 var = "Group", 
                 do.label =F, 
                 reduction.use = "tSNE", 
                 size = 0.5) +
  ggtitle("Group on tSNE") + tidydr::theme_dr() + 
  scale_color_manual('Group',values =c("#E64B35FF","#4DBBD5FF")) + 
  theme(plot.title = element_text(size = 14,hjust = 0.5),
        # legend.position = 'none',
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border =element_blank())

ggsave(filename = "Fig4B.pdf",plot = p,width = 4.5,height = 4)

### Fig4C
cluster.frequency.table <- table(spe_Mac$Sample,spe_Mac$celltype)
dat <-as.data.frame(proportions(cluster.frequency.table,1))
colnames(dat) <- c("Sample","Cluster","Freq")


sample_group=unique(as.data.frame(colData(spe_Mac))[,c("Sample","Group")])
rownames(sample_group)= sample_group[,1]
colnames(sample_group) <- c("Sample","Group")
diff.data <- merge(dat,sample_group,by="Sample")
diff.data <- as.data.frame(diff.data)
diff.data$Freq <-diff.data$Freq*100

stat.test <- diff.data %>% 
  dplyr::group_by(Cluster) %>% 
  wilcox_test(Freq~Group) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.0001,0.001,0.01,0.05,1),symbols = c('****','***','**','*','ns')) %>%
  add_xy_position(x='Cluster')


p<- ggplot(diff.data, aes(x = Cluster, y = Freq)) +
  geom_quasirandom(aes(color = Group), 
                   size = 1,width = 0.15, 
                   alpha = 0.4,
                   dodge.width =0.5,
                   method = "quasirandom") +
  
  stat_summary(aes(group = Group),
               fun.y =median,
               geom ="point",
               size = 1,
               col ="black",
               position = position_dodge(0.5)) +
  stat_summary(aes(group = Group),
               fun.y =median,
               fun.ymin =median,
               fun.ymax =median,
               geom ="crossbar",
               size = 0.3,
               width =0.5,
               col ="black",
               position = position_dodge(0.5)) +
  xlab('celltype')+
  ylab('Proportion (%)')+
  stat_pvalue_manual(stat.test,label = '{p.signif}',
                     tip.length = 0,coord.flip = F,
                     bracket.size = 0.8,vjust = 0,
                     size =4)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16)) +
  scale_color_manual(values = c("#E64B35FF","#4DBBD5FF")) 


ggsave("Fig4C-1.pdf", plot = p, width = 5, height =4,limitsize = FALSE)

celltype_list <- list(celltype=mycolor2,outline=c('outline'='black'))

pdf(paste0('Fig4C-2.pdf'),width = 12, height = 12)
plotCells( mask = masks,
           exprs_values = 'Exprs',
           object = spe_Mac, 
           cell_id = "ObjectNumber", 
           img_id = "Sample",
           colour_by = 'Celltype',
           outline_by = "outline",
           colour = celltype_list,
           margin = 25,
           thick = TRUE)
dev.off()

### Fig4D
library(mixOmics)
cluster.frequency.table <- table(spe_Mac$Sample,spe_Mac$celltype)
dat <-as.data.frame(proportions(cluster.frequency.table,1))
dat <- dcast(dat,Var1~Var2) # 长转宽
colnames(dat)[1] <- 'Sample'
sample_group=unique(as.data.frame(colData(spe))[,c("Sample","Group")])
rownames(sample_group)= sample_group[,1]
colnames(sample_group) <- c("Sample","Group")
diff.data <- merge(dat,sample_group,by="Sample")

# Compute MDS
mds <- diff.data[2:5] %>%
  dist() %>%            
  cmdscale() %>%        
  as_tibble()           
colnames(mds) <- c("Dim.1", "Dim.2")  
head(mds)

mds$Sample <- diff.data$Sample
mds$Group <- diff.data$Group


#使用ggplot2包绘图
p <-ggplot(data=mds,aes(x=Dim.1,y=Dim.2))+#指定数据、X轴、Y轴，颜色
  geom_point(aes(color = Group), shape = 19, size=3)+#绘制点图并设定大小
  #theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed", size = 1, color = 'grey50')+
  geom_hline(yintercept = 0,lty="dashed", size = 1, color = 'grey50')+#图中虚线
  theme_classic()+#主题设置
  #geom_text(aes(label=samples, y=NMDS2+0.03,x=NMDS1+0.03,vjust=0, color = samples),size=3.5, show.legend = F)+#添加数据点的标签
  stat_ellipse(data=mds,
               geom = "polygon",level=0.95,
               linetype = 2,size=0.5,
               aes(fill=Group),
               alpha=0.2)+
  scale_color_manual(values = c("#E64B35FF","#4DBBD5FF")) +
  scale_fill_manual(values = c("#E64B35FF","#4DBBD5FF"))+
  theme(axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14,angle=90),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        panel.grid=element_blank())
#geom_text_repel(aes(x=Dim.1,y=Dim.2,label=mds$Sample),size=4.5,max.overlaps = 200,force_pull = 0.5)


ggsave(filename = "Fig4D.pdf",p,width = 6,height = 5)

### Fig4E
somemarkers <- c('MHCI','MHCII','HLA-DR',"CD80","CD86",
                 "LAG-3","CD24","CD28","CTLA-4","CD47",
                 "MERTK","PD-1","PD-L1","Siglec-10","SirPa","TIM3","THBS1")

spe_T <- spe_Mac[,spe_Mac$Group%in%'CT']
somemarker_count <- as.data.frame(t(assay(spe_T, "data")))
somemarker_count$celltype <- colData(spe_T)[rownames(somemarker_count),]$celltype
somemarker_count$sample_id <- colData(spe_T)[rownames(somemarker_count),]$sample_id

dat1 <- somemarker_count %>% 
  group_by_at(c('sample_id',"celltype")) %>%
  summarise_if(is.numeric,'mean',na.rm=TRUE)%>%
  select_at(vars(one_of(somemarkers,'sample_id',"celltype")))%>%
  data.frame()

colnames(dat1)<- gsub('\\.','-',colnames(dat1))

diff.data <- melt(dat1,id.vars=c('celltype','sample_id'))
colnames(diff.data) <- c('celltype','sample_id','marker','Expr')

library(rstatix)
library(ggpubr)

df_p_val <- diff.data %>% 
  dplyr::group_by(marker) %>% 
  wilcox_test(formula = Expr ~ celltype) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.0001,0.001,0.01,0.05,1),symbols = c('****','***','**','*','ns')) %>% 
  add_xy_position(x='marker')

p<- ggplot(diff.data, aes(x = marker, y = Expr)) +
  geom_boxplot(aes(fill = celltype),position = position_dodge(width = 1),na.rm = T)+ 
  # geom_jitter(aes(fill = celltype),color='black',position = position_jitterdodge(1),size=0.5,alpha=0.3)+
  xlab('Marker')+
  ylab('Expression level')+
  stat_pvalue_manual(df_p_val,label = '{p.signif}',
                     tip.length = 0,coord.flip = T,hide.ns = T,
                     bracket.size = 0.8,vjust = 0,
                     size =4)+
  theme_classic()+ coord_flip()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16)) +
  scale_fill_manual(values = ggsci::pal_igv()(4)) + 
  geom_hline(yintercept=log2(10+1),color='red',
             linetype="dotdash",
             size=0.5) 

ggsave(filename = 'Fig4E.pdf',plot = p,width = 6,height = 9)

### Fig4F
dir.create('Fig4F')
for(i in unique(spe_Mac$celltype)){
  spe_Mac1 <- spe_Mac[,spe_Mac$celltype%in%i & spe_Mac$Group %in%'CT']
  somemarker_count <- as.data.frame(t(assay(spe_Mac1, "data")))
  somemarker_count <- somemarker_count[,somemarkers]
  somemarker_count$sample_id <- colData(spe)[rownames(somemarker_count),]$sample_id
  
  data1 <-somemarker_count%>% 
    group_by_at(c('sample_id')) %>%
    summarise_if(is.numeric,'mean',na.rm=TRUE)%>%
    select_at(vars(one_of(somemarkers,'sample_id')))%>%
    data.frame()
  
  colnames(data1) <- colnames(somemarker_count)
  data2 <- data1[,1:20]
  corr=round(cor(data2),3)
  p.mat=round(ggcorrplot::cor_pmat(data2),3)
  ggcorrplot(corr,type = 'lower')
  
  p<- ggcorrplot_fix(corr[,18:20], 
                     method = "square", 
                     type = "full", 
                     ggtheme = ggplot2::theme_void, 
                     title = "", 
                     show.legend = TRUE, 
                     legend.title = "Corr", 
                     show.diag = T, 
                     colors = c("#839EDB", "white", "#FF8D8D"), 
                     outline.color = "white", 
                     lab = F, 
                     lab_col = "black", 
                     lab_size = 3, 
                     p.mat = p.mat[,18:20], 
                     sig.level = 0.05, 
                     insig = c("pch"), 
                     pch = 8, 
                     pch.col = "black", 
                     pch.cex = 2, 
                     tl.cex = 12, 
                     tl.col = "black", 
                     tl.srt = 45, 
                     digits = 2) 
  
  ggsave(filename = paste0("Fig4F/marker_Expression.cor.",i,'.pdf'),plot = p,width = 6,height = 3)
}

### Fig5A
spe$Celltypes <- as.character(spe$celltype)
colData(spe)[colnames(spe_Mac),]$Celltypes <- as.character(spe_Mac$Celltypes)
library(ggcorrplot)
library(ggplot2)
library(ggthemes)

spe_tumour <- spe[,spe$Group %in%'CT']

###########
cluster.frequency.table <- table(spe_tumour$sample_id,spe_tumour$Celltypes)
dat <-as.data.frame(proportions(cluster.frequency.table,1))
dat <- dcast(dat,Var1~Var2)[,-1] # 长转宽

##相关性correlation
corr=round(cor(dat),3)
p.mat=round(ggcorrplot::cor_pmat(dat),3)

ggcorrplot(corr,type = 'lower')

p<- ggcorrplot_fix(corr, 
                   method = "square", 
                   type = "upper", 
                   ggtheme = ggplot2::theme_void, 
                   title = "", 
                   show.legend = TRUE, 
                   legend.title = "Corr", 
                   show.diag = T, 
                   colors = c("#839EDB", "white", "#FF8D8D"), 
                   outline.color = "white", 
                   lab = F, 
                   lab_col = "black", 
                   lab_size = 3, 
                   p.mat = p.mat, 
                   sig.level = 0.05, 
                   insig = c("pch"), 
                   pch = 8, 
                   pch.col = "black", 
                   pch.cex = 2, 
                   tl.cex = 12, 
                   tl.col = "black", 
                   tl.srt = 45, 
                   digits = 2) 

ggsave(filename = "Fig5A.pdf",plot = p,width = 6,height = 5)

### Fig5B
spe <- buildSpatialGraph(spe, img_id = "sample_id",
                         type = "knn",max_dist = 25,
                         k = 20)

####### Neighborhood aggregation
spe <- aggregateNeighbors(spe,
                          colPairName = "knn_interaction_graph",
                          aggregate_by = "metadata",
                          count_by = "Celltypes")

head(spe$aggregatedNeighbors)

##########
library(scales)
library(bluster)
library(BiocParallel)
library(ggplot2)
library('scico')

out1 <- testInteractions(spe, 
                         group_by = "sample_id",
                         label = "Celltypes", 
                         colPairName = "knn_interaction_graph",
                         method = "histocat",
                         BPPARAM = SerialParam(RNGseed = 221029))

out1$interaction <- ifelse(out1$sig & out1$sigval==1,TRUE,FALSE)
out1$sigval[which(is.na(out1$sigval))] <- 0

p <- Interaction_hetamap(object=spe,
                         het_col=scico(50, palette = 'vik'),
                         cell_order='Celltypes',
                         Interactions_file=out1,
                         group.var=c("sample_id","Group"))

pdf("Fig5B.pdf",width = 8,height = 8)
print(p)
dev.off()

### Fig5C
cluster.frequency.table <- table(spe_tumour$Sample,spe_tumour$Celltypes)
dat <-as.data.frame(proportions(cluster.frequency.table,1))
dat <- dcast(dat,Var1~Var2)[,-1] 
rownames(dat) <- rownames(cluster.frequency.table)

dat_plot <- dat
groupinfo <- as.data.frame(unique(colData(spe_tumour)[,c("Sample","Group")]))
rownames(groupinfo) <- groupinfo$Sample
dat_plot$Group <- groupinfo[rownames(dat_plot),]$Group

dir.create("Fig5C")

for(j in 1:(ncol(dat_plot)-1)){
  dir.create(paste0("Fig5C-1/",colnames(dat_plot)[j]))
  for(i in 1:(ncol(dat_plot)-1)){
    dat_plot1 <- dat_plot[,c(j,i)]
    colnames(dat_plot1) <- c('A',"B")
    p <- ggplot(data = dat_plot1, 
                mapping = aes(x = A, y = B)) +
      stat_density_2d(
        geom = "raster",
        aes(fill = after_stat(density)),
        contour = FALSE
      ) + scale_fill_gradientn(colours = colorRampPalette(c("#426ab3","#FFFFB3","#f15a22"))(256)) + 
      theme_classic()+
      geom_point() +
      geom_smooth(method = "lm", color = "black", size = 1,se = F) +
      stat_cor(method='pearson',
               label.x = max(dat_plot[,j])*0.18, label.y = max(dat_plot[,i])*0.9, #相关系数和P值位置，不设置也行
               label.sep = "\n",
      ) + 
      xlab(paste0(colnames(dat_plot)[j]," (%)")) + 
      ylab(paste0(colnames(dat_plot)[i]," (%)"))+
      theme(panel.grid = element_blank(),
            axis.title = element_text(size = 16,face = 'bold'),
            axis.text = element_text(size = 14))+ylim(0,NA)
    ggsave(filename = paste0("Fig5C-1/",colnames(dat_plot)[j],'/',colnames(dat_plot)[j],'_',colnames(dat_plot)[i],"_ggscatter.pdf"),plot = p,width = 4.5,height = 3.5)
  }
}

dir.create('Fig5C-2/M2_C1_Stromal')
spe$Celltype <- as.character(spe$Celltypes)
colData(spe)[!spe$Celltypes%in%c('M2-C1','Stromal'),]$Celltype <- 'Other' 
spe$Celltype <- factor(spe$Celltype,c('M2-C1','Stromal','Other'))

celltype_list <- list(Celltype=setNames(c("#5050FFFF","#FB8072",'#d3d7d4'),c('M2-C1','Stromal','Other')),out_line=c('outline'='black'))

### 
for(i in unique(spe$patient_id)) {
  plot.masks <- names(masks)[mcols(masks)$patient_id%in% i]
  ### 
  cur_spe <- spe[,colData(spe)$patient_id %in% c(i)]
  png(paste0('Fig5C-2/M2_C1_Stromal/',i,'.png'),width = 2500, height = 2500, res=300)
  plotCells(masks[plot.masks], object = cur_spe,
            img_id = "Sample", 
            outline_by='out_line',
            cell_id = "ObjectNumber",
            colour_by = "Celltype",
            missing_colour = "grey",
            exprs_values = "exprs",
            colour = celltype_list,
            margin = 25)
  dev.off()
  pdf(paste0('Fig5C-2/M2_C1_Stromal/',i,'.pdf'),width = 6, height = 6)
  plotCells(masks[plot.masks], object = cur_spe,
            img_id = "Sample", cell_id = "ObjectNumber",
            colour_by = "Celltype",
            outline_by='out_line',
            exprs_values = "exprs",
            missing_colour = "grey",
            colour = celltype_list,
            margin = 25)
  dev.off()
}


dir.create('Fig5C-2/M2_C2_CD8CTL')
spe$Celltype <- as.character(spe$Celltypes)
colData(spe)[!spe$Celltypes%in%c('M2-C2','CD8 CTL'),]$Celltype <- 'Other' 
spe$Celltype <- factor(spe$Celltype,c('M2-C2','CD8 CTL','Other'))

celltype_list <- list(Celltype=setNames(c("#CE3D32FF","#1e90ff",'#d3d7d4'),c('M2-C2','CD8 CTL','Other')),out_line=c('outline'='black'))

### 
for(i in unique(spe$patient_id)) {
  plot.masks <- names(masks)[mcols(masks)$patient_id%in% i]
  ### 
  cur_spe <- spe[,colData(spe)$patient_id %in% c(i)]
  png(paste0('Fig5C-2/M2_C2_CD8CTL/',i,'.png'),width = 2500, height = 2500, res=300)
  plotCells(masks[plot.masks], object = cur_spe,
            img_id = "Sample", 
            outline_by='out_line',
            cell_id = "ObjectNumber",
            colour_by = "Celltype",
            missing_colour = "grey",
            exprs_values = "exprs",
            colour = celltype_list,
            margin = 25)
  dev.off()
  pdf(paste0('Fig5C-2/M2_C2_CD8CTL/',i,'.pdf'),width = 6, height = 6)
  plotCells(masks[plot.masks], object = cur_spe,
            img_id = "Sample", cell_id = "ObjectNumber",
            colour_by = "Celltype",
            outline_by='out_line',
            exprs_values = "exprs",
            missing_colour = "grey",
            colour = celltype_list,
            margin = 25)
  dev.off()
}

### Fig5D
library(cytomapper)
library(imcRtools)
library(ggplot2)
library(ggraph)
dir.create("Fig5/D")
dir.create("Fig5/E")
#####
type <- somemarkers
type <- gsub('-','_',type)
type <- c(type ,'GZMB')

distance_data <- cbind(colData(spe)[,c("sample_id","Celltypes")],spatialCoords(spe))
distance_data1 <- CellDistance_min1cell(distance_data,celltypes="M2-C1",vars='Celltypes')
spe$distance_to_M2_C1 <- distance_data1[colnames(spe),]$distance

spe_T <- spe[,spe$Celltypes%in%c('CD4 T','CD4 CTL') & spe$Group%in%'CT']

dat2 <- data.frame()
for(i in 1:50) {
  somemarker_count <- as.data.frame(t(assay(spe_T, "data")))
  somemarker_count$distance_to_M2_C1 <- colData(spe)[rownames(somemarker_count),]$distance_to_M2_C1
  somemarker_count$sample_id <- colData(spe)[rownames(somemarker_count),]$sample_id
  somemarker_count <- somemarker_count[somemarker_count$distance_to_M2_C1>=(i-1)*2 & 
                                         somemarker_count$distance_to_M2_C1<i*2,]
  dat1 <- somemarker_count  %>%
    summarise_if(is.numeric,'median',na.rm=TRUE)%>%
    select_at(vars(one_of(type)))%>%
    data.frame()
  
  dat1$distance_to_M2_C1 <- i*2
  dat2 <- rbind(dat2,dat1)
}
dat2_M1 <- na.omit(dat2)
dat2_M1$Celltype <- 'M2-C1'
###########################
distance_data <- cbind(colData(spe)[,c("sample_id","Celltypes")],spatialCoords(spe))
distance_data1 <- CellDistance_min1cell(distance_data,celltypes="M2-C2",vars='Celltypes')
spe$distance_to_M2_C2 <- distance_data1[colnames(spe),]$distance

spe_T <- spe[,spe$Celltypes%in%c('CD4 T','CD4 CTL') & spe$Group%in%'CT']

dat2 <- data.frame()
for(i in 1:50) {
  somemarker_count <- as.data.frame(t(assay(spe_T, "data")))
  somemarker_count$distance_to_M2_C2 <- colData(spe)[rownames(somemarker_count),]$distance_to_M2_C2
  somemarker_count$sample_id <- colData(spe)[rownames(somemarker_count),]$sample_id
  somemarker_count <- somemarker_count[somemarker_count$distance_to_M2_C2>=(i-1)*2 & 
                                         somemarker_count$distance_to_M2_C2<i*2,]
  dat1 <- somemarker_count  %>%
    summarise_if(is.numeric,'median',na.rm=TRUE)%>%
    select_at(vars(one_of(type)))%>%
    data.frame()
  
  dat1$distance_to_M2_C2 <- i*2
  dat2 <- rbind(dat2,dat1)
}

dat2_M2 <- na.omit(dat2)
dat2_M2$Celltype <- 'M2-C2'
####################
distance_data <- cbind(colData(spe)[,c("sample_id","Celltypes")],spatialCoords(spe))
distance_data1 <- CellDistance_min1cell(distance_data,celltypes="M2-C3",vars='Celltypes')
spe$distance_to_M2_C3 <- distance_data1[colnames(spe),]$distance

spe_T <- spe[,spe$Celltypes%in%c('CD4 T','CD4 CTL') & spe$Group%in%'CT']

dat2 <- data.frame()
for(i in 1:50) {
  somemarker_count <- as.data.frame(t(assay(spe_T, "data")))
  somemarker_count$distance_to_M2_C3 <- colData(spe)[rownames(somemarker_count),]$distance_to_M2_C3
  somemarker_count$sample_id <- colData(spe)[rownames(somemarker_count),]$sample_id
  somemarker_count <- somemarker_count[somemarker_count$distance_to_M2_C3>=(i-1)*2 & 
                                         somemarker_count$distance_to_M2_C3<i*2,]
  dat1 <- somemarker_count  %>%
    summarise_if(is.numeric,'median',na.rm=TRUE)%>%
    select_at(vars(one_of(type)))%>%
    data.frame()
  
  dat1$distance_to_M2_C3 <- i*2
  dat2 <- rbind(dat2,dat1)
}

dat2_M3 <- na.omit(dat2)
dat2_M3$Celltype <- 'M2-C3'

dat1 <- table(spe$sample_id,spe$Celltypes)
spe_sub <- spe[,spe$sample_id%in%rownames(dat1)[dat1[,15]>0]]

distance_data <- cbind(colData(spe_sub)[,c("sample_id","Celltypes")],spatialCoords(spe_sub))
distance_data1 <- CellDistance_min1cell(distance_data,celltypes="M2-C4",vars='Celltypes')
spe_sub$distance_to_M2_C4 <- distance_data1[colnames(spe_sub),]$distance

spe_T <- spe_sub[,spe_sub$Celltypes%in%c('CD4 T','CD4 CTL') & spe_sub$Group%in%'CT']

dat2 <- data.frame()
for(i in 1:50) {
  somemarker_count <- as.data.frame(t(assay(spe_T, "data")))
  somemarker_count$distance_to_M2_C4 <- colData(spe_sub)[rownames(somemarker_count),]$distance_to_M2_C4
  somemarker_count$sample_id <- colData(spe_sub)[rownames(somemarker_count),]$sample_id
  somemarker_count <- somemarker_count[somemarker_count$distance_to_M2_C4>=(i-1)*2 & 
                                         somemarker_count$distance_to_M2_C4<i*2,]
  dat1 <- somemarker_count  %>%
    summarise_if(is.numeric,'median',na.rm=TRUE)%>%
    select_at(vars(one_of(type)))%>%
    data.frame()
  
  dat1$distance_to_M2_C4 <- i*2
  dat2 <- rbind(dat2,dat1)
}

dat2_M4 <- na.omit(dat2)
dat2_M4$Celltype <- 'M2-C4'
colnames(dat2_M1) <- gsub('distance_to_M2_C1','distance_to_M2',colnames(dat2_M1))
colnames(dat2_M2) <- gsub('distance_to_M2_C2','distance_to_M2',colnames(dat2_M2))
colnames(dat2_M3) <- gsub('distance_to_M2_C3','distance_to_M2',colnames(dat2_M3))
colnames(dat2_M4) <- gsub('distance_to_M2_C4','distance_to_M2',colnames(dat2_M4))

dat <- rbind(dat2_M1,dat2_M2,dat2_M3,dat2_M4)

p <- ggplot(data = dat, aes(x = distance_to_M2, y = PD_1,  color = Celltype)) + 
  geom_point(aes(shape=Celltype)) + geom_smooth(method = lm,se = F) +
  scale_color_manual(values = ggsci::pal_igv()(4))  +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  stat_cor(method = 'pearson', aes(x = distance_to_M2, y = PD_1,color = Celltype))

ggsave(filename = 'Fig5E/distance_to_M2_CD4T.cor_PD1.pdf',plot = p,width = 4.5,height = 3.5)

p <- ggplot(data = dat, aes(x = distance_to_M2, y = LAG_3,  color = Celltype)) + 
  geom_point(aes(shape=Celltype)) + geom_smooth(method = lm,se = F) +
  scale_color_manual(values = ggsci::pal_igv()(4))  +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  stat_cor(method = 'pearson', aes(x = distance_to_M2, y = LAG_3,color = Celltype))

ggsave(filename = 'Fig5E/distance_to_M2_CD4T.cor_LAG3.pdf',plot = p,width = 4.5,height = 3.5)

p <- ggplot(data = dat, aes(x = distance_to_M2, y = TIM3,  color = Celltype)) + 
  geom_point(aes(shape=Celltype)) + geom_smooth(method = lm,se = F) +
  scale_color_manual(values = ggsci::pal_igv()(4))  +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  stat_cor(method = 'pearson', aes(x = distance_to_M2, y = TIM3,color = Celltype))

ggsave(filename = 'Fig5E/distance_to_M2_CD4T.cor_TIM3.pdf',plot = p,width = 4.5,height = 3.5)

p <- ggplot(data = dat, aes(x = distance_to_M2, y = CTLA_4,  color = Celltype)) + 
  geom_point(aes(shape=Celltype)) + geom_smooth(method = lm,se = F) +
  scale_color_manual(values = ggsci::pal_igv()(4))  +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  stat_cor(method = 'pearson', aes(x = distance_to_M2, y = CTLA_4,color = Celltype))

ggsave(filename = 'Fig5E/distance_to_M2_CD4T.cor_CTLA_4.pdf',plot = p,width = 4.5,height = 3.5)

p <- ggplot(data = dat, aes(x = distance_to_M2, y = CD28,  color = Celltype)) + 
  geom_point(aes(shape=Celltype)) + geom_smooth(method = lm,se = F) +
  scale_color_manual(values = ggsci::pal_igv()(4))  +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  stat_cor(method = 'pearson', aes(x = distance_to_M2, y = CD28,color = Celltype),label.y.npc = 0.25)

ggsave(filename = 'Fig5E/distance_to_M2_CD4T.cor_CD28.pdf',plot = p,width = 4.5,height = 3.5)

p <- ggplot(data = dat, aes(x = distance_to_M2, y = GZMB,  color = Celltype)) + 
  geom_point(aes(shape=Celltype)) + geom_smooth(method = lm,se = F) +
  scale_color_manual(values = ggsci::pal_igv()(4))  +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  stat_cor(method = 'pearson', aes(x = distance_to_M2, y = GZMB,color = Celltype),label.y.npc = 0.25)

ggsave(filename = 'Fig5E/distance_to_M2_CD4T.cor_GZMB.pdf',plot = p,width = 4.5,height = 3.5)

################

distance_data <- cbind(colData(spe)[,c("sample_id","Celltypes")],spatialCoords(spe))
distance_data1 <- CellDistance_min1cell(distance_data,celltypes="M2-C1",vars='Celltypes')
spe$distance_to_M2_C1 <- distance_data1[colnames(spe),]$distance

spe_T <- spe[,spe$Celltypes%in%c('CD8 T','CD8 CTL') & spe$Group%in%'CT']

dat2 <- data.frame()
for(i in 1:50) {
  somemarker_count <- as.data.frame(t(assay(spe_T, "data")))
  somemarker_count$distance_to_M2_C1 <- colData(spe)[rownames(somemarker_count),]$distance_to_M2_C1
  somemarker_count$sample_id <- colData(spe)[rownames(somemarker_count),]$sample_id
  somemarker_count <- somemarker_count[somemarker_count$distance_to_M2_C1>=(i-1)*2 & 
                                         somemarker_count$distance_to_M2_C1<i*2,]
  dat1 <- somemarker_count  %>%
    summarise_if(is.numeric,'median',na.rm=TRUE)%>%
    select_at(vars(one_of(type)))%>%
    data.frame()
  
  dat1$distance_to_M2_C1 <- i*2
  dat2 <- rbind(dat2,dat1)
}
dat2_M1 <- na.omit(dat2)
dat2_M1$Celltype <- 'M2-C1'
###########################
distance_data <- cbind(colData(spe)[,c("sample_id","Celltypes")],spatialCoords(spe))
distance_data1 <- CellDistance_min1cell(distance_data,celltypes="M2-C2",vars='Celltypes')
spe$distance_to_M2_C2 <- distance_data1[colnames(spe),]$distance

spe_T <- spe[,spe$Celltypes%in%c('CD8 T','CD8 CTL') & spe$Group%in%'CT']

dat2 <- data.frame()
for(i in 1:50) {
  somemarker_count <- as.data.frame(t(assay(spe_T, "data")))
  somemarker_count$distance_to_M2_C2 <- colData(spe)[rownames(somemarker_count),]$distance_to_M2_C2
  somemarker_count$sample_id <- colData(spe)[rownames(somemarker_count),]$sample_id
  somemarker_count <- somemarker_count[somemarker_count$distance_to_M2_C2>=(i-1)*2 & 
                                         somemarker_count$distance_to_M2_C2<i*2,]
  dat1 <- somemarker_count  %>%
    summarise_if(is.numeric,'median',na.rm=TRUE)%>%
    select_at(vars(one_of(type)))%>%
    data.frame()
  
  dat1$distance_to_M2_C2 <- i*2
  dat2 <- rbind(dat2,dat1)
}

dat2_M2 <- na.omit(dat2)
dat2_M2$Celltype <- 'M2-C2'
####################
distance_data <- cbind(colData(spe)[,c("sample_id","Celltypes")],spatialCoords(spe))
distance_data1 <- CellDistance_min1cell(distance_data,celltypes="M2-C3",vars='Celltypes')
spe$distance_to_M2_C3 <- distance_data1[colnames(spe),]$distance

spe_T <- spe[,spe$Celltypes%in%c('CD8 T','CD8 CTL') & spe$Group%in%'CT']

dat2 <- data.frame()
for(i in 1:50) {
  somemarker_count <- as.data.frame(t(assay(spe_T, "data")))
  somemarker_count$distance_to_M2_C3 <- colData(spe)[rownames(somemarker_count),]$distance_to_M2_C3
  somemarker_count$sample_id <- colData(spe)[rownames(somemarker_count),]$sample_id
  somemarker_count <- somemarker_count[somemarker_count$distance_to_M2_C3>=(i-1)*2 & 
                                         somemarker_count$distance_to_M2_C3<i*2,]
  dat1 <- somemarker_count  %>%
    summarise_if(is.numeric,'median',na.rm=TRUE)%>%
    select_at(vars(one_of(type)))%>%
    data.frame()
  
  dat1$distance_to_M2_C3 <- i*2
  dat2 <- rbind(dat2,dat1)
}

dat2_M3 <- na.omit(dat2)
dat2_M3$Celltype <- 'M2-C3'

dat1 <- table(spe$sample_id,spe$Celltypes)
spe_sub <- spe[,spe$sample_id%in%rownames(dat1)[dat1[,15]>0]]

distance_data <- cbind(colData(spe_sub)[,c("sample_id","Celltypes")],spatialCoords(spe_sub))
distance_data1 <- CellDistance_min1cell(distance_data,celltypes="M2-C4",vars='Celltypes')
spe_sub$distance_to_M2_C4 <- distance_data1[colnames(spe_sub),]$distance

spe_T <- spe_sub[,spe_sub$Celltypes%in%c('CD8 T','CD8 CTL') & spe_sub$Group%in%'CT']

dat2 <- data.frame()
for(i in 1:50) {
  somemarker_count <- as.data.frame(t(assay(spe_T, "data")))
  somemarker_count$distance_to_M2_C4 <- colData(spe_sub)[rownames(somemarker_count),]$distance_to_M2_C4
  somemarker_count$sample_id <- colData(spe_sub)[rownames(somemarker_count),]$sample_id
  somemarker_count <- somemarker_count[somemarker_count$distance_to_M2_C4>=(i-1)*2 & 
                                         somemarker_count$distance_to_M2_C4<i*2,]
  dat1 <- somemarker_count  %>%
    summarise_if(is.numeric,'median',na.rm=TRUE)%>%
    select_at(vars(one_of(type)))%>%
    data.frame()
  
  dat1$distance_to_M2_C4 <- i*2
  dat2 <- rbind(dat2,dat1)
}

dat2_M4 <- na.omit(dat2)
dat2_M4$Celltype <- 'M2-C4'
colnames(dat2_M1) <- gsub('distance_to_M2_C1','distance_to_M2',colnames(dat2_M1))
colnames(dat2_M2) <- gsub('distance_to_M2_C2','distance_to_M2',colnames(dat2_M2))
colnames(dat2_M3) <- gsub('distance_to_M2_C3','distance_to_M2',colnames(dat2_M3))
colnames(dat2_M4) <- gsub('distance_to_M2_C4','distance_to_M2',colnames(dat2_M4))

dat <- rbind(dat2_M1,dat2_M2,dat2_M3,dat2_M4)

p <- ggplot(data = dat, aes(x = distance_to_M2, y = PD_1,  color = Celltype)) + 
  geom_point(aes(shape=Celltype)) + geom_smooth(method = lm,se = F) +
  scale_color_manual(values = ggsci::pal_igv()(4))  +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  stat_cor(method = 'pearson', aes(x = distance_to_M2, y = PD_1,color = Celltype))

ggsave(filename = 'Fig5D/distance_to_M2_CD8T.cor_PD1.pdf',plot = p,width = 4.5,height = 3.5)

p <- ggplot(data = dat, aes(x = distance_to_M2, y = LAG_3,  color = Celltype)) + 
  geom_point(aes(shape=Celltype)) + geom_smooth(method = lm,se = F) +
  scale_color_manual(values = ggsci::pal_igv()(4))  +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  stat_cor(method = 'pearson', aes(x = distance_to_M2, y = LAG_3,color = Celltype))

ggsave(filename = 'Fig5D/distance_to_M2_CD8T.cor_LAG3.pdf',plot = p,width = 4.5,height = 3.5)

p <- ggplot(data = dat, aes(x = distance_to_M2, y = TIM3,  color = Celltype)) + 
  geom_point(aes(shape=Celltype)) + geom_smooth(method = lm,se = F) +
  scale_color_manual(values = ggsci::pal_igv()(4))  +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  stat_cor(method = 'pearson', aes(x = distance_to_M2, y = TIM3,color = Celltype))

ggsave(filename = 'Fig5D/distance_to_M2_CD8T.cor_TIM3.pdf',plot = p,width = 4.5,height = 3.5)

p <- ggplot(data = dat, aes(x = distance_to_M2, y = CTLA_4,  color = Celltype)) + 
  geom_point(aes(shape=Celltype)) + geom_smooth(method = lm,se = F) +
  scale_color_manual(values = ggsci::pal_igv()(4))  +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  stat_cor(method = 'pearson', aes(x = distance_to_M2, y = CTLA_4,color = Celltype))

ggsave(filename = 'Fig5D/distance_to_M2_CD8T.cor_CTLA_4.pdf',plot = p,width = 4.5,height = 3.5)

p <- ggplot(data = dat, aes(x = distance_to_M2, y = CD28,  color = Celltype)) + 
  geom_point(aes(shape=Celltype)) + geom_smooth(method = lm,se = F) +
  scale_color_manual(values = ggsci::pal_igv()(4))  +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  stat_cor(method = 'pearson', aes(x = distance_to_M2, y = CD28,color = Celltype),label.y.npc = 0.2)

ggsave(filename = 'Fig5D/distance_to_M2_CD8T.cor_CD28.pdf',plot = p,width = 4.5,height = 3.5)

p <- ggplot(data = dat, aes(x = distance_to_M2, y = GZMB,  color = Celltype)) + 
  geom_point(aes(shape=Celltype)) + geom_smooth(method = lm,se = F) +
  scale_color_manual(values = ggsci::pal_igv()(4))  +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5)) + 
  stat_cor(method = 'pearson', aes(x = distance_to_M2, y = GZMB,color = Celltype),label.y.npc = 0.3)

ggsave(filename = 'Fig5D/distance_to_M2_CD8T.cor_GZMB.pdf',plot = p,width = 4.5,height = 3.5)

### Fig6A
PFS <- read.csv('PFS.csv',header = T)

cluster.frequency.table <- table(spe_Mac$Sample,spe_Mac$Celltypes)
dat <-as.data.frame(proportions(cluster.frequency.table,1))
colnames(dat) <- c("Sample","Cluster","Freq")
dat <- dat[dat$Sample %in% PFS$Sample,]

dat <- dcast(dat,Sample~Cluster)

df.final <- merge(dat,PFS,by='Sample')
df.final$time <- df.final$PFS.time/30
df.final$event <- df.final$PFS

for( i in colnames(df.final)[2:5]){
  p <- my_ggsurvplot(data =df.final,Var = i)
  
  ggsave(paste0('Fig6/','ggsurvplot_',i,'.pdf'),plot =p,height = 4,width = 4)
  
}

### Fig6B
cluster.frequency.table <- table(spe_Mac$Sample,spe_Mac$celltype)
dat <-as.data.frame(proportions(cluster.frequency.table,1))
colnames(dat) <- c("Sample","Cluster","Freq")


sample_group=unique(as.data.frame(colData(spe_Mac))[,c("Sample","RECIST")])
rownames(sample_group)= sample_group[,1]
colnames(sample_group) <- c("Sample","Group")
diff.data <- merge(dat,sample_group,by="Sample")
diff.data <- as.data.frame(diff.data)
diff.data$Freq <-diff.data$Freq*100
diff.data <- diff.data[diff.data$Group%in%c("HR","LR"),]

p<- ggboxplot(diff.data, x = 'Cluster', y = 'Freq',color ="Group",
              palette =c("#F39B7FFF","#8491B4FF"),
              add = "jitter", 
              # bxp.errorbar = T,
              ylab ='Proportion(%)',xlab = '') + 
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)) + 
  stat_compare_means(aes(group = Group),size=4.5,
                     label ="p.signif",label.y = 108,
                     # method = "t.test",
                     vjust=0.5) + theme(legend.position="right") +
  ylim(c(0,115))


ggsave("Fig6B.pdf", plot = p, width = 4.5, height =3.5,limitsize = FALSE)

### Fig6C
celltype_list <- list(celltype=mycolor2,outline=c('outline'='black'))


pdf('Fig6C',width = 12, height = 12)
plotCells( mask = masks,
           exprs_values = 'Exprs',
           object = cur_spe, 
           cell_id = "ObjectNumber", 
           img_id = "Sample",
           colour_by = 'celltype',
           outline_by = "outline",
           colour = celltype_list,
           margin = 25,
           thick = TRUE)
dev.off()

### Fig6D
library(mixOmics)
spe_sub <- spe_Mac[,spe_Mac$RECIST%in%c("HR","LR")]
cluster.frequency.table <- table(spe_sub$Sample,spe_sub$celltype)
dat <-as.data.frame(proportions(cluster.frequency.table,1))
dat <- dcast(dat,Var1~Var2) 
colnames(dat)[1] <- 'Sample'
sample_group=unique(as.data.frame(colData(spe))[,c("Sample","RECIST")])
rownames(sample_group)= sample_group[,1]
colnames(sample_group) <- c("Sample","Group")
diff.data <- merge(dat,sample_group,by="Sample")

# Compute MDS
mds <- diff.data[2:5] %>%
  dist() %>%             
  cmdscale() %>%        
  as_tibble()           
colnames(mds) <- c("Dim.1", "Dim.2")  
head(mds)

mds$Sample <- diff.data$Sample
mds$Group <- diff.data$Group


p <-ggplot(data=mds,aes(x=Dim.1,y=Dim.2))+
  geom_point(aes(color = Group), shape = 19, size=3)+
  #theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed", size = 1, color = 'grey50')+
  geom_hline(yintercept = 0,lty="dashed", size = 1, color = 'grey50')+
  theme_classic()+#
  #geom_text(aes(label=samples, y=NMDS2+0.03,x=NMDS1+0.03,vjust=0, color = samples),size=3.5, show.legend = F)+#添加数据点的标签
  stat_ellipse(data=mds,#
               geom = "polygon",level=0.95,#
               linetype = 2,size=0.5,#
               aes(fill=Group),
               alpha=0.2)+
  scale_color_manual(values = c("#F39B7FFF","#8491B4FF")) +
  scale_fill_manual(values = c("#F39B7FFF","#8491B4FF"))+
  theme(axis.title.x=element_text(size=14),#
        axis.title.y=element_text(size=14,angle=90),#
        axis.text.y=element_text(size=12),#
        axis.text.x=element_text(size=12),#
        panel.grid=element_blank())+ #
  geom_text_repel(aes(x=Dim.1,y=Dim.2,label=mds$Sample),size=4.5,max.overlaps = 200,force_pull = 0.5)


ggsave(filename = "Fig6D.pdf",p,width = 6,height = 5)

### Fig7A
distance_data <- cbind(colData(spe)[,c("sample_id","Celltypes")],spatialCoords(spe))
distance_data1 <- CellDistance_min1cell(distance_data,celltypes="Epithelial",vars='Celltypes')
spe$distance_to_Epi <- distance_data1[colnames(spe),]$distance
spe$distance_to_Epi_fix <- spe$distance_to_Epi
spe$distance_to_Epi_fix[spe$distance_to_Epi_fix>50] <- 50

plot.dat <- as.data.frame(colData(spe)[,c("distance_to_Epi","Celltypes","sample_id")])
colnames(plot.dat) <- c("distance_to_Epi","celltype","sample_id")
plot.dat <- plot.dat[plot.dat$celltype!='Epithelial'& plot.dat$distance_to_Epi <= 25,]

max_bin <- ceiling(max(plot.dat$distance_to_Epi/1))
plot.dat$binID <- 0
for(i in 1:max_bin){
  if(nrow(plot.dat[((i-1)*1) < plot.dat$distance_to_Epi & plot.dat$distance_to_Epi <= 1*i,])>0) {
    plot.dat[((i-1)*1) < plot.dat$distance_to_Epi & plot.dat$distance_to_Epi <= 1*i,]$binID <- i
  } else {next}
}

plot.dat1 <- c()
for(i in levels(plot.dat$celltype)[levels(plot.dat$celltype) %in% unique(plot.dat$celltype)]){
  plot.dat3 <- as.data.frame(matrix(0,nrow = max_bin,ncol = 6))
  plot.dat2 <- plot.dat %>% group_by(sample_id,binID) %>% dplyr::summarize(cell_count = sum(celltype==i),all_count=n(),cell_freq =cell_count/all_count*100) %>% as.data.frame() 
  plot.dat2$celltype <- i
  colnames(plot.dat3) <- colnames(plot.dat2)
  plot.dat3$celltype <- i
  plot.dat3$binID <- c(1:max_bin)
  for(j in levels(plot.dat2$sample_id)) {
    plot.dat3$sample_id <- j
    plot.dat2 <- rbind(plot.dat2,plot.dat3[setdiff(plot.dat3$binID,unique(plot.dat2[plot.dat2$sample_id==j,]$binID)),])
  }
  plot.dat1 <- rbind(plot.dat2,plot.dat1)
  plot.dat1 <- arrange(plot.dat1,sample_id,celltype,binID)
}

groupinfo <- unique(as.data.frame(colData(spe)[,c("sample_id","Group",'RECIST')]))
plot.dat1 <- merge(groupinfo,plot.dat1,by="sample_id")
plot.dat1 <- arrange(plot.dat1,sample_id,celltype,binID)
# plot.dat2 <- plot.dat1 %>% group_by(Group,binID,celltype) %>% dplyr::summarize(Freq=mean(cell_freq)) %>% as.data.frame() 

plot.dat1$celltype <- factor(plot.dat1$celltype,levels = levels(plot.dat$celltype))
plot.dat1 <- plot.dat1[plot.dat1$RECIST%in%c('HR','LR') & plot.dat1$celltype %in%levels(plot.dat1$celltype)[12:15],]


p <- ggplot(data = plot.dat1,aes(x=binID,y=cell_freq,color=celltype))+
  geom_smooth(method="loess",se=F) + theme_classic2() + theme(panel.grid = element_blank()) + 
  ggtitle("The Distance of Others to Epithelial") + theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~RECIST,scales = 'free_x') + coord_cartesian(ylim = c(0,20)) + scale_x_continuous(breaks = seq(0,35,3)) +
  scale_color_manual(values =ggsci::pal_igv()(8)) +
  xlab("distance (µm)") + ylab("Cell type fraction")  

ggsave("Fig7A.pdf", plot = p, width = 8, height =5,limitsize = FALSE)

### Fig7B
celltype_list <- list(Celltypes=setNames(c('#78331e',ggsci::pal_igv()(4)),levels(spe$Celltypes)[c(1,12:15)]),
                      distance_to_Epi_fix = rev(cmap(100)))

spe_sub <- spe[,spe$Celltypes%in%levels(spe$Celltypes)[c(1,12:15)]]

pdf(paste0('Fig7/distance/Epi/map/',i,'_with_celltype.pdf'),width = 6, height = 6)
plotCells(masks[unique(spe_sub$sample_id)], 
          object = spe_sub,
          img_id = "Sample", 
          cell_id = "ObjectNumber",
          colour_by = "distance_to_Epi_fix",
          legend = list(margin = 200),
          outline_by = "Celltypes",
          margin = 25,
          colour = celltype_list)
dev.off()

### Fig7C
spe_sub <- spe[,spe$RECIST%in%c("HR","LR")]
out2 <- out1[out1$group_by%in%unique(spe_sub )]
library('scico')
p <- Interaction_hetamap(object=spe_sub,
                         het_col=scico(50, palette = 'vik'),
                         cell_order='celltype',
                         Interactions_file=out2,
                         group.var=c("sample_id","Group"))

pdf("Fig7C.pdf",width = 6,height = 6)
print(p)
dev.off()

### Fig7D
mycolor1 <- c("#843900", "#FB8072", 
              "#7BAFDE","#1e90ff",
              "#BF0A3D","#EF8ECC","#F4800C","#33A02C",
              "#882E72",
              "#B2DF8A", "#55A1B1", 
              "#8DD3C7","#8DD3C7","#8DD3C7","#8DD3C7",
              "#A6761D","#d9d6c3","#999999")

neighborhood_dat <- as.data.frame(colPair(spe, "knn_interaction_graph"))
spe$cellID <- colnames(spe) 
neighborhood_dat$from_cellID <- spe$cellID[neighborhood_dat$form]
neighborhood_dat$to_cellID <- spe$cellID[neighborhood_dat$to]
neighborhood_dat$from_celltype <- spe$Celltypes[neighborhood_dat$form]
neighborhood_dat$to_celltype <- spe$Celltypes[neighborhood_dat$to]
spe$outline <- 'outline'

for(i in levels(spe$Celltypes)[12:15]){
  neighborhood_dat_sub <- neighborhood_dat[neighborhood_dat$from_celltype%in%i,]
  spe_sub <- spe[,unique(c(neighborhood_dat_sub$from_cellID,neighborhood_dat_sub$to_cellID))]
  celltype_list <- list(Celltypes=setNames(mycolor1,levels(spe$Celltypes)),outline=c('outline'='black'))
  pdf(paste0('Fig7D/',i,'.pdf'),width = 6, height = 6)
  plotCells(masks[unique(spe_sub$sample_id)], 
            object = spe_sub,
            img_id = "Sample", 
            cell_id = "ObjectNumber",
            colour_by = "Celltypes",
            outline_by='outline',
            exprs_values = "exprs",
            missing_colour = "grey",
            colour = celltype_list,
            margin = 25)
  dev.off()
}

### Fig7E
somemarkers <- c('Foxp3','a-SMA',somemarkers)
interaction_exprs <- data.frame()
for(i in levels(spe$Celltypes)[12:15]){
  neighborhood_dat_sub <- neighborhood_dat[neighborhood_dat$from_celltype%in%i,]
  spe_sub <- spe[,unique(c(neighborhood_dat_sub$from_cellID,neighborhood_dat_sub$to_cellID))]
  spe_sub <- spe_sub[,spe_sub$Group%in%'CT']
  somemarker_count <- as.data.frame(t(assay(spe_sub, "data")))
  somemarker_count <- somemarker_count[,somemarkers]
  somemarker_count$Sample <- colData(spe)[rownames(somemarker_count),]$Sample
  somemarker_count$Celltype <- colData(spe)[rownames(somemarker_count),]$Celltypes
  data1_C1 <-somemarker_count
  data1_C1$Term <- i
  interaction_exprs <- rbind(interaction_exprs,data1_C1)
}

interaction_exprs$Term <- gsub('M2-C','CN',interaction_exprs$Term)
dat2 <- melt(interaction_exprs,id.vars=c('Celltype','Sample','Term'))
colnames(dat2) <- c('Celltype','Sample','Term','marker','expr')

diff.data1 <- interaction_exprs[interaction_exprs$Celltype %in%'Hepatocyte',]

stat.test <- diff.data1 %>% 
  dplyr::group_by(marker) %>% 
  wilcox_test(expr~Term) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.0001,0.001,0.01,0.05,1),symbols = c('****','***','**','*','ns')) %>%
  add_xy_position(x='marker')


p<- ggviolin(data = diff.data1,
             x ='marker',
             y = 'expr',
             trim = F,
             color = "Term",
             width = 0.8,
             palette = ggsci::pal_lancet()(4),
             ylab = 'Expression level',
             add = c("boxplot"), #添加图 "dotplot", "jitter", "boxplot", "point"
             add.params = list( width = 0.2,#宽度
                                linetype = 1),#线型
             short.panel.labs = F) + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)) + 
  stat_pvalue_manual(stat.test,label = '{p.signif}',
                     tip.length = 0,
                     bracket.size = 0.8,vjust = 0,
                     size =4)+
  # ylim(c(0,8))+
  theme(legend.position="right",plot.title = element_text(size = 16,hjust = 0.5))+
  geom_hline(yintercept = log2(10+1),
             linetype="longdash",
             size=0.5,
             color = "#ed1941")

ggsave(filename = 'Fig7E.pdf',
       plot = p,height = 4.5,width = 8)

### Fig7F
interaction_exprs <- data.frame()
for(i in levels(spe$Celltypes)[12:15]){
  neighborhood_dat_sub <- neighborhood_dat[neighborhood_dat$from_celltype%in%i,]
  spe_sub <- spe[,unique(c(neighborhood_dat_sub$from_cellID,neighborhood_dat_sub$to_cellID))]
  spe_sub <- spe_sub[,spe_sub$Group%in%'CT']
  somemarker_count <- as.data.frame(t(assay(spe_sub, "data")))
  somemarker_count <- somemarker_count[,somemarkers]
  somemarker_count$Sample <- colData(spe)[rownames(somemarker_count),]$Sample
  somemarker_count$Celltype <- colData(spe)[rownames(somemarker_count),]$Celltypes
  data1_C1 <-somemarker_count
  data1_C1$Term <- i
  interaction_exprs <- rbind(interaction_exprs,data1_C1)
}

interaction_exprs$Term <- gsub('M2-C','CN',interaction_exprs$Term)
dat2 <- melt(interaction_exprs,id.vars=c('Celltype','Sample','Term'))
colnames(dat2) <- c('Celltype','Sample','Term','marker','expr')

diff.data1 <- interaction_exprs[interaction_exprs$Celltype %in%c('CD8 T','CD8 CTL'),]

stat.test <- diff.data1 %>% 
  dplyr::group_by(marker) %>% 
  wilcox_test(expr~Term) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.0001,0.001,0.01,0.05,1),symbols = c('****','***','**','*','ns')) %>%
  add_xy_position(x='marker')


p<- ggviolin(data = diff.data1,
             x ='marker',
             y = 'expr',
             trim = F,
             color = "Term",
             width = 0.8,
             palette = ggsci::pal_lancet()(4),
             ylab = 'Expression level',
             add = c("boxplot"), #添加图 "dotplot", "jitter", "boxplot", "point"
             add.params = list( width = 0.2,#宽度
                                linetype = 1),#线型
             short.panel.labs = F) + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)) + 
  stat_pvalue_manual(stat.test,label = '{p.signif}',
                     tip.length = 0,
                     bracket.size = 0.8,vjust = 0,
                     size =4)+
  # ylim(c(0,8))+
  theme(legend.position="right",plot.title = element_text(size = 16,hjust = 0.5))+
  geom_hline(yintercept = log2(10+1),
             linetype="longdash",
             size=0.5,
             color = "#ed1941")

ggsave(filename = 'Fig7F-1.pdf',
       plot = p,height = 4.5,width = 8)

diff.data1 <- interaction_exprs[interaction_exprs$Celltype %in%c('CD4 T','CD4 CTL'),]

stat.test <- diff.data1 %>% 
  dplyr::group_by(marker) %>% 
  wilcox_test(expr~Term) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.0001,0.001,0.01,0.05,1),symbols = c('****','***','**','*','ns')) %>%
  add_xy_position(x='marker')


p<- ggviolin(data = diff.data1,
             x ='marker',
             y = 'expr',
             trim = F,
             color = "Term",
             width = 0.8,
             palette = ggsci::pal_lancet()(4),
             ylab = 'Expression level',
             add = c("boxplot"), #添加图 "dotplot", "jitter", "boxplot", "point"
             add.params = list( width = 0.2,#宽度
                                linetype = 1),#线型
             short.panel.labs = F) + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)) + 
  stat_pvalue_manual(stat.test,label = '{p.signif}',
                     tip.length = 0,
                     bracket.size = 0.8,vjust = 0,
                     size =4)+
  # ylim(c(0,8))+
  theme(legend.position="right",plot.title = element_text(size = 16,hjust = 0.5))+
  geom_hline(yintercept = log2(10+1),
             linetype="longdash",
             size=0.5,
             color = "#ed1941")

ggsave(filename = 'Fig7F-2.pdf',
       plot = p,height = 4.5,width = 8)
