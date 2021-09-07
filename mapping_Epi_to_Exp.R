### load R packages
suppressWarnings(library(cisTopic))
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb<-TxDb.Mmusculus.UCSC.mm9.knownGene
library(ChIPseeker)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
### get cell name
Exp_k4 <- readRDS('Exp_k4.rds')
Exp_k27 <- readRDS('Exp_k27.rds')
cell.k4<-colnames(Exp_k4)
cell.k27<-colnames(Exp_k27)

### load cisTopic object
### this object is created by standard pipline of cisTopic
Epi_k27<-readRDS('cisTopic_result/Wnt_beads_H3K27me3_scSET-seq_Exp.rds')
Epi_k4<-readRDS('cisTopic_result/Wnt_beads_H3K4me3_scSET-seq_Exp.rds')

### function to mapping Barcode name to Epi_data
changeBC<-function(x,i7,i5,i7start=1,i7end=length(i7),i5start=1,i5end=length(i5)){
  for (i in i7start:i7end) {
    for (j in i5start:i5end) {
      if (i==i7start & j==i5start) {
        x<-gsub(x,pattern = i5[j],replacement = LETTERS[j])
        x<-gsub(x,pattern = i7[i],replacement = i)
      }  else {
        x<-gsub(x,pattern = i5[j],replacement = LETTERS[j])
        x<-gsub(x,pattern = i7[i],replacement = i)
      }
    }
  }  
  return(x)
}
changeEpiBC<-function(x){
  i5<-read.table('i5.txt')
  i5<-as.matrix(i5)
  i7<-read.table('i7.txt')
  i7<-as.matrix(i7)
  cell.names <- x@cell.names
  new.cell.names<-gsub(pattern ='.bam',replacement = '',cell.names)
  new.cell.names<-gsub(pattern ='_Epi',replacement = '',cell.names)
  new.cell.names<-changeBC(new.cell.names,i7,i5)
  x <- renameCells(x, new.cell.names)
  cell.type<-x@cell.names
  i7=as.character(c(1:12))
  index1 = seq(1,11,2) 
  index2 = seq(2,12,2)
  ### use last 2 character to define cell distance to Wnt beads
  for (i in index1) {
    cell.type[which(str_detect(str_sub(cell.type,-2,-1), pattern = i7[i]))] <- "Proximal"
  }
  for (i in index2) {
    cell.type[which(str_detect(str_sub(cell.type,-2,-1), pattern = i7[i]))] <- "Distal"
  }
  cellData<-data.frame(row.names = x@cell.names,cell.type)
  
  x <- addCellMetadata(x, cell.data = cellData)
  return(x)
}

#### mapping new metadata to cisTopic object (How to create cisTopic object: https://github.com/aertslab/cisTopic)
Epi_k27<-changeEpiBC(Epi_k27)
Epi_k27@count.matrix
Epi_k4<-changeEpiBC(Epi_k4)
Epi_k4@count.matrix
### subset cells by scSET-seq_Exp results
subset_counts <- Epi_k27@count.matrix[,cell.k27]
subset_metadata <- Epi_k27@cell.data[cell.k27,]
Epi_k27 <- createcisTopicObject(subset_counts)
Epi_k27 <- addCellMetadata(Epi_k27, subset_metadata)
subset_counts <- Epi_k4@count.matrix[,cell.k4]
subset_metadata <- Epi_k4@cell.data[cell.k4,]
Epi_k4 <- createcisTopicObject(subset_counts)
Epi_k4 <- addCellMetadata(Epi_k4, subset_metadata)

### select best model for Epi data
Epi_k27<- runWarpLDAModels(Epi_k27,topic = c(2,5,10:25,30,35,40),seed=987,nCores = 2,iterations = 500,addModels = F)
par(mfrow=c(3,3))
Epi_k27<-selectModel(Epi_k27,type='maximum')
Epi_k27<-selectModel(Epi_k27,type='perplexity')
Epi_k27<-selectModel(Epi_k27,type='derivative')
Epi_k27<-runUmap(Epi_k27,target = 'cell',seed = 123,method = 'Probability')
Epi_k4<- runWarpLDAModels(Epi_k4,topic = c(2,5,10:25,30,35,40),seed=987,nCores = 2,iterations = 500,addModels = F)
par(mfrow=c(3,3))
Epi_k4<-selectModel(Epi_k4,type='maximum')
Epi_k4<-selectModel(Epi_k4,type='perplexity')
Epi_k4<-selectModel(Epi_k4,type='derivative')
Epi_k4<-runUmap(Epi_k4,target = 'cell',seed = 123,method = 'Probability')
### save RDS
saveRDS(Epi_k27,'Epi_k27.rds')
Epi_k27<-readRDS('Epi_k27.rds')
saveRDS(Epi_k4,'Epi_k4.rds')
Epi_k4<-readRDS('Epi_k4.rds')

### load marker gene info 
wnt <- readRDS('wnt.rds')
### get top50 marker genes
wnt.markers <- FindAllMarkers(wnt, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
wnt.markers <- as.data.frame(wnt.markers %>% group_by(cluster))
wnt.markers<- as.data.frame(wnt.markers %>% group_by(cluster) %>% top_n(n =50, wt = avg_logFC))

### get top avg exp gene
Dista<-subset(wnt.markers,subset=cluster=='Dista')$gene
Proxi<-subset(wnt.markers,subset=cluster=='Proxi')$gene
Mix<-subset(wnt.markers,subset=cluster=='Mix')$gene

##get peaks closet to tss (distacne >3000)
peakk27<-readPeakFile('signature/H3K27me3_SET-seq_10^4_cells_Epi_peaks.broadPeak',header=F)
peakanno <- annotatePeak(peakk27,TxDb = txdb, tssRegion = c(-3000,3000),annoDb = 'org.Mm.eg.db')
Epi_k27_density<-as.data.frame(peakanno)
Epi_k27_density<-Epi_k27_density[,c(1,2,3,5,20,22)]
Epi_k27_density<-Epi_k27_density[abs(Epi_k27_density$distanceToTSS)>3000,]
names(Epi_k27_density)[1]<-'chr'
Epi_k27_density<-Epi_k27_density[,-c(4,5)]

peakk4<-readPeakFile('signature/H3K4me3_SET-seq_10^4_cells_Epi_peaks.narrowPeak',header=F)
peakanno <- annotatePeak(peakk4,TxDb = txdb, tssRegion = c(-3000,3000),annoDb = 'org.Mm.eg.db')
Epi_k4_density<-as.data.frame(peakanno)
Epi_k4_density<-Epi_k4_density[,c(1,2,3,5,20,22)]
Epi_k4_density<-Epi_k4_density[abs(Epi_k4_density$distanceToTSS)>3000,]
names(Epi_k4_density)[1]<-'chr'
Epi_k4_density<-Epi_k4_density[,-c(4,5)]

#load all peak in tss +/- 3kb (by bedtools intersect)
overlap<-read.table('signature/k27.overlap.bed',header = F,sep = '\t')
head(overlap)
overlap<-overlap[,c(1,2,3,16)]
head(overlap)
colnames(overlap)<-c('chr','start','end','SYMBOL')
head(overlap)

overlap<-read.table('signature/k4.overlap.bed',header = F,sep = '\t')
head(overlap)
overlap<-overlap[,c(1,2,3,16)]
head(overlap)
colnames(overlap)<-c('chr','start','end','SYMBOL')
head(overlap)

### merge two kinds of peaks
Epi_k27_density<-rbind(overlap,Epi_k27_density)
Epi_k27_density<-unique(Epi_k27_density)
Epi_k4_density<-rbind(overlap,Epi_k4_density)
Epi_k4_density<-unique(Epi_k4_density)
### get top exp peaks
Dista.use<-Epi_k27_density[Epi_k27_density$SYMBOL %in% Dista,]
Proxi.use<-Epi_k27_density[Epi_k27_density$SYMBOL %in% Proxi,]
Mix.use<-Epi_k27_density[Epi_k27_density$SYMBOL %in% Mix,]
write.table(Dista.use,'signature/Dista.k27peak',row.names = F,col.names = F,sep = '\t',quote = F)
write.table(Proxi.use,'signature/Proxi.k27peak',row.names = F,col.names = F,sep = '\t',quote = F)
write.table(Mix.use,'signature/Mix.k27peak',row.names = F,col.names = F,sep = '\t',quote = F)

Dista.use<-Epi_k4_density[Epi_k4_density$SYMBOL %in% Dista,]
Proxi.use<-Epi_k4_density[Epi_k4_density$SYMBOL %in% Proxi,]
Mix.use<-Epi_k4_density[Epi_k4_density$SYMBOL %in% Mix,]
write.table(Dista.use,'signature/Dista.k4peak',row.names = F,col.names = F,sep = '\t',quote = F)
write.table(Proxi.use,'signature/Proxi.k4peak',row.names = F,col.names = F,sep = '\t',quote = F)
write.table(Mix.use,'signature/Mix.k4peak',row.names = F,col.names = F,sep = '\t',quote = F)

### get Epi signals enrichments in each cell
library(plyr)
library(AUCell)
pred.matrix <- predictiveDistribution(Epi_k27)
path_to_signatures <- 'signature/'
Bulk_signatures <- paste(path_to_signatures, list.files(path_to_signatures,pattern = '.k27peak'), sep='')
labels  <- gsub('.k27peak', '', list.files(path_to_signatures,pattern = '.k27peak'))
Epi_k27 <- getSignaturesRegions(Epi_k27,signatures = Bulk_signatures,labels = labels,minOverlap = 0.4)
names(Epi_k27@signatures) <- labels
aucellRankings <- AUCell_buildRankings(pred.matrix, plot=T, verbose=FALSE)
Epi_k27 <- signatureCellEnrichment(Epi_k27, aucellRankings, selected.signatures='all', 
                                   aucMaxRank = 0.3*nrow(aucellRankings), plot=T,method='Umap')
signal<-Epi_k27@cell.data[,c(1,8,9,10)]
signal$sample<-rownames(signal)
saveRDS(signal,'scSET_Epi/k27signal.rds')

pred.matrix <- predictiveDistribution(Epi_k4)
path_to_signatures <- 'signature/'
Bulk_signatures <- paste(path_to_signatures, list.files(path_to_signatures,pattern = '.k4peak'), sep='')
labels  <- gsub('.k4peak', '', list.files(path_to_signatures,pattern = '.k4peak'))
Epi_k4 <- getSignaturesRegions(Epi_k4,signatures = Bulk_signatures,labels = labels,minOverlap = 0.4)
names(Epi_k4@signatures) <- labels
aucellRankings <- AUCell_buildRankings(pred.matrix, plot=T, verbose=FALSE)
Epi_k4 <- signatureCellEnrichment(Epi_k4, aucellRankings, selected.signatures='all', 
                                   aucMaxRank = 0.3*nrow(aucellRankings), plot=T,method='Umap')
signal<-Epi_k4@cell.data[,c(1,8,9,10)]
signal$sample<-rownames(signal)
saveRDS(signal,'scSET_Epi/k4signal.rds')


### mapping signals to use H3K4me3 data as example
signal <- readRDS('sc-scSET_Epi/k4signal.rds')
signal[,c('Mix','Proxi','Dista')][signal[,c('Mix','Proxi','Dista')]>0.3] =0.3 ### when do vlnplot and boxplot ,de annotation this, H3K27me3 signals use 0.5 as threshhold. 

Exp_k4@meta.data$sample<-rownames(Exp_k4@meta.data)
Exp_k4@meta.data<-merge(Exp_k4@meta.data,signal,by='sample')
rownames(Exp_k4@meta.data)<-Exp_k4@meta.data$sample

###draw Exp_Epi mapping UMAP clustering plots
saveRDS(Exp_k4,'EpiExp-k4.rds')
Exp_k4<-readRDS('EpiExp-k4.rds')
colorGradient <- colorRampPalette(c("blue","orange","red"))(30)
umap2<-DimPlot(Exp_k4,reduction = "umap",label=T,group.by = 'cell.type',pt.size = 3)
umapdata<-umap2$data
umapdata$sample<-rownames(umapdata)
umapdata<-merge(umapdata,signal,by='sample')
DimPlot(Exp_k4,reduction = "umap",label=T,pt.size = 3)
umap2

## Dista 
umap3<-ggplot(data=umapdata, aes(x=UMAP_1, y=UMAP_2, color=Dista))+geom_point(size=4)+
  scale_colour_gradientn(colours = c("grey","lightgoldenrod3","red"),breaks = c(0.1,0.2,0.3),limits=c(0,0.3))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border = element_blank())
umap3<-umap3+theme(legend.text = element_text(family = "Arial",size = 18,face = 'bold'),
                   legend.title = element_blank(),
                   axis.line = element_line(colour = "black",size = 2,lineend = "square"),
                   axis.title = element_text(family = "Arial",size = 20,face = 'bold'),
                   axis.text = element_text(family = "Arial",size = 20,face = 'bold',color = "black"),
                   axis.ticks.length = unit(.4,"lines"),
                   axis.ticks= element_line(size = 2))+ coord_fixed()
ggsave('k4_Dista.svg',plot = umap3,device = 'svg',path = 'final_plot',dpi=600)
### Proxi
umap4<-ggplot(data=umapdata, aes(x=UMAP_1, y=UMAP_2, color=Proxi))+geom_point(size=4)+
  scale_colour_gradientn(colours = c("grey","lightgoldenrod3","red"),breaks = c(0.1,0.2,0.3),limits=c(0,0.3))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border = element_blank())
umap4<-umap4+theme(legend.text = element_text(family = "Arial",size = 18,face = 'bold'),
                   #legend.title = element_text(size=10,face = 'bold'),
                   legend.title = element_blank(),
                   axis.line = element_line(colour = "black",size = 2,lineend = "square"),
                   axis.title = element_text(family = "Arial",size = 20,face = 'bold'),
                   axis.text = element_text(family = "Arial",size = 20,face = 'bold',color = "black"),
                   axis.ticks.length = unit(.4,"lines"),
                   axis.ticks= element_line(size = 2))+ coord_fixed()
ggsave('k4_Proxi.svg',plot = umap4,device = 'svg',path = 'final_plot',dpi=600)
### Mix
umap5<-ggplot(data=umapdata, aes(x=UMAP_1, y=UMAP_2, color=Mix))+geom_point(size=4)+
  scale_colour_gradientn(colours = c("grey","lightgoldenrod3","red"),breaks = c(0.1,0.2,0.3),limits=c(0,0.3))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border = element_blank())
umap5<-umap5+theme(legend.text = element_text(family = "Arial",size = 18,face = 'bold'),
                   #legend.title = element_text(size=10,face = 'bold'),
                   legend.title = element_blank(),
                   axis.line = element_line(colour = "black",size = 2,lineend = "square"),
                   axis.title = element_text(family = "Arial",size = 20,face = 'bold'),
                   axis.text = element_text(family = "Arial",size = 20,face = 'bold',color = "black"),
                   axis.ticks.length = unit(.4,"lines"),
                   axis.ticks= element_line(size = 2))+ coord_fixed()
ggsave('k4_Mix.svg',plot = umap5,device = 'svg',path = 'final_plot',dpi=600)

### draw vlnplot to quantify the signal distribution
## set comparison lists
my_comparisons <- list(c("Mix", "Dista"), c("Mix", "Proxi"), c("Dista", "Proxi"))
## remove outliers
idents<-data.frame(Idents(Exp_k4),names(Idents(Exp_k4)))
names(idents)<-c('cluster_name','sample')
Exp_k4@meta.data$sample<-rownames(Exp_k4@meta.data)
Exp_k4@meta.data<-merge(Exp_k4@meta.data,idents,by='sample')
rownames(Exp_k4@meta.data)<-Exp_k4@meta.data$sample
Exp_k4@meta.data<-Exp_k4@meta.data[,-1]
Exp_k4@meta.data
p<-ggplot(data = Exp_k4@meta.data,aes(x=cluster_name,y=Mix,fill=cluster_name)) +  
  geom_boxplot(size=2)  
out=layer_data(p)['outliers']
out=as.numeric(unlist(out))
pdata <- Exp_k4@meta.data
pdata <- pdata[! pdata$Mix %in% out,]
pobj<-Exp_k4
pobj@meta.data <- pdata
## draw vlnplots
p<-VlnPlot(pobj,features = c('Mix'),pt.size = 2,y.max = 0.75)+ 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",
                     face='bold',size=14,hide.ns = T,
                     bracket.size = 2,method = 't.test',p.adjust.methods='BH',vjust = 0.5,
                     label.y = c(0.45,0.55,0.65))
p$layers[[1]]$aes_params$size = 2 
p<- p+theme(axis.line = element_line(colour = "black",size = 2,lineend = "square"),
            axis.title = element_text(family = "Arial",size = 20,face = 'bold'),
            axis.text = element_text(family = "Arial",size = 20,face = 'bold',color = "black"),
            axis.text.x = element_text(angle = 0,hjust = 0.5),
            axis.ticks.length = unit(.4,"lines"),
            axis.ticks= element_line(size = 2))+guides(fill=F) + 
  scale_x_discrete(labels = c('Mix','Proxi','Dista'))
p
ggsave('k4_vln_Mix.svg',plot = p,device = 'svg',path = 'final_plot',dpi=600,height = 10,width = 7)
p<-ggplot(data = Exp_k4@meta.data,aes(x=cluster_name,y=Proxi,fill=cluster_name)) +  
  geom_boxplot(size=2)  
out=layer_data(p)['outliers']
out=as.numeric(unlist(out))
pdata <- Exp_k4@meta.data
pdata <- pdata[! pdata$Proxi %in% out,]
pobj<-Exp_k4
pobj@meta.data <- pdata
p<-VlnPlot(pobj,features = c('Proxi'),pt.size = 2,y.max = 0.75)+ 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",
                     face='bold',size=14,hide.ns = T,
                     bracket.size = 2,method = 't.test',p.adjust.methods='BH',vjust = 0.5,
                     label.y = c(0.45,0.55,0.65))
p$layers[[1]]$aes_params$size = 2
p<- p+theme(axis.line = element_line(colour = "black",size = 2,lineend = "square"),
            axis.title = element_text(family = "Arial",size = 20,face = 'bold'),
            axis.text = element_text(family = "Arial",size = 20,face = 'bold',color = "black"),
            axis.text.x = element_text(angle = 0,hjust = 0.5),
            axis.ticks.length = unit(.4,"lines"),
            axis.ticks= element_line(size = 2))+guides(fill=F) + 
  scale_x_discrete(labels = c('Mix','Proxi','Dista'))
p
ggsave('k4_vln_Proxi.svg',plot = p,device = 'svg',path = 'final_plot',dpi=600,height = 10,width = 7)
p<-ggplot(data = Exp_k4@meta.data,aes(x=cluster_name,y=Dista,fill=cluster_name)) +  
  geom_boxplot(size=2)  
out=layer_data(p)['outliers']
out=as.numeric(unlist(out))
pdata <- Exp_k4@meta.data
pdata <- pdata[! pdata$Dista %in% out,]
pobj<-Exp_k4
pobj@meta.data <- pdata
p<-VlnPlot(Exp_k4,features = c('Dista'),pt.size = 2,y.max = 0.75)+ 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",
                     face='bold',size=14,hide.ns = T,
                     bracket.size = 2,method = 't.test',p.adjust.methods='BH',vjust = 0.5,
                     label.y = c(0.45,0.55,0.65))
p$layers[[1]]$aes_params$size = 2
p<- p+theme(axis.line = element_line(colour = "black",size = 2,lineend = "square"),
            axis.title = element_text(family = "Arial",size = 20,face = 'bold'),
            axis.text = element_text(family = "Arial",size = 20,face = 'bold',color = "black"),
            axis.text.x = element_text(angle = 0,hjust = 0.5),
            axis.ticks.length = unit(.4,"lines"),
            axis.ticks= element_line(size = 2))+guides(fill=F) + 
  scale_x_discrete(labels = c('Mix','Proxi','Dista'))
p
ggsave('k4_vln_Dista.svg',plot = p,device = 'svg',path = 'final_plot',dpi=600,height = 10,width = 7)



