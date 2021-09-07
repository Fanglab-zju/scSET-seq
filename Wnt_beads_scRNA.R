library(Seurat)
library(dplyr)
library(stringr)
######## scRNA-seq analysis by seurat

##### load the result from zUMIs
Wnt_beads_H3K27me3_scSET_1<-readRDS('zUMIs_results/Wnt_beads_H3K27me3_scSET-seq_Exp_1.dgecounts.rds')
Wnt_beads_H3K27me3_scSET_2<-readRDS('zUMIs_results/Wnt_beads_H3K27me3_scSET-seq_Exp_2.dgecounts.rds')
Wnt_beads_H3K27me3_scSET_3<-readRDS('zUMIs_results/Wnt_beads_H3K27me3_scSET-seq_Exp_3.dgecounts.rds')
Wnt_beads_H3K27me3_scSET_4 <- readRDS("zUMIs_results/Wnt_beads_H3K27me3_scSET-seq_Exp_4.dgecounts.rds")
Wnt_beads_H3K27me3_scSET_5 <- readRDS("zUMIs_results/Wnt_beads_H3K27me3_scSET-seq_Exp_5.dgecounts.rds")
Wnt_beads_H3K27me3_scSET_6 <- readRDS("zUMIs_results/Wnt_beads_H3K27me3_scSET-seq_Exp_6.dgecounts.rds")
Wnt_beads_H3K27me3_scSET_7 <- readRDS("zUMIs_results/Wnt_beads_H3K27me3_scSET-seq_Exp_7.dgecounts.rds")
Wnt_beads_H3K27me3_scSET_8 <- readRDS("zUMIs_results/Wnt_beads_H3K27me3_scSET-seq_Exp_8.dgecounts.rds")
Wnt_beads_H3K27me3_scSET_9 <- readRDS("zUMIs_results/Wnt_beads_H3K27me3_scSET-seq_Exp_9.dgecounts.rds")
Wnt_beads_H3K27me3_scSET_10 <- readRDS("zUMIs_results/Wnt_beads_H3K27me3_scSET-seq_Exp_10.dgecounts.rds")
Wnt_beads_H3K4me3_scSET_1 <- readRDS("zUMIs_results/Wnt_beads_H3K4me3_scSET-seq_Exp_1.dgecounts.rds")
Wnt_beads_H3K4me3_scSET_2 <- readRDS("zUMIs_results/Wnt_beads_H3K4me3_scSET-seq_Exp_2.dgecounts.rds")
Wnt_beads_H3K4me3_scSET_3 <- readRDS("zUMIs_results/Wnt_beads_H3K4me3_scSET-seq_Exp_3.dgecounts.rds")
Wnt_beads_H3K4me3_scSET_4 <- readRDS("zUMIs_results/Wnt_beads_H3K4me3_scSET-seq_Exp_4.dgecounts.rds")
Wnt_beads_H3K4me3_scSET_5 <- readRDS("zUMIs_results/Wnt_beads_H3K4me3_scSET-seq_Exp_5.dgecounts.rds")
Wnt_beads_H3K4me3_scSET_6 <- readRDS("zUMIs_results/Wnt_beads_H3K4me3_scSET-seq_Exp_6.dgecounts.rds")

###get expression matrix
Wnt_beads_H3K27me3_scSET_1<-Wnt_beads_H3K27me3_scSET_1$readcount$inex$all
Wnt_beads_H3K27me3_scSET_2<-Wnt_beads_H3K27me3_scSET_2$readcount$inex$all
Wnt_beads_H3K27me3_scSET_3<-Wnt_beads_H3K27me3_scSET_3$readcount$inex$all
Wnt_beads_H3K27me3_scSET_4<-Wnt_beads_H3K27me3_scSET_4$readcount$inex$all
Wnt_beads_H3K27me3_scSET_5<-Wnt_beads_H3K27me3_scSET_5$readcount$inex$all
Wnt_beads_H3K27me3_scSET_6<-Wnt_beads_H3K27me3_scSET_6$readcount$inex$all
Wnt_beads_H3K27me3_scSET_7<-Wnt_beads_H3K27me3_scSET_7$readcount$inex$all
Wnt_beads_H3K27me3_scSET_8<-Wnt_beads_H3K27me3_scSET_8$readcount$inex$all
Wnt_beads_H3K27me3_scSET_9<-Wnt_beads_H3K27me3_scSET_9$readcount$inex$all
Wnt_beads_H3K27me3_scSET_10<-Wnt_beads_H3K27me3_scSET_10$readcount$inex$all

Wnt_beads_H3K4me3_scSET_1<-Wnt_beads_H3K4me3_scSET_1$readcount$inex$all
Wnt_beads_H3K4me3_scSET_2<-Wnt_beads_H3K4me3_scSET_2$readcount$inex$all
Wnt_beads_H3K4me3_scSET_3<-Wnt_beads_H3K4me3_scSET_3$readcount$inex$all
Wnt_beads_H3K4me3_scSET_4<-Wnt_beads_H3K4me3_scSET_4$readcount$inex$all
Wnt_beads_H3K4me3_scSET_5<-Wnt_beads_H3K4me3_scSET_5$readcount$inex$all
Wnt_beads_H3K4me3_scSET_6<-Wnt_beads_H3K4me3_scSET_6$readcount$inex$all
### create seurat object
Wnt_beads_H3K27me3_scSET_1 <-CreateSeuratObject(counts = Wnt_beads_H3K27me3_scSET_1, project = "Wnt_beads_H3K27me3_scSET-seq_Exp_1",assay = "RNA")
Wnt_beads_H3K27me3_scSET_2 <-CreateSeuratObject(counts = Wnt_beads_H3K27me3_scSET_2, project = "Wnt_beads_H3K27me3_scSET-seq_Exp_2",assay = "RNA")
Wnt_beads_H3K27me3_scSET_3 <-CreateSeuratObject(counts = Wnt_beads_H3K27me3_scSET_3, project = "Wnt_beads_H3K27me3_scSET-seq_Exp_3",assay = "RNA")
Wnt_beads_H3K27me3_scSET_4 <-CreateSeuratObject(counts = Wnt_beads_H3K27me3_scSET_4, project = "Wnt_beads_H3K27me3_scSET-seq_Exp_4",assay = "RNA")
Wnt_beads_H3K27me3_scSET_5 <-CreateSeuratObject(counts = Wnt_beads_H3K27me3_scSET_5, project = "Wnt_beads_H3K27me3_scSET-seq_Exp_5",assay = "RNA")
Wnt_beads_H3K27me3_scSET_6 <-CreateSeuratObject(counts = Wnt_beads_H3K27me3_scSET_6, project = "Wnt_beads_H3K27me3_scSET-seq_Exp_6",assay = "RNA")
Wnt_beads_H3K27me3_scSET_7 <-CreateSeuratObject(counts = Wnt_beads_H3K27me3_scSET_7, project = "Wnt_beads_H3K27me3_scSET-seq_Exp_7",assay = "RNA")
Wnt_beads_H3K27me3_scSET_8 <-CreateSeuratObject(counts = Wnt_beads_H3K27me3_scSET_8, project = "Wnt_beads_H3K27me3_scSET-seq_Exp_8",assay = "RNA")
Wnt_beads_H3K27me3_scSET_9 <-CreateSeuratObject(counts = Wnt_beads_H3K27me3_scSET_9, project = "Wnt_beads_H3K27me3_scSET-seq_Exp_9",assay = "RNA")
Wnt_beads_H3K27me3_scSET_10 <-CreateSeuratObject(counts = Wnt_beads_H3K27me3_scSET_10, project = "Wnt_beads_H3K27me3_scSET-seq_Exp_10",assay = "RNA")
Wnt_beads_H3K4me3_scSET_1 <-CreateSeuratObject(counts = Wnt_beads_H3K4me3_scSET_1, project = "Wnt_beads_H3K4me3_scSET-seq_Exp_1",assay = "RNA")
Wnt_beads_H3K4me3_scSET_2 <-CreateSeuratObject(counts = Wnt_beads_H3K4me3_scSET_2, project = "Wnt_beads_H3K4me3_scSET-seq_Exp_2",assay = "RNA")
Wnt_beads_H3K4me3_scSET_3 <-CreateSeuratObject(counts = Wnt_beads_H3K4me3_scSET_3, project = "Wnt_beads_H3K4me3_scSET-seq_Exp_3",assay = "RNA")
Wnt_beads_H3K4me3_scSET_4 <-CreateSeuratObject(counts = Wnt_beads_H3K4me3_scSET_4, project = "Wnt_beads_H3K4me3_scSET-seq_Exp_4",assay = "RNA")
Wnt_beads_H3K4me3_scSET_5 <-CreateSeuratObject(counts = Wnt_beads_H3K4me3_scSET_5, project = "Wnt_beads_H3K4me3_scSET-seq_Exp_5",assay = "RNA")
Wnt_beads_H3K4me3_scSET_6 <-CreateSeuratObject(counts = Wnt_beads_H3K4me3_scSET_6, project = "Wnt_beads_H3K4me3_scSET-seq_Exp_6",assay = "RNA")

###
### load barcode information
i5<-read.table('i5.txt')
i5<-as.matrix(i5)
i7<-read.table('i7.txt')
i7<-as.matrix(i7)

### build bacrode-beads mapping function
beadsBC<-function(x,i5=i5,i7=i7){
  metadata <- x@meta.data
  metadata$cells <- rownames(metadata)
  metadata$cell.type <- NA
  i7<-as.matrix(read.table("i7.txt"))
  index1 = seq(1,11,2) 
  index2 = seq(2,12,2)
  for (i in index1) {
    metadata$cell.type[which(str_detect(metadata$cells, pattern = i7[i]))] <- "Proximal"
  }
  for (i in index2) {
    metadata$cell.type[which(str_detect(metadata$cells, pattern = i7[i]))] <- "Distal"
  }
  x@meta.data <- metadata
  return(x)
}

### mapping bacorde to beads states
Wnt_beads_H3K27me3_scSET_4<-beadsBC(Wnt_beads_H3K27me3_scSET_4)
Wnt_beads_H3K27me3_scSET_5<-beadsBC(Wnt_beads_H3K27me3_scSET_5)
Wnt_beads_H3K27me3_scSET_6<-beadsBC(Wnt_beads_H3K27me3_scSET_6)
Wnt_beads_H3K27me3_scSET_7<-beadsBC(Wnt_beads_H3K27me3_scSET_7)
Wnt_beads_H3K27me3_scSET_8<-beadsBC(Wnt_beads_H3K27me3_scSET_8)
Wnt_beads_H3K27me3_scSET_9<-beadsBC(Wnt_beads_H3K27me3_scSET_9)
Wnt_beads_H3K27me3_scSET_10<-beadsBC(Wnt_beads_H3K27me3_scSET_10)
Wnt_beads_H3K27me3_scSET_1<-beadsBC(Wnt_beads_H3K27me3_scSET_1)
Wnt_beads_H3K27me3_scSET_2<-beadsBC(Wnt_beads_H3K27me3_scSET_2)
Wnt_beads_H3K27me3_scSET_3<-beadsBC(Wnt_beads_H3K27me3_scSET_3)
Wnt_beads_H3K4me3_scSET_1<-beadsBC(Wnt_beads_H3K4me3_scSET_1)
Wnt_beads_H3K4me3_scSET_2<-beadsBC(Wnt_beads_H3K4me3_scSET_2)
Wnt_beads_H3K4me3_scSET_3<-beadsBC(Wnt_beads_H3K4me3_scSET_3)
Wnt_beads_H3K4me3_scSET_4<-beadsBC(Wnt_beads_H3K4me3_scSET_4)
Wnt_beads_H3K4me3_scSET_5<-beadsBC(Wnt_beads_H3K4me3_scSET_5)
Wnt_beads_H3K4me3_scSET_6<-beadsBC(Wnt_beads_H3K4me3_scSET_6)

### building function mapping barcode to cell_name
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

### mapping barcode to cell_name
cell.names <- colnames(Wnt_beads_H3K27me3_scSET_1)
cell.names<-changeBC(cell.names,i7,i5)
Wnt_beads_H3K27me3_scSET_1<-RenameCells(Wnt_beads_H3K27me3_scSET_1,new.names = cell.names)
cell.names <- colnames(Wnt_beads_H3K27me3_scSET_2)
cell.names<-changeBC(cell.names,i7,i5)
Wnt_beads_H3K27me3_scSET_2<-RenameCells(Wnt_beads_H3K27me3_scSET_2,new.names = cell.names)
cell.names <- colnames(Wnt_beads_H3K27me3_scSET_3)
cell.names<-changeBC(cell.names,i7,i5)
Wnt_beads_H3K27me3_scSET_3<-RenameCells(Wnt_beads_H3K27me3_scSET_3,new.names = cell.names)
cell.names <- colnames(Wnt_beads_H3K27me3_scSET_4)
cell.names<-changeBC(cell.names,i7,i5)
Wnt_beads_H3K27me3_scSET_4<-RenameCells(Wnt_beads_H3K27me3_scSET_4,new.names = cell.names)
cell.names <- colnames(Wnt_beads_H3K27me3_scSET_5)
cell.names<-changeBC(cell.names,i7,i5)
Wnt_beads_H3K27me3_scSET_5<-RenameCells(Wnt_beads_H3K27me3_scSET_5,new.names = cell.names)
cell.names <- colnames(Wnt_beads_H3K27me3_scSET_6)
cell.names<-changeBC(cell.names,i7,i5)
Wnt_beads_H3K27me3_scSET_6<-RenameCells(Wnt_beads_H3K27me3_scSET_6,new.names = cell.names)
cell.names <- colnames(Wnt_beads_H3K27me3_scSET_7)
cell.names<-changeBC(cell.names,i7,i5)
Wnt_beads_H3K27me3_scSET_7<-RenameCells(Wnt_beads_H3K27me3_scSET_7,new.names = cell.names)
cell.names <- colnames(Wnt_beads_H3K27me3_scSET_8)
cell.names<-changeBC(cell.names,i7,i5)
Wnt_beads_H3K27me3_scSET_8<-RenameCells(Wnt_beads_H3K27me3_scSET_8,new.names = cell.names)
cell.names <- colnames(Wnt_beads_H3K27me3_scSET_9)
cell.names<-changeBC(cell.names,i7,i5)
Wnt_beads_H3K27me3_scSET_9<-RenameCells(Wnt_beads_H3K27me3_scSET_9,new.names = cell.names)
cell.names <- colnames(Wnt_beads_H3K27me3_scSET_10)
cell.names<-changeBC(cell.names,i7,i5)
Wnt_beads_H3K27me3_scSET_10<-RenameCells(Wnt_beads_H3K27me3_scSET_10,new.names = cell.names)
cell.names <- colnames(Wnt_beads_H3K4me3_scSET_1)
cell.names<-changeBC(cell.names,i7,i5)
Wnt_beads_H3K4me3_scSET_1<-RenameCells(Wnt_beads_H3K4me3_scSET_1,new.names = cell.names)
cell.names <- colnames(Wnt_beads_H3K4me3_scSET_2)
cell.names<-changeBC(cell.names,i7,i5)
Wnt_beads_H3K4me3_scSET_2<-RenameCells(Wnt_beads_H3K4me3_scSET_2,new.names = cell.names)
cell.names <- colnames(Wnt_beads_H3K4me3_scSET_3)
cell.names<-changeBC(cell.names,i7,i5)
Wnt_beads_H3K4me3_scSET_3<-RenameCells(Wnt_beads_H3K4me3_scSET_3,new.names = cell.names)
cell.names <- colnames(Wnt_beads_H3K4me3_scSET_4)
cell.names<-changeBC(cell.names,i7,i5)
Wnt_beads_H3K4me3_scSET_4<-RenameCells(Wnt_beads_H3K4me3_scSET_4,new.names = cell.names)
cell.names <- colnames(Wnt_beads_H3K4me3_scSET_5)
cell.names<-changeBC(cell.names,i7,i5)
Wnt_beads_H3K4me3_scSET_5<-RenameCells(Wnt_beads_H3K4me3_scSET_5,new.names = cell.names)
cell.names <- colnames(Wnt_beads_H3K4me3_scSET_6)
cell.names<-changeBC(cell.names,i7,i5)
Wnt_beads_H3K4me3_scSET_6<-RenameCells(Wnt_beads_H3K4me3_scSET_6,new.names = cell.names)

### merge samples by their sequencing batches
wntdata<-c(Wnt_beads_H3K27me3_scSET_5,Wnt_beads_H3K27me3_scSET_6)
cell.ids=c('Wnt_beads_H3K27me3_scSET-seq_Exp_4','Wnt_beads_H3K27me3_scSET-seq_Exp_5','Wnt_beads_H3K27me3_scSET-seq_Exp_6')
mix1<-merge(Wnt_beads_H3K27me3_scSET_4, y = wntdata, add.cell.ids = cell.ids, project = "wnt3a-mix1")
wntdata<-c(Wnt_beads_H3K27me3_scSET_8,Wnt_beads_H3K27me3_scSET_9,Wnt_beads_H3K27me3_scSET_10)
cell.ids=c('Wnt_beads_H3K27me3_scSET-seq_Exp_7','Wnt_beads_H3K27me3_scSET-seq_Exp_8','Wnt_beads_H3K27me3_scSET-seq_Exp_9','Wnt_beads_H3K27me3_scSET-seq_Exp_10')
mix2<-merge(Wnt_beads_H3K27me3_scSET_7, y = wntdata, add.cell.ids = cell.ids, project = "wnt3a-mix2")
wntdata<-c(Wnt_beads_H3K27me3_scSET_2,Wnt_beads_H3K27me3_scSET_3)
cell.ids=c('Wnt_beads_H3K27me3_scSET-seq_Exp_1','Wnt_beads_H3K27me3_scSET-seq_Exp_2','Wnt_beads_H3K27me3_scSET-seq_Exp_3')
mix3<-merge(Wnt_beads_H3K27me3_scSET_1, y = wntdata, add.cell.ids = cell.ids, project = "wnt3a-mix3")
wntdata<-c(Wnt_beads_H3K4me3_scSET_2,Wnt_beads_H3K4me3_scSET_3,Wnt_beads_H3K4me3_scSET_4)
cell.ids=c('Wnt_beads_H3K4me3_scSET-seq_Exp_1','Wnt_beads_H3K4me3_scSET-seq_Exp_2','Wnt_beads_H3K4me3_scSET-seq_Exp_3','Wnt_beads_H3K4me3_scSET-seq_Exp_4')
mix4<-merge(Wnt_beads_H3K4me3_scSET_1, y = wntdata, add.cell.ids = cell.ids, project = "wnt3a-mix4")
wntdata<-c(Wnt_beads_H3K4me3_scSET_6)
cell.ids=c('Wnt_beads_H3K4me3_scSET-seq_Exp_5','Wnt_beads_H3K4me3_scSET-seq_Exp_6')
mix5<-merge(Wnt_beads_H3K4me3_scSET_5, y = wntdata, add.cell.ids = cell.ids, project = "wnt3a-mix5")

### use CCA to integrate samples from different batches
scRNAlist<-c(mix1,mix2,mix3,mix4,mix5)
for (i in 1:length(scRNAlist)) {
  scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = nFeature_RNA > 1000 & nFeature_RNA  <5000 )
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
  scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]], selection.method = "vst",nfeatures = 2000)
}

scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist,anchor.features = 2000,k.filter = 28,dims = 1:28,k.score = 28)
wnt <- IntegrateData(anchorset = scRNA.anchors)

### get differential genes between Proximal_beads cells and Distal_beads cells
library(DESeq2)
proximal<-subset(wnt, subset = cell.type=='Proximal' )
distal<-subset(wnt, subset = cell.type=='Distal' )
proximal<-proximal@assays$RNA@counts
distal<-distal@assays$RNA@counts
mycounts<-cbind(proximal,distal)
ncol(proximal)
ncol(distal)
condition <- factor(c(rep("Proximal",278),rep("Distal",267)), levels = c("Proximal","Distal"))
colData <- data.frame(row.names=colnames(mycounts), condition)
colData
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds)
dds
saveRDS(dds,'dds.rds')
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="condition")


### filter differential genes
dds<-readRDS('dds.rds')
res = results(dds, contrast=c("condition", "Proximal", "Distal"))
diff_res <-subset(res, pvalue < 0.05 & abs(log2FoldChange) > 2)
diff<-rownames(diff_res)
all<-rownames(wnt)
### scale the sc-RNA exp counts
wnt <- ScaleData(wnt,features = all)
### do dimension reduction
wnt <- RunPCA(wnt,features = intersect(all,diff))
DimPlot(wnt, reduction = "pca",group.by = "cell.type")
DimPlot(wnt, reduction = "pca",group.by = "orig.ident")
ElbowPlot(wnt)
### do umap clustering
wnt <- FindNeighbors(wnt,dims = 1:4)
wnt <- FindClusters(wnt, resolution =0.2)
table(wnt@meta.data$RNA_snn_res.0.2)
wnt <- FindClusters(wnt, resolution =0.6)
wnt <- RunUMAP(object = wnt,dims = 1:4)

### results of clustering
p<-DimPlot(wnt,reduction = "umap",label=F,pt.size = 4)+
           theme(legend.text = element_blank(),
           axis.line = element_line(colour = "black",size = 2,lineend = "square"),
           axis.title = element_text(family = "Arial",size = 20,face = 'bold'),
           axis.text = element_text(family = "Arial",size = 20,face = 'bold',color = "black"),
           axis.ticks.length = unit(.4,"lines"),
           axis.ticks= element_line(size = 2))+ coord_fixed()
p
ggsave('all_cluster_no_define.svg',plot = p,device = 'svg',path = 'final_plot',dpi=600)

### show cell type clustering plot (beads distribution)
p<-DimPlot(wnt,reduction = "umap",label=F,group.by = 'cell.type',pt.size = 4,
           cols = c('deeppink3','#00BFC4'))+
           theme(legend.text = element_blank(),
           axis.line = element_line(colour = "black",size = 2,lineend = "square"),
           axis.title = element_text(family = "Arial",size = 20,face = 'bold'),
           axis.text = element_text(family = "Arial",size = 20,face = 'bold',color = "black"),
           axis.ticks.length = unit(.4,"lines"),
           axis.ticks= element_line(size = 2))+ coord_fixed()
p
ggsave('all_cluster_cell.type.svg',plot = p,device = 'svg',path = 'final_plot',dpi=600)

### show orig.ident from different batches
p<-DimPlot(wnt,reduction = "umap",label=F,group.by = 'orig.ident',pt.size = 4)
p<-p+theme(legend.text = element_blank(),
           axis.line = element_line(colour = "black",size = 2,lineend = "square"),
           axis.title = element_text(family = "Arial",size = 20,face = 'bold'),
           axis.text = element_text(family = "Arial",size = 20,face = 'bold',color = "black"),
           axis.ticks.length = unit(.4,"lines"),
           axis.ticks= element_line(size = 2))+ coord_fixed()
p
ggsave('all_cluster_orig.ident.svg',plot = p,device = 'svg',path = 'final_plot',dpi=600)

### rename cells 
new.cluster.ids<-c("Mix","Proxi","Dista","Mix","Mix")
names(new.cluster.ids) <- levels(wnt)
wnt<-RenameIdents(wnt,new.cluster.ids)
p<-DimPlot(wnt,reduction = "umap",label=F,pt.size = 4)
p<-p+theme(legend.text = element_blank(),
           axis.line = element_line(colour = "black",size = 2,lineend = "square"),
           axis.title = element_text(family = "Arial",size = 20,face = 'bold'),
           axis.text = element_text(family = "Arial",size = 20,face = 'bold',color = "black"),
           axis.ticks.length = unit(.4,"lines"),
           axis.ticks= element_line(size = 2))+ coord_fixed()
p
ggsave('all_cluster_define.svg',plot = p,device = 'svg',path = 'final_plot',dpi=600)
DimPlot(wnt,reduction = "umap",label=T,group.by = 'cell.type',pt.size = 3)
DimPlot(wnt,reduction = "umap",label=T,group.by = 'orig.ident',pt.size = 3)
saveRDS(wnt,'wnt.rds')
wnt<-readRDS('wnt.rds')
### get filtered cell metadata
wnt@meta.data$Cluster <- wnt@active.ident
wnt@meta.data$UMAP_1 <- as.data.frame(wnt@reductions$umap@cell.embeddings)$UMAP_1
wnt@meta.data$UMAP_2 <-  as.data.frame(wnt@reductions$umap@cell.embeddings)$UMAP_2

metadata <- wnt@meta.data
metadata
write.table(metadata,'table/wnt_metadata.txt',sep = '\t',quote = F,col.names = T,row.names = T)

### add antibody information to object
Antibody<-rownames(wnt@meta.data)
for (i in ncol(wnt)) {
  Antibody[which(str_detect(Antibody, pattern = 'K27'))] <- "H3K27me3"
}
for (i in ncol(wnt)) {
  Antibody[which(str_detect(Antibody, pattern = 'K4'))] <- "H3K4me3"
}
wnt@meta.data$Antibody<-Antibody

### show cell number in different clusters
### split wnt object into H3K27me3 exp object and H3K4me3 exp object
Exp_k4<-subset(wnt,subset= Antibody=='H3K4me3')
length(Exp_k4@active.ident[Exp_k4@active.ident=='Dista'])
length(Exp_k4@active.ident[Exp_k4@active.ident=='Proxi'])
length(Exp_k4@active.ident[Exp_k4@active.ident=='Mix'])

Exp_k27<-subset(wnt,subset= Antibody=='H3K27me3')
length(Exp_k27@active.ident[Exp_k27@active.ident=='Dista'])
length(Exp_k27@active.ident[Exp_k27@active.ident=='Proxi'])
length(Exp_k27@active.ident[Exp_k27@active.ident=='Mix'])

### rename cells_name for mapping with Epi_data
cell.k4<-colnames(Exp_k4)
cell.k27<-colnames(Exp_k27)
cell.k4<-gsub(cell.k4,replacement = '',pattern = '_Exp')
cell.k27<-gsub(cell.k27,replacement = '',pattern = '_Exp')
Exp_k4<-RenameCells(Exp_k4,new.names = cell.k4)
Exp_k27<-RenameCells(Exp_k27,new.names = cell.k27)

### show clustering results in H3K27me3 Exp and H3K4me3 Exp data
umap1<-DimPlot(Exp_k27,reduction = "umap",label=F,pt.size = 4)
umap1<-umap1+theme(legend.text = element_blank(),
                   axis.line = element_line(colour = "black",size = 2,lineend = "square"),
                   axis.title = element_text(family = "Arial",size = 20,face = 'bold'),
                   axis.text = element_text(family = "Arial",size = 20,face = 'bold',color = "black"),
                   axis.ticks.length = unit(.4,"lines"),
                   axis.ticks= element_line(size = 2))+ coord_fixed()
umap1
ggsave('k27_cluster_define.svg',plot = umap1,device = 'svg',path = 'final_plot',dp3=600)
umap2<-DimPlot(Exp_k27,reduction = "umap",label=F,group.by = 'cell.type',pt.size = 4,
               cols = c('deeppink3','#00BFC4'))
umap2<-umap2+theme(legend.text = element_blank(),
                   axis.line = element_line(colour = "black",size = 2,lineend = "square"),
                   axis.title = element_text(family = "Arial",size = 20,face = 'bold'),
                   axis.text = element_text(family = "Arial",size = 20,face = 'bold',color = "black"),
                   axis.ticks.length = unit(.4,"lines"),
                   axis.ticks= element_line(size = 2))+ coord_fixed()
umap2
ggsave('k27_cluster_cell.type.svg',plot = umap2,device = 'svg',path = 'final_plot',dpi=600 )


umap1<-DimPlot(Exp_k4,reduction = "umap",label=F,pt.size = 4)
umap1<-umap1+theme(legend.text = element_blank(),
                   axis.line = element_line(colour = "black",size = 2,lineend = "square"),
                   axis.title = element_text(family = "Arial",size = 20,face = 'bold'),
                   axis.text = element_text(family = "Arial",size = 20,face = 'bold',color = "black"),
                   axis.ticks.length = unit(.4,"lines"),
                   axis.ticks= element_line(size = 2)) + coord_fixed()
umap1
ggsave('k4_cluster_define.svg',plot = umap1,device = 'svg',path = 'final_plot',dpi=600)
umap2<-DimPlot(Exp_k4,reduction = "umap",label=F,group.by = 'cell.type',pt.size = 4,
               cols = c('deeppink3','#00BFC4'))
umap2<-umap2+theme(legend.text = element_blank(),
                   axis.line = element_line(colour = "black",size = 2,lineend = "square"),
                   axis.title = element_text(family = "Arial",size = 20,face = 'bold'),
                   axis.text = element_text(family = "Arial",size = 20,face = 'bold',color = "black"),
                   axis.ticks.length = unit(.4,"lines"),
                   axis.ticks= element_line(size = 2))+ coord_fixed()
umap2
ggsave('k4_cluster_cell.type.svg',plot = umap2,device = 'svg',path = 'final_plot',dpi=600)

### save RDS
saveRDS(Exp_k4,'Exp_k4.rds')


### get marker gene
wnt.markers <- FindAllMarkers(wnt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(wnt.markers,'table/Wnt_beads_marker.txt',sep = '\t',col.names = T,row.names = T,quote = F)
saveRDS(Exp_k27,'Exp_k27.rds')


