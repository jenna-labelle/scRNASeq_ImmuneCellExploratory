---
title: "R Notebook"
output: html_notebook
---
#Exploratory analysis of single cell RNASeq data
Data from

##Filtration, normalization, and Dimensionality reduction/clustering analyses performed (PCA, t-SNE)

Import packages
```{r}
library("org.Hs.eg.db") 
```

#Read in data, set up seurat object
```{r}
#Read in data
wd<- "C:/Users/jenna/OneDrive/Desktop/Seurat/"
dataInput<-read.table(file=paste(wd, "immune_control_expression_matrix.txt.gz",sep=""))
data<-dataInput[1:10000,]

#Convert to integer, if necessary
data_int<-apply(dataInput[,-1], 2, as.integer)
data<-as.data.frame(data_int)
rownames(data)<-dataInput$V1
data<-cbind(data, data)
data<-cbind(data, data)
colnames(data)<- paste("V", 1:ncol(data), sep ="")

#Set up seurat object
sdata<-CreateSeuratObject(counts=data, project="SC", min.cells=30)
sdata <- NormalizeData(sdata, verbose = FALSE)
```

#Visualize QC metrics, use for filtration
```{r}
VlnPlot(sdata, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
```
Will want to remove outlier cells


```{r}
sdata <- subset(sdata, subset = nFeature_RNA > 100)
VlnPlot(sdata, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
```

#Normalize
```{r}
sdata <- NormalizeData(sdata, normalization.method = "LogNormalize", scale.factor = 10000)
```

#Find variable features
```{r}
sdata <- FindVariableFeatures(sdata, selection.method = "vst", nfeatures = 1000)

#20 most variable features
top20 <- head(VariableFeatures(sdata), 20)

#Plot variable features, labeling top 20
plot <- VariableFeaturePlot(sdata)
LabelPoints(plot = plot, points = top20, repel = TRUE)

```

#Cell cycle: test for clustering based on likely cell cycle location of each cell
```{r}
#Read in list of cell cycle genes
cell_cycle <- read.csv(paste(wd, "CellCycleGenes.csv", sep=""))

#convert to gene symbols
cell_cycle$geneID <- unname(mapIds(org.Hs.eg.db, keys =cell_cycle$geneID, keytype = "ENSEMBL", column="SYMBOL"))

#Extract the G2/M genes
g2m_genes <- dplyr::filter(cell_cycle, phase == "G2/M") %>%
  pull(geneID) %>%
  as.character() 
  
#Extract the S genes
s_genes <- dplyr::filter(cell_cycle, phase == "S") %>%
  pull(geneID) %>%
  as.character() 

#Perform cell cycle scoring
sdata <- CellCycleScoring(sdata, g2m.features = g2m_genes, s.features = s_genes)

#Perform PCA and color by cell cycle phase
sdata = RunPCA(
  sdata,
  pc.genes = c(s_genes, g2m_genes),
  do.print = FALSE)

PCAPlot(sdata, 
        group.by= "Phase")
```
No clustering by cell cycle- don't need to include s.score and g2m.score in regression

#Dimensionality reduction
```{r}
#Apply linear transformation prior to PCA
sdata <- ScaleData(sdata, features = rownames(sdata))

#run pCA
sdata <- RunPCA(sdata, npcs = 30)

#plot PCA
DimPlot(sdata, reduction = "pca")
```
```{r}
#Plot heatmap of pcs
DimHeatmap(sdata, dims = 1:10, cells = 500, balanced = TRUE)
```
#Determine how many dimensions to use for further clustering
```{r}
#Using JackStraw resampling
sdata <- JackStraw(sdata, num.replicate = 100)
sdata <- ScoreJackStraw(sdata, dims = 1:20)

JackStrawPlot(sdata, dims = 1:20)
```
PCs 1-15 will be used for downstream clustering

#Perform additional clustering: UMAP and tSNE

#UMAP:
```{r}
sdata <- FindNeighbors(sdata, dims = 1:15)
sdata <- FindClusters(sdata, resolution = 0.5)

#Run UMAP clustering
sdata <- RunUMAP(sdata,  reduction = "pca", dims = 1:15)

DimPlot(sdata, reduction = "umap")
```
#tSNE clustering
```{r}
sdata <- RunTSNE(
  sdata,
  dims.use = 1:15,
  do.fast = TRUE)

#plot
DimPlot(sdata, reduction="tsne") 
```


#Find markers within each group
```{r}
# find all markers of cluster 1
cluster1.markers <- FindMarkers(sdata, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
```

```{r}
# find all markers distinguishing cluster 1 from cluster 0 
cluster1.markers <- FindMarkers(sdata, ident.1 = 1, ident.2 = 0, min.pct = 0.25)
head(cluster1.markers, n = 5)

```

#Look at expression for top 5 genes DE between groups 1 and 2
```{r}
VlnPlot(sdata, rownames(head(cluster1.markers, n = 5)))
```
#Visualize expression of top 5 DE genes in cluster 1 in cluster plot
```{r}
FeaturePlot(sdata,rownames(head(cluster1.markers, n = 5)))
```

#Find markers in all groups
```{r}
# find markers for every cluster compared to all remaining cells, report only the positive ones (top 2)
sdata.markers <- FindAllMarkers(sdata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sdata.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
```


#Rename clusters based on markers- if known
```{r}
new.cluster.ids <- c("Cell type 1", "Cell type 2")
names(new.cluster.ids) <- levels(sdata)
sdata <- RenameIdents(sdata, new.cluster.ids)
DimPlot(sdata, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```

