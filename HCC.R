.libPaths(c("~/SeuratV4",.libPaths()))
library(Seurat)
packageVersion("Seurat")

library(Seurat)
########fanjia数据#################
getwd()
setwd("D:/r/HCC_single_cell")
cellmetadata <- read.table("HCC_cell_metadata.txt",sep = "\t",header = T)
matrix <- read.table("HCC_log_tpm_expression_matrix.txt",sep = "\t",header = T)
rownames(matrix) <- matrix[,1]
matrix <- matrix[,-1]
head(matrix)
matrix[1:10,1:10]
cellmetadata <- cellmetadata[-1,]
head(cellmetadata)
rownames(cellmetadata)=cellmetadata[,1]
head(cellmetadata)
data=CreateSeuratObject(counts = matrix,project = "fanjia") 
fanjia <- data
View(data@meta.data)
a <- data@meta.data
a <- a[-1,]
data@meta.data <- a
View(data@meta.data)
save(data,file = "fanjia.rdata")
table(data@meta.data$orig.ident)
rm(data.big)

table(data@meta.data$orig.ident)
View(data@meta.data)

















############################GSE149614处理#####################
getwd()
setwd("D:/r/HCC_single_cell/GSE149614")
cellmetadata <- read.table("GSE149614_metadata.txt",sep = "\t",header = T)
matrix <- read.table("GSE149614_HCC.txt",sep = "\t",header = T)
cellmetadata <- read.table("HCC_cell_metadata.txt",sep = "\t",header = T)
matrix <- read.table("HCC_log_tpm_expression_matrix.txt",sep = "\t",header = T)
matrix <- matrix[-1,]
head(matrix)
matrix[1:10,1:10]

head(cellmetadata)
rownames(cellmetadata)=cellmetadata[,1]
cellmetadata <- cellmetadata[,-1]
head(cellmetadata)
data=CreateSeuratObject(counts = matrix,project = "GSE149614") 
View(data@meta.data)
save(data,file = "GSE149614.rdata")
table(data@meta.data$orig.ident)
############################GSE107747处理#####################
getwd()
setwd("D:/r/HCC_single_cell/GSE107747")
library(Seurat)
library(tidyverse)
library(reshape2)
library(SingleR)
library(BiocParallel)
library(patchwork)
library(SeuratWrappers)
library(batchelor)
library(SingleCellExperiment)
# rm(list = ls())
scRNAlist <- list()
samples=list.files("Merge1")
samples
dir <- file.path('Merge1',samples)
names(dir) <- samples
for(i in 1:length(dir)){
  print(i)
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, 
                                       min.cells = 3,min.features = 200)
}
scRNA2 <- merge(scRNAlist[[1]], y=c(scRNAlist[[2]]),
                add.cell.ids = names(scRNAlist))
dim(scRNA2)   #查看基因数和细胞总数
table(scRNA2@meta.data$orig.ident)  #查看每个样本的细胞数
save(scRNA2,file = 'GSE107747.rdata')
View(scRNA2@meta.data)
############################GSE162616处理#####################
getwd()
setwd("D:/r/HCC_single_cell/GSE162616")
library(Seurat)
library(tidyverse)
library(reshape2)
library(SingleR)
library(BiocParallel)
library(patchwork)
library(SeuratWrappers)
library(batchelor)
library(SingleCellExperiment)
rm(list = ls())
scRNAlist <- list()
samples=list.files("Merge1")
samples
dir <- file.path('Merge1',samples)
names(dir) <- samples
for(i in 1:length(dir)){
  print(i)
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, 
                                       min.cells = 3,min.features = 200)
}
scRNA2 <- merge(scRNAlist[[1]], y=c(scRNAlist[[2]],scRNAlist[[3]],scRNAlist[[4]],
                                    scRNAlist[[5]],scRNAlist[[6]],scRNAlist[[7]],
                                    scRNAlist[[8]],scRNAlist[[9]]),
                add.cell.ids = names(scRNAlist))
dim(scRNA2)   #查看基因数和细胞总数
table(scRNA2@meta.data$orig.ident)  #查看每个样本的细胞数
save(scRNA2,file = 'GSE162616.rdata')
View(scRNA2@meta.data)
############################GSE175793处理#####################
getwd()
setwd("D:/r/HCC_single_cell/GSE175793")
rm(list = ls())
scRNAlist <- list()
samples=list.files("Merge1")
samples
dir <- file.path('Merge1',samples)
names(dir) <- samples
for(i in 1:length(dir)){
  print(i)
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, 
                                       min.cells = 3,min.features = 200)
}
scRNA2 <- merge(scRNAlist[[1]], y=c(scRNAlist[[2]],scRNAlist[[3]],scRNAlist[[4]]),
                add.cell.ids = names(scRNAlist))
dim(scRNA2)   #查看基因数和细胞总数
table(scRNA2@meta.data$orig.ident)  #查看每个样本的细胞数
save(scRNA2,file = 'GSE175793.rdata')
View(scRNA2@meta.data)
############################GSE202642处理#####################
getwd()
setwd("D:/r/HCC_single_cell/GSE202642")
rm(list = ls())
scRNAlist <- list()
samples=list.files("Merge1")
samples
dir <- file.path('Merge1',samples)
names(dir) <- samples
counts <- Read10X(data.dir = dir[1])
scRNA2 <- CreateSeuratObject(counts,min.cells = 3,min.features = 200)
dim(scRNA2)   #查看基因数和细胞总数
table(scRNA2@meta.data$orig.ident)  #查看每个样本的细胞数
save(scRNA2,file = 'GSE202642.rdata')
View(scRNA2@meta.data)
###############单细胞数据集合并###############
library(Seurat)
table(GSE175793@meta.data$orig.ident)
rm(list = ls())
fanjia <- data
GSE107747 <-scRNA2
rm(data)
rm(scRNA2)
GSE149614 <- data
rm(data)
GSE162616 <- scRNA2
GSE175793 <- scRNA2
GSE101225 <- fanjia
data_merge <- merge(fanjia, 
                  y = c(GSE162616,GSE175793),
                  add.cell.ids = c("fanjia","GSE162616","GSE175793"), 
                  project = "data")
View(data_merge@meta.data)
unique(sapply(X = strsplit(colnames(data_merge), split = "_"), FUN = "[", 1))

######################降维聚类###################################
save(data_merge,file = "data_merge1.rdata")

#####################降维和聚类########################
.libPaths()
library(Seurat)
library(tidyverse)
library(reshape2)
library(SingleR)
library(celldex)
library(BiocParallel)
library(harmony)
library(Seurat)
library(tidyverse)
library(stringr)
library(patchwork)
library(clustree)
getwd()


#####################PCA聚类########################
??ScaleData
scRNA5 <- data_merge
scRNA5 <- ScaleData(scRNA5, do.center = T)
scRNA5 <- ScaleData(scRNA5, do.center = T)
scRNA5 <- FindVariableFeatures(scRNA5)

scRNA <- RunPCA(object = scRNA5,
                npcs = 50,
                rev.pca = FALSE,
                weight.by.var = TRUE,
                verbose = TRUE, 
                ndims.print = 1:5, 
                nfeatures.print = 30, 
                reduction.key = "PC_")
ElbowPlot(scRNA,
          ndims = 50) # 寻找对细胞差异贡献度较⼤的主成分




#热图的⾏是基因，列是细胞。dims = c(1:6)表示我只要展示前6
DimHeatmap(object = scRNA, 
           reduction = "pca", 
           dims = c(1:6),
           nfeatures = 30,
           disp.min = -2.5,
           balanced = TRUE,
           projected = FALSE,
           fast = TRUE,
           raster = TRUE,
           slot = "scale.data",
           combine = TRUE)


##查看批次效应
p1 <- DimPlot(object = scRNA, 
              reduction = "pca", 
              group.by = "orig.ident",
              dims = c(1,2),
              shuffle = TRUE,
              label = TRUE,
              label.size = 4,
              label.color = "black",
              label.box = TRUE,
              sizes.highlight = 1)
p1# 按照分组观察是否具有批次效应

p11 <- DimPlot(object = scRNA, 
               reduction = "pca", 
               group.by = "Group",
               dims = c(1,2),
               shuffle = TRUE,
               label = TRUE,
               label.size = 4,
               label.color = "black",
               label.box = TRUE,
               sizes.highlight = 1)
p11# 按照细胞周期阶段观察是否具有批次效应


p2 <- DimPlot(object = scRNA, 
              reduction = "pca",
              group.by = "Phase",
              dims = c(1,2),
              shuffle = TRUE,
              label = TRUE,
              label.size = 4,
              label.color = "black",
              label.box = TRUE,
              sizes.highlight = 1)
p2
p1+p2

###############tSNE分析##########################
scRNA <- RunTSNE(scRNA, 
                 dims = 1:20)
p3 <- DimPlot(object = scRNA, 
              reduction = "tsne", 
              group.by = "orig.ident",
              dims = c(1,2),
              shuffle = TRUE,
              label = TRUE,
              label.size = 4,
              label.color = "black",
              label.box = TRUE,
              sizes.highlight = 1)
p3
p33 <- DimPlot(object = scRNA, 
               reduction = "tsne", 
               group.by = "Group",
               dims = c(1,2),
               shuffle = TRUE,
               label = TRUE,
               label.size = 4,
               label.color = "black",
               label.box = TRUE,
               sizes.highlight = 1)
p33

p4 <- DimPlot(object = scRNA, 
              reduction = "tsne", 
              group.by = "Phase",
              dims = c(1,2),
              shuffle = TRUE,
              label = TRUE,
              label.size = 4,
              label.color = "black",
              label.box = TRUE,
              sizes.highlight = 1)
p4


####################用UMAP降维######################

scRNA <- RunUMAP(scRNA, 
                 dims = 1:20)
p5 <- DimPlot(object = scRNA, 
              reduction = "umap", 
              group.by = "orig.ident",
              dims = c(1,2),
              shuffle = TRUE,
              label = TRUE,
              label.size = 4,
              label.color = "black",
              label.box = TRUE,
              sizes.highlight = 1)
p5

p55 <- DimPlot(object = scRNA, 
               reduction = "umap", 
               group.by = "Group",
               dims = c(1,2),
               shuffle = TRUE,
               label = TRUE,
               label.size = 4,
               label.color = "black",
               label.box = TRUE,
               sizes.highlight = 1)
p55
p6 <- DimPlot(object = scRNA, 
              reduction = "umap", 
              group.by = "Phase",
              dims = c(1,2),
              shuffle = TRUE,
              label = TRUE,
              label.size = 4,
              label.color = "black",
              label.box = TRUE,
              sizes.highlight = 1)
p6

library(cowplot)
plot_grid(p2, p4, p6,
          labels = c("A", "B", "C"),
          nrow = 1,
          align = "h")
plot_grid(p1, p3, p5,
          labels = c("A", "B", "C"),
          nrow = 1,
          align = "h")
plot_grid(p11, p33, p55,
          labels = c("A", "B", "C"),
          nrow = 1,
          align = "h")

save(p1,p2,p3,p4,p5,p11,p33,p55,file = "UMAP.Rdata")
save(scRNA9_umaplog,file = "scRNA9_log")
save(scRNA8_tsnelog,file = "scRNA8_log")

FeaturePlot(scRNA9_log, 
            reduction = "umap", 
            features = "ACTA1")




#################单细胞聚类及可视化分析####################

#tSNE聚类
scRNA <- FindNeighbors(scRNA,
                       k.param = 20,
                       dims = 1:20)


scRNA <- FindClusters(scRNA,
                      resolution = 0.6,    #这个值需要自行调整设置越大，得到的cluster越多
                      method = "igraph", 
                      algorithm = 1, 
                      random.seed = 2022)
#tSNE可视化
DimPlot(object = scRNA,
        group.by = "seurat_clusters",
        reduction = "tsne",
        label = TRUE)


??FindClusters
#umap聚类
scRNA <- FindNeighbors(scRNA,
                       k.param = 20,
                       dims = 1:20)


scRNA <- FindClusters(scRNA,
                      resolution = 0.6,    #这个值需要自行调整设置越大，得到的cluster越多
                      method = "igraph", 
                      algorithm = 1, 
                      random.seed = 2022)
#umap可视化
DimPlot(object = scRNA,
        group.by = "seurat_clusters",
        reduction = "umap",
        label = TRUE)
DimPlot(object = scRNA,
        group.by = "seurat_clusters",
        reduction = "tsne",
        label = TRUE)
save(scRNA,file = "scRNA.rdata")














