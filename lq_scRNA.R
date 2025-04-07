getwd()
rm(list = ls())
setwd("/home/data/t030312/LY_work/LQ/scrna/")
.libPaths("~/R")
.libPaths("/home/data/refdir/Rlib")
library(Seurat)
library(tidyverse)
library(reshape2)
library(SingleR)
library(celldex)
library(BiocParallel)
rm(list = ls())
samples=list.files("Merge1")
samples
dir <- file.path('Merge1',samples)
names(dir) <- samples
#未分组合并方法
scRNAlist <- list()
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
save(scRNA2,file = 'merge.Rdata')


#NC组合并
samplesNC=list.files("NC")
samplesNC
dirNC <- file.path('NC',samples1)
names(dirNC) <- samplesNC
scRNAlistNC <- list()
for(i in 1:length(dirNC)){
  print(i)
  counts <- Read10X(data.dir = dirNC[i])
  scRNAlistNC[[i]] <- CreateSeuratObject(counts, # 表达矩阵，可以是稀疏矩阵，也可以是普通矩阵
                                         min.cells = 3,# 去除在⼩于3个细胞中表达的基因
                                         min.features = 200)# 去除只有200个以下基因表达的细胞
}
scRNANC <- merge(scRNAlistNC[[1]], y=c(scRNAlistNC[[2]], scRNAlistNC[[3]]),
                 add.cell.ids = names(scRNAlistNC))
dim(scRNANC)   #查看基因数和细胞总数
table(scRNANC@meta.data$orig.ident)  #查看每个样本的细胞数
save(scRNANC,file = 'scRNA_NC.orig.Rdata')
typeof(scRNA2) # "S4"
class(scRNA2) # "dgCMatrix"



# 计算稀疏矩阵的内存⼤⼩
sparse.size <- object.size(x = scRNA2)
sparse.size
# 查看稀疏矩阵数据
scRNA2@assays[["RNA"]][1:10, 1:4]

a <- scRNA2@assays[["RNA"]][3]

head(scRNA2@meta.data)





#Tu组合并
samplesTu=list.files("Tu")
samplesTu
dirTu <- file.path('Tu',samplesTu)
names(dirTu) <- samplesTu
scRNAlistTu <- list()
for(i in 1:length(dirTu)){
  print(i)
  counts <- Read10X(data.dir = dirTu[i])
  scRNAlistTu[[i]] <- CreateSeuratObject(counts, 
                                         min.cells = 3,min.features = 200)
}
scRNATu <- merge(scRNAlistTu[[1]], y=c(scRNAlistTu[[2]], scRNAlistTu[[3]]),
                 add.cell.ids = names(scRNAlistTu))
dim(scRNATu)   #查看基因数和细胞总数
table(scRNATu@meta.data$orig.ident)  #查看每个样本的细胞数
save(scRNATu,file = 'scRNA_Tu.orig.Rdata')
typeof(scRNATu) # "S4"
class(scRNATu) # "dgCMatrix"




########Seurat对象的处理########

#找表达矩阵
scRNA_counts <- scRNA2@assays[["RNA"]]@counts

#基因信息表就是矩阵的行名

#细胞信息表
Cell_information <- scRNA2@meta.data
View(scRNA2@meta.data)

sparse.size <- object.size(x = scRNA2)
sparse.size

scRNA2.orig <- scRNA2
save(scRNA2.orig, file = "scRNA2.orig.Rda")


scRNA2@meta.data



a <- scRNA2@meta.data %>%
  mutate(Group = if_else(str_sub(scRNA2@meta.data[,1],1,2) == "PB",
                         "PBS",
                         "TD1"))
scRNA2@meta.data <- a
scRNA2.orig <- scRNA2
save(scRNA2.orig, file = "scRNA2.orig.Rda")
save(scRNA2, file = "scRNA2.Rda")

scRNA2.orig@meta.data

#######单细胞数据预处理#######

#去除线粒体基因和极端值
load("scRNA2.Rda")

#当线粒体基因含量很⾼，⽽这时候细胞的转录本数⽬⼜⽐较少的话，
#我们就认为这个细胞可能发⽣了死亡，我们就要把这个细胞过滤掉

#找到线粒体基因

mito.genes <- str_subset(string = rownames(scRNA2),
                         pattern = "^mt-")
mito.genes


PercentageFeatureSet()  #可以用来富集分析

#用这种方法，可以求感兴趣的基因集，在每个细胞中的比例
scRNA2[["percent.mt"]] <- PercentageFeatureSet(scRNA2, 
                                               pattern = "^mt-")   #在细胞信息部分添加了线粒体比例

head(scRNA2[["percent.mt"]])


library(farver)
VlnPlot(object = scRNA2, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        group.by  = "Group",
        log = T,
        pt.size = 0) #表示点的大小

#####山峦图####
RidgePlot(object = scRNA2, 
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
          # log = T,
          ncol = 1,
          group.by  = "Group")


#####散点图########
p1 <- FeatureScatter(object = scRNA2, 
                     feature1 = "nCount_RNA", 
                     feature2 = "nFeature_RNA",
                     group.by = "Group")

p2 <- FeatureScatter(object = scRNA2, 
                     feature1 = "nFeature_RNA", 
                     feature2 = "percent.mt",
                     group.by = "Group")
p1 + p2


1

scRNA2

a1 <- scRNA2@meta.data[scRNA2@meta.data[["percent.mt"]]>25,]

scRNA3 <- subset(scRNA2,
                 subset = nFeature_RNA > 200 &
                   nFeature_RNA < 5000 &
                   percent.mt > -Inf & # 极端值
                   percent.mt < 10)

######数据的质控2: 消除细胞周期##########

#细胞周期打分
scRNA3 <- CellCycleScoring(scRNA3, 
                           s.features = cc.genes$s.genes, 
                           g2m.features = cc.genes$g2m.genes)
head(scRNA3@meta.data[ , 6:8])
library(ggplot2)
ggplot(data = scRNA3@meta.data, 
       aes(x = scRNA3$S.Score,
           y = scRNA3$G2M.Score,
           fill = factor(scRNA3$Phase))) +
  geom_point(shape = 21, 
             size = 5.5,
             alpha = 0.5) +
  ggtitle("Cell Cycle Score of Single Cell RNA") +
  theme(plot.title = element_text(hjust = 0.50,
                                  color="black",
                                  size = 15))
ggsave("Cell Cycle Score.pdf")
save(scRNA3, file = "scRNA3.Rda")
#消除细胞周期的影响
scRNA <- ScaleData(scRNA3, 
                   vars.to.regress = c("S.Score", "G2M.Score"),
                   features = rownames(scRNA3))
save(scRNA,file = "scRNA.Rda")

load("scRNA.Rda")
#此步骤耗时较长
view(scRNA@meta.data)
#提取数据预处理后的Tu组和NC组子集的三个样本
scRNA_T <- subset(scRNA, Group == "Tumor")
scRNA_N <- subset(scRNA, Group == "Control")


#######数据的批次矫正,方法一 anchor########


##先进行数据标准化,数据标准化和识别高变基因需要在不同的批次下面完成
##去除批间差需要合并后完成
#####多个Seurat对象的处理######
scRNA@meta.data[["orig.ident"]]

scRNA41 <- SplitObject(scRNA,
                       split.by = "orig.ident")
dim(scRNA41)
length(scRNA41)


for(i in 1:length(scRNA41)){
  scRNA41_1 <- scRNA41[[i]]
  # Log标准化
  scRNA41_1 <- NormalizeData(scRNA41_1,
                             normalization.method = "LogNormalize",
                             scale.factor = 100000)
  # 计算⾼度可变基因
  scRNA41_1 <- FindVariableFeatures(scRNA41_1, 
                                    selection.method = "vst", 
                                    nfeatures = 1000) 
  scRNA41[[i]] <- scRNA41_1    # 表示将处理好的Seurat对象⼜重新存回scRNA4中
}

####开始综合合并数据#####
#第一步是寻找锚定对象
# 多组数据整合
# 寻找锚定对象
AnchorSet <- FindIntegrationAnchors(object.list = scRNA41, 
                                    normalization.method = "LogNormalize") 
AnchorSet

save(AnchorSet,file = "AnchorSet1.Rdata")
scRNA5<- IntegrateData(anchorset = AnchorSet)  #根据锚定对象将数据整合到一起
DefaultAssay(scRNA5) <- "integrated"

save(scRNA5,file = "scRNA5.Rdata")

# scale去除批间差
scRNA5 <- ScaleData(scRNA5, 
                    features = rownames(scRNA5))

scRNA5_Anchor <- scRNA5
save(scRNA5_Anchor,file = "scRNA_Log.Rdata")



#######数据的批次矫正,方法二 SCT########

rm(list = ls())
load("scRNA3.Rda")

scRNA4 <- SplitObject(scRNA3,
                      split.by = "orig.ident")
view(scRNA@meta.data)

scRNA4_Cont <- scRNA4 [[1]]
scRNA4_Cont <- SCTransform(scRNA4_Cont, 
                           variable.features.n = 2000,
                           # vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
                           verbose = FALSE)
?SCTransform
#查看计算结果
HVFInfo(object = scRNA4_Cont)[1:15,1:3]
#可以看到2000个基因的计算结果列表
#提取top10的基因
top10 <- head(VariableFeatures(scRNA4_Cont), 10)

#进行可视化
plot1 <- VariableFeaturePlot(scRNA4_Cont)
LabelPoints(plot = plot1,
            points = top10,
            repel = TRUE)

#循环运行
##多个seruat
for(i in 1:length(scRNA4)){
  scRNA4[[i]] <- SCTransform(scRNA4[[i]],
                             variable.features.n = 2000,
                             vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
                             verbose = FALSE)
}
#数据整合
#选择用于整合的基因
features <- SelectIntegrationFeatures(object.list = scRNA4,
                                      nfeatures = 2000)
#准备整合
scRNA4 <- PrepSCTIntegration(object.list = scRNA4, 
                             anchor.features = features)
#寻找锚定基因
AnchorSet <- FindIntegrationAnchors(object.list = scRNA4, 
                                    reference = 1,
                                    normalization.method = "SCT", #参数修改
                                    anchor.features = features)
#整合数据
scRNA4_all <- IntegrateData(anchorset = AnchorSet, 
                            normalization.method = "SCT")
DefaultAssay(scRNA4_all) <- "integrated"

scRNA4_SCT <- scRNA4_all
save(scRNA4_SCT,file = "scRNA4_SCT.Rdata")


#####################降维和聚类########################

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
rm(list = ls())
load("scRNA6_SCT_final.Rdata")
dim(scRNA6_SCT)
load("scRNA41_SCT.Rdata")


#####################PCA聚类########################
??ScaleData
scRNA5 <- ScaleData(scRNA5_Anchor, do.center = T)
scRNA <- RunPCA(object = scRNA5,
                npcs = 50,
                rev.pca = FALSE,
                weight.by.var = TRUE,
                verbose = TRUE, #shuc
                ndims.print = 1:5, 
                nfeatures.print = 30, 
                reduction.key = "PC_")
ElbowPlot(scRNA,
          ndims = 50) # 寻找对细胞差异贡献度较⼤的主成分
save(scRNA,file = "scRNA.rdata")



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

p <- p1+p11+p2
p
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
p <- p3+p33+p4
p
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
plot_grid(p5, p55, p6,
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
save(scRNA,file = "scRNA.rdata")
#⽐如原来我们聚类后发现细胞虽然是不同cluster，但是他们在降维层⾯上分不开，
#那你就把这个值设置⼀下，看看能不能把那⼏个分不开的归
#为同⼀类细胞，然后再将这⼀类细胞提取出来重新跑前⾯流程细分亚型。如果这会你的细胞数⽬⾮常多，有上万
#个，那没必要⼀下⼦得到太多的Cluster对吧，可以先看看整体结果，然后再进⼀步进⾏细化分析，这时候你就可
#以将这个参数设置得⼩⼀些，得到的cluster数⽬也少⼀些，⽅便你后续再做其他分析。


table(scRNA13_SCT@meta.data$seurat_clusters)
table(scRNA11_SCT@meta.data$seurat_clusters)
save(scRNA13_SCT,file = "umap_all_cluster.Rda")
save(scRNA11_SCT,file = "tsne_all_cluster.Rda")


# 从seurat对象中提取身份和样本信息，以确定每个类群中每个样本的细胞数
n_cells <- FetchData(scRNA11_SCT, 
                     vars = c("ident","orig.ident"))%>%
  dplyr::count(ident,orig.ident)  %>%
  tidyr::spread(ident, n)

# 查看表格
View(n_cells)
#提取子集

DimPlot(subset(scRNA13_SCT, Group == "Tumor"), 
        reduction = "umap", # pca, umap, tsne
        group.by = "seurat_clusters",
        label = T)
DimPlot(scRNA13_SCT, 
        reduction = "umap", # pca, umap, tsne
        group.by = "Group",
        label = T)


#######细胞注释##########
rm(list = ls())
library(tidyverse)
library(reshape2)
library(patchwork)
library(Seurat)
load(file = "tsne_all_cluster.Rda")#scRNA11_SCT
load(file = "umap_all_cluster.Rda")#scRNA13_SCT
# 从CellMarker下载细胞注释表格
# CellMarker：http://bio-bigdata.hrbmu.edu.cn/CellMarker/index.jsp
dev.off()
dev.new()
cell_marker <- cell_marker_ly
?separate_rows
cell_marker <- cell_marker %>%
  separate_rows(cellMarker, sep = ",") %>%
  na.omit() %>%
  distinct() %>% 
  arrange(cellName) 
save(cell_marker,file = "cell_marker.csv")
cell_marker$cellMarker <- str_trim(cell_marker$cellMarker, "left") 
p1 <- DimPlot(scRNA11_SCT, 
              reduction = "tsne", # pca, umap, tsne
              group.by  = "seurat_clusters",
              label = T)
p2 <- DotPlot(scRNA11_SCT, 
              features = unique(cell_marker$cellMarker)) +
  theme(axis.text = element_text(size = 8,
                                 angle = 45,
                                 hjust = 1))
p2
p3 <- FeaturePlot(scRNA11_SCT, 
                  reduction = "tsne", 
                  features = c("COL1A1","COL1A2","COL3A1","SPARC","DCN","FSTL1",
                               "CFD","FNDC1","ABCA10","COL5A1","GSN","DPT"), 
                  label = TRUE) 
p3
p1 + p3
range(scRNA11_SCT@assays$integrated@scale.data)

##############Marker核密度图###########
#安装包
.libPaths("/home/data/t030312/R/x86_64-conda-linux-gnu-library/4.1")
.libPaths("/home/data/t030312/miniconda3/envs/myenv/lib/R/library")
.libPaths("/home/data/refdir/Rlib")
library(rgeos)
library(SeuratObject)
library(Seurat)
library(Nebulosa)
library(farver)
plot_density(scRNA13_SCT, "CD4")
plot_density(scRNA11_SCT, c("COL1A1","COL1A2","COL3A1","SPARC","DCN","FSTL1",
                            "CFD","FNDC1","ABCA10","COL5A1","GSN","DPT"))
#"CD79A","PECAM1","DCN","ACTA2","MS4A2","MYLK",
#"CD68","MYH11","PMP2","S100A9","RGS5","CD3E"

#11群为myocytes，5群为SMC
# 添加 cell type
Cell_type <- c("0" = "Fibroblast Cell",
               "1" = "Endothelial Cell",
               "2" = "Fibroblast Cell",
               "3" = "Pericyte",
               "4" = "Myocyte",
               "5" = "Myocyte",
               "6" = "Fibroblast Cell",
               "7" = "Lymphoid Cell",
               "8" = "Endothelial Cell",
               "9" = "Lymphoid Cell",
               "10"= "Fibroblast Cell",
               "11"= "Myocyte",                                                                                                                        
               "12"= "Myeloid Cell",
               "13"= "Myeloid Cell",
               "14"= "Myeloid Cell",
               "15"= "Fibroblast Cell",
               "16"= "Myofibroblast",
               "17"= "Neuron",
               "18"= "Myeloid Cell",
               "19"= "Fibroblast Cell",
               "20"= "Endothelial Cell",
               "21"= "Endothelial Cell",
               "22"= "Lymphoid Cell")
scRNA13_SCT[['cell_type']] <- unname(Cell_type[scRNA13_SCT@meta.data$seurat_clusters])
save(scRNA13_SCT,file = "scRNA13_SCT.Rdata")
DimPlot(scRNA13_SCT, 
        reduction = "umap", 
        group.by = "cell_type",
        label = TRUE, 
        pt.size = 0.2) + NoLegend()
library(tidyverse)
tsne = scRNA11_SCT@reductions$tsne@cell.embeddings %>%
  as.data.frame() %>% cbind(tx = scRNA11_SCT@meta.data$seurat_clusters)

## 用ggplot绘图
jpeg(file = "tSNE.by_cluster.tiff", width = 15, height = 10, units = "in", res = 600)

ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, color = tx)) + 
  geom_point(size = 0.5, alpha = 2) +
  scale_color_manual(values=c("0" = "green3",
                              "1" = "red", 
                              "2" = "green3",
                              "3" = "yellow2", 
                              "4" = "darkred", 
                              "5" = "darkred",
                              "6" = "green3",
                              "7" = "hotpink3",
                              "8" = "red",
                              "9" = "blue",
                              "10"= "green3",
                              "11"= "darkred",
                              "12"= "darkred",
                              "13"= "Aquamarine",
                              "14"= "YellowGreen",
                              "15"= "green3",
                              "16"= "Firebrick4",
                              "17"= "Tan3",
                              "18"= "Orange1",
                              "19"= "green3",
                              "20"= "red",
                              "21"= "Gold1",
                              "22"= "Tomato"))
dev.off() 
head(scRNA11_SCT@active.ident)

####改变active.ident

colnames(scRNA11_SCT@meta.data)
temp <- factor(scRNA11_SCT@meta.data$cell_type)
names(temp) <- rownames(scRNA11_SCT@meta.data)
scRNA11_SCT@active.ident <- temp
scRNA11_SCT@active.ident
scRNA11_SCT@meta.data

table(scRNA11_SCT@active.ident)
is.na(scRNA11_SCT@active.ident)
table(is.na(scRNA11_SCT@active.ident))
CDH15

colnames(scRNA11_SCT@meta.data)
################Vlnplot#########################
VlnPlot(scRNA11_SCT, 
        features = c("COL1A1","COL1A2","COL3A1","SPARC","DCN","FSTL1",
                     "CFD","FNDC1","ABCA10","COL5A1","GSN","DPT"),
        group.by = "cell_type",
        pt.size = 0,
        ncol = 3)
VlnPlot(scRNA11_SCT, 
        features = c("KLHL41","MYF6","NNMT","GPC1"),
        group.by = "cell_type",
        pt.size = 0,
        ncol = 2)
save(scRNA11_SCT, file = "scRNA11_tsne_all_cell.Rda")
view(scRNA13_SCT@meta.data)
top10 <- scRNA13_SCT
DoHeatmap(scRNA13_SCT, features = top10) + NoLegend()
View(scRNA13_SCT@meta.data)
################改变active.ident###########################
scRNA13_SCT@active.ident
colnames(scRNA13_SCT@meta.data)
temp <- factor(scRNA13_SCT@meta.data$cell_type)
names(temp) <- rownames(scRNA13_SCT@meta.data)
scRNA13_SCT@active.ident <- temp
scRNA13_SCT@active.ident
scRNA13_SCT@meta.data

table(scRNA11_SCT@active.ident)
is.na(scRNA11_SCT@active.ident)
table(is.na(scRNA11_SCT@active.ident))


colnames(scRNA11_SCT@meta.data)

#############寻找每个群的Marker Gene###################
#找每个cluster的marker
View(scRNA@meta.data)
Idents(scRNA) <- "seurat_clusters"
# remove.packages("tidyverse",.libPaths())
getwd()
setwd("/home/data/t030312/LY_work/LQ/scrna")
library(Seurat)
# ??JoinLayers
# scRNA <- JoinLayers(scRNA)
markers <- FindAllMarkers(scRNA,only.pos = TRUE,
                          min.pct = 0.1, # 过滤掉那些在25%以下细胞中检测到的基因
                          logfc.threshold = 0.25) 

head(markers)
getwd()
setwd("/home/data/t030312/LY_work/LQ/scrna")
library(readr)
write.csv(markers,file="markers.csv")
markers %>% group_by(cluster)%>%top_n(n=2,wt=avg_log2FC)
library(dplyr)
top10 <- markers%>%
  group_by(cluster)%>%
  top_n(n=10,wt=avg_log2FC)
write.csv(top10,file = "top10.csv")
head(top10)

#画热图
.libPaths("/home/data/refdir/Rlib")
library(farver)
DoHeatmap(scRNA,features = top10$gene,
          group.by = "seurat_clusters")+NoLegend()




tp <- scRNA13_SCT@meta.data
scRNA@active.ident
save(scRNA11_SCT,file = "scRNA11_SCT.Rdata")

t <- c("a","b","p","MSC","tumor")
t[i]
a <- FindMarkers(scRNA13_SCT,
                 ident.1 = 11,
                 #ident.2 = 20,
                 min.pct = 0.5, # 过滤掉那些在50%以下细胞中检测到的基因
                 logfc.threshold = log(1))
view(a)
write.csv(a,file = "cluster11.csv")
###SingleR注释###
.libPaths("/home/data/t030312/miniconda3/envs/myenv/lib/R/library")
.libPaths("/home/data/t030312/R/x86_64-conda-linux-gnu-library/4.1")
.libPaths("/home/data/refdir/Rlib")
library(SingleR)
library(celldex)
library(BiocParallel)
library(ggplot2)
library(reshape2)
library(ExperimentHub)
library(pheatmap)
library(viridis)
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se <- MouseRNAseqData()
help("Deprecated") 
# celldex: https://github.com/LTLA/celldex/tree/master/R
# ref <- BlueprintEncodeData() 
# ref <- MouseRNAseqData() 
# ref <- ImmGenData() 
# ref <- DatabaseImmuneCellExpressionData() 
# ref <- NovershternHematopoieticData() 
# ref <- MonacoImmuneData()
hpca.se <- ref
save(ref, file = "HumanPrimaryCellAtlasData.Rda")
meta=scRNA@meta.data 
head(meta)
plot1 <- DimPlot(scRNA, reduction = "umap", label = TRUE)
plot2<-DimPlot(scRNA, reduction = "tsne",
               label = TRUE)
plot1 + plot2
load("HumanPrimaryCellAtlasData.Rda")
#进行singleR注释
pbmc_for_SingleR <- GetAssayData(scRNA, slot="data") ##获取标准化矩阵
pbmc.hesc <- SingleR(test = pbmc_for_SingleR, ref = hpca.se, labels = hpca.se$label.main) #
pbmc.hesc
table(pbmc.hesc$labels,meta$seurat_clusters)
#seurat 和 singleR的table表
table(pbmc.hesc$labels,meta$seurat_clusters)


scRNA13_SCT@meta.data$labels <-pbmc.hesc$labels
print(DimPlot(scRNA13_SCT, group.by = c("seurat_clusters", "labels"),reduction = "umap"))
#基于scores within cells
print(plotScoreHeatmap(pbmc.hesc))


#基于 per-cell “deltas”诊断
plotDeltaDistribution(pbmc.hesc, ncol = 3)

#与cluster结果比较
tab <- table(label = pbmc.hesc$labels,cluster = meta$seurat_clusters)
pheatmap(log10(tab + 10))

#使用多个数据库注释
#使用BP和HPCA两个数据库综合注释，使用list函数读入多个数据库
scRNA13_SCT <- pbmcpbmc3.hesc <- SingleR(test = pbmc_for_SingleR, ref = list(BP=bpe.se, HPCA=hpca.se), 
                                         labels = list(bpe.se$label.main, hpca.se$label.main)) 
table(pbmc3.hesc$labels,meta$seurat_clusters)
pbmc3@meta.data$labels <-pbmc3.hesc$labels
print(DimPlot(pbmc3, group.by = c("seurat_clusters", "labels"),reduction = "umap"))

#umap
singler <- SingleR(scRNA13_SCT@assays$RNA@data,
                   ref = ref,
                   labels = ref$label.main,
                   clusters = scRNA13_SCT@meta.data$seurat_clusters,
                   fine.tune = TRUE,
                   BPPARAM = SnowParam(8))
#tsne
singler <- SingleR(scRNA11_SCT@assays$RNA@data,
                   ref = ref,
                   labels = ref$label.main,
                   clusters = scRNA13_SCT@meta.data$seurat_clusters,
                   fine.tune = TRUE)
#BPPARAM = SnowParam(8))
# 查看细胞类型
singler$pruned.labels
#绘制带cell label的tsne和umap图
cell_ids <- singler$pruned.labels
names(cell_ids) <- levels(scRNA11_SCT)
levels(scRNA11_SCT)
scRNA11_SCT <- RenameIdents(scRNA11_SCT, cell_ids)
levels(scRNA11_SCT)
TSNEPlot(object = scRNA11_SCT, pt.size = 0.5, label = TRUE)  
UMAPPlot(object = scRNA13_SCT, pt.size = 0.5, label = TRUE)
# 导出注释结果
clusterAnn <- singler %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  dplyr::select(id, labels)
write.csv(clusterAnn, "cluster_annotation.txt")
singler2 <- SingleR(test = srt3_all@assays$RNA@counts, 
                    ref = ref, 
                    labels = ref$label.main)
cellAnn <- singler2 %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  dplyr::select(id, labels)
write_tsv(cellAnn, "cell_annotation.txt")
# 将数据写⼊metadata⽂件
singleR_type <- srt3_all@active.ident
srt3_all[['singleR_type']] <- unname(singleR_type[rownames(srt3_all@meta.data)])


##############copyKAT#########################

library(Seurat)
library(copykat)
library(tidyverse)
rm(list = ls())
scRNA <- readRDS("scRNA.rds")
scRNA <- UpdateSeuratObject(scRNA)
scRNA_NC <- subset(scRNA, Group == "Control")
scRNA_Tu <- subset(scRNA, Group == "Tumor")
counts_Tu <- as.matrix(scRNA_Tu@assays$RNA@counts)
cnv_Tu <- copykat(rawmat=counts_Tu, ngene.chr=5, sam.name="SCLC", n.cores=54)
# ngene.chr参数是过滤细胞的一个标准，它要求被用来做CNV预测的细胞，一个染色体上至少有5个基因。
# sam.name定义样本名称 (sample name)，会给出来的文件加前缀
saveRDS(cnv_Tu, "cnv_Tu.rds")
counts_NC <- as.matrix(scRNA_NC@assays$RNA@counts)
cnv_NC <- copykat(rawmat=counts_NC, ngene.chr=5, sam.name="SCLC", n.cores=54)
# ngene.chr参数是过滤细胞的一个标准，它要求被用来做CNV预测的细胞，一个染色体上至少有5个基因。
# sam.name定义样本名称 (sample name)，会给出来的文件加前缀
saveRDS(cnv_NC, "cnv_NC.rds")
df <- scRNA13_SCT@meta.data
View(df)
df2 <- cnv_Tu$prediction
View(df2)
df3 <- cnv_NC$prediction
View(df3)
library(readr)
library(dplyr)

df_Tu <- cnv_Tu$prediction
df_Tu
scRNA_T <- subset(scRNA11_SCT, Group == "Tumor")
dim(scRNA_T)  #35422
dim(scRNA_N)  #33286
table(duplicated(df_Tu$cell.names))
table(df_Tu$copykat.pred)
df_Tu <- df_Tu[!duplicated(df_Tu$cell.names),]
rownames(df_Tu) <- df_Tu$cell.names
scRNA_Tu <- AddMetaData(scRNA_T, metadata = df_Tu)
meta_Tu <- scRNA_Tu$copykat.pred
view(meta_Tu)
DimPlot(scRNA_Tu, group.by = "copykat.pred")
DimPlot(scRNA_Tu,reduction = "tsne",group.by="copykat.pred",
        pt.size = 1) + scale_color_manual(values = c("red", "gray","gray"))
table(scRNA_Tu@meta.data$copykat.pred)
view(scRNA_Tu@meta.data)
Tu_cells <- subset(scRNA_Tu,copykat.pred =="aneuploid")

DimPlot(subset(scRNA_Tu, copykat.pred =="aneuploid"), 
        reduction = "tsne", # pca, umap, tsne
        group.by = "cell_type",
        label = T)



# SCLC_copykat_prediction.txt是上一步生成的结果
mallignant <- read.delim("SCLC_copykat_prediction.txt", row.names = 1)

# 原始count中有889个细胞，malignant只有824个，有65个细胞被ngene.chr=5过滤掉了
# 把细胞的良恶性信息加入metadata
scRNA <- AddMetaData(scRNA, metadata = mallignant)
p1 <- DimPlot(scRNA, group.by = "celltype", label = T) + ggsci::scale_color_lancet()
p2 <- DimPlot(scRNA, group.by = "copykat.pred") + scale_color_manual(values = c("red", "gray"))
pc <- p1 + p2
ggsave("pred_mallignant.pdf", pc, width = 12, height = 5)
scRNA$celltype2 <- scRNA$celltype
scRNA$celltype2[scRNA$copykat.pred=="aneuploid"] <- "mallignant_cell"
table(scRNA$celltype2)
head(cnv_Tu)
Tumor_cells <- subset(cnv_Tu@prediction,copykat.pred=="aneuploid")
save(scRNA13_SCT,file = "scRNA13_SCT.Rdata")
############细分群#############################
head(scRNA11_SCT@meta.data)
subset <- subset(scRNA11_SCT,cell_type=="Fibroblast cells")
head(subset)





















###################GO分析###################################
options(stringsAsFactors = F)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
scRNA11_SCT@meta.data
df1 <- subset(cluster_marker,cluster=="16")
genelist <- bitr(df1$gene,fromType = "SYMBOL",
                 toType = "ENTREZID",OrgDb = 'org.Hs.eg.db')
#BP生物过程
ego_BP <- enrichGO(gene = genelist$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = T)
write.csv(ego_BP,file = "GO_BP_C16.csv")
barplot(ego_BP,showCategory = 30,color = "pvalue")
#CC细胞组分
ego_CC <- enrichGO(gene = genelist$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = T)
write.csv(ego_CC,file = "GO_CC_C17.csv")
barplot(ego_CC,showCategory = 30,color = "pvalue")
#分子功能
ego_MF <- enrichGO(gene = genelist$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = T)
write.csv(ego_MF,file = "GO_MF_C8.csv")
barplot(ego_MF,showCategory = 30,color = "pvalue")



###############拟时序分析#####################
.libPaths("/home/data/t030312/R/x86_64-conda-linux-gnu-library/4.1")
.libPaths("/home/data/refdir/Rlib")
library(farver)
library(monocle)
library(dplyr)
library(patchwork)
library(Seurat)
library(BiocGenerics)
library(Matrix)
library(stats4)
library(splines)
library(DDRTree)
library(irlba)
#monocle构建CDS需要3个矩阵：expr.matrix、pd、featuredata
load("scRNA11_SCT.Rdata")
View(scRNA11_SCT@meta.data)
scRNA11_SCT@active.ident
sub_cellgroup <- subset(scRNA13_SCT,cell_type=="Fibroblast Cell")
table(sub_cellgroup@meta.data$cell_type)
data <- as(as.matrix(sub_cellgroup@assays$RNA@counts), 'sparseMatrix')
#构建featuredata，一般featuredata需要两个col，一个是gene_id,一个是gene_short_name,row对应counts的rownames
fData <- data.frame(gene_short_name = row.names(data),
                    row.names = row.names(data))
#head(gene_ann)
pd <- new("AnnotatedDataFrame",data=sub_cellgroup@meta.data)
fd <- new("AnnotatedDataFrame",data=fData)
#构建matrix
ct=as(as.matrix(sub_cellgroup@assays$RNA@counts), 'sparseMatrix')     #单细胞counts矩阵
#构建monocle对象
sc_cds <- newCellDataSet(data,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
sc_cds <- detectGenes(sc_cds, min_expr = 1)
sc_cds <- sc_cds[fData(sc_cds)$num_cells_expressed > 10, ]#计算每个基因在多少细胞中表达
cds <- sc_cds
#估计size factor和离散度，类似于seurat的数据归一化处理
#大的数据集使用稀疏矩阵，节省内存，加快运算
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 3)
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10))
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds) 
plot_pc_variance_explained(cds, return_all = F)
cds <- reduceDimension(cds, max_components = 2, num_dim = 20,
                       reduction_method = 'DDRTree', verbose = T)
cds <- clusterCells(cds, num_clusters = 4) 
plot_cell_clusters(cds, 1, 2)
plot_cell_clusters(cds, 1, 2,color_by = "State",
                   markers =c("CSRP2","POSTN","COL1A1","ACTA2","COL3A1","EDNRA","LRRC15","CFD","MYC"))
pData(cds)$Cluster=pData(cds)$celltype
cds <- clusterCells(cds, num_clusters = 6)
plot_cell_clusters(cds, 1, 2, color_by = "CellType",
                   markers = c("CSRP2","POSTN","COL1A1","ACTA2","COL3A1","EDNRA","LRRC15","CFD","MYC"))
plot_cell_trajectory(cds,color_by = "CellType")
plot_cell_trajectory(cds,color_by = "State")
plot_cell_trajectory(cds,color_by = "Pseudotime")
plot_cell_trajectory(cds,color_by = "CellType")+facet_wrap(~celltype,nrow = 1)
diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr = "~Cluster")
sig_genes <- subset(diff_test_res, qval < 0.01)
sig_genes=sig_genes[order(sig_genes$pval),]
head(sig_genes[,c("gene_short_name", "pval", "qval")] ) 
cg=as.character(head(sig_genes$gene_short_name)) 
#  挑选差异最显著的基因可视化
plot_genes_jitter(cds[cg,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )
cg2=as.character(tail(sig_genes$gene_short_name)) 
plot_genes_jitter(cds[cg2,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )
# 第一步: 挑选合适的基因
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
ordering_genes
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)
#第二步降维
cds <- reduceDimension(cds, max_components = 3,
                       method = 'DDRTree')
# 第三步: 对细胞进行排序
cds <- orderCells(cds)
#可视化细胞分化轨迹
plot_cell_trajectory(cds, color_by = "Cluster")
plot_genes_in_pseudotime(cds[cg,],
                         color_by = "Cluster")
library(ggsci)
plot_cell_trajectory(cds, color_by = "Cluster")  + scale_color_nejm()
plot_cell_trajectory(cds, color_by = "State")  + scale_color_npg()

sub_cellgroup$Group
plot_cell_trajectory(cds, color_by = "Group")
plot_cell_trajectory(cds, color_by = "Group") +
  facet_wrap(~State, nrow = 1)
plot_cell_trajectory(cds, color_by = "Group") +
  facet_wrap(~Group, nrow = 1)
pData(cds)$CCN1 = log2(exprs(cds)['CCN1',]+1)
diff_test_res$id <- rownames(diff_test_res)
#library(dplyr)
diff_test_res %>% arrange(qval) %>% head(50) %>% select(id) -> gene_to_cluster
diff_test_res
rownames(diff_test_res)[1:30]
diff_test_res
#gene_to_cluster <- gene_to_cluster$id
my_pseudotime_cluster <- plot_pseudotime_heatmap(cds[gene_to_cluster$id,],
                                                 num_clusters = 5,
                                                 show_rownames = TRUE)
BEAM_res <- BEAM(cds, branch_point = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
table(BEAM_res$qval < 1e-10)
library(BiocGenerics)
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res, qval < 1e-20)),],
                            branch_point = 1,
                            num_clusters = 4,
                            use_gene_short_name = TRUE,
                            show_rownames = TRUE)

##################GSVA################################
.libPaths("/home/data/refdir/Rlib")
.libPaths("/home/data/t030312/miniconda3/envs/myenv/lib/R/library")
.libPaths("/home/data/t030312/R/x86_64-conda-linux-gnu-library/4.1")
.libPaths("/home/data/t030312/miniconda3/envs/myenv/lib")
library(Seurat)
library(rgeos)
library(tidyverse)
library(CellChat)
library(msigdbr)
library(GSVA)
scRNA11_SCT@meta.data
subexpr <- subset(scRNA11_SCT,cell_type=="Fibroblast cells")
expr=as.matrix(subexpr@assays$RNA@data)
##通路基因集
msgdC2 <- msigdbr(species = "Homo sapiens",category = "C2",subcategory = "KEGG")
keggSet <- msgdC2%>%split(x=.$gene_symbol,f=.$gs_description)
kegg <- gsva(expr,gset.idx.list = keggSet,kcdf="Gaussian",method="zscore",
             parallel.sz=1)
library(limma)
## limma gsva通路活性评估
de_gsva <- function(exprSet,meta,compare = NULL){
  allDiff = list()
  design <- model.matrix(~0+factor(meta))
  colnames(design)=levels(factor(meta))
  rownames(design)=colnames(exprSet)
  fit <- lmFit(exprSet,design)
  if(length(unique(meta))==2){
    if(is.null(compare)){
      stop("there are 2 Groups,Please set  compare value")
    }
    contrast.matrix<-makeContrasts(contrasts = compare,levels = design)
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    fit2 <- eBayes(fit2)
    tempOutput = topTable(fit2,adjust='fdr', coef=1, number=Inf)
    allDiff[[compare]] = na.omit(tempOutput)
    
  }else if(length(unique(meta))>2){
    for(g in colnames(design)){
      fm = ""
      for(gother in colnames(design)[which(!colnames(design) %in% g)]){
        fm = paste0(fm,"+",gother)
      } 
      
      fm = paste0(g,"VsOthers = ",g,"-(",substring(fm,2),")/",ncol(design)-1)
      contrast.matrix <- makeContrasts(contrasts = fm,levels=design)
      fit2 <- contrasts.fit(fit, contrast.matrix) 
      fit2 <- eBayes(fit2)
      allDiff[[g]]=topTable(fit2,adjust='fdr',coef=1,number=Inf)
    }
  }else{
    stop("error only have one group")
  }
  
  return(allDiff)
}
meta <- subexpr@meta.data[,c("Group")]
head(subexpr@meta.data)
Diff =de_gsva(exprSet = kegg ,meta = meta,compare = "Control")
head(expr)
write.csv(expr,file = "GSVA.csv")
#可视化
idiff <-Diff[["Control"]]
head(idiff)
df <- data.frame(ID = rownames(idiff), score = idiff$t )
head(df)
df$group =sapply(1:nrow(idiff),function(x){
  if(idiff[x,"logFC"]>0 & idiff[x,"adj.P.Val"]<Padj_threshold){return("up")}
  else if(idiff[x,"logFC"]<0 & idiff[x,"adj.P.Val"]<Padj_threshold){return("down")}
  else{return("noSig")}
})
Padj_threshold
# 按照score排序
df$hjust = ifelse(df$score>0,1,0)
df$nudge_y = ifelse(df$score>0,-0.1,0.1)
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
limt = max(abs(df$score))
ggplot(sortdf, aes(ID, score,fill=group)) + 
  geom_bar(stat = 'identity',alpha = 0.7) + 
  scale_fill_manual(breaks=c("down","noSig","up"),
                    values = c("#008020","grey","#08519C"))+
  geom_text(data = df, aes(label = df$ID, y = df$nudge_y),
            nudge_x =0,nudge_y =0,hjust =df$hjust,
            size = tex.size)+
  labs(x = paste0(type," pathways"),
       y=paste0("t value of GSVA score\n",compare),
       title = title)+
  scale_y_continuous(limits=c(-limt,limt))+
  coord_flip() + 
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank())+
  theme(panel.border = element_rect(size = 0.6)
        #panel.border = element_blank()
  )+
  theme(plot.title = element_text(hjust = 0.5,size = 18),
        axis.text.y = element_blank(),
        axis.title = element_text(hjust = 0.5,size = 18),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = limt
  )



setwd("/home/data/t030312/LY_work/LQ/scrna")
View(scRNA@meta.data)
getwd()
setwd("/home/data/t030312/LY_work/LQ/scrna/B细胞分亚群")
rm(list = ls())
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
# slamf1
# 示例数据
scRNA <- scRNA_harmony
table(scRNA$seurat_clusters) #分组信息存在了Group
Idents(scRNA) <- "seurat_clusters"
# 1.提取数据
# 使用seurat自带的FetchData进行提取数据。注意！以下方式二选一，有所区别。方法一如果直接提取数据，画小提琴的时候会出现负值。使用方法二则不会。

interested_gene <- c("Slamf1") #这里面设置自己想要的基因名字

# 1. 提取基因表达数据
#使用FetchData直接提取，这和直接用seurat画图的数据是一样的
#方法一：直接提取
interested_gene <- c("Slamf1")
data_df <- Seurat::FetchData(scRNA,vars = interested_gene)
data_df <- as.data.frame(data_df)
#方法二：设置默认的数据对象，这样画出的没有负值
DefaultAssay(scRNA) <- "RNA"
interested_gene <- c("Slamf1")
data_df <- Seurat::FetchData(scRNA,vars = interested_gene)
data_df <- as.data.frame(data_df)

# 2. 整理数据为长格式
data_long <- data_df %>% 
  tibble::rownames_to_column("Cell") %>% 
  tidyr::pivot_longer(-Cell, names_to = "Gene", values_to = "Expression")
# 将分组信息（group名字）配对到表格里
data_long$seurat_clusters <- scRNA$seurat_clusters[data_long$Cell]
# 3. 使用ggplot2绘图
#--boxplot
plot_box <- ggplot(data_long, aes(x = Gene, y = Expression, fill = seurat_clusters)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(size=0.1) +   #size=散的大小,禁用这行代码就不显示散点
  labs(title = "Boxplot of gene expression", y = "Expression level") + 
  # scale_fill_manual(values = c("#00ABBD", "#FF9933")) + #分组配色
  theme_minimal() +
  theme(legend.position = "none")
print(plot_box)




aa <- as.character(scRNA_harmony@meta.data$seurat_clusters)    #得到所有细胞对应的细胞类型
head(aa)
        
class(aa)


names(aa) <- rownames(scRNA_harmony@meta.data$seurat_clusters)   #给aa中的每个元素命名，所命名为元素对应细胞名（Barcode）
head(aa)

class(aa)

#######改变基因active.ident############
scRNA_harmony@active.ident <- as.factor(aa)   #用新的细胞类型来替代原有的细胞类型，aa此时必须变成因子形式
head(scRNA_harmony@active.ident)   #查看是否变换成功，变换完成




library(Seurat)
human_data <- scRNA_harmony
Idents(scRNA_harmony) <- "Group"
DefaultAssay(human_data) <- "RNA"
all.markers  <- FindAllMarkers(human_data, 
                               only.pos = TRUE, 
                               min.pct = 0.1, 
                               logfc.threshold = 0.25)

###############基因表达量绘图################

gene <- c( "Ighd",
           "Pax5",
           "Bcl6",
           "Cxcr3",
           "Pdcd1lg2",
           "Ccr6",
           "Zeb2",
           "Bcl2",
           "Aim2",
           "Iglc1",
           "Tnfrsf13b",
           "Zbtb32",
           "Cd80",
           "Tnfrsf17",
           "Ighm",
           "Igha",
           "Il9r",
           "Spib",
           "Cd86",
           "Jchain",
           "Eaf2",
           "Ighg1",
           "Top2a")

# 选择需要的marker gene进行展示，
# 平均表达量使用seurat自带函数AverageExpression进行计算。
# 热图使用Complexheatmap做即可。


#计算平均表达量
gene_cell_exp <- AverageExpression(pbmc,
                                   features = gene,
                                   group.by = 'cell_type',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)

#complexheatmap作图
library(ComplexHeatmap)
#顶部细胞类型注释
df <- data.frame(colnames(gene_cell_exp))
colnames(df) <- 'clusters'
top_anno = HeatmapAnnotation(df = df,#细胞名/cluster
                             border = F,
                             show_annotation_name = T,
                             gp = gpar(col = 'black'),
                             col = list(clusters = c('Naïve B cell'="#FF34B3",
                                                  'Plasma'="#F6F5B4",
                                                  'sw-MBC'="#20B2AA",
                                                  'unsw-MBC'="#FFA500",
                                                  '4'="#ADFF2F",
                                                  '5'="#FF6A6A",
                                                  '6'="#9ECABE",
                                                  '7'="#CD5C5C",
                                                  '8'="#708090")))#颜色设置
cell_type_cols <- c(brewer.pal(9, "Set1"), "#FF34B3","#BC8F8F","#20B2AA",
                    "#00F5FF","#FFA500","#ADFF2F","#FF6A6A","#7FFFD4", 
                    "#AB82FF","#90EE90","#00CD00","#008B8B","#6495ED",
                    "#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00","#000000","#708090")  
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)
#数据标准化缩放一下
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
Heatmap(marker_exp,
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = T,
        column_title = NULL,
        heatmap_legend_param = list(
          title=' '), 
        col = colorRampPalette(c("#B22222","white","#00008B"))(100),
        border = 'black',
        rect_gp = gpar(col = "black", lwd = 1),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        top_annotation = top_anno)



#################批量GO分析####################
library(readr)
df <- markers
# rm(scRNA_harmony2)
table(scRNA@meta.data$seurat_clusters)
options(stringsAsFactors = F)
.libPaths("/home/data/t030312/miniconda3/envs/myenv/lib/R/library")
.libPaths("/home/data/refdir/Rlib")
#BiocManager::install("org.Hs.eg.db")
select_cluster <-  c(0:22)
getwd()
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(tidyverse)
for (i in 1:length(select_cluster)){
  select_cluster[i]
  print(select_cluster[i])
  df1 <- subset(df,cluster==select_cluster[i])
  genelist <- bitr(df1$gene,fromType = "SYMBOL",
                   toType = "ENTREZID",OrgDb = 'org.Mm.eg.db')
  #BP生物过程
  ego_BP <- enrichGO(gene = genelist$ENTREZID,
                     OrgDb = org.Mm.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     readable = T)
  tiff(filename = paste("./GO_BP_barplot/",select_cluster[i],".tiff",sep = ""),width = 1500,
       height = 1800,res = 200)
  name=barplot(ego_BP,showCategory = 20,color = "pvalue",
               title = paste("GO_analysis_of C",select_cluster[i],sep = ""))
  print(name)
  dev.off()
}
#CC细胞组分
for (i in 1:length(select_cluster)){
  select_cluster[i]
  print(select_cluster[i])
  df1 <- subset(df,cluster==select_cluster[i])
  genelist <- bitr(df1$gene,fromType = "SYMBOL",
                   toType = "ENTREZID",OrgDb = 'org.Hs.eg.db')
  #BP生物过程
  ego_CC <- enrichGO(gene = genelist$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     readable = T)
  tiff(filename = paste("./GO_CC_barplot/",select_cluster[i],".tiff",sep = ""),width = 1500,
       height = 1800,res = 200)
  name=dotplot(ego_CC,showCategory = 20,color = "pvalue")
  print(name)
  dev.off()
}

#分子功能
#MF细胞组分
for (i in 1:length(select_cluster)){
  select_cluster[i]
  print(select_cluster[i])
  df1 <- subset(df,cluster==select_cluster[i])
  genelist <- bitr(df1$gene,fromType = "SYMBOL",
                   toType = "ENTREZID",OrgDb = 'org.Hs.eg.db')
  #BP生物过程
  ego_MF <- enrichGO(gene = genelist$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.1,
                     qvalueCutoff = 0.1,
                     readable = T)
  tiff(filename = paste("./GO_MF_dotplot/",select_cluster[i],".tiff",sep = ""),width = 1500,
       height = 1800,res = 200)
  name=dotplot(ego_MF,showCategory = 20,color = "pvalue")
  print(name)
  dev.off()
}
###########GO分析#########################
df1 <- subset(markers,cluster=="22")
genelist <- bitr(df1$gene,fromType = "SYMBOL",
                 toType = "ENTREZID",OrgDb = 'org.Mm.eg.db')
#BP生物过程
ego_BP <- enrichGO(gene = genelist$ENTREZID,
                   OrgDb = org.Mm.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = T)
barplot(ego_BP,showCategory = 20,color = "pvalue")



#####################细胞数目统计图大群#################################
# rm(list = ls())
library(RColorBrewer)
table(scRNA@meta.data$Group)
scRNA
Idents(scRNA) <- "seurat_clusters"
table(scRNA$Group)#查看各组细胞数
prop.table(table(scRNA$seurat_clusters))
table(Idents(scRNA),scRNA$seurat_clusters)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(scRNA), scRNA$Group), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
sample_cols <- c(brewer.pal(9,"Set1"), "#FF34B3","#BC8F8F","#20B2AA",
                 "#00F5FF","#FFA500","#ADFF2F","#FF6A6A","#7FFFD4", 
                 "#AB82FF","#90EE90","#00CD00","#008B8B","#6495ED",
                 "#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00","#000000","#708090")  
ggplot(Cellratio)+ 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5)+
  theme_classic()+
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = sample_cols)+
  # scale_fill_brewer(palette = "Set2")+
  theme(panel.border = element_rect(fill=NA,color="black", 
                                    linewidth=0.5,linetype="solid"))


###################删除细胞语法##################
meta <- scRNA15filt@meta.data
scRNA15filt@meta.data <- meta
table(sc)
a <- meta %>%
  mutate(quality = case_when(Afp > 0 & cell_type == "B cell" ~ "A",
                             Afp > 0 & cell_type == "15" ~ "B",
                             Afp < 0 & seurat_clusters != "15" ~ "B"))
scRNA@meta.data <- a
                                                                                          
table(a$quality)
dim(scRNA@meta.data)
scRNAfilter <- subset(scRNA,quality%in%"B")
scRNA <- scRNAfilter
#########################提取单细胞基因表达量################
expr <- scRNA15filt@assays$RNA
gene_expression <- expr %>% .['Afp',] %>% as.data.frame() %>% t()
gene_expression <- as.data.frame(gene_expression)
colnames(gene_expression) <- 'Afp'
meta$Afp <- gene_expression$Afp
gene_expression$cell <- rownames(gene_expression)
gene_expression_sel <- gene_expression[which(gene_expression$PLAT > 0),]
sce_select <- sce[,rownames(gene_expression_sel)]
dim(gene_expression_sel) / dim(gene_expression)
######################李强细胞注释###########################
rm(list = ls())
table(pbmc@meta.data$seurat_clusters)
Cell_type <- c("0" = "Memory B cell",
               "1" = "Naïve B cell",
               "2" = "Memory B cell",
               "3" = "Plasma cell",
               "4" = "Memory B cell",
               "5" = "Memory B cell",
               "6" = "Naïve B cell",
               "7" = "Naïve B cell")
Cell_type <- c("0" = "Naïve B cell",
               "1" = "sw-MBC",
               "2" = "unsw-MBC",
               "3" = "sw-MBC",
               "4" = "Naïve B cell",
               "5" = "Plasma",
               "6" = "unsw-MBC",
               "7" = "sw-MBC")
table(pbmc@meta.data$cell_type)
pbmc[['cell_type']] <- unname(Cell_type[pbmc@meta.data$seurat_clusters])
library(ggplot2)
library(RColorBrewer)
cell_type_cols <- c(brewer.pal(9, "Set1"), "#FF34B3","#BC8F8F","#20B2AA",
                    "#00F5FF","#FFA500","#ADFF2F","#FF6A6A","#7FFFD4", 
                    "#AB82FF","#90EE90","#00CD00","#008B8B","#6495ED",
                    "#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00","#000000","#708090")  
# scRNA_umap1 <- subset(scRNA_umap,Group != c("PBMC of healthy donor"))
p1 <- DimPlot(pbmc, 
              label = T,
              reduction = "umap",
              group.by = "cell_type",
              # split.by = "Group",
              cols= cell_type_cols, #设置颜色
              pt.size = 0.5,#设置点的大小
              # ncol = 3,
              # label.size = 2,
              repel = T,
              raster = FALSE)+labs(x = "UMAP1", y = "UMAP2")+ 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))  #加边框
# dev.off()



p2 <- DimPlot(pbmc, 
              label = T,
              reduction = "tsne",
              group.by = "cell_type",
              # split.by = "Group",
              cols= cell_type_cols, #设置颜色
              pt.size = 0.5,#设置点的大小
              # ncol = 3,
              # label.size = 2,
              repel = T,
              raster = FALSE)+labs(x = "tSNE1", y = "tSNE2")+ 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p <- p1+p2
p
dev.off()

getwd()
setwd("/home/data/t030312/LY_work/LQ/scrna/B细胞分亚群")
scRNAB <- subset(scRNA,cell_type%in%"B cell")

save(pbmc,file = "pbmc.rdata")

library(harmony)
library(Seurat)
library(cowplot)
library(harmony)
load('data/pbmc_stim.RData') #加载矩阵数据
#在运行Harmony之前，创建一个Seurat对象并按照标准PCA进行分析。
# pbmc <- CreateSeuratObject(counts = cbind(stim.sparse, ctrl.sparse), project = "PBMC", min.cells = 5) %>%
#   Seurat::NormalizeData(verbose = FALSE) %>%
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
#   ScaleData(verbose = FALSE) %>%
#   RunPCA(pc.genes = pbmc@var.genes, npcs = 20, verbose = FALSE) #R语言中%>%的含义是什么呢，管道函数啦，就是把左件的值发送给右件的表达式，并作为右件表达式函数的第一个参数。
# pbmc@meta.data$stim <- c(rep("STIM", ncol(stim.sparse)),rep("CTRL", ncol(ctrl.sparse)))#赋值条件变量
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = scRNAB, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = scRNAB, features = "PC_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1,p2)
options(repr.plot.height = 2.5, repr.plot.width = 6)
scRNAB <- scRNAB %>%
  RunHarmony("orig.ident", plot_convergence = TRUE) #Harmony converged after 8 iterations
scRNAB@reductions$harmony#Harmory运行后的结果储存
# 使用Embeddings命令访问新的Harmony embeddings
harmony_embeddings <- Embeddings(scRNAB,'harmony')
harmony_embeddings[1:5, 1:5]
# 让我们查看确认数据集在Harmony运行之后的前两个维度中得到很好的整合。


options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = scRNAB, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = scRNAB, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1,p2)
# 
# Downstream analysis
# 
# 许多下游分析是在低维嵌入而不是基因表达上进行的。要使用校正后的Harmony embeddings而不是PC，
# 设置reduction ='harmony'。
scRNAB <- scRNAB %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5) %>%
  identity()
options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(scRNAB, reduction = "umap", group.by = "orig.ident", pt.size = .1, split.by = 'orig.ident')
# TSNE分析


scRNAB=RunTSNE(scRNAB,reduction = "harmony", dims = 1:20)
TSNEPlot(object = scRNAB, 
         pt.size = 0.5, 
         label = TRUE,
         split.by='orig.ident') 
# 两样本合并的TSNE和UMAP图


p1 <- DimPlot(scRNAB, reduction = "umap",pt.size = .1,label = TRUE)
p2 <- TSNEPlot(scRNAB, pt.size = .1, label = TRUE) 
p <- p1+p2
p
save(scRNAB,file = "scRNAB.rdata")
getwd()
setwd("/home/data/t030312/LY_work/LQ/scrna/B细胞分亚群")



library(Seurat)
data <- scRNAB
library(ggplot2)
library(RColorBrewer)
cell_type_cols <- c(brewer.pal(9, "Set1"), "#FF34B3","#BC8F8F","#20B2AA",
                    "#00F5FF","#FFA500","#ADFF2F","#FF6A6A","#7FFFD4", 
                    "#AB82FF","#90EE90","#00CD00","#008B8B","#6495ED",
                    "#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00","#000000","#708090")  
p1 <- DimPlot(scRNA123, 
        label = T,
        group.by = "cell_type",
        reduction = "umap",
        # split.by = "sample",
        cols= cell_type_cols, #设置颜色
        pt.size = 0.1,#设置点的大小
        #col = 3,
        label.size = 4,
        repel = T)+labs(x = "UMAP1", y = "UMAP2")+ 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))  #加边框

p2 <- DimPlot(scRNA123, 
              label = T,
              group.by = "cell_type",
              reduction = "tsne",
              # split.by = "sample",
              cols= cell_type_cols, #设置颜色
              pt.size = 0.1,#设置点的大小
              #col = 3,
              label.size = 4,
              repel = T)+labs(x = "UMAP1", y = "UMAP2")+ 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))  #加边框
p <- p1+p2
p



rm(list = ls())
getwd()
setwd("/home/data/t030312/LY_work/LQ/scrna")
###################删除细胞语法##################
meta <- scRNA@meta.data

a <- meta %>%
  mutate(quality = case_when(Afp > 0 & cell_type == "B cell" ~ "A",
                             Afp > 0 & cell_type == "15" ~ "B",
                             Afp < 0 & seurat_clusters != "15" ~ "B"))
scRNA@meta.data <- a

table(a$quality)
dim(scRNA@meta.data)
scRNAfilter <- subset(scRNA,quality%in%"B")
scRNA <- scRNAfilter
meta <- scRNA15filt@meta.data
scRNA15filt@meta.data <- meta
table(sc)
a <- meta %>%
  mutate(quality = case_when(Afp > 0 & cell_type == "B cell" ~ "A",
                             Afp > 0 & cell_type == "15" ~ "B",
                             Afp < 0 & seurat_clusters != "15" ~ "B"))
scRNA@meta.data <- a

table(a$quality)
dim(scRNA@meta.data)
scRNAfilter <- subset(scRNA,quality%in%"B")
scRNA <- scRNAfilter
#########################提取单细胞基因表达量################
expr <- scRNA@assays$RNA
gene_expression <- expr %>% .['Afp',] %>% as.data.frame() %>% t()
gene_expression <- as.data.frame(gene_expression)
colnames(gene_expression) <- 'Afp'
meta$Afp <- gene_expression$Afp

gene_expression <- expr %>% .['Alb',] %>% as.data.frame() %>% t()
gene_expression <- as.data.frame(gene_expression)
colnames(gene_expression) <- 'Alb'
meta$Alb <- gene_expression$Alb


gene_expression <- expr %>% .['Cd3g',] %>% as.data.frame() %>% t()
gene_expression <- as.data.frame(gene_expression)
colnames(gene_expression) <- 'Cd3g'
meta$Cd3g <- gene_expression$Cd3g



gene_expression <- expr %>% .['Cd45',] %>% as.data.frame() %>% t()
gene_expression <- as.data.frame(gene_expression)
colnames(gene_expression) <- 'Cd45'
meta$Cd45 <- gene_expression$Cd45
# gene_expression$cell <- rownames(gene_expression)
# gene_expression_sel <- gene_expression[which(gene_expression$PLAT > 0),]
# sce_select <- sce[,rownames(gene_expression_sel)]
# dim(gene_expression_sel) / dim(gene_expression)
getwd()
#############寻找每个群的Marker Gene###################
#找每个cluster的marker
View(pbmc@meta.data)
Idents(pbmc) <- "Group"

library(Seurat)
getwd()
setwd("/home/data/t030312/LY_work/LQ/scrna/B细胞分亚群")
getwd()
markers <- FindAllMarkers(pbmc,only.pos = TRUE,
                          min.pct = 0.1, # 过滤掉那些在25%以下细胞中检测到的基因
                          logfc.threshold = 0.25) 

head(markers)
getwd()
# setwd("/home/data/t030312/LY_work/LQ/scrna")
library(readr)
write.csv(markers,file="markers.csv")
markers %>% group_by(cluster)%>%top_n(n=2,wt=avg_log2FC)
library(dplyr)
top10 <- markers%>%
  group_by(cluster)%>%
  top_n(n=10,wt=avg_log2FC)
write.csv(top10,file = "top10.csv")
head(top10)

################2024-03-1-01分析##################
getwd()
setwd("/home/data/t030312/LY_work/LQ/scrna/B细胞分亚群")

Idents(pbmc) <- "cell_type"
source("DoheatmapPlot.R")
library(patchwork)
DefaultAssay(scRNAB) <- "RNA"

###############基因表达量绘图################

# gene <- c("Cd79a","Cd79b","Havcr1","Pdcd1","Tigit","Lag3","Havcr2",
#           # "Nt5e","Entpd1","Ctla4","Tgfb1","Tgfb2","Il10","Itgam",
#           "Ebi3","Cd1d1","Lta","Tim1","Tle3")

# 选择需要的marker gene进行展示，
# 平均表达量使用seurat自带函数AverageExpression进行计算。
# 热图使用Complexheatmap做即可。
gene_list <- c(             "Csf2rb",
                            "Stat6",
                            "Ciita",
                            "Irf1",
                            "Btla",
"Ifnar2",
"Myc",
"Place8",
"Cd9",
"Nt5e",
"Entpd1",
"Ctla4",
"Tgfb1",
"Il10",
"Bcl2",
"Bcl6",
"Hhex",
"Ebi3",
"Cd1d1",
"Lta",
"Gzmb",
"Tle3")
gene_list <- c("Slamf1")
pbmc1 <- subset(pbmc,Group%in%"TD1")
pbmc2 <- subset(pbmc,Group%in%"PBS")
#计算平均表达量
scRNAB <- pbmc
Idents(pbmc1) <- "cell_type"
Idents(pbmc2) <- "cell_type"
pbmc_PBS <- subset(pbmc,Group%in%"PBS")
pbmc_TD1 <- subset(pbmc,Group%in%"TD1")
gene_cell_exp <- AverageExpression(pbmc,
                                   features = gene_list,
                                   group.by = 'Group',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)

#complexheatmap作图
library(ComplexHeatmap)
#顶部细胞类型注释
df <- data.frame(colnames(gene_cell_exp))
colnames(df) <- 'clusters'
# top_anno = HeatmapAnnotation(df = df,#细胞名/cluster
#                              border = F,
#                              show_annotation_name = T,
#                              gp = gpar(col = 'black'),
#                              col = list(clusters = c('0'="#FF34B3",
#                                                      '1'="#F6F5B4",
#                                                      '2'="#20B2AA",
#                                                      '3'="#FFA500",
#                                                      '4'="#ADFF2F",
#                                                      '5'="#FF6A6A",
#                                                      '6'="#9ECABE",
#                                                      '7'="#8B008B",
#                                                      '8'="#FFA500")))
getwd()
# setwd("/home/data/t030312/LY_work/LQ/scrna/B细胞分亚群")
save(pbmc,file = "pbmc.rdata")
table(pbmc@meta.data$seurat_clusters)
top_anno = HeatmapAnnotation(df = df,#细胞名/cluster
                             border = F,
                             show_annotation_name = T,
                             gp = gpar(col = 'black'),
                             col = list(clusters = c('sw-MBC'="#FF34B3",
                                                     'Naïve B cell'="#F6F5B4",
                                                     'Plasma'="#20B2AA",
                                                     'unsw-MBC'="#FFA500")))#颜色设置
table(pbmc@meta.data$Group)
top_anno = HeatmapAnnotation(df = df,#细胞名/cluster
                             border = F,
                             show_annotation_name = T,
                             gp = gpar(col = 'black'),
                             col = list(clusters = c('PBS'="#FF34B3",
                                                     'TD1'="#F6F5B4")))#颜色设置
cell_type_cols <- c(brewer.pal(9, "Set1"), "#FF34B3","#BC8F8F","#20B2AA",
                    "#00F5FF","#FFA500","#ADFF2F","#FF6A6A","#7FFFD4", 
                    "#AB82FF","#90EE90","#00CD00","#008B8B","#6495ED",
                    "#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00","#000000","#708090")  
#数据标准化缩放一下
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
Heatmap(marker_exp,
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = T,
        column_title = NULL,
        heatmap_legend_param = list(
          title=' '), 
        col = colorRampPalette(c("#00008B","white","#B22222"))(100),
        border = 'black',
        rect_gp = gpar(col = "black", lwd = 1),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        top_annotation = top_anno)
FeaturePlot(scRNA2,
            features = c("Alp"),
            dims = c(1, 2),
            cells = NULL,
            pt.size = NULL,
            order = FALSE,
            min.cutoff = 0,
            max.cutoff = 8,
            reduction = NULL,
            split.by = NULL,
            shape.by = NULL,
            slot = "scale.data",
            # blend = FALSE,
            blend.threshold = 0.5,
            label = T,
            label.size = 3,
            repel = FALSE,
            # ncol = 2,
            coord.fixed = FALSE,
            by.col = TRUE,
            sort.cell = NULL,
            interactive = FALSE,
            combine = TRUE)

rm(scRNA456)
scRNA2 <- scRNA456
save(scRNA2,file = "scRNA2.rdata")
getwd()



#################B细胞分亚群###########################
#################harmony矫正批次########################
library(devtools)
# install_github("immunogenomics/harmony")
library(Seurat)
library(cowplot)
library(harmony)
scRNAB <- subset(scRNA2,cell_type%in%"B cell")
# load('data/pbmc_stim.RData') #加载矩阵数据
#在运行Harmony之前，创建一个Seurat对象并按照标准PCA进行分析。
pbmc <- scRNAB %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = scRNAB@var.genes, npcs = 20, verbose = FALSE) #R语言中%>%的含义是什么呢，管道函数啦，就是把左件的值发送给右件的表达式，并作为右件表达式函数的第一个参数。
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = pbmc, reduction = "pca", pt.size = .1, group.by = "Group")
p2 <- VlnPlot(object = pbmc, features = "PC_1", group.by = "Group", pt.size = .1)
plot_grid(p1,p2)
options(repr.plot.height = 2.5, repr.plot.width = 6)
pbmc <- pbmc %>%
  RunHarmony("Group", plot_convergence = TRUE) #Harmony converged after 8 iterations
pbmc@reductions$harmony
harmony_embeddings <- Embeddings(pbmc, 'harmony')
harmony_embeddings[1:5, 1:5]
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = pbmc, reduction = "harmony", pt.size = .1, group.by = "Group")
p2 <- VlnPlot(object = pbmc, features = "harmony_1", group.by = "Group", pt.size = .1)
plot_grid(p1,p2)
pbmc <- pbmc %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.05) %>%
  identity()
options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(pbmc, reduction = "umap", group.by = "Group", pt.size = .1, split.by = 'Group')
pbmc=RunTSNE(pbmc,reduction = "harmony", dims = 1:20)
TSNEPlot(object = pbmc, pt.size = 0.5, label = TRUE,split.by='Group') 

DimPlot(pbmc, reduction = "umap",pt.size = .1,  label = TRUE)
TSNEPlot(pbmc, pt.size = .1, label = TRUE) 


getwd()

######################B细胞注释#######################
rm(list = ls())
table(pbmc@meta.data$cell_type)

Cell_type <- c("0" = "B cell",
               "1" = "T cell",
               "2" = "Monocyte",
               "3" = "T cell",
               "4" = "T cell",
               "5" = "Malignant cell",
               "6" = "NK cell",
               "7" = "Neutrophil",
               "8" = "Monocyte",
               "9" = "NKT cell",
               "10"= "Malignant cell",
               "11"= "DC cell",                                                                                                                        
               "12"= "DC cell",
               "13"= "Macrophage",
               "14"= "DC cell",
               "15"= "B cell",
               "16"= "Monocyte",
               "17"= "Macrophage",
               "18"= "Mast cell",
               "19"= "Macrophage",
               "20"= "Neutrophil",
               "21"= "B cell",
               "22"= "Erythroid cell")

scRNA2[['cell_type']] <- unname(Cell_type[scRNA2@meta.data$seurat_clusters])
library(ggplot2)
library(RColorBrewer)
cell_type_cols <- c(brewer.pal(9, "Set1"), "#FF34B3","#BC8F8F","#20B2AA",
                    "#00F5FF","#FFA500","#ADFF2F","#FF6A6A","#7FFFD4", 
                    "#AB82FF","#90EE90","#00CD00","#008B8B","#6495ED",
                    "#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00","#000000","#708090")  
# scRNA_umap1 <- subset(scRNA_umap,Group != c("PBMC of healthy donor"))
p1 <- DimPlot(scRNA2, 
              label = T,
              reduction = "umap",
              group.by = "cell_type",
              # split.by = "Group",
              cols= cell_type_cols, #设置颜色
              pt.size = 0.5,#设置点的大小
              # ncol = 3,
              # label.size = 2,
              repel = T,
              raster = FALSE)+labs(x = "UMAP1", y = "UMAP2")+ 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))  #加边框
# dev.off()

pbmc

p2 <- DimPlot(scRNA2, 
              label = T,
              reduction = "tsne",
              group.by = "cell_type",
              # split.by = "Group",
              cols= cell_type_cols, #设置颜色
              pt.size = 0.5,#设置点的大小
              # ncol = 3,
              # label.size = 2,
              repel = T,
              raster = FALSE)+labs(x = "tSNE1", y = "tSNE2")+ 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p2
p <- p1+p2
p
dev.off()



FeaturePlot(scRNA2,
            features = c("Cd3g"),
            dims = c(1, 2),
            cells = NULL,
            pt.size = NULL,
            order = FALSE,
            min.cutoff = 0,
            max.cutoff = 8,
            reduction = NULL,
            split.by = NULL,
            shape.by = NULL,
            slot = "scale.data",
            # blend = FALSE,
            blend.threshold = 0.5,
            label = T,
            label.size = 3,
            repel = FALSE,
            # ncol = 2,
            coord.fixed = FALSE,
            by.col = TRUE,
            sort.cell = NULL,
            interactive = FALSE,
            combine = TRUE)







################T细胞分亚群#################################
rm(list = ls())
library(devtools)
# install_github("immunogenomics/harmony")
library(Seurat)
library(cowplot)
library(harmony)
scRNAT <- subset(scRNA2,cell_type%in%"T cell")
getwd()
setwd("/home/data/t030312/LY_work/LQ/scrna/CD8")
# load('data/pbmc_stim.RData') #加载矩阵数据
#在运行Harmony之前，创建一个Seurat对象并按照标准PCA进行分析。
scRNAT_harmony <- scRNAT %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = scRNAT@var.genes, npcs = 20, verbose = FALSE) #R语言中%>%的含义是什么呢，管道函数啦，就是把左件的值发送给右件的表达式，并作为右件表达式函数的第一个参数。
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = scRNAT_harmony, reduction = "pca", pt.size = .1, group.by = "Group")
p2 <- VlnPlot(object = scRNAT_harmony, features = "PC_1", group.by = "Group", pt.size = .1)
plot_grid(p1,p2)
options(repr.plot.height = 2.5, repr.plot.width = 6)
scRNAT_harmony <- scRNAT_harmony %>%
  RunHarmony("Group", plot_convergence = TRUE) #Harmony converged after 8 iterations
scRNAT_harmony@reductions$harmony
harmony_embeddings <- Embeddings(scRNAT_harmony, 'harmony')
harmony_embeddings[1:5, 1:5]
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = scRNAT_harmony, reduction = "harmony", pt.size = .1, group.by = "Group")
p2 <- VlnPlot(object = scRNAT_harmony, features = "harmony_1", group.by = "Group", pt.size = .1)
plot_grid(p1,p2)
scRNAT_harmony1 <- scRNAT_harmony %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.1) %>%
  identity()
options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(scRNAT_harmony1, reduction = "umap", group.by = "Group", pt.size = .1, split.by = 'Group')
scRNAT_harmony1=RunTSNE(scRNAT_harmony1,reduction = "harmony", dims = 1:20)
TSNEPlot(object = scRNAT_harmony1, pt.size = 0.5, label = TRUE,split.by='Group') 

DimPlot(scRNAT_harmony1, reduction = "umap",pt.size = .1,  label = TRUE)
TSNEPlot(scRNAT_harmony1, pt.size = .1, label = TRUE) 


getwd()
save(scRNAT_harmony1,file = "scRNAT.rdata")




#########################提取TD1表达量##########################
getwd()
rm(list = ls())
gene_list <- c("Ccr7",   
               "Lef1",   
               "Sell",   
               "Myc",   
               "S1pr1",   
               "Tcf7",   
               "Ifng",   
               "Bcl2",   
               "Cd44",   
               "Cxcr3",   
               "Ly6e",   
               "Pdcd1",              
               "Gzmk",   
               "Gzma",   
               "Klre1",   
               "Gzmb",   
               "Cx3cr1",   
               "Zeb2",   
               "Ctla4",   
               "Ly6a",   
               "Il7r",   
               "Cxcr4")

#计算平均表达量
scRNAB <- scRNAT_harmony1
gene_cell_exp <- AverageExpression(scRNAB,
                                   features = gene_list,
                                   group.by = 'seurat_clusters',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)

#complexheatmap作图
library(ComplexHeatmap)
#顶部细胞类型注释
df <- data.frame(colnames(gene_cell_exp))
colnames(df) <- 'clusters'
top_anno = HeatmapAnnotation(df = df,#细胞名/cluster
                             border = F,
                             show_annotation_name = T,
                             gp = gpar(col = 'black'),
                             col = list(clusters = c('0'="#FF34B3",
                                                     '1'="#F6F5B4",
                                                     '2'="#20B2AA",
                                                     '3'="#FFA500",
                                                     '4'="#ADFF2F",
                                                     '5'="#FF6A6A",
                                                     '6'="#9ECABE",
                                                     '7'="#8B008B",
                                                     '8'="#FFA500")))
table(pbmc@meta.data$seurat_clusters)
# top_anno = HeatmapAnnotation(df = df,#细胞名/cluster
#                              border = F,
#                              show_annotation_name = T,
#                              gp = gpar(col = 'black'),
#                              col = list(clusters = c('sw-MBC'="#FF34B3",
#                                                      'Naïve B cell'="#F6F5B4",
#                                                      'Plasma'="#20B2AA",
#                                                      'unsw-MBC'="#FFA500")))#颜色设置
cell_type_cols <- c(brewer.pal(9, "Set1"), "#FF34B3","#BC8F8F","#20B2AA",
                    "#00F5FF","#FFA500","#ADFF2F","#FF6A6A","#7FFFD4", 
                    "#AB82FF","#90EE90","#00CD00","#008B8B","#6495ED",
                    "#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00","#000000","#708090")  
#数据标准化缩放一下
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
Heatmap(marker_exp,
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = T,
        column_title = NULL,
        heatmap_legend_param = list(
          title=' '), 
        col = colorRampPalette(c("#00008B","white","#B22222"))(100),
        border = 'black',
        rect_gp = gpar(col = "black", lwd = 1),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        top_annotation = top_anno)
FeaturePlot(pbmc,
            features = c("Jchain"),
            dims = c(1, 2),
            cells = NULL,
            pt.size = NULL,
            order = FALSE,
            min.cutoff = 0,
            max.cutoff = 8,
            reduction = NULL,
            split.by = NULL,
            shape.by = NULL,
            slot = "scale.data",
            # blend = FALSE,
            blend.threshold = 0.5,
            label = T,
            label.size = 3,
            repel = FALSE,
            # ncol = 2,
            coord.fixed = FALSE,
            by.col = TRUE,
            sort.cell = NULL,
            interactive = FALSE,
            combine = TRUE)
table(scRNAT_harmony1@meta.data$seurat_clusters)
VlnPlot(pbmc,
        features = c("Jchain"),
        split.by = "Group")
Cell_type <- c("0" = "Naïve T cell",
               "1" = "Cytotoxic T cell",
               "2" = "Exhausted T cell",
               "3" = "Exhausted T cell",
               "4" = "Ctla4+ T cell",
               "5" = "Memory T cell")
scRNAT_harmony1[['cell_type']] <- unname(Cell_type[scRNAT_harmony1@meta.data$seurat_clusters])
library(ggplot2)
library(RColorBrewer)
cell_type_cols <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF",
                    "#999999", "#FF34B3","#BC8F8F","#20B2AA",
                    "#00F5FF","#FFA500","#ADFF2F","#FF6A6A","#7FFFD4", 
                    "#AB82FF","#90EE90","#00CD00","#008B8B","#6495ED",
                    "#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00","#000000","#708090") 
cell_type_cols <- c("#377EB8","#4DAF4A","#FF7F00","#E41A1C","#FFFF33","#A65628","#F781BF",
                    "#999999", "#FF34B3","#BC8F8F","#20B2AA",
                    "#00F5FF","#FFA500","#ADFF2F","#FF6A6A","#7FFFD4", 
                    "#AB82FF","#90EE90","#00CD00","#008B8B","#6495ED",
                    "#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00","#000000","#708090") 
# scRNA_umap1 <- subset(scRNA_umap,Group != c("PBMC of healthy donor"))
# brewer.pal(9, "Set1")
p1 <- DimPlot(scRNAT_harmony1, 
              label = T,
              reduction = "tsne",
              group.by = "cell_type",
              # split.by = "Group",
              cols= cell_type_cols, #设置颜色
              pt.size = 0.5,#设置点的大小
              # ncol = 3,
              # label.size = 2,
              repel = T,
              raster = FALSE)+labs(x = "UMAP1", y = "UMAP2")+ 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))  #加边框
# dev.off()

p1 

p2 <- DimPlot(scRNAT_harmony1, 
              label = T,
              reduction = "tsne",
              group.by = "cell_type",
              # split.by = "Group",
              cols= cell_type_cols, #设置颜色
              pt.size = 0.5,#设置点的大小
              # ncol = 3,
              # label.size = 2,
              repel = T,
              raster = FALSE)+labs(x = "tSNE1", y = "tSNE2")+ 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p2
p <- p1+p2
p
dev.off()
table(pbmc@meta.data$cell_type)
Idents(pbmc) <- "Group"
cluster.markers <- FindMarkers(pbmc,
                               ident.1 = "PBS",
                               ident.2 = "TD1")
getwd()
setwd("/home/data/t030312/LY_work/LQ/scrna/B细胞分亚群")
write.csv(cluster.markers,file = "deg.csv")



rm(list = ls())
# install.packages("devtools")
devtools::install_github("Simon-Leonard/FlexDotPlot")
1
library(FlexDotPlot)
Idents(scRNAT_harmony1) <- "Group"
features <- c("Cd5","Ccr2","Il7r","Xcl1","Ccl5")
dp=DotPlot(scRNAT_harmony1, features = features) + RotatedAxis()
dp


library(FlexDotPlot)
# Dot plot with shape type (and not size) controlled by "Percent Expressed" parameter 
dot_plot(dp$data[,c(3,4,1,2,5)], size_var = "pct.exp", col_var = "avg.exp.scaled", 
         size_legend = "Percent Expressed", col_legend = "Average Expression",
         x.lab.pos = "bottom", display_max_sizes = FALSE)




VlnPlot(scRNAT_harmony1,
        features = c("Cd5","Ccr2","Il7r","Xcl1","Ccl5"),
        pt.size = 0.1)




VlnPlot(pbmc,
        features = c("Ddit3","Egr1","Fcer1g","Flt3",
        "Ifi204","Nr4a3","Klrk1","Rgcc","Gclc","Eaf2"),
        pt.size = 0.1)

Idents(pbmc) <- "cell_type"
VlnPlot(pbmc,
        features = c("Slamf1"),
        pt.size = 0.5)
rm(list = ls())





