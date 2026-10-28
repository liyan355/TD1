.libPaths(c("~/SeuratV4",.libPaths()))
library(Seurat)
packageVersion("Seurat")


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
scRNA2 <- merge(scRNAlist[[1]], y=c(scRNAlist[[2]]),
                add.cell.ids = names(scRNAlist))
dim(scRNA2)   #查看基因数和细胞总数
table(scRNA2@meta.data$orig.ident)  #查看每个样本的细胞数
save(scRNA2,file = 'scRNA42.merge.Rdata')

#未分组合并方法
############使用fastMNN整合数据#######################################
scRNAlist <- lapply(scRNAlist,FUN = function(x)NormalizeData(x))
scRNAlist <- lapply(scRNAlist,FUN = function(x)FindVariableFeatures(x))
scRNA <- RunFastMNN(object.list = scRNAlist)
scRNA <- RunUMAP(scRNA,reduction = "mnn",dims = 1:30)
scRNA <- FindNeighbors(scRNA,reduction = "mnn",dims = 1:30)
scRNA <- FindClusters(scRNA)
p1 <- DimPlot(scRNA, group.by = "orig.ident", pt.size=0.1) + 
  ggtitle("Integrated by fastMNN")
p1
p2 <- DimPlot(scRNA, group.by="orig.ident", pt.size=0.1) + 
  ggtitle("No integrated")
p2
p = p1 + p2 + plot_layout(guides='collect')
ggsave('fastMNN.png', p, width=8, height=4)




########Seurat对象的处理########

#找表达矩阵
scRNA2_counts <- scRNA2@assays[["RNA"]]@counts

#基因信息表就是矩阵的行名

#细胞信息表
Cell_information <- scRNA2@meta.data
View(scRNA@meta.data)

sparse.size <- object.size(x = scRNA2)
sparse.size

scRNA2.orig <- scRNA2
save(scRNA2.orig, file = "scRNA2.orig.Rda")


scRNA2@meta.data
library(stringr)

a <- scRNA2@meta.data %>%
  mutate(Group = case_when(str_sub(scRNA2@meta.data[,1],1,2) == "KO" ~ 'KO',
                           str_sub(scRNA2@meta.data[,1],1,2) == "WT" ~ 'WT'))
a <- scRNA@meta.data %>%
  mutate(Group = if_else(str_sub(scRNA2@meta.data[,1], 3, 5) == "AF",
                         "Control",
                         "Tumor"))
scRNA2@meta.data <- a
scRNA2.orig <- scRNA2
save(scRNA2.orig, file = "scRNA2.orig.Rda")
save(scRNA2, file = "scRNA2.Rda")

View(scRNA2@meta.data)
library(dplyr)
library(tidyverse)
library(stringr)
#######单细胞数据预处理#######

#去除线粒体基因和极端值
load("scRNA2.Rda")

#当线粒体基因含量很⾼，⽽这时候细胞的转录本数⽬⼜⽐较少的话，
#我们就认为这个细胞可能发⽣了死亡，我们就要把这个细胞过滤掉

#找到线粒体基因
??str_subset
mito.genes <- str_subset(string = rownames(scRNA2),
                         pattern = "MT-")
mito.genes


PercentageFeatureSet()  #可以用来富集分析

#用这种方法，可以求感兴趣的基因集，在每个细胞中的比例
scRNA2[["percent.mt"]] <- PercentageFeatureSet(scRNA2, 
                                               pattern = "^MT-")   #在细胞信息部分添加了线粒体比例

head(scRNA2[["percent.mt"]])


library(farver)
VlnPlot(object = scRNA2, 
        features  = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
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




scRNA2

a1 <- scRNA2@meta.data[scRNA2@meta.data[["percent.mt"]]>25,]

scRNA3 <- subset(scRNA2,
                 subset = nFeature_RNA > 200 &
                   nFeature_RNA < 5000 &
                   percent.mt > -Inf & # 极端值
                   percent.mt < 25)

######数据的质控2: 消除细胞周期##########

#细胞周期打分
scRNA3 <- CellCycleScoring(scRNA2, 
                           s.features = cc.genes$s.genes, 
                           g2m.features = cc.genes$g2m.genes)
head(scRNA3@meta.data[6:8, 6:8])
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
load("scRNA.Rdata")






scRNA3 <- scRNA2

#######数据的批次矫正,方法一 anchor########


##先进行数据标准化,数据标准化和识别高变基因需要在不同的批次下面完成
##去除批间差需要合并后完成
#####多个Seurat对象的处理######
scRNA3@meta.data[["orig.ident"]]

scRNA4 <- SplitObject(scRNA2,
                      split.by = "orig.ident")
scRNA4

for(i in 1:length(scRNA4)){
  scRNA4_all <- scRNA4[[i]]
  # Log标准化
  scRNA4_all <- NormalizeData(scRNA4_all,
                              normalization.method = "LogNormalize",
                              scale.factor = 100000)
  # 计算⾼度可变基因
  scRNA4_all <- FindVariableFeatures(scRNA4_all, 
                                     selection.method = "vst", 
                                     nfeatures = 1500) 
  scRNA4[[i]] <- scRNA4_all    # 表示将处理好的Seurat对象⼜重新存回scRNA4中
}

####开始综合合并数据#####
#第一步是寻找锚定对象
# 多组数据整合
# 寻找锚定对象
AnchorSet <- FindIntegrationAnchors(object.list = scRNA4, 
                                    normalization.method = "LogNormalize") 
AnchorSet

save(AnchorSet,file = "AnchorSet1500.Rdata")
scRNA5<- IntegrateData(anchorset = AnchorSet)  #根据锚定对象将数据整合到一起
DefaultAssay(scRNA5) <- "integrated"

save(scRNA5,file = "scRNA5.Rdata")

# scale去除批间差
scRNA5 <- ScaleData(scRNA5, 
                    features = rownames(scRNA4))

scRNA5_Anchor <- scRNA5
save(scRNA5,file = "scRNA5.rdata")



#######数据的批次矫正,方法二 log########

rm(list = ls())
save("scRNA4.Rdata")
scRNA4 <- SplitObject(use,
                      split.by = "orig.ident")
view(scRNA@meta.data)

scRNA4_Cont <- scRNA4 [[1]]
scRNA4_Cont <- SCTransform(scRNA4_Cont, 
                           variable.features.n = 1000,
                           vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
                           verbose = FALSE)
?logransform
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
                             variable.features.n = 1000,
                             #vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
                             verbose = FALSE)
}
#数据整合
#选择用于整合的基因
features <- SelectIntegrationFeatures(object.list = scRNA4,
                                      nfeatures = 1000)
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

scRNA4_log <- scRNA4_all
save(scRNA4_log,file = "scRNA4_log.Rdata")


#####################降维和聚类########################
library(Seurat)
library(tidyverse)
library(reshape2)
library(SingleR)
library(celldex)
library(BiocParallel)
library(Seurat)
library(tidyverse)
library(stringr)
library(patchwork)
# rm(list = ls())
# load("scRNA6_log_finCl.Rdata")
# dim(scRNA6_log)
# load("scRNA41_log.Rdata")
# 1+1

#####################PCA聚类########################
data <- scRNA5
data <- FindVariableFeatures(scRNA,selection.method = "vst", nfeatures = 1500)
data <- ScaleData(data)
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


###############tSNE分析##########################
scRNA <- RunTSNE(scRNA,dims = 1:19)
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
p33 <- DimPlot(object = scRNA6_tsne, 
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

p4 <- DimPlot(object = scRNA6_tsne, 
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
save(scRNA,file = "scRNA.Rdata")
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

p55 <- DimPlot(object = scRNA_umap, 
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
p6 <- DimPlot(object = scRNA_umap, 
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

save(scRNA,file = "scRNA_fanjia.rdata")
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
DefaultAssay(scRNA_umap) <- "RNA"

#####################添加 cell type###############################
Cell_type <- c("0" = "NK",
               "1" = "Macrophage",
               "2" = "Neutrophil",
               "3" = "Monocyte",
               "4" = "Neutrophil",
               "5" = "CD8+ T",
               "6" = "pDC",
               "7" = "CD8+ T",
               "8" = "Macrophage",
               "9" = "Neutrophil",
               "10"= "Macrophage",
               "11"= "B",                                                                                                                        
               "12"= "pDC",
               "13"= "CD4+ T",
               "14"= "Macrophage",
               "15"= "DC",
               "16"= "cDC",
               "17"= "CD8+ T",
               "18"= "Malignant cell",
               "19"= "Erythrocyte",
               "20"= "Mast cell")
scRNA[['cell_type']] <- unname(Cell_type[scRNA@meta.data$seurat_clusters])
save(scRNA,file = "scRNA.rdata")
dev.off()
dev.new()
DimPlot(scRNA, 
        reduction = "umap", 
        group.by = "cell_type",
        label = TRUE, 
        pt.size = 0.2) 
cnv_Tu$prediction
VlnPlot(scRNA11_log, 
        features = c("NGF","G-CSF","NPY","CX3CL","CXCL12","CCL2","CCL5"),
        pt.size = 0,
        ncol = 2,
        group.by = "cell_type")


###################改变active.ident##############

colnames(use@meta.data)
temp <- factor(use@meta.data$cell_type)
names(temp) <- rownames(use@meta.data)
use@active.ident <- temp
use@active.ident
use@meta.data

table(use@active.ident)
is.na(use@active.ident)
table(is.na(use@active.ident))
library(magrittr)
#############寻找每个群的Marker Gene##################
library(magrittr)
#找每个cluster的marker
scRNA15_log <- PreplogFindMarkers(scRNA_g4, assay = "log", verbose = TRUE)
markers <- FindAllMarkers(scRNA,
                          only.pos = TRUE,
                          min.pct = 0.1,
                          logfc.threshold = 0.25)

markers %>% group_by(cluster)%>%top_n(n=2,wt=avg_log2FC)
top10 <- markers%>%
  group_by(cluster)%>%
  top_n(n=10,wt=avg_log2FC)
write.csv(top10,file = "markers_top10.csv")
head(top10)
write.csv(markers,file = "markers.csv")
#画热图
library(farver)
use <- FindVariableFeatures(use)
use <- ScaleData(use)
DoHeatmap(use,features = top10$gene)+NoLegend()
tp <- scRNA13_log@meta.data
scRNA13_log@active.ident
save(scRNA11_log,file = "scRNA11_log.Rdata")
t <- c("a","b","p","MSC","tumor")
t[i]
a <- FindMarkers(scRNA13_log,
                 ident.1 = 11,
                 #ident.2 = 20,
                 min.pct = 0.5, # 过滤掉那些在50%以下细胞中检测到的基因
                 logfc.threshold = log(1))
view(a)
write.csv(a,file = "cluster11.csv") 




df1 <- Tu_cells.markers
df2 <- NTu_cells.markers
df <- filter(df1,!df1$gene%in%df2$gene)
write.csv(df,file = "fib中肿瘤组和非肿瘤组差异基因.csv")



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
hpca.se <- HumanPrimaryCellAtlasData
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
meta=scRNA13_log@meta.data 
head(meta)
plot1 <- DimPlot(scRNA13_log, reduction = "umap", label = TRUE)
plot2<-DimPlot(scRNA11_log, reduction = "tsne",
               label = TRUE)
plot1 + plot2
load("HumanPrimaryCellAtlasData.Rda")
#进行singleR注释
pbmc_for_SingleR <- GetAssayData(scRNA13_log, slot="data") ##获取标准化矩阵
pbmc.hesc <- SingleR(test = pbmc_for_SingleR, ref = hpca.se, labels = hpca.se$label.main) #
pbmc.hesc
table(pbmc.hesc$labels,meta$seurat_clusters)
#seurat 和 singleR的table表
table(pbmc.hesc$labels,meta$seurat_clusters)


scRNA13_log@meta.data$labels <-pbmc.hesc$labels
print(DimPlot(scRNA13_log, group.by = c("seurat_clusters", "labels"),reduction = "umap"))
#基于scores within cells
print(plotScoreHeatmap(pbmc.hesc))


#基于 per-cell “deltas”诊断
plotDeltaDistribution(pbmc.hesc, ncol = 3)

#与cluster结果比较
tab <- table(label = pbmc.hesc$labels,cluster = meta$seurat_clusters)
pheatmap(log10(tab + 10))




#使用多个数据库注释
#使用BP和HPCA两个数据库综合注释，使用list函数读入多个数据库
scRNA13_log <- pbmcpbmc3.hesc <- SingleR(test = pbmc_for_SingleR, ref = list(BP=bpe.se, HPCA=hpca.se), 
                                         labels = list(bpe.se$label.main, hpca.se$label.main)) 
table(pbmc3.hesc$labels,meta$seurat_clusters)

pbmc3@meta.data$labels <-pbmc3.hesc$labels

print(DimPlot(pbmc3, group.by = c("seurat_clusters", "labels"),reduction = "umap"))





#umap
singler <- SingleR(scRNA13_log@assays$RNA@data,
                   ref = ref,
                   labels = ref$label.main,
                   clusters = scRNA13_log@meta.data$seurat_clusters,
                   fine.tune = TRUE,
                   BPPARAM = SnowParam(8))
#tsne
singler <- SingleR(scRNA11_log@assays$RNA@data,
                   ref = ref,
                   labels = ref$label.main,
                   clusters = scRNA13_log@meta.data$seurat_clusters,
                   fine.tune = TRUE)
#BPPARAM = SnowParam(8))
# 查看细胞类型
singler$pruned.labels
#绘制带cell label的tsne和umap图
cell_ids <- singler$pruned.labels
names(cell_ids) <- levels(scRNA11_log)
levels(scRNA11_log)
scRNA11_log <- RenameIdents(scRNA11_log, cell_ids)
levels(scRNA11_log)
TSNEPlot(object = scRNA11_log, pt.size = 0.5, label = TRUE)  
UMAPPlot(object = scRNA13_log, pt.size = 0.5, label = TRUE)
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
df <- scRNA13_log@meta.data
View(df)
df2 <- cnv_Tu$prediction
View(df2)
df3 <- cnv_NC$prediction
View(df3)
library(readr)
library(dplyr)

df_Tu <- cnv_Tu$prediction
df_Tu
scRNA_T <- subset(fib_tsne2, Group == "Tumor")
dim(scRNA_T)  #35422
dim(scRNA_N)  #33286
table(duplicated(df_Tu$cell.names))
table(df_Tu$copykat.pred)
df_Tu <- df_Tu[!duplicated(df_Tu$cell.names),]
rownames(df_Tu) <- df_Tu$cell.names
scRNA_Tu <- AddMetaData(scRNA_T, metadata = df_Tu)
meta_Tu <- scRNA_Tu$copykat.pred
head(meta_Tu)
DimPlot(scRNA_Tu, group.by = "copykat.pred",pt.size = 1)
DimPlot(scRNA_Tu,reduction = "tsne",group.by="copykat.pred",
        pt.size = 1) + scale_color_manual(values = c("red", "gray","gray"))
table(scRNA_Tu@meta.data$copykat.pred)
view(scRNA_Tu@meta.data)
Tu_cells <- subset(scRNA_Tu,copykat.pred =="aneuploid")
NTu_cells <- subset(scRNA_Tu,copykat.pred =="diploid")

DimPlot(subset(scRNA_Tu, copykat.pred =="aneuploid"), 
        reduction = "tsne", # pca, umap, tsne
        group.by = "cell_type",
        label = T)


df_NC <- cnv_NC$prediction
df_NC 
dim(scRNA_T)  #35422
dim(scRNA_N)  #33286
table(duplicated(df_NC$cell.names))
table(df_NC$copykat.pred)
df_NC <- df_NC[!duplicated(df_NC$cell.names),]
rownames(df_NC) <- df_NC$cell.names
scRNA_NC <- AddMetaData(scRNA11_log, metadata = df_NC)
meta_NC <- scRNA_NC$copykat.pred
view(meta_NC)
DimPlot(scRNA_NC, group.by = "copykat.pred")
DimPlot(scRNA_NC,reduction = "tsne",group.by="copykat.pred",
        pt.size = 1) + scale_color_manual(values = c("red","gray","gray","gray"))
table(scRNA_NC@meta.data$copykat.pred)
view(scRNA_NC@meta.data)
Tu_cells <- subset(scRNA_NC,copykat.pred =="aneuploid")

DimPlot(subset(scRNA_NC, copykat.pred =="aneuploid"), 
        reduction = "tsne", # pca, umap, tsne
        group.by = "cell_type",
        label = T)



# SCLC_copykat_prediction.txt是上一步生成的结果
mallignant <- read.delim("SCLC_copykat_prediction.txt", row.names = 1)


# 原始count中有889个细胞，malignant只有824个，有65个细胞被ngene.chr=5过滤掉了
# 把细胞的良恶性信息加入metadata
scRNA <- AddMetaData(scRNA11_log, metadata = mallignant)
p1 <- DimPlot(scRNA11_log, group.by = "cell_type", label = T) + ggsci::scale_color_lancet()
p2 <- DimPlot(scRNA11_log, group.by = "copykat.pred") + scale_color_manual(values = c("red", "gray"))
p1 
p2
ggsave("pred_mallignant.pdf", pc, width = 12, height = 5)
scRNA$celltype2 <- scRNA$celltype
scRNA$celltype2[scRNA$copykat.pred=="aneuploid"] <- "mallignant_cell"
table(scRNA$celltype2)
Tumor_cells <- subset(cnv_Tu$prediction,copykat.pred=="aneuploid")
head(cnv_Tu)
head(Tumor_cells)
dim(Tumor_cells)
###################改变active.ident##############

colnames(scRNA_Tu@meta.data)
temp <- factor(scRNA_Tu@meta.data$copykat.pred)
names(temp) <- rownames(scRNA_Tu@meta.data)
scRNA_Tu@active.ident <- temp
scRNA_Tu@active.ident
scRNA_Tu@meta.data

table(scRNA11_log@active.ident)
is.na(scRNA11_log@active.ident)
table(is.na(scRNA11_log@active.ident))



###################GO和KEGG分析##############################
#BP,CC和MF三种通路一起富集
library(Seurat)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
ego_ALL <- enrichGO(gene)

BiocManager::install("copykat")
library(devtools)
install.packages("devtools")
install_github("navinlabcode/copykat")
library(Seurat)
library(copykat)
library(tidyverse)
load("~/Project/SingleCell/Company_data/mt_scRNA3500_11_log.Rdata")
mt_scRNA3500_11_Tu_0261619 <- subset(mt_scRNA3500_11_log, seurat_clusters %in% c("0","2","6","16","19"))
dim(mt_scRNA3500_11_Tu_0261619@assays$RNA@counts)



#####################添加 cell type###############################
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
scRNA11_log[['cell_type']] <- unname(Cell_type[scRNA11_log@meta.data$seurat_clusters])

DimPlot(scRNA11_log, 
        reduction = "tsne", 
        group.by = "cell_type",
        label = TRUE, 
        pt.size = 0.2) + NoLegend()
###############拟时序分析#####################
.libPaths("/home/data/t030312/R/x86_64-conda-linux-gnu-library/4.1")
library(Matrix)
library(VGAM)
library(splines)
library(DDRTree)
library(irlba)
library(monocle)
library(dplyr)
library(patchwork)
library(Seurat)
#monocle构建CDS需要3个矩阵：expr.matrix、pd、featuredata
load("scRNA11_log.Rda")
scRNA11_log@meta.data$cell_type
sub_cellgroup <- subset(scRNA11_log,cell_type=="Fibroblast Cell")
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
sc_cds <- sc_cds[fData(sc_cds)$num_cells_expressed > 10, ]
cds <- sc_cds
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
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters = 4) 
plot_cell_clusters(cds, 1, 2)
plot_cell_clusters(cds, 1, 2, color = "CellType",
                   markers =c("CSRP2","POSTN","COL1A1","ACTA2","COL3A1","EDNRA","LRRC15","CFD","MYC"))
pData(cds)$Cluster=pData(cds)$celltype
cds <- clusterCells(cds, num_clusters = 3)
plot_cell_clusters(cds, 1, 2, color = "CellType",
                   markers = c("CSRP2","POSTN","COL1A1","ACTA2","COL3A1","EDNRA","LRRC15","CFD","MYC"))
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
sub_cellgroup$samples
plot_cell_trajectory(cds, color_by = "samples")
plot_cell_trajectory(cds, color_by = "samples") +
  facet_wrap(~State, nrow = 1)
plot_cell_trajectory(cds, color_by = "samples") +
  facet_wrap(~samples, nrow = 1)
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
BEAM_res <- BEAM(cds, branch_point = 1, cores = 8)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
table(BEAM_res$qval < 1e-10)
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res, qval < 1e-20)),],
                            branch_point = 1,
                            num_clusters = 4,
                            use_gene_short_name = TRUE,
                            show_rownames = TRUE)
#####################添加 cell type为了cellchat###############################
Cell_type <- c("0" = "Fibroblast cells 1",
               "1" = "Endothelial cells 1",
               "2" = "Fibroblast cells 2",
               "3" = "Pericytes",
               "4" = "Muscle cells 1",
               "5" = "Muscle cells 2",
               "6" = "Fibroblast cells 3",
               "7" = "NK cells",
               "8" = "Endothelial cells 2",
               "9" = "T cells",
               "10"= "Fibroblast cells 4",
               "11"= "Muscle cells 3",                                                                                                                        
               "12"= "Muscle cells 4",
               "13"= "Neutrophils",
               "14"= "Mast cells",
               "15"= "Fibroblast cells 5",
               "16"= "Myofibroblasts",
               "17"= "Neuron",
               "18"= "Myeloid cells",
               "19"= "Fibroblast cells 6",
               "20"= "Endothelial cells 3",
               "21"= "Endothelial cells 4",
               "22"= "B cells")
scRNA11_log[['cell_type']] <- unname(Cell_type[scRNA11_log@meta.data$seurat_clusters])

DimPlot(scRNA11_log, 
        reduction = "tsne", 
        group.by = "cell_type",
        label = TRUE, 
        pt.size = 0.2) + NoLegend()
###################细胞通讯cellchat#################################
rm(list = ls())
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
save("scRNA11_log.Rdata")
load("scRNA11_log.Rdata")
scRNA11_log@meta.data
Tu_log <- subset(scRNA11_log,Group=="Tumor")
table(Tu_log@meta.data$cell_type)
Tu_log@meta.data$orig.ident
Tu_log@meta.data$cell_type
cellchat <- createCellChat(Tu_log@assays$RNA@data,
                           meta = Tu_log@meta.data,
                           group.by = "cell_type")
summary(cellchat)
str(cellchat)
levels(cellchat@idents)
cellchat <- setIdent(cellchat,ident.use = "cell_type")
groupSize <- as.numeric(table(cellchat@idents))
groupSize
#导入配体受体数据库
CellChatDB <- CellChatDB.human
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
#dev.new()下一步不出图的时候运行
showDatabaseCategory(CellChatDB)
unique(CellChatDB$interaction$annotation)#查看可以选择的侧面，也就是上图左中的三种
#选择secreated signaling进行后续细胞互作分析
CellChatDB.use <- subsetDB(CellChatDB,search = "ECM-Receptor")
cellchat@DB <- CellChatDB.use
#在矩阵的所有基因中提取signaling gene，结果保存在data.signaling 里面
cellchat <- subsetData(cellchat)
#future::plan("multisession",workers=4)#可以不设置这一步
#multiprocess
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat,PPI.human)
#根据表达值推测细胞互作的概率（cellphoneDB是用平均表达值代表互作强度）
cellchat <- computeCommunProb(cellchat,raw.use=F,population.size=T)
cellchat <- filterCommunication(cellchat,min.cells=10)
df.net <- subsetCommunication(cellchat)
#write.csv(df.net,"net_lr.csv")
#推断信号通路水平的细胞通讯网络
cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat,slot.name = "netP")
#write.csv(df.netp,"net_pathway_Pt.csv")



#至此，统计分析部分已经完成
#细胞互作关系展示
#所有细胞群总体观，细胞互作数量与强度统计分析
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow=c(1,1),xpd=TRUE)
netVisual_circle(cellchat@net$count,vertex.weight = groupSize,weight.scale = T,
                 label.edge = F,title.name = "Number of Interactions")
netVisual_circle(cellchat@net$weight,vertex.weight = groupSize,weight.scale = T,
                 label.edge = F,title.name = "Number of Interactions")
#检查每种细胞发出的信号
mat <- cellchat@net$count
par(mfrow =c(2,3),xpd=T)
for (i in 1:nrow(mat)){
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize, weight.scale = T, arrow.width = 1,
                   arrow.size = 0.3, edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
mat <- cellchat@net$weight
par(mfrow=c(2,3),xpd=T)
for (i in 1:nrow(mat)){
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize, weight.scale = T, arrow.width = 1,
                   arrow.size = 0.3, edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
cellchat@netP$pathways
pathway.show <- c("RELN")
#[1] "COLLAGEN" "FN1"      "LAMININ"  "THBS"     "TENASCIN" "HSPG"     "AGRN"     "VTN"     
#[9] "RELN"
levels(cellchat@idents)
#层次图
vertex.receiver=c(1,4)
netVisual_aggregate(cellchat,signaling = pathway.show,
                    vertex.receiver = vertex.receiver)
par(mfrow=c(1,1))
#网络图
netVisual_aggregate(cellchat,signaling = pathway.show,layout = "circle")
#和弦图
netVisual_aggregate(cellchat,signaling = pathway.show,layout = "chord")
#热图
netVisual_heatmap(cellchat,signaling = pathway.show,color.heatmap = "Reds")
#配体-受体层级的可视化
#计算配体受体对选定信号通路的贡献值（看哪条信号通路对某通路贡献大）
netAnalysis_contribution(cellchat,signaling = pathway.show)
pairR. <- extractEnrichedLR(cellchat,signaling = pathway.show,
                                geneLR.return = F)#提取对WNT有贡献的所有配体
#层次图
LR.show <- pairR.[1,]
vertex.receiver <- c(1,2,3,4)
netVisual_individual(cellchat,signaling = pathway.show,pairLR.use = LR.show,
                     vertex.receiver = verte)
netVisual_individual(cellchat,signaling = pathway.show,pairLR.use = LR.show,
                     layout = "circle")
#和弦图
netVisual_individual(cellchat,signaling = pathway.show,pairLR.use = LR.show,
                     layout = "chord")
#自动批量保存每个信号通路的作用结果
pathway.show.all <- cellchat@netP$pathways
levels(cellchat@idents)
vertex.receiver=c(1:23)
dir.create("ECM_NC(6-26)")
setwd("ECM_NC(6-26)")
for(i in 1:length(pathway.show.all)){
  netVisual(cellchat,signaling = pathway.show.all[i],out.format = c("pdf"),
            vertex.receiver = vertex.receiver,layout = "circle")#绘制网络图
  gg <- netAnalysis_contribution(cellchat,signaling = pathway.show.all[i])
  ggsave(filename = paste0(pathway.show.all[i],"_L-R_contribution.pdf"),
         plot = gg,width = 10,height = 10,units = 'in',dpi = 300)
}
setwd("D:\\r\\Singlecell")
#气泡图
library(ggplot2)
levels(cellchat@idents)
p <- netVisual_bubble(cellchat,sources.use = c(19),
                      targets.use = c(4,11,13,17),remove.isolate = F)
p
#气泡图，指定配体受体或者指定信号通路
p1 <- netVisual_bubble(cellchat,sources.use = c(19),targets.use = c(11,13,18),
                 signaling = "AGRN")
p1               
#气泡图（指定配体受体或信号通路并指定细胞）
pairLR.use <- extractEnrichedLR(cellchat,signaling = "AGRN")
p3 <- netVisual_bubble(cellchat,sources.use = c(19),targets.use = c(11,13,18),
                 pairLR.use = pairLR.use,remove.isolate = T)
p3
#参与某条信号通路（如TGFb）的所有基因在细胞群中的表达情况展示（小提琴图和气泡图）
p4 <- plotGeneExpression(cellchat,signaling = "AGRN")
p4

#社会网络分析（通讯网络中的角色识别）
cellchat <- netAnalysis_computeCentrality(cellchat,slot.name ="netP")
pathway.show <- c("VTN")
#"COLLAGEN","FN1","LAMININ","THBS",
#"TENASCIN","HSPG","AGRN","VTN","RELN"
p5 <- netAnalysis_signalingRole_network(cellchat,signaling = pathway.show,
                                  width = 15,height = 6,font.size = 10)
#识别细胞的信号流模式
ht1 <- netAnalysis_signalingRole_heatmap(cellchat,pattern = "outgoing",font.size = 5)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat,pattern = "incoming",font.size = 5)
ht1
ht2
#非负矩阵分解（NMF）识别细胞的通讯模式
#信号输出细胞的模式识别
library(NMF)
library(ggalluvial)
selectK(cellchat,pattern = "incoming")
nPatterns=2
cellchat <- identifyCommunicationPatterns(cellchat,pattern = "incoming",k=nPatterns,
                                          width = 5,height = 14,font.size = 8)
cellchat
#河流图
netAnalysis_river(cellchat,pattern = "incoming")
#气泡图
netAnalysis_dot(cellchat,pattern = "incoming")
#incoming,outgoing
#在二维空间中可视化主要的发送者（源）和接受者（目标）
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat,signaling = "COLLAGEN")
gg1+gg2
#信号网络的流行学习与分类
#基于功能相似性的流行学习与分类
install.packages("umap-learn")
cellchat <- computeNetSimilarity(cellchat,type = "functional")
cellchat <- netEmbedding(cellchat,type = "functional")
Y
cellchat <- netClustering(cellchat,type = "functional")


###################GO分析###################################
options(stringsAsFactors = F)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
df1 <- subset(fib,cluster=="0")
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
write.csv(ego_BP,file = "GO_BP_fib0.csv")
dotplot(ego_BP,showCategory = 20,color = "pvalue")
#CC细胞组分
ego_CC <- enrichGO(gene = genelist$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = T)
write.csv(ego_CC,file = "GO_BP_fib0.csv")
dotplot(ego_CC,showCategory = 20,color = "pvalue")
#分子功能
ego_MF <- enrichGO(gene = genelist$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = T)
write.csv(ego_MF,file = "GO_BP_fib0.csv")
barplot(ego_MF,showCategory = 20,color = "pvalue")

#####################细胞数目统计图#################################
table(fib_umap2$orig.ident)#查看各组细胞数
prop.table(table(Idents(fib_umap2)))
table(Idents(fib_umap2), fib_umap2$orig.ident)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(fib_umap2), fib_umap2$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = allcolour)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))


