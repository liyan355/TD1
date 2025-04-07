####################GSVA分析###########################
getwd()
rm(list = ls())
View(scRNAb@meta.data)
table(scRNAb@meta.data$cancer)

.libPaths(c("~/SeuratV4",.libPaths()))
library(Seurat)
packageVersion("Seurat")
library(dplyr)
# unloadNamespace("SummarizedExperiment")
library(GSVA)
library(GSEABase)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(limma)
library(tidyverse)
library(Seurat)
library(pheatmap)
# devtools::install_github("xjsun1221/tinyarray")
library(tinyarray)
patient_list <- c("PCall_HCC10_P","PCall_HCC38B_T","PCall_HCC19_P",
                  "PCall_HCC24_T","PCall_HCC19_T","PCall_HCC41B_T","PCall_HCC40B_T","PCall_HCC18_P",
                  "PCall_HCC39B_P","PCall_HCC22_P","PCall_HCC21_P","PCall_HCC22_T","PCall_HCC21_T","PCall_HCC23_T","PCall_HCC18_T",
                  "PCall_HCC31_T","PCall_HCC32_T","PCall_HCC28_T","PCall_HCC36B_T","PCall_HCC36B_P")
scRNAp <- subset(scRNAb,orig.ident%in%patient_list)
table(scRNAp@meta.data$orig.ident)
# genesets <- getGmt('./c2.cp.v2023.2.Hs.symbols.gmt')
# scRNAb_LIHC
# immuneT <- subset(immune, celltype=="T cells")#提取我们需要分析的细胞类型
# immuneT <- as.matrix(scRNAb_LIHC@assays$RNA@counts)#提取count矩阵
# View(scRNAb_LIHC@meta.data)
# meta <- scRNAb_LIHC@meta.data[,c("orig.ident","patient","type","celltype","TD1","Breg")]#分组信息，为了后续作图
# 
msigdbr_species()
scRNA_LIHC <- subset(scRNAb,cancer=="LIHC")
# human[1:4,1:4]
save(scRNA_LIHC,file = "scRNA_LIHC.rdata")

human <- msigdbr(species = "Homo sapiens")
table(human[,1])
table(human$gs_subcat)
# genesets <- msigdbr(species = "Homo sapiens", category = "C5") 
# genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
# genesets <- split(genesets$gene_symbol, genesets$gs_name)
Idents(scRNAb) <- "TD1"
# expr <- AverageExpression(scRNAb, assays = "RNA", slot = "data")[[1]]
# expr <- expr[rowSums(expr)>0,]  #选取非零基因
# expr <- as.matrix(expr)
# head(expr)
# setwd("/home/wucheng/jianshu/function/data")
pbmc <-scRNA_LIHC
expr <- as.data.frame(pbmc@assays$RNA@data) #表达矩阵
meta <- pbmc@meta.data[,c("orig.ident","celltype","TD1")] #类别
m_df = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") #选取物种人类
msigdbr_list = split(x = m_df$gene_symbol, f = m_df$gs_name)
expr=as.matrix(expr) 

kegg <- gsva(expr, 
             msigdbr_list, 
             kcdf="Poisson") #gsva
save(kegg,file = "kegg.rdata")
kegg_df <- read.csv("GSVA.csv",row.names = 1)
pheatmap(kegg_df, 
         show_rownames=T, 
         show_colnames=F, 
         annotation_col=meta,
         fontsize_row=5, 
         filename='gsva_heatmap.png', 
         width=15, 
         height=12)#绘制热图
pheatmap(kegg, 
         show_rownames=F, 
         show_colnames=T, 
         # annotation_col=meta,
         filename='gsva_heatmap1.png', 
         width=15, 
         height=12)#绘制热图







###############取平均表达量作GSVA#############################



rm(list = ls())
scRNA <- scRNA_LIHC
DefaultAssay(scRNA) <- "RNA"
scRNA <- NormalizeData(scRNA)

# 选择基因集
genesets <- msigdbr(species = "Homo sapiens", category = "C5") 
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)
#Idents设置按什么来取组的表达均值（计算case control之间的均值也可以）
Idents(scRNA) <- "orig.ident" 
expr1 <- AverageExpression(scRNA_LIHC, assays = "RNA", slot = "data")[[1]]
# expr <- expr[rowSums(expr)>0,]  #选取非零基因
expr1 <- as.matrix(expr1)
# gsva默认开启全部线程计算
gsva.res <- gsva(expr1, genesets, method="ssgsea") 
saveRDS(gsva.res, "gsva.res.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_res.csv", row.names = F)




#####################################Addmodulescore打分绘制箱线图######################
library(tidyverse)
library(Matrix)
library(cowplot)
pbmc <- readRDS("pbmc.rds") 
rm(list = ls())

getwd()
setwd("/home/data/t030312/LY_work/scRNA_B_cells/GSVA")

Idents(scRNA_LIHC) <- "TD1"


markers <- FindAllMarkers(scRNA_LIHC,
                          only.pos = TRUE,
                          min.pct = 0.1, # 过滤掉那些在25%以下细胞中检测到的基因
                          logfc.threshold = 0.25) 
TD1_gene <- subset(markers,cluster=="SLAMF1+")
View(TD1_gene)

TD1_gene <- TD1_gene$gene
cd_features <- list(TD1_gene)


# scRNA_LIHC <- AddModuleScore(
#   object = scRNA_LIHC,
#   features = cd_features,
#   ctrl = 100, #默认值是100
#   name = 'CD_Features'
# )

scRNA_LIHC <- AddModuleScore(
  object = scRNA_LIHC,
  features = cd_features,
  ctrl = 100, #默认值是100
  name = 'CD_Features'
)
colnames(scRNA_LIHC@meta.data)

colnames(scRNA_LIHC@meta.data)[18] <- 'CD150_AddmoduleScore' 
VlnPlot(scRNA_LIHC,
        raster=FALSE,
        group.by = "orig.ident",
        features = 'CD150_AddmoduleScore')
Idents(scRNA_LIHC) <- "orig.ident"

dev.off()
View(scRNA_LIHC@meta.data)

dev.off()
Idents(scRNA_LIHC) <- "orig.ident"
table(Idents(scRNA_LIHC))



pc1 <- FetchData(object = scRNA_LIHC, 
                  vars =  c('orig.ident', 'CD150_AddmoduleScore'))
head(x = pc1)


library(dplyr)

avg_sd <- pc1 %>%
  group_by(orig.ident) %>%
  summarize(
    avg_score = mean(CD150_AddmoduleScore),
    sd_score = sd(CD150_AddmoduleScore)
  )

avg_sd


###转换数据矩阵############
ave_sd <- as.data.frame(t(avg_sd))

colnames(ave_sd) <- ave_sd[1,]

ave_sd <- ave_sd[-1,]

ave_sd <- ave_sd[-2,]
rownames(ave_sd) <- "CD150_AddmoduleScore"

#############rbind#######################
new_df <- rbind(ave_sd,kegg_df)




GO_data <- readxl::read_xlsx("孙老师通路.xlsx",col_names = F)

colnames(GO_data) <- "pathway_list"

new_df$pathway_list <- rownames(new_df)



df <- filter(new_df,new_df$pathway_list%in%GO_data$pathway_list)
ave_sd$pathway_list <- rownames(ave_sd)
DF <- rbind(ave_sd,df)


save(new_df,file = "new_df.rdata")
save(ave_sd,file = "ave_sd.rdata")
save(GO_data,file = "GO_data.rdata")
save(DF,file = "DF.rdata")
DF <- DF[,-39]

DF1 <- DF




DF <- as.data.frame(lapply(DF,
                           function(x) as.numeric(as.character(x))))

rownames(DF) <- rownames(DF1)
pheatmap(DF, 
         show_rownames=T, 
         show_colnames=T, 
         cluster_rows = T,
         # annotation_col=meta,
         filename='gsva_heatmap1.png', 
         width=15, 
         height=12)#绘制热图
pheatmap(DF, 
         show_rownames=T, 
         show_colnames=T, 
         cluster_rows = F,
         # annotation_col=meta,
         filename='gsva_heatmap2.png', 
         width=15, 
         
         
         
         height=12)#绘制热图


######################至此GSVA分析完毕################################









dim(kegg)
dim(scRNA_LIHC@meta.data)
VlnPlot(scRNA_LIHC,
        features = 'SLAMF1',
        raster=FALSE,
        cols = my36colors,
        group.by = "cancer")+
  geom_boxplot(width=.2,col="black",fill="white")
library(RColorBrewer)


my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')

p1 <- ggboxplot(scRNA_LIHC@meta.data, x="cancer", y="CD150_AddmoduleScore", width = 0.6, 
                color = "black",#轮廓颜色
                fill="cancer",#填充
                palette = my36colors,
                xlab = F, #不显示x轴的标签
                bxp.errorbar=T,#显示误差条
                bxp.errorbar.width=0.5, #误差条大小
                size=1, #箱型图边线的粗细
                outlier.shape=NA, #不显示outlier
                legend = "right")+
  scale_y_continuous(limits = c(-0.3,1),breaks = seq(-0.5,1,0.2))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        # axis.line = element_line(colour = "#000000",size = 2),
        # axis.text = element_text(colour = "#000000" ,size = 6),
        axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75), ##就是这里
        # axis.ticks = element_line(colour = "#000000" ,size = 2) ,
        # axis.ticks.length = unit(2,'mm'),
        # plot.margin = unit(c(0.5,0,0,0),"cm"),
        # axis.title.y = element_text(size = 27),
        # axis.title.x = element_blank(),
        # plot.title = element_text(hjust = 0.5),
        legend.position = "right")
p1
View(scRNA_LIHC@meta.data)
View(scRNA_LIHC@meta.data)
DefaultAssay(scRNA_LIHC) <- "RNA"
write.csv(TD1_gene,file = "TD1_GENE.csv")
getwd()
setwd("/home/data/t030312/LY_work/scRNA_B_cells/7-26")
table(scRNA_LIHC_cancer@meta.data$orig.ident)
scRNA_LIHC_cancer <- subset(scRNA_LIHC,type%in%"Cancer")
library(ggpubr)
p2 <- ggboxplot(scRNA_LIHC_cancer@meta.data, x="cancer", y="CD150_AddmoduleScore", width = 0.6, 
                order = c("ESCA","UCEC","LIHC","CESC","LC",
                          "BRCA","CTCL","HNSC","GBC",
                          "OV","THCA","THYM","GIST","COAD","BLCA",
                          "STAD","PAAD","RCC","NB","CHOL"),
                color = "black",#轮廓颜色
                fill="cancer",#填充
                palette = my36colors,
                xlab = F, #不显示x轴的标签
                bxp.errorbar=T,#显示误差条
                bxp.errorbar.width=0.2, #误差条大小
                size=1, #箱型图边粗细
                outlier.shape=NA, #不显示outlier
                legend = "right")+
  scale_y_continuous(limits = c(-0.3,1.0),breaks = seq(-0.5,1,0.2))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        # axis.line = element_line(colour = "#000000",size = 2),
        # axis.text = element_text(colour = "#000000" ,size = 6),
        axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75), ##就是这里
        # axis.ticks = element_line(colour = "#000000" ,size = 2) ,
        # axis.ticks.length = unit(2,'mm'),
        # plot.margin = unit(c(0.5,0,0,0),"cm"),
        # axis.title.y = element_text(size = 27),
        # axis.title.x = element_blank(),
        # plot.title = element_text(hjust = 0.5),
        legend.position = "right")
p2
save(scRNA_LIHC,file = "scRANb.rdata")
save(scRNA_LIHC_cancer,file = "scRNA_LIHC_cancer.rdata")

library(ggpubr)
p3 <- ggboxplot(scRNA_LIHC@meta.data, x="celltype", y="CD150_AddmoduleScore", width = 0.6, 
                # order = c("ESCA","GBC","UCEC","HNSC","CTCL","CESC","LIHC",
                #           "OV","LC","THYM","COAD","BLCA",
                #           "BRCA","THCA","NB","STAD","RCC","CHOL","GIST","PAAD"),
                color = "black",#轮廓颜色
                fill="celltype",#填充
                palette = my36colors,
                xlab = F, #不显示x轴的标签
                bxp.errorbar=T,#显示误差条
                bxp.errorbar.width=0.5, #误差条大小
                size=1, #箱型图边线的粗细
                outlier.shape=NA, #不显示outlier
                legend = "right")+
  scale_y_continuous(limits = c(-0.3,1.0),breaks = seq(-0.5,1,0.2))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        # axis.line = element_line(colour = "#000000",size = 2),
        # axis.text = element_text(colour = "#000000" ,size = 6),
        axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75), ##就是这里
        # axis.ticks = element_line(colour = "#000000" ,size = 2) ,
        # axis.ticks.length = unit(2,'mm'),
        # plot.margin = unit(c(0.5,0,0,0),"cm"),
        # axis.title.y = element_text(size = 27),
        # axis.title.x = element_blank(),
        # plot.title = element_text(hjust = 0.5),
        legend.position = "right")
p3
table()
p4 <- ggboxplot(scRNA_LIHC_cancer@meta.data, x="celltype", y="CD150_AddmoduleScore", width = 0.6, 
                order = c("B.04.MT1X+B","B.14.Plasmablast",
                          "B.13.Cycling_GCB","B.12.LMO2+LZ_GCB",
                          "B.10.ENO1+Pre_GCB","B.09.DUSP4+AtM",
                          "B.07.CCR7+ACB3","B.11.SUGCT+DZ_GCB",
                          "B.03.HSP+B","B.05.EGR1+ACB",
                          "B.15.Plasma cell",
                          "B.02.IFIT3+B","B.08.ITGB1+SwBm","B.06.NR4A2+ACB2","B.01.TCL1A+naiveB"),
                color = "black",#轮廓颜色
                fill="celltype",#填充
                palette = my36colors,
                xlab = F, #不显示x轴的标签
                bxp.errorbar=T,#显示误差条
                bxp.errorbar.width=0.2, #误差条大小
                size=1, #箱型图边线的粗细
                outlier.shape=NA, #不显示outlier
                legend = "right")+
  scale_y_continuous(limits = c(-0.3,1.0),breaks = seq(-0.5,1,0.2))+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        # axis.line = element_line(colour = "#000000",size = 2),
        # axis.text = element_text(colour = "#000000" ,size = 6),
        axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75), ##就是这里
        # axis.ticks = element_line(colour = "#000000" ,size = 2) ,
        # axis.ticks.length = unit(2,'mm'),
        # plot.margin = unit(c(0.5,0,0,0),"cm"),
        # axis.title.y = element_text(size = 27),
        # axis.title.x = element_blank(),
        # plot.title = element_text(hjust = 0.5),
        legend.position = "right")
p4
View(scRNA_LIHC@meta.data)
table(scRNA_LIHC@meta.data$celltype)


df <- read.csv("Heatmap_new1.CSV",row.names = 1)
library(pheatmap)
df1 <- df[-1,]
pheatmap(df1, 
         show_rownames=T, 
         show_colnames=T, 
         cluster_rows = F,
         cluster_cols = F,
         color = colorRampPalette(c("#53A85F","white","#E95C59"))(100), #表示热图颜色,(100)表示100个等级
         # annotation_col=meta,
         filename='gsva_heatmap3.png', 
         fontsize_number = 10, #表示热图上显示数字的字体大小
         number_color = "black",
         display_numbers = T, 
         border_color = NA, #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         width=14, 
         height=8)#绘制热图
