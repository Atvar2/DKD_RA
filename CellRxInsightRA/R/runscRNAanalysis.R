# The codes were used for performing analysis in diabetic nephropathy project,
# and the codes mainly basing on package Seurat v4.

suppressPackageStartupMessages(library(sceasy))
suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(progress))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(plot1cell))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggsignif))
suppressPackageStartupMessages(library("scales"))
suppressPackageStartupMessages(library(CellChat))


# maxsize only impacts the maximum amount of data that is passed to the new future process, not how much RAM it can use
options(stringsAsFactors = FALSE, future.globals.maxSize = 20000 * 1024^2)

# ================1.quality contrl for 10 x data ===============================

CTL1 <- Read10X("./KC3/outs/filtered_feature_bc_matrix/")
CTL2 <- Read10X("./KC2/outs/filtered_feature_bc_matrix/")
DK1<- Read10X("./Q8/outs/filtered_feature_bc_matrix/")
DK2 <- Read10X("./Q14/outs/filtered_feature_bc_matrix/")
OR_1 <- Read10X("./Q66/outs/filtered_feature_bc_matrix/")
OR_2 <- Read10X("./Q68/outs/filtered_feature_bc_matrix/")

CTL1 <- CreateSeuratObject(counts = CTL1, project = "CTL1", min.cells = 3, min.features = 200)
CTL2 <- CreateSeuratObject(counts = CTL2, project = "CTL2", min.cells = 3, min.features = 200)
DK1 <- CreateSeuratObject(counts = DK1, project = "DK1", min.cells = 3, min.features = 200)
DK2 <- CreateSeuratObject(counts = DK2, project = "DK2", min.cells = 3, min.features = 200)
OR_1 <- CreateSeuratObject(counts = OR_1, project = "OR_1", min.cells = 3, min.features = 200)
OR_2 <- CreateSeuratObject(counts = OR_2, project = "OR_2", min.cells = 3, min.features = 200)


CTL1$Type <- "CTL"
CTL2$Type <- "CTL"
DK1$Type <- "DK"
DK2$Type <- "DK"
OR_1$Type <- "OR"
OR_2$Type <- "OR"

CTL1[["percent.mt"]] <- PercentageFeatureSet(CTL1, pattern = "^mt-")
CTL2[["percent.mt"]] <- PercentageFeatureSet(CTL2, pattern = "^mt-")
DK1[["percent.mt"]] <- PercentageFeatureSet(DK1, pattern = "^mt-")
DK2[["percent.mt"]] <- PercentageFeatureSet(DK2, pattern = "^mt-")
OR_1[["percent.mt"]] <- PercentageFeatureSet(OR_1, pattern = "^mt-")
OR_2[["percent.mt"]] <- PercentageFeatureSet(OR_2, pattern = "^mt-")

pdf("0-feature.pdf")
VlnPlot(CTL1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
VlnPlot(CTL2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
VlnPlot(DK1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
VlnPlot(DK2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
VlnPlot(OR_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
VlnPlot(OR_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

CTL1 <- subset(CTL1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mt < 50)
CTL2 <- subset(CTL2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mt < 50)
DK1 <- subset(DK1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 500 & nCount_RNA < 30000 & percent.mt < 50)
DK2 <- subset(DK2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mt < 50)
OR_1 <- subset(OR_1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 500 & nCount_RNA < 20000 & percent.mt < 50)
OR_2 <- subset(OR_2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 500 & nCount_RNA < 25000 & percent.mt < 50)

pdf("1-feature.pdf")
VlnPlot(CTL1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
VlnPlot(CTL2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
VlnPlot(DK1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
VlnPlot(DK2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
VlnPlot(OR_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
VlnPlot(OR_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

CTL1 <- SCTransform(CTL1, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
CTL2 <- SCTransform(CTL2, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
DK1 <- SCTransform(DK1, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
DK2 <- SCTransform(DK2, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
OR_1 <- SCTransform(OR_1, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
OR_2 <- SCTransform(OR_2, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)

object_list = list(CTL1,CTL2,DK1,DK2,OR_1,OR_2)
selfeatures <- SelectIntegrationFeatures(object.list = object_list, nfeatures = 3000)
scc.list <- PrepSCTIntegration(object.list = object_list, anchor.features = selfeatures, verbose = FALSE)
scc.anchors <- FindIntegrationAnchors(object.list = scc.list, normalization.method = "SCT",anchor.features = selfeatures, verbose = FALSE)
scc_integrated <- IntegrateData(anchorset = scc.anchors, normalization.method = "SCT",verbose = FALSE)
pdf("RA-feature.pdf")
VlnPlot(scc_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
dev.off()
saveRDS(scc_integrated,file="RA_mt0.5_SCT_combine.sample_S2.rds")

# ================ 2.Annotation and cell type analysis =========================
# 2.1 Dimensionality reduction and clustreing
data<-readRDS("./RA_mt0.5_SCT_combine.sample_S2.rds")
i=45           #// Reduction
data <- ScaleData(data, verbose = FALSE)
DefaultAssay(data)<-"integrated"
data <- RunPCA(data, npcs = i, verbose = FALSE)
data <- RunUMAP(data, reduction = "pca", dims = 1:i)
#data <- RunTSNE(data,dims = 1:i)
data <- FindNeighbors(data, reduction = "pca", dims = 1:i)
data <- FindClusters(data, resolution = 0.5)
saveRDS(file="dataintegratedPCA45.rds",data)

# 2.2 Celltype annotation
# PT
gene<-c("Slc13a3", "Slc34a1","Gpx3","Dcxr","Haf4a","Slc22a12","Slc27a2","Lrp2","Atp11a")  # 小鼠
gene<-c("Keg1","G6pc", "Alpl","Slc13a1","Scin")
gene<-c("Slc22a7", "Fxyd2", "Hrsp12", "Slc17a3")  #PST
#DCT
gene<-c("Slc12a3","Pvalb","Emx1","Cwh43")
# CD-IC 集合管间介细胞
gene<-c("Atp6v0d2","Atp6v1g3", "Slc4a1","Aqp6", "Slc26a4","Hmx2","Tmem61")
# Endo
gene<-c("Nrp1","Kdr","Ehd3","Plat","Flt1","Eng","Ptprb")
# Pericytes and vascular smooth muscle  Peri 血管内皮
gene<-c("Vim","S100a4")
# DLH
gene<-c("Aqp1","Bst1","Proser2")
# Podo 足细胞
gene<-c("Nphs1","Nphs2","Cdkn1c","Bcam","Wt1","Itga3","Lamb2")
# ALH 亨勒升环
gene<-c("Slc12a1","Umod","Cldn8","Cldn10","Krt18","Krt8")
# CD-PC 集合管主体 / 上皮细胞
gene<-c("Aqp2","Hsd11b2","Tmem45b")
# CD-Trans 集合管移行细胞
gene<-c("Rhbg","Insrr","Parm1","Sec23b")
# Novel 新型细胞
gene<-c("Cdca3","Mki67","Slc27a2","Lrp2","Stmn1")
# Fib   成纤维细胞
gene<-c("S100a4","Plac8","Mmp2","Emilin1","Sfrp2")
gene<-c("S100a4","Plac8","Mmp2","Emilin1","Sfrp2","Hsp47","Col1a1", "Col1a2", "Col5a1", "Loxl1", "Lum", "Fbln1","Fbln2")   # fibroblast
#MyoFib 髓系成纤维细胞
gene<-c("Acta2","Pdgfrb")
#B      B淋巴细胞
gene<-c("Cd79a","Cd79b","Ly6d")
# T/NK
gene<-c("Cxcr6","Ltb","Il7r","Cd3d","Cd3e", "Gzma","Nkg7","Gnly")
# Macro 巨噬细胞
gene<-c("C1qa","C1qb","Ncf2","Lyz2","Alox5ap")
# Mono 单核细胞
gene<-c("Lyz","Cd14")
# Neutro        中性粒细胞
gene<-c("S100a8","S100a9")

data@meta.data <- data@meta.data %>% mutate(Type=case_when(orig.ident == "CTL1" ~ "Normal1",
                                                           orig.ident == "CTL2" ~ "Normal2",
                                                           orig.ident == "DK1" ~ "DKD1",
                                                           orig.ident == "DK2" ~ "DKD2",
                                                           orig.ident == "OR_1" ~ "DKD_RA1",
                                                           orig.ident ==  "OR_2" ~ "DKD_RA2")
)
data@meta.data <- data@meta.data %>% mutate(Groups=case_when(orig.ident == "CTL1" ~ "Normal",
                                                             orig.ident == "CTL2" ~ "Normal",
                                                             orig.ident == "DK1" ~ "DKD",
                                                             orig.ident == "DK2" ~ "DKD",
                                                             orig.ident == "OR_1" ~ "DKD-RA",
                                                             orig.ident == "OR_2" ~ "DKD-RA"
))

data@meta.ec.integrated$Celltype="PT"
data@meta.ec.integrated$Celltype[ec.integrated$seurat_clusters %in% c(0,1,2,3,4,5,6,8,9,11,15,27)]="PT"
data@meta.ec.integrated$Celltype[ec.integrated$seurat_clusters %in% c(7,24)]="ALH"
data@meta.ec.integrated$Celltype[ec.integrated$seurat_clusters %in% c(12,17)]="Macro"
data@meta.ec.integrated$Celltype[ec.integrated$seurat_clusters %in% c(10,14)]="DCT"
data@meta.ec.integrated$Celltype[ec.integrated$seurat_clusters %in% c(13)]="Endo"
data@meta.ec.integrated$Celltype[ec.integrated$seurat_clusters %in% c(18)]="T/NK"
data@meta.ec.integrated$Celltype[ec.integrated$seurat_clusters %in% c(20)]="B"
data@meta.ec.integrated$Celltype[ec.integrated$seurat_clusters %in% c(23)]="CD-PC"
data@meta.ec.integrated$Celltype[ec.integrated$seurat_clusters %in% c(21)]="Neutro"
Celltype[data$seurat_clusters %in% c(16,19)]="CD-IC"
data@meta.data$Celltype[data$seurat_clusters %in% c(25)]="Poto"
data@meta.data$Celltype[data$seurat_clusters %in% c(22)]="DLH"
data@meta.data$Celltype[data$seurat_clusters %in% c(26)]="Novel"
data$Celltype<-factor(data$Celltype,levels=c("PT","DLH","DCT","ALH","CD-PC","CD-IC","Poto","Endo","Macro","Neutro","T/NK","B","Novel"))
data$Groups<-factor(data$Groups,levels=c("Normal","DKD","DKD-RA"))
data$Type<-factor(data$Type,levels=c("Normal1", "Normal2", "DKD1", "DKD2","DKD_RA1","DKD_RA2"))
sceasy::convertFormat(data, from="seurat", to="anndata", outFile='DbdataAnno.h5ad',dtype="float32")

# 2.2 Cellular composition and gene differential expression analysis
Cellratio<-prop.table(table(data$orig.ident,data$seurat_clusters), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))

# Correlation analysis between public scRNA and bulk RNAseq and data in the research
pdata<-readRDS("./publick_Data/Seurat_object_mouseDB70.RDS") #// 已发表代谢
raw.data <- read.table("/szrmyy/Drwang/scRNA/chenjunhui/P00010_Drzhang/00.data/publick_Data/GSE107585_count.txt",sep="\t", header=T, row.names = 1)
main_tiss <- CreateSeuratObject(counts = raw.data)
metadata<-read.table("/szrmyy/Drwang/scRNA/chenjunhui/P00010_Drzhang/00.data/publick_Data/GSE107585_metadata.txt",sep="\t",header=T,row.names = 1)
metadata<-data.frame((t(metadata)))
main_tiss <- AddMetaData(object = main_tiss, metadata = metadata)
saveRDS(file="ScienceNomaleMouse.RDS", main_tiss)
data<-readRDS("dataintegratedPCA45Anno.rds")
sdata<-readRDS("ScienceNomaleMouse.RDS")
markers<-read.csv(file="2-Celltype_markers.csv", header=T)
bulkRNA<-read.csv(file="./Kidneydissect.csv",header = T, row.names = 1)

nonimmuData<-subset(data,Celltype %in% c("PT","DLH", "DCT", "ALH", "CD-PC", "CD-IC","Poto", "Endo","Novel"))
exprMax <- nonimmuData@assays$RNA@data
exprMax<-as.data.frame(exprMax)
interbulkGenes <- intersect(rownames(bulkRNA), markers$gene)
testcor <- cor(exprMax[interbulkGenes,],bulkRNA[interbulkGenes,])
metadataSorted <- nonimmuData@meta.data  %>% arrange(Celltype)
testcor <- testcor[row.names(metadataSorted),]
annotation_col = data.frame(
  celltype = factor(metadataSorted$Celltype)
)
rownames(annotation_col)=as.character(rownames(metadataSorted))
ann_colors=list(celltype=c(PT= "#0073C2FF", DLH="#EFC000FF", DCT="#868686FF", ALH="#CD534CFF", "CD-PC"="#7AA6DCFF", "CD-IC"="#003C67FF",
                           Poto="#8F7700FF", Endo="#3B3B3BFF",  Novel="#00A1D5FF"))  #Macro="#A73030FF", Neutro="#4A6990FF", "T/NK"="#374E55FF", B="#DF8F44FF")
pdf("1-bulkcorrelationship.pdf")
pheatmap(testcor,cluster_rows=F, cluster_cols=F,show_rownames =F, annotation_row=annotation_col,annotation_colors=ann_colors)
dev.off()

# // 已发表文章对应关系
DefaultAssay(data) <-"RNA"
data<-SCTransform(data)
sdata$seurat_clusters<-sdata$Cluster_Number
sdata@meta.data$Celltype="PT"
sdata@meta.data$Celltype[sdata$seurat_clusters %in% c(3)]="PT"
sdata@meta.data$Celltype[sdata$seurat_clusters %in% c(4)]="LOH"
sdata@meta.data$Celltype[sdata$seurat_clusters %in% c(11)]="Macro"
sdata@meta.data$Celltype[sdata$seurat_clusters %in% c(5)]="DCT"
sdata@meta.data$Celltype[sdata$seurat_clusters %in% c(1)]="Endo"
sdata@meta.data$Celltype[sdata$seurat_clusters %in% c(14,15)]="T/NK"
sdata@meta.data$Celltype[sdata$seurat_clusters %in% c(13)]="B"
sdata@meta.data$Celltype[sdata$seurat_clusters %in% c(6)]="CD-PC"
sdata@meta.data$Celltype[sdata$seurat_clusters %in% c(12)]="Neutro"
sdata@meta.data$Celltype[sdata$seurat_clusters %in% c(7)]="CD-IC"
sdata@meta.data$Celltype[sdata$seurat_clusters %in% c(8)]="CD-Trans"
sdata@meta.data$Celltype[sdata$seurat_clusters %in% c(2)]="Poto"
sdata@meta.data$Celltype[sdata$seurat_clusters %in% c(9, 16)]="Novel"
sdata@meta.data$Celltype[sdata$seurat_clusters %in% c(10)]="Fib"
sdata<- NormalizeData(object=sdata)
sdata<-ScaleData(object = sdata)
sdata<-SCTransform(sdata)
sdata <- FindVariableFeatures(sdata, selection.method = "vst", nfeatures = 2000)


sdata <- RunPCA(sdata, npcs = 20, verbose = FALSE)
sdata <- RunUMAP(sdata, reduction = "pca", dims = 1:20)
sdata <- RunTSNE(sdata,dims = 1:20)
sdata <- FindNeighbors(sdata, reduction = "pca", dims = 1:20)
sdata <- FindClusters(sdata, resolution = 0.5)
p2 <- DimPlot(sdata, reduction = "tsne", label = TRUE, repel = TRUE)+ggtitle(paste("PCA ",20,sep=""))+
  theme(plot.title = element_text(hjust = 0.5))
Idents(sdata)<-"Celltype"
saveRDS(sdata,file="ScienceNomaleMouse.RDS")
markers<-read.csv("2-Celltype_markers.csv", header = T, stringsAsFactors = F)
sdata<-readRDS("ScienceNomaleMouse.RDS")

data<-seuratdata
data@meta.data$tmpCt<-as.character(data@meta.data$Celltype)
data@meta.data$tmpCt[data@meta.data$tmpCt %in% c("DLH","ALH")]<-"LOH"
data@meta.data$tmpCt<-factor(data@meta.data$tmpCt, levels=c("PT","DCT","LOH","CD-PC","CD-IC","Poto","Endo","Macro","Neutro","T/NK","B","Novel"))
Idents(data)<-"tmpCt"
DefaultAssay(data)<-"RNA"
data<- NormalizeData(object=data)
dat<-AverageExpression(data)
dat<-dat$RNA
sdata<-subset(sdata,Celltype %in% c("PT", "DCT", "CD-PC", "B", "T/NK", "LOH", "Endo", "CD-IC", "Neutro", "Macro", "Novel", "Poto"))
sdata$Celltype<-factor(sdata$Celltype,levels=c("PT","DCT","LOH","CD-PC","CD-IC","Poto","Endo","Macro","Neutro","T/NK","B","Novel"))
DefaultAssay(sdata)<-"RNA"
Idents(sdata)<-"Celltype"
sdat<-AverageExpression(sdata)
sdat<-sdat$RNA
interGenes <- intersect(rownames(sdat), markers$gene)
ct2sdat<-cor(dat[interGenes,], sdat[interGenes,])
color_1 <- colorRampPalette(c("#133059", "#FFFFFE","#631120"))(50)
pdf("1-scienceCorrelationship.pdf")
corrplot(ct2sdat,method = "color",tl.col= "black",type="lower")    # c("circle", "square", "ellipse", "number", "shade", "color", "pie")
# corrplot.mixed(matrix, lower = "number", upper = "pie", tl.col = "black",lower.col = "black", number.cex = 1)
dev.off()

# gene differential expression analysis
markers <- FindAllMarkers(data, assay = "RNA",slot= "data", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <-subset(markers, p_val_adj<0.05)
markers <-read.csv(file="2-Celltype_markers.csv",header=T)

top10 <- markers  %>%  group_by(cluster)  %>%  top_n(n = 50, wt = avg_log2FC)
write.csv(markers,"2-Celltype_markers.csv",quote=F,row.names=F)
dat<-AverageExpression(data)
dat<-dat$RNA
celltype<-unique(Idents(data))
annotation_col = data.frame(
  celltype = factor(celltype)
)
rownames(annotation_col)=as.character(celltype)
ann_colors=list(celltype=c(PT= "#0073C2FF", DLH="#EFC000FF", DCT="#868686FF", ALH="#CD534CFF", "CD-PC"="#7AA6DCFF", "CD-IC"="#003C67FF",
                           Poto="#8F7700FF", Endo="#3B3B3BFF",  Novel="#00A1D5FF",Macro="#A73030FF", Neutro="#4A6990FF", "T/NK"="#374E55FF", B="#DF8F44FF", Novel="#00A1D5FF"))
pdf("2-celltype-markersTop50.pdf",height=10,width=8)
pheatmap(dat[top10$gene,],scale="row", cluster_rows = FALSE,cluster_cols = FALSE,annotation_col=annotation_col,annotation_colors=ann_colors,show_rownames = F,color=colorRampPalette(c("#133059", "#FFFFFE","#631120"))(1000),border_color=NA,angle_col=45)
dev.off()

Idents(data)<-"Groups"
DefaultAssay(data)<-"RNA"
Idents(data)<-"Celltype"
celltype = as.character(levels(data))
dat = data.frame(gene="",Difference=0,logFC=0,p=0,adj.P=0,Celltype="",State="")
for (i in celltype) {
  subdata=subset(data,idents=i)
  print(unique(subdata$Celltype))
  Idents(subdata)<-"Groups"
  marker<-FindMarkers(subdata,ident.1="DKD",ident.2="Normal",only.pos=F, min.pct=0.1,logfc.threshold=0.25)
  marker<- marker %>% mutate(Difference = pct.1-pct.2)
  submarker = data.frame(gene=row.names(marker),Difference=marker$Difference, logFC=marker$avg_log2FC,p=marker$p_val,adj.P=marker$p_val_adj,Celltype=i,State=ifelse(marker$p_val_adj<0.05,ifelse(marker$avg_log2FC>0.25,"Up",ifelse(marker$avg_log2FC< -0.25,"Down","No")),"No"))

  dat=rbind(dat,submarker)
}
dat<-dat[-1,]
write.csv(dat,"3-Celltype_DKDvsNormal.csv",quote=F,row.names=F)
dat<-read.csv(file="3-Celltype_DKD_RAvsDKD.csv",header=T)
dat$State<-factor(dat$State,levels=c("Down","No","Up"))
plots<-list(); j=0;

for(i in celltype){
  j=j+1;
  subdat<-subset(dat,Celltype==i)
  p<-ggplot(subdat,aes(x=Difference,y=logFC)) +
    geom_point(size=0.5, aes(color=State))+
    scale_color_manual(values=c(Down="blue",No="grey",Up="red"))+
    #geom_label_repel(data=subset(dat,State !="No"), aes(label=names), segment.size=0.25, size=2.5)+
    geom_vline(xintercept=0,linetype=2)+ geom_hline(yintercept = 0.25,linetype=2)+geom_hline(yintercept = -0.25,linetype=2)+
    theme_classic()+
    ggtitle(i)+
    theme(plot.title = element_text(hjust = 0.5),legend.position="none")
  plots[[j]]<-p
}

pdf("3-Celltypedegvioplot_DKD_RAvsDKD.pdf",width=10,height=8)
grid.arrange(grobs=plots, ncol=5)
dev.off()

Dkd<-read.csv(file="3-Celltype_DKDvsNormal.csv",header = T)
Dkdup<-subset(Dkd,State == "Up")
top2 <- Dkdup  %>%  group_by(Celltype)  %>%  top_n(n = 2, wt = logFC)
up_genes<-top2$gene
Dkddown<-subset(Dkd,State == "Down")
top2 <- Dkddown  %>%  group_by(Celltype)  %>%  top_n(n = 2, wt = logFC)
down_genes<-top2$gene
RA<-read.csv(file="3-Celltype_DKD_RAvsDKD.csv",header = T)
Dkdup<-subset(RA,State == "Up")
top2 <- Dkdup  %>%  group_by(Celltype)  %>%  top_n(n = 2, wt = logFC)
up_genes<-top2$gene
Dkddown<-subset(RA,State == "Down")
top2 <- Dkddown  %>%  group_by(Celltype)  %>%  top_n(n = 2, wt = logFC)
down_genes<-top2$gene


PlotGeneGroup2<-function(object, feature, splitby ){
  if (is.null(levels(object@meta.data[,splitby]))){
    object@meta.data[,splitby] <-factor(object@meta.data[,splitby], levels = names(table(object@meta.data[,splitby])))
  }
  dataplot<-DotPlot(object, features = feature, split.by =  splitby, cols = topo.colors(n=length((levels(object@meta.data[,splitby])))))  #material.heat
  dataplot<-dataplot$data
  dataplot$avg.exp<-scale(dataplot$avg.exp)
  dataplot$Cluster<-gsub( "_.*$", "", dataplot$id )
  dataplot$Disease<-gsub( ".*_", "", dataplot$id )
  dataplot$Disease<-factor(dataplot$Disease, levels = levels(object@meta.data[,splitby]))
  dataplot$Cluster<-factor(dataplot$Cluster, levels = levels(object))
  colnames(dataplot)[1:2]<-c('Avg.Exp', 'Pct.Exp')
  dotplot<-ggplot(dataplot, aes(y = Cluster, x = Disease)) +  geom_tile(fill="white", color="white") +
    geom_point(aes( colour=Avg.Exp, size =Pct.Exp))  +  scale_color_gradientn(colours  =  colorRampPalette(c('grey80','lemonchiffon1','indianred1','darkred'))(255)
    )+ scale_size(range = c(0, 10))+
    theme(axis.line = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(size = 16,hjust = 0.5, face = 'bold'),
          axis.text = element_text(size = 12),axis.title=element_text(size=8),legend.text=element_text(size=8),
          legend.title = element_text(size = 8),legend.position="right")+ylab("")+xlab("")+ggtitle(feature)
}

PlotMultiGeneGroup<-function(object, features, splitby){
  pb <- progress_bar$new(
    format = "  Ploting [:bar] :percent eta: :eta",
    clear = FALSE, total = length(features), width = 100)
  features=rev(features)
  plot_list<-list()
  for(i in 1:length(features)){
    pp<-PlotGeneGroup2(object = object, feature = features[i], splitby = splitby)
    plot_list[[i]]<-pp$data
    pb$tick()
    Sys.sleep(1 / length(features))
  }
  all_data<-do.call('rbind', plot_list)
  dotplot<-ggplot(all_data, aes(x = Disease, y = features.plot)) +  geom_tile(fill="white", color="white") +
    geom_point(aes( colour=Avg.Exp, size =Pct.Exp), alpha=0.9)  +  scale_color_gradientn(colours  =  colorRampPalette(c('grey80','lemonchiffon1','indianred1','darkred'))(255)
    )+ scale_size(range = c(0, 10))+
    theme(axis.line = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(size = 16,hjust = 0.5, face = 'bold'),
          axis.text = element_text(size = 12),axis.title=element_text(size=8),legend.text=element_text(size=8),
          legend.title = element_text(size = 8),legend.position="right",strip.text = element_text(size = 14,colour = 'yellow',face = 'bold'))+ylab("")+xlab("")+ggtitle('')+facet_wrap(~Cluster, ncol = length(levels(object)))
  print(dotplot)
}


dn_groups<-data  #// seurat 对象

p <- PlotMultiGeneGroup(object = dn_groups, features = up_genes, splitby = 'Groups')
g1 <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-t', g1$layout$name))
fills<-c(pal_jco()(10),pal_jama()(7))
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g1$grobs[[i]]$grobs[[1]]$childrenOrder))
  g1$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

png(filename =  '3-RAdotplot_mutliple1.png', width = 18, height = 8,units = 'in', res = 600)
grid.draw(g1)
dev.off()
pdf(file =  '3-RAdotplot_mutliple1.pdf', width = 18, height = 8)
grid.draw(g1)
dev.off()


p <- PlotMultiGeneGroup(object = dn_groups, features = down_genes, splitby = 'Groups')+
  scale_color_gradientn(colours  =  colorRampPalette(c('grey80','lightcyan','deepskyblue','blue4'))(255))

g2 <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-t', g2$layout$name))
fills<-c(pal_jco()(10),pal_jama()(7))
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
png(filename =  '3-RAdotplot_mutliple2.png', width = 18, height = 8,units = 'in', res = 600)
grid.draw(g2)
dev.off()
pdf(file =  '3-RAdotplot_mutliple2.pdf', width = 18, height = 8)
grid.draw(g2)
dev.off()

# deg number plot
deg<-read.csv(file="3-Celltype_DKD_RAvsDKD.csv",header=T)
deg<-subset(deg,State != "No")
deg<-as.data.frame(table(deg$Celltype))
deg<-deg %>% arrange(Freq)
deg$Var1<-factor(deg$Var1,levels=c(deg$Var1))
pdf("3-Celltype_DKD_RAvsDKDDegstat.pdf")
ggplot(data = deg,aes(x = Var1, y = Freq)) +
  geom_bar(stat = 'identity',
           fill = ifelse(deg$Freq>0,"#003C67FF" ,'gray40'), # 根据y值的正负设置颜色
           width = 0.9) +
  # 横纵坐标轴名称
  labs(x = 'Name', y = 'Value')+
  theme(panel.grid.major =element_blank(),  # 去除边框
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title =  element_text(size=12,face = "bold"),
        axis.text.x =   element_text(angle=90,  # 横坐标文字旋转九十度
                                     hjust = 0.5, # 调整横坐标文字位置
                                     size=15),  # 调整横坐标文字字号
        plot.margin=unit(rep(3,4),'lines'))+
  coord_flip()

number<-read.table(file="3-Celltype_DKD_RAvsDKD.stat.xls", sep="\t", header = T)
number<-number %>% arrange(desc(ctnumber))
number<-rbind(number[9,],number[-9,])
number$ctnumber<-factor(number$ctnumber,levels=c(number$ctnumber))
ggplot(data = number,aes(x = ctnumber, y = number)) +
  geom_bar(stat = 'identity',
           fill = ifelse(number$number>0,"#003C67FF" ,'gray40'), # 根据y值的正负设置颜色
           width = 0.9) +
  # 横纵坐标轴名称
  labs(x = 'Name', y = 'Value')+
  theme(panel.grid.major =element_blank(),  # 去除边框
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title =  element_text(size=12,face = "bold"),
        axis.text.x =   element_text(angle=90,  # 横坐标文字旋转九十度
                                     hjust = 0.5, # 调整横坐标文字位置
                                     size=15),  # 调整横坐标文字字号
        plot.margin=unit(rep(3,4),'lines'))+
  coord_flip()
dev.off()

# genes rescued by treatment
DKDdeg<-read.csv(file="3-Celltype_DKDvsNormal.csv",header=T)
DKDdeg<-subset(DKDdeg,State !="No")
RAdeg<-read.csv(file="3-Celltype_DKD_RAvsDKD.csv",header=T)
RAdeg<-subset(RAdeg,State !="No")
gene<-c("Lgfbp5","Ptpro", "Vegfa",
        "Dach1","Ddx1","Slc47a1",
        "Slc22a2","Dpep1","Slc34a1",
        "Slc7a9","Lrp2","Sypl2",
        "Slc6a13","A1cf","Umod",
        "Stc1","Mpped2","Prkag2",
        "Wdr37","Kcnq1","Wdr72",
        "Gatm","Nfkb1","Etv5","Skil","Cacna1s"
)


rescueGenes<-read.table(file="3-Celltype_DKD_rescuedGenes.xls",header=T, sep="\t")
numberplot<-as.data.frame(table(rescueGenes$Celltype))
numberplot<- numberplot %>% arrange(Freq)
numberplot$Var1<-factor(numberplot$Var1,levels=c(numberplot$Var1))
pdf("4-RescuedGenesbyRA.pdf")
ggplot(data = numberplot,aes(x = Var1, y = Freq)) +
  geom_bar(stat = 'identity',
           fill = ifelse(numberplot$Freq>0,"#003C67FF" ,'gray40'), # 根据y值的正负设置颜色
           width = 0.9) +
  # 横纵坐标轴名称
  labs(x = 'Name', y = 'Value')+
  theme(panel.grid.major =element_blank(),  # 去除边框
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title =  element_text(size=12,face = "bold"),
        axis.text.x =   element_text(angle=90,  # 横坐标文字旋转九十度
                                     hjust = 0.5, # 调整横坐标文字位置
                                     size=15),  # 调整横坐标文字字号
        plot.margin=unit(rep(3,4),'lines'))+
  coord_flip()
dev.off()

geneset<-read.table("./00.data/DKDgeneset.xls",sep="\t",header = T, stringsAsFactors = F)
geneset<-read.table("./00.data/OXIDATIVE_STRESS.txt",sep="\t",header = T, stringsAsFactors = F)
pdf("4-RescuedGenesbyRAheatmap.pdf")
rescueGenes<-read.table(file="3-Celltype_DKD_rescuedGenes.xls",header=T, sep="\t")

allrescued<-rescueGenes
allrescued<-data.frame(gene=allrescued$gene,Celltype=allrescued$Celltype,logFC=allrescued$logFC)
allrescued %>% arrange(Celltype,logFC)->allrescued
widedata3 <- spread(allrescued,key = "Celltype",value = "logFC")
widedata3[is.na(widedata3)] <- 0
rownames(widedata3)<-widedata3$gene
widedata3<-widedata3[,-1]
widedata3<-as.matrix(widedata3)
dscaled_mat = t(scale(t(widedata3)))
d = dist(dscaled_mat, method = 'euclidean')
tree = hclust(d, method = 'complete')
genelist <- row.names(dscaled_mat)[tree$order]
dscaled_mat<-dscaled_mat[genelist,celltypeName]
vpt<-geneset[geneset$gene %in% rescueGenes$gene,]
DKDgenes<-unique(vpt[vpt$term %in% c("TGF-β","Inflammation"),]$gene)
CelltypeOrder<-sort(table(rescueGenes$Celltype),decreasing = TRUE)
celltypeNumber<-as.numeric(unname(CelltypeOrder))
celltypeName<-names(CelltypeOrder)

col_fun = colorRamp2(c(-3,0,3), c("purple","white", "orange"))
ct<-colnames(dscaled_mat)
celltypeCol=c(PT= "#0073C2FF", DLH="#EFC000FF", DCT="#868686FF", ALH="#CD534CFF", "CD-PC"="#7AA6DCFF", "CD-IC"="#003C67FF",
              Endo="#3B3B3BFF",  Novel="#00A1D5FF",Macro="#A73030FF",  "T/NK"="#374E55FF")
column_ha = HeatmapAnnotation(Celltype=ct,
                              show_annotation_name = TRUE,
                              annotation_name_gp = gpar(fontsize = 7),
                              show_legend = TRUE,
                              col = list(Celltype = celltypeCol),bar1 = anno_barplot(celltypeNumber))
dpos<-which(rownames(dscaled_mat) %in% DKDgenes)
genelabel=rownames(dscaled_mat)[rownames(dscaled_mat) %in% DKDgenes]
dha = rowAnnotation(foo = anno_mark(at =dpos, labels = genelabel))
Heatmap(dscaled_mat, name = "mat",col=col_fun, show_row_names=F, cluster_rows = FALSE, cluster_columns = FALSE,right_annotation = dha, top_annotation =column_ha,
        row_names_side = "left", row_names_gp = gpar(fontsize = 4))
#热图2
RArescueGenes<-read.table(file="3-Celltype_RA_rescuedGenes.xls",header=T, sep="\t")
allrescued<-RArescueGenes
allrescued<-data.frame(gene=allrescued$gene,Celltype=allrescued$Celltype,logFC=allrescued$logFC)
allrescued %>% arrange(Celltype,logFC)->allrescued
widedata3 <- spread(allrescued,key = "Celltype",value = "logFC")
widedata3[is.na(widedata3)] <- 0
rownames(widedata3)<-widedata3$gene
widedata3<-widedata3[,-1]
widedata3<-as.matrix(widedata3)
dscaled_mat = t(scale(t(widedata3)))
dscaled_mat<-dscaled_mat[genelist,celltypeName]
col_fun = colorRamp2(c(-3,0,3), c("purple","white", "orange"))
ct<-colnames(dscaled_mat)
celltypeCol=c(PT= "#0073C2FF", DLH="#EFC000FF", DCT="#868686FF", ALH="#CD534CFF", "CD-PC"="#7AA6DCFF", "CD-IC"="#003C67FF",
              Endo="#3B3B3BFF",  Novel="#00A1D5FF",Macro="#A73030FF",  "T/NK"="#374E55FF")
column_ha = HeatmapAnnotation(Celltype=ct,
                              show_annotation_name = TRUE,
                              annotation_name_gp = gpar(fontsize = 7),
                              show_legend = TRUE,
                              col = list(Celltype = celltypeCol),bar1 = anno_barplot(celltypeNumber))
dpos<-which(rownames(dscaled_mat) %in% DKDgenes)
genelabel=rownames(dscaled_mat)[rownames(dscaled_mat) %in% DKDgenes]
dha = rowAnnotation(foo = anno_mark(at =dpos, labels = genelabel))
Heatmap(dscaled_mat, name = "mat",col=col_fun, show_row_names=F, cluster_rows = FALSE, cluster_columns = FALSE,right_annotation = dha, top_annotation =column_ha,
        row_names_side = "left", row_names_gp = gpar(fontsize = 4))
dev.off()

seuratdata<-readRDS("./01.anno/dataintegratedPCA45Anno.rds")
seuratdata$Celltype=factor(seuratdata$Celltype, levels = c("PT","DLH", "DCT","ALH","CD-PC", "CD-IC","Poto","Endo","Macro","Neutro", "T/NK", "B", "Novel"))
subsetData<-subset(seuratdata,Celltype %in% c("PT","DCT","ALH","CD-IC","Endo","CD-PC","Macro"))
prefix="OXIDATIVE_STRESS"
gmtfile="./GOBP_RESPONSE_TO_OXIDATIVE_STRESS.v2022.1.Hs.gmt"
genes<-read.gmt(gmtfile)
DefaultAssay(seuratdata) <- "RNA"
pgene<-capitalize(tolower(genes$gene))
gene=list(capitalize(tolower(genes$gene)))
Idents(seuratdata)<-"Celltype"
pbmc <- AddModuleScore(
  object = seuratdata,
  features = gene,
  ctrl = 100, #默认值是100
  assay="RNA",
  name = 'GLY'
)
colnames(pbmc@meta.data)[12] <- 'GLY'
data<- FetchData(pbmc, vars = c("Groups", "Celltype", "GLY"))
for (i in 1:dim(data)[1]){
  if (data[i,1] == "Normal"){
    data[i,1] = "Normal"
  }
  if (data[i,1] == "DKD"){
    data[i,1] = "DKD"
  }
  if (data[i,1] == "DKD_RA"){
    data[i,1] = "DKD_RA"
  }
}
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}


data2 <- data %>%
  group_by_at(.vars = c("Groups")) %>%
  mutate(GLY = case_when(TRUE ~ remove_outliers(GLY),TRUE ~ GLY))
data2 <- data2[,c(1,2,3)]
data2 <- melt(data2)

my_cols = c(pal_jco()(10),pal_jama()(4))
pdf(paste("3-",prefix,".Score.pdf",sep=""), width = 10, height = 4)
p <- ggplot(data2, aes(Groups, value))
p<-p+geom_boxplot( colour=rep(c("#CC33FFFF", "#6699FFFF","red"), 13) )+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_signif(comparisons = list(c("Normal", "DKD"),c("DKD","DKD_RA")),  map_signif_level = T,
              test = wilcox.test, size=0.5, textsize = 2, step_increase = 0.1) +
  xlab("") + ylab(prefix) +
  theme_bw() +
  facet_wrap( ~Celltype , ncol = 6 ) +  #按照细胞类型分列
  theme( strip.text = element_text(colour="white")  )  +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))


g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- my_cols
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)
dev.off()
# key genes
gene<-c("Acaa1a","Cyc1","Alas1","Atp5k","Ndufa9","Mpc1",
        "Kap", "Rps29","Serpina3g","Pnp",   #Ostree
        "Smad7","Id3","Tgfbr1","Pmepa1"
)
rescuedGenes<-read.table(file="3-Celltype_All_rescuedGenes.xls",header=T)
DKDup<-subset(rescueGenes, State!="No")
gene<-unique(DKDup$gene)
pathgene<-read.table(file="./genes.txt",header=T)
gene<-gene[(gene %in% pathgene$gene)]
gene<-as.character(gene)

png(filename =  '3-FocusrescuedGenesdotplot_multiple_migration.png', width = 10, height = 6,units = 'in', res = 300)
complex_dotplot_multiple(seu_obj = subsetData, features = gene,group = "Groups", celltypes = c("PT","DCT","ALH","CD-IC","Endo","CD-PC","Macro"))
dev.off()
png(filename =  'vlnplot_multiple_genes.png', width = 6, height = 6,units = 'in', res = 300)
complex_vlnplot_multiple(iri.integrated, features = c("Havcr1",  "Slc34a1", "Vcam1",   "Krt20"  , "Slc7a13", "Slc5a12"), celltypes = c("PTS1" ,   "PTS2"  ,  "PTS3"  ,  "NewPT1" , "NewPT2"), group = "Group", add.dot=T, pt.size=0.01, alpha=0.01, font.size = 10)
dev.off()


allrescuedGene<-read.table(file="3-Celltype_All_rescuedGenes.xls",sep="\t",header=T)
DKDup<-subset(allrescuedGene,RescuedType=="DKDUp")
DKDdown<-subset(allrescuedGene,RescuedType=="DKDDown")
subset(geneset,path=="NFKB2")           # Glycolysis / Gluconeogenesis, TGF-β, Apoptosis，AGE-RAGE，Tyrosine metabolism，Steroid hormone biosynthesis
subset(geneset,path=="ROS(GSEA)")
subset(geneset,path=="Inflammation")
subset(geneset,path=="Ros generation")
subset(geneset,path=="PI4K_AKT_mTOR")

genes<-allrescuedGene$gene
databases = c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021","GO_Biological_Process_2021","KEGG_2021_Mouse")
enriched <- enrichR::enrichr(genes, databases)
contrast_enrich <- NULL
for(database in databases) {
  temp <- enriched[database][[1]] %>%
    mutate(var = database)
  contrast_enrich <- rbind(contrast_enrich, temp)
}
contrast_enrich$contrast <- "AFTvsVehicle"
contrast_enrich$n <- length(genes)
write.table(file=paste("3-RescureGenesDEGs_BP_Enrichment_4Databases.xls"),contrast_enrich,sep="\t", quote = FALSE)

eg <- bitr(x, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")
ego_BP <- enrichGO(OrgDb="org.Mm.eg.db", gene = eg[,2], pvalueCutoff = 0.05, ont = "BP", readable=TRUE)
#'自定义输入gmt文件，至下载MsigGB数据库
gmtfile<-read.gmt("./00.data/gmt/REACTOME_kegg_hallmarker.gmt")
result <- enricher(x,TERM2GENE=gmtfile,pvalueCutoff = 0.05,pAdjustMethod = "BH")
write.table(file="3-RescuedGenesEnrichmentRKG.xls",data.frame(result ),sep="\t")

DegFun<-read.table(file="3-RescuedgenesFunctionplot.xls",sep="\t",header=T)
DegFun<-DegFun %>% arrange(GeneRatio)
DegFun$Description<-factor(DegFun$Description,levels=DegFun$Description)
pdf("3-RescuedgenesFunctionplot_Enrichment.pdf",width=7,height=8)
ggplot(DegFun,aes(GeneRatio,Description)) + geom_point(aes(size=Count,color=p.adjust)) + theme_classic() +
  scale_colour_gradientn(colors = c('#CD534CFF','#88b8d8')) +scale_size(limits=c(3,10))+ theme_minimal() +
  theme(axis.text.x=element_text(size=30,angle=90),axis.text.y = element_text(size = 30,color="black")) +
  theme_classic()
dev.off()


features<-list(subset(geneset,path=="Inflammation")$gene)
data<-AddModuleScore(
  data,
  features,
  pool = NULL,
  ctrl = 100, name = "Inflammation")
FeaturePlot(data, features = "Inflammation1", min.cutoff = "q10", max.cutoff = "q95",split.by = "Groups")
FeaturePlot(data, features = "Inflammation1", min.cutoff = "q10", max.cutoff = "q95")
VlnPlot(data, features="Inflammation1", group.by ="Celltype",split.by ="Groups",pt.size=0)

#pdf("4-LipidPFOAVSCONscore.pdf",width=6,height=6)
comparisons=list(c("Normal", "DKD"), c("DKD", "DKD-RA"))
p<-ggviolin(data@meta.data, x="Groups", y="Inflammation1", color = "Groups",palette = c('lancet'),add = "boxplot",add.params = list(fill="white"))+ stat_compare_means(comparisons = comparisons,label =  "p.signif")
p+facet_wrap( ~Celltype, nrow = 4 )+
  theme(strip.background = element_blank())
#dev.off()

# ================ 3.subtypes analysis ==========================
# 3.1 PT subtype
ec <- subset(data, idents="PT")
ec.list<-SplitObject(ec, split.by = 'orig.ident')
for (i in names(ec.list)) {
  ec.list[[i]] <- SCTransform(ec.list[[i]], verbose = F)
  print(paste(i, 'done!'))
}

ec.features <- SelectIntegrationFeatures(object.list = ec.list, nfeatures = 1000)
ec.list <- PrepSCTIntegration(object.list = ec.list, anchor.features = ec.features)
ec.anchors <- FindIntegrationAnchors(object.list = ec.list,
                                     anchor.features = ec.features, normalization.method = "SCT")
Sys.time()
ec.integrated <- IntegrateData(anchorset = ec.anchors)
Sys.time()
pc=30
ec.integrated<-ScaleData(ec.integrated)
ec.integrated <- RunPCA(object = ec.integrated, npcs = pc)
ec.integrated <- RunUMAP(object = ec.integrated, dims = 1:pc)
ec.integrated <- FindNeighbors(ec.integrated, dims = 1:pc)
ec.integrated <- FindClusters(ec.integrated, resolution = 0.6)
save(ec.integrated, file = 'ec.integrated.Rda')
ec.integrated@meta.data$Celltype="S1"
ec.integrated@meta.data$Celltype[ec.integrated$seurat_clusters %in% c(3)]="S3"
ec.integrated@meta.data$Celltype[ec.integrated$seurat_clusters %in% c(0,1)]="S2"
Idents(ec.integrated) <-"Celltype"

# 3.2 Other segments subtype
setwd("./01.anno/Othersegments")
ec<-subset(seuratdata, Celltype %in% c("DCT","ALH","CD-IC","CD_PC","DLH","Poto"))
ec.list<-SplitObject(ec, split.by = 'orig.ident')
for (i in names(ec.list)) {
  ec.list[[i]] <- SCTransform(ec.list[[i]], verbose = F)
  print(paste(i, 'done!'))
}

ec.features <- SelectIntegrationFeatures(object.list = ec.list, nfeatures = 1000)
ec.list <- PrepSCTIntegration(object.list = ec.list, anchor.features = ec.features)
ec.anchors <- FindIntegrationAnchors(object.list = ec.list,
                                     anchor.features = ec.features, normalization.method = "SCT")
Sys.time()
ec.integrated <- IntegrateData(anchorset = ec.anchors)
Sys.time()
bkec.integrated<-ec.integrated
pc=20
ec.integrated<-ScaleData(ec.integrated)
ec.integrated <- RunPCA(object = ec.integrated, npcs = pc)
ec.integrated <- RunUMAP(object = ec.integrated, dims = 1:pc)
ec.integrated <- FindNeighbors(ec.integrated, dims = 1:pc)
ec.integrated <- FindClusters(ec.integrated, resolution = 0.6)
save(ec.integrated, file = 'ec.integrated.Rda')
Sys.time()
DimPlot(ec.integrated,label=T)+ theme(plot.title = element_text(vjust = 0.5,hjust = 0.5,size = 15))
markers <- FindAllMarkers(ec.integrated, assay = "RNA",slot= "data", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <-subset(markers, p_val_adj<0.05)
top10 <- markers  %>%  group_by(cluster)  %>%  top_n(n = 5, wt = avg_log2FC)
write.table(file="5-Osmakers.xls",markers, sep="\t")
ec.integrated@meta.data$subType="DCT1"
ec.integrated@meta.data$subType[ec.integrated$seurat_clusters %in% c(1,3,4,6,15)]="DCT1"
ec.integrated@meta.data$subType[ec.integrated$seurat_clusters %in% c(2)]="DCT2"
ec.integrated@meta.data$subType[ec.integrated$seurat_clusters %in% c(10)]="B-IC"
ec.integrated@meta.data$subType[ec.integrated$seurat_clusters %in% c(18)]="A-IC"
ec.integrated@meta.data$subType[ec.integrated$seurat_clusters %in% c(14)]="CD-Trans"
ec.integrated@meta.data$subType[ec.integrated$seurat_clusters %in% c(8,9)]="CD-PC"
ec.integrated@meta.data$subType[ec.integrated$seurat_clusters %in% c(7,16,17)]="Podo"
ec.integrated@meta.data$subType[ec.integrated$seurat_clusters %in% c(0,5,12,13)]="ALH"
ec.integrated@meta.data$subType[ec.integrated$seurat_clusters %in% c(11)]="DLH"
ec.integrated$subType<-factor(ec.integrated$subType,levels=c("DLH","DCT1","DCT2","ALH","CD-PC","A-IC","CD-Trans","B-IC","Podo"))
Idents(ec.integrated)<-"subType"
saveRDS("5-ossubType.RDS",ec.integrated)

# 3.3 Macrophages analysis
Idents(data)<-"Celltype"
ec <- subset(data, idents="Macro")
ec.list<-SplitObject(ec, split.by = 'orig.ident')
for (i in names(ec.list)) {
  ec.list[[i]] <- SCTransform(ec.list[[i]], verbose = F)
  print(paste(i, 'done!'))
}

ec.features <- SelectIntegrationFeatures(object.list = ec.list, nfeatures = 1500)
ec.list <- PrepSCTIntegration(object.list = ec.list, anchor.features = ec.features)
ec.anchors <- FindIntegrationAnchors(object.list = ec.list,
                                     anchor.features = ec.features, normalization.method = "SCT")
Sys.time()
ec.integrated <- IntegrateData(anchorset = ec.anchors, k.weight = 60)
Sys.time()

pc=30
DefaultAssay(ec.integrated)="integrated"
ec.integrated<-ScaleData(ec.integrated)
ec.integrated <- RunPCA(object = ec.integrated, npcs = pc)
ec.integrated <- RunUMAP(object = ec.integrated, dims = 1:pc)
ec.integrated <- FindNeighbors(ec.integrated, dims = 1:pc)
ec.integrated <- FindClusters(ec.integrated, resolution = 0.6)
save(ec.integrated, file = 'MacrophageObject.Rda')
Sys.time()
ec.integrated<-subset(ec.integrated, seurat_clusters != 8)
ec.integrated@meta.data$subType[ec.integrated$seurat_clusters %in% c(0,1)]="Cd86_Macro"
ec.integrated@meta.data$subType[ec.integrated$seurat_clusters %in% c(3,5,6,7)]="S100a4_Macro"
ec.integrated@meta.data$subType[ec.integrated$seurat_clusters %in% c(2,4)]="Gpx3_Macro"
ec.integrated@meta.data$subType[ec.integrated$seurat_clusters %in% c(9)]="Cd163_Macro"
Idents(ec.integrated)<-"subType"
saveRDS(file="6-Macrophage_subtype.RDS",ec.integrated)

# 3.4 Lymphocytes  cell type analysis
seuratdata<-readRDS("./01.anno/dataintegratedPCA45Anno.rds")
ec <- subset(seuratdata, Celltype %in% c("T/NK","B"))
#ec<-FindVariableFeatures(ec, selection.method = "vst", nfeatures = 1500)
ec.list<-SplitObject(ec, split.by = 'Groups')
for (i in names(ec.list)) {
  ec.list[[i]] <- SCTransform(ec.list[[i]], verbose = F)
  print(paste(i, 'done!'))
}

ec.features <- SelectIntegrationFeatures(object.list = ec.list, nfeatures = 1000)
ec.list <- PrepSCTIntegration(object.list = ec.list, anchor.features = ec.features)
ec.anchors <- FindIntegrationAnchors(object.list = ec.list,
                                     anchor.features = ec.features, normalization.method = "SCT")
Sys.time()
ec.integrated <- IntegrateData(anchorset = ec.anchors,k.weight =50)
Sys.time()


pc=30
ec.integrated<-ScaleData(ec.integrated)
ec.integrated <- RunPCA(object = ec.integrated, npcs = pc)
ec.integrated <- RunUMAP(object = ec.integrated, dims = 1:pc)
ec.integrated <- FindNeighbors(ec.integrated, dims = 1:pc)
ec.integrated <- FindClusters(ec.integrated, resolution = 1.8)
save(ec.integrated, file = 'ec.integrated.Rda')
Sys.time()
ec.integrated@meta.data$subType="T_C1"
ec.integrated@meta.data$subType[ec.integrated$seurat_clusters %in% c(0,2,5)]="B_C1"
ec.integrated@meta.data$subType[ec.integrated$seurat_clusters %in% c(9)]="B_C2"
ec.integrated@meta.data$subType[ec.integrated$seurat_clusters %in% c(10)]="B_C3"
ec.integrated@meta.data$subType[ec.integrated$seurat_clusters %in% c(6,7)]="NK"
ec.integrated@meta.data$subType[ec.integrated$seurat_clusters %in% c(1,3)]="T_C1"
ec.integrated@meta.data$subType[ec.integrated$seurat_clusters %in% c(8)]="T_C2"
ec.integrated@meta.data$subType[ec.integrated$seurat_clusters %in% c(4,11)]="T_C3"
ec.integrated@meta.data$subType[ec.integrated$seurat_clusters %in% c(12)]="T_C4"

ec.integrated$subType<-factor(ec.integrated$subType,levels=c("T_C1","T_C2","T_C3","T_C4","B_C1","B_C2","B_C3","NK"))

Idents(ec.integrated)<-"subType"
saveRDS(file="ec.integrated.Rda",ec.integrated)

# ================ 4. Cellular crosstalk ==========================

data <- readRDS("./dataintegratedPCA45AnnosubType.rds")
Idents(data)<-"subType"
data <- subset(data, subType != "Macro")
meta<-data@meta.data
mat<-data@assays$RNA@data

Control_meta<-subset(meta,Groups=="Normal")
Control_meta$subType = droplevels(Control_meta$subType, exclude = setdiff(levels(Control_meta$subType),unique(Control_meta$subType)))
Treatment_meta<-subset(meta,Groups=="DKD")
Treatment_meta$subType = droplevels(Treatment_meta$subType, exclude = setdiff(levels(Treatment_meta$subType),unique(Treatment_meta$subType)))
RADKD_meta<-subset(meta,Groups=="DKD_RA")
RADKD_meta$subType = droplevels(RADKD_meta$subType, exclude = setdiff(levels(RADKD_meta$subType),unique(RADKD_meta$subType)))

Control_data<-data@assays$RNA@data[,as.character(row.names(Control_meta))]
Treatment_data<-data@assays$RNA@data[,as.character(row.names(Treatment_meta))]
RADKD_data <-data@assays$RNA@data[,as.character(row.names(RADKD_meta))]

cellchat <- createCellChat(object = Treatment_data, meta = Treatment_meta, group.by = "subType")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

pdf("DKD_CellChat.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
saveRDS(cellchat, "DKD_cellchat.RDS")


##Normal
cellchat <- createCellChat(object = Control_data, meta = Control_meta, group.by = "subType")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
#Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

pdf("Normal_CellChat.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
saveRDS(cellchat, "Normal_cellchat.RDS")


# RA DKD
cellchat <- createCellChat(object = RADKD_data, meta = RADKD_meta, group.by = "subType")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
#Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

pdf("DKD_RA_CellChat.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
saveRDS(cellchat, "DKD_RA_cellchat.RDS")

setwd("./02.Cellchat/RunCellchat")
cellchat.Normal <- readRDS("./Normal_cellchat.RDS")
cellchat.DKD <- readRDS("./DKD_cellchat.RDS")
cellchat.DKD_RA <- readRDS("./DKD_RA_cellchat.RDS")
object.list <- list(Normal =cellchat.Normal, DKD = cellchat.DKD, DKD_RA=cellchat.DKD_RA)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight")
gg1 + gg2


gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
pdf("1-NumberOfinteraction.pdf")
print(gg1)
dev.off()

cellchat.Normal <- updateCellChat(cellchat.Normal)
mat1 <- cellchat.Normal@net$weight
groupSize1 <- as.numeric(table(cellchat.Normal@idents))

cellchat.DKD <- updateCellChat(cellchat.DKD)
mat2 <- cellchat.DKD@net$weight
groupSize2 <- as.numeric(table(cellchat.DKD@idents))

cellchat.DKD_RA <- updateCellChat(cellchat.DKD_RA)
mat3 <- cellchat.DKD_RA@net$weight
groupSize3 <- as.numeric(table(cellchat.DKD_RA@idents))

tName <- c("Normal", "DKD","DKD_RA")
edgeMax<-15
for(i in c(21)){  #  c(1, 8, 14,15, 21
  pdf("1-NKnetVisual_circle.pdf",width = 10,height=5)
  par(mfrow = c(1,3), xpd=TRUE)
  matrix2 <- matrix(0, nrow = nrow(mat1), ncol = ncol(mat1), dimnames = dimnames(mat1))
  matrix2[i, ] <- mat1[i, ]
  netVisual_circle(matrix2, vertex.weight = groupSize1, weight.scale = T, edge.weight.max = max(mat1), edge.width.max =edgeMax, title.name = tName[1])

  matrix2 <- matrix(0, nrow = nrow(mat2), ncol = ncol(mat2), dimnames = dimnames(mat2))
  matrix2[i, ] <- mat2[i, ]
  netVisual_circle(matrix2, vertex.weight = groupSize2, weight.scale = T, edge.weight.max = max(mat2), edge.width.max =edgeMax, title.name = tName[2])

  matrix2 <- matrix(0, nrow = nrow(mat3), ncol = ncol(mat3), dimnames = dimnames(mat3))
  matrix2[20, ] <- mat3[20, ]
  netVisual_circle(matrix2, vertex.weight = groupSize3, weight.scale = T, edge.weight.max = max(mat3), edge.width.max =edgeMax, title.name = tName[3])  #NK 20
  dev.off()
}

pairLR.use <- read.table(file="FilterLR.txt", header = T)

#  数据库筛选
CellChatDB <- CellChatDB.mouse
CellChatDB$interaction
pathways.show.all <- cellchat.Normal@netP$pathways

# Circle plot show L-R pairs
pathways.show <- c("CXCL","TGFb","IL2", "TNF","IL1","CCL","CD86")
pairLR.CXCL <- extractEnrichedLR(cellchat.Normal, signaling = pathways.show, geneLR.return = FALSE)
pairLR.use$interaction_name <- factor(pairLR.use$interaction_name, levels = pairLR.use$interaction_name)
df$interaction_name_2<- factor(df$interaction_name_2, levels = c("Mif  - (Cd74+Cd44)","Cd86  - Ctla4","Cd274  - Pdcd1","Ccl2  - Ccr2","Ccl3  - Ccr5","Ccl5  - Ccr5","Ccl6  - Ccr2","H2-q10  - Cd8a","H2-q10  - Cd8b1","Lamb3  - Cd44"))

pdf("2-BuppleLRPaire.pdf",width = 12,height = 6)
p <- netVisual_bubble(cellchat,
                      sources.use = c("Cd86_Macro","S100a4_Macro","NK"),
                      targets.use = c("T_C1","T_C2","T_C3","B_C1","Neutro","Cd163_Macro","CD_PC"),return.data=F,
                      comparison = c(1, 2,3),pairLR.use =pairLR.use, angle.x = 45, color.heatmap = "Spectral",max.quantile=0.8)    # pairLR.use = pairLR.use, angle.x = 45
print(p)
dev.off()

# ================ 5. data visualization ==========================
# The following R packages were used for data visualization

#   SCP (version 0.4.8):           https://zhanghao-njmu.github.io/SCP/
#   ggplot2 (version 3.4.2):       https://github.com/tidyverse/ggplot2
#   pheatmap (version 1.0.12):      https://github.com/raivokolde/pheatmap
#   ComplexHeatmap (version 2.12.0): https://github.com/jokergoo/ComplexHeatmap
#   plot1cell (version 0.0.0.9000):           https://github.com/TheHumphreysLab/plot1cell

# ================ 6. Other instructions ==========================
# Gene sets, including those related to oxidative stress and inflammation, utilized for calculating module scores, were sourced from MGI (Mouse Genome Informatics).
# Genes associated with injury were obtained based on research papers (Kirita Y et al., 2020, PNAS; Wu H, 2022, Cell Metabolism).


