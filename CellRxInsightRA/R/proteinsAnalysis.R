# The code employed for conducting the analysis of proteins in the study.
#

suppressPackageStartupMessages(library(DEP))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ReactomePA))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(enrichR))
suppressPackageStartupMessages(library(venn))

contrastP<-"Model_vs_con"
contrast<-paste(contrastP,"_significant",sep="")

dF <- read.xlsx('./03.Verify/ZQ_DN_mouse_LFQ_ROS.xlsx', sheetIndex=1)
dF[,6:14][is.na(dF[,6:14])] <- 0
dF<-dF[rowSums( dF[,6:14]> 0) >= 5,]
data <- dF

data$Genes %>% duplicated() %>% any()
data %>% group_by(Genes) %>% summarise(frequency = n()) %>%
  arrange(desc(frequency)) %>% filter(frequency > 1)
data[is.na(data$Genes),]$Genes<-data$Protein.Ids[is.na(data$Genes)]
data_unique <- make_unique(data,"Genes","Protein.Ids",delim = ";")
Intensity_columns <- 6:14
experiment_design<-read.table(file="./experiment_design.txt",sep="\t",header=T)   #// design files should be provided
data_se <- make_se(data_unique,Intensity_columns,experiment_design)

# Data QC
pdf("1-plot_frequency.pdf")
plot_frequency(data_se)
dev.off()
data_filter <- filter_missval(data_se,thr = 0)
pdf("1-plot_numbers.pdf")
plot_numbers(data_filter)
dev.off()
pdf("1-plot_coverage.pdf")
plot_coverage(data_filter)
dev.off()
data_norm <- normalize_vsn(data_filter)
pdf("1-plot_normalization.pdf")
plot_normalization(data_filter,data_norm)
dev.off()
pdf("1-plot_missval.pdf")
plot_missval(data_filter)
dev.off()
pdf("1-plot_detect.pdf")
plot_detect(data_filter)
dev.off()

data_imp <- impute(data_norm,fun = "MinProb",q = 0.01)
pdf("1-plot_imputation.pdf")
plot_imputation(data_imp)
dev.off()
data_diff <- test_diff(data_imp, type = "control", control = "con")
dep <- add_rejections(data_diff,alpha = 1,lfc = 0)

pdf("1-plot_pca.pdf")
plot_pca(dep,x = 1,y = 2,n = 500,point_size = 4)
dev.off()
assay(dep)[is.na(assay(dep))]<-0
source("./plot_volcano.R")
pdf("1-DEPDistributio_plt0.05_logFCgt1.3.pdf")
plot_cor(dep,significant = TRUE,lower = 1,upper = -1,pal = rev("RdBu"))

plot_heatmap(dep,type = "centered",kmeans = TRUE,k = 6,col_limit = 4,
             show_row_names = TRUE,indicate = c("condition","replicate"))
plot_heatmap(dep,type = "contrast",kmeans = TRUE,k = 6,col_limit = 10,show_row_names = FALSE)
plot_volcano(dep,contrast = contrastP,label_size = 2,add_names = TRUE,topgenes = 10)
plot_single(dep, proteins = "ANAPC13", type = "centered")
dev.off()

data_results <- get_results(dep)
data_results %>% filter(significant) %>% nrow()
df_wide <- get_df_wide(dep)
df_long <- get_df_long(dep)
write.table(file="2-DeGprotein_modelVsother.xls",df_wide,sep="\t")

save(data_se, data_norm, data_imp, data_diff, dep, file = "data.RData")


#========================================== function enrichment  ===========================================

databases = c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021","GO_Biological_Process_2021")
gsea_results_per_contrast <- test_gsea(dep, databases = c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021","GO_Biological_Process_2021"))
write.table(gsea_results_per_contrast,file="3-test_plt0.05_logFCgt1.3gseaResult.xls",sep="\t",row.names = FALSE,quote = FALSE)

row_data <- rowData(dep, use.names = FALSE)
# Run background list
message("Background")
background <- gsub("[.].*", "", row_data$name)
backgroud<-gsub(".*[(]", "", background)
background_enriched <- enrichR::enrichr(background, databases)
df_background <- NULL
for(database in databases) {
  temp <- background_enriched[database][[1]] %>%
    mutate(var = database)
  df_background <- rbind(df_background, temp)
}
df_background$contrast <- "background"
df_background$n <- length(background)

OUT <- df_background %>%
  mutate(bg_IN = as.numeric(gsub("/.*", "", Overlap)),
         bg_OUT = n - bg_IN)  %>%
  dplyr::select(Term, bg_IN, bg_OUT)
write.table(df_background,file="2-AllinformationbackgroudEnrichment.xls",sep="\t",row.names = FALSE,quote = FALSE)

df <- row_data %>%
  as.data.frame() %>%
  dplyr::select(name, ends_with("_significant")) %>%
  mutate(name = gsub("[.].*", "", name))

df_enrich <- NULL
message(gsub("_significant", "", contrast))
significant <- df[df[[contrast]],]
genes <- significant$name
enriched <- enrichR::enrichr(genes, databases)

# Tidy output
contrast_enrich <- NULL
for(database in databases) {
  temp <- enriched[database][[1]] %>%
    mutate(var = database)
  contrast_enrich <- rbind(contrast_enrich, temp)
}
contrast_enrich$contrast <- contrast
contrast_enrich$n <- length(genes)

# Background correction   # //
cat("Background correction... ")
contrast_enrich <- contrast_enrich %>%
  mutate(IN = as.numeric(gsub("/.*", "", Overlap)),
         OUT = n - IN) %>%
  dplyr::select(-n) %>%
  left_join(OUT, by = "Term") %>%
  mutate(log_odds = log2((IN * bg_OUT) / (OUT * bg_IN)))
cat("Done.")


pdf("1-DEPDistributio_plt0.05_logFCgt1.3.pdf")  # 1-DEPDistributio_Distance.pdf
plot_cor(dep,significant = TRUE,lower = 1,upper = -1,pal = rev("RdBu"))

plot_heatmap(dep,type = "centered",kmeans = TRUE,k = 6,col_limit = 3,
             show_row_names = TRUE,indicate = c("condition","replicate"))
plot_volcano(dep,contrast = contrastP,label_size = 2,add_names = TRUE,topgenes = 10)
dev.off()
row_data <- rowData(dep, use.names = FALSE)
con_vs_Model<-row_data[row_data$con_vs_Model_significant,]
ROS_vs_Model<-row_data[row_data$ROS_vs_Model_significant ,]

pdf("3-Dep_volcano.pdf",width=6,height=7)
Compared <- "Model_vs_con"
plot_volcano(dep,contrast = Compared,label_size = 2,add_names = F,topgenes = 10)
Compared <- "ROS_vs_con"
plot_volcano(dep,contrast = Compared,label_size = 2,add_names = F,topgenes = 10)
dev.off()

cols_diff<-grep("_diff", colnames(row_data))
row_data[,cols_diff]
row_data$Model_vs_conState=ifelse(row_data$Model_vs_con_significant,ifelse(row_data$Model_vs_con_diff > 0.25,"Up",ifelse(row_data$Model_vs_con_diff< -0.25,"Down","No")),"No")
row_data$ROS_vs_ModelState=ifelse(row_data$ROS_vs_Model_significant,ifelse(row_data$ROS_vs_Model_diff > 0.25,"Up",ifelse(row_data$ROS_vs_Model_diff< -0.25,"Down","No")),"No")
write.table(file="3-DeGprotein_modelVsotherAddupDown.xls",data.frame(row_data),sep="\t")

# venn

venn_list <- list()
vennnames<-c("Model_vs_con","ROS_vs_con")
venn_list[[1]]<-Model_vs_con$name
venn_list[[2]]<-ROS_vs_con$name
names(venn_list) <- vennnames
pdf('3-Degvenn_Overlapgenes.pdf', width = 6, height = 6)
venn(venn_list, zcolor = 'style',cex=6)
dev.off()
comongenes <- Reduce(intersect,venn_list)  # overlap genes
write.table(file="3-Model_vs_con_ROS_vs_con_overlapgenes.xls",comongenes,sep="\t",quote = FALSE)

row_data<-read.table(file="3-DeGprotein_Addup_down.xls",header = T, sep="\t",stringsAsFactors = F)
upvenn_list <- list()
downvenn_list <- list()
upvennnames<-c("Model_vs_conUp","ROS_vs_conUp")
downvennnames<-c("Model_vs_conDown","ROS_vs_conDown")
Model_vs_conUp<-subset(row_data, Model_vs_conState== "Up")
Model_vs_conDown<-subset(row_data, Model_vs_conState== "Down")
ROS_vs_conUp<-subset(row_data, ROS_vs_conState== "Up")
ROS_vs_conDown<-subset(row_data, ROS_vs_conState== "Down")
upvenn_list[[1]]<-Model_vs_conUp$name
downvenn_list[[1]]<-Model_vs_conDown$name
upvenn_list[[2]]<-ROS_vs_conUp$name
downvenn_list[[2]]<-ROS_vs_conDown$name
names(upvenn_list) <- upvennnames
names(downvenn_list) <- downvennnames
pdf('3-Degvenn_Overlapgenes_UpDown.pdf', width = 6, height = 6)
#venn(venn_list, zcolor = 'style',cex=6)
venn(upvenn_list, zcolor = 'style',cex=6)
venn(downvenn_list, zcolor = 'style',cex=6)
dev.off()

upcomongenes <- Reduce(intersect,upvenn_list)  # overlap genes
downcomongenes <- Reduce(intersect,downvenn_list)
comongenes<-cbind(upcomongenes, downcomongenes)
write.table(file="3-Model_vs_con_ROS_vs_con_Up_overlapgenes.xls",upcomongenes,sep="\t",quote = FALSE)
write.table(file="3-Model_vs_con_ROS_vs_con_Down_overlapgenes.xls",downcomongenes,sep="\t",quote = FALSE)
write.table(file="3-DeGprotein_Addup_down.xls",row_data,sep="\t")

# 富集分析
databases = c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021","GO_Biological_Process_2021","KEGG_2021_Human")
Compared <- "Model_vs_con_ROS_vs_con"
genes<-comongenes
enriched <- enrichR::enrichr(genes, databases)
contrast_enrich <- NULL
for(database in databases) {
  temp <- enriched[database][[1]] %>%
    mutate(var = database)
  contrast_enrich <- rbind(contrast_enrich, temp)
}
contrast_enrich$contrast <- Compared
contrast_enrich$n <- length(genes)
write.table(file=paste("3-","Model_vs_con_ROS_vs_con_EnrichmentOverlap.xls",sep=""),contrast_enrich,sep="\t", quote = FALSE)

#  Model_vs_conUp, Model_vs_conDown; ROS_vs_conUp, ROS_vs_conDown
ROS_modelUp<-subset(row_data, ROS_vs_ModelState =="Up")$name
ROS_modelDown<-subset(row_data, ROS_vs_ModelState =="Down")$name
databases = c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021","GO_Biological_Process_2021","KEGG_2021_Mouse")
Compared <- "ROS_modelDownDeg"
genes<-ROS_modelDown
enriched <- enrichR::enrichr(genes, databases)
contrast_enrich <- NULL
for(database in databases) {
  temp <- enriched[database][[1]] %>%
    mutate(var = database)
  contrast_enrich <- rbind(contrast_enrich, temp)
}
contrast_enrich$contrast <- Compared
contrast_enrich$n <- length(genes)
write.table(file=paste("3-",Compared,"Enrichment.xls",sep=""),contrast_enrich,sep="\t", quote = FALSE)

