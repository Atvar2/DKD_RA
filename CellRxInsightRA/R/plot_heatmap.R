# /* centered: 每个值减去行均值，中心化，contrast展现的是倍数关系。clustering_distance：聚类方式

library(dplyr)
library(ComplexHeatmap)
dep
indicate=c("condition","replicate")
type="centered"
k=6
row_font_size=6
col_font_size = 10
col_limit=1
row_data <- rowData(dep, use.names = FALSE)
col_data <- colData(dep) %>%
  as.data.frame()

# 判断数据的完整性
if(any(!c("label", "condition", "replicate") %in% colnames(col_data))) {
  stop(paste0("'label', 'condition' and/or 'replicate' columns are not present in '",
              deparse(substitute(dep)), "'"),
       call. = FALSE)
}
if(length(grep("_diff", colnames(row_data))) < 1) {
  stop(paste0("'[contrast]_diff' columns are not present in '",
              deparse(substitute(dep)),
              "'.\nRun test_diff() to obtain the required columns."),
       call. = FALSE)
}
if(!"significant" %in% colnames(row_data)) {
  stop(paste0("'significant' column is not present in '",
              deparse(substitute(dep)),
              "'.\nRun add_rejections() to obtain the required column."),
       call. = FALSE)
}

# 获取注释信息
get_annotation <- function(dep, indicate) {
  assertthat::assert_that(
    inherits(dep, "SummarizedExperiment"),
    is.character(indicate))

  # Check indicate columns
  col_data <- colData(dep) %>%
    as.data.frame()
  columns <- colnames(col_data)
  if(all(!indicate %in% columns)) {
    stop("'",
         paste0(indicate, collapse = "' and/or '"),
         "' column(s) is/are not present in ",
         deparse(substitute(dep)),
         ".\nValid columns are: '",
         paste(columns, collapse = "', '"),
         "'.",
         call. = FALSE)
  }
  if(any(!indicate %in% columns)) {
    indicate <- indicate[indicate %in% columns]
    warning("Only used the following indicate column(s): '",
            paste0(indicate, collapse = "', '"),
            "'")
  }

  # Get annotation
  anno <- dplyr::select(col_data, indicate)

  # Annotation color
  names <- colnames(anno)
  anno_col <- vector(mode="list", length=length(names))
  names(anno_col) <- names
  for(i in names) {
    var = anno[[i]] %>% unique() %>% sort()
    if(length(var) == 1)
      cols <- c("black")
    if(length(var) == 2)
      cols <- c("orangered", "cornflowerblue")
    if(length(var) < 7 & length(var) > 2)
      cols <- RColorBrewer::brewer.pal(length(var), "Pastel1")
    if(length(var) > 7)
      cols <- RColorBrewer::brewer.pal(length(var), "Set3")
    names(cols) <- var
    anno_col[[i]] <-  cols
  }

  # HeatmapAnnotation object Heatmap注释信息
  HeatmapAnnotation(df = anno,
                    col = anno_col,
                    show_annotation_name = TRUE)
}
ha1 <- get_annotation(dep, indicate) # !is.null(indicate) & type == "centered")

filtered <- dep[row_data$significant, ]   # filtered data, significant
#filtered<-dep
# Check for missing values
#  中心化， 数据减去均值
if(type == "centered") {
  rowData(filtered)$mean <- rowMeans(assay(filtered), na.rm = TRUE)
  df <- assay(filtered) - rowData(filtered, use.names = FALSE)$mean
}

set.seed(1)
df_kmeans <- kmeans(df, k)
if(type == "centered") {
  # Order the k-means clusters according to the maximum fold change
  # in all samples averaged over the proteins in the cluster
  order <- data.frame(df) %>%
    cbind(., cluster = df_kmeans$cluster) %>%
    mutate(row = apply(.[, seq_len(ncol(.) - 1)], 1, function(x) max(x))) %>%
    group_by(cluster) %>%
    summarize(index = sum(row)/n()) %>%
    arrange(desc(index)) %>%
    pull(cluster) %>%
    match(seq_len(k), .)
  df_kmeans$cluster <- order[df_kmeans$cluster]
}


# 是否聚类自动判断，即是否有一列
if(ncol(df) == 1) {
  col_clust = FALSE
} else {
  col_clust = TRUE
}
if(nrow(df) == 1) {
  row_clust = FALSE
} else {
  row_clust = TRUE
}

if(clustering_distance == "gower") {
  clustering_distance <- function(x) {
    dist <- cluster::daisy(x, metric = "gower")
    dist[is.na(dist)] <- max(dist, na.rm = TRUE)
    return(dist)
  }
}

# Legend info
legend <- ifelse(type == "contrast",
                 "log2 Fold change",
                 "log2 Centered intensity")

anno<-data.frame(Samples = c("Con1", "Con2", "Con3","Model1","Model2","Model3", "ROS1","ROS2","ROS3" ))
anno_col<-list(type = c("a" =  "red", "b" = "blue"))
ha1<-ComplexHeatmap::HeatmapAnnotation(df = anno,
                  #col = anno_col,
                  show_annotation_name = TRUE)

focusGenes <- read.table(file="./Focusinflammatory_Response.txt", header = F)
genes<-focusGenes$V1
mat<-test[genes,c("Con1", "Con2", "Con3","Model1","Model2","Model3", "ROS1","ROS2","ROS3")]
mat[is.na(mat)] <- 1
ht1 = ComplexHeatmap::Heatmap(mat,
              col = circlize::colorRamp2(
                seq(-col_limit, col_limit, (col_limit/5)),
                rev(RColorBrewer::brewer.pal(11, "RdBu"))),
              cluster_rows = T, cluster_columns = F,#width = ncol(mat)*unit(5, "mm"), height = nrow(mat)*unit(5, "mm"),
              row_names_side = "left",
              show_row_names = T,
              column_names_side = "top",
              heatmap_legend_param = list(color_bar = "continuous",
                                          legend_direction = "horizontal",
                                          legend_width = unit(5, "cm"),
                                          title_position = "lefttop"),
              name = legend,
              row_names_gp = gpar(fontsize = row_font_size),
              column_names_gp = gpar(fontsize = col_font_size),
              top_annotation = ha1)
pdf("3-DepheatmapInflammatory_Response.pdf",width = 6,height = 10)
ComplexHeatmap::draw(ht1, heatmap_legend_side = "top")
dev.off()
