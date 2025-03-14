plot_volcano <- function(dep, contrast, label_size = 3,topgenes=10,
  add_names = TRUE, adjusted = FALSE, plot = TRUE) {
  # Show error if inputs are not the required classes
  if(is.integer(label_size)) label_size <- as.numeric(label_size)
  assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
    is.character(contrast),
    length(contrast) == 1,
    is.numeric(label_size),
    length(label_size) == 1,
    is.logical(add_names),
    length(add_names) == 1,
    is.logical(adjusted),
    length(adjusted) == 1,
    is.logical(plot),
    length(plot) == 1)

  row_data <- rowData(dep, use.names = FALSE)

  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(row_data))) {
    stop(paste0("'name' and/or 'ID' columns are not present in '",
                deparse(substitute(dep)),
                "'.\nRun make_unique() to obtain required columns."),
         call. = FALSE)
  }
  if(length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '",
                deparse(substitute(dep)),
                "'.\nRun test_diff() to obtain the required columns."),
         call. = FALSE)
  }
  if(length(grep("_significant", colnames(row_data))) < 1) {
    stop(paste0("'[contrast]_significant' columns are not present in '",
                deparse(substitute(dep)),
                "'.\nRun add_rejections() to obtain the required columns."),
         call. = FALSE)
  }

  # Show error if an unvalid contrast is given
  if (length(grep(paste(contrast, "_diff", sep = ""),
                  colnames(row_data))) == 0) {
    valid_cntrsts <- row_data %>%
      data.frame() %>%
      select(ends_with("_diff")) %>%
      colnames(.) %>%
      gsub("_diff", "", .)
    valid_cntrsts_msg <- paste0("Valid contrasts are: '",
                                paste0(valid_cntrsts, collapse = "', '"),
                                "'")
    stop("Not a valid contrast, please run `plot_volcano()` with a valid contrast as argument\n",
         valid_cntrsts_msg,
         call. = FALSE)
  }

  # Generate a data.frame containing all info for the volcano plot
  diff <- grep(paste(contrast, "_diff", sep = ""),
               colnames(row_data))
  if(adjusted) {
    p_values <- grep(paste(contrast, "_p.adj", sep = ""),
                     colnames(row_data))
  } else {
    p_values <- grep(paste(contrast, "_p.val", sep = ""),
                     colnames(row_data))
  }
  signif <- grep(paste(contrast, "_significant", sep = ""),
                 colnames(row_data))
  df <- data.frame(x = row_data[, diff],
                   y = -log10(row_data[, p_values]),
                   significant = row_data[, signif],
                   name = row_data$name) %>%
    filter(!is.na(significant)) %>%
    arrange(significant)

  name1 <- gsub("_vs_.*", "", contrast)
  name2 <- gsub(".*_vs_", "", contrast)
  
  df$color <- ifelse(df$y>1.3 & abs(df$x)>= 0.26,ifelse(df$x > 0.26,'red','blue'),'gray')
  Topupgenes<-(df[df$x>0 & df$significant==TRUE,]  %>% arrange(desc(x)) %>% head(topgenes))$name
  Topdowngenes<-(df[df$x<0 & df$significant==TRUE,]  %>% arrange(x) %>% head(topgenes))$name
  labelgenes<-c(Topupgenes,Topdowngenes)
  df$lable<-df$name %in% labelgenes
  # Plot volcano with or without labels
  p <- ggplot(df, aes(x, y)) +
    geom_vline(xintercept = 0) +
    geom_point(aes(col = color)) +
    geom_text(data = data.frame(), aes(x = c(Inf, -Inf),
                                       y = c(-Inf, -Inf),
                                       hjust = c(1, 0),
                                       vjust = c(-1, -1),
                                       label = c(name1, name2),
                                       size = 5,
                                       fontface = "bold")) +
    labs(title = contrast,
      x = expression(log[2]~"Fold change")) +
    theme_DEP1() +
    theme(legend.position = "none") +
    scale_color_manual(values = c(red = "red", blue = "blue",gray="gray"))   # // 上调和下调基因先颜色显示
  if (add_names) {
    p <- p + ggrepel::geom_text_repel(data = filter(df, lable),
                                      aes(label = name), max.overlaps=10,
                                      size = label_size,                 #// label_size=3
                                      box.padding = unit(0.1, 'lines'),
                                      point.padding = unit(0.1, 'lines'),
                                      segment.size = 0.5)
  }
  if(adjusted) {
    p <- p + labs(y = expression(-log[10]~"Adjusted p-value"))
  } else {
    p <- p + labs(y = expression(-log[10]~"P-value"))
  }
  if(plot) {
    return(p)
  } else {
    df <- df %>%
      select(name, x, y, significant) %>%
      arrange(desc(x))
    colnames(df)[c(1,2,3)] <- c("protein", "log2_fold_change", "p_value_-log10")
    if(adjusted) {
      colnames(df)[3] <- "adjusted_p_value_-log10"
    }
    return(df)
  }
}
