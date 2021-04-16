## ----options, include = FALSE-------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----bioc_install, eval = FALSE-----------------------------------------------
#  packages.bioc <- c("preprocessCore",
#                     "IRanges",
#                     "GenomicRanges",
#                     "Rsamtools")
#  new.packages.bioc <- packages.bioc[!(packages.bioc %in% installed.packages()[,"Package"])]
#  if(length(new.packages.bioc)){
#    if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#    BiocManager::install(new.packages.bioc)
#  }

## ----cran_install, eval = FALSE-----------------------------------------------
#  install.packages("noisyr")

## ----github_install, eval = FALSE---------------------------------------------
#  packages.cran <- c("utils",
#                     "grDevices",
#                     "tibble",
#                     "dplyr",
#                     "magrittr",
#                     "ggplot2",
#                     "philentropy",
#                     "doParallel",
#                     "foreach")
#  new.packages.cran <- packages.cran[!(packages.cran %in% installed.packages()[,"Package"])]
#  if(length(new.packages.cran))
#    install.packages(new.packages.cran)
#  
#  if (!requireNamespace("devtools", quietly = TRUE))
#    install.packages("devtools")
#  devtools::install_github("Core-Bioinformatics/noisyR")

## ----setup--------------------------------------------------------------------
library(noisyr)

## ----read---------------------------------------------------------------------
counts.in <- system.file("extdata", "counts_raw.csv", package = "noisyr")
df <- read.csv(counts.in, row.names = 1)
str(df)
head(df)

## ----convert------------------------------------------------------------------
expression.matrix <- noisyr::cast_matrix_to_numeric(df)

## ----noisyr_counts------------------------------------------------------------
expression.matrix.denoised.standard <- noisyr::noisyr(
  approach.for.similarity.calculation = "counts", 
  expression.matrix = df
)

## ----rm_noise_out-------------------------------------------------------------
head(expression.matrix.denoised.standard)
apply(expression.matrix.denoised.standard, 2, min)

## ----runCM--------------------------------------------------------------------
expression.summary <- noisyr::calculate_expression_similarity_counts(
  expression.matrix = expression.matrix, 
  similarity.measure = "correlation_pearson"
)
str(expression.summary)

## ----dist_metrics-------------------------------------------------------------
noisyr::get_methods_correlation_distance()

## ----window_opt---------------------------------------------------------------
noisyr::optimise_window_length(
  expression.matrix = expression.matrix,
  similarity.measure = "correlation_pearson"
)

## ----simple_plot, warning = FALSE---------------------------------------------
plotlist <- noisyr::plot_expression_similarity(
  expression.summary = expression.summary)
plotlist[[1]]

## ----combined_plot, warning = FALSE-------------------------------------------
plotdf.line <- tibble::tibble()
for(i in 1:4){
  lineid <- i * 2 - 1
  plotdf.line <- rbind(
    plotdf.line, 
    dplyr::mutate(plotlist[[lineid]]$data,
                  Sample=colnames(expression.matrix)[i]))
}

ggplot2::ggplot(plotdf.line) +
    ggplot2::theme_minimal() + 
    ggplot2::geom_line(ggplot2::aes(x=x, y=y, colour=Sample)) +
    ggplot2::geom_smooth(ggplot2::aes(x,y,colour=Sample), method="loess",
                         formula= y ~ x, span=0.1) +
    ggplot2::ylim(0:1) +
    ggplot2::xlab("log2(expression)") +
    ggplot2::ylab("Pearson correlation") +
    ggplot2::geom_hline(yintercept=0.25, color="black")

## ----calc_thr-----------------------------------------------------------------
noise.thresholds <- noisyr::calculate_noise_threshold(expression = expression.summary)
noise.thresholds

## ----get_thr_methods----------------------------------------------------------
noisyr::get_methods_calculate_noise_threshold()

## ----calc_thr_range-----------------------------------------------------------
similarity.threshold.sequence <- seq(0.2, 0.3, by=0.01)
stats.table <- noisyr::calculate_noise_threshold_method_statistics(
  expression = expression.summary,
  similarity.threshold.sequence = similarity.threshold.sequence
)
row.min.coef.var <- which.min(stats.table$noise.threshold.coefficient.of.variation)
# adjust column names for printing
colnames(stats.table) <- c("approach", "method", "corr.thr", "min", "mean", "coef.var", "max", "all")
stats.table[row.min.coef.var, 1:7]
dplyr::filter(stats.table, round(corr.thr, 2) == 0.21)[, 1:7]
dplyr::filter(stats.table, method == "loess10_smoothing")[, 1:7]

## ----calc_thr_adjusted--------------------------------------------------------
noise.thresholds <- noisyr::calculate_noise_threshold(
  expression = expression.summary,
  similarity.threshold = 0.21,
  method.chosen = "Line_plot-loess10_smoothing"
)
noise.thresholds

## ----rm_noise_matrix----------------------------------------------------------
expression.matrix.denoised <- noisyr::remove_noise_from_matrix(
  expression.matrix = expression.matrix,
  noise.thresholds = noise.thresholds)
str(expression.matrix.denoised)

## ----rm_noise_param-----------------------------------------------------------
expression.matrix.denoised.fixed <- noisyr::remove_noise_from_matrix(
  expression.matrix = expression.matrix, 
  noise.thresholds = mean(noise.thresholds))
nrow(expression.matrix.denoised); nrow(expression.matrix.denoised.fixed)

## ----pipeline_comparison------------------------------------------------------
expression.matrix.denoised.full.pipeline <- noisyr::noisyr(
  approach.for.similarity.calculation = "counts", 
  expression.matrix = expression.matrix,
  similarity.measure = "correlation_pearson",
  optimise.window.length.logical = FALSE,
  similarity.threshold.sequence = seq(0.2, 0.3, by=0.01),
  method.chosen.sequence = get_methods_calculate_noise_threshold()
)
identical(expression.matrix.denoised,
          expression.matrix.denoised.full.pipeline)

## ----edgeR, eval = FALSE------------------------------------------------------
#  DE_edgeR = function(expression.matrix, metadata){
#    # load or install edgeR
#    if(!require(edgeR)){
#      if (!requireNamespace("BiocManager", quietly = TRUE))
#        install.packages("BiocManager")
#      BiocManager::install("edgeR")
#    }
#  
#    # create metadata
#    metadata <- data.frame(id = colnames(expression.matrix),
#                           timepoint = c("0h", "0h", "12h", "12h"))
#  
#    # quantile normalise
#    expression.matrix.normalised <-
#      preprocessCore::normalize.quantiles(expression.matrix)
#    rownames(expression.matrix.normalised) <- base::rownames(expression.matrix)
#    colnames(expression.matrix.normalised) <- base::colnames(expression.matrix)
#  
#    # process using edgeR
#    expression.matrix.for.de <- round(expression.matrix.normalised)
#    expression.matrix.for.de <-
#      expression.matrix.for.de[apply(expression.matrix.for.de, 1, sum) > 0, ]
#    design <- model.matrix(~ 0 + metadata$timepoint)
#    edger <- DGEList(counts = expression.matrix.for.de)
#    edger <- estimateDisp(edger, design)
#    edger.fit <- glmFit(edger, design)
#    edger.lrt <- glmLRT(edger.fit, contrast=c(-1, 1))
#  
#    # extract results
#    res <- topTags(edger.lrt, n = Inf)$table
#    res$DE <- res$FDR < 0.05 & abs(res$logFC) > 1
#  
#    # make volcano plot
#    print(
#      ggplot2::ggplot(res) +
#        ggplot2::theme_minimal() +
#        ggplot2::geom_point(ggplot2::aes(x=logFC,
#                                         y=-log10(FDR),
#                                         colour=DE),
#                            show.legend=FALSE) +
#        ggplot2::scale_color_manual(values=c("black", "red")) +
#        ggplot2::lims(x=c(-12, 12), y=c(0, 120))
#    )
#  
#    return(res)
#  
#  }
#  
#  results.raw <- DE_edgeR(expression.matrix)
#  results.denoised <- DE_edgeR(expression.matrix.denoised)

## ----volcanos, echo = FALSE, out.width = "100%"-------------------------------
knitr::include_graphics("volcanos.png") 

## ----session------------------------------------------------------------------
sessionInfo()

