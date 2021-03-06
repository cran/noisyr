% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_distance_matrices_transcript.R
\name{calculate_distance_matrices_transcript}
\alias{calculate_distance_matrices_transcript}
\title{Calcualte the distance matrices using the BAM files}
\usage{
calculate_distance_matrices_transcript(
  path.bams = ".",
  bams = NULL,
  path.gtf = list.files(".", pattern = "\\\\.g[tf]f$"),
  genes = NULL,
  expression.matrix = NULL,
  subsample.genes = FALSE,
  make.index = FALSE,
  uniqueOnly = TRUE,
  mapq.unique = 255,
  slack = 200,
  method = "correlation_pearson",
  save.image.every.1000 = FALSE,
  ncores = 1
)
}
\arguments{
\item{path.bams, bams}{either a path to the directory where the BAM files are
or a vector of paths to each individual file; if a path is specified,
it extracts all files that end in .bam; looks in the working directory by default}

\item{path.gtf}{the path to the gtf/gff annotation file;
if unspecified, looks for one in the working directory}

\item{genes}{a tibble of the genes extracted from the gtf file;
this is meant for speed if the output of cast_gtf_to_genes() is already generated}

\item{expression.matrix}{expression matrix; not necessary but is used to filter the
gtf to fewer entries and for subsampling if subsample.genes=TRUE;
if not provided, raw read counts are extracted from the BAM files}

\item{subsample.genes}{logical, whether to subsample low abundance genes to decrease
computational time; the first minimum of the distribution of abundances is calculated,
and genes lower than it are subsampled to match the number of genes higher than it;
a plot is generated to show that minimum}

\item{make.index}{whether a BAM index should be generated; if this is FALSE (the default)
and no index exists, the function will exit with an error; the index needs to have
the same name as each BAM file, but ending with .bam.bai}

\item{uniqueOnly}{whether only uniquely mapped reads should contribute to the profile;
default is TRUE}

\item{mapq.unique}{The values of the mapping quality field in the BAM file that corresponds
to uniquely mapped reads; by default, values of 255 are used as these correspond to
the most popular aligners, but an adjustment might be needed;
the mapq scores should be as follows: 255 for STAR, 60 for hisat2,
255 for bowtie in -k mode, 40 for bowtie2 default, 50 for tophat}

\item{slack}{slack needs to be >=readLength, adjust for efficiency; the default is 200,
as it is higher than most modern sequencing experiments}

\item{method}{one of the distance metrics to be used, defaults to pearson correlation
CURRENTLY ONLY CORRELATIONS!}

\item{save.image.every.1000}{whether to save a workspace image after every 1000 exons
are processed; default is FALSE}

\item{ncores}{Number of cores to use for parallel computation; needs doParallel installed;
defaults to sequential computation}
}
\value{
A list with three elements: the first element is the expression matrix,
        as supplied or calculated; the other two are the abundance and distance matrix;
        they have the same # of columns as the expression matrix,
        and as many rows as genes processed.
}
\description{
This function generates an average correlation/distance coefficient
for every exon present in the BAM files. This is done by calculating
the point-to-point correlation/distance of the distribution of reads
across the transcript of each exon and comparing it across samples.
}
\examples{
bams <- rep(system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE), 2)
genes <- data.frame("id" = 1:2,
                    "gene_id" = c("gene1", "gene2"),
                    "seqid" = c("seq1", "seq2"),
                    "start" = 1,
                    "end" = 1600)
obj <- calculate_distance_matrices_transcript(
  bams = bams,
  genes = genes,
  mapq.unique = 99
)
}
