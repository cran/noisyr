% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_noise_matrix.R
\name{remove_noise_matrix}
\alias{remove_noise_matrix}
\title{Function to remove the noisy reads from the expression matrix}
\usage{
remove_noise_matrix(
  expression.matrix,
  abn.thresh,
  add.thresh = TRUE,
  average.thresh = TRUE,
  remove.noisy.features = TRUE,
  export.csv = NULL
)
}
\arguments{
\item{expression.matrix}{the expression matrix}

\item{abn.thresh}{a vector of abundance thresholds; must be the same length
as the number of columns of the expression matrix}

\item{add.thresh}{whether to add the noise threshold to all values in the expression matrix
(default), or set entries below the threshold to the threshold}

\item{average.thresh}{if TRUE (default), uses tthe average of the vector of thresholds across all samples;
if FALSE, uses the thresholds as supplied}

\item{remove.noisy.features}{logical, whether rows of the expression matrix that are
fully under the noise threshold should be removed (default TRUE)}

\item{export.csv}{option to write the matrix into a csv after the noise removal;
should be NULL or the name of the output file.}
}
\value{
Returns a matrix of the same dims as the expression matrix,  with the noise removed.
This matrix has no entries remaining below the noise threshold.
}
\description{
This function is used to remove the noisy reads from the expression matrix.
It uses as input a vector of abundance thresholds;
all entries below the noise threshold are replaced with the noise threshold.
}
\examples{
remove_noise_matrix(
    expression.matrix = matrix(1:100, ncol=5),
    abn.thresh=c(5,30,45,62,83))
}
