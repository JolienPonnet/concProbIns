\name{concProbContKMeans}

\alias{concProbContKMeans}

\title{Concordance probability based on k-means approximation for severity data.}

\description{Flexible calculation function of the concordance probability based on the k-means approximation for severity data.}

\usage{concProbContKMeans(inputDT, nu = 0, nClus, letsTime = TRUE, log = FALSE)}

\arguments{
  \item{inputDT}{data.table object consisting of the columns 'observed' and 'predicted'. This is an obligatory argument, without default value.}
  \item{nu}{Numeric vector of length 1 referring to the value that the observations have to differ at least from one another in order to be taken into account by the algorithm. Default value is 0.}
  \item{nClus}{Numeric integer indicating the number of clusters.}
  \item{log}{Logical vector of length 1 indicating which selection criterion is used for comparison. If TRUE, pairs are compaired when log(obs_i) - log(obs_j) < nu, if FALSE,  pairs are compaired when obs_i - obs_j < nu. Default value is FALSE.}
  \item{letsTime}{Logical vector of length 1 indicating whether or not the computing time of this function should be measured. Default value is TRUE.}
}

\value{
  \code{concProbContKMeans} returns a list that contains the following items:
    \item{concProbGlobal}{the concordance probability, obtained following the supplied values for the arguments}
    \item{time}{if \code{letsTime} is TRUE, the computing time of this function}
}

\references{Ponnet, J., Van Oirbeek, R., Verdonck, T. (2021). Concordance Probability for Insurance Pricing Models in a Big Data Setting.}

\author{Jolien Ponnet, Robin Van Oirbeek and Tim Verdonck}

\examples{
  library(concProbIns)
  library(data.table)

  nSamples <- 10000
  preds <- runif(n=nSamples, min=3500, max=5500)
  inputDT <- as.data.table(data.frame(observed = rgamma(nSamples, shape = preds/2000,
  scale = 2000), predicted = preds))
  concProbContKMeans(inputDT, nClus = 10)$concProbGlobal
  concProbContKMeans(inputDT, nu = 0.4, nClus = 20, log = TRUE)
  concProbContKMeans(inputDT, nu = 0.1, nClus = 30, letsTime = FALSE)
}
