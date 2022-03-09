\name{concProbGridLowHigh}

\alias{concProbGridLowHigh}

\title{Concordance probability based on marginal approximation for frequency data.}

\description{Flexible calculation function of the concordance probability based on the marginal approximation for frequency data. The splitpoints can be different for both considered groups.
}

\usage{concProbGridLowHigh(inputDT, lowCat, highCat, quantSplitsLow, quantSplitsHigh, gamma,
fine = FALSE, letsTime = TRUE)}

\arguments{
  \item{inputDT}{data.table object consisting of the columns 'observed' and 'predicted'. This is an obligatory argument, without default value.}
  \item{lowCat}{Numeric vector of length 1 or 2 indicating the lowest category of the event of interest (positive integer, containing 0) of the concordance probability, as defined by the 'observed' column of the 'inputDT' argument. If the numeric vector is of length 1, then just this corresponding value defines the lowest category. If the numeric vector is of length 2, then the corresponding values are the bounds of the inclusive interval. This is an obligatory argument, without default value.}
  \item{highCat}{Numeric vector of length 1 or 2 indicating the highest category of the event of interest (positive integer, containing 0) of the concordance probability, as defined by the 'observed' column of the 'inputDT' argument. If the numeric vector is of length 1, then just this corresponding value defines the highest category. If the numeric vector is of length 2, then the corresponding values are the bounds of the inclusive interval. This is an obligatory argument, without default value.}
  \item{quantSplitsLow}{Numeric vector representing the quantiles to define the boundary values of the grid for the lowest category.}
  \item{quantSplitsHigh}{Numeric vector representing the quantiles to define the boundary values of the grid for the highest category.}
  \item{gamma}{Maximum value that two observations their exposure can differ, such that their predictions are still comparable to determine the concordance probability.}
  \item{fine}{Logical vector of length 1. If TRUE, the for loop goes over the group with the biggest number of unique exposures (after clustering) or over the lowest category group when both groups have the same number of unique exposures (after clustering). If FALSE, the for loop goes over the group with the smallest number of unique exposures (after clustering) or over the highest category group when both groups have the same number of unique exposures (after clustering). Default value is FALSE.}
  \item{letsTime}{Logical vector of length 1 indicating whether or not the computing time of this function should be measured. Default value is TRUE.}
}

\value{
  \code{concProbGridLowHigh} returns a list that contains the following items:
    \item{concProbGlobal}{the concordance probability, obtained following the supplied values for the arguments}
  \item{uniqExp}{the unique exposures over which the for loop is done according to the value of \code{fine}}
  \item{concProbs}{the local concordance probabilities, one for each exposure in \code{uniqExp}}
  \item{compPairs}{the number of comparable pairs, one for each exposure in \code{uniqExp}}
  \item{time}{if \code{letsTime} is TRUE, the computing time of this function}
}

\references{Ponnet, J., Van Oirbeek, R., Verdonck, T. (2021). Concordance Probability for Insurance Pricing Models in a Big Data Setting.}

\author{Jolien Ponnet, Robin Van Oirbeek and Tim Verdonck}

\examples{
  library(concProbIns)
  library(data.table)

  nSamples <- 10000
  obsValues <- rpois(nSamples, lambda = 1)
  inputDT <- as.data.table(data.frame(exposure = runif(nSamples), observed = obsValues,
  predicted = abs(obsValues + rnorm(nSamples))))
  concProbGridLowHigh(inputDT, 0, 1, quantSplitsLow = seq(0.02, 0.98, 0.02),
  quantSplitsHigh = seq(0.02, 0.98, 0.02), gamma = 0.05, letsTime = FALSE)$concProbGlobal
  concProbGridLowHigh(inputDT, 0, 1, quantSplitsLow = seq(0.1, 0.9, 0.1),
  quantSplitsHigh = seq(0.01, 0.99, 0.01), gamma = 0.1, fine = TRUE)
  concProbGridLowHigh(inputDT, 1, 2, quantSplitsLow = seq(0.2, 0.8, 0.2),
  quantSplitsHigh = seq(0.1, 0.9, 0.1), gamma = 0.1, fine = TRUE)$compPairs
  concProbGridLowHigh(inputDT, c(0,1), c(2,4), quantSplitsLow = seq(0.01, 0.99, 0.01),
  quantSplitsHigh = seq(0.02, 0.98, 0.02), gamma = 0.15, letsTime = FALSE)$concProbs
}
