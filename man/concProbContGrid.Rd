\name{concProbContGrid}

\alias{concProbContGrid}

\title{Concordance probability based on marginal approximation for severity data.}

\description{Flexible calculation function of the concordance probability based on the marginal approximation for severity data.}

\usage{concProbContGrid(inputDT, nu = 0, quantSplits, direction = 'left', letsTime = TRUE,
log = FALSE)}

\arguments{
  \item{inputDT}{data.table object consisting of the columns 'observed' and 'predicted'. This is an obligatory argument, without default value.}
  \item{nu}{Numeric vector of length 1 referring to the value that the observations have to differ at least from one another in order to be taken into account by the algorithm. Default value is 0.}
  \item{quantSplits}{Numeric vector representing the quantiles to define the boundary values of the grid.}
  \item{direction}{Character vector of length 1, with as possible values 'left', 'right', 'up', 'down' or 'all'. It is the direction that is used to make sure that each region-to-region comparison contributes once only. Default value is 'left'.}
  \item{log}{Logical vector of length 1 indicating which selection criterion is used for comparison. If TRUE, pairs are compaired when log(obs_i) - log(obs_j) < nu, if FALSE,  pairs are compaired when obs_i - obs_j < nu. Default value is FALSE.}
  \item{letsTime}{Logical vector of length 1 indicating whether or not the computing time of this function should be measured. Default value is TRUE.}
}

\value{
  \code{concProbContGrid} returns a list that contains the following items:
    \item{concProbGlobal}{the concordance probability, obtained following the supplied values for the arguments}
  \item{time}{if \code{letsTime} is TRUE, the computing time of this function}
  \item{concProbMat}{list that consists of two matrices; \code{concMat}, a matrix representing the number of concordant pairs and \code{discMat}, a matrix representing the number of discordant pairs}
  \item{splits}{vector with the values of the splits that are used to construct the grid}
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
  concProbContGrid(inputDT, quantSplits = seq(0.1,0.9,0.1))$concProbGlobal
  concProbContGrid(inputDT, nu = 0.4, quantSplits = seq(0.1,0.9,0.01), log = TRUE)
  concProbContGrid(inputDT, nu = 0.2, quantSplits = seq(0.1,0.9,0.1))$concProbMat
  concProbContGrid(inputDT, nu = 0.1, quantSplits = seq(0.1,0.9,0.1), direction ='up',
  letsTime = FALSE)
}
