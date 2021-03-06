#' A common problem faced by journal reviewers and authors is the question of
#' whether the results of a replication study are consistent with the original
#' published study. One solution to this problem is to examine the effect size
#' from the original study and generate the range of effect sizes that could
#' reasonably be obtained (due to random sampling) in a replication attempt
#' (i.e., calculate a prediction interval).This package has functions that
#' calculate the prediction interval for the correlation (i.e., r),
#' standardized mean difference (i.e., d-value), and mean.
#'
#'\tabular{ll}{
#'Package: \tab predictionInterval\cr
#'Type: \tab Package\cr
#'Version: \tab 1.0.0\cr
#'Date: \tab 2016-08-19\cr
#'License: \tab MIT License + file LICENSE\cr
#'}
#'
#'\code{\link{pi.r}} creates a prediction interval for a correlation (i.e., \emph{r} )\cr
#'\code{\link{pi.d}} creates a prediction interval for a standardized mean difference (i.e., \emph{d} )\cr
#'\code{\link{pi.m}} creates a prediction interval for a mean (i.e., \emph{M} )\cr
#'
#'\code{\link{pi.r.demo}} demonstrates PI capture percentage for a correlation (i.e., \emph{r} )\cr
#'\code{\link{pi.d.demo}} demonstrates PI capture percentage for a standardized mean difference (i.e., \emph{d} )\cr
#'\code{\link{pi.m.demo}} demonstrates PI capture percentage for a mean (i.e., \emph{M} )\cr
#'
#'@name predictionInterval-package
#'@aliases predictionInterval
#'@docType package
#'@title Prediction Interval Functions
#'@author
#'\tabular{ll}{
#'Author: \tab David J. Stanley \email{dstanley@@uoguelph.ca}\cr
#'Maintainer: \tab David J. Stanley \email{dstanley@@uoguelph.ca}
#'}
#'@references
#'Spence, J.R. & Stanley, D.J.(in prep). Prediction Interval: What to expect when you're expecting a replication. \cr\cr\cr
#'Also: \cr\cr
#'Cumming, G. & Maillardet, R. (2006). Confidence intervals and replication: where will the next mean fall? \emph{Psychological Methods, 11(3)}, 217-227. \cr\cr
#'Estes, W.K. (1997). On the communication of information by displays of standard error and confidence intervals. \emph{Psychonomic Bulleting & Review, 4(3)}, 330-341. \cr\cr
#'Zou, G.Y. (2007). Toward using a confidence intervals to compare correlations. \emph{Psychological Methods, 12(4)}, 399-413. \cr
#'@keywords package
#'@examples
#' pi.r(r=.35,n=100,rep.n=200)
#' pi.d(d=.65,n1=50,n2=50,rep.n1=100,rep.n2=100)
#' pi.m(M=2.53,SD=1.02,n=40,rep.n=80)
#' @import "ggplot2"
#' @importFrom "MASS" "mvrnorm"
#' @importFrom "MBESS" "ci.smd"
#' @importFrom "stats" "cor" "qt" "rnorm" "var"
NULL

