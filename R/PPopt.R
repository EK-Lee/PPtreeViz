#' PP optimization using various PP indices
#' 
#' Find the q-dim optimal projection using various projectin pursuit indices with class information
#' @usage PPopt(origclass,origdata,q=1,PPmethod="LDA",weight=TRUE,r=1,lambda=0.1,
#'              energy=0,cooling=0.999,TOL=0.0001,maxiter = 50000)
#' @param origclass class information
#' @param origdata A training data  without class information
#' @param q dimension of projection matrix
#' @param PPmethod method for projection pursuit; "LDA", "PDA", "Lr", "GINI", and "ENTROPY"
#' @param weight weight flag in LDA, PDA and Lr index
#' @param r r in Lr index
#' @param lambda lambda in PDA index
#' @param energy energy parameter
#' @param cooling cooling parameter
#' @param TOL tolerance
#' @param maxiter number of maximum iteration
#' @return indexbest optimal value of PP index value
#' @return projbest optimal q-dim projection
#' @return origclass class information
#' @return origdata A training data  without class information
#' @references Lee, EK., Cook, D., Klinke, S., and Lumley, T.(2005) 
#' Projection Pursuit for exploratory supervised classification, 
#' Journal of Computational and Graphical statistics, 14(4):831-846.
#' @export
#' @keywords projection pursuit
#' @examples
#' data(iris)
#' PP.proj.result <- PPopt(iris[,5],as.matrix(iris[,1:4]))
#' PP.proj.result

PPopt <- function(origclass, origdata, q=1, PPmethod = "LDA", weight = TRUE, r = 1, lambda = 0.1, energy = 0, cooling = 0.999, TOL = 0.0001, maxiter = 50000) {
    .Call('PPtreeViz_PPopt', PACKAGE = 'PPtreeViz', origclass, origdata, q, PPmethod, weight, r, lambda, energy, cooling, TOL, maxiter)
}
