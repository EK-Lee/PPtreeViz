#' PDA PP index 
#' 
#' Calculate PDA projection pursuit index value 
#' @usage PDAindex(projclass,projdata,weight=TRUE,lambda=0.1)
#' @param projclass class information
#' @param projdata projected data  
#' @param weight weight flag in PDA
#' @param lambda lambda in PDA index
#' @references Lee, EK., Cook, D., Klinke, S., and Lumley, T.(2005) 
#' Projection Pursuit for exploratory supervised classification, 
#' Journal of Computational and Graphical statistics, 14(4):831-846.
#' @export
#' @keywords projection pursuit
#' @examples
#' data(iris)
#' PDAindex(iris[,5],as.matrix(iris[,1:4]))
#' 

PDAindex <- function(projclass, projdata, weight = TRUE, lambda = 0.1) {
    .Call('PPtreeViz_PDAindex', PACKAGE = 'PPtreeViz', projclass, projdata, weight, lambda)
}
