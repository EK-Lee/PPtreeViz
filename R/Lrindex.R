#' Calculate Lr projection pursuit index value
#' @usage Lrindex(projclass,projdata,weight=TRUE,r=1)
#' @param projclass class information
#' @param projdata A training data  without class information
#' @param weight weight flag in LDA, PDA and Lr index
#' @param r r in Lr index
#' @references Lee, EK., Cook, D., Klinke, S., and Lumley, T.(2005) 
#' Projection Pursuit for exploratory supervised classification, 
#' Journal of Computational and Graphical statistics, 14(4):831-846.
#' @export
#' @keywords projection pursuit
#' @examples
#' data(iris)
#' Lrindex(iris[,5],as.matrix(iris[,1:4]))
#' 

Lrindex <- function(projclass, projdata, weight = TRUE, r = 1) {
    .Call('PPtreeViz_Lrindex', PACKAGE = 'PPtreeViz', projclass, projdata, weight, r)
}
