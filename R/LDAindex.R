#' LDA PP index
#' 
#' Calculate LDA projection pursuit index value
#' @usage LDAindex(projclass,projdata,weight=TRUE)
#' @param projclass class information
#' @param projdata projected data  
#' @param weight weight flag in LDA
#' @references Lee, EK., Cook, D., Klinke, S., and Lumley, T.(2005) 
#' Projection Pursuit for exploratory supervised classification, 
#' Journal of Computational and Graphical statistics, 14(4):831-846.
#' @export
#' @keywords projection pursuit
#' @examples
#' data(iris)
#' LDAindex(iris[,5],as.matrix(iris[,1:4]))
#' 

LDAindex <- function(projclass, projdata, weight = TRUE) {
    .Call('PPtreeViz_LDAindex', PACKAGE = 'PPtreeViz', projclass, projdata, weight)
}

