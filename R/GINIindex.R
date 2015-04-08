#' GINI PP index
#' 
#' Calculate GINI projection pursuit index value
#' @usage GINIindex1D(projclass,projdata)
#' @param projclass class information
#' @param projdata projected data  
#' @references Lee, EK., Cook, D., Klinke, S., and Lumley, T.(2005) 
#' Projection Pursuit for exploratory supervised classification, 
#' Journal of Computational and Graphical statistics, 14(4):831-846.
#' @export
#' @keywords projection pursuit
#' @examples
#' data(iris)
#' GINIindex1D(iris[,5],iris[,1])

GINIindex1D <- function(projclass, projdata) {
    .Call('PPtreeViz_GINIindex1D', PACKAGE = 'PPtreeViz', projclass, projdata)
}
