#' ENTROPY PP index 
#' 
#' Calculate ENTROPY projection pursuit index value
#' @usage ENTROPYindex1D(projclass, projdata)
#' @param projclass class information
#' @param projdata projected data  
#' @export
#' @keywords projection pursuit
#' @examples
#' data(iris)
#' ENTROPYindex1D(iris[,5],iris[,1])
ENTROPYindex1D <- function(projclass, projdata) {
    .Call('PPtreeViz_ENTROPYindex1D', PACKAGE = 'PPtreeViz', projclass, projdata)
}
