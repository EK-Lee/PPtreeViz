#' Find the q-dim optimal projection using LDA projectin pursuit index
#' @usage LDAopt(origclass, origdata,  q=1, weight = TRUE,...) 
#' @param origclass class information
#' @param origdata data without class information
#' @param q dimension of projection matrix
#' @param weight weight flag using in LDA index
#' @param ... arguments to be passed to methods
#' @return indexbest optimal value of LDA index value
#' @return projbest optimal q-dim projection
#' @references Lee, EK., Cook, D., Klinke, S., and Lumley, T.(2005) 
#' Projection Pursuit for exploratory supervised classification, 
#' Journal of Computational and Graphical statistics, 14(4):831-846.
#' @export
#' @keywords projection pursuit
#' @examples
#' data(iris)
#' LDA.proj.result <- LDAopt(iris[,5],iris[,1:4])
#' LDA.proj.result
LDAopt<- function (origclass, origdata, q=1, weight=TRUE, ...) 
{ 
   origdata<-as.matrix(origdata)
   class.table <- table(origclass)
   class.name<-names(class.table)
   p<-ncol(origdata)
   mean.g<-matrix(apply(origdata,2,function(x) 
                  tapply(x,origclass,mean, na.rm=TRUE)),ncol=p)
   mean.all<-matrix(apply(origdata,2,mean),ncol=p)
                   
   B<-matrix(0,ncol=p,nrow=p)
   W<-matrix(0,ncol=p,nrow=p)
   g<-length(class.table)
   n<-length(origclass)
   
   for(i in 1:g)
   {  sel.id<-which(origclass==class.name[i])
      temp.m1<-mean.g[i,]-mean.all
      temp.m2<-origdata[sel.id,]-
                matrix(1,length(sel.id),ncol=1)%*%mean.g[i,,drop=FALSE]
      
      gn1<-ifelse(weight,length(sel.id),n/g)
      B <- B+ gn1*t(temp.m1)%*%temp.m1
      W <- W+ gn1*t(temp.m2)%*%temp.m2
   }
   
   opt<-eigen(MASS::ginv(W)%*%B)  
   optVector<-matrix(as.numeric(opt$vectors[,1:q]),ncol=q)
   proj.data<-origdata%*%optVector
   optindex<-LDAindex(origclass,proj.data,weight)  
   
   return(list(indexbest=optindex,projbest=optVector))
}
