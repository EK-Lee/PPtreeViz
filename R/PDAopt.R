#' Find the q-dim optimal projection using PDA index for classification
#' @usage PDAopt(origclass, origdata,q=1, weight=TRUE,lambda=0.1, ...) 
#' @param origdata A training data  without class information
#' @param origclass class information
#' @param q dimension of projection matrix
#' @param weight weight flag using in PDA index
#' @param lambda lambda in PDA index
#' @param ... arguments to be passed to methods
#' @return indexbest optimal value of PDA index value
#' @return projbest optimal q-dim projection
#' @references Lee, EK., Cook, D.(2010) 
#' A projection pursuit index for large p small n data, 
#' Statistics and Computing, 20:381-392.
#' @export
#' @keywords projection pursuit
#' @examples
#' data(iris)
#' PDA.proj.result <- PDAopt(iris[,5],iris[,1:4],weight=TRUE,q=2,lambda=0.1)
#' PDA.proj.result
PDAopt <- function(origclass,origdata,q=1,weight=TRUE,lambda=0.1,...) 
{ 
   origdata<-as.matrix(origdata)
   class.table<-table(origclass)
   g<-length(class.table)
   class.name<-names(class.table)
   p<-ncol(origdata)
   n<-nrow(origdata)
   mean.g<-matrix(apply(origdata,2,function(x) 
                  tapply(x,origclass,mean, na.rm=TRUE)),ncol=p)
   mean.all<-matrix(apply(origdata,2,mean),ncol=p)
                   
   B<-matrix(0,ncol=p,nrow=p)
   W<-matrix(0,ncol=p,nrow=p)
  
   for(i in 1:length(class.table)){
      sel.id<-which(origclass==class.name[i])
      temp.m1<-mean.g[i,]-mean.all
      temp.m2<-origdata[sel.id,]-
                 matrix(1,length(sel.id),ncol=1)%*%mean.g[i,,drop=FALSE]
      gn1<-ifelse(weight,length(sel.id),n/g)
      B <- B+ gn1*t(temp.m1)%*%temp.m1
      W <- W+ gn1*t(temp.m2)%*%temp.m2      
   }
   
   I<-diag(rep(1,p))
   W.t<-(1-lambda)*W+lambda*n*I
   opt<-eigen(MASS::ginv(W.t)%*%B)
  
   optVector<-matrix(as.numeric(opt$vectors[,1:q]),ncol=q)
   proj.data<-origdata%*%optVector
   optindex<-PDAindex(origclass,proj.data,weight,lambda)  
   return(list(indexbest=optindex,projbest=optVector))
}
  
    
  