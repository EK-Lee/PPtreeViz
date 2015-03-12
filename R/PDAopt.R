#' Find the 1D-optimal projection using PDA index for classification
#' @usage PDA.opt.1D(i.class, i.data, weight=TRUE,lambda, ...) 
#' @param i.data A training data  without class information
#' @param i.class class information
#' @param weight weight flag using in PDA index
#' @param lambda lambda in PDA index
#' @return index.best optimal value of LDA index value
#' @return proj.best optimal 1D projection
#' @references Lee, EK., Cook, D.(2010) 
#' A projection pursuit index for large p small n data, 
#' Statistics and Computing, 20:381-392.
#' @export
#' @keywords projection pursuit
#' @examples
#' data(iris)
#' PDA.proj.result <- PDA.opt.1D(iris[,5],iris[,1:4],lambda=0.1)
#' PDA.proj.result
PDA.opt.1D<- function(i.class,i.data,lambda,weight, ...) 
{ 
   i.data<-as.matrix(i.data)
   class.table<-table(i.class)
   class.name<-names(class.table)
   p<-ncol(i.data)
   n<-nrow(i.data)
   mean.g<-matrix(apply(i.data,2,function(x) 
                  tapply(x,i.class,mean, na.rm=TRUE)),ncol=p)
   mean.all<-matrix(apply(i.data,2,mean),ncol=p)
                   
   B<-matrix(0,ncol=p,nrow=p)
   W<-matrix(0,ncol=p,nrow=p)
  
   for(i in 1:length(class.table)){
      sel.id<-which(i.class==class.name[i])
      temp.m1<-mean.g[i,]-mean.all
      temp.m2<-i.data[sel.id,]-
                 matrix(1,length(sel.id),ncol=1)%*%mean.g[i,,drop=FALSE]
      if(weight){
         B<-B+length(sel.id)*t(temp.m1)%*%temp.m1
         W<-W+length(sel.id)*t(temp.m2)%*%temp.m2
      } else{
         B<-B+t(temp.m1)%*%temp.m1
         W<-W+t(temp.m2)%*%temp.m2
      }
   }
   I<-diag(rep(1,p))
   W.t<-(1-lambda)*W+lambda*n*I
   opt<-eigen(solve(W.t)%*%B)
  
   optVector<-opt$vectors[,1]
   proj.data<-i.data%*%optVector
   optindex<-PDAindex1(as.numeric(as.factor(as.character(i.class))),proj.data,lambda)  
   return(list(index.best=optindex,proj.best=optVector))
}


  
    
  