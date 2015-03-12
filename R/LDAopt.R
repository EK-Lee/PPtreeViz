#' Find the 1D-optimal projection using LDA index for classification
#' @usage LDA.opt.1D(i.class, i.data, weight = TRUE, ...) 
#' @param i.class class information
#' @param i.data data without class information
#' @param weight weight flag using in LDA index
#' @return index.best optimal value of LDA index value
#' @return proj.best optimal 1D projection
#' @references Lee, EK., Cook, D., Klinke, S., and Lumley, T.(2005) 
#' Projection Pursuit for exploratory supervised classification, 
#' Journal of Computational and Graphical statistics, 14(4):831-846.
#' @export
#' @keywords projection pursuit
#' @examples
#' data(iris)
#' LDA.proj.result <- LDA.opt.1D(iris[,5],iris[,1:4])
#' LDA.proj.result
LDA.opt.1D<- function (i.class, i.data, weight, ...) 
{ 
   i.data<-as.matrix(i.data)
   class.table <- table(i.class)
   class.name<-names(class.table)
   p<-ncol(i.data)
   mean.g<-matrix(apply(i.data,2,function(x) 
                  tapply(x,i.class,mean, na.rm=TRUE)),ncol=p)
   mean.all<-matrix(apply(i.data,2,mean),ncol=p)
                   
   B<-matrix(0,ncol=p,nrow=p)
   W<-matrix(0,ncol=p,nrow=p)
  
   for(i in 1:length(class.table))
   {  sel.id<-which(i.class==class.name[i])
      temp.m1<-mean.g[i,]-mean.all
      temp.m2<-i.data[sel.id,]-
                matrix(1,length(sel.id),ncol=1)%*%mean.g[i,,drop=FALSE]
      if(weight){
         B <- B+ length(sel.id)*t(temp.m1)%*%temp.m1
         W <- W+ length(sel.id)*t(temp.m2)%*%temp.m2
      } else{
         B <- B+ t(temp.m1)%*%temp.m1
         W <- W+ t(temp.m2)%*%temp.m2
      }
        
   }
   opt<-eigen(MASS::ginv(W)%*%B)  

   optVector<-matrix(as.numeric(opt$vectors[,1]),ncol=1)
   proj.data<-i.data%*%optVector
   optindex<-LDAindex1(as.numeric(as.factor(as.character(i.class))),proj.data)  
   return(list(index.best=optindex,proj.best=optVector))
}
