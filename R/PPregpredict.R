#' predict PPtree
#' 
#' Predict class for the test set and calculate prediction error 
#' After finding tree structure, predict class for the test set and calculate prediction error.
#' @usage PPreg.predict(PPtreeregOBJ,test.dataX,Rule=1,trueY=NULL,...)
#' @param PPtreeregOBJ PPtreereg object
#' @param test.dataX  the test dataset
#' @param Rule split rule 1: mean of two group means 
#'                        2: weighted mean of two group means - weight with group size
#'                        3: weighted mean of two group means - weight with group sd
#'                        4: weighted mean of two group means - weight with group se
#'                        5: mean of two group medians 
#'                        6: weighted mean of two group medians - weight with group size
#'                        7: weighted mean of two group median - weight with group IQR
#'                        8: weighted mean of two group median - weight with group IQR and size
#'                        9: cutoff that minimize error rates in each node                      
#' @param trueY true Y values of test dataset if available
#' @param ... arguments to be passed to methods
#' @return predict.Y predicted values
#' @return predict.MSE MSE of the predicted values
#' @export
#' @keywords tree
#' @examples
#' data(mtcars)
#' n <- nrow(mtcars)
#' tot <- c(1:n)
#' n.train <- round(n*0.8)
#' train <- sample(tot,n.train)
#' test <- tot[-train]
#' Tree.result <- PP.Tree.reg(mtcars[train,1],mtcars[train,-1],
#'                            final.rule=1,DEPTH=2)
#' PPreg.predict(Tree.result,mtcars[test,-1],1,mtcars[test,1])
PPreg.predict<-function(PPtreeregOBJ,test.dataX,Rule=1,trueY=NULL,...) {

   test.dataX<-as.matrix(test.dataX)
   trueY<-c(trueY)

   PP.Classification<-function(Tree.Struct,test.class.index,IOindex,
                                  test.class,id,rep){
      if(Tree.Struct[id,4]==0){
         i.class<-test.class
         i.class[i.class>0]<-1
         i.class<-1-i.class
         test.class<-test.class+IOindex*i.class*Tree.Struct[id,3]
         return(list(test.class=test.class,rep=rep))
      } else {  
         IOindexL<-IOindex*test.class.index[rep,]
         IOindexR<-IOindex*(1-test.class.index[rep,])
         rep<-rep+1
         a<-PP.Classification(Tree.Struct,test.class.index,IOindexL,
                                  test.class,Tree.Struct[id,2],rep)
         test.class<-a$test.class
         rep<-a$rep;
         a<-PP.Classification(Tree.Struct,test.class.index,IOindexR,
                              test.class,Tree.Struct[id,3],rep)
         test.class<-a$test.class
         rep<-a$rep
      }
      list(test.class=test.class,rep=rep)
   }
    
   PP.Class.index<-function(class.temp,test.class.index,test.data,
                            Tree.Struct,Alpha.Keep,C.Keep,id,Rule) {
      class.temp<-as.integer(class.temp)
      if(Tree.Struct[id,2]==0){
         return(list(test.class.index=test.class.index,class.temp=class.temp))
      } else {
         t.class<-class.temp 
         t.n<-length(t.class[t.class==0])
         t.index<-sort.list(t.class)
         if(t.n)
            t.index<-sort(t.index[-(1:t.n)])
         t.data<-test.data[t.index,]
         id.proj<-Tree.Struct[id,4]            
         proj.test<-as.matrix(test.data)%*%as.matrix(Alpha.Keep[id.proj,])
         proj.test<-as.double(proj.test)
         class.temp<-t(proj.test<C.Keep[id.proj,Rule]) 
         test.class.index<-rbind(test.class.index, class.temp)
         a<-PP.Class.index(class.temp,test.class.index,test.data,
                           Tree.Struct,Alpha.Keep,C.Keep,
                           Tree.Struct[id,2], Rule)
         test.class.index<-a$test.class.index
         a<-PP.Class.index(1-class.temp,test.class.index,test.data,
                             Tree.Struct, Alpha.Keep,C.Keep,
                             Tree.Struct[id,3],Rule)
         test.class.index<-a$test.class.index;
      }
      list(test.class.index=test.class.index,class.temp=class.temp)
   }
   Tree.result<-PPtreeregOBJ$Tree.result
   n<-nrow(test.dataX)
   class.temp<-rep(1, n)
   test.class.index<-NULL
   temp<-PP.Class.index(class.temp,test.class.index,test.dataX,
                        Tree.result$Tree.Struct,Tree.result$projbest.node,
                        Tree.result$splitCutoff.node,1,Rule)
   test.class<-rep(0,n)
   IOindex<-rep(1,n)
   rep<-1
   temp<-PP.Classification(Tree.result$Tree.Struct,temp$test.class.index,
                           IOindex,test.class,1,1)
   if(is.null(PPtreeregOBJ$coef.G)){
      predict.Y<-PPtreeregOBJ$mean.G[temp$test.class]
   } else{
      gt<-table(temp$test.class)
      predict.Y<-rep(0,length(temp$test.class))
      for(i in 1:length(gt)){
         sel.id<-which(temp$test.class==i)
         proj.data<-as.matrix(cbind(rep(1,nrow(test.dataX)),test.dataX))%*%
           matrix(PPtreeregOBJ$coef.G[i,])
         predict.Y[sel.id]<-proj.data[sel.id,1]
      }
   }
   if(!is.null(trueY)){
     predict.MSE<-mean(abs(trueY-predict.Y)^2)
   } else {
      predict.MSE<-NA
   }  
   list(predict.Y=predict.Y,predict.MSE=predict.MSE)
}

