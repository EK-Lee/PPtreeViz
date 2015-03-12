#' Finding the 1D-optimal projection with PDA index
#' 
#' Find the 1D-optimal projection using PDA index for classification
#' @usage print.PPtree(PPtreeOBJ,coef.print=FALSE,cutoff.print=FALSE)
#' @param PPtreeOBJ PPtree object
#' @param coef.print print projection coefficient in each node when TRUE
#' @param cutoff.print print cutoff values in each node when TRUE
#' @references Lee, E.-K., Cook, D. (2009) A Projection Pursuit Index for Large p Small n Data, 
#' Statistics and Computing,  
#' \url{http://www.springerlink.com/content/g47n0n342761838m/#?p=d2ff5a7b69eb45ef8abf7ef3aba69557&pi=3}.
#' @export
#' @keywords tree
#' @examples
#' data(iris)
#' Tree.result <- LDA.Tree.viz(iris[,5],iris[,1:4])
#' Tree.result
#' print(Tree.result,coef.print=TRUE,cutoff.print=TRUE)
print.PPtree<-function(PPtreeOBJ,coef.print=FALSE,cutoff.print=FALSE,verbose=TRUE){
   TS<-PPtreeOBJ$Tree.Struct
   Alpha<-PPtreeOBJ$Alpha.Keep
   cut.off<-PPtreeOBJ$C.Keep
   gName<-names(table(PPtreeOBJ$orig.class))

   pastemake<-function(k,arg,sep.arg="") {
      temp<-""
      for(i in 1:k)
         temp<-paste(temp,arg,sep=sep.arg)
      return(temp)
   }

   TreePrint<-"1) root"
   i<-1;  flag.L<-rep(FALSE,nrow(TS))
   keep.track<-1
   depth.track<-0
   depth<-0
   while(sum(flag.L)!=nrow(TS)){
      if(!flag.L[i]) {                    
         if(TS[i,2] == 0) {
            flag.L[i]<-TRUE
            n.temp<-length(TreePrint)
            tempp<-strsplit(TreePrint[n.temp],") ")[[1]]
            temp.L<-paste(tempp[1],")*",tempp[2],sep="")
            temp.L<- paste(temp.L, gName[TS[i,3]],sep="  ->  ")
            TreePrint<-TreePrint[-n.temp]
            id.l<-length(keep.track)-1
            i <- keep.track[id.l]
            depth <- depth -1
         } else if(!flag.L[TS[i,2]]) {
            depth <- depth +1
            emptyspace<-pastemake(depth,"   ")
            temp.L<- paste(emptyspace,TS[i,2],")  proj",TS[i,4],"*X < cut",TS[i,4],sep="")
            i<-TS[TS[i,2],1]   
         } else {
            depth <- depth +1
            emptyspace<-pastemake(depth,"   ")          
            temp.L<- paste(emptyspace,TS[i,3],")  proj",TS[i,4],"*X >= cut",TS[i,4],sep="")
            flag.L[i]<-TRUE
            i<-TS[TS[i,3],1]
         } 
         keep.track <- c(keep.track,i)
         depth.track<-c(depth.track,depth)
         TreePrint<-c(TreePrint,temp.L)
      } else {
         id.l<-id.l-1
         i<-keep.track[id.l]
         depth <- depth.track[id.l]
      }
   }
   colnames(Alpha)<-colnames(PPtreeOBJ$orig.data)
   rownames(Alpha)<-paste("proj",1:nrow(Alpha),sep="")

   colnames(cut.off)<-paste("Rule",1:ncol(cut.off),sep="")
   rownames(cut.off)<-paste("cut",1:nrow(cut.off),sep="")

   TreePrint.output<-paste("=====================================",
                           "\nProjection Pursuit Tree result",
                           "\n=====================================\n")
   for(i in 1:length(TreePrint))
      TreePrint.output<-paste(TreePrint.output,TreePrint[i],sep="\n")
   TreePrint.output<-paste(TreePrint.output,"\n",sep="")
   sample.data.X<-PPtreeOBJ$orig.data
   sample.data.class<-PPtreeOBJ$orig.class

   error.rate<-matrix(c(PP.classify.viz(sample.data.X, sample.data.class, 
                                 PPtreeOBJ,Rule=1)$predict.error,
                 PP.classify.viz(sample.data.X, sample.data.class, 
                                 PPtreeOBJ,Rule=2)$predict.error,
                 PP.classify.viz(sample.data.X, sample.data.class, 
                                 PPtreeOBJ,Rule=3)$predict.error,
                 PP.classify.viz(sample.data.X, sample.data.class, 
                                 PPtreeOBJ,Rule=4)$predict.error,
                 PP.classify.viz(sample.data.X, sample.data.class, 
                                 PPtreeOBJ,Rule=5)$predict.error,
                 PP.classify.viz(sample.data.X, sample.data.class, 
                                 PPtreeOBJ,Rule=6)$predict.error)/nrow(sample.data.X),nrow=1)
                 
   colnames(error.rate)<-colnames(cut.off)
   rownames(error.rate)<-c("error.rate")
   if(verbose){
      cat(TreePrint.output)
      if(coef.print){
         cat(#"\n-------------------------------------",
             "\nProjection Coefficient in each node",
             "\n-------------------------------------\n")
         print(Alpha)
      }
      if(cutoff.print){
         cat(#"\n-------------------------------------",
             "\nCutoff values of each node",
             "\n-------------------------------------\n")
         print(cut.off)
      }    
      cat(#"\n-------------------------------------",
          "\nError rates of various cutoff values",
          "\n-------------------------------------\n")
      print(error.rate)
   }    
      return(invisible(TreePrint)) 
   }
