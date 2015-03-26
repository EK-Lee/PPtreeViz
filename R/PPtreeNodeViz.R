#' Finding PP tree structure using LDA index
#' 
#' Find tree structure using linear discriminant in each split.
#' @usage PPtreeNode.Viz(PPtreeOBJ,node.id,Rule)
#' @param PPtreeOBJ pptree object
#' @param node.id node ID
#' @param Rule cutoff rule
#' @references Lee, YD, Cook, D., Park JW, and Lee, EK(2013) 
#' PPtree: Projection pursuit classification tree, 
#' Electronic Journal of Statistics, 7:1369-1386.
#' @export
#' @keywords tree
#' @examples
#' data(iris)
#' Tree.result <- LDA.Tree(iris[,5],iris[,1:4])
#' Tree.result
#' PPtreeNode.Viz(Tree.result,1,1)
PPtreeNode.Viz<-function(PPtreeOBJ,node.id,Rule){
   #library(ggplot2)
   searchGroup<-function(node.id,TS,gName){
      flag<-TRUE
      sel.id<-TS[node.id,2:3]
      sel.group<-NULL
      i<-1
      while(sel.id[i]!= 0 && i <length(sel.id)){
         if(TS[sel.id[i],2]!=0){
            sel.id<-c(sel.id,TS[sel.id[i],2:3])
         }   
         if(TS[sel.id[i+1],2]!=0){
            sel.id<-c(sel.id,TS[sel.id[i+1],2:3])
         }
         i<-i+2
      }
      return(gName[sort(TS[sel.id[which(TS[sel.id,2]==0)],3])])
   }
   
   TS<-PPtreeOBJ$Tree.Struct
   Alpha<-PPtreeOBJ$projbest.node
   cut.off<-PPtreeOBJ$splitCutoff.node
   origdata<-PPtreeOBJ$origdata
   origclass<-PPtreeOBJ$origclass
   gName<-names(table(origclass))
   if(TS[node.id,2]!=0){
      selG<-searchGroup(node.id,TS,gName)
      sel.id<-NULL
      for(i in 1:length(selG)){
         sel.id<-c(sel.id,which(origclass==selG[i])) 
      }
      proj.data<-c(as.matrix(origdata)%*%as.matrix(Alpha[TS[node.id,4],]))[sel.id]
      proj.class<-origclass[sel.id]
      ..density.. <- NULL
      plot.data<-data.frame(proj.data = proj.data,origclass=proj.class)
      p1<- ggplot(plot.data, aes(x = proj.data,group=origclass))+
                geom_histogram( aes(y = ..density.., fill = origclass))+
                geom_vline(xintercept=cut.off[TS[node.id,4],Rule],linetype="longdash",lwd=1,col=2)
      vID <-1:ncol(origdata)
      coef.data<-data.frame(vID = vID,coef=Alpha[TS[node.id,4],])
      p2<-ggplot(coef.data,aes(x=vID,y=coef))+geom_bar(stat="identity",width=0.1)+
         geom_hline(yintercept=0) + geom_hline(yintercept=c(-1,1)*1/ncol(origdata),col=2,linetype="dashed") +
        xlab("variable ID")+ggtitle(paste("Node",node.id,sep=" "))
      gridExtra::grid.arrange(p2, p1,nrow=1)
   } else{
      sel.id<-which(origclass==gName[TS[node.id,3]])
      find.i<-which((TS[,2]==node.id | TS[,3]==node.id )& TS[,2]!=0)
      proj.data<-c(as.matrix(origdata)%*%as.matrix(Alpha[find.i,]))[sel.id]
      plot.data<-data.frame(proj.data)
      p1<-ggplot(plot.data,aes(x=proj.data))+
           geom_histogram(aes(y = ..density..),fill="gray40")+
           ggtitle(paste("Node",node.id,": ",gName[TS[node.id,3]],sep=""))  
      gridExtra::grid.arrange(p1,nrow=1)     
   }
}