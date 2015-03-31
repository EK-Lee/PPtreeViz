#' Visualization of each node in PPtree
#' 
#' Explore PPtree with different Rules in each split.
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
   p<-ncol(origdata)
   gName<-names(table(origclass))
   if(TS[node.id,2]!=0){
      selG<-searchGroup(node.id,TS,gName)
      sel.id<-NULL
      for(i in 1:length(selG)){
         sel.id<-c(sel.id,which(origclass==selG[i])) 
      }
      proj.data<-c(as.matrix(origdata)%*%as.matrix(Alpha[TS[node.id,4],]))[sel.id]
      proj.class<-origclass[sel.id]
      plot.data<-data.frame(proj.data = proj.data,origclass=proj.class)
      p1<- ggplot(plot.data, aes(x = proj.data,group=origclass))+
                geom_histogram( aes(fill = origclass),position="stack")+
                geom_vline(xintercept=cut.off[TS[node.id,4],Rule],linetype="longdash",lwd=1,col=2)
      vID <-1:p
      coef.data<-data.frame(vID = vID,coef=Alpha[TS[node.id,4],])
      bin.width<-ifelse(p>100,1,0.1)
      y.max <-max(c(abs(coef.data$coef),1/sqrt(p)))
      
      p2<-ggplot(coef.data,aes(x=vID,y=coef))+geom_bar(stat="identity",width=bin.width)+
         geom_hline(yintercept=0) + geom_hline(yintercept=c(-1,1)*1/sqrt(ncol(origdata)),col=2,linetype="dashed") +
        xlab("variable ID")+ggtitle(paste("Node",node.id,sep=" "))+ylim(-y.max,y.max)
      gridExtra::grid.arrange(p2, p1,nrow=1)
   } else{
      #sel.id<-which(origclass==gName[TS[node.id,3]])
      sel.id<-which(PP.classify(PPtreeOBJ,origdata,Rule)$predict.class==gName[TS[node.id,3]])
      find.i<-TS[which((TS[,2]==node.id | TS[,3]==node.id )& TS[,2]!=0),4]
      proj.data<-c(as.matrix(origdata)%*%as.matrix(Alpha[find.i,]))[sel.id]
      proj.class<-origclass[sel.id]
      plot.data<-data.frame(proj.data=proj.data,proj.class=proj.class)
      
      p1<-ggplot(plot.data,aes(x=proj.data,group=proj.class))+
           geom_histogram(aes(fill=proj.class),position="stack")+
           ggtitle(paste("Node",node.id,": ",gName[TS[node.id,3]],sep=""))  
      gridExtra::grid.arrange(p1,nrow=1)     
   }
}