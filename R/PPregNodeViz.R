#' Visualization of each node in PPtree
#' 
#' Explore PPtree with different Rules in each split.
#' @usage PPregNode.Viz(PPtreeregOBJ,node.id,Rule)
#' @param PPtreeregOBJ  PPregtree object
#' @param node.id node ID
#' @param Rule cutoff rule
#' @export
#' @keywords tree
#' @examples
#' data(mtcars)
#' Tree.result <- PP.Tree.reg(mtcars[,1],mtcars[,-1],final.rule=1,DEPTH=2)
#' plot(Tree.result)
#' PPregNode.Viz(Tree.result,1,1)
PPregNode.Viz<-function(PPtreeregOBJ,node.id,Rule){
  
   searchGroup<-function(node.id,TS,gName){
      flag<-TRUE
      sel.id<-TS[node.id,2:3]
      sel.group<-NULL
      i<-1
      while((sel.id[i]!=0)&&(i<length(sel.id))){
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
   
   final.search<-function(PPtreeobj,node.id,direction){  
      TS<-PPtreeobj$Tree.Struct
      leaf.group<-NULL
      if(direction=="left"){
         keep.id<-TS[node.id,2]
         i<-1
         while(i<=length(keep.id)){
            if(TS[keep.id[i],2]==0){
               leaf.group<-c(leaf.group,TS[keep.id[i],3])
               i<-i+1
            } else{  
               keep.id<-c(keep.id,TS[keep.id[i],2:3])
               i<-i+1
            }                               
         }  
      } else if(direction=="right"){
         keep.id<-TS[node.id,3]
         i<-1
         while(i<=length(keep.id)){
            if(TS[keep.id[i],2]==0){
               leaf.group<-c(leaf.group,TS[keep.id[i],3])               
               i<-i+1
            } else{  
               keep.id<-c(keep.id,TS[keep.id[i],2:3])
               i<-i+1
            }                               
         }  
      }
      return(leaf.group)
   }

   PPtreeOBJ<-PPtreeregOBJ$Tree.result
   TS<-PPtreeOBJ$Tree.Struct
   Alpha<-PPtreeOBJ$projbest.node
   cut.off<-PPtreeOBJ$splitCutoff.node
   origdata<-PPtreeOBJ$origdata
   origclass<-PPtreeOBJ$origclass
   p<-ncol(origdata)
   gName<-names(table(origclass))
   if(node.id==0){
      Y<-PPtreeregOBJ$origY
      proj.data<-rep(1,length(Y))
      predY<-PPreg.predict(PPtreeregOBJ,origdata)$predict.Y
      pred.data<-predY
      plot.data<-data.frame(pred.data,Y)

      p1<-ggplot(plot.data,aes(x=pred.data,y=Y))+ 
          geom_point(position="jitter")+
          xlab("Fitted Y")+
          ylab("Observed Y")+
          geom_abline(col=2)+
          ggtitle("Fitted vs. Observed")  
      gridExtra::grid.arrange(p1,nrow=1)      
   } else if(TS[node.id,2]!=0){
      selG<-searchGroup(node.id,TS,gName)
      sel.id<-NULL
      for(i in 1:length(selG)){
         sel.id<-c(sel.id,which(origclass==selG[i])) 
      }
      proj.data<-c(as.matrix(origdata)%*%
                     as.matrix(Alpha[TS[node.id,4],]))[sel.id]
      proj.class<-factor(round(PPtreeregOBJ$mean.G[origclass[sel.id]],3))
      Y<-PPtreeregOBJ$origY[sel.id]
      plot.data<-data.frame(proj.data,origclass=proj.class,Y)
      min.X<-min(proj.data)
      max.X<-max(proj.data)
      left.group<-final.search(PPtreeOBJ,node.id,"left")
      right.group<-final.search(PPtreeOBJ,node.id,"right")   
      x1<-c(rep(c(min.X,cut.off[TS[node.id,4],Rule]),
                                         length(left.group)),
                                 rep(c(cut.off[TS[node.id,4],Rule],max.X),
                                         length(right.group)))
      y1<-rep(PPtreeregOBJ$mean.G[c(left.group,right.group)],each=2)
      line.data<-data.frame(x1,y1,origclass=factor(round(PPtreeregOBJ$mean.G[
                                 rep(c(left.group,right.group),each=2)],3)))
      p1<- ggplot()+
           geom_point(data = plot.data,aes(x = proj.data,y=Y,color=origclass))+
           geom_vline(xintercept=cut.off[TS[node.id,4],Rule],
                      linetype="longdash",lwd=1,col=2)+
           geom_line(data = line.data,aes(x=x1,y=y1,color=origclass),lwd=1)
      vID <-1:p
      coef.data<-data.frame(vID = vID,coef=Alpha[TS[node.id,4],])
      bin.width<-ifelse(p>100,1,0.1)
      y.max <-max(c(abs(coef.data$coef),1/sqrt(p)))     
      p2<-ggplot(coef.data,aes(x=vID,y=coef))+
          geom_segment(aes(yend=0,xend=vID,width=0.1))+        
          geom_hline(yintercept=0)+ 
          geom_hline(yintercept=c(-1,1)*1/sqrt(ncol(origdata)),
                     col=2,linetype="dashed")+        
          xlab("variable ID")+
          ggtitle(paste("Node",node.id,sep=" "))+
          ylim(-y.max,y.max)
      gridExtra::grid.arrange(p2,p1,nrow=1)
   } else{
      sel.id<-which(PP.classify(PPtreeOBJ,origdata,Rule)$predict.class==
                      gName[TS[node.id,3]])
      Yorig<-PPtreeregOBJ$origY
      proj.data<-rep(1,length(Y))
      proj.data<-c(as.matrix(origdata)%*%
                     as.matrix(Alpha[which((TS[,2]!=0&TS[,3]==node.id)|
                                             TS[,2]==node.id),]))      
#      predY<-PPreg.predict(PPtreeregOBJ,origdata)$predict.Y
      coef<-lm(Y[sel.id]~proj.data[sel.id])$coef
      predY<-proj.data
      pred.data<-predY[sel.id];Y=Yorig[sel.id]
      plot.data1<-data.frame(pred.data,Y)
      pred.data<-predY[-sel.id];Y=Yorig[-sel.id]
      plot.data2<-data.frame(pred.data,Y)     
      p1<- ggplot()+
                geom_point(data=plot.data2,
                           aes(x=pred.data,y=Y),col="grey50",size=1.3)+
                geom_point(data=plot.data1,aes(x=pred.data,y=Y))+
                geom_hline(yintercept=PPtreeregOBJ$mean.G[TS[node.id,3]],
                           col=4,linetype="dashed")+
#                geom_abline(intercept=0,slope=1,lwd=1,col=2)+       
                geom_abline(intercept=coef[1],slope=coef[2],col=2,linetype="dashed")+       
                xlab("Projected X")+ylab("Observed Y")+
                ggtitle(paste("Node",node.id,": ",
                              round(PPtreeregOBJ$mean.G[TS[node.id,3]],3),
                             "(",round(PPtreeregOBJ$sd.G[TS[node.id,3]],3),")",
                             sep=""))   
      predY<-PPreg.predict(PPtreeregOBJ,origdata)$predict.Y
      pred.data<-predY[sel.id];Y=Yorig[sel.id]
      plot.data3<-data.frame(pred.data,Y)
      pred.data<-predY[-sel.id];Y=Yorig[-sel.id]
      plot.data4<-data.frame(pred.data,Y)     

      p2<- ggplot()+
                geom_point(data=plot.data4,
                           aes(x=pred.data,y=Y),col="grey50",size=1.3)+
                geom_point(data=plot.data3,aes(x=pred.data,y=Y))+
                geom_hline(yintercept=PPtreeregOBJ$mean.G[TS[node.id,3]],
                           col=4,linetype="dashed")+
                geom_abline(intercept=0,slope=1,lwd=1,col=2)+       
                xlab("Fitted Y")+ylab("Observed Y")+
                ggtitle(paste("Node",node.id,": ",
                              round(PPtreeregOBJ$mean.G[TS[node.id,3]],3),
                             "(",round(PPtreeregOBJ$sd.G[TS[node.id,3]],3),")",
                             sep="")) 
      gridExtra::grid.arrange(p1,p2,nrow=1)     
   }
}