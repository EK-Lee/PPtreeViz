#' Draw Huber plot for 2D data with various PP index
#' 
#' Draw Huber plot for 2D data with various PP index
#' @usage Huber.plot(origdata2D,origclass,PPmethod,weight,r,lambda,
#'            opt.proj,UserDefFtn,...)
#' @param origdata2D 2D data for Huber plot
#' @param origclass class information
#' @param PPmethod method for projection pursuit; "LDA", "PDA", "Lr", "GINI", and "ENTROPY"
#' @param weight weight flag in LDA, PDA and Lr index
#' @param r r in Lr index
#' @param lambda lambda in PDA index 
#' @param opt.proj flag to show the best projection in the plot
#' @param UserDefFtn  User defined index function
#' @param ...  arguments to be passed to methods
#' @references Lee, EK., Cook, D., Klinke, S., and Lumley, T.(2005) 
#' Projection Pursuit for exploratory supervised classification, 
#' Journal of Computational and Graphical statistics, 14(4):831-846.
#' @export
#' @keywords projection pursuit
#' @examples
#' data(iris)
#' Huber.plot(iris[,1:2],iris[,5],PPmethod="LDA")
Huber.plot<-function(origdata2D,origclass,PPmethod="LDA",
                     weight=TRUE,r=1,lambda=0.5,opt.proj=TRUE,
                     UserDefFtn=NULL,...){  
   index<-NULL
   best.proj<-NULL
   best.index<-0
   origdata2D<-as.matrix(origdata2D)
   for(i in 0:360){   
      theta<-pi/180*i
      proj.data<-matrix(cos(theta)*origdata2D[,1]+sin(theta)*origdata2D[,2])
      proj<-matrix(c(cos(theta),sin(theta)),ncol=1)
      if(PPmethod=="LDA"){
         newindex<-LDAindex(origclass,origdata2D,proj=proj,weight=weight)
      } else if(PPmethod=="PDA"){
         newindex<-PDAindex(origclass,origdata2D,proj,weight=weight,
                            lambda=lambda)
      } else if(PPmethod=="Lr"){
         newindex<-Lrindex(origclass,origdata2D,proj,weight=weight,r=r)          
      } else if(PPmethod=="GINI"){
         newindex<-GINIindex1D(origclass,origdata2D,proj)          
      } else if(PPmethod=="ENTROPY"){
        newindex<-ENTROPYindex1D(origclass,origdata2D,proj)
      } else if(PPmethod=="UserDef"){
        newindex<-UserDefFtn(proj.data)
      } 
      index<-c(index,newindex)
   }
   sel.index<-which(index[1:360]>signif(max(index),6)-1.0E-6)
   theta.best.all<-pi/180*(sel.index-1)
   theta.best<-theta.best.all[1]
   proj.data.best<-matrix(cos(theta.best)*origdata2D[,1]+
                            sin(theta.best)*origdata2D[,2])
   index.best<-max(index)
   range<-round(max(index)-min(index),5)
   if(range==0){
      PPindex<-rep(4,length(index))
   } else {
      PPindex<-(index-min(index))/range*2+3
   }
   data.circle<-NULL
   data.index<-NULL
   for(i in 1:361){  
      theta<-pi/180*(i-1)
      data.index<-rbind(data.index,
                 c(PPindex[i]*cos(theta),PPindex[i]*sin(theta)))
      data.circle<-rbind(data.circle,c(4*cos(theta),4*sin(theta)))
   }
   orig.scaled<-apply(origdata2D,2,function(x) (x-mean(x))/diff(range(x))*3.5)
   data.cX<-data.circle[,1]
   data.cY<-data.circle[,2]
   data.X<-data.index[,1]
   data.Y<-data.index[,2]
   plot.data<-data.frame(data.cX,data.cY,data.X,data.Y)
   x<-orig.scaled[,1]
   y<-orig.scaled[,2]
   group<-origclass
   point.data<-data.frame(x,y,group)
   P1<-ggplot(data=plot.data,aes(x=data.X,y=data.Y))+
         geom_path()+
         geom_path(aes(x=data.cX,y=data.cY),linetype="dashed")+
         geom_point(data=point.data,aes(x=x,y=y,color=group,shape=group))+
         scale_x_continuous(breaks=NULL)+
         scale_y_continuous(breaks=NULL)+
         xlab("")+ylab("")+coord_fixed()+
         theme_bw()+theme(panel.border=element_blank())
   if(opt.proj){
      P1<-P1+geom_abline(intercept=0,slope=sin(theta.best)/cos(theta.best))
      if(length(theta.best.all)>1)
         for(i in 2:length(theta.best.all))
            P1<-P1+
                geom_abline(intercept=0,          
                            slope=sin(theta.best.all[i])/cos(theta.best.all[i]),
                            linetype="dashed")
   }    
   best.proj.data<-proj.data.best
   group<-origclass
   hist.data<-data.frame(best.proj.data,group)
   P2<-ggplot(data=hist.data,aes(x=best.proj.data,group=group))+
       geom_histogram(aes(fill=group),position="stack")
   gridExtra::grid.arrange(P1,P2,nrow=1)
}