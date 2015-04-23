#' Plot PP regression tree
#' 
#' plot projection pursuit regression tree structure
#' @usage plot.PPtreereg(x,font.size=17,width.size=1,...)
#' @param x regPPtree object
#' @param font.size font size of plot
#' @param width.size size of eclipse in each node.
#' @param ... arguments to be passed to methods
#' @export
#' @keywords tree
#' @aliases plot
#' @examples
#' data(mtcar)
#' Tree.result <- PP.Tree.reg(mtcars[,1],mtcars[,-1],final.rule=1,DEPTH=2)
#' plot(Tree.result)
plot.PPtreereg<-function(x,font.size=17,width.size=1,...){
   
   PPregtreeOBJ<-x
   PPtreeobj<-PPregtreeOBJ$Tree.result
   class(PPtreeobj)<-append(class(PPtreeobj),"PPtree")
   plotPPregtree<-function(PPtreeobj,node.id,xlim,ylim){
      TS<-PPtreeobj$Tree.Struct
      if(TS[node.id,2]==0) {
         x<-xlim[1]+0.5  
         y<-ylim[2]-1     
         Final.Node.V<-viewport(x=unit(x,"native"),
                                y=unit(y,"native"),
                                width=unit(1,"native"), 
                                height=unit(1,"native")-unit(2,"lines"),
                                just=c("center","top"))
         pushViewport(Final.Node.V)
         node.terminal.PPtree(PPtreeobj,node.id) 
         upViewport()
         return(NULL)
      }
      nl<-n.final(PPtreeobj,node.id,"left")
      nr<-n.final(PPtreeobj,node.id,"right")
      x0<-xlim[1]+nl 
      y0<-max(ylim)-1 
      lf<-ifelse(TS[TS[node.id,2],2]==0,0.5,
                   n.final(PPtreeobj,TS[node.id,2],"right"))  
      rf<-ifelse(TS[TS[node.id,3],2]==0,0.5,
                   n.final(PPtreeobj,TS[node.id,3],"left")) 
      x1l<-x0-lf
      x1r<-x0+rf
      y1<-y0-1
      grid.lines(x=unit(c(x0,x1l),"native"),y=unit(c(y0,y1),"native"))
      grid.lines(x=unit(c(x0,x1r),"native"),y=unit(c(y0,y1),"native")) 
      node.V<-viewport(x=unit(x0,"native"),
                       y=unit(y0,"native"),
                       width=unit(1,"native"),
                       height=unit(1,"native")-unit(1,"lines"))
      pushViewport(node.V)
      node.inner.PPtree(PPtreeobj,node.id)
      upViewport()
      ylpos<-y0-0.6
      yrpos<-y0-0.45
      xlpos<-x0-(x0-x1l)*0.6 
      xrpos<-x0-(x0-x1r)*0.45 
      LeftEdge.V<-viewport(x=unit(xlpos,"native"),
                           y=unit(ylpos,"native"),
                           width=unit(xlpos-xrpos,"native"),
                           height=unit(1,"lines")*1.2)
      pushViewport(LeftEdge.V)
      edge.lable.PPtree(PPtreeobj,node.id,left=TRUE)
      upViewport()
      RightEdge.V<-viewport(x=unit(xrpos,"native"),
                            y=unit(yrpos,"native"),
                            width=unit(xlpos-xrpos,"native"),
                            height=unit(1,"lines"))
      pushViewport(RightEdge.V) 
      edge.lable.PPtree(PPtreeobj,node.id,left=FALSE)
      upViewport()     
      plotPPregtree(PPtreeobj,TS[node.id,2],c(xlim[1],x0),c(1,y1+1))
      plotPPregtree(PPtreeobj,TS[node.id,3],c(x0,xlim[2]),c(1,y1+1))
   }  

   n.final<-function(PPtreeobj,node.id,direction){  
      TS<-PPtreeobj$Tree.Struct
      n.leaf<-0
      if(direction=="left"){
         keep.id<-TS[node.id,2]
         i<-1
         while(i<=length(keep.id)){
            if(TS[keep.id[i],2]==0){
               n.leaf<-n.leaf+1
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
               n.leaf<-n.leaf+1
               i<-i+1
            } else{  
               keep.id<-c(keep.id,TS[keep.id[i],2:3])
               i<-i+1
            }                               
         }  
      }
      return(n.leaf)
   }

   edge.lable.PPtree<- function(PPtreeobj,node.id,left=TRUE){   
      TS<-PPtreeobj$Tree.Struct
      if(left){
         text.t<-paste("< cut",TS[node.id,4],sep="")
         grid.rect(gp=gpar(fill="white",lty=1,col="grey95"), 
                   width=unit(width.size,"strwidth",text.t)*1.2)
         grid.text(text.t,just="center",gp=gpar(fontsize=font.size))
      } else{
         text.t<-paste(">= cut",TS[node.id,4],sep="")
         grid.rect(gp=gpar(fill="white",lty=1,col="grey95"), 
                   width=unit(width.size,"strwidth",text.t)*1.2)
         grid.text(text.t,just="center",gp=gpar(fontsize=font.size))
      } 
   }

   node.inner.PPtree <- function(PPtreeobj,node.id){   
      TS<-PPtreeobj$Tree.Struct
      PPtreeobj
      class(PPtreeobj)<-append(class(PPtreeobj),"PPtree")
      PS<-print(PPtreeobj,verbose=FALSE)
      label1<-rep(NA,length(PS))
      label2<-rep(NA,length(PS))
      ID<-rep(NA,length(PS))
      final.group<-rep(NA,length(PS))
      temp<-strsplit(PS,"->")
      for(i in 1:length(temp)){
         t<-strsplit(temp[[i]][1],")")
         ID[i]<-as.numeric(t[[1]][1])
         tt<-strsplit(t[[1]][2]," ")[[1]]
         tt<-tt[tt!="" & tt!="*"]
         label1[i]<-tt[1]
         if(tt[1]!="root")
            label2[i]<-paste(tt[2],tt[3])
         if(length(temp[[i]])==2)
            final.group[i]<-temp[[i]][2]
      }
      label.t<-paste("proj",TS[node.id,4]," * X",sep="")
      Inner.Node.V<-viewport(x=unit(0.5,"npc"),
                             y=unit(0.5,"npc"),
                             width=unit(width.size*1.5,"strwidth",label.t), 
                             height=unit(width.size*2,"lines"))
      pushViewport(Inner.Node.V)
      xell<-c(seq(0,0.2,by=0.01),
              seq(0.2,0.8,by=0.05),
  	          seq(0.8,1,by=0.01))
	    yell<-sqrt(xell*(1-xell))
      grid.polygon(x=unit(c(xell,rev(xell)),"npc"),
                   y=unit(c(yell,-yell)+0.5,"npc"),
                   gp=gpar(fill="white"))
      grid.text(label.t,y=0.3,gp=gpar(fontsize=font.size),
                just=c("center","bottom"))
      Inner.Node.Id.V<-viewport(x=unit(0.5,"npc"), 
                                y=unit(1,"npc"),
	                              width=max(unit(1,"lines"), 
                                          unit(1.2,"strwidth",
                                               as.character(node.id))),
	                              height=max(unit(1,"lines"), 
                                           unit(1.2,"strheight",
                                                as.character(node.id))),
                                just=c("center","center"),
                                gp=gpar(fontsize=font.size))
      pushViewport(Inner.Node.Id.V)
      grid.rect(gp=gpar(fill="white",lty="solid",fontsize=font.size))
      grid.text(node.id,gp=gpar(fontsize=font.size))
      popViewport()
      upViewport()
   }    

   node.terminal.PPtree<- function(PPtreeobj,node.id){
      TS<-PPtreeobj$Tree.Struct
      gName1<-as.character(round(PPregtreeOBJ$mean.G,3))
      gName2<-as.character(round(PPregtreeOBJ$sd.G,3))      
      gN<-paste(gName1[TS[node.id,3]],"\n(",gName2[TS[node.id,3]],")",sep="")
      temp<-strsplit(as.character(gN),split="")[[1]]
      gN.width<-length(temp)      
      set.unit<-length(sum(tolower(temp)!=temp)*0.65+
                         sum(tolower(temp)==temp)*0.5)/gN.width
      Terminal.Node.V<-viewport(x=unit(0.5,"npc"),   
                                y=unit(0.8,"npc"),   
                                height=unit(1,"lines")*3,
                                width=unit(0.3,"lines")*(gN.width),
  	                            just=c("center","top"))
      pushViewport(Terminal.Node.V )	
      grid.rect(gp=gpar(fill="lightgray"))
      grid.text(y=0.05,gN,gp=gpar(fontsize=font.size),
                just=c("center","bottom"))
      Terminal.Node.Id.V<-viewport(x=unit(0.5,"npc"), 
                                   y=unit(1,"npc"),
                                   width=max(unit(1,"lines"), 
                                             unit(1.2,"strwidth",
                                                  as.character(node.id))),
                                   height=max(unit(1,"lines"), 
                                              unit(1.2,"strheight",
                                                   as.character(node.id))),
                                   just=c("center","center"),
                                   gp=gpar(fontsize=font.size))
      pushViewport(Terminal.Node.Id.V)
      grid.rect(gp=gpar(fill="lightgray",lty="solid",fontsize=font.size))
      grid.text(node.id,gp=gpar(fontsize=font.size))
      popViewport()
      upViewport()    
   }

   calc.depth<-function(PPtreeobj){
      TS<-PPtreeobj$Tree.Struct
      i<-1;  
      flag.L<-rep(FALSE,nrow(TS))
      keep.track<-1
      depth.track<-0
      depth<-0
      while(sum(flag.L)!=nrow(TS)){
         if(!flag.L[i]) {                    
            if(TS[i,2] == 0) {
               flag.L[i]<-TRUE
               id.l<-length(keep.track)-1
               i<-keep.track[id.l]
               depth<-depth -1
            } else if(!flag.L[TS[i,2]]) {
               depth<-depth +1
               i<-TS[TS[i,2],1]   
            } else {
               depth<-depth +1
               flag.L[i]<-TRUE
               i<-TS[TS[i,3],1]
            } 
            keep.track<-c(keep.track,i)
            depth.track<-c(depth.track,depth)
         } else {
            id.l<-id.l-1
            i<-keep.track[id.l]
            depth<-depth.track[id.l]
         }
      }
      depth<-max(depth.track)+2
      return(depth)
   }

   nx<-length(PPregtreeOBJ$mean.G)
   ny<-calc.depth(PPtreeobj)
   tnex<-1
   node.id<-1
   grid.newpage()
   PPtree.Main.V<-viewport(layout=grid.layout(3,3, 
                                    heights=unit(c(3,1,1),
                                                 c("lines","null","lines")),
     			                          widths=unit(c(1,1,1),
                                                 c("lines","null","lines"))))       
   pushViewport(PPtree.Main.V)
   PPtree.title.V<-viewport(layout.pos.col=2,layout.pos.row=1)
   pushViewport(PPtree.title.V)
   grid.text(y=unit(1,"lines"),"Projection Pursuit Regression Tree",
             just="center")
   upViewport()
   PPtree.Tree.V<-viewport(layout.pos.col=2,layout.pos.row=2, 
   			                   xscale=c(0,nx),yscale=c(0,ny+1))
   pushViewport(PPtree.Tree.V)
   plotPPregtree(PPtreeobj,1,c(0,nx),ylim=c(1,ny+1))   
}
