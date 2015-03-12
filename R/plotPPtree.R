#' Plot PPtree
#' @usage plot.PPtree(PPtreeOBJ)
#' @param PPtreeOBJ PPtree object
#' @references Lee, E.-K., Cook, D. (2009) A Projection Pursuit Index for Large p Small n Data, 
#' Statistics and Computing,  
#' \url{http://www.springerlink.com/content/g47n0n342761838m/#?p=d2ff5a7b69eb45ef8abf7ef3aba69557&pi=3}.
#' @export
#' @keywords tree
#' @examples
#' data(iris)
#' Tree.result <- LDA.Tree.viz(iris[,5],iris[,1:4])
#' Tree.result
#' plot(Tree.result)
plot.PPtree<-function(PPtreeobj){
   library(grid)
   plotPPtree<-function(PPtreeobj,node.id,xlim,
                         ylim,nx = nx, ny = ny,tnex = tnex){
      TS<-PPtreeobj$Tree.Struct

      if(TS[node.id,2]==0) {
         x<-xlim[1] + diff(xlim)/2
         y<-ylim[1] + 0.5       
         tn_vp <- viewport(x = unit(x, "native"),
                     y = unit(y, "native") - unit(0.5, "lines"),
                     width = unit(1, "native"), 
                     height = unit(tnex, "native") - unit(1, "lines"),
                     just = c("center", "top"),
                     name = paste("Node", node.id, sep = ""))
         pushViewport(tn_vp)
         node_terminal.PP(PPtreeobj,node.id) 
         upViewport()
         return(NULL)
      }

      nl <- nleaf(PPtreeobj,node.id,"left")
      nr <- nleaf(PPtreeobj,node.id,"right")

      x0 <- xlim[1] + (nl / (nl + nr)) * diff(xlim)
      y0 <- max(ylim)

      if(TS[TS[node.id,2],2]==0){
         lf <- 1/2
      } else {
         lf <- nleaf(PPtreeobj,TS[node.id,2],"left")  / 
                      (nleaf(PPtreeobj,TS[node.id,2],"left") + 
                       nleaf(PPtreeobj,TS[node.id,2],"right"))
      }

      if(TS[TS[node.id,3],2]==0){
         rf <- 1/2
      } else {
         rf <- nleaf(PPtreeobj,TS[node.id,3],"left") / 
                   (nleaf(PPtreeobj,TS[node.id,3],"left") + 
                    nleaf(PPtreeobj,TS[node.id,3],"right"))      
      }

      x1l <- xlim[1] + (x0 - xlim[1]) * lf
      x1r <- x0 + (xlim[2] - x0) * rf
     
      y1l <- y1r <- y0 - 1

      grid.lines(x = unit(c(x0, x1l), "native"), 
                 y = unit(c(y0, y1l), "native"))
      grid.lines(x = unit(c(x0, x1r), "native"), 
                 y = unit(c(y0, y1r), "native"))
  
      in_vp <- viewport(x = unit(x0, "native"),
                        y = unit(y0, "native"),
                        width = unit(1, "native"),
                        height = unit(1, "native") - unit(1, "lines"), 
                        name = paste("Node", node.id, sep = ""))
      pushViewport(in_vp)

      node_inner.PP(PPtreeobj,node.id)
      upViewport()

      y1lr <- max(y1l, y1r)
      ypos <- y0 - (y0 - y1lr) * 0.5
      xlpos <- x0 - (x0 - x1l) * 0.5 * (y0 - y1lr)/(y0 - y1l)
      xrpos <- x0 - (x0 - x1r) * 0.5 * (y0 - y1lr)/(y0 - y1r)

      lsp_vp <- viewport(x = unit(xlpos, "native"),
                         y = unit(ypos, "native"),
                         width = unit(xlpos - xrpos, "native"),
                         height = unit(1, "lines"), 
                         name =  paste("lEdge", node.id, sep = ""))
      pushViewport(lsp_vp)

      edge_simple.PP(PPtreeobj,node.id, left = TRUE)
      upViewport()

      rsp_vp <- viewport(x = unit(xrpos, "native"),
                         y = unit(ypos, "native"),
                         width = unit(xlpos - xrpos, "native"),
                         height = unit(1, "lines"),
                         name =  paste("rEdge", node.id, sep = ""))
      pushViewport(rsp_vp) 

      edge_simple.PP(PPtreeobj,node.id, left = FALSE)
      upViewport()     

      plotPPtree(PPtreeobj,TS[node.id,2],c(xlim[1], x0),c(y1l, 1),nx,ny,tnex)
      plotPPtree(PPtreeobj,TS[node.id,3],c(x0, xlim[2]),c(y1r, 1),nx,ny,tnex)
   }    

   nleaf<-function(PPtreeobj,node.id,direction){  
      TS<-PPtreeobj$Tree.Struct
      n.leaf<-0
      flag<-TRUE
      if(direction=="left"){
         while(flag){
            if(TS[node.id,2]==0){
               flag<-FALSE
            } else {
               n.leaf <- n.leaf+1
               node.id<-TS[node.id,2]
            }
         }
      } else if(direction=="right"){
         while(flag){
            if(TS[node.id,2]==0){
               flag<-FALSE
            } else {
               n.leaf <- n.leaf+1
               node.id<-TS[node.id,3]
            }
         }
      } else {
         print("direction should be left or right")
      }
      return(n.leaf+1)
   }

   edge_simple.PP<- function(PPtreeobj,node.id,left=TRUE){   
      TS<-PPtreeobj$Tree.Struct
      if(left){
         text.t<-paste("< cut",TS[node.id,4],sep="")
         grid.rect(gp = gpar(fill = "white", col = 0), 
                   width = unit(1, "strwidth", text.t))
         grid.text(text.t, just = "center")
      } else{
         text.t<-paste(">= cut",TS[node.id,4],sep="")
         grid.rect(gp = gpar(fill = "white", col = 0), 
                  width = unit(1, "strwidth", text.t))
         grid.text(text.t, just = "center")
      } 
   }

   node_inner.PP <- function(PPtreeobj,node.id){   
      TS<-PPtreeobj$Tree.Struct
      PS<-print(PPtreeobj,verbose=FALSE)
      label1<-rep(NA,length(PS))
      label2<-rep(NA,length(PS))
      ID<-rep(NA,length(PS))
      final.group<-rep(NA,length(PS))
      temp<-strsplit(PS,"->")
      for(i in 1:length(temp)){
         t<-strsplit(temp[[i]][1],")" )
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
      node_vp <- viewport(x = unit(0.5, "npc"),
                        y = unit(0.5, "npc"),
                        width = unit(1, "strwidth", label.t) * 1.3, 
                        height = unit(3, "lines"),
                  name = paste("node_inner", node.id, sep = ""))
      pushViewport(node_vp)

      xell <- c(seq(0, 0.2, by = 0.01),
	              seq(0.2, 0.8, by = 0.05),
		            seq(0.8, 1, by = 0.01))
	    yell <- sqrt(xell * (1-xell))
	
      grid.polygon(x = unit(c(xell, rev(xell)), "npc"),
                   y = unit(c(yell, -yell)+0.5, "npc"),
                   gp = gpar(fill = "white"))
      grid.text(label.t, y = unit(1.2, "lines"))

      nodeIDvp <- viewport(x = unit(0.5, "npc"), 
                           y = unit(1, "npc"),
	                         width = max(unit(1, "lines"), 
                                       unit(1.3, "strwidth", 
	                                        as.character(node.id))),
	                         height = max(unit(1, "lines"), 
                                        unit(1.3, "strheight", 
	                                         as.character(node.id))))
      pushViewport(nodeIDvp)
      grid.rect(gp = gpar(fill = "white"))
      grid.text(node.id)
      popViewport()
      upViewport()
   }    

   node_terminal.PP <- function(PPtreeobj,node.id){
      TS<-PPtreeobj$Tree.Struct
      gName <- names(table(PPtreeobj$orig.class))
      gN <-gName[TS[node.id,3]]
      node_vp <- viewport(x = unit(0.5, "npc"),   
                          y = unit(0.5, "npc"),   
                          width = unit(1, "strwidth", gN) * 1.2,
                          height = unit(2, "lines"),
		                      name = paste("node_terminal", node.id, sep = ""))
      pushViewport(node_vp)	
      grid.rect(gp = gpar(fill = "lightgray"))
      grid.text(y = unit(1, "lines"), gN)

      nodeIDvp <- viewport(x = unit(0.5, "npc"), 
                           y = unit(1, "npc"),
                           width = max(unit(1, "lines"), 
                                unit(1.3, "strwidth", as.character(node.id))),
	                         height = max(unit(1, "lines"), 
                                unit(1.3, "strheight", as.character(node.id))))
      pushViewport(nodeIDvp)
      grid.rect(gp = gpar(fill = "lightgray", lty = "solid"))
      grid.text(node.id)
      popViewport()
      upViewport()    
   }

   calc.depth<-function(PPtreeobj){
      TS<-PPtreeobj$Tree.Struct
      i<-1;  flag.L<-rep(FALSE,nrow(TS))
      keep.track<-1
      depth.track<-0
      depth<-0
      while(sum(flag.L)!=nrow(TS)){
         if(!flag.L[i]) {                    
            if(TS[i,2] == 0) {
               flag.L[i]<-TRUE
               id.l<-length(keep.track)-1
               i <- keep.track[id.l]
               depth <- depth -1
            } else if(!flag.L[TS[i,2]]) {
               depth <- depth +1
               i<-TS[TS[i,2],1]   
            } else {
               depth <- depth +1
               flag.L[i]<-TRUE
               i<-TS[TS[i,3],1]
            } 
            keep.track <- c(keep.track,i)
            depth.track<-c(depth.track,depth)
         } else {
            id.l<-id.l-1
            i<-keep.track[id.l]
            depth <- depth.track[id.l]
         }
      }
      depth<-max(depth.track)+2
      return(depth)
   }

   nx <- length(table(PPtreeobj$orig.class))
   ny <- calc.depth(PPtreeobj)
   tnex<-2
   node.id<-1
   main<-"PPtree"

   grid.newpage()

   root_vp <- viewport(layout = grid.layout(3, 3, 
                          heights = unit(c(ifelse(is.null(main), 0, 3), 1, 1), 
                                         c("lines", "null", "lines")),
    			                widths = unit(c(1, 1, 1), c("lines", "null", "lines"))), 
    			             name = "root")       
   pushViewport(root_vp)
   main_vp <- viewport(layout.pos.col = 2, layout.pos.row = 1,name = "main")
   pushViewport(main_vp)
   grid.text(y = unit(1, "lines"), main, just = "center")
   upViewport()

   tree_vp <- viewport(layout.pos.col = 2, layout.pos.row = 2, 
    			xscale = c(0, nx), yscale = c(0, ny + (tnex - 1)), 
                        name = "tree")
   pushViewport(tree_vp)
   plotPPtree(PPtreeobj,1, c(0, nx), ylim = c(0, ny - 0.5 + (tnex - 1)), nx, ny, tnex)
}
