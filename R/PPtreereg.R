#' Projection pursuit regression tree with various PP indices
#' 
#' Find tree structure using linear discriminant(LD) in each split.
#' @usage PP.Tree.reg(origY,origX,final.rule=1,DEPTH=NULL,Rr=1,PPmethod="LDA",weight=TRUE,lambda=0.1,r=1,TOL.CV=0.1,selP=NULL,energy=0,maxiter=500,...) 
#' @param origY dependent data vector
#' @param origX independent data matrix without dependent variable
#' @param final.rule rule to calculate the final node value
#' @param DEPTH depth of the projection pursuit regression tree
#' @param Rr cutoff rule in each node
#' @param PPmethod method for projection pursuit; "LDA", "PDA", "Lr", "GINI", and "ENTROPY"
#' @param weight weight flag in LDA, PDA and Lr index
#' @param lambda lambda in PDA index
#' @param r r in Lr index
#' @param TOL.CV CV limit for the final node
#' @param selP number of variables for the final node in Method 5
#' @param energy energy parameter
#' @param maxiter number of maximum iteration 
#' @param ... arguments to be passed to methods
#' @return Tree.result projection pursuit regression tree result with PPtreeclass object format
#' @return MSE mean squared error of the final tree
#' @return mean.G  means of the observations in the final node 
#' @return sd.G standard deviations of the observations in the final node. 
#' @return coef.G regression coefficients for Method 3, 4 and 5
#' @return origY original dependent variable vector
#' @references ...
#' @export
#' @keywords tree  
#' @examples
#' data(mtcars)
#' Tree.result <- PP.Tree.reg(mtcars[,1],mtcars[,-1],final.rule=1,DEPTH=2,PPmethod="LDA")
#' Tree.result
PP.Tree.reg<-function(origY,origX,final.rule=1,DEPTH=NULL,Rr=1,PPmethod="LDA",
                      weight=TRUE,lambda=0.1,r=1,TOL.CV=0.1,selP=NULL,
                      energy=0,maxiter=500,...){   
   Find.proj<-function(origclass,origdata,PPmethod="LDA",weight=TRUE,
                       lambda=0.1,r=1,...){
      n<-nrow(origdata)
      p<-ncol(origdata)
      g<-table(origclass)
      g.name<-as.numeric(factor(names(g)))
      G<-length(g)
      origclass<-as.numeric(factor(origclass))
      if(PPmethod=="LDA"){
         indexbest<-LDAindex(origclass,as.matrix(origdata),weight=weight);
      } else if(PPmethod=="PDA"){
         indexbest<-PDAindex(origclass,as.matrix(origdata),
                             weight=weight,lambda=lambda);
      } else if(PPmethod=="Lr"){
         indexbest<-Lrindex(origclass,as.matrix(origdata),
                            weight=weight,r=r);
      } else if(PPmethod=="GINI"){
         indexbest<-0;
         for(i in 1:p){
            tempdata<-origdata[,i];
            tempindex<-GINIindex1D(origclass,as.matrix(tempdata));  
            if(indexbest<tempindex)
               indexbest<-tempindex;
         }  
      } else if(PPmethod=="ENTROPY"){
         indexbest<-0;
         for(i in 1:p){
            tempdata<-origdata[,i];
            tempindex<-ENTROPYindex1D(origclass,as.matrix(tempdata));  
            if(indexbest<tempindex)
               indexbest<-tempindex;
         }  
      }       
      energy<-ifelse(energy==0,1-indexbest,energy)
      energy.temp<-1-indexbest
      TOL.temp<-energy.temp/1000000        
      if(PPmethod=="LDA"){
         a.proj.best<-LDAopt(as.numeric(as.factor(origclass)),
                             as.matrix(origdata),weight,q=1)$projbest
      } else if(PPmethod=="PDA"){
         a.proj.best<-PDAopt(as.numeric(as.factor(origclass)),
                               as.matrix(origdata),weight,
                               q=1,lambda=lambda)$projbest          
      } else {
         a.proj.best<-PPopt(as.numeric(as.factor(origclass)),
                            as.matrix(origdata),PPmethod=PPmethod,
                            r=r,q=1,energy=energy,cooling=0.999,
                            TOL=TOL.temp)$projbest          
      }
      proj.data<-as.matrix(origdata)%*%a.proj.best   
      if(diff(tapply(proj.data,origclass,mean))<0)
         a.proj.best<--a.proj.best
      proj.data<-as.matrix(origdata)%*%a.proj.best          
      class<-origclass
      m.LR<-tapply(proj.data,class,mean)
      sd.LR<-tapply(proj.data, class, function(x) 
                                         ifelse(length(x)>1,sd(x),0))
      IQR.LR<-tapply(proj.data, class, function(x) 
                                         ifelse(length(x)>1,IQR(x),0))
      median.LR<-tapply(proj.data, class, median)
      n.LR<-table(class)        
      c1<-(m.LR[1]+m.LR[2])/2
      c2<-(m.LR[1]*n.LR[2]+m.LR[2]*n.LR[1])/sum(n.LR)
      c3<-ifelse(sum(sd.LR==0)!=0,c1,(m.LR[1]*sd.LR[2]+m.LR[2]*sd.LR[1])/
                                      sum(sd.LR))
      c4<-ifelse(sum(sd.LR==0)!=0,c2,(m.LR[1]*sd.LR[2]/sqrt(n.LR[2])+ 
                                      m.LR[2]*sd.LR[1]/sqrt(n.LR[1]))/
                           (sd.LR[1]/sqrt(n.LR[1])+sd.LR[2]/sqrt(n.LR[2])))
      c5<-(median.LR[1]+median.LR[2])/2  
      c6<-(median.LR[1]*n.LR[2]+median.LR[2]*n.LR[1])/sum(n.LR)    
      c7<-ifelse(sum(IQR.LR==0)!=0,c5,(median.LR[1]*IQR.LR[2]+
                                         median.LR[2]*IQR.LR[1])/sum(IQR.LR))      
      c8<-ifelse(sum(IQR.LR==0)!=0,c6,(median.LR[1]*(IQR.LR[2]/sqrt(n.LR[2]))+
                                     median.LR[2]*(IQR.LR[1]/sqrt(n.LR[1])))/
                         ((IQR.LR[1]/sqrt(n.LR[1]))+(IQR.LR[2]/sqrt(n.LR[2]))))
      sel.proj<-sort(proj.data[which(proj.data>quantile(proj.data,prob=0.25)&
                                    proj.data<quantile(proj.data,prob=0.75))])
      sel.n<-length(sel.proj)
      temp.cut<-matrix((sel.proj[2:sel.n]+sel.proj[1:(sel.n-1)])/2,ncol=1)
      c9<-sel.proj[sort.list(apply(temp.cut,1,function(x) 
                                { temp<-table(class,proj.data>x[1]);
                                  return(prod(temp[,1])+prod(temp[,2]))}))[1]]        
      C<-c(c1, c2, c3, c4,c5,c6,c7,c8,c9)
      if(PPmethod=="LDA"){
         Index<-LDAindex(as.numeric(as.factor(class)),as.matrix(proj.data),
                         weight=weight)
      } else if(PPmethod=="PDA"){
         Index<-PDAindex(as.numeric(as.factor(class)),as.matrix(proj.data),
                         weight=weight,lambda=lambda)          
      } else if(PPmethod=="Lr"){
         Index<-Lrindex(as.numeric(as.factor(class)),as.matrix(proj.data),
                        weight=weight,r=r)          
      } else if(PPmethod=="GINI"){
         Index<-GINIindex1D(as.numeric(as.factor(class)),as.matrix(proj.data))          
      } else if(PPmethod=="ENTROPY"){
         Index<-ENTROPYindex1D(as.numeric(as.factor(class)),
                               as.matrix(proj.data))          
      } 
      Alpha<-t(a.proj.best)
      IOindexR<-NULL
      IOindexL<-NULL
      IOindexL<-class==1
      IOindexR<-class==2
      list(Index=Index,Alpha=Alpha,C=C,IOindexL=IOindexL,IOindexR=IOindexR)
   }

   Tree.construct<-function(origY,origdata,Tree.Struct,id,depth,
                             rep,rep1,rep2,projbest.node,splitCutoff.node,
                             G,cut.class.keep,TOL,...) {
      maxG<-500
      cut.class<-median(origY)
      origclass<-ifelse(origY<cut.class,1,2)
      if(length(table(origclass))==1){
         origclass<- ifelse(origY<=cut.class,1,2)
      } 
      n<-nrow(origdata)
      p<-ncol(origdata)
      g<-table(origclass)
      if(length(Tree.Struct)==0) {
         Tree.Struct<-matrix(1:(2*maxG-1),ncol=1)
         Tree.Struct<-cbind(Tree.Struct,0,0,0,0,0,0)
      }
      CV.test<-ifelse(sum(g < max((p/4),5))!=0,TOL.CV/2,
                    ifelse(abs(mean(origY))<1.0e-5,sd(origY),
                           sd(origY)/abs(mean(origY))))
      if(is.null(DEPTH)&(length(table(origclass))<=1|CV.test<TOL.CV)){
         G<-G+1
         Tree.Struct[id, 3]<-G
         Tree.Struct[id,7]<-sd(origY)
         Tree.Struct[id,6]<-mean(origY)
         depth<-depth
         cut.class.keep<-sort(cut.class.keep)
         return(list(Tree.Struct=Tree.Struct,projbest.node=projbest.node,
                     depth=depth,splitCutoff.node=splitCutoff.node,rep=rep,
                     rep1=rep1,rep2=rep2,G=G,cut.class.keep=cut.class.keep))
      } else if(ifelse(!is.null(DEPTH),depth>=DEPTH,FALSE)) {
         G<-G+1
         Tree.Struct[id, 3]<-G
         Tree.Struct[id,7]<-sd(origY)
         Tree.Struct[id,6]<-mean(origY)
         cut.class.keep<-sort(cut.class.keep)
         return(list(Tree.Struct=Tree.Struct,projbest.node=projbest.node, 
                     depth=depth,splitCutoff.node=splitCutoff.node,rep=rep,
                     rep1=rep1,rep2=rep2,G=G,cut.class.keep=cut.class.keep))
      } else {
         cut.class.keep<-c(cut.class.keep,cut.class)
         depth<-depth+1;
         Tree.Struct[id,2]<-rep1
         rep1<-rep1+1
         Tree.Struct[id,3]<-rep1
         rep1<-rep1+1
         Tree.Struct[id,4]<-rep2
         rep2<-rep2+1
         a<-Find.proj(origclass,origdata,PPmethod,weight,lambda)
         splitCutoff.node<-rbind(splitCutoff.node,a$C)
         Tree.Struct[id,5]<-a$Index
         projbest.node<-rbind(projbest.node,a$Alpha)
         t.data<-origdata[a$IOindexL,,drop=FALSE]
         t.Y<-origY[a$IOindexL]   
           
         b<-Tree.construct(t.Y,as.matrix(t.data),Tree.Struct,
                           Tree.Struct[id,2],depth,rep,rep1,rep2,
                           projbest.node,splitCutoff.node,G,
                           cut.class.keep,TOL.CV)
         G<-b$G
         cut.class.keep<-b$cut.class.keep
         Tree.Struct<-b$Tree.Struct
         projbest.node<-b$projbest.node
         splitCutoff.node<-b$splitCutoff.node
         rep<-b$rep
         rep1<-b$rep1
         rep2<-b$rep2
         t.data<-origdata[a$IOindexR,,drop=FALSE]
         t.Y<-origY[a$IOindexR]
         n<-nrow(t.data)
         b<-Tree.construct(t.Y,as.matrix(t.data),Tree.Struct, 
                           Tree.Struct[id,3],depth,rep,rep1,rep2, 
                           projbest.node,splitCutoff.node,G,
                           cut.class.keep,TOL.CV)
         Tree.Struct<-b$Tree.Struct
         depth<-b$depth
         cut.class.keep<-b$cut.class.keep 
         projbest.node<-b$projbest.node
         splitCutoff.node<-b$splitCutoff.node
         rep<-b$rep
         rep1<-b$rep1
         rep2<-b$rep2
         G<-b$G           
      }   
      cut.class.keep<-sort(cut.class.keep)
      return(list(Tree.Struct=Tree.Struct,projbest.node=projbest.node, 
                  splitCutoff.node=splitCutoff.node,depth=depth,rep=rep,
                  rep1=rep1,rep2=rep2,G=G,cut.class.keep=cut.class.keep))
   }

   origX<-as.matrix(origX)  
   splitCutoff.node<-NULL
   projbest.node<-NULL
   Tree.Struct<-NULL
   id<-1
   rep1<-2
   rep2<-1
   rep<-1
   Tree.final<-Tree.construct(origY,origX,Tree.Struct,id,0,rep,rep1,rep2, 
                              projbest.node,splitCutoff.node,0,NULL,TOL.CV)
   Tree.Struct<-Tree.final$Tree.Struct
   Tree.Struct<-Tree.Struct[-which(Tree.Struct[,3]==0),,drop=FALSE]
   origclass<-rep(0,length(origY))
   g<-length(Tree.final$cut.class.keep)+1 
   for(i in 1:(g-1)){
      sel.id<-which(origclass==0&origY<Tree.final$cut.class.keep[i])
      if(length(sel.id)==0)
         sel.id<-which(origclass==0&origY<=Tree.final$cut.class.keep[i])
      origclass[sel.id]<-i
   }    
   origclass[origclass==0]<-g   
   colnames(Tree.Struct)<-c("id","L.node.ID","R.F.node.ID",
                            "Coef.ID","Index","mean","sd")    
   projbest.node<-Tree.final$projbest.node
   splitCutoff.node<-Tree.final$splitCutoff.node     
   treeobj<-list(Tree.Struct=Tree.Struct,projbest.node=projbest.node, 
                 splitCutoff.node=splitCutoff.node,origclass=origclass,
                 origdata= origX)
   origclass<-factor(origclass,levels=1:g)
   if(final.rule==2){
      mean.G<-tapply(origY,origclass,median)
   } else{
      mean.G<-tapply(origY,origclass,mean)
   }  
   sd.G<-tapply(origY,origclass,sd) 
   predict.class<-origclass
   TS<-Tree.Struct
   p<-ncol(origX)
   selP<-ifelse(is.null(selP),max(min(p,2),round(p*0.2)),selP)
   g<-length(table(predict.class))
   if(final.rule<=2){
      predict.Y<-mean.G[as.numeric(predict.class)]
      coef.G<-NULL
   } else{
      coef.G<-matrix(0,ncol=(1+ncol(origX)),nrow=g)
      predict.Y<-rep(0,length(predict.class))
      for(i in 1:g){
         sel.id<-which(predict.class==i) 
         if(length(sel.id)!=0){
            cor.list<-sort.list(abs(suppressWarnings(
                           cor(origY[sel.id],origX[sel.id,]))),decreasing=TRUE)
            if(final.rule==3){
               str.id<-TS[c(which(TS[,2]==which(TS[,2]==0 &TS[,3]==i)),
                     which(TS[,2]!=0&TS[,3]==which(TS[,2]==0 &TS[,3]==i))),4]
               temp.G<-projbest.node[str.id,]
               proj.data<-as.matrix(origX)%*%matrix(temp.G)
               sel.data<-data.frame(Y=origY[sel.id],X=proj.data[sel.id,1])
               temp.coef<-coef(lm(Y~.,data=sel.data))        
               temp.G<-c(temp.coef[1],temp.G*temp.coef[2])
            } else if(final.rule==4){
               sel.data<-data.frame(Y=origY[sel.id],X=origX[sel.id,cor.list])
               temp.G<-rep(0,ncol(origX)+1)
               temp.G[c(1,cor.list+1)]<-coef(lm(Y~.,data=sel.data))                  
            } else if(final.rule==5){
               sel.data<-data.frame(Y=origY[sel.id],
                                    X=origX[sel.id,cor.list[1:selP]])
               temp.G<-rep(0,ncol(origX)+1)
               temp.G[c(1,cor.list[1:selP]+1)]<-coef(lm(Y~.,data=sel.data))        
            }
         } else{
            temp.G<-rep(0,ncol(coef.G))
         }
         if(!is.null(temp.G)) 
            temp.G[is.na(temp.G)]<-0
         predict.Y[sel.id]<-as.matrix(cbind(rep(1,length(sel.id)),
                                         origX[sel.id,]))%*%matrix(temp.G)[,1]
         coef.G[i,]<-temp.G
      }
      colnames(coef.G)<-c("intercept",colnames(origX))
   }   
   MSE<-mean((origY-predict.Y)^2)
   colnames(splitCutoff.node)<-paste("Rule",1:9,sep="")
   treeobj<-list(Tree.Struct=Tree.Struct,projbest.node=projbest.node, 
                 splitCutoff.node=splitCutoff.node,origclass=origclass,
                 origdata=origX)
   class(treeobj)<-append(class(treeobj),"PPtreeclass")
   regtreeobj<-list(Tree.result=treeobj,MSE=MSE,mean.G=mean.G,
                    sd.G=sd.G,coef.G=coef.G,origY=origY)
   class(regtreeobj)<-append(class(regtreeobj),"PPtreereg") 
   return(regtreeobj)
}