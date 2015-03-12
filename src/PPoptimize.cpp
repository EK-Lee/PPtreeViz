#include <RcppArmadillo.h> 
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp;

List VecSort(NumericVector IDdata,IntegerVector Auxdata) {
   NumericVector sortX = clone(IDdata);
   IntegerVector sortY = clone(Auxdata);
   
   int n = sortX.size();
   
   for(int i=0; i<(n-1); i++){
      for(int j=(i+1); j<n; j++){
          if(sortX(j) < sortX(i)){
            double temp = sortX(i);
            sortX(i) = sortX(j);
            sortX(j)=temp;
            int tempA = sortY(i);
            sortY(i) = sortY(j);
            sortY(j) = tempA;
          }
      }
   }
   
   return List::create(_["sortID"]=sortX,_["sortAux"]=sortY);
}

// [[Rcpp::export]]
NumericMatrix NormalizeD(NumericMatrix rawdata){

   int n=rawdata.nrow(),p=rawdata.ncol();
   
   NumericMatrix normdata(n,p);
   NumericVector vmean(p),vs(p);
   
   for(int k=0;k<p;k++){
      for(int i=0; i<n; i++){
        vmean(k) += rawdata(i,k);
        vs(k) += pow(rawdata(i,k),2);
     }
     vmean(k) /=n;
     vs(k) = pow((vs(k) - n*pow(vmean(k),2))/(n-1),0.5);
   }
   for(int k=0;k<p;k++){
      for(int i=0; i<n; i++){
        normdata(i,k) = (rawdata(i,k)-vmean(k))/vs(k);
     }
   }
   return normdata;
 
}


NumericMatrix NormalizeProj(NumericMatrix proj){

   int p=proj.nrow(),q=proj.ncol();
   
   NumericMatrix normproj(p,q);
   double ss=0;
   for(int i=0; i<p; i++){
       ss += proj(i,0)*proj(i,0);
   }
   for(int i=0; i<p; i++){
       normproj(i,0) = proj(i,0)/sqrt(ss);
   }  
   
   if(q>1){
      for(int j=1; j<q; j++){
        double temp1 = 0,temp2 = 0;
        for(int i=0; i<p; i++){
            temp1 += normproj(i,j)*normproj(i,j-1);
            temp2 += normproj(i,j-1)*normproj(i,j-1);
        }
        for(int i=0; i<p; i++){
           normproj(i,j) = normproj(i,j)-normproj(i,j-1)*temp1/temp2; 
        }
      }
   }
   return normproj;
 
}
// [[Rcpp::export]]
NumericVector NormalizeProjV(NumericVector proj){

   int p=proj.size();
   
   NumericVector normproj(p);
   
      double ss=0,s=0;
      for(int i=0; i<p; i++){
         ss += proj(i)*proj(i);
         s += proj(i)/p;
      }
      for(int i=0; i<p; i++){
         normproj(i) = (proj(i)-s)/sqrt(ss);
      }   
   return normproj;
 
}
// [[Rcpp::export]]
NumericVector NormalizeProjV1(NumericVector proj){

   int p=proj.size();
   
   NumericVector normproj(p);
   
      double ss=0,s=0;
      for(int i=0; i<p; i++){
         ss += proj(i)*proj(i);
      }
      for(int i=0; i<p; i++){
         normproj(i) = proj(i)/sqrt(ss);
      }   
   return normproj;
 
}
// [[Rcpp::export]]
int FindMaxID(NumericVector X){

   int p=X.size();
   
      double maxtemp = X(0); int maxID=0;
      for(int i=1; i<p; i++){
         if(maxtemp < X(i)){
            maxtemp = X(i);
            maxID = i;
           }
      }
  
   return maxID;
 
}


// [[Rcpp::export]]

double LDAindex1(IntegerVector projclass, NumericMatrix projdata){

   double index;
   int n=projdata.nrow(),p=projdata.ncol();
   Environment base("package:base");
   Function table = base["table"];
   NumericVector gn=table(projclass);
   int g=gn.size();  
   NumericMatrix W(p,p),WB(p,p),gsum(p,g);
   NumericVector allmean(p);
   
/*   projdata = NormalizeD(projdata);*/
   for(int i=0; i<n; i++){
      for(int k=0;k<p;k++){
        allmean(k) += projdata(i,k)/n;
        gsum(k,(projclass(i)-1)) += projdata(i,k);
     }
   }
   
   for(int i=0; i<n; i++){
         int l = projclass[i]-1;
         for(int j1=0;j1<p;j1++){
            for(int j2=0; j2<=j1; j2++) {
                W(j1,j2) += (projdata(i,j1)-gsum(j1,l)/gn(l))*(projdata(i,j2)-gsum(j2,l)/gn(l));
                W(j2,j1) = W(j1,j2);
                double temp = (projdata(i,j1)-gsum(j1,l)/gn(l))*(projdata(i,j2)-gsum(j2,l)/gn(l))+
                               (gsum(j1,l)/gn(l)-allmean(j1))*(gsum(j2,l)/gn(l)-allmean(j2));
                WB(j1,j2) += temp;
                WB(j2,j1) = WB(j1,j2);
            }
         }
   }
   Function det = base["det"];
   index =1.0-as<double>(det(wrap(W)))/as<double>(det(wrap(WB)));
   return index;
}


// [[Rcpp::export]]
double LDAindex2(IntegerVector projclass, NumericMatrix projdata){

   double index;
   int n=projdata.nrow(),p=projdata.ncol();
   Environment base("package:base");
   Function table = base["table"];
   NumericVector gn=table(projclass);
   int g=gn.size();   
   NumericMatrix W(p,p),WB(p,p),gsum(p,g);
   NumericVector allmean(p);
   
 /*  projdata = NormalizeD(projdata);*/
   for(int i=0; i<n; i++){
      for(int k=0;k<p;k++){
        allmean(k) += projdata(i,k)/n;
        gsum(k,(projclass(i)-1)) += projdata(i,k);
     }
   }
   double gn1 = n/g;
   for(int i=0; i<n; i++){
         int l = projclass[i]-1;
         for(int j1=0;j1<p;j1++){
            for(int j2=0; j2<=j1; j2++) {
                W(j1,j2) += (projdata(i,j1)-gsum(j1,l)/gn(l))*(projdata(i,j2)-gsum(j2,l)/gn(l))/gn(l);
                W(j2,j1) = W(j1,j2);
                double temp = (projdata(i,j1)-gsum(j1,l)/gn(l))*(projdata(i,j2)-gsum(j2,l)/gn(l))/gn(l)+
                               ((gsum(j1,l)/gn(l)-allmean(j1))*(gsum(j2,l)/gn(l)-allmean(j2)))/gn(l);
                WB(j1,j2) += temp;
                WB(j2,j1) = WB(j1,j2);
            }
         }
   }
   Function det = base["det"];
   index =1.0-as<double>(det(wrap(W)))/as<double>(det(wrap(WB)));
    
   return index;
}


// [[Rcpp::export]]
double Lpindex1(IntegerVector projclass, NumericMatrix projdata,int r){
 
   double index,B=0,W=0;
   int n=projdata.nrow(),p=projdata.ncol();
   Environment base("package:base");
   Function table = base["table"];
   NumericVector gn=table(projclass);
   int g=gn.size();   
   NumericMatrix gsum(p,g);
   NumericVector allmean(p);
   
   projdata = NormalizeD(projdata);
   for(int i=0; i<n; i++){
      for(int k=0;k<p;k++){
        allmean(k) = allmean(k)+ projdata(i,k)/n;
        gsum(k,(projclass(i)-1)) += projdata(i,k);
     }
   }
   for(int i=0; i<n; i++){
         int l = projclass(i)-1;
         for(int j=0;j<p;j++){
                W = W+ pow(pow((projdata(i,j)-gsum(j,l)/gn(l)),2*r),0.5);
                B = B+ pow(pow((gsum(j,l)/gn(l)-allmean(j)),2*r),0.5);
        }
   }
   
    W = pow(W,1/r); 
    double WB = pow(W+B,1/r);
    index = 1-W/WB;
   return index;
}

// [[Rcpp::export]]
double Lpindex2(IntegerVector projclass, NumericMatrix projdata,int r){
 
   double index,B=0,W=0;
   int n=projdata.nrow(),p=projdata.ncol();
   Environment base("package:base");
   Function table = base["table"];
   NumericVector gn=table(projclass);
   int g=gn.size();   
   NumericMatrix gsum(p,g);
   NumericVector allmean(p);
   
   projdata = NormalizeD(projdata);
   for(int i=0; i<n; i++){
      for(int k=0;k<p;k++){
        allmean(k) = allmean(k)+ projdata(i,k)/n;
        gsum(k,(projclass(i)-1)) += projdata(i,k);
     }
   }
   for(int i=0; i<n; i++){
         int l = projclass(i)-1;
         for(int j=0;j<p;j++){
                W = W+ pow(pow((projdata(i,j)-gsum(j,l)/gn(l)),2*r),0.5)/gn(l);
                B = B+ pow(pow((gsum(j,l)/gn(l)-allmean(j)),2*r),0.5)/gn(l);

        }
   }
   
    W = pow(W,1/r); 
    double WB = pow(W+B,1/r);
    index = 1-W/WB;
   return index;
}

// [[Rcpp::export]]
double PDAindex1(IntegerVector projclass, NumericMatrix projdata,double lambda){

   double index;
   int n=projdata.nrow(),p=projdata.ncol();
   Environment base("package:base");
   Function table = base["table"];
   NumericVector gn=table(projclass);
   int g=gn.size();   
   NumericMatrix W(p,p),WB(p,p),gsum(p,g),normdata(n,p);
   
   normdata = NormalizeD(projdata);
   for(int i=0; i<n; i++){
      for(int k=0;k<p;k++){
        gsum(k,(projclass(i)-1)) += normdata(i,k);
     }
   }
   
   for(int i=0; i<n; i++){
         int l = projclass[i]-1;
         for(int j1=0;j1<p;j1++){
            for(int j2=0; j2<=j1; j2++) {
                W(j1,j2) += (1-lambda)*(normdata(i,j1)-gsum(j1,l)/gn(l))*(normdata(i,j2)-gsum(j2,l)/gn(l));
                
                W(j2,j1) = W(j1,j2);
                double temp = (1-lambda)*(normdata(i,j1)-gsum(j1,l)/gn(l))*(normdata(i,j2)-gsum(j2,l)/gn(l))+
                               (gsum(j1,l)/gn(l))*(gsum(j2,l)/gn(l));
                WB(j1,j2) += temp;
                WB(j2,j1) = WB(j1,j2);
            }
         }
   }
   for(int j=0; j<p; j++){
      W(j,j) = W(j,j)+ n*lambda;
      WB(j,j) = WB(j,j)+ n*lambda;  
   }
   Function det = base["det"];
   index =1.0-as<double>(det(wrap(W)))/as<double>(det(wrap(WB)));
    
   return index;
}


// [[Rcpp::export]]
double PDAindex2(IntegerVector projclass, NumericMatrix projdata,double lambda){

   double index;
   int n=projdata.nrow(),p=projdata.ncol();
   Environment base("package:base");
   Function table = base["table"];
   NumericVector gn=table(projclass);
   int g=gn.size();   
   NumericMatrix W(p,p),WB(p,p),gsum(p,g),normdata(n,p);

   
   normdata = NormalizeD(projdata);
   for(int i=0; i<n; i++){
      for(int k=0;k<p;k++){
        gsum(k,(projclass(i)-1)) += normdata(i,k);
     }
   }
      double gn1 = n/g;
   for(int i=0; i<n; i++){
         int l = projclass[i]-1;
         for(int j1=0;j1<p;j1++){
            for(int j2=0; j2<=j1; j2++) {
                W(j1,j2) += (1-lambda)*(normdata(i,j1)-gsum(j1,l)/gn(l))*(normdata(i,j2)-gsum(j2,l)/gn(l))/gn(l)*gn1;                
                W(j2,j1) = W(j1,j2);
                double temp = (1-lambda)*(normdata(i,j1)-gsum(j1,l)/gn(l))*(normdata(i,j2)-gsum(j2,l)/gn(l))/gn(l)*gn1+
                               (1-lambda)*(gsum(j1,l)/gn(l))*(gsum(j2,l)/gn(l))/gn(l)*gn1;
                WB(j1,j2) += temp;
                WB(j2,j1) = WB(j1,j2);
            }
         }
   }
   for(int j=0; j<p; j++){
      W(j,j) = W(j,j)+ n*lambda;
      WB(j,j) = WB(j,j)+ n*lambda;  
   }

   Function det = base["det"];
   index =1.0-as<double>(det(wrap(W)))/as<double>(det(wrap(WB)));
    
   return index;
}



// [[Rcpp::export]]
double GINIindex1D(IntegerVector projclass, NumericVector projdata){
 
   Environment base("package:base");
   Function table = base["table"];
   NumericVector gn=table(projclass);
   int g=gn.size(), n = projclass.size();
   double n1,n2;
   double index=0,tempindex;
   List VecSortdata =  VecSort(projdata,projclass);
   NumericVector sortdata = as<NumericVector>(VecSortdata["sortID"]);
   IntegerVector sortclass = as<IntegerVector>(VecSortdata["sortAux"]);
   IntegerVector part1,part2,temptable1,temptable2;

   for(int j=0; j<g; j++){
      index += (gn(j)/n)*(1-gn(j)/n);
   }  
   
   for(int i=1; i<(n-1); i++){  
       part1 = sortclass[sortdata<=sortdata[i]];
       part2 = sortclass[sortdata>sortdata[i]];
       n1 = part1.size();
       n2 = part2.size();
       temptable1 = table(part1); int g1 = temptable1.size();  tempindex=0;
       temptable2 = table(part2); int g2 = temptable2.size(); 
       for(int j=0; j<g1; j++){
         tempindex += ((n1)/n)*(temptable1(j)/n1)*(1-temptable1(j)/n1);
       }
       for(int j=0; j<g2; j++){         
         tempindex += ((n2)/n)*(temptable2(j)/n2)*(1-temptable2(j)/n2);        
       } 
      if (tempindex < index) index = tempindex; 
            
   }   
   return 1-index*2;
}


// [[Rcpp::export]]

double ENTROPYindex1D(IntegerVector projclass, NumericVector projdata){
 
   Environment base("package:base");
   Function table = base["table"];
   NumericVector gn=table(projclass);
   int g=gn.size(), n = projclass.size();
   double n1,n2;
   double index=0,tempindex;
   List VecSortdata =  VecSort(projdata,projclass);
   NumericVector sortdata = as<NumericVector>(VecSortdata["sortID"]);
   IntegerVector sortclass = as<IntegerVector>(VecSortdata["sortAux"]);
   IntegerVector part1,part2,temptable1,temptable2;
   double maxindex = 0;

   for(int j=0; j<g; j++){
      if(gn(j)!=0)
      { index -= (gn(j)/n)*log(gn(j)/n);
      }
    }  
   maxindex = -log(1.0/g);
   for(int i=1; i<(n-1); i++){  
       part1 = sortclass[sortdata<=sortdata[i]];
       part2 = sortclass[sortdata>sortdata[i]];
       n1 = part1.size();
       n2 = part2.size();
       temptable1 = table(part1); int g1 = temptable1.size();  tempindex=0;
       temptable2 = table(part2); int g2 = temptable2.size(); 
       for(int j=0; j<g1; j++){
         if(temptable1(j)!=0)
         {  tempindex -= ((n1)/n)*(temptable1(j)/n1)*log(temptable1(j)/n1);
         }
       }
       for(int j=0; j<g2; j++){         
         if(temptable2(j)!=0)
         {  tempindex -= ((n2)/n)*(temptable2(j)/n2)*log(temptable2(j)/n2);        
         }
       } 
      if (tempindex < index) index = tempindex; 
            
   }   
   return (maxindex-index)/maxindex;
}

// [[Rcpp::export]]
List PPoptimizeqD(IntegerVector origclass, NumericMatrix origdata,std::string method,int q,
             int r=1,double lambda=0,double TOL=0.00001, int maxiter=5000){

   int n=origdata.nrow(),p=origdata.ncol();
   Environment base("package:base");
   Function table = base["table"];
   
   NumericMatrix projbest(p,q);
   for(int k=0; k<q; k++){
         projbest(_,k)= rnorm(p);
      }
   NumericMatrix projdata(n,q);
   projbest = NormalizeProj(projbest);  
   for(int i=0; i<n; i++){
      for(int k=0; k<q; k++){
         projdata(i,k)=0;
         for(int j=0; j<p; j++){
            projdata(i,k) += origdata(i,j)*projbest(j,k);
         }
      }
   }
  
  double indexbest,newindex;
  if(method=="LDA1"){
   indexbest = LDAindex1(origclass,projdata);
  } else if(method=="LDA2"){
   indexbest = LDAindex2(origclass,projdata);
  } else if(method=="Lp"){
   indexbest = Lpindex1(origclass,projdata,r);
  } else if(method=="PDA1"){
   indexbest = PDAindex1(origclass,projdata,lambda);
  }  else if(method=="PDA2"){
   indexbest = PDAindex2(origclass,projdata,lambda);
  }  else if(method=="GINI"){
   indexbest = GINIindex1D(origclass,projdata);
  }  else if(method=="ENTROPY"){
   indexbest = ENTROPYindex1D(origclass,projdata);
  } 
   double temp=1;
   int kk=0;
   
   NumericVector indexkeep(maxiter);
   while(temp > TOL && kk < maxiter){
      NumericMatrix projnew(p,q);
      for(int k=0; k<q; k++){
         projnew(_,k)= rnorm(p);
      }
      projnew = NormalizeProj(projnew);  
        
      for(int i=0; i<n; i++){
         for(int k=0; k<q; k++){
            projdata(i,k)=0;
            for(int j=0; j<p; j++){
               projdata(i,k) += origdata(i,j)*projnew(j,k);
            }
         }
      }

      if(method=="LDA1"){
         newindex = LDAindex1(origclass,projdata);
      } else if(method=="LDA2"){
         newindex = LDAindex2(origclass,projdata);
      } else if(method=="Lp"){
         newindex = Lpindex1(origclass,projdata,r);
      } else if(method=="PDA1"){
         newindex = PDAindex1(origclass,projdata,lambda);
      } else if(method=="PDA2"){
         newindex = PDAindex2(origclass,projdata,lambda);
      } else if(method=="GINI"){
         newindex = GINIindex1D(origclass,projdata);
      } else if(method=="ENTROPY"){
         newindex = ENTROPYindex1D(origclass,projdata);
      } 
 
      if(newindex > indexbest){
          temp = newindex - indexbest;
          indexbest = newindex; 
          projbest = projnew;     

      }
          indexkeep[kk]=newindex;      
      kk++;    
   }
         indexkeep[kk-1]=100;
   return Rcpp::List::create(Rcpp::Named("indexkeep") = indexkeep,
                          Rcpp::Named("optproj") = projbest,
                          Rcpp::Named("optindex") = indexbest);
}

// [[Rcpp::export]]
List PPoptimizeAnnealqD(IntegerVector origclass, NumericMatrix origdata,std::string method,int q,
                        int r=1,double lambda=0, double TOL=0.0001, int maxiter=50000,double energy=0.01, double cooling=0.999){

   int n=origdata.nrow(),p=origdata.ncol();
   Environment base("package:base");
   Function table = base["table"];
   double cool=cooling;
   GetRNGstate();
   NumericMatrix projbest(p,q);
   for(int k=0; k<q; k++){
         projbest(_,k)= rnorm(p);
      }
   projbest = NormalizeProj(projbest);  
   NumericMatrix projdata(n,q);  
  
   for(int i=0; i<n; i++){
      for(int k=0; k<q; k++){
         projdata(i,k)=0;
         for(int j=0; j<p; j++){
            projdata(i,k) += origdata(i,j)*projbest(j,k);
         }
      }
   }
  
   double indexbest,newindex;
   if(method=="LDA1"){
      indexbest = LDAindex1(origclass,projdata);
   } else if(method=="LDA2"){
      indexbest = LDAindex2(origclass,projdata);
   } else if(method=="Lp"){
      indexbest = Lpindex1(origclass,projdata,r);
   } else if(method=="PDA1"){
      indexbest = PDAindex1(origclass,projdata,lambda);
   } else if(method=="PDA2"){
      indexbest = PDAindex2(origclass,projdata,lambda);
   } 
   
   double temp=1;
   int kk=0;
   
   NumericVector indexkeep(maxiter);
   double diff = 100, keepdiff=diff;
   while(fabs(diff)>TOL && kk < maxiter){
      double tempp = energy/log(kk+2);
      if(kk>1000) {
      temp = temp*cool;
      } else {
      temp = temp*cool*cool;
      } 
      
      NumericMatrix projnew(p,q);
      for(int k=0; k<q; k++){
         projnew(_,k)= temp*rnorm(p)+projbest(_,k);/*/(kk+1)+kk/(kk+1)*projbest(_,k); */
      }
      projnew = NormalizeProj(projnew);  
        
      for(int i=0; i<n; i++){
         for(int k=0; k<q; k++){
            projdata(i,k)=0;
            for(int j=0; j<p; j++){
               projdata(i,k) += origdata(i,j)*projnew(j,k);
            }
         }
      }
      if(method=="LDA1"){
         newindex = LDAindex1(origclass,projdata);
      } else if(method=="LDA2"){
         newindex = LDAindex2(origclass,projdata);
      } else if(method=="Lp"){
         newindex = Lpindex1(origclass,projdata,r);
      } else if(method=="PDA1"){
         newindex = PDAindex1(origclass,projdata,lambda);
      } else if(method=="PDA2"){
         newindex = PDAindex2(origclass,projdata,lambda);
      } 
      
      NumericVector prob = runif(1);
      double difft = newindex - indexbest;
      double e = exp(difft/tempp);
      if( e>1){
          for(int i=0; i<p; i++){
             for(int j=0; j<q; j++){
                projbest(i,j) = projnew(i,j);
             }
          }
          indexbest = newindex;    
          diff = difft;
      kk++;  indexkeep[kk]=newindex; 
      } else if(prob[0] < e && difft>energy){
          for(int i=0; i<p; i++){
             for(int j=0; j<q; j++){
                projbest(i,j) = projnew(i,j);
             }
          }
          indexbest = newindex;    
          diff = difft;
      kk++;  indexkeep[kk]=newindex; 
      }
/*      printf("kk=%d : tempp = %f, prob = %f, difft = %f, e = %f\n, TOL = %f, abs.diff = %f, maxiter=%d\n",kk,tempp,prob[0],difft,e,TOL,fabs(diff),maxiter);               */
   }
/*         indexkeep[kk-1]=100;*/
   PutRNGstate();
   return Rcpp::List::create(Rcpp::Named("indexkeep") = indexkeep,
                          Rcpp::Named("optproj") = projbest,
                          Rcpp::Named("optindex") = indexbest);
}


// [[Rcpp::export]]
List PPoptimizeAnnealorig(IntegerVector origclass, NumericMatrix origdata,std::string method,int q,
                        int r=1,double lambda=0, double TOL=0.0001, int maxiter=50000,double energy=0.01, double cooling=0.999,
                        double tempstart=1){

   int n=origdata.nrow(),p=origdata.ncol();
   Environment base("package:base");
   Function table = base["table"];
   double cool=cooling;
   double tempend = 0.001,temp = tempstart;
   int maxproj = 1000;
   GetRNGstate();
   NumericMatrix projbest(p,q);
   for(int k=0; k<q; k++){
         projbest(_,k)= rnorm(p);
      }
   projbest = NormalizeProj(projbest);  
   NumericMatrix projdata(n,q);  
  
   for(int i=0; i<n; i++){
      for(int k=0; k<q; k++){
         projdata(i,k)=0;
         for(int j=0; j<p; j++){
            projdata(i,k) += origdata(i,j)*projbest(j,k);
         }
      }
   }
  
   double indexbest,newindex;
   if(method=="LDA1"){
      indexbest = LDAindex1(origclass,projdata);
   } else if(method=="LDA2"){
      indexbest = LDAindex2(origclass,projdata);
   } else if(method=="Lp"){
      indexbest = Lpindex1(origclass,projdata,r);
   } else if(method=="PDA1"){
      indexbest = PDAindex1(origclass,projdata,lambda);
   } else if(method=="PDA2"){
      indexbest = PDAindex2(origclass,projdata,lambda);
   } 
   
   
   int kk=1;
   
   NumericVector indexkeep(maxiter);
   double diff = 100, keepdiff=diff;
   while((temp>0.001 || fabs(diff)>(energy/1000000)) && kk <= maxiter){
      double tempp = energy/log(kk+1)/10000;
      temp = temp*cool;
      NumericMatrix projnew(p,q);
      for(int k=0; k<q; k++){
         projnew(_,k)= temp*rnorm(p)+projbest(_,k);/*/(kk+1)+kk/(kk+1)*projbest(_,k); */
      }
      projnew = NormalizeProj(projnew);  
        
      for(int i=0; i<n; i++){
         for(int k=0; k<q; k++){
            projdata(i,k)=0;
            for(int j=0; j<p; j++){
               projdata(i,k) += origdata(i,j)*projnew(j,k);
            }
         }
      }
      if(method=="LDA1"){
         newindex = LDAindex1(origclass,projdata);
      } else if(method=="LDA2"){
         newindex = LDAindex2(origclass,projdata);
      } else if(method=="Lp"){
         newindex = Lpindex1(origclass,projdata,r);
      } else if(method=="PDA1"){
         newindex = PDAindex1(origclass,projdata,lambda);
      } else if(method=="PDA2"){
         newindex = PDAindex2(origclass,projdata,lambda);
      } 
      
      NumericVector prob = runif(1);
      double diff = newindex - indexbest;
      double e = exp(diff/tempp);
      if(prob[0] < e){
          for(int i=0; i<p; i++){
             for(int j=0; j<q; j++){
                projbest(i,j) = projnew(i,j);
             }
          }
          indexbest = newindex;    
          kk++;  indexkeep[kk]=newindex; 
      }
 /*     printf("kk=%d : indexbest = %f, \ntemp = %f, prob = %f, difft = %f, e = %f\n, TOL = %f, abs.diff = %f, maxiter=%d\n",
              kk,indexbest,temp,prob[0],diff,e,TOL,fabs(diff),maxiter);       */        
   }
/*         indexkeep[kk-1]=100;*/
   PutRNGstate();
   return Rcpp::List::create(Rcpp::Named("indexkeep") = indexkeep,
                          Rcpp::Named("optproj") = projbest,
                          Rcpp::Named("optindex") = indexbest);
}

