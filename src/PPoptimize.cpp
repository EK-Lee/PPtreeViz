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
            temp1 += proj(i,j)*normproj(i,j-1);
            temp2 += normproj(i,j-1)*normproj(i,j-1);
        }
        for(int i=0; i<p; i++){
           normproj(i,j) = proj(i,j)-normproj(i,j-1)*temp1/temp2; 
        }
      }
   }
   return normproj;
 
}

// [[Rcpp::export]]

double LDAindex(IntegerVector projclass, NumericMatrix projdata, bool weight=true){

   double index;
   int n=projdata.nrow(),p=projdata.ncol();
   Environment base("package:base");
   Function table=base["table"];
   NumericVector gn=table(projclass);
   int g=gn.size();  
   NumericMatrix W(p,p),WB(p,p),gsum(p,g);
   NumericVector allmean(p);
   
   for(int i=0;i<n;i++){
      for(int k=0;k<p;k++){
         allmean(k) += projdata(i,k)/n;
         gsum(k,(projclass(i)-1)) += projdata(i,k);
      }
   }

   for(int i=0; i<n; i++){
      int l = projclass[i]-1;
      double gn1;
      if(weight){
         gn1 = gn(l);
      } else{
         gn1 = n/g; 
      }
      for(int j1=0;j1<p;j1++){
         for(int j2=0;j2<=j1;j2++) {
            W(j1,j2) += ((projdata(i,j1)-gsum(j1,l)/gn(l))*
                         (projdata(i,j2)-gsum(j2,l)/gn(l)))/gn(l)*gn1;
            W(j2,j1) = W(j1,j2);
            double temp = ((projdata(i,j1)-gsum(j1,l)/gn(l))*
                           (projdata(i,j2)-gsum(j2,l)/gn(l))+
                              (gsum(j1,l)/gn(l)-allmean(j1))*
                              (gsum(j2,l)/gn(l)-allmean(j2)))/gn(l)*gn1;
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

double PDAindex(IntegerVector projclass, NumericMatrix projdata,bool weight=true,double lambda=0.1){

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
      double gn1;
      if(weight){
         gn1 = gn(l);
      } else{
         gn1 = n/g; 
      }
      for(int j1=0;j1<p;j1++){
         for(int j2=0; j2<=j1; j2++) {
            W(j1,j2) += (1-lambda)*(normdata(i,j1)-gsum(j1,l)/gn(l))*
                                   (normdata(i,j2)-gsum(j2,l)/gn(l))/gn(l)*gn1;    
            W(j2,j1) = W(j1,j2);
            double temp = ((1-lambda)*(normdata(i,j1)-gsum(j1,l)/gn(l))*
                                     (normdata(i,j2)-gsum(j2,l)/gn(l))+
                               (gsum(j1,l)/gn(l))*(gsum(j2,l)/gn(l)))/gn(l)*gn1; 
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
double Lrindex(IntegerVector projclass, NumericMatrix projdata,bool weight=true,int r=1){
 
   double index,B=0,W=0;
   int n=projdata.nrow(),p=projdata.ncol();
   Environment base("package:base");
   Function table = base["table"];
   NumericVector gn=table(projclass);
   int g=gn.size();   
   NumericMatrix gsum(p,g);
   NumericVector allmean(p);
   
/*  projdata = NormalizeD(projdata); */
   for(int i=0; i<n; i++){
      for(int k=0;k<p;k++){
         allmean(k) = allmean(k)+ projdata(i,k)/n;
         gsum(k,(projclass(i)-1)) += projdata(i,k);
      }
   }
   for(int i=0; i<n; i++){
      int l = projclass(i)-1;
      double gn1;
      if(weight){
         gn1 = gn(l);
      } else{
         gn1 = n/g; 
      }         
      for(int j=0;j<p;j++){
         W = W+ pow(pow((projdata(i,j)-gsum(j,l)/gn(l)),2*r),0.5)/gn(l)*gn1;
         B = B+ pow(pow((gsum(j,l)/gn(l)-allmean(j)),2*r),0.5)/gn(l)*gn1;
      }
   }
   
   W = pow(W,1.0/r); 
   double WB = pow(W+B,1.0/r);
   index = 1-W/WB;
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

   for(int i=1; i<(n-1); i++){  
      part1 = sortclass[sortdata<=sortdata[i]];
      part2 = sortclass[sortdata>sortdata[i]];
      n1 = part1.size();
      n2 = part2.size();
      temptable1 = table(part1); int g1 = temptable1.size();  tempindex=0;
      temptable2 = table(part2); int g2 = temptable2.size(); 
      for(int j=0; j<g1; j++){
         tempindex += ((n1)/n)*(temptable1(j)/n1)*(1.0-temptable1(j)/n1);
      }
      for(int j=0; j<g2; j++){         
         tempindex += ((n2)/n)*(temptable2(j)/n2)*(1.0-temptable2(j)/n2);        
      } 
      tempindex = g-1.0-g*tempindex;
      if(tempindex > index) 
         index = tempindex;             
   }  

   return index;
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
 
   for(int i=1; i<(n-1); i++){  
      part1 = sortclass[sortdata<=sortdata[i]];
      part2 = sortclass[sortdata>sortdata[i]];
      n1 = part1.size();
      n2 = part2.size();
      temptable1 = table(part1); int g1 = temptable1.size();  tempindex=0;
      temptable2 = table(part2); int g2 = temptable2.size(); 
      for(int j=0; j<g1; j++){
         if(temptable1(j)!=0)
         {  tempindex += ((n1)/n)*(temptable1(j)/n1)*log(temptable1(j)/n1);
         }
      }
      for(int j=0; j<g2; j++){         
         if(temptable2(j)!=0)
         {  tempindex += ((n2)/n)*(temptable2(j)/n2)*log(temptable2(j)/n2);        
         }
      } 
      double maxI = log(2)-log(g);
      if((g/2)*2 != g){
        maxI = -0.5*log((g*g-1.0)/4.0)+1.0/(2.0*g)*log((g-1.0)/(g+1.0));
      }
      
      tempindex = (1+tempindex/log(g))/(1+maxI/log(g));
      if (tempindex > index) index = tempindex;            
   }   

   return index;
}

// [[Rcpp::export]]
List PPopt(IntegerVector origclass, NumericMatrix origdata,int q=1, std::string PPmethod="LDA", 
           bool weight = true, int r=1,double lambda=0.1, 
           double energy=0, double cooling=0.999, double TOL=0.0001, int maxiter=50000){

   int n=origdata.nrow(),p=origdata.ncol();
   Environment base("package:base");
   Function table = base["table"];
   GetRNGstate();
   NumericMatrix projbest(p,q);
   double indexbest=0,newindex=0;
   if((PPmethod=="GINI" || PPmethod=="ENTROPY") && q>1)
   {  /*printf("GINI and ENTROPY PPmethod is only for 1D projection. Use q=1\n"); */
      return Rcpp::List::create(Rcpp::Named("indexbest") = 0,
                             Rcpp::Named("projbest") = 0);
   } else {  
      if(PPmethod=="LDA"){
         indexbest = LDAindex(origclass,origdata,weight);
      } else if(PPmethod=="Lr"){
         indexbest = Lrindex(origclass,origdata,weight,r);
      } else if(PPmethod=="PDA"){
         indexbest = PDAindex(origclass,origdata,weight,lambda);  
      } else if(PPmethod=="GINI"){
         NumericVector tempdata;
         double tempindex=0;
         for(int i=0; i<p; i++){
            tempdata = origdata(_,i);
            tempindex = GINIindex1D(origclass,tempdata);  
            if(indexbest<tempindex)
                indexbest = tempindex;
         }        
      } else if(PPmethod=="ENTROPY"){
         NumericVector tempdata;
         double tempindex=0;
         for(int i=0; i<p; i++){
            tempdata = origdata(_,i);
            tempindex = ENTROPYindex1D(origclass,tempdata);  
            if(indexbest<tempindex)
                indexbest = tempindex;
         } 
      }   
   
      if(energy==0)
         energy = 1-indexbest;
   
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
  
      if(PPmethod=="LDA"){
         indexbest = LDAindex(origclass,projdata,weight);
      } else if(PPmethod=="Lr"){
         indexbest = Lrindex(origclass,projdata,weight,r);
      } else if(PPmethod=="PDA"){
         indexbest = PDAindex(origclass,projdata,weight,lambda);
      } else if(PPmethod=="GINI"){
         NumericVector tempdata;
         tempdata = projdata(_,0);
         indexbest = GINIindex1D(origclass,tempdata);        
      } else if(PPmethod=="ENTROPY"){
         NumericVector tempdata;
         tempdata = projdata(_,0);
         indexbest = ENTROPYindex1D(origclass,tempdata);        
      }
   
      double temp=1;
      int kk=0;
   
      double diff = 100;
      while(fabs(diff)>TOL && kk < maxiter){
         double tempp = energy/log(kk+2);
         if(kk>1000) {
            temp = temp*cooling;
         } else {
            temp = temp*cooling*cooling;
         }   
      
         NumericMatrix projnew(p,q);
         for(int k=0; k<q; k++){
            projnew(_,k)= temp*rnorm(p)+projbest(_,k);
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
      
         if(PPmethod=="LDA"){
            newindex = LDAindex(origclass,projdata,weight);
         } else if(PPmethod=="Lr"){
            newindex = Lrindex(origclass,projdata,weight,r);
         } else if(PPmethod=="PDA"){
            newindex = PDAindex(origclass,projdata,weight,lambda);
         } else if(PPmethod=="GINI"){
            NumericVector tempdata;
            tempdata = projdata(_,0);
            newindex = GINIindex1D(origclass,tempdata);        
         } else if(PPmethod=="ENTROPY"){
            NumericVector tempdata;
            tempdata = projdata(_,0);
            newindex = ENTROPYindex1D(origclass,tempdata);        
         }
      
         NumericVector prob = runif(1);
         double difft = newindex - indexbest;
         double e = exp(difft/tempp);
         if(e>1){
            for(int i=0; i<p; i++){
               for(int j=0; j<q; j++){
                  projbest(i,j) = projnew(i,j);
               }
            }
            indexbest = newindex;    
            diff = difft;
         } else if(prob[0] < e && difft>energy){
            for(int i=0; i<p; i++){
               for(int j=0; j<q; j++){
                  projbest(i,j) = projnew(i,j);
               }
            }
            indexbest = newindex;    
            diff = difft; 
         }
         kk++;
      }
      PutRNGstate();
      return Rcpp::List::create(Rcpp::Named("indexbest") = indexbest,
                             Rcpp::Named("projbest") = projbest,
                             Rcpp::Named("origclass") = origclass,
                             Rcpp::Named("origdata") = origdata);
   }                           
}

