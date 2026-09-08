#include <iostream>
#include <iomanip>
#include <string.h>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
std::vector<float> GetMahalanobis(int dimN, int intdimTest,
                             int dimgene, std::vector<float> &specinpt,
                             std::vector<float> &Siginpt,
                             std::vector<float> &dist,std::vector<float> &meansinp,
                             std::vector<float> &colsv, std::vector<float> &covinpt1,
                             std::vector<float> &covinpt2,
                             std::vector<float> &covinpt3,
                             std::vector<float> &covinpt4,
                             std::vector<float> &sigtest
                             ) {
  vector<float> distances = dist;
  vector<float> MeansVec = meansinp;
  vector<float> columnVec = colsv;
  vector<float> SignatureCellsTest = sigtest;
  vector<float> Covariance_matrix = covinpt1;
  vector<float> InvCovariance_matrix = covinpt2;
  vector<float> Covariance_matrix_cp = covinpt3;
  vector<float> InvCovariance_matrix_cp = covinpt4;
  vector<float> PivotVec = meansinp;
  vector<float> specificCells = specinpt;
  vector<float> SignatureCells = Siginpt;

  cout << "Function called: " << intdimTest << ", " << dimgene << ", " << (dimN+1);
  for(int z=0; z<intdimTest; z++){
  for(int q=0; q<(dimN+1); q++){
    
    for(int g=0; g<dimgene; g++){
      if(q==dimN){
        SignatureCellsTest[g*(dimN+1)+q]=SignatureCells[intdimTest*g + z];
      }
      else{
        SignatureCellsTest[g*(dimN+1)+q]=specificCells[dimN*g + q];
      }
    }
  }
  
  for(int n=0; n<(dimN); n++){
    for(int m=0; m<(dimgene); m++){
      PivotVec[m]=0;
      MeansVec[m] += SignatureCellsTest[m*(dimN+1)+n]/(dimN);
      for(int s=0; s<(dimgene); s++){
        if(s==m){
          InvCovariance_matrix[m*dimgene + s]=1;
        }
        else{
          InvCovariance_matrix[m*dimgene + s]=0;
        }
      }
    }
  }
  
  /*gene by gene matrix*/
  for(int i=0; i<dimgene; i++){
    for(int j=0; j<dimgene; j++){
      float covsum=0;
      for(int k=0; k<(dimN); k++){
        
        covsum = covsum + (SignatureCellsTest[i*(dimN+1)+k]-MeansVec[i])*(SignatureCellsTest[j*(dimN+1)+k]-MeansVec[j]);
        
      }
      Covariance_matrix[i*(dimgene)+j] = covsum/(dimN);
    }
  }
  
  for(int h1=0; h1<dimgene; h1++){
    for(int k1=0; k1<dimgene; k1++){
      
      Covariance_matrix_cp[h1*dimgene + k1]=Covariance_matrix[h1*dimgene + k1];
      InvCovariance_matrix_cp[h1*dimgene + k1]=InvCovariance_matrix[h1*dimgene + k1];
      
      //std::cout << "Covariance_matrix_cp[h1*dimgene + k1]: " << Covariance_matrix_cp[h1*dimgene + k1] << std::endl;
      //std::cout << "InvCovariance_matrix_cp[h1*dimgene + k1]: " << InvCovariance_matrix_cp[h1*dimgene + k1] << std::endl;
    }
  }
  
//Start Gauss Jordan elimination
 //First make the matrix upper triangular
 int h = 0;
    int k=0;
    while((h<dimgene) & (k<dimgene)){
      /*for(int f=h; f<dimgene; f++){
       PivotVec[f]=Covariance_matrix[f*dimgene + k];
       }*/
      
      /* float i_result=PivotVec[h]; //larger element pointer initially pointed to first element
       int i_max=0;
       if(i_result==0){
       for (int i = 1; i < dimgene; i++)
       {
       if (PivotVec[i]!=0)
       {
       i_max = i; //updating the pointer to maximum 
       i_result = PivotVec[i];
       std::cout << "One non zero imax: " << i_max << std::endl;
       }
       }
       
       for(int s=0; s<dimgene; s++){
       Covariance_matrix[i_max*dimgene + s]=Covariance_matrix_cp[h*dimgene + s];
       Covariance_matrix[h*dimgene + s]=Covariance_matrix_cp[i_max*dimgene + s];
       InvCovariance_matrix[i_max*dimgene + s]=InvCovariance_matrix_cp[h*dimgene + s];
       InvCovariance_matrix[h*dimgene + s]=InvCovariance_matrix_cp[i_max*dimgene + s];
       }
       }*/
      
      ////Change the current row to have a leading 1 by dividing by it
      double denom_divide=Covariance_matrix[h*dimgene + k];
      for(int d=0; d<dimgene; d++){
        InvCovariance_matrix[h*dimgene + d] = InvCovariance_matrix[h*dimgene + d]/denom_divide;
        if(d==k){
          Covariance_matrix[h*dimgene + d] = 1;
        }
        else{
          if(denom_divide!=0){
            
            Covariance_matrix[h*dimgene + d] = Covariance_matrix[h*dimgene + d]/denom_divide;
          }
        }
        
        
      }
      
      ////Change the current row to have a leading 1 by dividing by it      
      for(int i=(h+1); i<dimgene; i++){
        double valpre=Covariance_matrix[i*dimgene + k];
        for(int d=0; d<dimgene; d++){
          //if(Covariance_matrix[h*dimgene + d]!=0){
          InvCovariance_matrix[i*dimgene + d] = InvCovariance_matrix[i*dimgene + d] - valpre*InvCovariance_matrix[h*dimgene + d];
        }
        for(int d=0; d<dimgene; d++){
          Covariance_matrix[i*dimgene + d] = Covariance_matrix[i*dimgene + d] - valpre*Covariance_matrix[h*dimgene + d];
          // }
        }
        
      }
      
      
      /* for(int d=0; d<dimgene; d++){
       PivotVec[d]=0.0;
       for(int h1=0; h1<dimgene; h1++){
       Covariance_matrix_cp[h1*dimgene + d]=Covariance_matrix[h1*dimgene + d];
       InvCovariance_matrix_cp[h1*dimgene + d]=InvCovariance_matrix[h1*dimgene + d];
       
       }
      }*/
      h += 1;
      k += 1;
    }
    //Next do the reverse to create a diagonal matrix
    int p = 0;
    while((p<dimgene)){
      int h1=dimgene-(p+1);
      if(h1>0){
        for(int m=1; m<(h1+1); m++){
          int i=h1-m;
          double valpretwo=Covariance_matrix[i*dimgene + (dimgene-(p+1))];
          for(int s=0; s<dimgene; s++){
            i=h1-m;
            int d=dimgene-(s+1);
            InvCovariance_matrix[i*dimgene + d] = InvCovariance_matrix[i*dimgene + d] - valpretwo*InvCovariance_matrix[h1*dimgene + d];
          }
          for(int s=0; s<dimgene; s++){
            i=h1-m;
            int d=dimgene-(s+1);
            Covariance_matrix[i*dimgene + d] = Covariance_matrix[i*dimgene + d] - valpretwo*Covariance_matrix[h1*dimgene + d];
          }
        }
      }
      p += 1;
    }
  
 // for(int h1=0; h1<dimgene; h1++){
   // for(int k1=0; k1<dimgene; k1++){
      
      //std::cout << "Covariance_matrix_cp[h1*dimgene + k1]: " << Covariance_matrix[h1*dimgene + k1] << std::endl;
      //std::cout << "InvCovariance_matrix_cp[h1*dimgene + k1]: " << InvCovariance_matrix[h1*dimgene + k1] << std::endl;
    //}
//  }
  

  
  for(int k=0; k<dimgene; k++){
    columnVec[k]=0;
    for(int i=0; i<dimgene; i++){
      columnVec[k]=columnVec[k]+(SignatureCells[intdimTest*i + z]-MeansVec[i])*InvCovariance_matrix[i*dimgene + k];
    }
  }
  
  for(int k=0; k<dimgene; k++){
    distances[z]=distances[z] + columnVec[k]*(SignatureCells[intdimTest*k + z]-MeansVec[k]);
  }
  
  distances[z]=sqrt(distances[z]);
 /* cout << "distances[z]: " << distances[z];*/
  
  for(int h=0; h<dimgene; h++){
    PivotVec[h]=0;
    MeansVec[h] = 0;
    for(int g=0; g<dimgene; g++){
      Covariance_matrix[h*dimgene + g] = 0;
      if(h==g){
        InvCovariance_matrix[h*dimgene + g] = 1;
      }
      else{
      InvCovariance_matrix[h*dimgene + g] = 0;
      }
    }
    
  }
  
}
  

  MeansVec.clear();
  columnVec.clear();
  SignatureCellsTest.clear();
  Covariance_matrix.clear();
  InvCovariance_matrix.clear();
  Covariance_matrix_cp.clear();
  InvCovariance_matrix_cp.clear();
  PivotVec.clear();
  specificCells.clear();
  SignatureCells.clear();

  return(distances);
}

