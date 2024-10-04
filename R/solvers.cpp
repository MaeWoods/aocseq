// -*- C++ -*-
//===-------------------------- __GetMahalanobis.cpp ----------------------------------===//
//
// Part of the traceseq Project, under the MIT license.
// See https://opensource.org/license/mit/ for license information.
// SPDX short identifier: MIT
//
//===----------------------------------------------------------------------===//


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
vector<float> GetMahalanobis(int dimN, int intdimTest, 
                             int dimgene, std::vector<float> &specinpt, 
                             std::vector<float> &Siginpt,
                             std::vector<float> &dist,std::vector<float> &meansinp,
                             std::vector<float> &colsv, std::vector<float> &covinpt,
                             std::vector<float> &sigtest
                             ) {
  vector<float> distances = dist;
  vector<float> MeansVec = meansinp;
  vector<float> columnVec = colsv;
  vector<float> SignatureCellsTest = sigtest;
  vector<float> Covariance_matrix = covinpt;
  vector<float> InvCovariance_matrix = covinpt;
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
  
  for(int n=0; n<(dimN+1); n++){
    for(int m=0; m<(dimgene); m++){
      MeansVec[m] += SignatureCellsTest[m*(dimN+1)+n]/(dimN+1);
      for(int s=0; s<(dimgene); s++){
        if(s==m){
          InvCovariance_matrix[m*dimgene + s]=1;
        }
      }
    }
  }
  
  /*gene by gene matrix*/
  for(int i=0; i<dimgene; i++){
    for(int j=0; j<dimgene; j++){
      double covsum=0;
      for(int k=0; k<(dimN + 1); k++){
        
        covsum = covsum + (SignatureCellsTest[i*(dimN+1)+k]-MeansVec[i])*(SignatureCellsTest[j*(dimN+1)+k]-MeansVec[j]);
        
      }
      Covariance_matrix[i*(dimgene)+j] = covsum/(dimN);
    }
  }
  
  
 
  int h = 0;
  int k=0;
  while((h<dimgene) & (k<dimgene)){
    for(int f=h; f<dimgene; f++){
      PivotVec[f]=Covariance_matrix[f*dimgene + k];
    }
    
    int i_result=PivotVec[0]; //larger element pointer initially pointed to first element
    int i_max=0;
    for (int i = 1; i < dimgene; i++)
    {
      if (PivotVec[i] > i_result)
      {
        i_max = i; //updating the pointer to maximum 
        i_result = PivotVec[i];
      }
    }

    for(int s=0; s<dimgene; s++){
      Covariance_matrix[i_max*dimgene + s]=Covariance_matrix[h*dimgene + s];
      Covariance_matrix[h*dimgene + s]=Covariance_matrix[i_max*dimgene + s];
      InvCovariance_matrix[i_max*dimgene + s]=InvCovariance_matrix[h*dimgene + s];
      InvCovariance_matrix[h*dimgene + s]=InvCovariance_matrix[i_max*dimgene + s];
    }
    float f =0;
    float f1 = 0;
    if(Covariance_matrix[h*dimgene + k]>0){
      f = 1/Covariance_matrix[h*dimgene + k];
    }
    if(InvCovariance_matrix[h*dimgene + k]>0){
      f1 = 1/InvCovariance_matrix[h*dimgene + k];
    }
    
    for(int d=k; d<dimgene; d++){
      Covariance_matrix[h*dimgene + d] = Covariance_matrix[h*dimgene + d]*f;
      InvCovariance_matrix[h*dimgene + d] = InvCovariance_matrix[h*dimgene + d]*f1;
    }
    
    for(int i=0; i<dimgene; i++){
      if(i!=h){
        
        if(Covariance_matrix[i*dimgene + k]!=0){
          for(int d=0; d<dimgene; d++){
            Covariance_matrix[i*dimgene + d] = Covariance_matrix[i*dimgene + d] - Covariance_matrix[i*dimgene + d]*Covariance_matrix[h*dimgene + d];
            InvCovariance_matrix[i*dimgene + d] = InvCovariance_matrix[i*dimgene + d] - InvCovariance_matrix[i*dimgene + d]*InvCovariance_matrix[h*dimgene + d];
          }
        }
      }
    }
    
    
    for(int d=0; d<dimgene; d++){
      PivotVec[d]=0.0;
      for(int g=0; g<dimgene; g++){
        Covariance_matrix[d*dimgene + g]=Covariance_matrix[d*dimgene + g];
      }
    }
    h += 1;
    k += 1;
  }
  
  
/*cout << "z: " << z;*/
  
  for(int k=0; k<dimgene; k++){
    columnVec[k]=0;
    for(int i=0; i<dimgene; i++){
      columnVec[k]=columnVec[k]+(SignatureCells[intdimTest*k + z]-MeansVec[i])*InvCovariance_matrix[i*dimgene + k];
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
  PivotVec.clear();
  specificCells.clear();
  SignatureCells.clear();
  PivotVec.clear();
  

  return(distances);
}
