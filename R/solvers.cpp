// -*- C++ -*-
//===-------------------------- __solvers.cpp ----------------------------------===//
//
// Part of the aocseq Project, under the BSD2 license.
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
  
  for(int n=0; n<(dimN+1); n++){
    for(int m=0; m<(dimgene); m++){
      PivotVec[m]=0;
      MeansVec[m] += SignatureCellsTest[m*(dimN+1)+n]/(dimN+1);
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
      for(int k=0; k<(dimN + 1); k++){
        
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
    for(int f=h; f<dimgene; f++){
      PivotVec[f]=Covariance_matrix[f*dimgene + k];
    }
    
    float i_result=PivotVec[h]; //larger element pointer initially pointed to first element
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
    }
    
    for(int d=0; d<dimgene; d++){
      if((Covariance_matrix[h*dimgene + k]!=0)&(Covariance_matrix[h*dimgene + d]!=0)&(InvCovariance_matrix[h*dimgene + d]!=0)){
      if(d==k){
        Covariance_matrix[h*dimgene + d] = 1;
          InvCovariance_matrix[h*dimgene + d] = InvCovariance_matrix[h*dimgene + d]/Covariance_matrix[h*dimgene + k];
      }
      else{
      Covariance_matrix[h*dimgene + d] = Covariance_matrix[h*dimgene + d]/Covariance_matrix[h*dimgene + k];
      InvCovariance_matrix[h*dimgene + d] = InvCovariance_matrix[h*dimgene + d]/Covariance_matrix[h*dimgene + k];
      }
      }
     
    }
    
    for(int i=(h+1); i<dimgene; i++){
          for(int d=0; d<dimgene; d++){
            if(Covariance_matrix[h*dimgene + d]!=0){
            Covariance_matrix[i*dimgene + d] = Covariance_matrix[i*dimgene + d] - Covariance_matrix[i*dimgene + d]*Covariance_matrix[h*dimgene + d];
            InvCovariance_matrix[i*dimgene + d] = InvCovariance_matrix[i*dimgene + d] - Covariance_matrix[i*dimgene + d]*Covariance_matrix[h*dimgene + d];
            }
          }
  
    }
    
    
    for(int d=0; d<dimgene; d++){
      PivotVec[d]=0.0;
      for(int h1=0; h1<dimgene; h1++){
          Covariance_matrix_cp[h1*dimgene + d]=Covariance_matrix[h1*dimgene + d];
          InvCovariance_matrix_cp[h1*dimgene + d]=InvCovariance_matrix[h1*dimgene + d];
        
      }
    }
    h += 1;
    k += 1;
  }
  //Next do the reverse to create a diagonal matrix
  int p = 0;
  while((p<dimgene)){
    int h1=dimgene-(p+1);
    if(h1>0){
    for(int m=1; m<(h1+1); m++){
      for(int s=0; s<dimgene; s++){
        int i=h1-m;
        int d=dimgene-(s+1);
        Covariance_matrix[i*dimgene + d] = Covariance_matrix[i*dimgene + d] - Covariance_matrix[i*dimgene + d]*Covariance_matrix[h1*dimgene + d];
        InvCovariance_matrix[i*dimgene + d] = InvCovariance_matrix[i*dimgene + d] - Covariance_matrix[i*dimgene + d]*Covariance_matrix[h1*dimgene + d];
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

// -*- C++ -*-
//===-------------------------- __isoForest(...) ----------------------------------===//
//data is an array of the data for a single gene (gene kurtosis) - this is the gene with maximum kurtosis across all cells
//dataid keeps track of the data points being updated and split
//datainpt is the original data
//currheight is the current height of the tree
//currheight_vec is the height of all the data points
//dataidright keeps track of the data points being updated after they have been split
//dataidleft keeps track of the data points being updated after they have been split
//===----------------------------------------------------------------------===//

// [[Rcpp::export]]
vector<float> tree(std::vector<float> dataf,std::vector<float> &datainptf, 
                   std::vector<float> dataidf, std::vector<float> &currentheightvecf,
                   std::vector<float> &kurtosisvecf,
                   int dimN, int dimSplit, int currheight, 
                   int maxheight, int genekurtosis, int ngenes,
                   std::vector<float> dataidrightf,std::vector<float> dataidleftf,
                   std::vector<float> vecleftf,std::vector<float> vecrightf){
  
  std::vector<float> dat(dimN);
  std::vector<float> datID(dimN);
  std::vector<float> vL(dimN);
  std::vector<float> vR(dimN);
  std::vector<float> v_IL(dimN);
  std::vector<float> v_IR(dimN);
  for(int q=0; q<dimN; q++){
    dat[q]=dataf[q];
    datID[q]=dataidf[q];
    vL[q]=vecleftf[q];
    vR[q]=vecrightf[q];
    v_IL[q]=dataidleftf[q];
    v_IR[q]=dataidrightf[q];
  }
  
  if(currheight==0){
    float qv = 0.0;
  for(int q=0; q<dimN; q++){
    currentheightvecf[q]=(-1);
    dataidf[q] = qv;
    dataf[q]=datainptf[genekurtosis*dimN + q];
    qv += 1.0;
  }
  }
  else{
    ///find the gene with maximum kurtosis 
    for(int g=0; g<ngenes; g++){
      float fourth_moment_fill=0;
      float fourth_moment=0;
      float mean=0;
      float standard_dev_fill=0;
      float standard_dev=0;
      
      for(int q=0; q<dimN; q++){
        if(datID[q]!=(-1)){
        mean += (datainptf[g*(dimN)+q])/dimSplit;
        }
      }
      for(int q=0; q<dimN; q++){
        if(datID[q]!=(-1)){
        fourth_moment_fill += (datainptf[g*(dimN)+q]-mean)/dimSplit;
        standard_dev_fill += ((datainptf[g*(dimN)+q]-mean)*(datainptf[g*(dimN)+q]-mean))/dimSplit;
        }
      }
      
      fourth_moment=fourth_moment_fill*fourth_moment_fill*fourth_moment_fill*fourth_moment_fill;
      standard_dev = standard_dev_fill*standard_dev_fill;
      if(standard_dev==0){
        kurtosisvecf[g]=0; 
      }
      else{
      kurtosisvecf[g]=(fourth_moment/standard_dev);
      }
    }
    float maxkurtosisf =0;
    for(int j=0; j<ngenes; j++){
      if(j==0){
        maxkurtosisf=kurtosisvecf[j];
      }
      if(kurtosisvecf[j]>maxkurtosisf){
        maxkurtosisf=kurtosisvecf[j];
      }
    }
    
    for(int j=0; j<ngenes; j++){
      if(kurtosisvecf[j]==maxkurtosisf){
        genekurtosis=j;
      }
    }

      for(int q=0; q<dimN; q++){
      if(datID[q]!=(-1)){
      dat[q]=datainptf[genekurtosis*dimN + q];
      }
    }
    
  }

  int outout=0;
  for(int q=0; q<dimN; q++){
    if(currentheightvecf[q]==(-1)){
      outout=1;
    }
  }
  if(outout==0){
    std::cout << "Error: all cells have been assigned a distance" << std::endl;
  }
  
  if((outout==0)|(currheight == (maxheight))|(dimSplit==1)){
      for(int g=0; g<dimN; g++){
        if(datID[g]!=(-1)){
        currentheightvecf[g]=currheight;
        }
      }
      
      for(int g=0; g<dimN; g++){
        dat[g]=datainptf[genekurtosis*dimN + g];
      }
    }
    else{
      float mindata=0,maxdata=0;
      int counter_j=0;
      for(int j=0; j<dimN; j++){
        if(datID[j]!=(-1)){
        if(counter_j==0){
          mindata=dat[j];
          maxdata=dat[j];
        }
        if(dataf[j]<mindata){
          mindata=dat[j];
        }
        if(dataf[j]>maxdata){
          maxdata=dat[j];
        }
        counter_j += 1;
        }
      }
      
      float split_value = mindata+R::runif(0,1)*maxdata;
      int counter_left=0,counter_right=0;
      for(int j=0; j<dimN; j++){
        vL[j]= dat[j];
        vR[j]= dat[j];
        v_IL[j]= -1;
        v_IR[j]= -1;
      }
      if(mindata==maxdata){
        int split_valueL=ceil(dimSplit/2);
        for(int j=0; j<dimN; j++){
          if(datID[j]!=(-1)){
            if(counter_left<=(split_valueL-1)){
              vL[j] = dat[j];
              v_IL[j]=j;
              counter_left += 1;
            }
            else{
              vR[j] = dat[j];
              v_IR[j]=j;
              counter_right += 1; 
            }
          }
        }
      }
      else{
      float qv=0;
      for(int j=0; j<dimN; j++){
        if(datID[j]!=(-1)){
        if(dat[j]<=split_value){
          
          vL[j] = dat[j];
          v_IL[j]=qv;
          counter_left += 1;
          
          }
        
        else if(dat[j]>split_value){
          vR[j] = dat[j];
          v_IR[j]=qv;
         counter_right += 1; 
          
        }
        }
        qv += 1.0;
      }
      }
      
      if(counter_left!=0){
        for(int j=0; j<dimN; j++){
          dat[j]=vL[j];
          datID[j]=v_IL[j];
        }
        vecleftf = tree(dat,datainptf,datID,currentheightvecf,kurtosisvecf,dimN,counter_left,currheight+1,maxheight,genekurtosis,ngenes,v_IR,v_IL,vL,vR);
      }
      if(counter_right!=0){
        for(int j=0; j<dimN; j++){
          dat[j]=vR[j];
          datID[j]=v_IR[j];
        }
        vecrightf = tree(dat,datainptf,datID,currentheightvecf,kurtosisvecf,dimN,counter_right,currheight+1,maxheight,genekurtosis,ngenes,v_IR,v_IL,vL,vR);
      }
  
  }
   dat.clear(); 
    datID.clear(); 
    vL.clear(); 
    vR.clear(); 
    v_IL.clear(); 
    v_IR.clear(); 
    
  return(currentheightvecf);
}

// [[Rcpp::export]]
vector<float> isoForest(int num_trees, int ngenes, 
                             int dimN, int dimClone, std::vector<float> &distinpt, 
                             std::vector<float> &smallvecinpt,std::vector<float> &smallvecinptcp,
                             std::vector<float> &dist, std::vector<float> &heightinpt,
                             std::vector<float> &currheightinpt,std::vector<float> &kurtosisinpt,
                             std::vector<float> &dataids_Rinpt,std::vector<float> &dataids_Linpt,
                             std::vector<float> &Vec_Rinpt,std::vector<float> &Vec_Linpt,float maxkurtosis,
                             std::vector<float> &Cell_inpt, std::vector<float> &Height_output
){
  vector<float> distances = dist;
  vector<float> testdf = distinpt;
  vector<float> avg_height = heightinpt;
  vector<float> currentheightvec = currheightinpt;
  vector<float> kurtosisvec = kurtosisinpt;
  vector<float> dataids_R = dataids_Rinpt;
  vector<float> dataids_L = dataids_Linpt;
  vector<float> small_vec = smallvecinpt;
  vector<float> small_veccp = smallvecinptcp;
  vector<float> Vec_R = Vec_Rinpt;
  vector<float> Vec_L = Vec_Linpt;
  vector<float> C_I = Cell_inpt;
  vector<float> H_O = Height_output;

  int maxheight=30;
  int genekurtosis=-1;
  
  for(int r=0; r<dimClone; r++){
    distances = dist;
    testdf = distinpt;
    avg_height = heightinpt;
    currentheightvec = currheightinpt;
    kurtosisvec = kurtosisinpt;
    dataids_R = dataids_Rinpt;
    dataids_L = dataids_Linpt;
    small_vec = smallvecinpt;
    small_veccp = smallvecinptcp;
    Vec_R = Vec_Rinpt;
    Vec_L = Vec_Linpt;
    C_I = Cell_inpt;
    H_O = Height_output;
    
    if(r>0){
      for(int g=0; g<ngenes; g++){
      for(int q=0; q<dimN; q++){
        
        if(q==(dimN-1)){
        testdf[g*(dimN)+q]=C_I[g*(dimClone)+r];
        }
        else{
        testdf[g*(dimN)+q]=distinpt[g*(dimN)+q];
        }
        
      }
    }
  }
 
      /*///////////////////////////////////////////////////////*/
      /*End of creating the matrix of unique values ///////////*/
      /*Now implement the kurtosis of all rows across          */
      /*all genes for all cells                                */
      /*///////////////////////////////////////////////////////*/
      /*///////////////////////////////////////////////////////*/
      
      for(int g=0; g<ngenes; g++){
        float fourth_moment_fill=0;
        float fourth_moment=0;
        float mean=0;
        float standard_dev_fill=0;
        float standard_dev=0;
        
        for(int q=0; q<dimN; q++){
          mean += (testdf[g*(dimN)+q])/dimN;
        }
        for(int q=0; q<dimN; q++){
          fourth_moment_fill += (testdf[g*(dimN)+q]-mean)/dimN;
          standard_dev_fill += ((testdf[g*(dimN)+q]-mean)*(testdf[g*(dimN)+q]-mean))/dimN;
        }
        
        fourth_moment=fourth_moment_fill*fourth_moment_fill*fourth_moment_fill*fourth_moment_fill;
        standard_dev = standard_dev_fill*standard_dev_fill;
        kurtosisvec[g]=(fourth_moment/standard_dev);
      }
      
      for(int j=0; j<ngenes; j++){
          if(j==0){
            maxkurtosis=kurtosisvec[j];
          }
          if(kurtosisvec[j]>maxkurtosis){
            maxkurtosis=kurtosisvec[j];
          }
      }
      
      for(int j=0; j<ngenes; j++){
        if(kurtosisvec[j]==maxkurtosis){
          genekurtosis=j;
        }
      }

      /*///////////////////////////////////////////////////////*/
      /*// Use gene with maximum kurtosis for isolation tree //*/
      /*///////////////////////////////////////////////////////*/
      
      
      for(int q=0; q<dimN; q++){
      currentheightvec[q]=-1;
      }
        
      for(int treecounter=0; treecounter<num_trees; treecounter++){
        int current_height_in=0;
        currentheightvec=tree(small_vec,testdf,small_veccp,currentheightvec,kurtosisvec,dimN,dimN,current_height_in,maxheight,genekurtosis,ngenes,dataids_R, dataids_L,Vec_R,Vec_L);

      for(int q=0; q<dimN; q++){
      avg_height[q] += currentheightvec[q]/num_trees;
      }
        
      }
      
      float c=2*(log(dimN-1)+0.5772156649) - (2.0*(log(dimN-1)/(log(dimN))));
      H_O[r]=pow(2,((-1*avg_height[(dimN-1)])/c));

    for(int q=0; q<dimN; q++){
    if(avg_height[q]==(-1)){
      cout << "Some data points were never sampled. Increase num_trees or subsample_count." << endl;
      break;
    }
    }
  
    distances.clear();
    testdf.clear();
    currentheightvec.clear();
    kurtosisvec.clear();
    dataids_R.clear();
    dataids_L.clear();
    small_vec.clear();
    small_veccp.clear();
    Vec_L.clear();
    Vec_R.clear();
    avg_height.clear();
    C_I.clear();
    
}

  return(H_O);
}
