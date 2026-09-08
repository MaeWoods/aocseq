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
std::vector<float> isoForest(int num_trees, int ngenes,
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

// -*- C++ -*-
//===-------------------------- __tree(...) ----------------------------------===//
//Recursive function that is called repeatedly until the algorithm is complete.
//===----------------------------------------------------------------------===//

// [[Rcpp::export]]
std::vector<float> tree(std::vector<float> dataf,std::vector<float> &datainptf,
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
    datID[q] = qv;
    dat[q]=datainptf[genekurtosis*dimN + q];
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
        if(dat[j]<mindata){
          mindata=dat[j];
        }
        if(dat[j]>maxdata){
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

