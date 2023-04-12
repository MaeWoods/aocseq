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
int GetOptimalMarkers(int UniqueGenes, 
     int Nmarkersin, int Nsetsin, std::vector<int> &Genedatain) {

  vector<int> Genedata = Genedatain;
  
  std::vector<vector<vector<int> > > sets;
  std::vector<vector<vector<int> > > markergenes;
  
  std::vector<vector<vector<int> >> markergenesins;
  std::vector<vector<vector<int> >> setsins;
  markergenesins.resize(1);
  markergenesins[0].resize(1);
  markergenesins[0][0].resize(1);
  setsins.resize(1);
  setsins[0].resize(1);
  setsins[0][0].resize(1);
  
  int counter_h=0;
  int counterh2=0;
  int matchedsets=0;
  int m=0;
  int q=0;
  int g=0;
  int s=0;
  int b=0;
  int matchB=0;
  int counterh3=0;

    for(s=0; s<(Nsetsin*Nmarkersin); s++){
      for(g=0; g<UniqueGenes; g++){
        Genedata[g*(Nsetsin*Nmarkersin) + s] = Genedatain[g*(Nsetsin*Nmarkersin) + s];
      }
    }
    
        
        sets.resize(1);
        sets[0].resize(1);
        sets[0][0].resize(1);
        markergenes.resize(1);
        markergenes[0].resize(1);
        markergenes[0][0].resize(1);
        sets[0][0][0]=0;
        markergenes[0][0][0]=0;
        for(g=0; g<UniqueGenes; g++){
        
        for(s=0; s<(Nsetsin*Nmarkersin); s++){
          counterh2=0;
          if(g==0){
          
          if(Genedata[g*(Nsetsin*Nmarkersin) + s]>0){
            
            if(counter_h==0){
              sets[0][0][0]=s;
              markergenes[0][0][0]=g;
              counter_h += 1;
            }
            else{
              sets[0][0].push_back(s);
              markergenes[0][0].push_back(g);
            }
            
          }
          }
          else{
            
              if(Genedata[g*(Nsetsin*Nmarkersin) + s]>0){
                
                for(q=0; q<sets[0].size(); q++){
                  
                  matchedsets=0;
                  for(m=0; m<sets[0][q].size(); m++){
                    if(Genedata[g*(Nsetsin*Nmarkersin) + sets[0][q][m]]>0){
                      matchedsets += 1;
                    }
                    
                  }
                  
                  if(matchedsets==(int)sets[0][q].size()){
                    
                    markergenesins[0][0][0] = g;
                    
                    markergenes[0][q].reserve(markergenes[0][q].size() + markergenesins[0][0].size());
                    markergenes[0][q].insert(markergenes[0][q].end(), markergenesins[0][0].begin(), markergenesins[0][0].end());
                    
                  }
                  else{
                    if(s==((Nsetsin*Nmarkersin)-1)){
                      matchB=0;
                      int checkb=0;
                      int counterbhprev=0;
                      counterbhprev = floor(sets[0][q][0]/Nsetsin);
                    for(m=0; m<sets[0][q].size(); m++){
                      if(counterbhprev!=floor(sets[0][q][m]/Nsetsin)){
                        checkb =0;
                      }
                 
                     
                      counterbhprev = floor(sets[0][q][m]/Nsetsin);
                    int step=0;
                    
                      for(b=0;b<Nsetsin;b++){
                      step=Nmarkersin+b*Nsetsin;
                      if(((sets[0][q][m]>=(step-Nsetsin))&&(checkb==0)&&(sets[0][q][m]<step))&&(Genedata[g*(Nsetsin*Nmarkersin) + sets[0][q][m]]>0)){
                        matchB += 1;
                        checkb = 1;
                      }
                      }
                      }
                    
                    //cout << "matchB: " << matchB << " Nsetsin: " << Nsetsin << " gamma: " << gamma << "counterh3: " << counterh3 <<   std::endl;

                    if(matchB==Nsetsin){
                      counterh3=0;
                      for(m=0; m<sets[0][q].size(); m++){

                          if(Genedata[g*(Nsetsin*Nmarkersin) + sets[0][q][m]]>0){
                          if(counterh3==0){
                            setsins[0][0][0] = sets[0][q][m];
                            markergenesins[0][0][0] = g;
                            
                            markergenes[0].reserve(markergenes[0].size() + markergenesins[0].size());
                            markergenes[0].insert(markergenes[0].end(), markergenesins[0].begin(), markergenesins[0].end());
                            sets[0].reserve(sets[0].size() + setsins[0].size());
                            sets[0].insert(sets[0].end(), setsins[0].begin(), setsins[0].end());
                            counterh3 += 1;
                          }
                          else{
                            sets[0][(sets[0].size())-1].push_back(sets[0][q][m]);
                            for(int n=0; n<markergenes[0][q].size(); n++){
                            markergenes[0][(markergenes[0].size())-1].push_back(markergenes[0][q][n]);
                            }
                          }
                          
                          }
                        }
                      
                    }
              
                    }
 
                  }
                  
                  
                }
                
                 if(counterh2==0){
                
                  setsins[0][0][0] = s;
                  markergenesins[0][0][0] = g;
                  
                  markergenes[0].reserve(markergenes[0].size() + markergenesins[0].size());
                  markergenes[0].insert(markergenes[0].end(), markergenesins[0].begin(), markergenesins[0].end());
                  
                  sets[0].reserve(sets[0].size() + setsins[0].size());
                  sets[0].insert(sets[0].end(), setsins[0].begin(), setsins[0].end());
                  counterh2 += 1;
                  }
                  else{
                  sets[0][(sets[0].size())-1].push_back(s);
                  markergenes[0][(markergenes[0].size())-1].push_back(g);
                 
                  }
                
              }
              
          }
          }
          
        }
        
        
        for(q=0; q<sets[0].size(); q++){
          //std::cout << "set: " << q << ": ";
          for(m=0; m<sets[0][q].size(); m++){
            
           // std::cout << sets[0][q][m] << ", ";

          }
          //std::cout << " " << std::endl;
        }
        
        
        for(q=0; q<markergenes[0].size(); q++){
          //std::cout << "markergene: " << q << ": ";
          for(m=0; m<markergenes[0][q].size(); m++){
            
          //std::cout << markergenes[0][q][m] << ", ";
            
          }
          //std::cout << " " << std::endl;
        }
        
        std::vector<int> maxgeneintersection;
        maxgeneintersection.resize(markergenes[0].size());
        
        for(q=0; q<markergenes[0].size(); q++){
        
        sort(markergenes[0][q].begin(), markergenes[0][q].end());
        vector<int>::iterator it;
        it = unique(markergenes[0][q].begin(), markergenes[0][q].end());  
        
        markergenes[0][q].resize(distance(markergenes[0][q].begin(),it)); 
        
        maxgeneintersection[q] = markergenes[0][q].size();
        //std::cout << "Size of q: " << maxgeneintersection[q] << std::endl;
        //std::cout << "markergenes unique: " << q << ": ";
        for(m=0; m<markergenes[0][q].size(); m++){
          
          //std::cout << markergenes[0][q][m] << ", ";
          
        }
        //std::cout << " " << std::endl;
        
        }
        
        int maxgene = std::max_element(maxgeneintersection.begin(),maxgeneintersection.end()) - maxgeneintersection.begin();
        std::cout << "Maximal intersection completed, maxgene: " << maxgene << std::endl;
        
         // std::cout << "set maxgene: " << maxgene << ": ";
          for(m=0; m<sets[0][maxgene].size(); m++){
            
           // std::cout << sets[0][maxgene][m] << ", ";
            
          }
          //std::cout << " " << std::endl;
        
          
        
  markergenesins.clear();
  setsins.clear();
  Genedata.clear();
  sets.clear();
  markergenes.clear();
  
/*delete &counter_h;
delete &counterh2;
delete &matchedsets;
delete &m;
delete &q;
delete &g;
delete &s;*/

  return(maxgene);
}