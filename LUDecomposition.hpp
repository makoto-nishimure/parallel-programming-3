#include <vector>
#include <array>
#include "mat.hpp"

using namespace std;

void LUDecompositionBlock(std::vector<std::vector<double> >& matrix, int n){
  int kk, k, i, j, w;

  if(n >= 100){
    w = n/100;
  }else{
    w = 1;
  }
  for(kk = 0; kk < n-w; kk+=w){
    for(k=kk;k < (kk + w);k++){
       for(j=k+1;j<n;j++){
         matrix[k][j] /= matrix[k][k];
       }
       for(i=k+1;i<kk+w;i++){
         for(j=k+1;j<n;j++){
           matrix[i][j] -= matrix[i][k] * matrix[k][j];
         }
       }
       for(i=kk+w;i<n;i++){
         for(j=k+1;j<kk+w;j++){
           matrix[i][j] -= matrix[i][k] * matrix[k][j];
         }
       }
     }
     for(i=kk+w;i<n;i++){
       for(k=kk;k < (kk + w);k++){
         for(j=kk+w;j<n;j++){
           matrix[i][j] -= matrix[i][k] * matrix[k][j];
         }
       }
     }
   }

   for(k=kk; k<n-1;k++)
    {
      for(j=k+1;j<n;j++){
        matrix[k][j] /= matrix[k][k];
      }
      for(i=k+1;i<n;i++)
        for(j=k+1;j<n;j++)
          matrix[i][j] -= matrix[i][k] * matrix[k][j];
    }

}

void LUDecompositionBlock(std::array<std::array<double, SIZE>, SIZE>& matrix, int n){
  int kk, k, i, j, w;

  if(n >= 100){
    w = n/100;
  }else{
    w = 1;
  }
  for(kk = 0; kk < n-w; kk+=w){
    for(k=kk;k < (kk + w);k++){
       for(j=k+1;j<n;j++){
         matrix[k][j] /= matrix[k][k];
       }
       for(i=k+1;i<kk+w;i++){
         for(j=k+1;j<n;j++){
           matrix[i][j] -= matrix[i][k] * matrix[k][j];
         }
       }
       for(i=kk+w;i<n;i++){
         for(j=k+1;j<kk+w;j++){
           matrix[i][j] -= matrix[i][k] * matrix[k][j];
         }
       }
     }
     for(i=kk+w;i<n;i++){
       for(k=kk;k < (kk + w);k++){
         for(j=kk+w;j<n;j++){
           matrix[i][j] -= matrix[i][k] * matrix[k][j];
         }
       }
     }
   }

   for(k=kk; k<n-1;k++)
    {
      for(j=k+1;j<n;j++){
        matrix[k][j] /= matrix[k][k];
      }
      for(i=k+1;i<n;i++)
        for(j=k+1;j<n;j++)
          matrix[i][j] -= matrix[i][k] * matrix[k][j];
    }

}

void LUDecompositionBlock(mat_t matrix, int n){
  int kk, k, i, j, w;

  if(n >= 100){
    w = n/100;
  }else{
    w = 1;
  }
  for(kk = 0; kk < n-w; kk+=w){
    for(k=kk;k < (kk + w);k++){
       for(j=k+1;j<n;j++){
         matrix[k][j] /= matrix[k][k];
       }
       for(i=k+1;i<kk+w;i++){
         for(j=k+1;j<n;j++){
           matrix[i][j] -= matrix[i][k] * matrix[k][j];
         }
       }
       for(i=kk+w;i<n;i++){
         for(j=k+1;j<kk+w;j++){
           matrix[i][j] -= matrix[i][k] * matrix[k][j];
         }
       }
     }
     for(i=kk+w;i<n;i++){
       for(k=kk;k < (kk + w);k++){
         for(j=kk+w;j<n;j++){
           matrix[i][j] -= matrix[i][k] * matrix[k][j];
         }
       }
     }
   }

   for(k=kk; k<n-1;k++)
    {
      for(j=k+1;j<n;j++){
        matrix[k][j] /= matrix[k][k];
      }
      for(i=k+1;i<n;i++)
        for(j=k+1;j<n;j++)
          matrix[i][j] -= matrix[i][k] * matrix[k][j];
    }

}

void naiveLU(vector<vector<double> >& matrix, int n){
  for(int k=0;k<n-1;k++)
   {
     for(int j=k+1;j<n;j++){
       matrix[k][j] /= matrix[k][k];
     }
     for(int i=k+1;i<n;i++)
       for(int j=k+1;j<n;j++)
         matrix[i][j] -= matrix[i][k] * matrix[k][j];
   }
}
