#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


//calculate the value of ||Xj'Y/n||
double* CrossProduct(double *x, double *y,int begin,int end, int n)
{
  //error checking
  if(end<begin)
  {
    printf("inside CrossProduct: end should not be smaller than begin!\n");
    abort();
  }
  double * val=Calloc(end-begin+1, double);
  
  for(int j=begin;j<=end;j++)
  {
    val[j-begin]=0;
    for(int i=0;i<n;i++)
    {
      val[j-begin]+=x[j*n+i]*y[i];
    }
    val[j-begin]=val[j-begin]/n;
    //printf(" j:%d, val:%f\n",j,val[j]);
  }  
  return val;
}

double CrossProductL2Norm(double *x, double *y,int begin,int end, int n)
{
  double * val=CrossProduct(x,y,begin,end,n);
  double sumSq=0;
  for(int i=0;i<end-begin+1;i++)
  {
    sumSq+=val[i]*val[i];
    //printf("val:%f\n",val[i]);
  }
  sumSq=sqrt(sumSq);
  //printf(" begin:%d, end:%d, sumsq:%f\n",begin,end,sumSq);
  Free(val);
  return sumSq;
}
//find the max ||Xj'Y/n|| for max lamda
SEXP MaxProduct(SEXP x_, SEXP y_, SEXP groupInfo_)
{
  
  //declaration  
  int p=nrows(groupInfo_);  
  int n=nrows(y_);  
  double *x=REAL(x_);  
  double *y=REAL(y_);  
  double *groupInfo=REAL(groupInfo_);  
  int begin=-1,end=-1;
  double maxVal=0;
  double val=0;
  
  for(int j=0;j<p;j++)
  {
    begin=end+1;
    end+=groupInfo[j];
    val= CrossProductL2Norm(x,y,begin,end,n);
    //printf("inside maxproduct: j:%d val:%f\n",j,val);
    if(maxVal<val)
    {
      maxVal=val;
    }     
  }  
  SEXP res;
  PROTECT(res = allocVector(REALSXP, 1));
  REAL(res)[0]=maxVal;
  UNPROTECT(1);
  return res;
  
}