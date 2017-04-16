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