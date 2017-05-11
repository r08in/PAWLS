#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#define false 0
#define true 1
typedef int bool; 

double VectorProduct2(double *x, double *y)
{
  int n=sizeof(y);
  double val=0;
  for(int i=0;i<n;i++)
  {
    val+=x[i]*y[i];
  }  
  return val;
}

double UpdateSoftThreshold(double z,double lambda)
{
  if(z>lambda)
  {
    return (z-lambda);
  }
  else if(z+lambda<0)
  {
    return(z+lambda); 
  }
  else
  {
    return 0;
  }
}

SEXP CleanupG2(double *r, double *betaPre, double *gamPre, double * shift,  double *lam1, double *lam2,
              SEXP beta_, SEXP Gam_, SEXP loss_, SEXP iter_) 
{
  Free(r);
  Free(betaPre);
  Free(gamPre);
  Free(shift);
  Free(lam1);
  Free(lam2);
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 4));
  SET_VECTOR_ELT(res, 0, beta_);
  SET_VECTOR_ELT(res, 1, Gam_);
  SET_VECTOR_ELT(res, 2, loss_);
  SET_VECTOR_ELT(res, 3, iter_);
  UNPROTECT(5);
  return(res);
}


SEXP INNERPAMLS( SEXP X_, SEXP Y_, SEXP Lambda1_, SEXP Lambda2_,
               SEXP Beta0_, SEXP Gam0_, SEXP Delta_, SEXP MaxIter_, 
               SEXP Intercept_, SEXP StarBeta_, SEXP StarGam_ )
{
  ////printf("begin\n");
  //data convert
  double *x=REAL(X_);
  double *y=REAL(Y_);
  double *lambda1=REAL(Lambda1_);
  double *lambda2=REAL(Lambda2_);
  double *beta0=REAL(Beta0_);
  double *Gam0=REAL(Gam0_);
  double *starBeta=NULL;
  if(StarBeta_!=NULL)
  {
    starBeta=REAL(StarBeta_);
  }
  double *starGam=NULL;
  if(StarGam_!=NULL)
  {
    starGam=REAL(StarGam_);
  }
  double delta=REAL(Delta_)[0];
  double maxIter = REAL(MaxIter_)[0];
  int intercept=REAL(Intercept_)[0];
  
  ////printf("pass data convert\n");
  
  //data declare
  int n=nrows(X_);
  int m=ncols(X_);
  int L1=length(Lambda1_);
  int L2=length(Lambda2_);
  int lstart1=0, lstart2=0;
  
  //data return
  SEXP res_, beta_, Gam_, loss_, iter_;
  PROTECT(beta_ = allocVector(REALSXP, L1*L2*m));
  PROTECT(Gam_ = allocVector(REALSXP, L1*L2*n));
  PROTECT(loss_ = allocVector(REALSXP, L1*L2));
  PROTECT(iter_ = allocVector(REALSXP, L1*L2)); 
  double * beta=REAL(beta_);
  double * Gam=REAL(Gam_);
  double *loss=REAL(loss_);
  double * iter=REAL(iter_);
  
  double *betaPre = Calloc(m, double);
  double *gamPre=Calloc(n, double);
  double *r=Calloc(n, double);
  
  //initial
  if(StarBeta_==NULL)
  {
     for(int i=0;i<m;i++)
   {
     betaPre[i]=0;
   }
  }
  else
  {
    for(int i=0;i<m;i++)
   {
     betaPre[i]=starBeta[i];
   }
  }
  if(StarGam_==NULL)
  {
    for(int i=0;i<n;i++)
   {
     gamPre[i]=0;
   }
  }
  else
  {
    for(int i=0;i<n;i++)
   {
     gamPre[i]=starGam[i];
   }
  }
  
  for(int i=0;i<L1*L2;i++)
  {
    loss[i]=iter[i]=0;
  }
  
  if(StarBeta_==NULL)
  {
    for(int i=0;i<n;i++)
    {
      r[i]=y[i];
    }
  }
  else
  {
    for(int i=0;i<n;i++)
   {
     double temp=0;
     for(int j=0;j<m;j++)
     {
       temp+=x[j*n+i]*betaPre[j];
     }
     r[i]=y[i]-temp;
   }     
  }
  
  
  //temp
  double *shift=Calloc(n+m, double);
  double *lam1=Calloc(n, double);
  double *lam2=Calloc(m, double);
  
  //test
 // FILE *f = fopen("debug.txt", "a");
  //fprintf(f,"begin:\n");
  ////printf("enter interation for each lambda1\n");
  //interation for each lambda1
  for(int l1=lstart1;l1<L1;l1++)
  {
    //initial
    for(int i=0;i<n+m;i++)
    {
      shift[i]=0;
    }
    
    for(int i=0;i<n;i++)
    {
      
      lam1[i]=lambda1[l1]/fabs(Gam0[i])*n;
    }
    
     //printf("enter interation for each lambda2\n");
    //iteration for each lambda2
    for(int l2=lstart2;l2<L2;l2++)
    {
      //fprintf(f,"\nlam2\n");
      for(int i=0;i<m;i++)
      {
        lam2[i]=lambda2[l2]/fabs(beta0[i]);
        //fprintf(f,"%f ",lam2[i]);
      }
      
      if(intercept==true)
      {
        lam2[0]=0;
        //lam2[1]=0;//for test
      }
       //printf("enter iteration for all covariates\n");
      //iteration for all covariates
      while(iter[l2*L1+l1]<maxIter)
      {
        iter[l2*L1+l1]+=1;
        
         //////printf("enter iteration for each beta\n");
        //iteration for each beta
        for(int j=0;j<m;j++)
        {
          //(1)calculate zj 
          double zj=0;
          for(int i=0;i<n;i++)
          {
            zj+=x[j*n+i]* (r[i] - gamPre[i]);
          }
          zj=zj/n + betaPre[j];
          ////fprintf(f,"\n zj: %f \n",zj);
          //(2)update betaj
          beta[j*L1*L2+l2*L1+l1]=UpdateSoftThreshold(zj,lam2[j]);
          //fprintf(f,"\nbeta(%d,%d,%d)=%f ",l1,l2,j,beta[j*L1*L2+l2*L1+l1]);
          //(3)update r
          shift[j]=beta[j*L1*L2+l2*L1+l1]-betaPre[j];
          //fprintf(f,"\n shift %d: %f \n",j,shift[j]);
          for(int i=0;i<n;i++)
          {
            r[i]-=x[j*n+i]*shift[j];
          }
        }
        //fprintf(f,"\n");
        //printf("enter update Gam\n");
        //update Gam
        for(int i=0;i<n;i++)
        {
          Gam[i*L1*L2+l2*L1+l1]=UpdateSoftThreshold(r[i],lam1[i]);
        }
  
        for(int i=0;i<n;i++)
        {
          shift[m+i]=Gam[i*L1*L2+l2*L1+l1]-gamPre[i];
          //fprintf(f,"\n shift %d: %f \n",i+m,shift[m+i]);
        }
        
        //update betaPre and gamPre for next iteration
        for(int i=0;i<m;i++)
        {
          betaPre[i]=beta[i*L1*L2+l2*L1+l1];
        }
        for(int i=0;i<n;i++)
        {
          gamPre[i]=Gam[i*L1*L2+l2*L1+l1];
        }
        
        //Check for convergence
        if(VectorProduct2(shift,shift)<delta)
        {
          break;
        }
        //fprintf(f,"\n VectorProduct2= %f ",VectorProduct2(shift,shift));
        
      }//end for the inner loop
      //printf("exit inner loopw\n");
      //compute square of loss
      double temp=0;
      for(int i=0;i<n;i++)
      {
        temp+=(r[i]-gamPre[i]) * (r[i]-gamPre[i])  ;
      }
      loss[l2*L1+l1]=temp;
      
      
    }//end iteration for each lambda2 fixed lambda1
    
  }//end iteration for each lambda1
  //fclose(f);
  //clean and return
  res_=CleanupG2(r, betaPre, gamPre, shift,lam1, lam2,
            beta_, Gam_, loss_,iter_);          
  return res_;       
}

