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

//calculate the value of ||Xj'Y\n||
double VectorProduct(double *x, double *y)
{
  int n=sizeof(y);
  double val=0;
  for(int i=0;i<n;i++)
  {
    val+=x[i]*y[i];
  }  
  return val;
}

double UpdateBeta(double z,double lambda,double c)
{
  if(z>lambda)
  {
    return (z-lambda)/c;
  }
  else if(z+lambda<0)
  {
    return(z+lambda)/c;
  }
  else
  {
    return 0;
  }
}

SEXP CleanupG(double *r, double *betaPre, double *wPre, double * shift, 
              double *c, double *lam1, double *lam2,
              SEXP beta_, SEXP w_, SEXP loss_, SEXP wloss_, SEXP iter_) 
{
  Free(r);
  Free(betaPre);
  Free(wPre);
  Free(shift);
  Free(c);
  Free(lam1);
  Free(lam2);
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 5));
  SET_VECTOR_ELT(res, 0, beta_);
  SET_VECTOR_ELT(res, 1, w_);
  SET_VECTOR_ELT(res, 2, wloss_);
  SET_VECTOR_ELT(res, 3, loss_);
  SET_VECTOR_ELT(res, 4, iter_);
  UNPROTECT(6);
  return(res);
}


SEXP INNERREG( SEXP X_, SEXP Y_, SEXP Penalty1_, SEXP Penalty2_, SEXP Lambda1_, SEXP Lambda2_,
               SEXP Beta0_, SEXP W0_, SEXP Delta_, SEXP MaxIter_, 
               SEXP Intercept_, SEXP StarBeta_, SEXP StarW_ )
{
  ////printf("begin\n");
  //data convert
  double *x=REAL(X_);
  double *y=REAL(Y_);
  const char *penalty1 = CHAR(STRING_ELT(Penalty1_, 0));
  const char *penalty2 = CHAR(STRING_ELT(Penalty2_, 0));
  double *lambda1=REAL(Lambda1_);
  double *lambda2=REAL(Lambda2_);
  double *beta0=REAL(Beta0_);
  double *w0=REAL(W0_);
  double *starBeta=NULL;
  if(StarBeta_!=NULL)
  {
    starBeta=REAL(StarBeta_);
  }
  double *starW=NULL;
  if(StarW_!=NULL)
  {
    starW=REAL(StarW_);
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
  SEXP res_, beta_, w_, loss_, wloss_, iter_;
  PROTECT(beta_ = allocVector(REALSXP, L1*L2*m));
  PROTECT(w_ = allocVector(REALSXP, L1*L2*n));
  PROTECT(loss_ = allocVector(REALSXP, L1*L2));
  PROTECT(wloss_ = allocVector(REALSXP, L1*L2));
  PROTECT(iter_ = allocVector(REALSXP, L1*L2)); 
  double * beta=REAL(beta_);
  double * w=REAL(w_);
  double *loss=REAL(loss_);
  double *wloss=REAL(wloss_);
  double * iter=REAL(iter_);
  
  double *betaPre = Calloc(m, double);
  double *wPre=Calloc(n, double);
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
  if(StarW_==NULL)
  {
    for(int i=0;i<n;i++)
   {
     wPre[i]=1;
   }
  }
  else
  {
    for(int i=0;i<n;i++)
   {
     wPre[i]=starW[i];
   }
  }
  
  for(int i=0;i<L1*L2;i++)
  {
    loss[i]=wloss[i]=iter[i]=0;
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
  double *c=Calloc(m, double);
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
    for(int i=0;i<m;i++)
    {
      c[i]=0;
    }
    loss[0*n+l1]=VectorProduct(y,y); //initial loss[l1,1]
    
    if(strcmp(penalty1,"log")==0)
    {
      
      for(int i=0;i<n;i++)
      {
        lam1[i]=sqrt(lambda1[l1]/fabs(log(w0[i]))*n) ;//init sqrt(lambda1/fabs(log(w0))n)
        
        
      }
    }
    else if(strcmp(penalty1,"1-w0")==0)//1-w0
    {
      //fprintf(f,"lam1\n");
      for(int i=0;i<n;i++)
      {
        lam1[i]=lambda1[l1]/fabs(1-w0[i])*n;//init sqrt(lambda1/fabs(log(w0))n)
        //fprintf(f,"lam1=%f lambda1=%f w0=%f ",lam1[i],lambda1[l1],w0[i]);
      }
      
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
        //calculate coefficient c
        for(int j=0;j<m;j++)
        {
          //fprintf(f,"\n cj \n");
          c[j]=0;
          for(int i=0;i<n;i++)
          {
            c[j]+=(x[j*n+i]*wPre[i])*(x[j*n+i]*wPre[i]);
            
          }
          c[j]=c[j]/n;
          //fprintf(f,"%f ",c[j]);
        }
        
         //////printf("enter iteration for each beta\n");
        //iteration for each beta
        for(int j=0;j<m;j++)
        {
          //(1)calculate zj 
          double zj=0;
          for(int i=0;i<n;i++)
          {
            zj+=x[j*n+i]*wPre[i]*wPre[i]*r[i];
          }
          zj=zj/n+c[j]*betaPre[j];
          ////fprintf(f,"\n zj: %f \n",zj);
          //(2)update betaj
          if(strcmp(penalty2,"LASSO")==0)
          {
            beta[j*L1*L2+l2*L1+l1]=UpdateBeta(zj,lam2[j],c[j]);
          }
          else if(strcmp(penalty2,"RIDGE")==0)
          {
            beta[j*L1*L2+l2*L1+l1]=zj/(c[j]+lam2[j]);
          }
          //fprintf(f,"\nbeta(%d,%d,%d)=%f ",l1,l2,j,beta[j*L1*L2+l2*L1+l1]);
          //(3)update r
          shift[j]=beta[j*L1*L2+l2*L1+l1]-betaPre[j];
          for(int i=0;i<n;i++)
          {
            r[i]-=x[j*n+i]*shift[j];
          }
        }
        //fprintf(f,"\n");
        //printf("enter update w\n");
        //update w
        if(strcmp(penalty1,"log")==0)
        {
          double fabsr=0;
          for(int i=0;i<n;i++)
          {
            fabsr=fabs(r[i]);
            if(fabsr>lam1[i])
            {
              w[i*L1*L2+l2*L1+l1]=lam1[i]/fabsr;
            }
            else 
            {
              w[i*L1*L2+l2*L1+l1]=1;
            }
          }
        }       
        else if(strcmp(penalty1,"1-w0")==0)//1-w0
        {
          //fprintf(f,"w:\n" );
          double sqr=0;
          for(int i=0;i<n;i++)
          {
            sqr=r[i]*r[i];
            if(sqr>lam1[i])
            {
              w[i*L1*L2+l2*L1+l1]=lam1[i]/sqr;
            }
            else 
            {
              w[i*L1*L2+l2*L1+l1]=1;
            }
            //fprintf(f,"w(%d,%d,%d)=%f ",l1,l2,i,w[i*L1*L2+l2*L1+l1]);
            //fprintf(f,"\n");
          }
        }
        /*
         else if (strcmp(penalty1,"null")==0)
        {
          //make sure all ws are 1
          printf("it is null\n");
        }*/
        
        
        for(int i=0;i<n;i++)
        {
          shift[m+i]=w[i*L1*L2+l2*L1+l1]-wPre[i];
        }
        
        //update betaPre and wPre for next iteration
        for(int i=0;i<m;i++)
        {
          betaPre[i]=beta[i*L1*L2+l2*L1+l1];
        }
        for(int i=0;i<n;i++)
        {
          wPre[i]=w[i*L1*L2+l2*L1+l1];
        }
        
        //Check for convergence
        if(VectorProduct(shift,shift)<delta)
        {
          break;
        }
        
      }//end for the inner loop
      //printf("exit inner loopw\n");
      //compute square of loss
      loss[l2*L1+l1]=VectorProduct(r,r);
      double temp=0;
      for(int i=0;i<n;i++)
      {
        temp+=r[i]*wPre[i]*r[i]*wPre[i];
      }
      wloss[l2*L1+l1]=temp;
      
      //update for next lambda
      /*
      for(int i=0;i<m;i++)
      {
         betaPre[i]=0;
      }
      for(int i=0;i<n;i++)
      {
         wPre[i]=1;
         r[i]=y[i];
      }*/
      
    }//end iteration for each lambda2 fixed lambda1
    
  }//end iteration for each lambda1
  //fclose(f);
  //clean and return
  res_=CleanupG(r, betaPre, wPre, shift, c, lam1, lam2,
            beta_, w_, loss_, wloss_, iter_);          
  return res_;       
}

