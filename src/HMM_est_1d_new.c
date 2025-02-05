#include <stdio.h> 
#include <stdlib.h>   
#include <math.h>   

#define cte 0.3989423

   double maxi(double a, double b)
   {
      if (a > b)
         return a;
      else
         return b;
   }

double sumx( double *x, int n)
    {
      int i;
      double sum = 0.0;
   
      for(i=0;i<n;i++)
         sum += x[i];
   
      return sum;
   
   }

   double meanx( double *x, int n)
   
   {
       
      return sumx(x,n)/((double)n);
   
   }


   double meanx2( double *x, double mx, int n)
     {
      int i;
      double z, sum = 0.0;
   
      for(i=0;i<n;i++)
      {
         z = x[i]-mx;
         sum += z*z;
      }
      return sum/((double)n-1.0);
   
   }


double Phi(double x)    /*Distribution function of the standard Gaussian r.v.      */
   {
      /* A&S formula 7.1.26 */

      double p  =  0.3275911;
      double a1 =  0.254829592;
      double a2 = -0.284496736;
      double a3 =  1.421413741;
      double a4 = -1.453152027;
      double a5 =  1.061405429;



      int sign = 1;
      if (x < 0)
         sign = -1;
      x = fabs(x)/sqrt(2.0);


      double t = 1.0/(1.0 + p*x);
      double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

      return 0.5*(1.0 + sign*y);
   }


double Sn_1d( double *U, double N) /** U is a sequence of pseudo-observations of uniform */
   {
      double sum1, sum2, sum3;
      int i,j;
   
   
   
      sum1 = 0.0; 
      for(i=0;i<N;i++)
      {
         for(j=0;j<i;j++)
            sum1 += 1.0 - maxi(U[i],U[j]);
      
      }
   
      sum1 = 2.0*sum1/((double)N);
   
      sum2 =0.0;
   
      for(i=0;i<N;i++)
         sum2 += 1.0 - U[i];
   
   
   
   
      sum2 = sum2/((double) N);
   
      sum3 =0.0;
   
      for(i=0;i<N;i++)
         sum3 += 1.0 - U[i]*U[i];
   
   
   
      return sum1 + sum2 - sum3 + ((double)N)/3.0;
   
   }  



   void est_Q(int *reg, int nstates, int N, double **hatQ)
   {
      int i, k, l;
      double sum2;
   
      for(k=0;k<nstates;k++)
      {
         sum2 = 0.0;
      
         for(l=0;l<nstates;l++)
            hatQ[k][l] = 0.0;
      
         for(i=0;i<N-1;i++)
         { 
            if( reg[i] == k)
            { 
               sum2 += 1.0;
               for(l=0;l<nstates;l++)
               {
                  if(reg[i+1] == l)
                     hatQ[k][l] += 1.0;
               }
            }
         }
      
      
         for(l=0;l<nstates;l++)
            hatQ[k][l] = hatQ[k][l]/sum2;
      }
   
   }






  
   void hmm_init_1d(double *R, int n, int m, double *mu, double *sigma, double **Q)
      {
      int i, j, k, n0;
      double *x, mx;
   
      x = (double *)calloc(n,sizeof(double));
   
   
      n0 = floor(1.0*n/(1.0*m));
   
   
   
      for(k=0;k<m;k++)
      {
         for(i=0; i< n0; i++)
            {
                x[i] = R[k*n0+i];
            }
         mx = meanx(x,n0);
         mu[k] = mx;
         sigma[k] = sqrt(meanx2(x,mx,n0));
       
      }
    
      free(x);
   
      for(i=0;i<m;i++)
      {
         for(j=0;j<m;j++)
            Q[i][j] = 1.0/((double) m);
      
      
      }
   }

   void weights(double **Q, int nstates, double *eta, double *W)
   {
      int i,j;
      double sum;
   
   
      for(i=0;i<nstates;i++)
      {
         sum = 0.0;
         for(j=0;j<nstates;j++) 
            sum += eta[j]*Q[j][i];
      
         W[i] =   sum;
      
      }
   
   }

 

   void E_step_1d(double *R, double *mu, double *sigma, double **Q, int N, int m, double **eta, double **g, double **qbar, double **lambda, double ***Lambda, double *L)
      {
   
      int i, j, k;
      double sum, *W, z;
   
      W = (double *)calloc(m, sizeof(double));
     
   
     
      for(i=0;i<N;i++)
      {
         for(k=0;k<m;k++)
         {     
            z = (R[i] - mu[k])/sigma[k];
            g[k][i] = maxi(exp(-0.5*z*z)/sigma[k],1E-10);
         
          
         }
      }
   
  
      for(k=0;k<m;k++)
      {
         eta[0][k] = 1.0/((double) m);
         qbar[N][k] = 1.0/((double) m);
      }
   
   
   
   
      for(i=1;i<=N;i++)
      {
         weights(Q,m,eta[i-1],W);
         sum = 0.0;
      
         for(j=0;j<m;j++)
            sum += W[j]*g[j][i-1];
        L[i-1] = sum;
      
         for(j=0;j<m;j++)
            eta[i][j] = W[j]*g[j][i-1]/sum;
      
         for(k=0;k<m;k++)
         { 
            sum = 0.0;
         
            for(j=0;j<m;j++)
               sum += qbar[N+1-i][j]*Q[k][j]*g[j][N-i];
         
            qbar[N-i][k] = sum;
         }
      
         sum = 0.0;
         for(k=0;k<m;k++) 
            sum += qbar[N-i][k];
      
         for(k=0;k<m;k++) 
            qbar[N-i][k] = qbar[N-i][k]/sum;
      
      }

         
       /**************  lambda(t,j) ***************/
      for(i=1;i<=N;i++)
      {
         sum = 0.0;
         for(k=0;k<m;k++)
         {
            lambda[i][k] = eta[i][k]*qbar[i][k];
            sum += lambda[i][k];
         }
         for(k=0;k<m;k++)
            lambda[i][k] = lambda[i][k]/sum;            
      }
   
      /**************  Lambda(t,j,k) ***************/
      for(i=1;i<N;i++)
      {
      
         sum = 0.0;
         for(j=0;j<m;j++)
         {
            for(k=0;k<m;k++)
            {
               Lambda[i][j][k] = Q[j][k]*eta[i][j]*qbar[i+1][k]*g[k][i];
               sum += Lambda[i][j][k];
            }
         }
      
         for(j=0;j<m;j++)
         {
            for(k=0;k<m;k++)
               Lambda[i][j][k] = Lambda[i][j][k]/sum;
         }
      }
   
      for(j=0;j<m;j++)
      {
         for(k=0;k<m;k++)
            Lambda[N][j][k] = lambda[N][j]*Q[j][k];
      }
   
   
   
   
    /******************************************************************/
   
      free(W); 
      
   }


  
   void M_step_1d(double *R, int N, int m, double **eta, double **g, double **qbar, double **lambda, double ***Lambda, double *nu, double *mu, double *sigma, double **Q)
   
   {
   
      int i, j, k;
      double sum;
   
   
      /**** sum of lambda's  **************/
   
      for(k=0;k<m;k++)
      {
         sum = 0.0;
      
         for(i=0;i<N;i++)
            sum += lambda[i+1][k];
      
         nu[k] = sum;
      }
   
   
   /***********  mu ***************/
   
      for(k=0;k<m;k++)
      {
         sum = 0.0;
         for(i=0;i<N;i++)
            sum += R[i]*lambda[i+1][k];
      
         mu[k] = sum/nu[k];
      
      
      }
   
   
   /**** A ************************************/
      for(k=0;k<m;k++)
      {
         sum = 0.0;
         for(i=0;i<N;i++)
            sum += (R[i]-mu[k])*(R[i]-mu[k])*lambda[i+1][k];
      
         sigma[k] = sqrt(sum/nu[k]);
      
      
      
      
      
      }
   
   
   /********  Q  *******************************************/
   
      for(j=0;j<m;j++)
      {
         for(k=0;k<m;k++)
         {
            sum = 0.0;
            for(i=1;i<=N;i++)
               sum += Lambda[i][j][k];
         
            Q[j][k] = sum/nu[j];
         }
      
      }
   
   
   /************** nu again!! *************************/
   
      for(k=0;k<m;k++)
         nu[k] = nu[k]/((double) N);
   
   
   }


	
/**********************   Estimation of 1-d HMM  **********************/	

void   est_hmm_1d_new(double *R, int *pnstates, int *pN, double *peps, int *pmax_iter, double *mu, double *sigma, double *nu, double *Qvec, double *etavec, double *lambdavec, double *U, double *cvm, double *L) 
{
  
  int i, j, k, l;
  double *mu0, **g, **qbar, **lambda, ***Lambda, sum1, sum2, U0,sum,W, **Q, **eta;
  int  nstates, N, max_iter;
  double eps;
  

  nstates = (int) pnstates[0];
  eps = (double) peps[0];
  N = (int) pN[0];
  max_iter = (int) pmax_iter[0];
  
 
  Q  = (double **)calloc(nstates, sizeof(double *));
  for(j=0;j< nstates;j++)
    Q[j] = (double *)calloc(nstates, sizeof(double));
  
  
  
  eta    = (double **)calloc((N+1), sizeof(double *));
  for(k=0;k<=N;k++)
    eta[k] = (double *)calloc(nstates, sizeof(double));
  
  
  
  
  for(i=0;i<nstates*nstates;i++)
  {
    Qvec[i]=0.0;
  }
  
  for(i=0;i<nstates;i++)
  {
    nu[i]=0.0;
    mu[i]=0.0;
    sigma[i]=0.0;
  }
  
  g = (double **)calloc(nstates, sizeof(double *));
  for(j=0;j< nstates;j++)
    g[j] = (double *)calloc(N, sizeof(double));
  
  lambda = (double **)calloc((N+1),sizeof(double *));
  for(i=0;i<=N;i++)
    lambda[i] = (double *)calloc(nstates, sizeof(double));
  
  qbar = (double **)calloc((N+1), sizeof(double *));
  for(i=0;i<=N;i++)
    qbar[i] = (double *)calloc(nstates, sizeof(double));
  
  Lambda = (double ***)calloc((N+1), sizeof(double **));
  for(i=0;i<=N;i++)
    Lambda[i] = (double **)calloc(nstates, sizeof(double *));
  
  for(i=0;i<=N;i++)
  {
    for(k=0;k<nstates;k++)
      Lambda[i][k] = (double *)calloc(nstates, sizeof(double));
  }
  
  mu0 = (double *)calloc(nstates, sizeof(double));
  
  
  
  
 
  
  
  
  hmm_init_1d(R,N, nstates, mu0, sigma, Q);
  
 
  

 
  
  for(i=0;i<100;i++)
  {
    E_step_1d(R,mu0,sigma,Q,N,nstates,eta,g,qbar,lambda,Lambda,L);
  
    M_step_1d(R,N,nstates,eta,g,qbar,lambda,Lambda,nu,mu0,sigma,Q);
  }
  
 
  
 

  for(k=0;k<nstates;k++)
    mu[k] = mu0[k];
  
  
  
  for(i=0;i<max_iter;i++)
  {
    E_step_1d(R,mu0,sigma,Q,N,nstates,eta,g,qbar,lambda,Lambda,L);
    
    M_step_1d(R,N,nstates,eta,g,qbar,lambda,Lambda,nu,mu,sigma,Q);
    
    
    sum1 = 0.0; sum2=0.0;
    for(k=0;k<nstates;k++)
    { 
      sum1 += fabs(mu0[k]);
      sum2 += fabs(mu[k]-mu0[k]);
    }
    if( sum2 < nstates*eps*sum1)
      break;
    
    for(k=0;k<nstates;k++)
      mu0[k] = mu[k];
    
  }
  
  E_step_1d(R,mu,sigma,Q,N,nstates,eta,g,qbar,lambda,Lambda,L);
  
 
 /* etavec and lambdavec*/

 i = 0;
  for(k=0;k<nstates;k++)
  {
    for(j=1;j<=N;j++)
    {
      etavec[i] = eta[j][k];
      lambdavec[i] = lambda[j][k];
      
      i ++;
    }
  }
  
 /* Qvec */
  
   i = 0;
  for(k=0;k<nstates;k++)
  {
    for(j=0;j<nstates;j++)
    {
      Qvec[i] = Q[j][k];
      
      i ++;
    }
  }
  
  
  
  
  for(i=0;i<N;i++)
  {
    sum = 0.0;
    
    for(k=0;k<nstates;k++)
    {
      W = 0.0;
      for(l=0;l<nstates;l++)
        W += eta[i][l]*Q[l][k];
      
      U0 = Phi((R[i]-mu[k])/sigma[k]);
      
      sum += W*U0;
    }
    
    U[i]=sum;
    
    
    
  }
  
  cvm[0] = Sn_1d(U,N);
  
  
  
   	
  for(j=0;j< nstates;j++)
    free(g[j]); 
  
  free(g);
  
  
  for(i=0;i<=N;i++)
    free(lambda[i]);
  
  free(lambda);
  
  
  for(i=0;i<=N;i++)
    free(qbar[i]); 
  
  free(qbar);
  
  
  
  for(i=0;i<=N;i++)
  {
    for(k=0;k<nstates;k++)
      free(Lambda[i][k]);
  }
  
  for(i=0;i<=N;i++)
    free(Lambda[i]);
  
  free(Lambda);
  
  free(mu0);
  
}
