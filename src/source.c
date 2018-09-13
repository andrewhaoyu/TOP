#include "./source.h"

/* Function to check the stopping criteria */
double checkStop(delta, delta0, n)
double *delta, *delta0;
int n;
{
  int i;
  double maxv=-9999.9, temp, rerror=-9999.9, *pd, *pd0, eps=1.0e-8, d1, d0;
  double reldiff, absdiff;

  /* Get max value */
  for (i=0, pd0=delta0; i<n; i++, pd0++) {
    temp = fabs(*pd0);
    if (temp > maxv) maxv = temp;
  }
  if (maxv < 0.1) maxv = 0.1;
  for (i=0, pd=delta, pd0=delta0; i<n; i++, pd++, pd0++) {
    temp = fabs(*pd - *pd0)/maxv;
    if (temp > rerror) rerror = temp;
  }
  
  /*
  for (i=0, pd=delta, pd0=delta0; i<n; i++, pd++, pd0++) {
    d1      = *pd;
    d0      = *pd0;
    absdiff = fabs(d1 - d0);
    reldiff = absdiff/(fabs(d0) + eps);
    temp    = fmin(absdiff, reldiff);
    if (temp > rerror) rerror = temp;
  }
  */

  return(rerror);

} /* END: checkStop */



/*calculate the weighted matrix for MLE*/
/* W is a long vector with length N*(M+1)*M/2*/
/* since W is sysmetric , we only need to record the lower triangle*/
/* j<==i */
/* use the sysmetric to fill in the other side */
void Weighted_W(double *p, double *W,int N,int M){
    int i, j, k, NM, NMi, NMj, Nj, Ni, NMiNi;
    double *pti, *ptj, *pWi, *pWj, tmp;

    NM = M*N;

    /*for(int i=0;i<M;i++){
      for(int j=0;j<M;j++){
        if(i==j){
          for(int k=0;k<N;k++){
            W[NM*i+N*j+k] = p[N*i+k]-p[N*i+k]*p[N*i+k];
          }
        }else{
          for(int k=0;k<N;k++){
            W[NM*i+N*j+k] = -p[N*i+k]*p[N*j+k];
          }
        }

      }
    }*/

      for(i=0; i<M; i++){
        NMi = NM*i;
        Ni  = N*i;
        for(j=0; j<i; j++){
          NMj = NM*j;
          Nj  = N*j;
          pWi = &W[NMi+Nj];
          pWj = &W[NMj+Ni];
          for(k=0, pti=&p[Ni], ptj=&p[Nj]; k<N; k++, pti++, ptj++, pWi++, pWj++){
            tmp  = - *pti * *ptj;
            *pWi = tmp;
            *pWj = tmp; /* use the sysmetric of W */
          }
        }

        NMiNi = NMi + Ni;
        for(k=0, pti=&p[Ni], pWi=&W[NMiNi]; k<N; k++, pti++, pWi++){
          tmp  = *pti;
          *pWi = tmp - tmp*tmp;
        }
      }

} /* END: Weighted_W */


double LogLikelihood(double *Y, double *ret_p, int N, int M){
 
  double *Y_temp, *p_temp, loglike=0.0, p_o, p_sum, *pY, *pp;
  int i, j, ind, k;

  Y_temp = dVec_alloc(M,0,0.0);
  p_temp = dVec_alloc(M,0,0.0);

  for(i=0;i<N;i++){
    p_o   = 0.0;
    ind   = 0;
    p_sum = 0.0;
    for(j=0, pp=p_temp, pY=Y_temp; j<M; j++, pp++, pY++) {
      k   = i+N*j;
      *pp = ret_p[k];
      *pY = Y[k];

      if(*pY > ALMOST_ZERO){
        p_o += *pp;
        ind  = 1;
      }
    }
    
    if (!ind){
      for(j=0, pp=p_temp;j<M; j++, pp++) p_sum += *pp;
      p_o = 1.0 - p_sum;
    }

    loglike += log(p_o);
  }
  free(Y_temp);
  free(p_temp);
  return(loglike);

} /* END: LogLikelihood */


void myconvert(int *n, double *p, double *snpvalue){
  int i;
  for(i=0; i<*n; i++){
    snpvalue[i] = p[(3*i+1)]+p[(3*i+2)]*2;
  }
  return;
}


/* Function to check for non-finite values */
int all_finite(vec, n)
double *vec;
int n;
{
  int i;
  double *pd;

  for (i=0, pd=vec; i<n; i++, pd++) {
    if (!R_FINITE(*pd)) return(0);
  }

  return(1);

} /* END: all_finite */

