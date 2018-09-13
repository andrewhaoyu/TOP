#include "./source.h"


/* Function to compute lxx = xx * beta = xx * z * delta */
static void get_lxx(Z, Znr, Znc, delta, X, Xnr, Xnc, M,out)
double **Z, *delta, **X, *out;
int Xnr, Xnc, Znr, Znc, M;
{
  double  sum,*beta, *pb, **prow, *px;
  int i, j, k, row, brow, bstart;
  beta  = dVec_alloc(Znr, 0, 0.0);

  /* Compute beta = Z*delta */
  X_y(Z, Znr, Znc, delta, beta);

  /* Compute out = XX*beta */
  row = 0;
  for (i=0; i<M; i++) {
      bstart = i*Xnc;
      for (j=0, prow=X; j<Xnr; j++, prow++) {
        sum  = 0.0;
        brow = bstart;
        for (k=0, px=*prow, pb=&beta[brow]; k<Xnc; k++, px++, pb++) sum += *px * *pb;  
        out[row] = sum;
        row++;
      }
  }
  free(beta);

} /* END: get_lxx */

/* Function to compute pxx */
static void get_pxx(lxx, N, M, out)
double *lxx, *out;
int N, M;
{
  int i, j, col, NM;
  double sum, *po, *pl;

  NM = N*M;
  for (i=0, po=out, pl=lxx; i<NM; i++, po++, pl++) *po = exp(*pl);

  /* Scale */
  for (i=0; i<N; i++) {
        sum = 0.0;
        col = i;
        for (j=0; j<M; j++) {
          sum += out[col];
          col += N;
        }
        col = i;
        sum = sum + 1.0;
        for (j=0; j<M; j++) {
          out[col] = out[col]/sum;
          col += N;
        }
  }

} /* END: get_pxx */

/* Function to compute W_y = yy-pxx+W%*%lxx */
static void get_Wy(Y, lxx, Pxx, N, M,W, out)
double *Y;  /* Nsub*M vector of outcomes */
int N, M;
double *lxx, *out, *Pxx, *W;
{
  int NM, k,i, j, t, NMik, Nj;
  double sum, *pY, *po, *pP;

  NM = N*M;

  for(t=0, po=out, pY=Y, pP=Pxx; t<NM; t++, po++, pY++, pP++){
      i    = t/N; /*t/N round to the nearst smaller integer */
      k    = t%N;
      NMik = NM*i + k;
      sum  = 0.0;
      for(j=0; j<M; j++){
        Nj   = N*j;
        sum += lxx[Nj+k]*W[NMik+Nj];
      }
      *po = *pY - *pP + sum;
  }


  /* for (row=0, p1=Pxx, p2=out, py=Y; row<NM; row++, p1++, p2++, py++) {
    prow = *p1;
    sum  = (prow-prow*prow)*lxx[row];
    ii   = row + N;
    jj   = row - N;
    for (i=2; i<MP1; i++) {
      if (ii < NMP1) {
        sum += -prow*Pxx[ii]*lxx[ii];
        ii   = ii + N;
      }
      if (jj > -1) {
        sum += -prow*Pxx[jj]*lxx[jj];
        jj   = jj - N;
      }
    }
    *p2 = *py - prow + sum;
  } */

} /* END: get_Wy */


  /* Function to compute delta */
static void get_delta(INV, nc, tXXZ, N, M, W_y, out)
double **INV, **tXXZ, *W_y, *out;
int nc, N, M;
{
  /* delta <- solve(Infor_M,t(xxz)%*%W_y) */
    int i, NM;
  double *p1, **p2, *tXXZWy;

  NM     = N*M;
  tXXZWy = dVec_alloc(nc, 0, 0.0);

  /* Compute t(xxz)%*%w_y */
    for (i=0, p1=tXXZWy; i<nc; i++, p1++) {
      *p1 = dotProd(tXXZ[i], W_y, NM);
    }

  /* Compute delta */
    for (i=0, p1=out, p2=INV; i<nc; i++, p1++, p2++) *p1 = dotProd(*p2, tXXZWy, nc);

    free(tXXZWy);

} /* END: get_delta */
 
  /*************************************************************************************/
  /*************************************************************************************/
  /*************************************************************************************/


/* Function Mvpoly */
void Mvpoly(double *deltai,int Nparm,double* Y,
                     double **X, double **Z,int Znr,int Znc,
                     int N, int M, int Ncov, int Niter,double tol,
                     int DEBUG,int* ret_rc,double* ret_delta,double**Info,
                     double *ret_p,double *lxx,double *W,double *beta,
                     double* w_y,double **XmWXm,double **Inv,double **tXXZ,double *XX){
    int iter, NM, rc, conv=0;
    double  rerror,*delta0;

    NM      = N*M;
    *ret_rc = 1;

      /* Allocate memory for matrix of covariates and Z map matrix */
      delta0   = dVec_alloc(Nparm, 0, 0.0);
      /* Copy initial estimates to delta0 */
      copy_dVec(delta0, deltai, Nparm);
      

        /*print_dMat(tXXZ,(Ncov-1)*Ncatp1+M,M*N,tXXZ);*/
        if (DEBUG) Rprintf("Begin Weighted Least Sqaure Maximization\n");
        for (iter=1; iter<=Niter; iter++) {
          if (DEBUG) Rprintf("Weighted Least Sqaure Iteration: %d\n", iter);

          if (DEBUG) Rprintf("Compute lxx\n");
          get_lxx(Z, Znr, Znc, delta0, X, N, Ncov, M,lxx);
          if (DEBUG) Rprintf("Compute pxx\n");
          get_pxx(lxx, N, M, ret_p);
          if (DEBUG) Rprintf("Compute weighted matrix\n");
          Weighted_W(ret_p, W, N, M);
          if (DEBUG) Rprintf("Compute XmWXm matrix\n");
          get_XmWXm(XX,X,W, M,N, Ncov,XmWXm);
          /*print_dMat(XmWXm,Znr,Znr,"XmWXm");*/
            if (DEBUG) Rprintf("Compute information matrix\n");

          /*get_Info(X, N, Ncov, M, Z, Znr, Znc, pxx, Info);*/

            QuadXKX(Z,XmWXm, Znr,Znc, Info);
          /* print_dMat(XmWXm,Znr,Znr,"Info");*/

            if (DEBUG) Rprintf("Compute covariance matrix\n");
          rc = cov_inv(Info, Znc, Inv);
          if (rc) {
            Rprintf("ERROR computing inverse of information matrix\n");
            free(delta0);
            return;
          }

          if (DEBUG) Rprintf("Compute W_y\n");
          get_Wy(Y, lxx, ret_p, N, M,W,w_y);
          if (DEBUG) Rprintf("Compute delta\n");

          get_delta(Inv, Nparm, tXXZ, N, M, w_y, ret_delta);

          if (DEBUG) print_dVec(ret_delta, Nparm, "delta");

          /* Check for non-finite values */
            if (!all_finite(ret_delta, Nparm)) {
              Rprintf("ERROR: Weighted Least Square Algorithm not converging, parameters have non-finite values\n");
              free(delta0);
              return;
            }

          if (DEBUG) Rprintf("Check Weighted Least Square stopping criteria\n");
          rerror = checkStop(ret_delta, delta0, Nparm);
          if (rerror <= tol) {
            conv = 1;
            break;
          }

          /* Update delta0 */
          if (DEBUG) Rprintf("Update parameters\n");
          copy_dVec(delta0, ret_delta, Nparm);
        }

        free(delta0);

        if (!conv) {
          Rprintf("ERROR: algorithm did not converge\n");
          return;
        }
        *ret_rc = 0;

  return;

} /* END: Mvpoly */

