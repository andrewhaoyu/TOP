#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <R.h>
#include <R_ext/Memory.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define ALMOST_ZERO 1e-16
#define NUM_ZERO 1e-100
#define ERROR_SINGULAR_MATRIX 1
#define CHECK_MEM(obj) if (obj == NULL) {Rprintf("ERROR: allocating memory \n"); error("1");}

static void print_dVec(vec, n, name)
double *vec;
int n;
char name[10];
{
  int i;
  printf("%s \n", name);
  for (i=0; i<n; i++) {
    printf(" %g ", vec[i]);
  }
  printf("\n \n");
}
static void print_iVec(vec, n, name)
int *vec;
int n;
char name[10];
{
  int i;
  printf("%s \n", name);
  for (i=0; i<n; i++) {
    printf(" %d ", vec[i]);
  }
  printf("\n \n");
}


static void print_dMat(mat, nr, nc, name)
double **mat;
int nr, nc;
char name[10];
{
  int i, j;
  printf("%s \n", name);
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) printf(" %g ", mat[i][j]);
    printf("\n");
  }
  printf("\n \n");
}

/* Function to allocate memory for a double vector */
  static double * dVec_alloc(n, initFlag, initVal)
int n, initFlag;
double initVal;
{
  int i;
  double *ret, *p;

  ret = (double *) malloc(n*sizeof(double));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }

  return(ret);

} /* END: dVec_alloc */

  /* Function to allocate a double matrix */
  static double ** dMat_alloc(nrow, ncol, initFlag, initVal)
int nrow, ncol, initFlag;
double initVal;
{
  double **mat, **ptr;
  int i;

  mat = (double **) malloc(nrow*sizeof(double *));
  CHECK_MEM(mat);
  for (i=0, ptr=mat; i<nrow; i++, ptr++) *ptr = dVec_alloc(ncol, initFlag, initVal);

  return(mat);

} /* END: dMat_alloc */

  /* Function to free a matrix */
  static void matrix_free(x, n)
void **x;
int n;
{
  int i;
  for (i=0; i<n; i++) {
    if (x[i]) free(x[i]);
  }
  free(x);

} /* END: matrix_free */

  /* Function to check for non-finite values */
  static int all_finite(vec, n)
double *vec;
int n;
{
  int i;

  for (i=0; i<n; i++) {
    if (!R_FINITE(vec[i])) return(0);
  }

  return(1);

} /* END: all_finite */

  /* Multiply to matrices (matrix mult)  */
  static void matrixMult(m1, m1_nr, m1_nc, m2, m2_nc, ret)
double **m1, **m2, **ret;
int m1_nr, m1_nc, m2_nc;
{
  int i, j, k;
  double sum;

  for (i=0; i<m1_nr; i++) {
    for (j=0; j<m2_nc; j++) {
      sum = 0.;
      for (k=0; k<m1_nc; k++) sum += m1[i][k]*m2[k][j];
      ret[i][j] = sum;
    }
  }

} /* END: matrixMult */


  /* two matrix minus  */
  static void MatrixMinus(double **mat1,double**mat2,int nr,int nc,double **ret){
    for(int i=0;i<nr;i++){
      for(int j=0;j<nc;j++){
        ret[i][j] = mat1[i][j]-mat2[i][j];
      }
    }
  }/* END: matrixminus */


  /* two vector minus */
  static void VectorMinus(double *vec1,double*vec2,int N,double *ret){
    for(int i=0;i<N;i++){
      ret[i] = vec1[i]-vec2[i];
    }
  }



/* Function to compute Xy for matrix X and vector */
  static void X_y(X, nr, nc, y, ret)
double **X, *y, *ret;
int nr, nc;
{
  int i, j;
  double sum, *p, *pret, *px;

  for (i=0, pret=ret; i<nr; i++, pret++) {
    sum  = 0.0;
    for (j=0, p=y, px=X[i]; j<nc; j++, p++, px++) {
      sum += *px * *p;
    }
    *pret = sum;
  }

} /* END: X_y */

  /* Function for dot product of two vectors */
  static double dotProd(v1, v2, n)
double *v1, *v2;
int n;
{
  int i;
  double sum=0.0;

  for (i=0; i<n; i++) sum += v1[i]*v2[i];

  return(sum);

} /* END: dotProd */

  /* Function to compute Z */
  /*static void get_Z(Z_design, M, Ncatp1, Ncov, out)
double **Z_design, **out;
int M, Ncatp1, Ncov;
{
  int i, j, k, row=0, col;
  double *vec;


  Get the column for diag(P) (X) z_design, loop over rows of z_design
  for(i = 0;i<M;i++){
    row = 0+ i*Ncov;
    col = i;
    out[row][col] = 1.0;
  }
  for(i = 0;i<M;i++){
    vec = Z_design[i];
    for(j =0; j<(Ncov-1);j++){
      row = j+1+i*(Ncov);
      col = M+j*Ncatp1;
      for(k=0;k<Ncatp1;k++){
        out[row][col] = vec[k];
        col++;
      }
    }

  }

}*/

  /* Get the column for diag(P) (X) z_design, loop over rows of z_design */
  /*for(i=0; i<M; i++) {
    vec = Z_design[i];
    for (j=0; j<Ncov; j++) {
      col = j*Ncatp1;
      for (k=0; k<Ncatp1; k++) {
        out[row][col] = vec[k];
        col++;
      }
      row++;
    }
  }

}  END: get_Z */

  /* Function to get a column of XXZ matrix */
  static void get_XXZ_col(col, X, Xnr, Xnc, M, Z, dtempvec, out)
double **X, **Z, *out, *dtempvec; /* dtempvec length M*Xnc */
  int col, Xnr, Xnc, M;
{
  int  i, j, k, zstart;
  double sum, *pt, *prow, *pout;

  /* Copy column of Z into dtempvec */
    for (i=0, pt=dtempvec; i<Xnc*M; i++, pt++) *pt = Z[i][col];

    pout = out;
    for (i=0; i<M; i++) {
      zstart = i*Xnc;
      for (j=0; j<Xnr; j++) {
        sum  = 0.0;
        for (k=0, prow=X[j], pt=&dtempvec[zstart]; k<Xnc; k++, prow++, pt++) {
          sum += *prow * *pt;
        }
        *pout++ = sum;
      }
    }

} /* END: get_XXZ_col */

  /* Function to compute t(XXZ) */
  static void get_tXXZ(X, Xnr, Xnc, M, Znc, Z, out)
double **X, **Z, **out;  /* out must have dim Zallnc x Xnr*M */
  int Xnr, Xnc, M,Znc;
{
  int i;
  double *temp;

  temp = dVec_alloc(Xnc*M, 0, 0.0);

  for (i=0; i< Znc; i++) {
    get_XXZ_col(i, X, Xnr, Xnc, M, Z, temp, out[i]);
  }

  free(temp);

} /* END: get_tXXZ */

  /* Function to compute lxx = xx * beta = xx * z * delta */
  static void get_lxx(Z, Znr, Znc, delta, X, Xnr, Xnc, M,out)
double **Z, *delta, **X, *out;
int Xnr, Xnc, Znr, Znc, M;
{
  double  sum,*beta;
  int i, j, k, row, brow, bstart;
  beta  = dVec_alloc(Znr, 0, 0.0);

  /* Compute beta = Z*delta */
  X_y(Z, Znr, Znc, delta, beta);

  /* Compute out = XX*beta */
    row = 0;
    for (i=0; i<M; i++) {
      bstart = i*Xnc;
      for (j=0; j<Xnr; j++) {
        sum  = 0.0;
        brow = bstart;
        for (k=0; k<Xnc; k++) {
          sum += X[j][k]*beta[brow];
          brow++;
        }
        out[row] = sum;
        row++;
      }
    }



} /* END: get_lxx */

  /* Function to compute pxx */
  static void get_pxx(lxx, N, M, out)
double *lxx, *out;
int N, M;
{
  int i, j, col, NM;
  double sum;

  NM = N*M;
  for (i=0; i<NM; i++) out[i] = exp(lxx[i]);

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

  /* Function to compute vec1*W*vec2 */
  static double v1Wv2(p, N, M, vec1, vec2)
double *p, *vec1, *vec2;  /* p is stored as a vector, out must be of length NM */
  int N, M;
{
  int i, ii, jj, NM, MP1, row, NMP1;
  double sum, prow, *p1, *pv2, *pv1, ret;

  NM   = N*M;
  MP1  = M + 1;
  NMP1 = NM + 1;

  ret = 0.0;
  for (row=0, p1=p, pv2=vec2, pv1=vec1; row<NM; row++, p1++, pv2++, pv1++) {
    prow = *p1;
    sum  = (prow-prow*prow)* *pv2;
    ii   = row + N;
    jj   = row - N;
    for (i=2; i<MP1; i++) {
      if (ii < NMP1) {
        sum += -prow*p[ii]*vec2[ii];
        ii   = ii + N;
      }
      if (jj > -1) {
        sum += -prow*p[jj]*vec2[jj];
        jj   = jj - N;
      }
    }
    ret += *pv1 * sum;
  }

  return(ret);

} /* END: v1Wv2 */



  /* Function to compute W_y = yy-pxx+W%*%lxx */
  static void get_Wy(Y, lxx, Pxx, N, M,W, out)
double *Y;  /* Nsub*M vector of outcomes */
  int N, M;
double *lxx, *out, *Pxx, *W;
{
  int  NM, k,i;
  double sum;

  NM   = N*M;

  for(int t =0; t <NM; t++){
    i = t/N; /*t/N round to the nearst smaller integer */
      k = t%N;
      sum =0.0;
      for(int j=0;j<M;j++){
        sum += lxx[N*j+k]*W[i*NM+j*N+k];
      }
      out[t] = Y[t]-Pxx[t]+sum;
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


  /* fill the info matrix to the result*/
  /* the Info matrix is sysmetric */
  static void fill_SysMat(Mat,Vec,Nr)
double **Mat,*Vec;
int Nr;
{
  int i,j;
  for(j=0;j<Nr;j++){
    for(i=0;i<(j+1);i++){
      Vec[j*Nr+i] = Mat[i][j];
      Vec[i*Nr+j] = Mat[i][j];
    }
  }


}

/*end fill_SysMat*/

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


  /* Function to check the stopping criteria */
  static double checkStop(delta, delta0, n)
double *delta, *delta0;
int n;
{
  int i;
  double maxv=-9999.9, temp, rerror;

  /* Get max value */
    for (i=0; i<n; i++) {
      temp = fabs(delta0[i]);
      if (temp > maxv) maxv = temp;
    }

  if (maxv < 0.1) maxv = 0.1;

  rerror = -9999.9;
  for (i=0; i<n; i++) {
    temp = fabs(delta[i] - delta0[i])/maxv;
    if (temp > rerror) rerror = temp;
  }

  return(rerror);

} /* END: checkStop */

  /* Function to fill in a matrix from a vector (by column) */
  static void fillMat(vec, nr, nc, addInt, out)
double *vec, **out;
int nr, nc, addInt;
{
  int i, j, col=0, ii;

  if (addInt) {
    /* Intercept for first column */
      for (i=0; i<nr; i++) out[i][0] = 1.0;
      col = 1;
  }

  ii = 0;
  for (j=0; j<nc; j++) {
    for (i=0; i<nr; i++) {
      out[i][col] = vec[ii];
      ii++;
    }
    col++;
  }

} /* END: fillMat */

  /***********************************************************************************/
  /* For computing inverse */
  /***********************************************************************************/

  /* Factor a symmetric positive definite matrix */
  static int symPosFactor(mat, n, ret, retdiag)
double **mat, **ret, *retdiag;
int n;
{
  int i, j, k;
  double sum, save, *ptr;

  /* Copy mat to ret */
    for (i=0; i<n; i++) {
      for (j=0; j<n; j++) ret[i][j] = mat[i][j];
    }

  for (i=0, ptr=retdiag; i<n; i++, ptr++) {
    for (j=i; j<n; j++) {
      sum = ret[i][j];
      for (k=i-1; k>-1; k--) sum -= ret[i][k]*ret[j][k];
      if (i == j) {
        if (sum < NUM_ZERO) return(ERROR_SINGULAR_MATRIX);
        save = sqrt(sum);
        *ptr = save;
      } else {
        ret[j][i] = sum/save;
      }
    }
  }

  /* Zero out the diagonal and above */
    for (i=0; i<n; i++) {
      for (j=i; j<n; j++) ret[i][j] = 0.;
    }

  return(0);

} /* END: symPosFactor */

  /* Invert a factor */
  static void symPosFacInv(L, diag, n, ret)
double **L, *diag;
int n;
double **ret;
{
  int i, j, k;
  double sum;

  /* Copy L to ret */
    for (i=0; i<n; i++) {
      for (j=0; j<n; j++) ret[i][j] = L[i][j];
    }

  for (i=0; i<n; i++) {
    ret[i][i] = 1./diag[i];
    for (j=i+1; j<n; j++) {
      sum = 0.;
      for (k=i; k<=j-1; k++) sum -= ret[j][k]*ret[k][i];

      ret[j][i] = sum/diag[j];
    }
  }

} /* END: symPosFacInv */

  /* Transpose of a square matrix */
  static void matTranspose(mat, n, ret)
double **mat;
int n;
double **ret;
{
  int i, j;

  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) ret[j][i] = mat[i][j];
  }

} /* END: rmatTranspose */

  /* Inverse of symmetric positive definite matrix */
  static int symPosMatInv(mat, n, ret)
double **mat;
int n;
double **ret;
{
  double **L, **Linv, *diag;
  int i;

  L = dMat_alloc(n, n, 1, 0.0);
  diag = (double *) dVec_alloc(n, 1, 0.0);

  /* Get cholesky lower triangle */
    i = symPosFactor(mat, n, L, diag);
    if (i) {
      free(diag);
      matrix_free((void **) L, n);
      return(i);
    }

    /* Get the lower inverse */
      Linv = dMat_alloc(n, n, 1, 0.0);
      symPosFacInv(L, diag, n, Linv);
      free(diag);

      /* Get the transpose of Linv */
        matTranspose(Linv, n, L);

      /* Inverse is t(L^-1)(L^-1) */
        matrixMult(L, n, n, Linv, n, ret);

      matrix_free((void **) L, n);
      matrix_free((void **) Linv, n);

      return(0);

} /* END: symPosMatInv */




  /*calculate the weighted matrix for MLE*/
  /* W is a long vector with length N*(M+1)*M/2*/
  /* since W is sysmetric , we only need to record the lower triangle*/
  /* j<==i */
  /* use the sysmetric to fill in the other side */
  static void Weighted_W(double *p, double *W,int N,int M){
    int NM;
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
      for(int i=0;i<M;i++){
        for(int j=0;j<(i);j++){
          for(int k=0;k<N;k++){
            W[NM*i+N*j+k] = -p[N*i+k]*p[N*j+k];
            W[NM*j+N*i+k] = W[NM*i+N*j+k]; /* use the sysmetric of W */

          }

        }
      }
    for(int i=0;i<M;i++){
      for(int k=0;k<N;k++){
        W[NM*i+N*i+k] = p[N*i+k]-p[N*i+k]*p[N*i+k];
      }
    }

  }

/* calculate the vector combination XX*/
  /* X = (X1,..Xp) */
  /* We save X1*X1T for later calculation of XmWXm */
  /* XX is saved as a N*p^2 vector */
  /* Use the sysmetric of XX property */
  static void get_XX(double **X,int Ncov,int N,double *XX){
    for(int i=0;i<Ncov;i++){
      for(int j=0;j<(i+1);j++){
        for(int k=0;k<N;k++){
          XX[i*N*Ncov+j*N+k] = X[k][i]*X[k][j];
          XX[j*N*Ncov+i*N+k] = XX[i*N*Ncov+j*N+k];
        }
      }
    }

  }



/*calculate the quadratic form  XmWXm */
  /* Since Xm is are sparse, we take the advantage of this property */
  /*XmWXm can be decomposed into M*M blcck, each block is p*p matrix*/
  /* i j represent the row and column for M*M block */
  /* k,l represent the row and column within the p*p matrix*/
  /* p*p matrix can be got by X*(W with N*N block)*X */
  /* Since W with N*N block is diagnonal matrix */
  /* X*(W with N*N block)*X is weigthed qudractic sum, t loop used for that */
  /* Use the sysmetric property of XmWXm */

  static void get_XmWXm(double *XX,double **X,double *W, int M,int N, int p,double **ret){
    int NM = N*M;
    double *Wtemp;
    Wtemp   = dVec_alloc(N,0, 0.0);
    double sum =0.0;

    for(int i=0;i<M;i++){
      for(int j=0;j<(i+1);j++){
        for(int v=0;v<N;v++){   /*fillin Wtemp for W with N*N block diagnals*/
            Wtemp[v] = W[NM*i+N*j+v];
        }
        for(int k=0;k<p;k++){
          for(int l=0;l<(k+1);l++){
            sum =0.0;
            for(int t=0;t<N;t++){
              sum += XX[k*N*p+l*N+t]*Wtemp[t];
            }
            ret[i*p+k][j*p+l] = sum;
            ret[i*p+l][j*p+k] = sum;
            ret[j*p+k][i*p+l] = sum;
            ret[j*p+l][i*p+k] = sum;
          }
        }

      }
    }
    free(Wtemp);
  }

/* Function for quadractic computation X^tkX */
  static void QuadXtKX(double **X,double ** K, int Xnr,int Xnc, double** ret){
    double sum = 0.0;

    for(int i=0;i<Xnc;i++){
      for(int j=0;j<(i+1);j++){
        /*the ith row and jth column of the ret*/
          /*the ith column X transpose times K times the jth column of the X */
          sum = 0.0;
          /* One vector times K times one Vector*/
            for(int k=0;k<Xnr;k++){
              for(int l=0;l<Xnr;l++){

                sum += X[k][i]*K[k][l]*X[l][j];
              }
            }
          ret[i][j] = sum;
          /* ret is sysmetric */
            ret[j][i] = sum;

      }
    }
  }

/* Function for quadractic computation X^tkX */
  static void QuadXKXt(double **X,double ** K, int Xnr,int Xnc, double** ret){
    double sum = 0.0;

    for(int i=0;i<Xnr;i++){
      for(int j=0;j<(i+1);j++){
        /*the ith row and jth column of the ret*/
          /*the ith row X transpose times K times the jth row of the X */
          sum = 0.0;
          /* One vector times K times one Vector*/
            for(int k=0;k<Xnc;k++){
              for(int l=0;l<Xnc;l++){

                sum += X[i][k]*K[k][l]*X[j][l];
              }
            }
          ret[i][j] = sum;
          /* ret is sysmetric */
            ret[j][i] = sum;

      }
    }
  }




/* Function for getting the observed weighted matrix  W_com-W_com|mis*/
  static void Get_ObservedW(double *W,double *W_obs,int N,int M, double *Y){
    double * W_mis  ;
    double **Info_mis, **XmWmisXm;
    W_mis = dVec_alloc((M*M*N),0,0.0);
    Weighted_W(Y, W_mis, N, M);
    VectorMinus(W,W_mis,N*M*M,W_obs);
    free(W_mis);
  }














/* Function for getting the observed information matrix */
  static void Get_ObservedInfo(int M,int N,double **Info_obs,int DEBUG,
                               double *XX,double **X,int Ncov,int Znr,int Znc,double **Z,double*W_obs){
    double **XmWobsXm;
    XmWobsXm    = dMat_alloc(Znr,Znr,0,0.0);
  /*  if (DEBUG) Rprintf("Compute XmWobsXm matrix\n");*/
    get_XmWXm(XX,X,W_obs, M,N, Ncov,XmWobsXm);
    /*if (DEBUG) Rprintf("Compute mis information matrix\n");*/
    /*get_Info(X, N, Ncov, M, Z, Znr, Znc, pxx, Info);*/
      QuadXtKX(Z,XmWobsXm, Znr,Znc, Info_obs);
    matrix_free((void **)XmWobsXm,Znr);

  }


/* calculate WX is a matrix with N*M row and M*p column*/
  /* WX is decomposed into M*M block, each block has N*p column*/
  /* Decompose W matrix into M*M block. each block is a diagnonal matrix with N*N */
  /* i represent the row for the M*M block */
  /* j represent the column for the M*M block */
  /* k represent 1 to Ncov column of X */
  /* v represent 1 to N element of X kth column element */
  /* W is sysmetric */

  static void get_WX(int M,int N,int Ncov,double **X, double *W,double **WX){
    int NM =N*M;
    double *Wtemp;
    Wtemp = dVec_alloc(N,0,0.0);
    for(int i=0;i<M;i++){
      for(int j=0;j<(i+1);j++){
        for(int v=0;v<N;v++){
          Wtemp[v] = W[NM*i+N*j+v];
        }

        for(int k=0;k<Ncov;k++){
          for(int v=0;v<N;v++){
            /*WX[N*i+v][Ncov*j+k] = */
              WX[N*j+v][Ncov*i+k] =  WX[N*i+v][Ncov*j+k]=Wtemp[v]*X[v][k];

          }
        }
      }
    }
    free(Wtemp);
  }

/* fill matrix into vector */
  /* vec was ordered by column */
  static void fill_vec(double **mat,int nr, int nc, double *ret){
    for(int j=0;j<nc;j++){
      for(int i=0;i<nr;i++){
        ret[j*nr+i] = mat[i][j];
      }
    }
  }

/* get the matrix WXZ */
  static void get_WXZ(int M,int N,int Ncov,double **X, double *W,
                      double **WX,double **WXZ,int Nparm,double **Z,int Znc,
                      double *WXZ_vec,double *WX_vec){
    /*Rprintf("Get WX");*/
      get_WX(M,N,Ncov,X, W, WX);
    fill_vec(WX,N*M,M*Ncov,WX_vec);
    /*Rprintf("Get WXZ");*/
      matrixMult(WX, N*M, M*Ncov, Z, Nparm, WXZ);
    fill_vec(WXZ,N*M,Znc,WXZ_vec);

  }



/*calculate the weighted matrix for MLE*/

  /*************************************************************************************/
  /*************************************************************************************/
  /*************************************************************************************/


  /* Function Mvpoly */
  static void Mvpoly(double *deltai,int Nparm,double* Y,
                     double **X, double **Z,int Znr,int Znc,
                     int N, int M, int Ncov, int Niter,double tolMStep,
                     int DEBUG,int* ret_rc,double* ret_delta,double**Info,
                     double *ret_p,double *lxx,double *W,double *beta,
                     double* w_y,double **XmWXm,double **Inv,double **tXXZ,double *XX){
    int i, iter, NM, rc, conv=0;
    double  rerror,*delta0;
    NM      = N*M;
    /*printf("Ncat\tNcatp1\tNcov0\tNcov\tZnr\Znc\n");
    printf("%i\t%i\t%i\t%i\t%i\t%i\n",Ncat,Ncatp1,Ncov0,Ncov,Znr,Znc);*/

      /* Allocate memory for matrix of covariates and Z map matrix */
      delta0   = dVec_alloc(Nparm, 0, 0.0);
      /* Copy initial estimates to delta0 */
        for (i=0; i<Nparm; i++) delta0[i] = deltai[i];
        /*print_dMat(tXXZ,(Ncov-1)*Ncatp1+M,M*N,tXXZ);*/
          /*  if (DEBUG) Rprintf("Begin Weighted Least Sqaure Maximization\n");*/
        for (iter=1; iter<=Niter; iter++) {
          /*if (DEBUG) Rprintf("Weighted Least Sqaure Iteration: %d\n", iter);*/

          /*if (DEBUG) Rprintf("Compute lxx\n");*/
          get_lxx(Z, Znr, Znc, delta0, X, N, Ncov, M,lxx);
          /*if (DEBUG) Rprintf("Compute pxx\n");*/
          get_pxx(lxx, N, M, ret_p);
          /*if (DEBUG) Rprintf("Compute weighted matrix\n");*/
          Weighted_W(ret_p, W, N, M);
          /*if (DEBUG) Rprintf("Compute XmWXm matrix\n");*/
          get_XmWXm(XX,X,W, M,N, Ncov,XmWXm);
          /*if (DEBUG) Rprintf("Compute information matrix\n");*/
          /*get_Info(X, N, Ncov, M, Z, Znr, Znc, pxx, Info);*/
            QuadXtKX(Z,XmWXm, Znr,Znc, Info);
          /*if (DEBUG) Rprintf("Compute covariance matrix\n");*/
          rc = cov_inv(Info, Znc, Inv);
          if (rc) {
            Rprintf("ERROR computing inverse of information matrix\n");
            error("ERROR");
          }
          /*if (DEBUG) Rprintf("Compute W_y\n");*/
          get_Wy(Y, lxx, ret_p, N, M,W,w_y);
          /*if (DEBUG) Rprintf("Compute delta\n");*/

          get_delta(Inv, Nparm, tXXZ, N, M, w_y, ret_delta);

          /*if (DEBUG) print_dVec(ret_delta, Nparm, "delta");*/

          /* Check for non-finite values */
            if (!all_finite(ret_delta, Nparm)) {
              Rprintf("ERROR: Weighted Least Square Algorithm not converging, parameters have non-finite values\n");
              error("ERROR");
            }

          /*if (DEBUG) Rprintf("Check Weighted Least Square stopping criteria\n");*/
          rerror = checkStop(ret_delta, delta0, Nparm);
          if (rerror <= tolMStep) {
            conv = 1;
            break;
          }

          /* Update delta0 */
            /*if (DEBUG) Rprintf("Update parameters\n");*/
          for (i=0; i<Nparm; i++) delta0[i] = ret_delta[i];
        }


        free(delta0);

        if (!conv) {
          Rprintf("ERROR: algorithm did not converge\n");
          error("ERROR");
        }
        *ret_rc = 0;


  } /* END: Mvpoly */


  static void EstepFitting(double *missing_vec, double **missing_Mat,
                           double *Y, double **X,double*beta,
                           int missing_number,int M,int N,int NCOV){
    double sum=0.0;
    int temp;
    int i;
    for(int t=0;t<missing_number;t++){
      i = missing_vec[t]-1;
      sum = 0.0;
      for(int j=0;j<M;j++){
        if(missing_Mat[t][j]==1){
          temp = j*N+i;
          Y[temp] = 0.0;
          for(int k=0;k<NCOV;k++){
            Y[temp] +=  X[i][k]*beta[j*NCOV+k]; /* get X[i,]*beta[,j] if j subtype is potentail true */
          }
          Y[temp] = exp(Y[temp]);
          sum += Y[temp];
        }
      }
      for(int j=0;j<M;j++){
        if(missing_Mat[t][j]==1){
          Y[j*N+i] = Y[j*N+i]/sum;
        }
      }
    }

  }

static void Free_Mem(double * XX,double **tXXZ,int Nparm,double**X,int N,
                     double *delta0,double**Z,int M,
                     double** XmWXm,int Znr,int Znc,
                     double *lxx,double **Info_obs, double* ret_info,
                     double **Inv, double *w_y, double *W,double *beta,
                     double **missing_Mat,int DEBUG, double **Info,double **Info_obs_inv
                     , int missing_number){
  /*if (DEBUG) Rprintf("Free memory\n");
  if (DEBUG) Rprintf("Free XX\n");*/
  free(XX);
  /*if (DEBUG) Rprintf("Free tXXZ\n");*/
  matrix_free((void **)tXXZ, Nparm);

  /*if (DEBUG) Rprintf("Free X\n");*/
  matrix_free((void **)X, N);
  /*if (DEBUG) Rprintf("Free delta0\n");*/
  free(delta0);
  /*if (DEBUG) Rprintf("Free Z\n");*/
  matrix_free((void **)Z, Znr);
  /*if (DEBUG) Rprintf("Free XmWXm\n");*/
  matrix_free((void**) XmWXm, Znr);
  /*if (DEBUG) Rprintf("Free lxx\n");*/
  free(lxx);
  /*if (DEBUG) Rprintf("Begin fill Mat\n");*/
  fill_SysMat(Info_obs,ret_info,Nparm);
  /*if (DEBUG) Rprintf("Free Info_obs\n");*/
  matrix_free((void **)Info_obs, Znc);
  /*if (DEBUG) Rprintf("Free Inv\n");*/
  matrix_free((void **)Inv, Znc);
  /*if (DEBUG) Rprintf("Free W\n");*/
  free(W);
  /*if (DEBUG) Rprintf("Free W_y\n");*/
  free(w_y);
  /*if (DEBUG) Rprintf("Free beta\n");*/
  free(beta);
  /*if (DEBUG) Rprintf("Free missing_Mat\n");*/
  matrix_free((void**)missing_Mat,missing_number);
  /*if (DEBUG) Rprintf("Free Info\n");*/
  matrix_free((void **)Info,Znc);
  /*if (DEBUG) Rprintf("Free Info_obs_inv\n");*/
  matrix_free((void **)Info_obs_inv,Znc);

}

void EMStepScoreTestSupportMixedModel(deltai, pNparm, Y, Xvec, ZallVec,Zallnr,Zallnc, pN, pM, pNcov, pNiter, ptol,ptolMaxstep,
                        pDEBUG, ret_rc, ret_delta,ret_info,ret_p,missing_vec,
                        missing_Mat_vec,pmissing_number,ret_Inv_info_vec,YminusP,W_obs, WXZ_vec,WX_vec)
double *deltai, *Y, *Xvec, *ptol, *ret_delta,*ret_info,*ret_p,*ZallVec,*missing_vec,
*missing_Mat_vec,*ptolMaxstep, *ret_Inv_info_vec, *YminusP, *W_obs;
int *pNparm, *pN, *pM, *pNcov, *pNiter, *ret_rc, *pDEBUG,*Zallnr,*Zallnc,*pmissing_number;

{
  int i, Niter, M, N, Ncov0, Ncov, iter, Znr, Znc, NM, rc, conv=0;
  int Nparm, DEBUG;
  double tol, **X, **Z_design, *delta0, **Z, rerror, **XmWXm;
  double *w_y, **Inv, **Info,*lxx, **tXXZ;
  double *beta;
  double **missing_Mat;
  int missing_number;
  double tolMaxstep;
  double *XX;
  double **Info_obs;
  double *W;
  double **Info_obs_inv;

  *ret_rc = 1;
  DEBUG   = *pDEBUG;
  Nparm   = *pNparm;
  Niter   = *pNiter;
  N       = *pN;
  M       = *pM;
  tol     = *ptol;
  tolMaxstep = *ptolMaxstep;
  Ncov0   = *pNcov;
  Ncov    = Ncov0 + 1;  /* Allow for intercept */
    Znr     = *Zallnr;
    Znc     = *Zallnc;
    NM      = N*M;
    missing_number = *pmissing_number;
    /*printf("Ncat\tNcatp1\tNcov0\tNcov\tZnr\Znc\n");
    printf("%i\t%i\t%i\t%i\t%i\t%i\n",Ncat,Ncatp1,Ncov0,Ncov,Znr,Znc);*/
      if (Nparm != Znc) error("Nparm != Znc\n");

    /* Allocate memory for matrix of covariates and Z map matrix */
      /*if (DEBUG) Rprintf("Allocate memory\n");*/


    tXXZ     = dMat_alloc(Nparm, NM, 0, 0.0);

    X        = dMat_alloc(N, Ncov, 0, 0.0); /* Xvec doesn't have intercept*/
    delta0   = dVec_alloc(Nparm, 0, 0.0);
    Z        = dMat_alloc(Znr, Znc, 0, 0.0); /* Initialize to 0 */
    lxx      = dVec_alloc(NM, 0, 0.0);
    XmWXm    = dMat_alloc(Znr,Znr,0,0.0);
    Info     = dMat_alloc(Nparm, Nparm, 0, 0.0);
    Info_obs = dMat_alloc(Nparm,Nparm,0,0.0);
    Inv      = dMat_alloc(Nparm, Nparm, 0, 0.0);
    w_y      = dVec_alloc(NM, 0, 0.0);
    W      = dVec_alloc(NM*M, 0, 0.0);
    beta  = dVec_alloc(Znr, 0, 0.0);
    missing_Mat = dMat_alloc(missing_number,M,0,0.0);
    XX = dVec_alloc((N*Ncov*Ncov),0,0.0);
    Info_obs_inv = dMat_alloc(Nparm,Nparm,0,0.0);


    /* Copy initial estimates to delta0 */
    /*if (DEBUG) Rprintf("Copy data\n");*/
    for (i=0; i<Nparm; i++) delta0[i] = deltai[i];
    fillMat(Xvec, N, Ncov0, 1, X);
    fillMat(missing_Mat_vec,missing_number,M,0,missing_Mat);
    /* fillMat(Zvec, M, Ncatp1, 0, Z_design);*/
    /*if (DEBUG) Rprintf("Get the matrix Z\n");*/
    fillMat(ZallVec,Znr,Znc,0,Z);
    /* Get the matrix t(XXZ) */
    /*if (DEBUG) Rprintf("Get the matrix t(XXZ)\n");*/
    get_tXXZ(X, N, Ncov, M, Znc, Z, tXXZ);
    /*if(DEBUG) Rprintf("Get beta\n");
    if(DEBUG) Rprintf("Get Matrix XX\n");*/
    get_XX(X,Ncov,N,XX);
    /* Compute beta = Z*delta */
    X_y(Z, Znr, Znc, delta0, beta);
    /*First EM Step */
    /* since first Estep is outside of C function,add first M step here */
    Mvpoly(deltai,Nparm, Y,
    X, Z,Znr,Znc,
    N, M,Ncov, Niter, tolMaxstep,
    DEBUG,ret_rc,ret_delta,Info,
    ret_p,lxx,W,beta,
    w_y,XmWXm,Inv,tXXZ,XX);
    /* Check the convergence for first EM round */
    if (!all_finite(ret_delta, Nparm)) {
    Rprintf("ERROR: EM algorithm not converging, parameters have non-finite values\n");
    error("ERROR");
    }

    /*if (DEBUG) Rprintf("Check EM algorithm stopping criteria\n");*/
    rerror = checkStop(ret_delta, delta0, Nparm);
    /* Check the convergence for first round EM */
    /* If converged, then we go out of the EM iteration and return result*/
    /* If not we keep iterating */
    if (rerror <= tol) {
    conv = 1;
    /*if (DEBUG) Rprintf("Get Observed Weighted Matrix\n");*/
    Get_ObservedW(W,W_obs ,N, M, Y);
    /*if (DEBUG) Rprintf("Get Observed Information Matrix Matrix\n");*/
    Get_ObservedInfo(M,N,Info_obs, DEBUG,
    XX,X,Ncov,Znr,Znc,Z,W_obs);
    /*if (DEBUG) Rprintf("Get Observed Information Matrix Matrix Inverse\n");*/
    cov_inv(Info_obs,Znc,Info_obs_inv);
    /*if (DEBUG) Rprintf("Get Y-P\n");*/
    VectorMinus(Y,ret_p,NM,YminusP);
    fill_SysMat(Info_obs_inv,ret_Inv_info_vec,Znc);


    Free_Mem(XX,tXXZ,Nparm,X,N,
    delta0,Z,M,
    XmWXm, Znr,Znc,
    lxx,Info_obs,ret_info,
    Inv, w_y, W,beta,
    missing_Mat,DEBUG,Info,Info_obs_inv,missing_number);

    if (!conv) {
    Rprintf("ERROR: algorithm did not converge\n");
    error("ERROR");
    }
    *ret_rc = 0;

    return;
    }else{

    /* Update the parameters for the first E step */
    /*if (DEBUG) Rprintf("Update parameters\n");*/
    for (i=0; i<Nparm; i++) delta0[i] = ret_delta[i];



    /*Begins the second EM round if not converged*/
    for(iter=2; iter<=Niter; iter++) {
    /*if(DEBUG) Rprintf("EM round: %d\n",iter);*/
    X_y(Z, Znr, Znc, delta0, beta);
    EstepFitting(missing_vec,missing_Mat,Y, X,beta,missing_number, M,N,Ncov);
    Mvpoly(delta0,Nparm, Y,
    X, Z,Znr,Znc,
    N, M,Ncov, Niter, tolMaxstep,
    DEBUG,ret_rc,ret_delta,Info,
    ret_p,lxx,W,beta,
    w_y,XmWXm,Inv,tXXZ,XX);


    /* Check for non-finite values */
    if (!all_finite(ret_delta, Nparm)) {
    Rprintf("ERROR: EM algorithm not converging, parameters have non-finite values\n");
    error("ERROR");
    }

    /*if (DEBUG) Rprintf("Check EM algorithm stopping criteria\n");*/
    rerror = checkStop(ret_delta, delta0, Nparm);
    if (rerror <= tol) {
    conv = 1;
    break;
    }

    /* Update delta0 */
    /*if (DEBUG) Rprintf("Update parameters\n");*/
    for (i=0; i<Nparm; i++) delta0[i] = ret_delta[i];
    }
    }

    /*if (DEBUG) Rprintf("Get Observed Weighted Matrix\n");*/
    Get_ObservedW(W,W_obs ,N, M, Y);
    /*if (DEBUG) Rprintf("Get Observed Information Matrix Matrix\n");*/
    Get_ObservedInfo(M,N,Info_obs, DEBUG,
    XX,X,Ncov,Znr,Znc,Z,W_obs);
    /*if (DEBUG) Rprintf("Get Observed Information Matrix Matrix Inverse\n");*/
    cov_inv(Info_obs,Znc,Info_obs_inv);
    /*if (DEBUG) Rprintf("Get Y-P\n");*/
    VectorMinus(Y,ret_p,NM,YminusP);
    fill_SysMat(Info_obs_inv,ret_Inv_info_vec,Znc);



    Free_Mem(XX,tXXZ,Nparm,X,N,
    delta0,Z,M,
    XmWXm, Znr,Znc,
    lxx,Info_obs,ret_info,
    Inv, w_y, W,beta,
    missing_Mat,DEBUG,Info,Info_obs_inv,missing_number);

    /*if (DEBUG) Rprintf("Mem_free Finished\n");*/
    if (!conv) {
    Rprintf("ERROR: algorithm did not converge\n");
    error("ERROR");
    }
    *ret_rc = 0;

    return;
} /* END: Mvpoly */

