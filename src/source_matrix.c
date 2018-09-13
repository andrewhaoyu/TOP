#include "./source.h"

void print_dVec(vec, n, name)
double *vec;
int n;
char name[10];
{
  int i;
  Rprintf("%s \n", name);
  for (i=0; i<n; i++) {
    Rprintf(" %g ", vec[i]);
  }
  Rprintf("\n \n");
}
void print_iVec(vec, n, name)
int *vec;
int n;
char name[10];
{
  int i;
  Rprintf("%s \n", name);
  for (i=0; i<n; i++) {
    Rprintf(" %d ", vec[i]);
  }
  Rprintf("\n \n");
}

void print_dMat(mat, nr, nc, name)
double **mat;
int nr, nc;
char name[10];
{
  int i, j;
  Rprintf("%s \n", name);
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) Rprintf(" %g ", mat[i][j]);
    Rprintf("\n");
  }
  Rprintf("\n \n");
}


/* Function to allocate memory for a double vector */
double * dVec_alloc(n, initFlag, initVal)
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
double ** dMat_alloc(nrow, ncol, initFlag, initVal)
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

/* Function to allocate memory for a integer vector */
int * iVec_alloc(n, initFlag, initVal)
int n, initFlag, initVal;
{
  int i, *ret, *p;

  ret = (int *) malloc(n*sizeof(int));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }

  return(ret);

} /* END: iVec_alloc */

  /* Function to allocate a integer matrix */
int ** iMat_alloc(nrow, ncol, initFlag, initVal)
int nrow, ncol, initFlag, initVal;
{
  int i, **mat, **ptr;

  mat = (int **) malloc(nrow*sizeof(int *));
  CHECK_MEM(mat);
  for (i=0, ptr=mat; i<nrow; i++, ptr++) *ptr = iVec_alloc(ncol, initFlag, initVal);

  return(mat);

} /* END: iMat_alloc */

/* Function to free a matrix */
void matrix_free(x, n)
void **x;
int n;
{
  int i;
  for (i=0; i<n; i++) {
    if (x[i]) free(x[i]);
  }
  free(x);

} /* END: matrix_free */

void copy_dVec(v1, v2, n)
double *v1, *v2;
int n;
{
  int i;
  double *p1, *p2;
  
  for (i=0, p1=v1, p2=v2; i<n; i++, p1++, p2++) *p1 = *p2;

} /* END: copy_dVec */

/* Function to fill in a matrix from a vector (by column) */
void fillMat(vec, nr, nc, addInt, out)
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

/* Function to fill in a matrix from a vector (by column) */
void filliMat(vec, nr, nc, addInt, out)
int *vec, **out;
int nr, nc, addInt;
{
  int i, j, col=0, ii;

  if (addInt) {
    /* Intercept for first column */
      for (i=0; i<nr; i++) out[i][0] = 1;
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

} /* END: filliMat */

/* fill the info matrix to the result*/
/* the Info matrix is sysmetric */
void fill_SysMat(Mat,Vec,Nr)
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


} /*end fill_SysMat*/



/* fill the info matrix to the result*/
/* the Info matrix is sysmetric */
void fill_SysMat_to_vec(Mat,Vec,Nr)
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

} /*end fill_SysMat_to_vec*/

/* Multiply to matrices (matrix mult)  */
void matrixMult(m1, m1_nr, m1_nc, m2, m2_nc, ret)
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

/* two matrix minues  */
void matrixminus(double **mat1,double**mat2,int nr,int nc,double **ret){

    int i, j;
    double **pr1, **pr2, *pc1, *pc2, **row, *col;

    for(i=0, pr1=mat1, pr2=mat2, row=ret; i<nr; i++, pr1++, pr2++, row++){
      for(j=0, pc1=*pr1, pc2=*pr2, col=*row; j<nc; j++, pc1++, pc2++, col++){
        *col = *pc1 - *pc2;
      }
    }

}/* END: matrixminus */

/* two vector minus */
void VectorMinus(double *vec1, double*vec2, int N, double *ret){
  int i;
  double *p1, *p2, *p3;

  for(i=0, p1=vec1, p2=vec2, p3=ret; i<N; i++, p1++, p2++, p3++) *p3 = *p1 - *p2;

} /* END: VectorMinus */

/* Function for dot product of two vectors */
double dotProd(v1, v2, n)
double *v1, *v2;
int n;
{
  int i;
  double sum=0.0, *p1, *p2;

  for (i=0, p1=v1, p2=v2; i<n; i++, p1++, p2++) sum += *p1 * *p2;

  return(sum);

} /* END: dotProd */

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
void get_tXXZ(X, Xnr, Xnc, M, Znc, Z, out)
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

/* calculate the vector combination XX*/
/* X is a vector */
/* We save X1*X1T for later calculation of XmWXm */
/* XX is saved as a N*p^2 vector */
/* Use the sysmetric of XX property */
void get_XX_vec(double *X, int N, double *XX){
  int k;
  double *px, *pxx;

  for(k=0, px=X, pxx=XX; k<N; k++, px++, pxx++) *pxx = *px * *px;

} /* get_XX_vec */

/* calculate the vector combination XX*/
/* X = (X1,..Xp) */
/* We save X1*X1T for later calculation of XmWXm */
/* XX is saved as a N*p^2 vector */
/* Use the sysmetric of XX property */
void get_XX(double **X,int Ncov,int N,double *XX){
    int i, j, k, iNNcov, jNNcov, NNcov, iN, jN, iNNcovjN, jNNcoviN;
    double tmp, *p1, *p2, **prow, *vec;

    NNcov = N*Ncov;

    for(i=0; i<Ncov; i++){
      iNNcov = i*NNcov;
      iN     = i*N;
      for(j=0; j<i+1; j++){
        jNNcov   = j*NNcov;
        jN       = j*N;
        iNNcovjN = iNNcov + jN;
        jNNcoviN = jNNcov + iN;
        for(k=0, p1=&XX[iNNcovjN], p2=&XX[jNNcoviN], prow=X; k<N; k++, p1++, p2++, prow++){
          vec = *prow;
          tmp = vec[i]*vec[j];
          *p1 = tmp;
          *p2 = tmp;
        }
      }
    }

} /* END: get_XX */

/* Function to compute Xy for matrix X and vector */
void X_y(X, nr, nc, y, ret)
double **X, *y, *ret;
int nr, nc;
{
  int i, j;
  double sum, *p, *pret, *px, **prow;

  for (i=0, pret=ret, prow=X; i<nr; i++, pret++, prow++) {
    sum  = 0.0;
    for (j=0, p=y, px=*prow; j<nc; j++, p++, px++) sum += *px * *p;
    *pret = sum;
  }

} /* END: X_y */

/*calculate the quadratic form  XmWXm */
/* Since Xm is are sparse, we take the advantage of this property */
/*XmWXm can be decomposed into M*M blcck, each block is p*p matrix*/
/* i j represent the row and column for M*M block */
/* k,l represent the row and column within the p*p matrix*/
/* p*p matrix can be got by X*(W with N*N block)*X */
/* Since W with N*N block is diagnonal matrix */
/* X*(W with N*N block)*X is weigthed qudractic sum, t loop used for that */
/* Use the sysmetric property of XmWXm */
void get_XmWXm(double *XX,double **X,double *W, int M,int N, int p, double **ret){

    int NM=N*M, i, j, v, k, t, l, NMi, NMiNj, kNp, kNplN, ip, jp, ipk, ipl, jpk, jpl;
    double *Wtemp, sum=0.0, *pW2, *pXX;

    for(i=0; i<M; i++){
      NMi = NM*i;
      ip  = i*p;
      for(j=0; j<i+1; j++){
        /*fill in Wtemp for W with N*N block diagnals*/
        NMiNj = NMi + N*j;
        Wtemp = &W[NMiNj];
        jp    = j*p;
        
        for(k=0; k<p; k++){
          ipk = ip + k;
          jpk = jp + k;
          kNp = k*N*p;
          for(l=0; l<k+1; l++){
            ipl    = ip + l;
            jpl    = jp + l;
            sum    = 0.0;
            kNplN  = kNp + l*N;
            for(t=0, pXX=&XX[kNplN], pW2=Wtemp; t<N; t++, pXX++, pW2++) sum += *pXX * *pW2;
            
            ret[ipk][jpl] = sum;
            ret[ipl][jpk] = sum;
            ret[jpk][ipl] = sum;
            ret[jpl][ipk] = sum;
          }
        }
      }
    }

} /* END: get_XmWXm */

/*calculate the quadratic form  XmWXm */
/* Since Xm is are sparse, we take the advantage of this property */
/*XmWXm can be decomposed into M*M blcck, each block is p*p matrix*/
/* i j represent the row and column for M*M block */
/* Since W with N*N block is diagnonal matrix */
/* X*(W with N*N block)*X is weigthed qudractic sum, t loop used for that */
/* Use the sysmetric property of XmWXm */
void get_XmWXm_vec(double *XX,double *W, int M,int N, double **ret){
    int NM = N*M;
    double *Wtemp;
    Wtemp   = dVec_alloc(N,0, 0.0);
    double sum =0.0;

    for(int i=0;i<M;i++){
      for(int j=0;j<(i+1);j++){
        for(int v=0;v<N;v++){   /*fillin Wtemp for W with N*N block diagnals*/
            Wtemp[v] = W[NM*i+N*j+v];
        }

        sum =0.0;
        for(int t=0;t<N;t++){
          sum += XX[t]*Wtemp[t];
        }
        ret[i][j] = sum;

        ret[j][i] = sum;

      }
    }
    free(Wtemp);

} /* END: get_XmWXm_vec */

/* Function for quadractic computation X^tkX */
void QuadXKX(double **X,double ** K, int Xnr,int Xnc, double** ret){
    double sum, Xki, *pK;
    int i, j, k, l;    

    for(i=0; i<Xnc; i++){
      for(j=0; j<i+1; j++){
        /*the ith row and jth column of the ret*/
        /*the ith column X transpose times K times the jth column of the X */
        sum = 0.0;
        /* One vector times K times one Vector*/
        for(k=0; k<Xnr; k++){
          Xki  = X[k][i];
          for(l=0, pK=K[k]; l<Xnr; l++, pK++) sum += Xki * *pK * X[l][j];    
        }
        ret[i][j] = sum;
        /* ret is sysmetric */
        ret[j][i] = sum;
      }
    }

} /* END: QuadXKX */

/* Function for quadractic computation X^tkX */
void QuadXtKX(double **X,double ** K, int Xnr,int Xnc, double** ret){
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

} /* END: QuadXtKX */

/* Function for quadractic computation XkXt */
void QuadXKXt(double **X,double ** K, int Xnr,int Xnc, double** ret){
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

} /* END: QuadXKXt */

/**transform a matrix **/
void transform_x(double **x, int nr, int nc, double** ret){
    for(int i=0;i <nr; i++){
      for(int j=0;j <nc; j++){
        ret[j][i] = x[i][j];
      }
    }

} /* END: transform_x */




/* fill the info matrix to the result*/
void fill_Info(Info,ret_info,Nparm)
double **Info,*ret_info;
int Nparm;
{
  int i,j;
  for(j=0; j<Nparm; j++){
    for(i=0; i< Nparm; i++){
      ret_info[i*Nparm+j] = Info[i][j];
    }
  }


} /*end fill_Info*/



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

/* Function to compute the inverse of a covariance matrix */
int cov_inv(cov, n, inv)
double **cov;
int n;
double **inv; /* Returned inverse */
{
  double cc, a, b, d;
  int ret;

  switch (n) {
    case 0:
      Rprintf("\nERROR: dimension of covariance matrix is 0\n");
      error("1");
    case 1:
      a = cov[0][0];
      if (fabs(a) < ALMOST_ZERO) return(ERROR_SINGULAR_MATRIX);
      inv[0][0] = 1.0/a;
      break;
      case 2:
        a  = cov[0][0];
        b  = cov[0][1];
        d  = cov[1][1];
        cc = a*d - b*b;
        if (fabs(cc) < ALMOST_ZERO) return(ERROR_SINGULAR_MATRIX);
        cc = 1.0/cc;
        inv[0][0] = d*cc;
        inv[0][1] = -b*cc;
        inv[1][0] = -b*cc;
        inv[1][1] = a*cc;
        break;
        default:
          ret = symPosMatInv(cov, n, inv);
          if (ret) return(ret);
          break;
  } /* END: switch */

    return(0);

} /* END: cov_inv */
