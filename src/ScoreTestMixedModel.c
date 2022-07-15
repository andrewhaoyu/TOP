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

 

  /**transform a matrix **/

  static void transform_x(double **x, int nr, int nc, double** ret){
    for(int i=0;i <nr; i++){
      for(int j=0;j <nc; j++){
        ret[j][i] = x[i][j];
      }
    }
  }

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

 

 



  /* fill the info matrix to the result*/
  /* the Info matrix is sysmetric */
  static void fill_SysMat_to_vec(Mat,Vec,Nr)
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

/*end fill_SysMat_to_vec*/

  

 

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

 



  

/* calculate the vector combination XX*/
  /* X is a vector */
  /* We save X1*X1T for later calculation of XmWXm */
  /* XX is saved as a N*p^2 vector */
  /* Use the sysmetric of XX property */
  static void get_XX_vec(double *X,int N,double *XX){

    for(int k=0;k<N;k++){
      XX[k] = X[k]*X[k];
    }
  }



/*calculate the quadratic form  XmWXm */
  /* Since Xm is are sparse, we take the advantage of this property */
  /*XmWXm can be decomposed into M*M blcck, each block is p*p matrix*/
  /* i j represent the row and column for M*M block */
  /* Since W with N*N block is diagnonal matrix */
  /* X*(W with N*N block)*X is weigthed qudractic sum, t loop used for that */
  /* Use the sysmetric property of XmWXm */

  static void get_XmWXm_vec(double *XX,double *W, int M,int N, double **ret){
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

/* Function for quadractic computation XkXt */
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












/* Function for getting the observed information matrix */
  static void Get_ObservedInfo(int M,int N,double **Info_obs,int DEBUG,
                               double *X,int Ncov,int Znr,int Znc,double **Z,double*W_obs){
    double **XmWobsXm;
    double *XX;
    XmWobsXm    = dMat_alloc(Znr,Znr,0,0.0);
    XX = dVec_alloc(N,0,0.0);
    /*if (DEBUG) Rprintf("get XX\n");*/
    get_XX_vec(X,N,XX);
    /*if (DEBUG) Rprintf("Compute XmWobsXm matrix\n");*/
    get_XmWXm_vec(XX,W_obs,M, N,XmWobsXm);
    /*if (DEBUG) Rprintf("Compute information matrix\n");*/
    /*get_Info(X, N, Ncov, M, Z, Znr, Znc, pxx, Info);*/
      QuadXtKX(Z,XmWobsXm, Znr,Znc, Info_obs);
    matrix_free((void **)XmWobsXm,Znr);
    free(XX);

  }








/* this function is used to compute out the XtYminusP*/
  /* Xt is the transpose of knroneker product of X and diag M */
  /* XtYminusP is a vector with M length*/
  static void get_XtYminusP(double *x_intere,double *YminusP, double *XtYminusP,
                            int N, int M){
    for(int i=0;i< M;i++){
      XtYminusP[i] =0;
      for(int k=0;k<N;k++){
        XtYminusP[i] += x_intere[k]*YminusP[N*i+k];
      }
    }


  }

/* Compute Score */
  /* Score = t(z_intere) %*% XtYminusP */
  static void get_Score(double ** z_intere,
                        double *XtYminusP,
                        int z_intere_nr,
                        int z_intere_nc,
                        double *score){
    double sum = 0.0;
    for(int i=0;i<z_intere_nc;i++){
      sum = 0.0;
      for(int j=0;j<z_intere_nr;j++){
        sum += z_intere[j][i] * XtYminusP[j];
      }
      score[i] = sum;
    }
  }



/*function to compute the t(x_intere_M)%*%W%*%X */
/* x_intere_M is sparse */
/* the output would be a M* (M*Ncov_c0) matrix */
/* i represent the ith row */
/* j represent the jth block in the 1*p block */
 /*in total there are M*M block*/
/* k represent kth column with the ith row, jth 1*p block */
static void get_tx_intere_W_X(double* x_intere,double *W,double **X,
                              int M,int N, int Ncov_c0,double **ret){
  double sum =0.0;
  double* Wtemp;
  int NM = N*M;
  int p = Ncov_c0;

  Wtemp = dVec_alloc(N,0,0.0);

  for(int i=0;i<M; i++){
    for(int j=0;j<M;j++){

      for(int l=0;l <N; l++){
        Wtemp[l] = W[NM*i+N*j+l];
      }


      for(int k=0;k<p;k++){
        sum =0.0;
        for(int l=0;l<N;l++){
          sum += x_intere[l]*Wtemp[l]*X[l][k];
        }
        ret[i][j*p+k] = sum;
      }
      }
  }

  free(Wtemp);
  }

static void get_tz_intere_tx_intere_W_X_zc(double** tz_intere,
                                           double** tx_intere_W_X,
                                           double **zc,
                                           double **ret,
                                           int z_intere_nc,
                                           int z_intere_nr,
                                           int zc_nr,
                                           int zc_nc){
 double ** tz_intere_tx_intere_W_X;
  tz_intere_tx_intere_W_X = dMat_alloc(z_intere_nc,zc_nr,0,0.0);

  matrixMult(tz_intere, z_intere_nc, z_intere_nr, tx_intere_W_X,
             zc_nr, tz_intere_tx_intere_W_X);

  matrixMult(tz_intere_tx_intere_W_X,z_intere_nc,
             zc_nr,zc,zc_nc,ret);


  matrix_free((void**)tz_intere_tx_intere_W_X,z_intere_nc);



}



void ScoreTestMixedModel( double *x_intere ,
                double *z_intere_vec,
                double *inv_info_vec,
                double *YminusP,
                double *W_obs,
                double *score,
                double *efficient_info_vec,
                int *pzc_nr,
                int *pzc_nc,
                int *pz_intere_nr,
                int *pz_intere_nc,
                int *pnparm_intere,
                int *pM,
                int *pN,
                int *pDEBUG,
                double *info_complete_vec,
                double *info_lost_vec,
                double *x_vec,
                double *zc_vec)

{
  int zc_nr = *pzc_nr;
  int zc_nc = *pzc_nc;
  int z_intere_nr = *pz_intere_nr;
  int z_intere_nc = *pz_intere_nc;
  int M = *pM;
  int N = *pN;
  int DEBUG = *pDEBUG;
  int Ncov = 1;
  int Ncov_c0 = (zc_nr/M);
  int Ncov_c = (Ncov_c0-1);

  double **z_intere;
  double **Inv_info;
  double **efficient_info;
  double * XtYminusP;



  double **info_lost;
  double **info_complete;
  double **X;
  double **tx_intere_W_X;
  double **zc;
  double **tz_intere;
  double **tz_intere_tx_intere_W_X_zc;
  /*if (DEBUG) Rprintf("Allocate memory\n");
  if (DEBUG) Rprintf("Allocate z_intere\n");*/
  z_intere = dMat_alloc(z_intere_nr,z_intere_nc,0,0.0);
  /*if (DEBUG) Rprintf("Finish z_intere\n");*/

  zc = dMat_alloc(zc_nr,zc_nc,0,0.0);

  tx_intere_W_X = dMat_alloc(M,M*Ncov_c0,0,0.0);

  X        = dMat_alloc(N, Ncov_c0,0, 0.0);

  tz_intere = dMat_alloc(z_intere_nc,z_intere_nr,0,0.0);
  tz_intere_tx_intere_W_X_zc = dMat_alloc(z_intere_nc,zc_nc,0,0.0);

  fillMat(x_vec, N, Ncov_c, 1, X);
  fillMat(zc_vec,zc_nr,zc_nc,0,zc);

  Inv_info = dMat_alloc(zc_nc,zc_nc,0,0.0);
  /*if (DEBUG) Rprintf("Finish inv_info\n");
  if (DEBUG) Rprintf("Allocate efficient_info\n");*/
  /*efficient_info = dMat_alloc(nparm_intere,nparm_intere,0,0.0);*/
  /*if (DEBUG) Rprintf("Allocate XtYminusP\n");*/
  XtYminusP = dVec_alloc(M,0,0.0);




  /*if (DEBUG) Rprintf("Allocate info_lost\n");*/
  info_lost = dMat_alloc(z_intere_nc,z_intere_nc,0,0.0);
  /*if (DEBUG) Rprintf("Allocate info_complete\n");*/
  info_complete = dMat_alloc(z_intere_nc,z_intere_nc,0,0.0);
  /*if (DEBUG) Rprintf("Allocate efficient_info\n");*/
  efficient_info = dMat_alloc(z_intere_nc,z_intere_nc,0,0.0);
  /*if (DEBUG) Rprintf("Fill in Matrix\n");*/
  fillMat(z_intere_vec,z_intere_nr,z_intere_nc,0,z_intere);
  fillMat(inv_info_vec,zc_nc,zc_nc,0,Inv_info);
  /*if (DEBUG) Rprintf("Get Xt(Y-P)\n");*/
  get_XtYminusP(x_intere,YminusP, XtYminusP,N,  M);
  /*if (DEBUG) Rprintf("Get Score\n");*/
  get_Score(z_intere, XtYminusP, z_intere_nr, z_intere_nc,score);

  transform_x(z_intere,z_intere_nr,z_intere_nc,tz_intere);

  get_tx_intere_W_X(x_intere,W_obs,X,
                     M,N,Ncov_c0,tx_intere_W_X);
  get_tz_intere_tx_intere_W_X_zc(tz_intere,
                                 tx_intere_W_X,
                                 zc,
                                 tz_intere_tx_intere_W_X_zc,
                                 z_intere_nc,
                                 z_intere_nr,
                                 zc_nr,
                                 zc_nc);



  /*if (DEBUG) Rprintf("Get info_lost\n");*/
  QuadXKXt(tz_intere_tx_intere_W_X_zc,Inv_info,z_intere_nc,zc_nc,info_lost);
  /*if (DEBUG) Rprintf("Get info_complete\n");*/
  Get_ObservedInfo(M,N,info_complete, DEBUG,
                   x_intere,Ncov,z_intere_nr,z_intere_nc,z_intere,W_obs);
  /*if (DEBUG) Rprintf("Get efficient information matrix\n");*/
  MatrixMinus(info_complete,info_lost, z_intere_nc,z_intere_nc, efficient_info);
  /*if (DEBUG) Rprintf("fill in efficient information matrix\n");*/
  fill_SysMat_to_vec(efficient_info,efficient_info_vec,z_intere_nc);
  /*if (DEBUG) Rprintf("fill in complete information matrix\n");*/
  fill_SysMat_to_vec(info_complete,info_complete_vec,z_intere_nc);
  /*if (DEBUG) Rprintf("fill in lost information matrix\n");*/
  fill_SysMat_to_vec(info_lost,info_lost_vec,z_intere_nc);



  /*if (DEBUG) Rprintf("Free Memory\n");*/
  matrix_free((void**)z_intere,z_intere_nr);
  matrix_free((void**)Inv_info,zc_nc);
  matrix_free((void**)efficient_info,z_intere_nc);
  free(XtYminusP);
  matrix_free((void**)info_lost,z_intere_nc);
  matrix_free((void**)info_complete,z_intere_nc);
  matrix_free((void**)X,N);
  matrix_free((void**)tx_intere_W_X,M);
  matrix_free((void**)zc ,zc_nr);
  matrix_free((void**)tz_intere ,z_intere_nc);
  matrix_free((void**)tz_intere_tx_intere_W_X_zc,z_intere_nc);


} /* END: ScoreTest */

