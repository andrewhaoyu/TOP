#include "./source.h"

/* this function is used to compute out the XtYminusP*/
  /* Xt is the transpose of knroneker product of X and diag M */
  /* XtYminusP is a vector with M length*/
static void get_XtYminusP(double *x_intere,double *YminusP, double *XtYminusP,
                            int N, int M){
    int NM = N*M;
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


/*function to compute the t(x_intere_M)%*%WXZ */
  /* x_intere_M is sparse */
  /* the output would be a M*zc_nc matrix */
  /* i represent the ith row */
  /* j represent the jth column */
  /* k represent kth element for a N length vector in WXZ */
static void get_tx_intereWXZ(double* x_intere,double **WXZ,int N,int M,int zc_nc,
                               double ** tx_intereWXZ){
    double sum =0.0;
    for(int i=0;i<M;i++){
      for(int j=0;j<zc_nc;j++){
        sum = 0.0;
        for(int k=0;k<N;k++){
          sum += WXZ[i*N+k][j]*x_intere[k];
        }
        tx_intereWXZ[i][j] = sum;
      }

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

/* Function for getting the observed information matrix */
static void Get_ObservedInfo(int M,int N,double **Info_obs,int DEBUG,
                               double *X,int Ncov,int Znr,int Znc,double **Z,double*W_obs){
    double **XmWobsXm;
    double *XX;

    XmWobsXm = dMat_alloc(Znr,Znr,0,0.0);
    XX       = dVec_alloc(N,0,0.0);
    if (DEBUG) Rprintf("get XX\n");
    get_XX_vec(X,N,XX);
    if (DEBUG) Rprintf("Compute XmWobsXm matrix\n");
    get_XmWXm_vec(XX,W_obs,M, N,XmWobsXm);
    QuadXtKX(Z,XmWobsXm, Znr,Znc, Info_obs);
    matrix_free((void **)XmWobsXm,Znr);
    free(XX);

}

void ScoreTest(double *x_intere,
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
                double *zc_vec) {

  int zc_nr        = *pzc_nr;
  int zc_nc        = *pzc_nc;
  int z_intere_nr  = *pz_intere_nr;
  int z_intere_nc  = *pz_intere_nc;
  int nparm_intere = *pnparm_intere;
  int M            = *pM;
  int N            = *pN;
  int DEBUG        = *pDEBUG;
  int Ncov         = 1;
  int Ncov_c0      = zc_nr/M;
  int Ncov_c       = Ncov_c0-1;

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

  if (DEBUG) Rprintf("Allocate memory\n");
  z_intere                   = dMat_alloc(z_intere_nr,z_intere_nc,0,0.0);
  zc                         = dMat_alloc(zc_nr,zc_nc,0,0.0);
  tx_intere_W_X              = dMat_alloc(M,M*Ncov_c0,0,0.0);
  X                          = dMat_alloc(N, Ncov_c0,0, 0.0);
  tz_intere                  = dMat_alloc(z_intere_nc,z_intere_nr,0,0.0);
  tz_intere_tx_intere_W_X_zc = dMat_alloc(z_intere_nc,zc_nc,0,0.0);
  Inv_info                   = dMat_alloc(zc_nc,zc_nc,0,0.0);
  XtYminusP                  = dVec_alloc(M,0,0.0);  
  info_lost                  = dMat_alloc(z_intere_nc,z_intere_nc,0,0.0);
  info_complete              = dMat_alloc(z_intere_nc,z_intere_nc,0,0.0);
  efficient_info             = dMat_alloc(z_intere_nc,z_intere_nc,0,0.0);

  if (DEBUG) Rprintf("Fill in Matrix\n");
  fillMat(x_vec, N, Ncov_c, 1, X);
  fillMat(zc_vec,zc_nr,zc_nc,0,zc);
  fillMat(z_intere_vec,z_intere_nr,z_intere_nc,0,z_intere);
  fillMat(inv_info_vec,zc_nc,zc_nc,0,Inv_info);

  if (DEBUG) Rprintf("Get Xt(Y-P)\n");
  get_XtYminusP(x_intere,YminusP, XtYminusP,N,  M);

  if (DEBUG) Rprintf("Get Score\n");
  get_Score(z_intere, XtYminusP, z_intere_nr, z_intere_nc,score);

  if (DEBUG) Rprintf("transform_x\n");
  transform_x(z_intere,z_intere_nr,z_intere_nc,tz_intere);

  if (DEBUG) Rprintf("get_tx_intere_W_X\n");
  get_tx_intere_W_X(x_intere,W_obs,X,M,N,Ncov_c0,tx_intere_W_X);

  if (DEBUG) Rprintf("get_tz_intere_tx_intere_W_X_zc\n");
  get_tz_intere_tx_intere_W_X_zc(tz_intere, tx_intere_W_X, zc,\
    tz_intere_tx_intere_W_X_zc, z_intere_nc, z_intere_nr, zc_nr, zc_nc);

  if (DEBUG) Rprintf("Get info_lost\n");
  QuadXKXt(tz_intere_tx_intere_W_X_zc,Inv_info,z_intere_nc,zc_nc,info_lost);

  if (DEBUG) Rprintf("Get info_complete\n");
  Get_ObservedInfo(M,N,info_complete, DEBUG,\
                   x_intere,Ncov,z_intere_nr,z_intere_nc,z_intere,W_obs);

  if (DEBUG) Rprintf("Get efficient information matrix\n");
  matrixminus(info_complete,info_lost, z_intere_nc,z_intere_nc, efficient_info);

  if (DEBUG) Rprintf("fill in efficient information matrix\n");
  fill_SysMat_to_vec(efficient_info,efficient_info_vec,z_intere_nc);

  if (DEBUG) Rprintf("fill in complete information matrix\n");
  fill_SysMat_to_vec(info_complete,info_complete_vec,z_intere_nc);

  if (DEBUG) Rprintf("fill in lost information matrix\n");
  fill_SysMat_to_vec(info_lost,info_lost_vec,z_intere_nc);

  if (DEBUG) Rprintf("Free Memory\n");
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

  return;

} /* END: ScoreTest */

