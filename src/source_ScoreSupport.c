#include "./source.h"

/* Function for getting the observed weighted matrix  W_com-W_com|mis*/
static void Get_ObservedW_EM(double *W,double *W_obs,int N,int M, double *Y){
    double * W_mis  ;
    double **Info_mis, **XmWmisXm;
    W_mis = dVec_alloc((M*M*N),0,0.0);
    Weighted_W(Y, W_mis, N, M);
    VectorMinus(W,W_mis,N*M*M,W_obs);
    free(W_mis);
}

/* Function for getting the observed information matrix */
static void Get_ObservedInfo_EM(int M,int N,double **Info_obs,int DEBUG,
                               double *XX,double **X,int Ncov,int Znr,int Znc,double **Z,double*W_obs){
    double **XmWobsXm;

    XmWobsXm = dMat_alloc(Znr,Znr,0,0.0);

    if (DEBUG) Rprintf("Compute XmWobsXm matrix\n");
    get_XmWXm(XX,X,W_obs, M,N, Ncov,XmWobsXm);
    if (DEBUG) Rprintf("Compute mis information matrix\n");
    /*get_Info(X, N, Ncov, M, Z, Znr, Znc, pxx, Info);*/
    QuadXtKX(Z,XmWobsXm, Znr,Znc, Info_obs);
    matrix_free((void **)XmWobsXm,Znr);

}

static void EstepFitting(int *missing_vec, int **missing_Mat,
                           double *Y, double **X, double *beta,
                           int missing_number,int M,int N,int NCOV){
    double sum, sum2, *pbeta, *px, dval;
    int temp, i, t, j, k, flag1, *ptri, **prow, *vec, *pvec, jNCOV;
    
    for(t=0, ptri=missing_vec, prow=missing_Mat; t<missing_number; t++, ptri++, prow++){
      i   = *ptri - 1;
      sum = 0.0;
      vec = *prow;
      for(j=0, pvec=vec; j<M; j++, pvec++){
        if(*pvec){
          temp = j*N+i;
          /*Y[temp] = 0.0;*/
          /* get X[i,]*beta[,j] if j subtype is potentail true */
          jNCOV = j*NCOV;
          sum2  = 0.0;
          for(k=0, px=X[i], pbeta=&beta[jNCOV]; k<NCOV; k++, px++, pbeta++) sum2 += *px * *pbeta; 
          dval    = exp(sum2);
          Y[temp] = dval;
          sum    += dval;
        }
      }
      for(j=0, pvec=vec; j<M; j++, pvec++){
        if(*pvec) Y[j*N+i] /= sum;
      }
    }

} /* END: EstepFitting */



void ScoreSupport(deltai, pNparm, Y, Xvec, ZallVec,Zallnr,Zallnc, pN, pM, pNcov, pNiter, ptol,ptolMaxstep,
                        pDEBUG, ret_rc, ret_delta,ret_info,ret_p,missing_vec,
                        missing_Mat_vec,pmissing_number,ret_Inv_info_vec,YminusP,W_obs)
double *deltai, *Y, *Xvec, *ptol, *ret_delta,*ret_info,*ret_p,*ZallVec, *ptolMaxstep, *ret_Inv_info_vec, *YminusP, *W_obs;
int *pNparm, *pN, *pM, *pNcov, *pNiter, *ret_rc, *pDEBUG,*Zallnr,*Zallnc,*pmissing_number;
int *missing_Mat_vec, *missing_vec;

{
  int i, Niter, M, N, Ncov0, Ncov, iter=1, Znr, Znc, NM, rc, conv=0;
  int Nparm, DEBUG, mvconv=0;
  double tol, **X, **Z_design, *delta0, **Z, rerror, **XmWXm;
  double *w_y, **Inv, **Info,*lxx, **tXXZ;
  double *beta;
  int **missing_Mat = NULL;
  int missing_number;
  double tolMaxstep;
  double *XX;
  double **Info_obs;
  double *W;
  double **Info_obs_inv;

  *ret_rc        = 1;
  DEBUG          = *pDEBUG;
  Nparm          = *pNparm;
  Niter          = *pNiter;
  N              = *pN;
  M              = *pM;
  tol            = *ptol;
  tolMaxstep     = *ptolMaxstep;
  Ncov0          = *pNcov;
  Ncov           = Ncov0 + 1;  /* Allow for intercept */
  Znr            = *Zallnr;
  Znc            = *Zallnc;
  NM             = N*M;
  missing_number = *pmissing_number;
  
  if (Nparm != Znc) error("Nparm != Znc\n");

  if (DEBUG) Rprintf("Allocate memory\n");
  tXXZ         = dMat_alloc(Nparm, NM, 0, 0.0);
  X            = dMat_alloc(N, Ncov, 0, 0.0); /* Xvec doesn't have intercept*/
  delta0       = dVec_alloc(Nparm, 0, 0.0);
  Z            = dMat_alloc(Znr, Znc, 0, 0.0); /* Initialize to 0 */
  lxx          = dVec_alloc(NM, 0, 0.0);
  XmWXm        = dMat_alloc(Znr,Znr,0,0.0);
  Info         = dMat_alloc(Nparm, Nparm, 0, 0.0);
  Info_obs     = dMat_alloc(Nparm,Nparm,0,0.0);
  Inv          = dMat_alloc(Nparm, Nparm, 0, 0.0);
  w_y          = dVec_alloc(NM, 0, 0.0);
  W            = dVec_alloc(NM*M, 0, 0.0);
  beta         = dVec_alloc(Znr, 0, 0.0);
  XX           = dVec_alloc((N*Ncov*Ncov),0,0.0);
  Info_obs_inv = dMat_alloc(Nparm,Nparm,0,0.0);
  if (missing_number) missing_Mat  = iMat_alloc(missing_number,M,0,0);

  /* Copy initial estimates to delta0 */
  if (DEBUG) Rprintf("Copy data\n");
  for (i=0; i<Nparm; i++) delta0[i] = deltai[i];
  fillMat(Xvec, N, Ncov0, 1, X);
  if (missing_number) filliMat(missing_Mat_vec,missing_number,M,0,missing_Mat);
  
  if (DEBUG) Rprintf("Get the matrix Z\n");
  fillMat(ZallVec,Znr,Znc,0,Z);

  /* Get the matrix t(XXZ) */
  if (DEBUG) Rprintf("Get the matrix t(XXZ)\n");
  get_tXXZ(X, N, Ncov, M, Znc, Z, tXXZ);

  if(DEBUG) Rprintf("Get Matrix XX\n");
  get_XX(X,Ncov,N,XX);

  /* Compute beta = Z*delta */
  if(DEBUG) Rprintf("Get beta\n");
  X_y(Z, Znr, Znc, delta0, beta);

  /*First EM Step */

  /* Change this later. We need to make sure The E step is consistent with delta0 */
  if (missing_number) EstepFitting(missing_vec,missing_Mat,Y, X,beta,missing_number, M,N,Ncov);

  /* since first Estep is outside of C function,add first M step here */
  Mvpoly(deltai,Nparm, Y, X, Z,Znr,Znc, N, M,Ncov, Niter, tolMaxstep,\
    DEBUG,ret_rc,ret_delta,Info, ret_p,lxx,W,beta, w_y,XmWXm,Inv,tXXZ,XX);

  /* Check the convergence for first EM round */
  if (!all_finite(ret_delta, Nparm)) {
    Rprintf("ERROR: EM algorithm not converging, parameters have non-finite values\n");
  } 

  /* Check the return code from Mvpoly */
  if (*ret_rc) {
    mvconv = 0;
  } else {
    mvconv = 1;
    if (!missing_number) conv = 1;
  }

  if (mvconv && missing_number) {
    if (DEBUG) Rprintf("Check EM algorithm stopping criteria\n");
    rerror = checkStop(ret_delta, delta0, Nparm);

    /* Check the convergence for first round EM */
    /* If converged, then we go out of the EM iteration and return result*/
    /* If not we keep iterating */
    if (rerror <= tol) conv = 1;

    if (!conv) {

      /* Update the parameters for the first E step */
      if (DEBUG) Rprintf("Update parameters\n");
      for (i=0; i<Nparm; i++) delta0[i] = ret_delta[i];

      /*Begins the second EM round if not converged*/
      for(iter=2; iter<=Niter; iter++) {
        if(DEBUG) Rprintf("EM round: %d\n",iter);
        X_y(Z, Znr, Znc, delta0, beta);

        if (missing_number) {
          if(DEBUG) Rprintf("Calling EstepFitting\n");
          EstepFitting(missing_vec,missing_Mat,Y, X,beta,missing_number, M,N,Ncov);
        }

        Mvpoly(delta0,Nparm, Y, X, Z,Znr,Znc,N, M,Ncov, Niter, tolMaxstep,\
             DEBUG,ret_rc,ret_delta,Info,ret_p,lxx,W,beta,w_y,XmWXm,Inv,tXXZ,XX);

        /* Check the return code from Mvpoly */
        if (*ret_rc) break; /* Did not converge */

        /* Check for non-finite values */
        if (!all_finite(ret_delta, Nparm)) {
          Rprintf("ERROR: EM algorithm not converging, parameters have non-finite values\n");
          break;
        }

        if (DEBUG) Rprintf("Check EM algorithm stopping criteria\n");
        rerror = checkStop(ret_delta, delta0, Nparm);
        if (rerror <= tol) {
          conv = 1;
          break;
        }

        /* Update delta0 */
        if (DEBUG) Rprintf("Update parameters\n");
        for (i=0; i<Nparm; i++) delta0[i] = ret_delta[i];
      }
    }
  } /* END: if (mvconv) */

  if (conv) {
      if (DEBUG) Rprintf("Get Observed Weighted Matrix\n");
      Get_ObservedW_EM(W,W_obs ,N, M, Y);
      if (DEBUG) Rprintf("Get Observed Information Matrix Matrix\n");
      Get_ObservedInfo_EM(M,N,Info_obs, DEBUG, XX,X,Ncov,Znr,Znc,Z,W_obs);
      if (DEBUG) Rprintf("Get Observed Information Matrix Matrix Inverse\n");
      cov_inv(Info_obs,Znc,Info_obs_inv);
      if (DEBUG) Rprintf("Get Y-P\n");
      VectorMinus(Y,ret_p,NM,YminusP);
      fill_SysMat(Info_obs_inv,ret_Inv_info_vec,Znc);
  }

  if (DEBUG) Rprintf("Free memory\n");
  matrix_free((void **)tXXZ, Nparm);
  matrix_free((void **)X, N);
  free(delta0);
  matrix_free((void **)Z, Znr);
  free(lxx);
  matrix_free((void**) XmWXm, Znr);
  matrix_free((void **)Info, Nparm);
  matrix_free((void **)Info_obs, Nparm);
  matrix_free((void **)Inv, Nparm);
  free(w_y);
  free(W);
  free(beta);
  if (missing_number) matrix_free((void**)missing_Mat,missing_number);
  free(XX);
  matrix_free((void **)Info_obs_inv, Nparm);

  if (!conv) {
    Rprintf("ERROR: algorithm did not converge\n");
    error("ERROR");
  } else {
    if (1) Rprintf("algorithm converged in %d iterations\n", iter);
  }
  *ret_rc = 0;

  return;

} /* END: ScoreSupport */

