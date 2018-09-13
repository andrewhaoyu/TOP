#ifndef __SOURCE_H__
#define __SOURCE_H__


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

/* source_matrix.c */
double * dVec_alloc(int, int, double);
double ** dMat_alloc(int, int, int, double);
int * iVec_alloc(int, int, int);
int ** iMat_alloc(int, int, int, int);
void matrix_free(void **, int);
void copy_dVec(double *, double *, int);
void fillMat(double *, int, int, int, double **);
void filliMat(int *, int, int, int, int **);
void fill_SysMat_to_vec(double **, double *, int);
void fill_SysMat(double **, double *, int);
void get_tXXZ(double **, int, int, int, int, double **, double **);
void get_XX(double **,int,int,double *);
void get_XX_vec(double *, int, double *);
void X_y(double **, int, int, double *, double *);
void get_XmWXm(double *,double **,double *, int,int, int,double **);
void get_XmWXm_vec(double *,double *, int, int, double **);
void matrixMult(double **, int, int, double **, int, double **);
void matrixminus(double **,double**,int,int,double **);
void VectorMinus(double *,double*,int,double *);
double dotProd(double *, double *, int);
int cov_inv(double **, int, double **);
void QuadXKX(double **, double **, int, int, double**);
void QuadXtKX(double **,double **, int, int, double**);
void QuadXKXt(double **,double **, int, int, double**);
void print_dVec(double *, int, char *);
void print_iVec(int *, int, char *);
void print_dMat(double **, int, int, char *);
void transform_x(double **, int, int, double **);
void fill_Info(double **, double *, int);

/* source_MVpoly.c */
void Mvpoly(double *, int, double*, double **, double **, int, int, int, int, int, int, double,\
            int, int*, double*, double**, double *, double *, double *, double *,\
            double*, double **, double **, double **, double *);

/* source.c */
void Weighted_W(double *, double *,int,int);
double checkStop(double *, double *, int);
double LogLikelihood(double *, double *, int, int);
void myconvert(int *, double *, double *);
int all_finite(double *, int);







#endif
