#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <R.h>
#include <R_ext/Memory.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
void convert(int *n,double *p, double *snpvalue){
  for(int i = 0; i<(*n); i++){
    snpvalue[i] = p[(3*i+1)]+p[(3*i+2)]*2;
  }
  return;
}
