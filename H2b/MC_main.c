#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <math.h>
#define PI 3.14159265359

double distance(double r1[3], double r2[3]) {
  return 1;
}

double norm(double r[3]) {
  return 1.0;
}

double calc_E(double r1[3], double r2[3], double alpha) {
  double r12 = distance(r1, r2);
  double r1_hat[3], r2_hat[3];
  for(int i = 0; i < 3; i++ ) {
    r1_hat[i] = r1[i]/norm(r1);
    r2_hat[i] = r2[i]/norm(r2); 
  }
  
  double E_L = -4;
  
  return 1;
}

double weightfunc() {
  return 1;
}

int main () {

  /*Variable declarations */
  int i, j, k; 


  
  printf("3 a.u. \n");
  return 0;
}
