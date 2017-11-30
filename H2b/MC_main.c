#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <math.h>
#define PI 3.14159265359

double distance(double r1[3], double r2[3]) {
  double dist = sqrt(pow(r2[0]-r1[0],2)+ pow(r2[1]-r1[1],2)+pow(r2[2]-r1[2],2));
  return dist;
}

double norm(double r[3]) {
  double norm = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
  return sqrt(norm);
}

double calc_E(double r1[3], double r2[3], double alpha) {
  double r12 = distance(r1, r2);
  double r1_hat[3], r2_hat[3], r_hat[3], r[3];
  double scal; 
  for(int i = 0; i < 3; i++ ) {
    r1_hat[i] = r1[i]/norm(r1);
    r2_hat[i] = r2[i]/norm(r2);
    r_hat[i] = r1_hat[i] - r2_hat[i];
    r[i] = r1[i] - r2[i];
    scal = r_hat[i]*r[i]; 
  }
  
  double E_L = -4 + scal/(r12*pow(1+alpha*r12,2)) - 1/(r12*pow(1+alpha*r12,3)) - 1/(4*pow(1+alpha*r12,4)) + 1/r12;
  
  return E_L;
}

double weightfunc(double r1[3], double r2[3], double alpha) {
  double r12 = distance(r1, r2);
  double psi_T = exp(-2*norm(r1))*exp(-2*norm(r2))*exp(0.5*r12/(1+alpha*r12));
  return abs(psi_T) * abs(psi_T);
}

int main () {

  /*Variable declarations */
  int i, j, k; 


  
  printf("3 a.u. \n");
  return 0;
}
