#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <math.h>

/*Potential for harmonix oscillator */
double potential(double x) {
  double potential = x*x/2;
  return potential; 
}

double weight_func(double x, double E_T, double timestep) {
  double p = potential(x);
  double W = exp(-timestep*(p - E_T));
  return W; 
}

int main() {
  printf("Spindelhora \n"); 
}
