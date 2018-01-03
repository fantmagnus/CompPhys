#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <math.h>

/*Potential for harmonix oscillator */
double potential(double x) {
  double potential = x*x/2;
  return potential; 
}

/*Weight function to determine if walkers are to be added or removed */
double weight_func(double x, double E_T, double timestep) {
  double p = potential(x);
  double W = exp(-timestep*(p - E_T));
  return W; 
}

int main() {
  /*Initializing variables */
  int N_0 = 300;
  int N = N_0;
  int timestep = 0.1;
  int nbr_of_timesteps = 10000;
  double E_T = 2;
  double E_0 = 1/2;

  const gsl_rng_type*o;
  gsl_rng *q;
  gsl_rng_env_setup();
  o = gsl_rng_degault;
  q = gsl_rng_alloc(o);
  gsl_rng_set(q, time(NULL));
  printf("");
    
}
