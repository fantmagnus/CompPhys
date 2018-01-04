#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <math.h>

/*Potential for harmonix oscillator */
double potential(double x) {
  double potential = x*x/2;
  //printf("%f\n", potential);
  return potential; 
}

/*Weight function to determine if walkers are to be added or removed */
double weight_func(double x, double E_T, double timestep) {
  double p = potential(x);
  double W = exp(-timestep*(p - E_T));
  printf("%f\n", W);
  return W; 
}

/* Branching function */
void make_branch(double walkers[], int *N, int index, double W,  gsl_rng *q) {
  double x = walkers[index];
  double u = gsl_rng_uniform(q);
  int m = (W + u);
  if (m == 0){
    walkers[index] = walkers[*N-1];
    *N = *N - 1;
    double* new_walkers = realloc(walkers, *N*sizeof(double));
    walkers = new_walkers;
  }
  else if (m > 1) {
    *N  = *N + m-1;
    double* new_walkers = realloc(walkers, *N*sizeof(double));
    for(int i = 0; i < m; i++) {
      new_walkers[*N-i-1] = x;
    }
    walkers = new_walkers;
  }

  
}


int main() {
  /*Initializing variables */
  int i,j,k;
  int N_0 = 300;
  int N = N_0;
  int cur_N;
  double timestep = 0.1;
  int nbr_of_timesteps = 10000;
  double E_T = 2;
  double E_0 = 1/2;
  double W;
  /* Initialize random number */
  const gsl_rng_type*o;
  gsl_rng *q;
  gsl_rng_env_setup();
  o = gsl_rng_default;
  q = gsl_rng_alloc(o);
  gsl_rng_set(q, time(NULL));
  /*---------------------------*/
  double* walkers = malloc(N_0 * sizeof(double));
  /* Initialize walkers */
  for (i = 0; i < N_0; i++){
    double u = gsl_rng_uniform(q);
    u = u; // Set range of random intial position
    if (i > N_0/2) { // Put every second walker at negative position
      u = -u;
    }
    walkers[i] = u;
    //printf("%f\n", u);
  }
  /*--------------------*/

  cur_N = N_0; 
  for (i = 0; i < cur_N; i++){
    W = weight_func(walkers[i], E_T, timestep);
    make_branch(walkers, &N, i, W, q);
  }

  printf("%d\n",N);
  


    
}
