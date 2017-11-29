#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <math.h>
#define PI 3.14159265359

double weightfunc(double x[3]){
  double b = -3.0/2.0;
  return pow(PI,b)*exp(-x[0]*x[0]-x[1]*x[1]-x[2]*x[2]);
}

double integrand(double x[3]){
  return x[0]*x[0] + x[0]*x[0]*x[1]*x[1] + x[0]*x[0]*x[1]*x[1]*x[2]*x[2]; 
}

void trial_state(double x_t[3], double x_m[3], double delta, gsl_rng *q){
  double u;
  for (int i = 0; i < 3; i++){
    u = gsl_rng_uniform(q);
    x_t[i] = x_m[i] + delta*(u-0.5);
    
  }
}

void next_state(double x_t[3], double x_m[3], int *count, gsl_rng *q){
  double p_t = weightfunc(x_t);
  double p_m = weightfunc(x_m);
  double xi = gsl_rng_uniform(q);
  if (p_t/p_m > xi){
    for (int i = 0; i < 3; i++){
      x_m[i] = x_t[i];
    }
    *count = *count+1;
  }
}




int main()  {

  double x_t[3];
  double x_m[3];
  double delta = 2.05;
  int count=0;
  double I_sum=0;
  double I;
  
  /* Initial state */
  x_m[0] = 0;
  x_m[1] = 0;
  x_m[2] = 0;

  
  int N = 1e7; // Number of steps
  /* Generate random numbers */
  const gsl_rng_type*T;
  gsl_rng *q;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  q = gsl_rng_alloc(T);
  gsl_rng_set(q, time(NULL));
  /*________________________*/

  /* Throw away states in begining */
  for (int i = 0; i < 100; i++){
    trial_state(x_t, x_m, delta, q);
    next_state(x_t, x_m, &count, q);
  }

    for (int i = 0; i < N; i++){
    trial_state(x_t, x_m, delta, q);
    next_state(x_t, x_m, &count, q);
    I_sum += integrand(x_m);
  }

    I=I_sum/N;
    printf("The stepping percentage is %f%%.\n",(double)count/N);
    printf("The result is I =  %f \n",I);
  return 0; 
}
