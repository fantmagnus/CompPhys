#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#define N 10000 

int main()  {
  //declarate a rng
  double u;
  double sum_I;
  const gsl_rng_type*T;
  gsl_rng *q;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  q = gsl_rng_alloc(T);
  gsl_rng_set(q, time(NULL));

  double *I = malloc(N * sizeof(double));

  sum_I = 0; 

  for(int i = 0; i < N; i++ ) {
    u = gsl_rng_uniform(q);
    I[i] = u*(1-u);
    sum_I  += I[i]; 
  }
  printf("%.5f\n", sum_I/N);
  return 0; 
}
