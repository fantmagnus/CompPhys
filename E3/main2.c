#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <math.h>
#define N 10000
#define PI 3.14159265359

int main()  {
  //declarate a rng
  double u;
  double sum_I;
  double x;
  const gsl_rng_type*T;
  gsl_rng *q;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  q = gsl_rng_alloc(T);
  gsl_rng_set(q, time(NULL));

  FILE *file1;

  double *I = malloc(N * sizeof(double));
  double *eta = malloc(N * sizeof(double));
  sum_I = 0; 

  for(int i = 0; i < N; i++ ) {
    u = gsl_rng_uniform(q);
    /* inverse to CDF for sin(pi*x) */
    x = 2*asin(sqrt(u))/PI;
    I[i] = x*(1-x)/(0.5*PI*sin(PI*x));
    sum_I  += I[i];
    eta[i] = x;
  }
  printf("%.5f\n", sum_I/N);

  file1 = fopen("generated_points.dat","w");
  for (int i = 0; i < N; i++) {
    fprintf(file1,"%.5f\n",eta[i]);
  }
  fclose(file1);

  free(I);
  free(eta);
  return 0; 
}


