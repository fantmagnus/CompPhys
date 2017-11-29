#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <math.h>
#define PI 3.14159265359

int main()  {
  //declarate a rng
  double u;
  double sum_I;
  double var_g;
  double x;
  int nbr_of_iterations[4];
  int N = 1e4;
  double I_i;
  const gsl_rng_type*T;
  gsl_rng *q;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  q = gsl_rng_alloc(T);
  gsl_rng_set(q, time(NULL));
  nbr_of_iterations[0]=1e1;
  nbr_of_iterations[1]=1e2;
  nbr_of_iterations[2]=1e3;
  nbr_of_iterations[3]=1e4;
  FILE *file1;

  double *eta = malloc(N * sizeof(double));

  sum_I = 0;
  var_g = 0;
  for (int k = 0; k < 4; k++){
    N = nbr_of_iterations[k];
    sum_I = 0;
    var_g = 0;
    for(int i = 0; i < N; i++ ) {
      u = gsl_rng_uniform(q);
      /* inverse to CDF for sin(pi*x) */
      x = 2*asin(sqrt(u))/PI;
      I_i = x*(1-x)/(0.5*PI*sin(PI*x));
      sum_I  += I_i;
      var_g += I_i*I_i/N;
      if (k==3){
	eta[i] = x;
      }
    }
    var_g -= (sum_I/N) * (sum_I/N);
  
    printf("N = %d,\t I_N = %.5f,\t var_g = %.5f\n ", N, sum_I/N, var_g);
  }

  file1 = fopen("generated_points.dat","w");
  for (int i = 0; i < N; i++) {
    fprintf(file1,"%.5f\n",eta[i]);
  }
  fclose(file1);
;
  free(eta);
  return 0; 
}
