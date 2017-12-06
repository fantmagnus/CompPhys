#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <math.h>
#define PI 3.14159265359
#define k_B 1.38064852e23 // m2 kg s-2 K-1

double calc_a(double x, double v, double f, double eta, gsl_rng *q) {
  double omega_0 = f*2*PI;
    double u = gsl_rng_uniform(q);
    double xi = u;
    double a = -eta * v - omega_0*omega_0*x + xi;
    return a;
}
/*Transformation according to Box-Muller */
void gen_gaussian(double *G1, double *G2, gsl_rng *q){
  double U = gsl_rng_uniform(q);
  double V = gsl_rng_uniform(q);
  *G1 = sqrt(-2*log(U))*cos(2*PI*V);
  *G2 = sqrt(-2*log(U))*sin(2*PI*V);  
}

int main()  {
  int i,j;
  double v_T;
  double x_T;
  double v_tilde = 0;
  double a;
  double timestep = 0.1; // in ms
  double f_0 = 3e3; //Frequency of 3 kHz
  double tau = 48.5e-3; // 48.5e-3 for case A, 147.3e-3 for case B
  double eta = 1/tau;
  double c_0 = exp(-eta * timestep);
  double sqrt_c_0 = sqrt(c_0);
  double m = 1;
  double T = 295; //K 
  double v_th = sqrt(k_B*T/m );
  double current_time = 0; 
  int N = 1e4;
  int progress = 0;
  double G1, G2;
  /* Generate random numbers */
  const gsl_rng_type*o;
  gsl_rng *q;
  gsl_rng_env_setup();
  o = gsl_rng_default;
  q = gsl_rng_alloc(o);
  gsl_rng_set(q, time(NULL));
  /*________________________*/


  double *x = malloc(N * sizeof(double));
  double *v = malloc(N * sizeof(double));
  
  FILE* x_file = fopen("x_data.dat","w");
  FILE* v_file = fopen("v_data.dat","w");
  
  x[0] = 0.1;
  x_T = x[0];
  v[0] = 2.0;
  v_T = v[0];
  a = calc_a(x_T, v_T, f_0, eta, q);

  for(i = 1; i < N; i ++) {
    gen_gaussian(&G1, &G2, q);
    /*Verlet velocity algorithm*/
    /*Velocity v(t+dt/2) and displacement q(t+dt) */
    v_tilde = 0.5*a*timestep + sqrt_c_0*v_T + v_th*sqrt(1-c_0)*G1;

    x_T = x_T + v_tilde*timestep;
    x[i] = x_T;
	
    /*Calculate the acceleration for v(t+dt) */
    a = calc_a(x_T, v_tilde, f_0, eta, q);
    /*Velocity v(t+dt) */
    
    v_T = 0.5*sqrt_c_0*a*timestep + sqrt_c_0*v_tilde + v_th*sqrt(1-c_0)*G2;
    v[i] = v_T;

    
    if (i%(N/100) == 0) {
      printf("\rProgress: %d %% done\n",progress++); // Print progress of main loop
      fflush(stdout);
    }
    /* ---------------------------End of main loop------------------------------------ */
  }

  for(i = 0; i < N; i++) {
    current_time = (double)i*timestep;
    fprintf(x_file, "%i \t %.6f \n", current_time , x[i]);
    fprintf(v_file, "%i \t %.6f \n", current_time , v[i]);
  }
  
  return 0; 
}
