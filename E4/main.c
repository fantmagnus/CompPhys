#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <math.h>
#define PI 3.14159265359
#define k_B 1.38064852e23 // m2 kg s-2 K-1

double calc_a(double x, double v, double f, double eta, double xi) {
  double omega_0 = f*2*PI;
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
  double timestep = 0.0001; // in ms
  double f_0 = 3; //Frequency of 3 kHz (1/ms)
  double tau = 48.5e-3; // 48.5e-3 ms for case A, 147.3e-3 ms for case B
  double eta = 1/tau;
  double c_0 = exp(-eta * timestep);
  double sqrt_c_0 = sqrt(c_0);
  double R = 1.395e-3; // mm
  double rho = 2.65e3; // g/mm^3
  double m = rho*3*PI*R*R*R/3; // g
  double T = 295; //K 
  double v_th = sqrt(k_B*T/m );
  double current_time; 
  int N = 1e5;
  int progress = 0;
  double G1, G2;
  double G3, G4;
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
  FILE* G_file = fopen("G_data.dat","w");
  
  x[0] = 0.1;
  x_T = x[0];
  v[0] = 2.0;
  v_T = v[0];
  gen_gaussian(&G3, &G4, q);
  a = calc_a(x_T, v_T, f_0, eta, G3);
  printf("The initial acceleration is %e mm/(ms)^2\n",a);
  for(i = 1; i < N; i ++) {
    gen_gaussian(&G1, &G2, q);
    gen_gaussian(&G3, &G4, q);
    /*Verlet velocity algorithm*/
    /*Velocity v(t+dt/2) and displacement q(t+dt) */
    v_tilde = 0.5*a*timestep + sqrt_c_0*v_T + v_th*sqrt(1-c_0)*G1;

    x_T = x_T + v_tilde*timestep;
    x[i] = x_T;
	
    /*Calculate the acceleration for v(t+dt) */
    a = calc_a(x_T, v_tilde, f_0, eta, G3);
    /*Velocity v(t+dt) */
    
    v_T = 0.5*sqrt_c_0*a*timestep + sqrt_c_0*v_tilde + v_th*sqrt(1-c_0)*G2;
    v[i] = v_T;

    fprintf(G_file, "%e\n", G1);
    
    if ((i-1)%((N-1)/100) == 0) {
      printf("\rProgress: %d %% done",progress++); // Print progress of main loop
      fflush(stdout);
    }
  }
  /* ---------------------------End of main loop--------------------------- */
  printf("\n");
  for(i = 0; i < N; i++) {
    current_time = i*timestep;
    fprintf(x_file, "%e \t %.6f \n", current_time , x[i]);
    fprintf(v_file, "%e \t %.6f \n", current_time , v[i]);
  }
  fclose(x_file);
  fclose(v_file);
  fclose(G_file);

  free(x);
  free(v);
  
  return 0; 
}
