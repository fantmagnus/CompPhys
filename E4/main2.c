#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <math.h>
#define PI 3.14159265359
#define k_B 1.38064852e-17 // m2 kg s-2 K-1

double calc_a(double x, double f) {
  double a = -4*PI*PI*f*f*x;
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
  int nbr_of_particles = 3000;
  double v_T[nbr_of_particles];
  double x_T[nbr_of_particles];
  double v_tilde = 0;
  double timestep = 0.0001; // in ms
  double f_0 = 3; //Frequency of 3 kHz (1/ms)
  double tau = 48.5e-3; // 48.5e-3 ms for case A, 147.3e-3 ms for case B
  double eta = 1/tau;
  double c_0 = exp(-eta * timestep);
  double sqrt_c_0 = sqrt(c_0);
  double R = 1.395; // mm
  double rho = 2.65e-15; // kg/mum^3
  double m = rho*3*PI*R*R*R/3; // g
  double T = 295; //K 
  double v_th = sqrt(k_B*T/m );
  double current_time; 
  int N = 1e5;
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


  double *mean_x = malloc(N * sizeof(double));
  double *mean_v = malloc(N * sizeof(double));
  double mean_x_sq = 0;
  double mean_v_sq = 0;

  double *var_x = malloc(N * sizeof(double)); 
  double *var_v = malloc(N * sizeof(double));
  
  double (*x)[5] = malloc(sizeof(double[N][5]));
  double (*v)[5] = malloc(sizeof(double[N][5]));
  double *a = malloc(nbr_of_particles * sizeof(double));
  
  FILE* x_file = fopen("x_data.dat","w");
  FILE* v_file = fopen("v_data.dat","w");

  
 
  for (j = 0; j < nbr_of_particles; j++) {
    if (j < 5){
      v[0][j] = 2.0;
      x[0][j] = 0.1;
    }
    x_T[j] = x[0][0];
    v_T[j] = v[0][0];
    mean_x[0] += x_T[j]/nbr_of_particles; 
    mean_v[0] += v_T[j]/nbr_of_particles;
    a[j] = calc_a(x_T[j], f_0);
  }
  var_x[0] = 0;
  var_v[0] = 0;
  
  for (i = 1; i < N; i++){
      for(j = 0; j < nbr_of_particles; j ++) {    
      gen_gaussian(&G1, &G2, q);
      /*Verlet velocity algorithm*/
      /*Velocity v(t+dt/2) and displacement q(t+dt) */
      v_tilde = 0.5*a[j]*timestep + sqrt_c_0*v_T[j] + v_th*sqrt(1-c_0)*G1;

      x_T[j] = x_T[j] + v_tilde*timestep;
	
      /*Calculate the acceleration for v(t+dt) */
      a[j] = calc_a(x_T[j], f_0);
      /*Velocity v(t+dt) */
    
      v_T[j] = 0.5*sqrt_c_0*a[j]*timestep + sqrt_c_0*v_tilde + v_th*sqrt(1-c_0)*G2;

      mean_x[i] += x_T[j]/nbr_of_particles; 
      mean_v[i] += v_T[j]/nbr_of_particles;
      mean_x_sq += x_T[j] * x_T[j]/nbr_of_particles; 
      mean_v_sq += v_T[j] * v_T[j]/nbr_of_particles;

      if (j < 5) {
	x[i][j] = x_T[j];
	v[i][j] = v_T[j];
      }
    }
      var_x[i] = mean_x_sq - mean_x[i]*mean_x[i];
      var_v[i] = mean_v_sq - mean_v[i]*mean_v[i];
      mean_x_sq = 0;
      mean_v_sq = 0;
      if ((i+1)%((N-1)/100) == 0){
	printf("\rProgress: Simulation %d %% complete  ",progress++); // Print progress of main loop
	fflush(stdout);
      }
      
  }
  printf("\n");


  for(i = 0; i < N; i++) {
    current_time = i*timestep;
    fprintf(x_file, "%e\t", current_time);
    fprintf(v_file, "%e\t", current_time);
    for (j = 0; j < 5; j++){ 
      fprintf(x_file, "%e \t", x[i][j]);
      fprintf(v_file, "%e \t", v[i][j]);
    }
    fprintf(x_file, "%e\t %e\n", mean_x[i], var_x[i]);
    fprintf(v_file, "%e\t %e\n", mean_v[i], var_v[i]);
  }
  fclose(x_file);
  fclose(v_file);

  free(x);
  free(v);
  free(var_x);
  free(mean_x);
  free(var_v);
  free(mean_v);
  return 0; 
}
