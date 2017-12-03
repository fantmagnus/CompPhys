#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <math.h>
#define PI 3.14159265359

double distance(double r1[3], double r2[3]) {
  double dist = sqrt(pow(r2[0]-r1[0],2)+ pow(r2[1]-r1[1],2)+pow(r2[2]-r1[2],2));
  return dist;
}

double norm(double r[3]) {
  double norm = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
  return sqrt(norm);
}

double calc_E(double r1[3], double r2[3], double alpha) {
  double r12 = distance(r1, r2);
  double r1_hat[3], r2_hat[3], r_hat[3], r[3];
  double scal=0; 
  for(int i = 0; i < 3; i++ ) {
    r1_hat[i] = r1[i]/norm(r1);
    r2_hat[i] = r2[i]/norm(r2);
    r_hat[i] = r1_hat[i] - r2_hat[i];
    r[i] = r1[i] - r2[i];
    scal += r_hat[i]*r[i]; 
  }
  
  double E_L = -4 + scal/(r12*pow(1+alpha*r12,2)) - 1/(r12*pow(1+alpha*r12,3)) - 1/(4*pow(1+alpha*r12,4)) + 1/r12;
  return E_L;
}

double weightfunc(double r1[3], double r2[3], double alpha) {
  double r12 = distance(r1, r2);
  double psi_t = exp(-2*norm(r1))*exp(-2*norm(r2))*exp(0.5*r12/(1+alpha*r12));
  double abs_psi_t = sqrt(psi_t * psi_t);
  return abs_psi_t * abs_psi_t;
}

void next_state(double r1[3], double r2[3], double alpha, double delta, int *count, gsl_rng *q){
  /* create trial step */
  double u1, u2; 
  double r1_t[3];
  double r2_t[3];
  for (int i = 0; i < 3; i++) {
    u1 = gsl_rng_uniform(q);
    u2 = gsl_rng_uniform(q);
    r1_t[i] = r1[i]+delta*(u1-0.5);
    r2_t[i] = r2[i]+delta*(u2-0.5);
  }
  /* Calulate relative probabilities */
  double p_m = weightfunc(r1, r2, alpha);
  double p_t = weightfunc(r1_t, r2_t, alpha);
  /* Accept or reject new state */
  double xi = gsl_rng_uniform(q);
  if (p_t/p_m > xi) { // if true accept
    for (int i = 0; i < 3; i++){
      r1[i] = r1_t[i];
      r2[i] = r2_t[i];
    }
    *count = *count+1;
  }
}

double calc_x(double r1[3], double r2[3]) {
  double scal=0;
  for (int i = 0; i < 3; i++){
    scal += r1[i] * r2[i];
  }
  double theta = acos( scal / (norm(r1) * norm(r2)));
  double x = cos(theta);
  return x;
}

int main () {

  /*Variable declarations */
  int i,k;
  int N = 1e7; 
  int count = 0;
  int progress = 0;
  int nbr_skipped_states = 1e7;
  double alpha = 0.1;
  double delta = 1;
  double r1[3], r2[3];
  double E_sum = 0;
  int mod_factor = 1000000;
  int k_span = 300;
  double mean_E=0;
  double mean_sq_E=0;
  double E_ik=0; 
  double Phi[k_span];
  double *E = malloc((N) * sizeof (double));
  
  
  /* gsl random number setup */ 
  const gsl_rng_type*T;
  gsl_rng *q;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  q = gsl_rng_alloc(T);
  gsl_rng_set(q, time(NULL));
  /* Files */
  FILE *dist_file = fopen("dist.dat","w");
  FILE *corr_file = fopen("corr.dat","w");
  FILE *corrfunc_file = fopen("phi.dat","w");

  /* Set initial state */
  // May be redunadant


  
  /* Throw away first states */
  for (i = 1; i < nbr_skipped_states+1; i++){
    next_state(r1, r2, alpha, delta, &count, q);
    if (i%mod_factor == 0){
      printf("Local energy after %d steps is %.5f\n",i,calc_E(r1,r2,alpha)); 
    }
  }
  printf("%.4f\t %.4f\t %.4f\t %.4f\n",r1[0],r1[1],r1[2],norm(r1));
  printf("%.4f\t %.4f\t %.4f\t %.4f\n",r2[0],r2[1],r2[2],norm(r2));
  count = 0;
  
  /* Step in configuration space */
  for (i = 0; i < N; i++){
    /* Metropolis algorithm */
    next_state(r1, r2, alpha, delta, &count, q);
    E[i] = calc_E(r1, r2, alpha);
    E_sum += E[i];
    fprintf(dist_file,"%.4f\n%.4f\n", norm(r1), norm(r2));
    fprintf(corr_file,"%.4f\n", calc_x(r1,r2));

    if (i%((N-1)/100) == 0){
      printf("\rProgress: %d %% done",progress++); // Print progress of main loop
      fflush(stdout);
    }
  }

  /* Error estimate */
  /* Calculate mean and variance */
    for(i = 0; i < N; i++) {
      mean_E += E[i]/N;
      mean_sq_E += E[i]*E[i]/N; 
  }
    /* Calculate correlation function */
    for(k = 0; k < k_span; k++) {
      for(i = 0; i < N-k; i++) {
	E_ik += E[i+k]*E[i];
      }
      E_ik = E_ik/(N-k);
      Phi[k] = (E_ik - mean_E*mean_E)/(mean_sq_E - mean_E*mean_E);
      fprintf(corrfunc_file, "%i \t %.6f \n", k, Phi[k]);
    }
    fclose(corrfunc_file);
    
    
  printf("Distance 1: %.4f\n",norm(r1));
  printf("Distance 2: %.4f\n",norm(r2));
  printf("E should be about -3 a.u. \n");
  printf("E is %f a.u. \n", E_sum/N);
  printf("The stepping procentage is %.0f%%\n",100*(double)count/N);

  /* Free memory */
  free(E);
  
  return 0;
}
