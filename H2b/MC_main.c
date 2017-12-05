#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <math.h>
#define PI 3.14159265359

/* Returns the distance between two points in space */
double distance(double r1[3], double r2[3]) {
  double dist = sqrt(pow(r2[0]-r1[0],2)+ pow(r2[1]-r1[1],2)+pow(r2[2]-r1[2],2));
  return dist;
}

/* Returns the length of a vector */
double norm(double r[3]) {
  double norm = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
  return sqrt(norm);
}

/* Returns local energy for a certain configuration and alpha */
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

/* Returns value of the trial function for a certain configuration and alpha */
double weightfunc(double r1[3], double r2[3], double alpha) {
  double r12 = distance(r1, r2);
  double psi_t = exp(-2*norm(r1))*exp(-2*norm(r2))*exp(0.5*r12/(1+alpha*r12));
  double abs_psi_t = sqrt(psi_t * psi_t);
  return abs_psi_t * abs_psi_t;
}

/* Determines the next configuration */
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

/* Returns cosine of angle between two vectors */
double calc_x(double r1[3], double r2[3]) {
  double scal=0;
  for (int i = 0; i < 3; i++){
    scal += r1[i] * r2[i];
  }
  double theta = acos( scal / (norm(r1) * norm(r2)));
  double x = cos(theta);
  return x;
}

/*Return derivative of logarithm of the trial function for a certain configuration and alpha */
double dPhi(double r1[3], double r2[3], double alpha) {
  double r12 = distance(r1, r2);
  double dPhi = -r12*r12*pow(1+alpha*r12, -2.0)/2;
  return dPhi; 
}


int main () {
  /*Variable declarations */
  int i,j,k;
  int dyn = 0; //For dynamic alpha
  int N = 1e6; 
  int count = 0;
  int progress = 0;
  int nbr_skipped_states = 1e5;
  double alpha = 0.145;
  int n_alpha = 1;//100;
  double delta = 1;
  double r1[3], r2[3];
  double E_sum = 0;
  int mod_factor = 1000;
  int k_span = 600;
  double mean_E=0;
  double mean_sq_E=0;
  double E_ik=0; 
  double Phi[k_span];
  double *E = malloc((N) * sizeof (double));
  int B, b;
  int j_span;
  int n_B = 500;
  double s[n_B];
  double mean_F=0;
  double mean_sq_F=0;
  int n_s=12;
  double var_E, var_I;
  double inv_Hess = 1;
  double beta = 0.8;
  double meanf_E;
  double meanf_div;
  double meanf_Ediv;
  double div_sum = 0;
  double Ediv_sum = 0; 
 
  
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
  FILE *energy_file = fopen("E.dat","w");
  FILE *block_file = fopen("S.dat","w");
  /* Set initial state */
  // May be redunadant

  for (j=0; j < n_alpha; j++) { // for different alphas
    if (n_alpha > 1) {
      alpha = 0.05+0.20/(n_alpha-1)*j;
    }

    /* Reset */
    E_sum = 0;
    progress=0;
    printf("\r");
    mean_E = 0;
    meanf_E = 0;
    meanf_Ediv = 0;
    meanf_div = 0;
    mean_sq_E = 0;
    
    /* Throw away first states */
    for (i = 1; i < nbr_skipped_states+1; i++){
      next_state(r1, r2, alpha, delta, &count, q);
      fprintf(energy_file,"%.7f\n", calc_E(r1,r2,alpha));
      if (i%mod_factor == 0){
	//printf("Local energy after %d steps is %.5f\n",i,calc_E(r1,r2,alpha)); 
      }
    }

    count = 0;
  
    /* Step in configuration space */
    for (i = 0; i < N; i++){
      /* Metropolis algorithm */
      next_state(r1, r2, alpha, delta, &count, q);
      E[i] = calc_E(r1, r2, alpha);
      E_sum += E[i];
      div_sum += dPhi(r1, r2, alpha);
      Ediv_sum += E[i]*dPhi(r1, r2, alpha);
      if (n_alpha == 1){ // For fixed alpha simulations
	fprintf(dist_file,"%.4f\n%.4f\n", norm(r1), norm(r2));
	fprintf(corr_file,"%.4f\n", calc_x(r1,r2));
	fprintf(energy_file,"%.7f\n", E[i]);
      }
      if (i%((N-1)/100) == 0){
	printf("\rProgress: Simulation %d of %d %d %% done  ",j+1,n_alpha,progress++); // Print progress of main loop
	fflush(stdout);
      }
      mean_E += E[i]/N;
      mean_sq_E += E[i]*E[i]/N;
      if(dyn > 0) {
	meanf_E = E_sum/(i+1);
	meanf_div = div_sum/(i+1);
	meanf_Ediv = Ediv_sum/(i+1);
	alpha = alpha - inv_Hess*pow((i+1),-beta)*2*(meanf_Ediv - meanf_div*meanf_E);
	
      }
    }
   /* Calculate mean and variance */
    var_E = mean_sq_E - mean_E*mean_E;
    var_I = n_s*var_E/N;
    if (n_alpha > 1){
      fprintf(energy_file,"%.4f\t %.7f\t %.7f\n",alpha, E_sum/N, var_I);
    }
  }
  
  fclose(energy_file);
  fclose(corr_file);
  fclose(dist_file);
  /* Error estimate */

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

  /* Calcualte statistical inefficiency with block averaging */
  for (b = 1; b < n_B; b++){
    B = b*10; 
    j_span = (int)N/B;
    double *F = malloc(j_span * sizeof (double));
    for (j = 0; j < j_span; j++){
      for (i = 0; i < B; i++){
	F[j] += E[j*B + i];
      }
      F[j] = F[j]/B;
      mean_F += F[j]/j_span;
      mean_sq_F += F[j]*F[j]/j_span;
    }
    s[b] = B*(mean_sq_F - mean_F*mean_F)/(mean_sq_E - mean_E*mean_E);
    fprintf(block_file, "%i \t %.6f \n", B, s[b]);
    free(F);
    mean_F = 0;
    mean_sq_F = 0;
  }
  fclose(block_file);
 
    
  printf("\n");
  printf("Distance 1: %.4f\n",norm(r1));
  printf("Distance 2: %.4f\n",norm(r2));
  printf("E should be about -3 a.u. \n");
  printf("E is %f a.u. \n", E_sum/N);
  printf("Alpha is %f. \n", alpha);
  printf("The stepping percentage is %.0f%%\n",100*(double)count/N);
  printf("The standard deviation is %f. \n", sqrt(var_I));
  printf("The angle between the electrons is %f. \n",acos(calc_x(r1,r2))/PI*180);
  /* Free memory */
  free(E);
  
  return 0;
}
