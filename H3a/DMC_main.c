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
  return W; 
}

/* Branching function */
int make_branch(double walkers[], int m[], int *N, double E_T, double timestep, gsl_rng *q) {
  int M = 0;
  for (int i = 0; i < *N; i++) {
    double u = gsl_rng_uniform(q);
    m[i] = (int)(weight_func(walkers[i], E_T, timestep) + u);  // Calculate number of walkers to replace walker att position walkers[i]
    M = M + m[i];  // Calculate new number of walkers
  }
  return M;
}

void place_bin(double bins[], double walkers[], int nbr_of_bins, double x_min, double x_max, int N) {
  double size_of_bin = (x_max-x_min)/nbr_of_bins; 
  for(int i = 0; i < N; i++) { // Place al walkers in a bin
    int bool = 1;
    int index = 0;
    double x = x_min;
    if(walkers[i] > x_max || walkers[i] < x_min) { // If walkers are outside the interval
      bool = 0; 
    }
    while (bool == 1) { // Until walker position is placed in a bin
      if(walkers[i] > x && walkers[i] < x + (double)size_of_bin) {
	bins[index] += 1;
	bool = 0; 
      }
      index += 1;
      x = x + size_of_bin;
    }
  }
}
/* Calculate trial energy */
double calc_E(double E_T, int N, int N_0, double alpha, double timestep) {
  double Nk = (double)N;
  double N0 = (double)N_0;
  return E_T-alpha/timestep*log(Nk/N0);
}
/* Makes a diffusive step with one walker */
void diffusive_step(double walkers[], int index, double timestep, gsl_rng *q) {
  double G = gsl_ran_ugaussian(q); // Gaussian rn
  walkers[index] = walkers[index] + sqrt(timestep)*G;
}

int main() {
  /*Initializing variables */
  int i,j,k;
  int N_0 = 300;
  int N = N_0;
  int M;
  int index;
  int progress = 0;
  double timestep = 0.01;
  double alpha = 5e-5;
  int nbr_of_timesteps = 20000;
  int nbr_of_eq_steps = 0; // Number of equilibration steps
  double E_T = 1; // Initial value for trial energy
  double E_0 = 0.5; // Correct ground-state energy
  /* Set interval for wave-function and number of bins*/
  double x_min = -5;
  double x_max = 5;
  int nbr_of_bins = 1000;
  /*-------------------------------*/
  double size_of_bin = (x_max-x_min)/nbr_of_bins;
  double bin_sum = 0;
  /* Initialize random number */
  const gsl_rng_type*o;
  gsl_rng *q;
  gsl_rng_env_setup();
  o = gsl_rng_default;
  q = gsl_rng_alloc(o);
  gsl_rng_set(q, time(NULL));
  /*---------------------------*/
  double* walkers = malloc(N_0 * sizeof(double));
  double* bins = malloc(nbr_of_bins * sizeof(double));
  /* Initialize walkers */
  for (i = 0; i < N_0; i++){
    double G = gsl_ran_ugaussian(q);
    walkers[i] = G; 
  }
  /* Files */
  FILE* energy_file = fopen("E.dat","w");
  FILE* walkers_file = fopen("M.dat","w");
  FILE* dist_file = fopen("dist.dat","w");
  FILE* bin_file = fopen("bin.dat","w");
  
  N = N_0;
  /*----------------------- Main loop-------------------------- */
  for (j = 0; j < nbr_of_timesteps; j++) {
    E_T = calc_E(E_T, N, N_0, alpha, timestep);
    if (j > nbr_of_eq_steps - 1) {
      fprintf(energy_file,"%f\t %f\n", E_T, timestep*j); // Save energy 
      place_bin(bins, walkers, nbr_of_bins, x_min, x_max, N); // Generate histogram 
    }
    for (i = 0; i < N; i++){ // For all walkers
      diffusive_step(walkers, i, timestep, q);
    }
    /* Branching */
    double* temp_walkers = malloc(N * sizeof(double));
    for (i = 0; i < N; i++){
      temp_walkers[i] = walkers[i];
    }
    int* m = malloc(N * sizeof(double));
    M = make_branch(walkers, m, &N, E_T, timestep, q);
    walkers = realloc(walkers, M * sizeof(double));
    index = 0;
    for (i = 0; i < N; i++) { // Fill new walker array 
      for(k = 0; k < m[i]; k++) {
	walkers[index] = temp_walkers[i];
	index = index + 1;
      }
    }
    fprintf(walkers_file, "%d\t %f\n", N, timestep*j);
    N = M;
    /* Print progress of main loop */
    if (j%((nbr_of_timesteps)/100) == 0){
      printf("\rProgress: Simulation %d %% done ",progress++); // Print progress of main loop
      fflush(stdout);
    }
    /* Free allocated arrays */
    free(temp_walkers);
    free(m);
  }
  /*--------------------------------------------------------*/

  for (i = 0; i < nbr_of_bins; i++) {
    bins[i] = bins[i]/(nbr_of_timesteps); // Time averaging of wave-function
    bin_sum +=  bins[i];    
  }
  
  for (i = 0; i < nbr_of_bins; i++) {
    bins[i] = bins[i]/(bin_sum*size_of_bin); // Normalization of wave-function
    bins[i] = bins[i] * bins[i]; // Calculate density from wave-function
  }
  
  for (i = 0; i < N; i++){
    fprintf(dist_file,"%f\n", walkers[i]);
  }

  /* Save density */
  for (i = 0; i < nbr_of_bins; i++){
    fprintf(bin_file,"%f\t %f\n", x_min+i*size_of_bin+size_of_bin/2, bins[i] );
  }
  
  printf("E_T = %f\n", E_T);
  printf("It should be %f\n", E_0);
  if (N == 0) { // If all walkers were eliminated
    printf("WALKER HOLOCAUST! N = 0\n");
  }
    else {
    printf("There are %d walkers.\n", N);
  }
}
