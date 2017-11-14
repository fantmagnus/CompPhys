
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include "initfcc.h"
#include "alpotential.h"

#define PI 3.141592653589
#define k_B 0.000086173303
#define nbr_of_timesteps_eq 10000 /* nbr_of_timesteps+1 = power of 2, for best speed */
#define nbr_of_timesteps 100
#define nbr_of_atoms 256

void temp_scale(double v[][3], double n, double T_eq, double tau_T, double T, double timestep)
{
  double alpha_T;
  alpha_T = 1 + timestep/tau_T*(T_eq - T)/T;
  for(int i = 0; i < n; i ++) {
    for( int j = 0; j < 3; j++ ) {
      v[i][j] = alpha_T * v[i][j];  
    }
  }
}

void press_scale(double positions[][3], double P_eq, double p, double tau_P, double kappa_T, double timestep, double n, double tot_volume)
{
  double alpha_P = 0;
  alpha_P = 1 - kappa_T * timestep * (P_eq-p)/tau_P;
  double exponent = 1/3;
  for(int i = 0; i < n; i ++) {
    for( int j = 0; j < 3; j++ ) {
      positions[i][j] = pow(alpha_P,exponent) * positions[i][j];  
    }
  }
  tot_volume = tot_volume * alpha_P;
}
/* Main program */
int main()
{ 
  // Declaration of variables
  double timestep;
  double timestep_sq,current_time;
  double m;
  int n_cell; // Number of atoms in unit cell
  double cell_length;
  double tot_length;
  double T_eq;
  double tau_T;
  double tau_P;
  double W;
  double tot_volume;
  double P_eq;
  double kappa_T;
  double nbr_of_cells;
  double mean_T;
  double mean_P;
  double mean_E_kin; 

  /* declare file variable */
  FILE *energy_file;
  FILE *temp_file;
  FILE *press_file; 

  /* displacement, velocity and acceleration */
 	   
  /* Allocating memory for large vectors */
  double (*positions)[3] = malloc(sizeof(double[nbr_of_atoms][3]));
  double (*v)[3] = malloc(sizeof(double[nbr_of_atoms][3]));
  double (*F)[3] = malloc(sizeof(double[nbr_of_atoms][3]));
  double *E_pot = malloc(nbr_of_timesteps * sizeof(double));
  double *E_kin_eq = malloc(nbr_of_timesteps_eq * sizeof(double));
  double *E_kin = malloc(nbr_of_timesteps * sizeof(double));
  double *temp = malloc(nbr_of_timesteps_eq * sizeof(double));
  double *press = malloc(nbr_of_timesteps_eq * sizeof(double));
  
  

  /* Set variables */
  timestep = 0.001;
  timestep_sq = timestep * timestep;
  cell_length = 4.045;
  n_cell = 4;
  m = 27*1.0364*0.0001;
  tot_length = cell_length*n_cell;
  tot_volume = tot_length * tot_length * tot_length;
  T_eq = 500 + 273.15;
  P_eq = 6.2415 * 0.0000001; // 1 atm är la typ 100 kPa som är typ detta i eV/Å^3 
  tau_T = timestep * nbr_of_timesteps_eq * 0.01;
  tau_P = timestep * nbr_of_timesteps_eq * 0.1;
  kappa_T = 21.7576;
  nbr_of_cells = nbr_of_atoms / n_cell;
  mean_E_kin = 0;
 

  /* Initialize lattice*/
  init_fcc(positions, n_cell, cell_length);
  
  // Initialization for lattice and random deviation in positions
  double u;
  const  gsl_rng_type *T; /*  static  info  about  rngs */
  gsl_rng *q; /* rng  instance  */
  gsl_rng_env_setup (); /*  setup  the  rngs */
  T = gsl_rng_default; /*  specify  default  rng */
  q = gsl_rng_alloc(T); /*  allocate  default  rng */
  gsl_rng_set(q,time(NULL)); /*  Initialize  rng */
  u = gsl_rng_uniform(q); /*  generate  random  number
			      (repeatable) */
  for(int i = 0; i < nbr_of_atoms; i++)
    {
      for(int j = 0; j < 3; j++)
	{
	  positions[i][j] += cell_length/20*gsl_rng_uniform(q)-cell_length/20*gsl_rng_uniform(q);
	  v[i][j] = 0;
	}
    }
  
  /* Calculate initial accelerations (forces)  based on initial displacements */
  get_forces_AL(F, positions, tot_length, nbr_of_atoms); 
  E_pot[0] = get_energy_AL(positions, tot_length, nbr_of_atoms);

  /* Equilibration of the system */
  for (int i = 1; i < nbr_of_timesteps_eq + 1; i++) {
    for (int x = 0; x < 3; x++) { // For three space directions 
      for (int j = 0; j < nbr_of_atoms; j++) {
	v[j][x] += timestep * 0.5 * F[j][x]/m; // v(t+0.5*dt)
	positions[j][x] += timestep * v[j][x];  // q(t+dt)
      }      
    }
    // a(t+dt)
     get_forces_AL(F, positions, tot_length, nbr_of_atoms);
    for(int j = 0; j < nbr_of_atoms; j++) {
      for(int x = 0; x < 3; x++){
	v[j][x] += timestep * 0.5 * F[j][x]/m; // v(t+dt)
	E_kin_eq[i] += m * v[j][x] * v[j][x]*0.5;
      }
    }
    temp[i] = 2 / (3 * nbr_of_atoms * k_B) * E_kin_eq[i];
    temp_scale(v, nbr_of_atoms, T_eq, tau_T, temp[i], timestep);
    W = get_virial_AL(positions, cell_length, nbr_of_atoms);
    press[i] = (W + nbr_of_atoms * k_B * temp[i]) / (tot_volume);
    printf("%e\n",press[i]);
    press_scale(positions, P_eq, press[i], tau_P, kappa_T, timestep, nbr_of_atoms, tot_volume);
  }
  /* Reset the matrixes */
  E_kin[0] = E_kin_eq[nbr_of_timesteps_eq - 1];
  E_pot[0] = get_energy_AL(positions, tot_length, nbr_of_atoms);
  
  W = 0; 

  /* Simulation */
  /* timesteps according to velocity Verlet algorithm */
  for (int i = 1; i < nbr_of_timesteps + 1; i++) {
    for (int x = 0; x < 3; x++) { // For three space directions 
      for (int j = 0; j < nbr_of_atoms; j++) {
	v[j][x] += timestep * 0.5 * F[j][x]/m; // v(t+0.5*dt)
	positions[j][x] += timestep * v[j][x];  // q(t+dt)
      }      
    }
    // a(t+dt)
    get_forces_AL(F, positions, tot_length, nbr_of_atoms);
    for(int j = 0; j < nbr_of_atoms; j++) {
      for(int x = 0; x < 3; x++){
	v[j][x] += timestep * 0.5 * F[j][x]/m; // v(t+dt)
	E_kin[i] += m * v[j][x] * v[j][x]*0.5;
      }
    }
    E_pot[i] = get_energy_AL(positions, tot_length, nbr_of_atoms);
    W += get_virial_AL(positions, cell_length, nbr_of_atoms);
  }

  for(int i = 0; i < nbr_of_timesteps; i++){
    mean_E_kin += E_kin[i]/nbr_of_timesteps;
  }

  mean_T = 2*mean_E_kin/(3*k_B*nbr_of_atoms);
  printf("T =  %f K\n", mean_T); 
  mean_P = (nbr_of_atoms*k_B*mean_T + W/nbr_of_timesteps)/(tot_volume);
  printf("P = %f\n", mean_P);
  
  /* Print temperature, pressure and and energy data to output file */
  energy_file = fopen("energy.dat","w");
  temp_file = fopen("temp.dat","w");
  press_file = fopen("press.dat","w");
  for (int i = 0; i < nbr_of_timesteps + 1; i++) {
    current_time = i * timestep;
    fprintf(energy_file, "%.4f \t %e \t %e \t %e \n", current_time, E_pot[i], E_kin[i], E_pot[i]+E_kin[i]);
  }
  fclose(energy_file);
  for (int i = 0; i < nbr_of_timesteps_eq + 1; i++) {
  current_time = i * timestep;
  fprintf(temp_file, "%.4f \t %e \n", current_time, temp[i]);
  fprintf(press_file, "%.4f \t %e \n", current_time, press[i]);
  }
  fclose(temp_file);
  fclose(press_file); 
  /* Free allocated memory */
  gsl_rng_free(q); 
  free(positions); 
  free(v);
  free(E_pot);
  free(E_kin);
  free(E_kin_eq);
  free(F);
  free(temp);
  free(press);

}
