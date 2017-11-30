#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include "func.h"
#include "alpotential.h"

#define PI 3.141592653589
#define nbr_of_timesteps 32769 //32769 //44000 //Should be power of 2 + 1 for improved simulation time
#define nbr_of_eq_timesteps 20000 //20000 //30000
#define nbr_of_atoms 256

/* Equilibration routine for the temperature and pressure */
double alpha_T(double T_eq, double T, double tau_T, double timestep) {
  double alpha_T = 1 + timestep/tau_T*(T_eq-T)/T;
  return alpha_T; 
}

double alpha_P(double kappa_T, double P_eq, double P, double tau_P, double timestep) {
  double alpha_P = 1 - kappa_T*timestep/tau_P*(P_eq-P);
  return alpha_P; 
}



/* Main program */
int main()
{ 
  // Declaration of variables
  int n_cell;
  int i,j,k; 
  double timestep;
  double timestep_sq,current_time;
  double m;
  double cell_length;
  double tot_length;
  double T_eq;
  double tau_T;
  double tau_P; 
  double k_B;
  double P_eq;
  double kappa;
  double W;
  int simul = nbr_of_timesteps-nbr_of_eq_timesteps;
  double mean_P=0;
  double C_v;
  int progress=0;
  double T_mean=0;
  double E_kin_mean=0;
  double E_pot_mean=0;
  double E_kin_var;
  double E_kin_sqr_mean=0;
  int b = 4; // Temperature is raised with b*delta_T
  
  
  /* declare file variable */
  FILE *energy_file;
  FILE *temp_file;
  FILE *press_file;
  FILE *traj_file;
  FILE *radial_file;
  FILE *pos_file;
  FILE *vel_file;
 
  /* Allocating memory for large vectors */
  double (*positions)[3] = malloc(sizeof(double[nbr_of_atoms][3]));
  double (*v)[3] = malloc(sizeof(double[nbr_of_atoms][3]));
  double (*F)[3] = malloc(sizeof(double[nbr_of_atoms][3]));
  double *E_pot = malloc(nbr_of_timesteps * sizeof(double));
  double *E_kin = malloc(nbr_of_timesteps * sizeof(double));
  double *Temp = malloc(nbr_of_timesteps * sizeof(double));
  double *Press = malloc(nbr_of_timesteps * sizeof(double));
  double *Trajectory = malloc(nbr_of_timesteps * sizeof(double));
  /* Set variables */
  timestep = 0.001;
  timestep_sq = timestep * timestep;
  cell_length = 4.045;
  n_cell = 4;
  m = 27*1.0364*0.0001;
  tot_length  =  cell_length*n_cell;
  tau_T = 0.5;
  tau_P = 0.1; 
  T_eq = 500+273.15+b*20; // Task 6
  P_eq = 6.324e-7; 
  k_B = 8.617330e-5;
  kappa = 2.219;
  // tot_length =  17.021875; // From equilibration run (liquid)
  tot_length =  16.357553; // From equilibration run (solid)
  int nbr_r_steps = 350;
  int bin = 0;
  double dr = 0;
  double r = 0;
  double N_ideal;
  double N_rk[nbr_r_steps*2];
  double g[nbr_r_steps];
  double sxi, syi, szi, sxjk, syjk, szjk;
  double *sx = malloc(nbr_of_atoms * sizeof (double));
  double *sy = malloc(nbr_of_atoms * sizeof (double));
  double *sz = malloc(nbr_of_atoms * sizeof (double));
  
  
  /*  Initialize lattice and a small displacement*/
  init_fcc(positions, n_cell, cell_length);

  
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
  /*________________________________________________*/

  /* Read in velocities and positions from equilibration run */
    FILE *myFile1;
    FILE *myFile2;
    myFile1 = fopen("pos.dat", "r");
    myFile2 = fopen("vel.dat", "r");
    //read file into array
    for (int i = 0; i < nbr_of_atoms ; i++)
      {
	for (int j = 0; j < 3; j++){
	  fscanf(myFile1, "%lf", &positions[i][j]);
	  fscanf(myFile2, "%lf", &v[i][j]);
	}
      }
  /*_____________________________*/

  /* Calc initial acceleration */ 
  printf("Initial cell length is %f\n", tot_length); 
  get_forces_AL(F, positions,tot_length, nbr_of_atoms);
  
  /* --------------------------------Main simulation loop--------------------------------- */
  for(i = 0; i < nbr_of_timesteps + 1; i ++) {
      /*Verlet velocity algorithm*/
      /*Velocity v(t+dt/2) and displacement q(t+dt) */
      for(j = 0; j < nbr_of_atoms; j++){
	for(k = 0; k < 3; k++){
	  v[j][k] += timestep * 0.5 * F[j][k]/m;
	}
      }
      for(j = 0; j < nbr_of_atoms; j++){
	for(k = 0; k < 3; k++){
	  positions[j][k] += timestep * v[j][k];
	}
      }
      /*Calculate the acceleration for v(t+dt) */
      get_forces_AL(F, positions,tot_length, nbr_of_atoms);
      /*Velocity v(t+dt) */
      for(j = 0; j < nbr_of_atoms; j++){
	for(k = 0; k < 3; k++){
	  v[j][k] += timestep * 0.5 * F[j][k]/m;
	  E_kin[i] += m/2.0*v[j][k]*v[j][k];
	}
      }
      /* Calculate energies, temperature and pressure */
      E_pot[i] = get_energy_AL(positions, tot_length, nbr_of_atoms);
      Temp[i] = 2/(3*nbr_of_atoms*k_B) * E_kin[i];
      W = get_virial_AL(positions, tot_length, nbr_of_atoms);
      Press[i] = (Temp[i]*nbr_of_atoms*k_B+W)/(tot_length*tot_length*tot_length);

      /* Keep the temperature extra high for some steps to melt system */
      //      if(i > nbr_of_eq_timesteps/4){  
    //	T_eq=700+273.15+20*b;
      // }
      /*         Equilibration phase           */
      if(i < nbr_of_eq_timesteps) {
	for(j = 0; j < nbr_of_atoms; j++){
	  for(k = 0; k < 3; k++){
	    v[j][k] *= sqrt(alpha_T(T_eq, Temp[i], tau_T, timestep));
	    positions[j][k] *= cbrt(alpha_P(kappa, P_eq, Press[i], tau_P, timestep));
	  }
	}
       tot_length *= cbrt(alpha_P(kappa, P_eq, Press[i], tau_P, timestep));
       //printf("%f\n", tot_length); //This is used to check if the pressure is stable, i.e. the total length oscillates around the same value
      }
       /*________________________________________*/
       
   /*         Production phase           */
      if( i > nbr_of_timesteps-simul)
	{
	  int k = i - nbr_of_eq_timesteps; 
	  for(int j = 0; j < 3; j ++) {
	    Trajectory[k] += sqrt(positions[4][j]*positions[4][j]); // Save trajectory for atom 4 (Task 3,4)
	  }
	}
      /*Production for Task 7 with periodic boundary conditions and calculations for g(r)*/
      if(i > nbr_of_timesteps-simul)
	{
	for(j = 0; j < nbr_of_atoms; j++){
	  sx[j] = positions[j][0] / tot_length;
	  sy[j] = positions[j][1] / tot_length;
	  sz[j] = positions[j][2] / tot_length;
	}
	dr = tot_length/(2*nbr_r_steps);			
	for(j = 0; j < nbr_of_atoms - 1; j++){
	  sxi = sx[j] - floor(sx[j]);
	  syi = sy[j] - floor(sy[j]);
	  szi = sz[j] - floor(sz[j]);							
	  for(k = j + 1; k < nbr_of_atoms; k++){
	    /* Periodically translate atom j to positive quadrants and calculate distance to it. */
	    sxjk = sxi - (sx[k] - floor(sx[k]));
	    syjk = syi - (sy[k] - floor(sy[k]));
	    szjk = szi - (sz[k] - floor(sz[k]));
			  
	    /* Periodic boundary conditions. */
	    sxjk = sxjk - (int)floor(sxjk + 0.5);
	    syjk = syjk - (int)floor(syjk + 0.5);
	    szjk = szjk - (int)floor(szjk + 0.5);
					
	    r = tot_length*sqrt(sxjk*sxjk + syjk*syjk + szjk*szjk);
	    bin = ceil(r/dr) - 1;
	    N_rk[bin] += 2;
	  }	
	}
	/*_______________________________________________________________________________*/
	/* Production for Tasks 1,3,5,6 */
	E_kin_mean += E_kin[i];
	E_pot_mean += E_pot[i];
	T_mean += Temp[i];
	E_kin_sqr_mean += E_kin[i] * E_kin[i];
	/*_______________________________*/
	
	}
      if (i%(nbr_of_timesteps/100) == 0){
	printf("\rProgress: %d %% done\n",progress++); // Print progress of main loop
	fflush(stdout);
      }
  }
  /* ---------------------------End of main loop------------------------------------ */

  /* Calculate mean temperature, pressure, energy, C_v and print it */
  T_mean = T_mean/simul;
  E_kin_var = E_kin_sqr_mean/simul-(E_kin_mean/simul)*(E_kin_mean/simul);
  
  printf("Final cell length is %f\n", tot_length);
  printf("Average energy is %.10f\n", E_kin_mean/simul+E_pot_mean/simul);

  for(i = 0; i<simul; i++)
    {
      mean_P += Press[nbr_of_timesteps-simul+i];
    }
  printf("Average pressure is %.10f\n", mean_P/simul);

  printf("Average temperature is %f\n", T_mean);

  C_v = 3.0*nbr_of_atoms*k_B/2.0 / (1- 2.0/(3.0*nbr_of_atoms*k_B*k_B*T_mean*T_mean)*E_kin_var);
  printf("The heat capacity is %f\n", C_v);
  /*_______________________________________________________________*/

 /* Calculate radial distribution function */ 
  for (k = 1; k < nbr_r_steps; k++)
    {
      N_ideal = (nbr_of_atoms-1)/(tot_length*tot_length*tot_length) * 4*PI/3 * (3*k*k - 3*k + 1) * dr * dr * dr;  // Calculate ideal spatial distribution
      g[k] = N_rk[k]/(nbr_of_atoms*simul*N_ideal); //Calculate the radial distribution
      if (dr*k < 2.1)
	{
	  g[k] = 0; 
	}
    }
  /*_____________________________________*/
  
	
  

  /* data to output file */
  energy_file = fopen("energy.dat","w");
  for (int i = 0; i < nbr_of_timesteps + 1; i++) {
    current_time = i * timestep;
    fprintf(energy_file, "%.4f \t %e \t %e \t %e \n", current_time, E_pot[i], E_kin[i], E_pot[i]+E_kin[i]);
  }
  fclose(energy_file);

  temp_file = fopen("temp.dat","w"); // Temperature
  for (int i = 0; i < nbr_of_timesteps + 1; i++) {
    current_time = i * timestep;
    fprintf(temp_file, "%.4f \t %e\n", current_time, Temp[i]);
  }
  fclose(temp_file);

  press_file = fopen("press.dat","w"); // Pressure
  for (int i = 0; i < nbr_of_timesteps + 1; i++) {
    current_time = i * timestep;
    fprintf(press_file, "%.4f \t %e\n", current_time, Press[i]);
  }
  fclose(press_file);

  traj_file = fopen("traj.dat","w"); // Trajectory
  for (int i = 0; i < nbr_of_timesteps + 1; i++) {
    current_time = i * timestep;
    fprintf(press_file, "%.4f \t %e\n", current_time, Trajectory[i]);
  }
  fclose(traj_file);

 
  radial_file = fopen("radial.dat","w"); // Radial distribution function
  for (int i = 0; i < nbr_r_steps ; i++) {
    current_time = i * timestep;
    fprintf(press_file, "%.4f \t %e\n", dr*i, g[i]);
  }
  fclose(radial_file);
  
  //pos_file = fopen("pos.dat","w");   // Final positions, done only for equilibration run
  for (int i = 0; i < nbr_of_atoms ; i++) {
    for (int j = 0; j < 3; j++){
      // fprintf(pos_file, "%.4f \t", positions[i][j]);
    }
    //fprintf(pos_file, "\n");
  }
  //fclose(pos_file);

  //vel_file = fopen("vel.dat","w");  // Final velocities, done only for equilibration run
  for (int i = 0; i < nbr_of_atoms ; i++) {
    for (int j = 0; j < 3; j++){
      //fprintf(pos_file, "%.4f \t", v[i][j]);
    }
    //fprintf(vel_file, "\n");
  }
  //fclose(vel_file); /* Free allocated memory */
  /*________________________________________________________*/

  /* Free allocated memory */
  gsl_rng_free(q); 
  free(positions); 
  free(v);
  free(F);
  free(E_pot);
  free(E_kin);
  free(Temp);
  free(Press);
  free(Trajectory);
  free(sx); free(sy); free(sz);

  return 0;
}
