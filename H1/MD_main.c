/*
 MD_main.c
 
 Created by Anders Lindman on 2013-10-31.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include "initfcc.h"
#include "alpotential.h"

#define PI 3.141592653589
#define nbr_of_timesteps 10000 /* nbr_of_timesteps+1 = power of 2, for best speed */
#define nbr_of_atoms 256

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

  /* declare file variable */
  FILE *energy_file;

  /* displacement, velocity and acceleration */
 	   
  /* Allocating memory for large vectors */
  double (*positions)[3] = malloc(sizeof(double[nbr_of_atoms][3]));
  double (*v)[3] = malloc(sizeof(double[nbr_of_atoms][3]));
  double (*F)[3] = malloc(sizeof(double[nbr_of_atoms][3]));
  double *E_pot = malloc(nbr_of_timesteps * sizeof(double));
  double *E_kin = malloc(nbr_of_timesteps * sizeof(double));

  

  /* Set variables */
  timestep = 0.001;
  timestep_sq = timestep * timestep;
  cell_length = 4.045;
  n_cell = 4;
  m = 27*1.0364*0.0001;
  tot_length  =  cell_length*n_cell;

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
  /* Set initial displacements and velocites */

  /* Calculate initial accelerations (forces)  based on initial displacements */
  get_forces_AL(F, positions, tot_length, nbr_of_atoms); 
  
  /* timesteps according to velocity Verlet algorithm */
  for (int i = 0; i < nbr_of_timesteps + 1; i++) {
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
  }
  
  /* Print displacement and energy data to output file */
  energy_file = fopen("energy.dat","w");
  for (int i = 0; i < nbr_of_timesteps + 1; i++) {
    current_time = i * timestep;
    fprintf(energy_file, "%.4f \t %e \t %e \t %e \n", current_time, E_pot[i], E_kin[i], E_pot[i]+E_kin[i]);
  }
  fclose(energy_file);
   
  /* Free allocated memory */
  gsl_rng_free(q); 
  free(positions); 
  free(v);
  free(E_pot);
  free(E_kin);
  free(F);
  /* 
     Function that generates a fcc lattice in units of [Å]. Nc is the number of 
     primitive cells in each direction and a0 is the lattice parameter. The
     positions of all the atoms are stored in pos which should be a matrix of the
     size N x 3, where N is the number of atoms. The first, second and third column
     correspond to the x,y and z coordinate respectively.
    */
    /*
     init_fcc(pos, Nc, a0);
    */
    
    /* 
     Function that calculates the potential energy in units of [eV]. pos should be
     a matrix containing the positions of all the atoms, L is the length of the 
     supercell and N is the number of atoms.
    */
    /*
     double energy;
     energy = get_energy_AL(pos, L, N);
    */
    
    /* 
     Function that calculates the virial in units of [eV]. pos should be a matrix
     containing the positions of all the atoms, L is the length of the supercell 
     and N is the number of atoms.
    */
    /*
     double virial;
     virial = get_virial_AL(pos, L, N);
    */
    
    /*
     Function that calculates the forces on all atoms in units of [eV/Å]. the 
     forces are stored in f which should be a matrix of size N x 3, where N is the
     number of atoms and column 1,2 and 3 correspond to the x,y and z component of
     the force resepctively . pos should be a matrix containing the positions of 
     all the atoms, L is the length of the supercell and N is the number of atoms.
    */
    /*
     get_forces_AL(f,pos, L, N);
    */    
}
