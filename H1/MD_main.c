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
#include "func.h"

#define PI 3.141592653589
#define nbr_of_timesteps 1000 /* nbr_of_timesteps+1 = power of 2, for best speed */
#define nbr_of_atoms 256

/* Main program */
int main()
{ 
  // Declaration of variables
  double timestep;
  double timestep_sq,current_time;
  double m;
  double kappa;
  int n_cell; // Number of atoms in unit cell
  double cell_length;

  /* declare file variable */
  FILE *energy_file;

  /* displacement, velocity and acceleration */
  double disp_x[nbr_of_atoms];
  double disp_y[nbr_of_atoms];
  double disp_z[nbr_of_atoms];
  double v[nbr_of_atoms][3];	   
  double a_x[nbr_of_atoms];
  double a_y[nbr_of_atoms];
  double a_z[nbr_of_atoms];
  /* Allocating memory for large vectors */
  double (*positions)[3] = malloc(sizeof(double[nbr_of_atoms][3]));
  double (*positions_0)[3] = malloc(sizeof(double[nbr_of_atoms][3]));
  double (*a)[3] = malloc(sizeof(double[nbr_of_atoms][3]));
  double *E_pot = malloc(nbr_of_timesteps * sizeof(double));
  double *E_kin = malloc(nbr_of_timesteps * sizeof(double));

  

  /* Set variables */
  timestep = 0.0001;
  kappa = 99.86; /*Hopefully correct in the correct units */
  timestep_sq = timestep * timestep;
  cell_length = 4.03;
  n_cell = 4;
  m = 1;

  /* Initialize lattice*/
  init_fcc(positions, n_cell, cell_length);
  
  /* Calculate initial displacements */
  for (int i = 0; i < nbr_of_atoms; i++) {
    disp_x[i] = positions[i][0]-positions_0[i][0];
    disp_y[i] = positions[i][1]-positions_0[i][1];
    disp_z[i] = positions[i][2]-positions_0[i][2];
  }
  
  
  
  // Initialiasation for lattice and random deviation in positions
  init_fcc(positions, n_cell, cell_length); // Initialize lattice
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
	}
    }
  /* Set initial displacements and velocites */

  /* Calculate initial accelerations based on initial displacements */
  calc_acc(a_x, disp_x, m, kappa, nbr_of_atoms);
  calc_acc(a_y, disp_y, m, kappa, nbr_of_atoms);
  calc_acc(a_z, disp_z, m, kappa, nbr_of_atoms);

  
  /* Calculate initial energies */
  E_pot[0] = get_energy_AL(positions, cell_length, nbr_of_atoms);
  
  /* timesteps according to velocity Verlet algorithm */
  for (int i = 1; i < nbr_of_timesteps + 1; i++) {
    for (int x = 0; x < 3; x++) { // For three space directions 
      for (int j = 0; j < nbr_of_atoms; j++) {
	a[j][0] = a_x[j];
	a[j][1] = a_y[j];
	a[j][2] = a_z[j];
	v[j][x] += timestep * 0.5 * a[j][x]; // v(t+dt)
	positions[j][x] += timestep * v[j][x];  // q(t+dt)
      }
      calc_acc(a_x, disp_x, m, kappa, nbr_of_atoms); // a(t+dt)
      calc_acc(a_y, disp_y, m, kappa, nbr_of_atoms); 
      calc_acc(a_z, disp_z, m, kappa, nbr_of_atoms); 
      
      // a(t+dt)
      for(int j = 0; j < nbr_of_atoms; j++) {
	v[j][x] += timestep * 0.5 * a[j][x]; // v(t+dt)
	E_kin[i] += m * v[j][x] * v[j][x]*0.5; 
      }
      E_pot[i] = get_energy_AL(positions, cell_length, nbr_of_atoms);
    }
  }
  /* Print displacement and energy data to output file */
  energy_file = fopen("energy.dat","w");
  for (int i = 0; i < nbr_of_timesteps + 1; i++) {
    current_time = i * timestep;
    fprintf(energy_file, "%.4f \t %e \t %e \t %e \n", current_time, E_pot[i], E_kin[i], E_pot[i]+E_kin[i]);
  }
  fclose(energy_file);
   
  /* Free allocated memory */
  gsl_rng_free(q); /*  deallocate  rng */
  free(positions); 
  free(positions_0);
  free(E_pot);
  free(a);
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
