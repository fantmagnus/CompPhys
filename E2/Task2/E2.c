/*
 E2.c
 
 Created by Anders Lindman on 2014-11-04.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "func.h"
#define PI 3.141592653589
#define nbr_of_timesteps 25000 /* nbr_of_timesteps+1 = power of 2, for best speed */


/* Main program */
int main()
{
    int nbr_of_particles = 32 ;
    double factor;
    double trans_matrix[nbr_of_particles][nbr_of_particles];
    double timestep;
    double timestep_sq,current_time;
    double m;
    double kappa;
    double E_0;
    double sum_Q;
    double sum_P;
    double omega_k;
    
    /* declare file variable */
    FILE *modes_file;

    /* displacement, velocity and acceleration */
    double q[nbr_of_particles];
    double v[nbr_of_particles];
    double a[nbr_of_particles];
    double Q[nbr_of_particles];
    double P[nbr_of_particles];
    /* Allocating memory for large vectors */
    /* displacements for writing to file */
    double (*q_i)[nbr_of_particles] = malloc(sizeof (double[nbr_of_timesteps+1][nbr_of_particles]));
    double (*v_i)[nbr_of_particles] = malloc(sizeof (double[nbr_of_timesteps+1][nbr_of_particles]));
    double (*E)[nbr_of_particles] = malloc(sizeof (double[nbr_of_timesteps+1][nbr_of_particles]));

    /* Set variables */
    timestep = 0.1;
    m = 1.0;
    kappa = 1.0;
    timestep_sq = timestep * timestep;
    E_0 = nbr_of_particles;
    
    /* Construction of tranformation matrix */
    factor = 1 / ((double) nbr_of_particles + 1);
    for (int i = 0; i < nbr_of_particles; i++) {
        for (int j = 0; j < nbr_of_particles; j++) {
            trans_matrix[i][j] = sqrt(2 * factor) * sin((j + 1) * (i + 1) * PI * factor);
        }
    }
    /* Initial conditions */
    /* Set initial displacements and velocites */
    for (int j = 0; j < nbr_of_particles; j++) {
      q[j] = 0;
      v[j] = sqrt(2*E_0)*sqrt(2/(m*(nbr_of_particles+1)))*sin(j*PI/(nbr_of_particles+1));
      q_i[0][j] = 0;
      v_i[0][j] = v[j];
    }
    /* Initial acceleration and energy */
    calc_acc(a, q, m, kappa, nbr_of_particles);
    /* Transformation to normal modes Q from displacements q.  */
      sum_Q = 0;
      sum_P = 0;
      for (int j = 0; j < nbr_of_particles; j++){
	sum_Q += q[j] * trans_matrix[0][j];
	sum_P += v[j] * trans_matrix[0][j];
      }
      Q[0] = sum_Q;
      P[0] = sum_P;
      omega_k = 2 * sqrt(kappa/m);
      E[0][0] = 0.5 * (P[0]*P[0] + omega_k*omega_k*Q[0]*Q[0]);

    /*___________________________________________________________*/
    /* timesteps according to velocity Verlet algorithm */
    for (int t = 1; t < nbr_of_timesteps + 1; t++) {
      /* v(t+dt/2) */
      for (int j = 0; j < nbr_of_particles; j++) {
	v[j] += timestep * 0.5 * a[j];
      } 

      /* q(t+dt) */
      for (int j = 0; j < nbr_of_particles; j++) {
	q[j] += timestep * v[j];
	q_i[t][j] = q[j];
      }

      /* a(t+dt) */
      calc_acc(a, q, m, kappa, nbr_of_particles);

      /* v(t+dt) */
      for (int j = 0; j < nbr_of_particles; j++) {
	v[j] += timestep * 0.5 * a[j];
	v_i[t][j] = v[j];
      } 

      /* Transformation to normal modes Q from displacements q.  */
      for (int k = 0; k < nbr_of_particles; k++){
        sum_Q = 0;
	sum_P = 0;
        for (int j = 0; j < nbr_of_particles; j++){
	  sum_Q += q[j] * trans_matrix[k][j];
	  sum_P += v[j] * trans_matrix[k][j];
        }
        Q[k] = sum_Q;
	P[k] = sum_P;
	omega_k = 2 * sqrt(kappa/m) * sin(k * PI)/(2 * (nbr_of_particles+1));
	E[t][k] = 0.5 * (P[k]*P[k] + omega_k*omega_k*Q[k]*Q[k]);
      }      
    }    
    /*_____________________________________________________________________*/ 

    /* Print displacement data to output file */
    modes_file = fopen("modeenergies.dat","w");

    for (int i = 0; i < nbr_of_timesteps + 1; i++) {
      current_time = i * timestep;
      fprintf(modes_file, "%.4f \t", current_time);
	for (int j = 0; j < nbr_of_particles; j++) {
	  fprintf(modes_file, "%e \t", E[i][j]);	
      }
      fprintf(modes_file, "\n");
    }
    fclose(modes_file);

    /* Free allocated memory */
    free(q_i); 
    free(v_i);
    free(E);
    return 0;
}
