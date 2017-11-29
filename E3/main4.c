#include <stdio.h>
#include <stdlib.h>

int main()
{
  int i, N;
  FILE *in_file;
  
  N  = 1e6; /* The number of lines in MC.txt. */
  double *data = malloc((N) * sizeof (double));
  
  /* Read data from file. */
  in_file = fopen("MC.txt","r");
  for (i=0; i<N; i++) {
    fscanf(in_file,"%lf",&data[i]);
  }
  fclose(in_file);

  /* Variables and FILES */
  int j, k;
  int k_span = 300;
  double mean_f=0;
  double mean_sq_f=0;
  double f_ik=0; 
  double Phi[k_span];
  int B, b;
  int j_span;
  int n_B = 500;
  double s[n_B];
  double mean_F=0;
  double mean_sq_F=0;

  FILE *corrfile = fopen("phi.dat","w");
  FILE *blockfile = fopen("S.dat","w");

  for(i = 0; i < N; i++) {
    mean_f += data[i]/N;
    mean_sq_f += data[i]*data[i]/N; 
  }

  for(k = 0; k < k_span; k++) {
    for(i = 0; i < N-k; i++) {
      f_ik += data[i+k]*data[i];
    }
    f_ik = f_ik/(N-k);
    Phi[k] = (f_ik - mean_f*mean_f)/(mean_sq_f - mean_f*mean_f);
    fprintf(corrfile, "%i \t %.6f \n", k, Phi[k]);
  }
  fclose(corrfile);
  /*_______________________________________________*/
  
  for (b = 1; b < n_B; b++){
    B = b*10; 
    j_span = (int)N/B;
    printf("%d\n", j_span);
    double *F = malloc(j_span * sizeof (double));
    for (j = 0; j < j_span; j++){
      for (i = 0; i < B; i++){
	F[j] += data[j*B + i];
      }
      F[j] = F[j]/B;
      mean_F += F[j]/j_span;
      mean_sq_F += F[j]*F[j]/j_span;
    }
    printf("%d\n", j_span);
    s[b] = B*(mean_sq_F - mean_F*mean_F)/(mean_sq_f - mean_f*mean_f);
    fprintf(blockfile, "%i \t %.6f \n", B, s[b]);
    free(F);
    mean_F = 0;
    mean_sq_F = 0;
  }
  fclose(blockfile);
  free(data);
  
}
