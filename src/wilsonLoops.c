#include <complex.h>
#include <stdio.h>
#include "gauge.h"
#include <time.h>


void one_x_one(double complex U[], int pos[], int mu, int nu, int Nt, int Ns, double complex plaq[NC][NC]){
    // Calculate the value of the plaquette
    // in the mu-nu plane at position pos.

    int U_shape[DIM] = {Nt, Ns, Ns, Ns, ND, NC, NC}; // Shape of the lattice
    
    //double complex plaq[3][3]; // Plaquette matrix
    
    double complex U_mu[3][3];    // Working matrices for calculations
    double complex U_nu[3][3];
    double complex U_mu_dag[3][3];
    double complex U_nu_dag[3][3];
    
    // Get position of U_mu on the lattice
    int U_mu_pos[DIM] = {pos[0], pos[1], pos[2], pos[3], mu, 0, 0}; // Start SU(3) indices at zero
    int U_mu_idx = idx(U_mu_pos, U_shape, DIM); // Calculate address in memory of U_mu    
    // Populate U_mu with value of U_mu
    construct_3x3(U_mu, U, U_mu_idx);
    // up to here is fine
    // Get neighbour in the mu-direction
    int next_pos[ND];
    memcpy(next_pos, pos, sizeof(int)*4);  // Explicitly copy pos
    next_pos[mu] += 1;
    if(mu == 0){ next_pos[mu] %= Nt; }else{ next_pos[mu] %= Ns; } // Enforce periodicity
    
    // Get position of U_nu(x+mu) on the lattice
    int U_nu_pos[DIM] = {next_pos[0], next_pos[1], next_pos[2], next_pos[3], nu, 0, 0};
    int U_nu_idx = idx(U_nu_pos, U_shape, DIM); // Calculate address in memory of U_nu
    // Populate U_nu with value of U_nu(x+mu)
    construct_3x3(U_nu, U, U_nu_idx);
    // Calculate U_mu (plaq) * U_nu(x+mu) (plaq)
    MultiplyMat(plaq, U_mu, U_nu);
    
    // Get neighbour in the nu-direction
    memcpy(next_pos, pos, sizeof(int)*4);
    next_pos[nu] += 1;
    if(nu == 0){ next_pos[nu] %= Nt; }else{ next_pos[nu] %= Ns; } // Enforce periodicity
    
    // Get position of U_mu^\dag(x+nu)
    int U_mu_dag_pos[DIM] = {next_pos[0], next_pos[1], next_pos[2], next_pos[3], mu, 0, 0};
    int U_mu_dag_idx = idx(U_mu_dag_pos, U_shape, DIM); // Calculate address in memory of U_mu^\dag

    // Populate U_mu_dag with value of U_mu^\dag(x+nu)
    construct_3x3(U_mu_dag, U, U_mu_dag_idx);
    ConjTranspose(U_mu_dag);
    
    // Get position of U_nu^\dag
    int U_nu_dag_pos[DIM] = {pos[0], pos[1], pos[2], pos[3], nu, 0, 0};
    int U_nu_dag_idx = idx(U_nu_dag_pos, U_shape, DIM); // Calculate address in memory of U_nu^\dag

    // Populate U_nu_dag with value of U_nu and transpose
    construct_3x3(U_nu_dag, U, U_mu_dag_idx);
    ConjTranspose(U_nu_dag);
    
    // Temporarily store U_mu^\dag(x+nu) * U_nu^\dag in U_mu
    
    MultiplyMat(U_mu, U_mu_dag, U_nu_dag);
    
    // Calculate full plaquette;
    MultiplyMat(plaq, plaq, U_mu);

      double trace = 0;
    
      for(int i=0; i<NC; i++){
	trace += U_mu[i][i];
      }
      //printf("%f\n", creal(trace));

    
    // Return the real part of the trace
    /*
        double trace = 0;
    
    for(int i=0; i<NC; i++){
        trace += plaq[i][i];
    }
    
    return creal(trace);
    */
    // return plaq;
}

int main(){
  // does stuff
  int NS = 24;
  int NT = 8;
  // Lattice volume
  int nmemb = NT*NS*NS*NS*ND*NC*NC; 
  int size = sizeof(double complex);
  // Allocate memory to store the lattice
  double complex *U = (double complex*) malloc(nmemb * size);
  readGauge_C(NS, NT, "../conf/Gen2_8x24_gfAr0.C", U);

  double complex plaq[3][3]; // Plaquette matrix

  int pos[4] = {1, 1, 1, 1};
  
  one_x_one(U, pos, 1, 2, NT, NS, plaq);

  double trace = 0;
    
  for(int i=0; i<NC; i++){
    trace += plaq[i][i];
  }
  printf("%f", creal(trace));


  int S_nP = 0; // The number of spatial plaquettes
  int T_nP = 0; // The number of temporal plaquettes
  
  double S_sumReTrP = 0; // Sum of real trace of spatial plaquettes
  double T_sumReTrP = 0; // Sum of real trace of temporal plaquettes
  
  clock_t start = clock();
  
  // Loop over all sites
  for(int t=0; t<NT; t++){
    for(int i=0; i<NS; i++){
      for(int j=0; j<NS; j++){
	for(int k=0; k<NS; k++){
	  
	  int pos[4] = {t, i, j, k};
          
	  // Loop over temporal Lorentz indices
	  int mu = 0;
          
	  for(int nu=1; nu<ND; nu++){
	    one_x_one(U, pos, mu, nu, NT, NS, plaq);
	    trace = 0;
	    for(int i=0; i<NC; i++){
	      trace += plaq[i][i];
	    }
	    T_sumReTrP += trace; // plaquette(U, pos, mu, nu, Nt, Ns);
	    T_nP += 1;
	  }
          
	  // Loop over spatial Lorentz indices
	  for(int mu=1; mu<ND; mu++){
	    for(int nu=mu+1; nu<ND; nu++){
	      one_x_one(U, pos, mu, nu, NT, NS, plaq);
	      trace = 0;
	      for(int i=0; i<NC; i++){
		trace += plaq[i][i];
	      }
	      S_sumReTrP += trace; //plaquette(U, pos, mu, nu, Nt, Ns);
	      S_nP += 1;
	    }
	  }
	}
      }
        }
  }

    clock_t end = clock();

    printf("Type:\t\tsumReTrP:\tnP:\tAvg:\n");
    printf("whole\t\t%lf\t%d\t%lf\n", S_sumReTrP+T_sumReTrP,S_nP+T_nP,(S_sumReTrP+T_sumReTrP)/(S_nP+T_nP));
    printf("spatial\t\t%lf\t%d\t%lf\n", S_sumReTrP,S_nP,S_sumReTrP/S_nP);
    printf("temporal\t%lf\t%d\t%lf\n", T_sumReTrP,T_nP,T_sumReTrP/T_nP);
    printf("Total execution time: %fs\n",(float)(end-start)/CLOCKS_PER_SEC);
  
}
