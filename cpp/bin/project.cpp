#include <complex>
#include <cstdio>
// #include <stdio.h>
#include <ctime>
// #include <time.h>

#include "../src/gauge.h"
#include "../src/wilsonLoops.h"

using namespace std;

int main() {
  // does stuff
  int NS = 32;
  int NT = 32;
  // Lattice volume
  int nmemb = NT * NS * NS * NS * ND * NC * NC;
  int size = sizeof(dc);
  // Allocate memory to store the lattice
  dc* U = new dc[nmemb * size];  //(dc*)malloc(nmemb * size);
  // printf("%s\n", "Successfully allocated U memory");
  //  readGauge_C(NS, NT, "conf/Gen2_8x24_gfAr0.C", U);
  // readGauge_C(NS, NT, "../conf/Gen2P_k278480_128x48_Tune_053n2.C", U);
  readGauge_C(NS, NT, "/cosma5/data/do009/dc-bign2/conf/Gen2l_32x32n954.C", U);

  // dc plaq[3][3];  // Plaquette matrix
  // Allocate the path arrays with 6 paths of length 2
  int* fPath = allocatePaths(6, 2);
  int* bPath = allocatePaths(6, 2);
  // Set up the paths
  one_x_one_Paths(fPath, bPath);

  int RTPathSize = calculatePathSize(20, 1);
  printf("%s %d \n", "The path length in 1x1 is ", RTPathSize);
  int* fPathRT = allocatePaths(3, RTPathSize);
  int* bPathRT = allocatePaths(3, RTPathSize);
  // Set up the paths
  wRT(20, 1, fPathRT, bPathRT, RTPathSize);

  int S_nP = 0;  // The number of spatial plaquettes
  int T_nP = 0;  // The number of temporal plaquettes

  double S_sumReTrP = 0;  // Sum of real trace of spatial plaquettes
  double T_sumReTrP = 0;  // Sum of real trace of temporal plaquettes

  double sumReTrWRT = 0;  // alternative cac

  clock_t start = clock();

  
  int pos[4];
  int t, i, j, k;
  int mu, nu;

  // Loop over all sites
#pragma omp parallel for shared(U, fPathRT, bPathRT) reduction(+:sumReTrWRT,T_sumReTrP,S_sumReTrP,T_nP,S_nP) private(pos,t,i,j,k,mu,nu) default(none), firstprivate(NT,NS,RTPathSize)
    for (t = 0; t < NT; t++) {
      for (i = 0; i < NS; i++) {
	for (j = 0; j < NS; j++) {
	  for (k = 0; k < NS; k++) {
	  pos[0]=t;
	  pos[1]=i;
	  pos[2]=j;
	  pos[3]=k;

          sumReTrWRT += tracePathWilsonLoop(U, pos, NT, NS, 3, RTPathSize,
                                            fPathRT, bPathRT);
          // Loop over temporal Lorentz indices
          mu = 0;
          for (nu = 1; nu < ND; nu++) {
            T_sumReTrP += one_x_one(U, pos, mu, nu, NT, NS);
            T_nP += 1;
          }

          // Loop over spatial Lorentz indices
          for (mu = 1; mu < ND; mu++) {
            for (nu = mu + 1; nu < ND; nu++) {
              S_sumReTrP += one_x_one(U, pos, mu, nu, NT, NS);
              S_nP += 1;
            }
          }
	
        }
      }
    }
  }

  clock_t end = clock();

  printf("%s %d \n", "The path length in 1x1 is ", RTPathSize);
  printf("\n\n WRT calc %lf \n", sumReTrWRT);
  printf("Type:\t\tsumReTrP:\tnP:\tAvg:\n");
  printf("whole\t\t%lf\t%d\t%lf\n", S_sumReTrP + T_sumReTrP, S_nP + T_nP,
         (S_sumReTrP + T_sumReTrP) / (S_nP + T_nP));
  printf("spatial\t\t%lf\t%d\t%lf\n", S_sumReTrP, S_nP, S_sumReTrP / S_nP);
  printf("temporal\t%lf\t%d\t%lf\n", T_sumReTrP, T_nP, T_sumReTrP / T_nP);
  printf("Total execution time: %fs\n", (float)(end - start) / CLOCKS_PER_SEC);
  // Free the path arrays
  free(fPath);
  free(bPath);
  return 0;
}
