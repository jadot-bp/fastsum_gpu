#include <complex.h>

#include <cstdio>
// #include <stdio.h>
// #include <ctime>
#include <time.h>

#include "../src/gauge.h"
#include "../src/wilsonLoops.h"

using namespace std;

int main() {
  // does stuff
  int NS = 32;
  int NT = 64;
  // Lattice volume
  int nmemb = NT * NS * NS * NS * ND * NC * NC;
  int size = sizeof(double complex);
  // Allocate memory to store the lattice
  double complex* U = (double complex*)malloc(nmemb * size);
  // printf("%s\n", "Successfully allocated U memory");
  //  readGauge_C(NS, NT, "conf/Gen2_8x24_gfAr0.C", U);
  // readGauge_C(NS, NT, "../conf/Gen2P_k278480_128x48_Tune_053n2.C", U);
  readGauge_C(NS, NT, "../conf/Gen2l_64x32n100.C", U);

  double complex plaq[3][3];  // Plaquette matrix
  // Allocate the path arrays with 6 paths of length 2
  int* fPath = allocatePaths(6, 2);
  int* bPath = allocatePaths(6, 2);
  // Set up the paths
  one_x_one_Paths(fPath, bPath);

  int RTPathSize = calculatePathSize(1, 1);
  printf("%s %d \n", "The path length in 1x1 is ", RTPathSize);
  int* fPathRT = allocatePaths(3, RTPathSize);
  int* bPathRT = allocatePaths(3, RTPathSize);
  // Set up the paths
  wRT(1, 1, fPathRT, bPathRT, RTPathSize);

  int S_nP = 0;  // The number of spatial plaquettes
  int T_nP = 0;  // The number of temporal plaquettes

  double S_sumReTrP = 0;  // Sum of real trace of spatial plaquettes
  double T_sumReTrP = 0;  // Sum of real trace of temporal plaquettes

  double sumReTrP = 0;  // alternative cac

  double sumReTrWRT = 0;  // alternative cac

  clock_t start = clock();

  // Loop over all sites
  for (int t = 0; t < NT; t++) {
    for (int i = 0; i < NS; i++) {
      for (int j = 0; j < NS; j++) {
        for (int k = 0; k < NS; k++) {
          // for (int t = 0; t < 1; t++) {
          //   for (int i = 0; i < 1; i++) {
          //     for (int j = 0; j < 1; j++) {
          //       for (int k = 0; k < 1; k++) {
          int pos[4] = {t, i, j, k};

          // Alternative plaq calc
          sumReTrP +=
              tracePathWilsonLoop(U, pos, NT, NS, 6, 2, fPath, bPath, plaq);
          printf("%s", "\n");
          sumReTrWRT += tracePathWilsonLoop(U, pos, NT, NS, 3, RTPathSize,
                                            fPathRT, bPathRT, plaq);

          // Loop over temporal Lorentz indices
          int mu = 0;

          for (int nu = 1; nu < ND; nu++) {
            one_x_one(U, pos, mu, nu, NT, NS, plaq);
            double trace = 0;
            for (int i = 0; i < NC; i++) {
              trace += plaq[i][i];
            }
            T_sumReTrP += trace;  // plaquette(U, pos, mu, nu, Nt, Ns);
            T_nP += 1;
          }

          // Loop over spatial Lorentz indices
          for (int mu = 1; mu < ND; mu++) {
            for (int nu = mu + 1; nu < ND; nu++) {
              one_x_one(U, pos, mu, nu, NT, NS, plaq);
              double trace = 0;
              for (int i = 0; i < NC; i++) {
                trace += plaq[i][i];
              }
              S_sumReTrP += trace;  // plaquette(U, pos, mu, nu, Nt, Ns);
              S_nP += 1;
            }
          }
        }
      }
    }
  }

  clock_t end = clock();

  printf("\n\n alt plaq calc %lf \n", sumReTrP);
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
