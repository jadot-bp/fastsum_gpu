#include <omp.h>

#include <complex>
#include <cstdio>

#include "../src/gauge.h"
#include "../src/wilsonLoops.h"

using namespace std;

wc* convertComplexDoubleToFloat(std::complex<double>* U, int size) {
  wc* U_float = new wc[size];
  for (int i = 0; i < size; ++i) {
    U_float[i] = wc(U[i].real(), U[i].imag());
  }
  return U_float;
}

int main() {
  // does stuff
  int NS = 32;
  int NT = 64;
  // Lattice volume
  int nmemb = NT * NS * NS * NS * ND * NC * NC;
  int size = sizeof(dc);
  // Allocate memory to store the lattice
  dc* U = new dc[nmemb * size];  //(dc*)malloc(nmemb * size);
  // printf("%s\n", "Successfully allocated U memory");
  //  readGauge_C(NS, NT, "conf/Gen2_8x24_gfAr0.C", U);
  // readGauge_C(NS, NT, "../conf/Gen2P_k278480_128x48_Tune_053n2.C", U);
  readGauge_C(NS, NT, "../conf/Gen2l_64x32n74.C", U);
  // readGauge_C(NS, NT, "../conf/Gen2_8x24_gfAr0.C", U);

  // convert to single precision
  wc* Uf = convertComplexDoubleToFloat(U, nmemb);

  double start = omp_get_wtime();

  int* fPathAr[NS - 3];
  int* bPathAr[NS - 3];

  int RTPathSize[NS - 3];

  float sumReTrWRT[NS - 3];
  // init to zero

#pragma omp parallel default(none) \
    shared(sumReTrWRT, Uf, fPathAr, bPathAr, RTPathSize) firstprivate(NS, NT)
#pragma omp single
  {
#pragma omp taskloop
    for (int R = 0; R < NS - 3; R++) {
      // sumReTrWRT[R] = 0.0;
      RTPathSize[R] = calculatePathSize(R + 1, 1);
      printf("The path length in %dx1 is %d\n", R + 1, RTPathSize[R]);
      fPathAr[R] = allocatePaths(3, RTPathSize[R]);
      bPathAr[R] = allocatePaths(3, RTPathSize[R]);
      // Set up the paths
      wRT(R + 1, 1, fPathAr[R], bPathAr[R], RTPathSize[R]);
    }
  }

#pragma omp taskwait
  printf("Have finished making all the paths \n");
#pragma omp single
  {
    for (int R = 0; R < NS - 3; ++R) {
#pragma omp task shared(sumReTrWRT, Uf, fPathAr, bPathAr, RTPathSize) \
    firstprivate(NS, NT, R)
      {
        // Loop over all sites
        float result = 0.0;
        int pos[4];
        int t, i, j, k;
#pragma omp target teams distribute parallel for reduction(+:result) collapse(4) private(pos) map(to:Uf[:nmemb],fPathAr[R][:3 * RTPathSize[R]],bPathAr[R][:3 * RTPathSize[R]],RTPathSize[R],NS,NT,R)
        for (t = 0; t < NT; t++) {
          for (i = 0; i < NS; i++) {
            for (j = 0; j < NS; j++) {
              for (k = 0; k < NS; k++) {
                pos[0] = t;
                pos[1] = i;
                pos[2] = j;
                pos[3] = k;
                result += tracePathWilsonLoop(Uf, pos, NT, NS, 3, RTPathSize[R],
                                              fPathAr[R], bPathAr[R]);
                // printf("%d %d %d %d %f\n", t, i, j, k, result);
              }
            }
          }
        }
        sumReTrWRT[R] = result;
        // printf("R %d %f %f\n", R, sumReTrWRT[R], result);
      }
    }
  }
#pragma omp taskwait
  {
    double end = omp_get_wtime();
    printf("\n\n\n");
    for (int R = 0; R < NS - 3; R++) {
      printf("The path length in %dx1 is %d for %f \n", R + 1, RTPathSize[R],
             sumReTrWRT[R]);
    }
    printf("Total execution time: %lfs\n", (float)(end - start));
  }
  return 0;
}
