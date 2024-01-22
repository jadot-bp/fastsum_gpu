#include <complex>
#include <cstdio>
//#include <ctime>
#include <chrono>

#include<ranges>
#include<execution>
#include<algorithm>
#include<numeric>


#include "../src/gauge.h"
#include "../src/wilsonLoops.h"

using namespace std;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;

double pathSum(int NS, int NT, int RTPathSize, int* fPathRT, int* bPathRT, dc U[]){
  // Does the loop over all paths
  // We construct an enumerate of the cartesian product
  // The enumerate is required for the transform_reduce
  auto myRange = std::views::enumerate(std::views::cartesian_product(
								     std::ranges::views::iota(0,NT),
								     std::ranges::views::iota(0,NS),
								     std::ranges::views::iota(0,NS),
								     std::ranges::views::iota(0,NS)
								     ));


  // mySum holds all of the sum_{paths} at each  (t, x, y, z)
  std::vector<double> mySum(NT * NS * NS * NS);
  return std::transform_reduce(std::execution::par_unseq, myRange.begin(), myRange.end(), 0., std::plus{}, [NS, NT, RTPathSize, fPathRT, bPathRT, U, mySum = mySum.data()] (auto index){
      // a lot here
      // execution method for parallel or not, etc (seq, unseq, par, par_unseq)
      // start iter
      // end iter
      // start value for sum
      // what operation doing (sum)
      // capturing most variables by value (which is often a pointer anyway)
      // caputuring mySum by reference
      // index is a tuple of <enum, <t, i, j, k>>
      // Unwrap it in two stages
      auto [ii, tijk] = index;
      auto [t, i, j, k] = tijk;
      // assemble position
      int pos[4] = {t, i, j, k};
      //printf("%d
      // do work
      mySum[ii] += tracePathWilsonLoop(U, pos, NT, NS, 3, RTPathSize,
				      fPathRT, bPathRT);
     
      return mySum[ii];
   });
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
  readGauge_C(NS, NT, "/cosma5/data/do009/dc-bign2/conf/Gen2l_64x32n325.C", U);


  const auto start = high_resolution_clock::now();

  int* fPathAr[NS - 3];
  int* bPathAr[NS - 3];

  int RTPathSize[NS - 3];

  double sumReTrWRT[NS - 3] = {};  // init to zero


  for (int R=0; R< NS - 3; R++)
    {
      RTPathSize[R] = calculatePathSize(R + 1, 1);
      printf("The path length in %dx1 is %d\n", R+1, RTPathSize[R]);
      fPathAr[R] = allocatePaths(3, RTPathSize[R]);
      bPathAr[R] = allocatePaths(3, RTPathSize[R]);
      // Set up the paths
      wRT(R + 1, 1, fPathAr[R], bPathAr[R], RTPathSize[R]);
      double mySum = pathSum(NS, NT, RTPathSize[R], fPathAr[R], bPathAr[R], U);
      sumReTrWRT[R] = mySum;
    }


  const auto end = high_resolution_clock::now();

  printf("\n\n\n");
  for (int R=0; R<NS - 3; R++){
    printf("The path length in %dx1 is %d for %lf \n", R+1, RTPathSize[R], sumReTrWRT[R]);
  }
  printf("Total execution time: %lfs\n", duration_cast<duration<double>>(end - start).count());
  return 0;
}
