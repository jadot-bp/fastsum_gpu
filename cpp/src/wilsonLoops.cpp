#include <cmath>  // pow
#include <complex>
// #include <stdlib.h>  // ???
#include <iostream>

#include "gauge.h"

using namespace std;

int* allocatePaths(int N, int M) {
  // Allocate memory to hold the paths
  // Return pointer to first element
  int* paths = new int[N * M];  //(int*)malloc(N * M * sizeof(int));
  if (paths == NULL) {
    cout << "Memory not allocated." << endl;
    exit(0);
  }
  return paths;
}

void one_x_one_Paths(int* fPath, int* bPath) {
  // Put the 6 paths for the whole plaquette
  // into the path variable
  // Space-space paths
  // mu=1, nu = 2
  int f12[] = {1, 2};
  int b21[] = {2, 1};
  // mu=1, nu = 3
  int f13[] = {1, 3};
  int b31[] = {3, 1};
  // mu=2, nu = 3
  int f23[] = {2, 3};
  int b32[] = {3, 2};
  // temporal paths mu = 0
  // nu = 1
  int f01[] = {0, 1};
  int b10[] = {1, 0};
  // nu = 2
  int f02[] = {0, 2};
  int b20[] = {2, 0};
  // nu = 3
  int f03[] = {0, 3};
  int b30[] = {3, 0};
  // Now put it into paths
  for (int ii = 0; ii < 2; ii++) {
    // spatial
    fPath[ii] = f12[ii];
    fPath[ii + 1 * 2] = f13[ii];
    fPath[ii + 2 * 2] = f23[ii];
    // temporal
    fPath[ii + 3 * 2] = f01[ii];
    fPath[ii + 4 * 2] = f02[ii];
    fPath[ii + 5 * 2] = f03[ii];
    // backwards
    bPath[ii] = b21[ii];
    bPath[ii + 1 * 2] = b31[ii];
    bPath[ii + 2 * 2] = b32[ii];
    // temporal
    bPath[ii + 3 * 2] = b10[ii];
    bPath[ii + 4 * 2] = b20[ii];
    bPath[ii + 5 * 2] = b30[ii];
  }
}

int calculatePathSize(double maxR, double maxT) {
  /*
    Calculates the number of steps or links
    in the forward (backward) path
    hence half the total number of steps
    Only to do for one mu, as the are all the same
  */
  int mu;
  double t0, t1;
  int stepCount = 0;

  for (mu = 1; mu < 2; mu++) {
    double x[3] = {0.0, 0.0, 0.0};
    double y[3] = {0.0, 0.0, 0.0};
    stepCount = 0;
    while (sqrt(pow(x[0] - y[0], 2) + pow(x[1] - y[1], 2) +
                pow(x[2] - y[2], 2)) < maxR) {
      stepCount++;
      y[mu - 1] = y[mu - 1] + 1.0;
    }

    t0 = 0.0;
    t1 = 0.0;
    while (fabs(t1 - t0) < maxT) {
      stepCount++;
      t1 = t1 + 1.0;
    }
  }
  return stepCount;
}

void wRT(double maxR, double maxT, int* fPaths, int* bPaths, int pathSize) {
  /*
    Constructs the paths for loops of size maxR * maxT
    // Note that FF corresponds to backwards loops that get conjugated
    // and BF corresponds to forward loops that don't
    // even though that's not what the forward/backward below say
   */
  int mu;
  int* FF;
  int stepCount;
  int BStepCount;
  double t0, t1;
  int* BF;
  double x[3] = {0.0, 0.0, 0.0};
  double y[3];
  // FF = (int*)malloc(sizeof(int) * pathSize);
  FF = new int[pathSize];
  // BF = (int*)malloc(sizeof(int) * pathSize);
  BF = new int[pathSize];
  // Loop over each direction
  // mu = 1
  mu = 1;
  y[0] = 0.0;
  y[1] = 0.0;
  y[2] = 0.0;
  stepCount = 0;
  // Here we set the mu into the list for the forward path
  while (sqrt(pow(x[0] - y[0], 2) + pow(x[1] - y[1], 2) + pow(x[2] - y[2], 2)) <
         maxR) {
    FF[stepCount] = mu;
    stepCount++;
    y[mu - 1] = y[mu - 1] + 1.0;
  }
  // Similarly for the T steps for this forward and backwards path
  BStepCount = 0;
  t0 = 0.0;
  t1 = 0.0;
  while (fabs(t1 - t0) < maxT) {
    FF[stepCount] = 0;
    BF[BStepCount] = 0;
    stepCount++;
    BStepCount++;
    t1 = t1 + 1.0;
  }
  y[0] = 0.0;
  y[1] = 0.0;
  y[2] = 0.0;
  while (sqrt(pow(x[0] - y[0], 2) + pow(x[1] - y[1], 2) + pow(x[2] - y[2], 2)) <
         maxR) {
    // Here we set the mu into the list for the backward path
    BF[BStepCount] = mu;
    BStepCount++;
    y[mu - 1] = y[mu - 1] + 1.0;
  }
  // printf("mu = %d \n", mu);
  // Now put them into the overall list of paths
  for (int ii = 0; ii < pathSize; ii++) {
    bPaths[ii + 0 * pathSize] = FF[ii];
    // printf("ii %d ind %d FF %d \n", ii, ii + 0 * pathSize, FF[ii]);
  }
  for (int ii = 0; ii < pathSize; ii++) {
    fPaths[ii + 0 * pathSize] = BF[ii];
    // printf("ii %d ind %d BF %d \n", ii, ii + 0 * pathSize, BF[ii]);
  }

  // mu = 2
  mu = 2;
  y[0] = 0.0;
  y[1] = 0.0;
  y[2] = 0.0;
  stepCount = 0;
  // Here we set the mu into the list for the forward path
  while (sqrt(pow(x[0] - y[0], 2) + pow(x[1] - y[1], 2) + pow(x[2] - y[2], 2)) <
         maxR) {
    FF[stepCount] = mu;
    stepCount++;
    y[mu - 1] = y[mu - 1] + 1.0;
  }
  // Similarly for the T steps for this forward and backwards path
  BStepCount = 0;
  t0 = 0.0;
  t1 = 0.0;
  while (fabs(t1 - t0) < maxT) {
    FF[stepCount] = 0;
    BF[BStepCount] = 0;
    stepCount++;
    BStepCount++;
    t1 = t1 + 1.0;
  }
  y[0] = 0.0;
  y[1] = 0.0;
  y[2] = 0.0;
  while (sqrt(pow(x[0] - y[0], 2) + pow(x[1] - y[1], 2) + pow(x[2] - y[2], 2)) <
         maxR) {
    // Here we set the mu into the list for the backward path
    BF[BStepCount] = mu;
    BStepCount++;
    y[mu - 1] = y[mu - 1] + 1.0;
  }
  // Now put them into the overall list of paths
  // printf("mu = %d \n", mu);
  for (int ii = 0; ii < pathSize; ii++) {
    bPaths[ii + 1 * pathSize] = FF[ii];
    // printf("ii %d ind %d FF %d \n", ii, ii + 1 * pathSize, FF[ii]);
  }
  for (int ii = 0; ii < pathSize; ii++) {
    fPaths[ii + 1 * pathSize] = BF[ii];
    // printf("ii %d ind %d BF %d \n", ii, ii + 1 * pathSize, BF[ii]);
  }

  // mu = 3
  mu = 3;
  y[0] = 0.0;
  y[1] = 0.0;
  y[2] = 0.0;
  stepCount = 0;
  // Here we set the mu into the list for the forward path
  while (sqrt(pow(x[0] - y[0], 2) + pow(x[1] - y[1], 2) + pow(x[2] - y[2], 2)) <
         maxR) {
    FF[stepCount] = mu;
    stepCount++;
    y[mu - 1] = y[mu - 1] + 1.0;
  }
  // Similarly for the T steps for this forward and backwards path
  BStepCount = 0;
  t0 = 0.0;
  t1 = 0.0;
  while (fabs(t1 - t0) < maxT) {
    FF[stepCount] = 0;
    BF[BStepCount] = 0;
    stepCount++;
    BStepCount++;
    t1 = t1 + 1.0;
  }
  y[0] = 0.0;
  y[1] = 0.0;
  y[2] = 0.0;
  while (sqrt(pow(x[0] - y[0], 2) + pow(x[1] - y[1], 2) + pow(x[2] - y[2], 2)) <
         maxR) {
    // Here we set the mu into the list for the backward path
    BF[BStepCount] = mu;
    BStepCount++;
    y[mu - 1] = y[mu - 1] + 1.0;
  }
  // Now put them into the overall list of paths
  // printf("mu = %d \n", mu);
  for (int ii = 0; ii < pathSize; ii++) {
    bPaths[ii + 2 * pathSize] = FF[ii];
    // printf("ii %d ind %d FF %d \n", ii, ii + 2 * pathSize, FF[ii]);
  }
  for (int ii = 0; ii < pathSize; ii++) {
    fPaths[ii + 2 * pathSize] = BF[ii];
    // printf("ii %d ind %d BF %d \n", ii, ii + 2 * pathSize, BF[ii]);
  }

  for (int pp = 0; pp < 3 * pathSize; pp += pathSize) {
    // printf("Path %d: \n", pp / pathSize + 1);
    for (int jj = 0; jj < pathSize; jj++) {
      // printf("link %d: \n", jj + 1);
      // printf("pp %d, jj %d, F %d, B %d \n", pp, jj, fPaths[pp + jj],
      // bPaths[pp + jj]);
    }
  }
  //  exit(1);
}

#pragma omp declare target
float tracePathWilsonLoop(wc U[], int pos[], int Nt, int Ns, int NPath,
                          int pathLength, int* fPath, int* bPath) {
  // Follows the wilson loop path from a fPath and bPath
  // Starting at pos
  // Returns the sum of the trace of the 3x3 link SU(3) matrices
  // Output
  float trace = 0;
  // Working matrices for calculations
  wc U_fLeft[3][3];
  wc U_fRight[3][3];
  wc U_bLeft[3][3];
  wc U_bRight[3][3];
  wc plaq[3][3];
  int fPos[4];
  int bPos[4];
  // Shape of the lattice
  int U_shape[DIM] = {Nt, Ns, Ns, Ns, ND, NC, NC};
  // Loop over all of the paths
  for (int pp = 0; pp < NPath * pathLength; pp += pathLength) {
    // printf("Path %d: \n", pp / pathLength + 1);
    // Copy position
    for (int ii = 0; ii < 4; ii++) {
      fPos[ii] = pos[ii];
      bPos[ii] = pos[ii];
    }

    // Get starting position
    // Start SU(3) indices at zero
    // printf("FOR x, y, z, t, mu, %d %d %d %d %d\n", fPos[0], fPos[1], fPos[2],
    // fPos[3], fPath[pp + 0]);
    // printf("BAC x, y, z, t, mu, %d %d %d %d %d\n",
    // bPos[0], bPos[1], bPos[2], bPos[3], bPath[pp + 0]);
    int U_f_pos[DIM] = {fPos[0],       fPos[1], fPos[2], fPos[3],
                        fPath[pp + 0], 0,       0};
    int U_b_pos[DIM] = {bPos[0],       bPos[1], bPos[2], bPos[3],
                        bPath[pp + 0], 0,       0};
    // Calculate address in memory of U_f, U_b
    int U_f_idx = idx(U_f_pos, U_shape, DIM);
    int U_b_idx = idx(U_b_pos, U_shape, DIM);
    // Populate U_[f,b]Left with value
    construct_3x3(U_fLeft, U, U_f_idx);
    construct_3x3(U_bLeft, U, U_b_idx);
    // printf("3x3 %f\n", real(U_fLeft[1][1]));
    // printf("idx %d\n", U_f_idx);
    // printf("\n blah %d \n\n", fPath[pp + 0]);
    // Update position
    fPos[fPath[pp + 0]] = pos[fPath[pp + 0]] + 1;
    bPos[bPath[pp + 0]] = pos[bPath[pp + 0]] + 1;
    // Enforce periodicity
    if (fPath[pp + 0] == 0) {
      fPos[fPath[pp + 0]] %= Nt;
    } else {
      fPos[fPath[pp + 0]] %= Ns;
    }
    if (bPath[pp + 0] == 0) {
      bPos[bPath[pp + 0]] %= Nt;
    } else {
      bPos[bPath[pp + 0]] %= Ns;
    }
    // Loop over each step in the path
    // Start at 1 as have already got the first link
    for (int jj = 1; jj < pathLength; jj++) {
      // printf("link %d: \n", jj + 1);
      // Get the next link
      // Start SU(3) indices at zero
      // printf("FOR x, y, z, t, mu, %d %d %d %d %d\n", fPos[0], fPos[1],
      //     fPos[2], fPos[3], fPath[pp + jj]);
      // printf("BAC x, y, z, t, mu, %d %d %d %d %d\n", bPos[0], bPos[1],
      // bPos[2], bPos[3], bPath[pp + jj]);
      int U_f_pos[DIM] = {fPos[0],        fPos[1], fPos[2], fPos[3],
                          fPath[pp + jj], 0,       0};
      int U_b_pos[DIM] = {bPos[0],        bPos[1], bPos[2], bPos[3],
                          bPath[pp + jj], 0,       0};
      // Calculate address in memory of U_f, U_b
      int U_f_idx = idx(U_f_pos, U_shape, DIM);
      int U_b_idx = idx(U_b_pos, U_shape, DIM);
      // Populate U_[f,b]Right with value
      construct_3x3(U_fRight, U, U_f_idx);
      construct_3x3(U_bRight, U, U_b_idx);
      // Do the matrix multiplication and store in left
      // Format is out, left, right
      MultiplyMat(U_fLeft, U_fLeft, U_fRight);
      MultiplyMat(U_bLeft, U_bLeft, U_bRight);
      // Update position
      fPos[fPath[pp + jj]] = fPos[fPath[pp + jj]] + 1;
      bPos[fPath[pp + jj]] = bPos[fPath[pp + jj]] + 1;
      // Enforce periodicity
      if (fPath[pp + jj] == 0) {
        fPos[fPath[pp + jj]] %= Nt;
      } else {
        fPos[fPath[pp + jj]] %= Ns;
      }
      if (bPath[pp + jj] == 0) {
        bPos[bPath[pp + jj]] %= Nt;
      } else {
        bPos[bPath[pp + jj]] %= Ns;
      }
    }
    // Have finished the path
    // So now multiply forward * conj(back)
    ConjTranspose(U_bLeft);
    MultiplyMat(plaq, U_fLeft, U_bLeft);
    // printf("3x3 %f %f %f\n", real(plaq[1][1]), real(U_fLeft[1][1]),
    // real(U_bLeft[1][1]));
    for (int i = 0; i < NC; i++) {
      trace += real(plaq[i][i]);
      // printf("trace %f\n", trace);
    }
  }
  // return summed trace
  // printf("trace %f\n", trace);
  return trace;
}
#pragma omp end declare target
/*
double one_x_one(dc U[], int pos[], int mu, int nu, int Nt, int Ns) {
  // Calculate the value of the plaquette
  // in the mu-nu plane at position pos.

  int U_shape[DIM] = {Nt, Ns, Ns, Ns, ND, NC, NC};  // Shape of the lattice

  // double complex plaq[3][3]; // Plaquette matrix

  dc U_mu[3][3];  // Working matrices for calculations
  dc U_nu[3][3];
  dc U_mu_dag[3][3];
  dc U_nu_dag[3][3];
  dc plaq[3][3];

  // Get position of U_mu on the lattice
  int U_mu_pos[DIM] = {pos[0], pos[1], pos[2], pos[3],
                       mu,     0,      0};  // Start SU(3) indices at zero
  int U_mu_idx =
      idx(U_mu_pos, U_shape, DIM);  // Calculate address in memory of U_mu

  // Populate U_mu with value of U_mu
  construct_3x3(U_mu, U, U_mu_idx);

  // Get neighbour in the mu-direction
  int next_pos[ND];
  memcpy(next_pos, pos, sizeof(int) * 4);  // Explicitly copy pos
  next_pos[mu] += 1;
  if (mu == 0) {
    next_pos[mu] %= Nt;
  } else {
    next_pos[mu] %= Ns;
  }  // Enforce periodicity

  // Get position of U_nu(x+mu) on the lattice
  int U_nu_pos[DIM] = {
      next_pos[0], next_pos[1], next_pos[2], next_pos[3], nu, 0, 0};
  int U_nu_idx =
      idx(U_nu_pos, U_shape, DIM);  // Calculate address in memory of U_nu

  // Populate U_nu with value of U_nu(x+mu)
  construct_3x3(U_nu, U, U_nu_idx);

  // Calculate U_mu (plaq) * U_nu(x+mu) (work)
  MultiplyMat(plaq, U_mu, U_nu);

  // Get neighbour in the nu-direction
  memcpy(next_pos, pos, sizeof(int) * 4);
  next_pos[nu] += 1;
  if (nu == 0) {
    next_pos[nu] %= Nt;
  } else {
    next_pos[nu] %= Ns;
  }  // Enforce periodicity

  // Get position of U_mu^\dag(x+nu)
  int U_mu_dag_pos[DIM] = {
      next_pos[0], next_pos[1], next_pos[2], next_pos[3], mu, 0, 0};
  int U_mu_dag_idx = idx(U_mu_dag_pos, U_shape,
                         DIM);  // Calculate address in memory of U_mu^\dag

  // Populate U_mu_dag with value of U_mu^\dag(x+nu)
  construct_3x3(U_mu_dag, U, U_mu_dag_idx);
  ConjTranspose(U_mu_dag);

  // Get position of U_nu^\dag
  int U_nu_dag_pos[DIM] = {pos[0], pos[1], pos[2], pos[3], nu, 0, 0};
  int U_nu_dag_idx = idx(U_nu_dag_pos, U_shape,
                         DIM);  // Calculate address in memory of U_nu^\dag

  // Populate U_nu_dag with value of U_nu and transpose
  construct_3x3(U_nu_dag, U, U_nu_dag_idx);
  ConjTranspose(U_nu_dag);

  // Temporarily store U_mu^\dag(x+nu) * U_nu^\dag in U_mu

  MultiplyMat(U_mu, U_mu_dag, U_nu_dag);

  // Calculate full plaquette;
  MultiplyMat(plaq, plaq, U_mu);

  double trace = 0;
  for (int i = 0; i < NC; i++) {
    trace += real(plaq[i][i]);
  }
  return trace;
}
*/
