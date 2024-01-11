#include <complex.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define ND 4
#define NC 3
#define DIM 7  // Rank of the gauge field array Nt x Ns^3 x Nd x Nc^2

double plaquette(double complex U[], int pos[], int mu, int nu, int Nt, int Ns);
void MultiplyMat(double complex M[NC][NC], double complex A[NC][NC],
                 double complex B[NC][NC]);
void ConjTranspose(double complex M[NC][NC]);
int idx(int pos[], int shape[], int Nd);
void construct_3x3(double complex M[3][3], double complex U[], int idx);

int main() {
  // Calculates the average value of the trace of the spatial
  // and temporal plaquettes.
  // Assumes Nd=4 and Nc=3

  int Nt = 64;
  int Ns = 32;

  // Open pointer to binary file containing data
  FILE *fptr = fopen("Gen2l_64x32n100.C", "r");
  // FILE *fptr = fopen("testfile", "r");

  if (fptr == NULL) {
    fprintf(stderr, "Couldn't open %s: %s\n", "Gen2l_64x32n100.C",
            strerror(errno));
    exit(1);
  }

  int nmemb = Nt * Ns * Ns * Ns * ND * NC * NC;  // Lattice volume
  int size = sizeof(double complex);

  // Allocate memory to store the lattice
  double complex *U = (double complex *)malloc(nmemb * size);

  // Read the gaugefields into U
  fread(U, size, nmemb, fptr);
  fclose(fptr);

  int S_nP = 0;  // The number of spatial plaquettes
  int T_nP = 0;  // The number of temporal plaquettes

  double S_sumReTrP = 0;  // Sum of real trace of spatial plaquettes
  double T_sumReTrP = 0;  // Sum of real trace of temporal plaquettes

  clock_t start = clock();

  // Loop over all sites
  for (int t = 0; t < Nt; t++) {
    for (int i = 0; i < Ns; i++) {
      for (int j = 0; j < Ns; j++) {
        for (int k = 0; k < Ns; k++) {
          int pos[4] = {t, i, j, k};

          // Loop over temporal Lorentz indices
          int mu = 0;

          for (int nu = 1; nu < ND; nu++) {
            T_sumReTrP += plaquette(U, pos, mu, nu, Nt, Ns);
            T_nP += 1;
          }

          // Loop over spatial Lorentz indices
          for (int mu = 1; mu < ND; mu++) {
            for (int nu = mu + 1; nu < ND; nu++) {
              S_sumReTrP += plaquette(U, pos, mu, nu, Nt, Ns);
              S_nP += 1;
            }
          }
        }
      }
    }
  }

  clock_t end = clock();

  printf("Type:\t\tsumReTrP:\tnP:\tAvg:\n");
  printf("whole\t\t%lf\t%d\t%lf\n", S_sumReTrP + T_sumReTrP, S_nP + T_nP,
         (S_sumReTrP + T_sumReTrP) / (S_nP + T_nP));
  printf("spatial\t\t%lf\t%d\t%lf\n", S_sumReTrP, S_nP, S_sumReTrP / S_nP);
  printf("temporal\t%lf\t%d\t%lf\n", T_sumReTrP, T_nP, T_sumReTrP / T_nP);
  printf("Total execution time: %fs\n", (float)(end - start) / CLOCKS_PER_SEC);

  return EXIT_SUCCESS;
};

double plaquette(double complex U[], int pos[], int mu, int nu, int Nt,
                 int Ns) {
  // Calculate the value of the real part of the trace of the plaquette
  // in the mu-nu plane at position pos.

  int U_shape[DIM] = {Nt, Ns, Ns, Ns, ND, NC, NC};  // Shape of the lattice

  double complex plaq[3][3];  // Plaquette matrix

  double complex U_mu[3][3];  // Working matrices for calculations
  double complex U_nu[3][3];
  double complex U_mu_dag[3][3];
  double complex U_nu_dag[3][3];

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

  // Return the real part of the trace

  double trace = 0;

  for (int i = 0; i < NC; i++) {
    trace += plaq[i][i];
  }

  return creal(trace);
}

void MultiplyMat(double complex M[NC][NC], double complex A[NC][NC],
                 double complex B[NC][NC]) {
  // Multiply two matrices A and B, and return the product in M.
  // Assumes matrices are 3x3 double complex.

  int N = 3;

  // construct work matrix to prevent overloading if A or B equals M

  double complex work[N][N];

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      work[i][j] = 0;
      for (int k = 0; k < N; k++) {
        work[i][j] += A[i][k] * B[k][j];
      }
    }
  }

  // Populate M with result
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      M[i][j] = work[i][j];
    }
  }
};

void ConjTranspose(double complex M[NC][NC]) {
  // Calculates the conjugate transpose of M.

  int N = 3;

  // Create working matrix to avoid overloading M

  double complex work[N][N];

  // Calculate the conjugate transpose and store in work
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      work[j][i] = conj(M[i][j]);
    }
  }

  // Return work as M
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      M[i][j] = work[i][j];
    }
  }
};

int idx(int pos[], int shape[], int Nd) {
  // Calculate the row-major address in  memory at position pos for given shape.
  //
  // Args:
  //
  // pos[]   : Position within the array
  // shape[] : Shape (dimensionality) of the array
  // Nd      : Number of dimensions of the array.

  int idx = 0;

  for (int k = 0; k < Nd; k++) {
    int prod_dim = 1;

    for (int l = k + 1; l < Nd; l++) {
      prod_dim *= shape[l];
    }
    idx += prod_dim * pos[k];
  }
  return idx;
};

void construct_3x3(double complex M[3][3], double complex U[], int idx) {
  // Constructs a 3x3 matrix M from address idx in memory U[]
  //
  // The SU(3) matrices are stored as a stream of 9 complex doubles
  // starting at address idx.

  for (int c1 = 0; c1 < NC; c1++) {
    for (int c2 = 0; c2 < NC; c2++) {
      M[c1][c2] = U[idx];

      idx += 1;  // Get next SU(3) index in memory
    }
  }
};
