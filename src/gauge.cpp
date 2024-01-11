#include "gauge.h"  //

#include <complex.h>  // For complex numbers
#include <errno.h>    // the errno in opening gaugefield
#include <math.h>     // Not sure what this is for

#include <iostream>  // For FILE
// #include <stdlib.h>   // or this ?
#include <string.h>  // for strerror
using namespace std;
/*
# define DIM 7 // Rank of the gauge field array Nt x Ns^3 x Nd x Nc^2
# define NC 3  // NC = 3 only
# define ND 4  // ND = 4 only

void MultiplyMat(double complex MM[NC][NC], double complex left[NC][NC], double
complex right[NC][NC]);

void readGauge_C(int NS, int NT, const char* filename, double complex* U);
int idx(int pos[], int shape[], int Nd);

void construct_3x3(double complex M[3][3], double complex U[], int idx);
*/

void MultiplyMat(double complex MM[NC][NC], double complex left[NC][NC],
                 double complex right[NC][NC]) {
  // Multiply two matrices left and right, and return the product in MM.
  // Assumes matrices are 3x3 double complex.
  // Explicitly write out the maths for it
  //
  // Args:
  //
  // MM    : output 3x3 double matrix M = left * right
  // left  : output 3x3 double matrix
  // right : output 3x3 double matrix

  // construct work matrix to prevent overloading if A or B equals M
  double complex work[NC][NC];

  // !# first index
  work[0][0] = left[0][0] * right[0][0] + left[0][1] * right[1][0] +
               left[0][2] * right[2][0];
  work[1][0] = left[1][0] * right[0][0] + left[1][1] * right[1][0] +
               left[1][2] * right[2][0];
  work[2][0] = left[2][0] * right[0][0] + left[2][1] * right[1][0] +
               left[2][2] * right[2][0];
  // !# second index
  work[0][1] = left[0][0] * right[0][1] + left[0][1] * right[1][1] +
               left[0][2] * right[2][1];
  work[1][1] = left[1][0] * right[0][1] + left[1][1] * right[1][1] +
               left[1][2] * right[2][1];
  work[2][1] = left[2][0] * right[0][1] + left[2][1] * right[1][1] +
               left[2][2] * right[2][1];
  // !# third index
  work[0][2] = left[0][0] * right[0][2] + left[0][1] * right[1][2] +
               left[0][2] * right[2][2];
  work[1][2] = left[1][0] * right[0][2] + left[1][1] * right[1][2] +
               left[1][2] * right[2][2];
  work[2][2] = left[2][0] * right[0][2] + left[2][1] * right[1][2] +
               left[2][2] * right[2][2];
  /*
  for (int i=0; i < NC; i++) {
    for (int j=0; j < NC; j++) {
      work[i][j] = 0;
      for (int k=0; k < NC; k++) {
        work[i][j] += left[i][k]*right[k][j];
      }
    }
  }
  */
  // Populate M with result
  for (int i = 0; i < NC; i++) {
    for (int j = 0; j < NC; j++) {
      MM[i][j] = 0;
      MM[i][j] = work[i][j];
    }
  }
}

void readGauge_C(int NS, int NT, const char* filename, double complex* U) {
  // Reads a gaugefield in 'C' format from dumping numpy array from lyncs in
  // python Args:
  //
  // NS       : spatial dim
  // NT       : temporal dim
  // filename : full path to gaugefield
  // U        : memory for lattice data
  FILE* fptr = fopen(filename, "r");
  if (fptr == NULL) {
    fprintf(stderr, "Couldn't open %s: %s\n", filename, strerror(errno));
    exit(1);
  }
  int nmemb = NT * NS * NS * NS * ND * NC * NC;  // Lattice volume
  int size = sizeof(double complex);
  // Read the gaugefields into U
  fread(U, size, nmemb, fptr);
  fclose(fptr);
}

int idx(int pos[], int shape[], int Nd) {
  // Calculate the row-major address in  memory at
  // position pos for given shape.
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
}

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
}

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
}

/*
int main() {
  // does stuff i guess
  // will be removed later
  int NS = 24;
  int NT = 8;
  // Lattice volume
  int nmemb = NT*NS*NS*NS*ND*NC*NC;
  int size = sizeof(double complex);
  // Allocate memory to store the lattice
  double complex *U = (double complex*) malloc(nmemb * size);
  readGauge_C(NS, NT, "../conf/Gen2_8x24_gfAr0.C", U);
}
*/
