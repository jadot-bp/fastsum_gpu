// #include <stdlib.h>  // or this ?
#include <complex.h>  // For complex numbers
#include <math.h>     // Not sure what this is for
// #include <stdio.h>   // For FILE
#include <errno.h>   // the errno in opening gaugefield
#include <string.h>  // for strerror

#include <iostream>

#define DIM 7  // Rank of the gauge field array Nt x Ns^3 x Nd x Nc^2
#define NC 3   // NC = 3 only
#define ND 4   // ND = 4 only

void MultiplyMat(double complex MM[NC][NC], double complex left[NC][NC],
                 double complex right[NC][NC]);  // 3x3 only

void readGauge_C(int NS, int NT, const char* filename, double complex* U);
int idx(int pos[], int shape[], int Nd);

void construct_3x3(double complex M[3][3], double complex U[], int idx);
void ConjTranspose(double complex M[NC][NC]);
