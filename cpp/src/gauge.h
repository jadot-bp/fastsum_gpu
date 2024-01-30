// #include <stdlib.h>  // or this ?
#include <cmath>    // Not sure what this is for
#include <complex>  // For complex numbers
#include <cstdio>
// #include <stdio.h>   // For FILE
#include <cerrno>   // the errno in opening gaugefield
#include <cstring>  // for strerror
#include <iostream>
typedef std::complex<double> dc;

typedef std::complex<float> wc;

#define DIM 7  // Rank of the gauge field array Nt x Ns^3 x Nd x Nc^2
#define NC 3   // NC = 3 only
#define ND 4   // ND = 4 only

void MultiplyMat(wc MM[NC][NC], wc left[NC][NC],
                 wc right[NC][NC]);  // 3x3 only

void readGauge_C(int NS, int NT, const char* filename, dc* U);
int idx(int pos[], int shape[], int Nd);

void construct_3x3(wc M[3][3], wc U[], int idx);
void ConjTranspose(wc M[NC][NC]);
