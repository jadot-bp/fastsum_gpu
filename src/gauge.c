#include <stdlib.h>  // or this ?
#include <complex.h> // For complex numbers
#include <math.h>    // Not sure what this is for
#include <stdio.h>   // For FILE
#include <string.h>  // for strerror
#include <errno.h>   // the errno in opening gaugefield


# define DIM 7 // Rank of the gauge field array Nt x Ns^3 x Nd x Nc^2
# define NC 3  // NC = 3 only
# define ND 4  // ND = 4 only

void MultiplyMat(double complex MM[NC][NC], double complex left[NC][NC], double complex right[NC][NC]);

void readGauge_C(int NS, int NT, const char* filename, double complex* U);
int idx(int pos[], int shape[], int Nd);


void MultiplyMat(double complex MM[NC][NC], double complex left[NC][NC], double complex right[NC][NC]){
  // Multiply two matrices left and right, and return the product in MM.
  // Assumes matrices are 3x3 double complex.
  // Explicitly write out the maths for it
  //
  // Args:
  //
  // MM    : output 3x3 double matrix M = left * right
  // left  : output 3x3 double matrix
  // right : output 3x3 double matrix

  // !# first index
  MM[0][0] = left[0][0] * right[0][0] + left[0][1] * right[1][0] + left[0][2] * right[2][0];
  MM[1][0] = left[1][0] * right[0][0] + left[1][1] * right[1][0] + left[1][2] * right[2][0];
  MM[2][0] = left[2][0] * right[0][0] + left[2][1] * right[1][0] + left[2][2] * right[2][0];
  // !# second index
  MM[0][1] = left[0][0] * right[0][1] + left[0][1] * right[1][1] + left[0][2] * right[2][1];
  MM[1][1] = left[1][0] * right[0][1] + left[1][1] * right[1][1] + left[1][2] * right[2][1];
  MM[2][1] = left[2][0] * right[0][1] + left[2][1] * right[1][1] + left[2][2] * right[2][1];
  // !# third index
  MM[0][2] = left[0][0] * right[0][2] + left[0][1] * right[1][2] + left[0][2] * right[2][2];
  MM[1][2] = left[1][0] * right[0][2] + left[1][1] * right[1][2] + left[1][2] * right[2][2];
  MM[2][2] = left[2][0] * right[0][2] + left[2][1] * right[1][2] + left[2][2] * right[2][2];
}

void readGauge_C(int NS, int NT, const char* filename, double complex* U) {
  // Reads a gaugefield in 'C' format from dumping numpy array from lyncs in python
  // Args:
  //
  // NS       : spatial dim
  // NT       : temporal dim
  // filename : full path to gaugefield
  // U        : memory for lattice data
  FILE *fptr = fopen(filename,"r");
  if( fptr == NULL ) {
    fprintf(stderr, "Couldn't open %s: %s\n", filename, strerror(errno));
    exit(1);
  }
  int nmemb = NT*NS*NS*NS*ND*NC*NC; // Lattice volume
  int size = sizeof(double complex);
  // Read the gaugefields into U
  fread(U, size, nmemb, fptr);
  fclose(fptr);
}


int idx(int pos[], int shape[], int Nd){
    // Calculate the row-major address in  memory at position pos for given shape.
    //
    // Args:
    //
    // pos[]   : Position within the array
    // shape[] : Shape (dimensionality) of the array
    // Nd      : Number of dimensions of the array.
    int idx = 0;
    for(int k=0; k < Nd; k++){
        int prod_dim = 1;
        for(int l=k+1; l < Nd; l++){
            prod_dim *= shape[l];
        }
        idx += prod_dim*pos[k];
    }
    return idx;
};


int main(){
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
