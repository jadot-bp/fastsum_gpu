#include <cerrno>  // the errno in opening gaugefield
#include <complex>
#include <cstdio>    // C stdlib
#include <iostream>  // For FILE
#include <iostream>
#include <sycl/sycl.hpp>

//#include "../src/gauge.h"

typedef std::complex<double> dc;

const int NS = 24;
const int NT = 8;
const int ND = 4;  // Assuming ND is 4
const int NC = 3;  // Assuming NC is 3
const int nmemb = NT * NS * NS * NS * ND * NC * NC;
int size = sizeof(dc);
dc* U = new dc[nmemb * size];

const int NPath = 3;
const int pathLength = 2;
int* fPath = new int[NPath * pathLength];
int* bPath = new int[NPath * pathLength];

void one_x_one_Paths(int* fPa, int* bPa) {
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
    fPa[ii] = f12[ii];
    fPa[ii + 1 * 2] = f13[ii];
    fPa[ii + 2 * 2] = f23[ii];
    // temporal
    fPa[ii + 3 * 2] = f01[ii];
    fPa[ii + 4 * 2] = f02[ii];
    fPa[ii + 5 * 2] = f03[ii];
    // backwards
    bPa[ii] = b21[ii];
    bPa[ii + 1 * 2] = b31[ii];
    bPa[ii + 2 * 2] = b32[ii];
    // temporal
    bPa[ii + 3 * 2] = b10[ii];
    bPa[ii + 4 * 2] = b20[ii];
    bPa[ii + 5 * 2] = b30[ii];
  }
}

void PMultiplyMat(dc MM[NC][NC], dc left[NC][NC], dc right[NC][NC]) {
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
  dc work[NC][NC];

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

void PConjTranspose(dc M[NC][NC]) {
  // Calculates the conjugate transpose of M.

  int N = 3;

  // Create working matrix to avoid overloading M

  dc work[3][3];

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

int Pidx(int pos[], int shape[], int Nd) {
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

template <typename ComplexAccessor2D, typename ComplexAccessor1D>
void Pconstruct_3x3(ComplexAccessor2D M, ComplexAccessor1D U, int idx) {
  // Constructs a 3x3 matrix M from address idx in memory U[]
  // The SU(3) matrices are stored as a stream of 9 complex doubles
  // starting at address idx.

  for (int c1 = 0; c1 < NC; c1++) {
    for (int c2 = 0; c2 < NC; c2++) {
      M[c1][c2] = U[idx];
      idx += 1;  // Get next SU(3) index in memory
    }
  }
}

/*
void Pconstruct_3x3(dc M[3][3], dc U[], int idx) {
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
*/

// double PtracePathWilsonLoop(dc U[], int pos[], int Nt, int Ns, int NPath,
template <typename ComplexAccessor, typename IntAccessor>
double PtracePathWilsonLoop(ComplexAccessor U, int pos[], int Nt, int Ns,
                            int NPath, int pathLength, IntAccessor fPath,
                            IntAccessor bPath) {
  // Follows the wilson loop path from a fPath and bPath
  // Starting at pos
  // Returns the sum of the trace of the 3x3 link SU(3) matrices
  // Output
  double trace = 0;
  // Working matrices for calculations
  dc U_fLeft[3][3];
  dc U_fRight[3][3];
  dc U_bLeft[3][3];
  dc U_bRight[3][3];
  dc plaq[3][3];
  int fPos[4];
  int bPos[4];
  // Shape of the lattice
  int U_shape[7] = {Nt, Ns, Ns, Ns, ND, NC, NC};
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
    // fPos[3], fPath[pp + 0]); printf("BAC x, y, z, t, mu, %d %d %d %d %d\n",
    // bPos[0], bPos[1], bPos[2], bPos[3], bPath[pp + 0]);
    int U_f_pos[7] = {fPos[0], fPos[1], fPos[2], fPos[3], fPath[pp + 0], 0, 0};
    int U_b_pos[7] = {bPos[0], bPos[1], bPos[2], bPos[3], bPath[pp + 0], 0, 0};
    // Calculate address in memory of U_f, U_b
    int U_f_Pidx = Pidx(U_f_pos, U_shape, 7);
    int U_b_Pidx = Pidx(U_b_pos, U_shape, 7);
    // Populate U_[f,b]Left with value
    Pconstruct_3x3(U_fLeft, U, U_f_Pidx);
    Pconstruct_3x3(U_bLeft, U, U_b_Pidx);
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
      // fPos[2], fPos[3], fPath[pp + jj]); printf("BAC x, y, z, t, mu, %d %d %d
      // %d %d\n", bPos[0], bPos[1], bPos[2], bPos[3], bPath[pp + jj]);
      int U_f_pos[7] = {fPos[0],        fPos[1], fPos[2], fPos[3],
                        fPath[pp + jj], 0,       0};
      int U_b_pos[7] = {bPos[0],        bPos[1], bPos[2], bPos[3],
                        bPath[pp + jj], 0,       0};
      // Calculate address in memory of U_f, U_b
      int U_f_Pidx = Pidx(U_f_pos, U_shape, 7);
      int U_b_Pidx = Pidx(U_b_pos, U_shape, 7);
      // Populate U_[f,b]Right with value
      Pconstruct_3x3(U_fRight, U, U_f_Pidx);
      Pconstruct_3x3(U_bRight, U, U_b_Pidx);
      // Do the matrix multiplication and store in left
      // Format is out, left, right
      PMultiplyMat(U_fLeft, U_fLeft, U_fRight);
      PMultiplyMat(U_bLeft, U_bLeft, U_bRight);
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
    PConjTranspose(U_bLeft);
    PMultiplyMat(plaq, U_fLeft, U_bLeft);
    for (int i = 0; i < NC; i++) {
      trace += real(plaq[i][i]);
    }
  }
  // return summed trace
  return trace;
}

void iterateOverPos() {
  sycl::device myDevice;
  for (const auto& dev : sycl::device::get_devices()) {
    if (dev.is_cpu()) {
      myDevice = dev;
      break;
    }
  }

  sycl::queue queue(myDevice);

  std::cout << "Device: "
            << queue.get_device().get_info<sycl::info::device::name>() << "\n";

  sycl::range<3> dimensions(NS, NS, NS);
  double sum = 0;

  // Create buffers for the data
  {
    sycl::buffer<dc, 1> U_buf(U, sycl::range<1>(nmemb));
    sycl::buffer<int, 1> fPath_buf(fPath, sycl::range<1>(NPath * pathLength));
    sycl::buffer<int, 1> bPath_buf(bPath, sycl::range<1>(NPath * pathLength));
    sycl::buffer<double> buf_sum(&sum, sycl::range(1));

    queue.submit([&](sycl::handler& cgh) {
      auto U_acc = U_buf.get_access<sycl::access::mode::read>(cgh);
      auto fPath_acc = fPath_buf.get_access<sycl::access::mode::read>(cgh);
      auto bPath_acc = bPath_buf.get_access<sycl::access::mode::read>(cgh);

      auto sumReduction = reduction(buf_sum, cgh, std::plus<>());
      cgh.parallel_for(dimensions, sumReduction,
                       [=](sycl::item<3> item, auto& temp) {
                         int pos[4];
                         // pos[0] = idx[0];
                         // pos[1] = idx[1];
                         // pos[2] = idx[2];
                         pos[0] = 0;
                         pos[1] = 0;
                         pos[2] = 0;
                         for (pos[3] = 0; pos[3] < NT; pos[3]++) {
                           // Call PtracePathWilsonLoop with pos and U here
                           // double tempSum = PtracePathWilsonLoop(U_acc, pos,
                           // NT, NS, NPath, pathLength, fPath_acc, bPath_acc);
                           temp.combine(1.0);
                         }
                       });
    });
    // Ensure all operations are finished before we delete
    queue.wait_and_throw();
  };
  printf("Done %lf", sum);
}

int main() {
  FILE* fptr = fopen("../conf/Gen2_8x24_gfAr0.C", "r");
  // Read the gaugefields into U
  int result = fread(U, size, nmemb, fptr);
  int res = fclose(fptr);
  printf("have read gaugefield\n");
  // one_x_one_Paths(fPath, bPath);
  iterateOverPos();
  delete[] U;
  delete[] fPath;
  delete[] bPath;
  return 0;
}
