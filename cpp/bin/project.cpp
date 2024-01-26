#include <complex>
#include <cstdio>
// #include <stdio.h>
#include <ctime>
// #include <time.h>
#include <sycl/sycl.hpp>

#include "../src/gauge.h"
#include "../src/wilsonLoops.h"

using namespace std;

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

double PtracePathWilsonLoop(dc U[], int pos[], int Nt, int Ns, int NPath,
                            int pathLength, int* fPath, int* bPath) {
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
    // fPos[3], fPath[pp + 0]); printf("BAC x, y, z, t, mu, %d %d %d %d %d\n",
    // bPos[0], bPos[1], bPos[2], bPos[3], bPath[pp + 0]);
    int U_f_pos[DIM] = {fPos[0],       fPos[1], fPos[2], fPos[3],
                        fPath[pp + 0], 0,       0};
    int U_b_pos[DIM] = {bPos[0],       bPos[1], bPos[2], bPos[3],
                        bPath[pp + 0], 0,       0};
    // Calculate address in memory of U_f, U_b
    int U_f_Pidx = Pidx(U_f_pos, U_shape, DIM);
    int U_b_Pidx = Pidx(U_b_pos, U_shape, DIM);
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
      int U_f_pos[DIM] = {fPos[0],        fPos[1], fPos[2], fPos[3],
                          fPath[pp + jj], 0,       0};
      int U_b_pos[DIM] = {bPos[0],        bPos[1], bPos[2], bPos[3],
                          bPath[pp + jj], 0,       0};
      // Calculate address in memory of U_f, U_b
      int U_f_Pidx = Pidx(U_f_pos, U_shape, DIM);
      int U_b_Pidx = Pidx(U_b_pos, U_shape, DIM);
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

int main() {
  // does stuff
  int NS = 24;
  int NT = 8;
  // Lattice volume
  int nmemb = NT * NS * NS * NS * ND * NC * NC;
  int size = sizeof(dc);
  // Allocate memory to store the lattice
  // dc* U = new dc[nmemb * size];  //(dc*)malloc(nmemb * size);

  // create a device queue
  // sycl::queue q{sycl::cpu_selector{}};

  sycl::device cpu_device;
  for (const auto& dev : sycl::device::get_devices()) {
    if (dev.is_cpu()) {
      cpu_device = dev;
      break;
    }
  }

  // dc* U = static_cast<dc*>(malloc_host(nmemb * size, q));
  // dc* U = static_cast<dc*>(malloc_shared(nmemb * size, q));
  dc* U = new dc[nmemb * size];

  // printf("%s\n", "Successfully allocated U memory");
  //  readGauge_C(NS, NT, "conf/Gen2_8x24_gfAr0.C", U);
  // readGauge_C(NS, NT, "../conf/Gen2P_k278480_128x48_Tune_053n2.C", U);
  readGauge_C(NS, NT, "../conf/Gen2_8x24_gfAr0.C", U);
  // q.wait();
  dc plaq[3][3];  // Plaquette matrix
  // Allocate the path arrays with 6 paths of length 2
  int* fPath = allocatePaths(6, 2);
  int* bPath = allocatePaths(6, 2);
  // Set up the paths
  one_x_one_Paths(fPath, bPath);

  int RTPathSize = calculatePathSize(1, 1);
  printf("%s %d \n", "The path length in 1x1 is ", RTPathSize);
  int* fPathRT = allocatePaths(3, RTPathSize);
  int* bPathRT = allocatePaths(3, RTPathSize);
  // int* paths = new int[N * M];  //(int*)malloc(N * M * sizeof(int));
  // int* fPathRT = static_cast<int*>(malloc_shared(3 * RTPathSize *
  // sizeof(int), q)); int* bPathRT = static_cast<int*>(malloc_shared(3 *
  // RTPathSize * sizeof(int), q));
  // Set up the paths
  wRT(1, 1, fPathRT, bPathRT, RTPathSize);

  int S_nP = 0;  // The number of spatial plaquettes
  int T_nP = 0;  // The number of temporal plaquettes

  double S_sumReTrP = 0;  // Sum of real trace of spatial plaquettes
  double T_sumReTrP = 0;  // Sum of real trace of temporal plaquettes

  double sumReTrP = 0;  // alternative cac

  double sumReTrWRT = 0;  // alternative cac

  clock_t start = clock();

  // sycl::queue q{sycl::cpu_selector_v(0)};
  // sycl::queue q;
  sycl::queue q(cpu_device);

  std::cout << "Device: " << q.get_device().get_info<sycl::info::device::name>()
            << "\n";

  /*
  dc* U_device = sycl::malloc_device<dc>(nmemb * size, q);
  printf("HERE1\n");
  q.memcpy(U_device, U, nmemb * size);
  //q.wait();
  printf("HERE2\n");
  */

  double sum = 0;
  {
    // create buffers for data and sum
    sycl::buffer buf_sum(&sum, sycl::range(1));

    // Create buffers for the data
    sycl::buffer<dc> U_buf(U, sycl::range<1>(sizeof(U) / sizeof(dc)));
    //    sycl::buffer<int> pos_buf(pos, sycl::range<1>(sizeof(pos) /
    //    sizeof(int)));
    sycl::buffer<int> fPath_buf(fPathRT,
                                sycl::range<1>(sizeof(fPathRT) / sizeof(int)));
    sycl::buffer<int> bPath_buf(bPathRT,
                                sycl::range<1>(sizeof(bPathRT) / sizeof(int)));

    q.submit([&](sycl::handler& h) {
      sycl::stream out(1024, 256, h);

      auto U_acc = U_buf.get_access<sycl::access::mode::read_write>(h);
      // auto pos_acc = pos_buf.get_access<sycl::access::mode::read>(h);
      auto fPath_acc = fPath_buf.get_access<sycl::access::mode::read>(h);
      auto bPath_acc = bPath_buf.get_access<sycl::access::mode::read>(h);

      h.parallel_for(
          sycl::range<3>(NT, NS, NS), reduction(buf_sum, h, plus<>()),
          [=](sycl::item<3> item, auto& temp) {
            // auto i = item.get_id();
            int R1 = item.get_id()[0];
            int R2 = item.get_id()[1];
            int R3 = item.get_id()[2];
            // out << R1 << " " << R2 << " " << R3 << " " << sycl::endl;
            dc plaq[3][3];
            for (int k = 0; k < NS; k++) {
              int pos[4] = {R1, R2, R3, k};
              double blah = 5.0;
              blah = PtracePathWilsonLoop(U_acc.get_pointer(), pos, NT, NS, 6,
                                          2, fPath_acc.get_pointer(),
                                          bPath_acc.get_pointer());
              // one_x_one(U_device, pos, 2, 1, 64, 32, plaq);

              int U_shape[DIM] = {NT, NS, NS, NS,
                                  ND, NC, NC};  // Shape of the lattice
              // Get position of U_mu on the lattice
              int U_mu_pos[DIM] = {
                  pos[0], pos[1], pos[2], pos[3],
                  1,      0,      0};  // Start SU(3) indices at zero
              int U_mu_idx = Pidx(U_mu_pos, U_shape,
                                  DIM);  // Calculate address in memory of U_mu
              // Pconstruct_3x3(plaq, U, U_mu_idx);
              // double blah = real(U[U_mu_idx]);
              temp.combine(blah);
            }
          });
    });
  }

  printf("The summed number was %f", sum);

  /*

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
          sumReTrP += tracePathWilsonLoop(U, pos, NT, NS, 6, 2, fPath, bPath);
          printf("%s", "\n");
          sumReTrWRT += tracePathWilsonLoop(U, pos, NT, NS, 3, RTPathSize,
                                            fPathRT, bPathRT);

          // Loop over temporal Lorentz indices
          int mu = 0;

          for (int nu = 1; nu < ND; nu++) {
            one_x_one(U, pos, mu, nu, NT, NS, plaq);
            double trace = 0;
            for (int i = 0; i < NC; i++) {
              trace += real(plaq[i][i]);
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
                trace += real(plaq[i][i]);
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
  */
  return 0;
}
