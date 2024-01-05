#include <complex.h>

#include "gauge.h"

void one_x_one(double complex U[], int pos[], int mu, int nu, int Nt, int Ns,
               double complex plaq[3][3]) {
  // Calculate the value of the plaquette
  // in the mu-nu plane at position pos.

  int U_shape[DIM] = {Nt, Ns, Ns, Ns, ND, NC, NC};  // Shape of the lattice

  // double complex plaq[3][3]; // Plaquette matrix

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
}
