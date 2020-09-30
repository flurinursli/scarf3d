/*
Copyright (c) 2020, Eidgenoessische Technische Hochschule Zurich, ETHZ.

Written by:
Walter Imperatori (walter.imperatori@sed.ethz.ch)

All rights reserved.

This file is part of SCARF3D, version: 2.4

SCARF3D is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SCARF3D is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SCARF3D.  If not, see <https://www.gnu.org/licenses/>.
*/

/*
 Purpose:
   To provide a sample C++ driver program for the SCARF3D library

 Revisions:
     Date                    Description of change
     ====                    =====================
   04/05/20                  original version
*/

#include <mpi.h>
#include <iostream>
#include <iomanip>
#include "scarf3d.h"

// FORTRAN subroutine prototype
extern "C"{
  extern void sample_mesh(const int* nd, const int* rank, const int* ntasks, const int n[], int fs[], int fe[]);
}

// timing functions
inline void watch_start(double* t){MPI_Barrier(MPI_COMM_WORLD); *t = MPI_Wtime();}
inline void watch_stop(double* t){MPI_Barrier(MPI_COMM_WORLD); *t = MPI_Wtime() - *t;}


int main(){

  // initialise MPI
  MPI_Init(nullptr, nullptr);

  // get number of tasks
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // get rank of each task
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (world_rank == 0){
    std::cout << ""                                                                    << std::endl;
    std::cout << "This program will compute small 2D and 3D random fields"             << std::endl;
    std::cout << "using the FIM (and SRM - if selected at compile-time) algorithm(s)." << std::endl;
    std::cout << "Test should complete in a few minutes. Output can be inspected"      << std::endl;
    std::cout << "with Matlab/Octave script 'scarf3d.m' and 'vizme.m'."                << std::endl;
    std::cout << ""                                                                    << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // Set parameters

  // number of grid points in model
  const int n[3] = {500, 500, 500};

  // grid-step
  const real ds = 50.;
  const real dh = 50.;

  // autocorrelation function (0=von karman/exponential, 1=gaussian)
  const int acf = 0;

  // correlation length
  const real cl[3] = {2000., 2000., 2000.};

  // standard deviation
  const real sigma = 0.05;

  // seed number
  const int seed = 1235;

  // hurst exponent
  const real hurst = 0.25;

  // number of POI
  const int npoi = 2;

  const real taper = 2500;
  const real mute  = 500;

  // rotation angles
  const real alpha = 30;
  const real beta  = -10;

  const int nwriters = 2;

  // other variables used throughout the code
  double tictoc;
  int fs[3], fe[3], dims[3];
  int npts, nd;
  long c;
  real stats[8];
  real *x, *y, *z, *field, *poi;

  // set formatting
  std::cout << std::right << std::fixed << std::boolalpha;
  std::cout.precision(5);

  // ================================================================================================================================
  // --------------------------------------------------------------------------------------------------------------------------------
  // create sample cartesian mesh: 2D case
  // --------------------------------------------------------------------------------------------------------------------------------
  // ================================================================================================================================

  nd = 2;

  sample_mesh(&nd, &world_rank, &world_size, n, fs, fe);

  for (int i = 0; i < nd; i++){
    dims[i] = (fe[i] - fs[i] + 1);
  }

  x     = new real[dims[0]*dims[1]];
  y     = new real[dims[0]*dims[1]];
  field = new real[dims[0]*dims[1]];

  for (int j = 0; j < dims[1]; j++){
    for (int i = 0; i < dims[0]; i++){
      c    = (j * dims[0]) + i;
      x[c] = (i + fs[0] - 1) * ds;
      y[c] = (j + fs[1] - 1) * ds;
    }
  }

  npts = dims[0] * dims[1];

  poi = new real[nd * npoi];

  poi[0] = 400 * ds;
  poi[1] = 250 * ds;
  poi[2] = 200 * ds;
  poi[3] = 150 * ds;

  // ================================================================================================================================
  // --------------------------------------------------------------------------------------------------------------------------------
  // tests FIM algorithm
  // --------------------------------------------------------------------------------------------------------------------------------
  // ================================================================================================================================

  if (world_rank == 0){
    std::cout << ""                                                  << std::endl;
    std::cout << "*************************************************" << std::endl;
    std::cout << "****** FIM algorithm, 2D, cartesian mesh *******" << std::endl;
  }

  {

    Scarf3D::options.hurst = hurst;
    Scarf3D::options.alpha = alpha;
    Scarf3D::options.npoi  = npoi;
    Scarf3D::options.poi   = poi;
    Scarf3D::options.taper = taper;
    Scarf3D::options.mute  = mute;

    Scarf3D::Initialize<fim> S(nd, fs, fe, ds, acf, cl, sigma);

    watch_start(&tictoc);
    S.execute(seed, field, stats);
    watch_stop(&tictoc);

    if (world_rank == 0){
      std::cout << ""                                                                                       << std::endl;
      std::cout << "Statistics for current simulation"                                                      << std::endl;
      std::cout << "*************************************************"                                      << std::endl;
      std::cout << "Elapsed time                  |" << std::setw(12) << (float) tictoc   << " sec" << " |" << std::endl;
      std::cout << "   + spectrum                 |" << std::setw(12) << (float) stats[4] << " sec" << " |" << std::endl;
      std::cout << "   + symmetry                 |" << std::setw(12) << (float) stats[5] << " sec" << " |" << std::endl;
      std::cout << "   + FFT                      |" << std::setw(12) << (float) stats[6] << " sec" << " |" << std::endl;
      std::cout << "   + interpolation            |" << std::setw(12) << (float) stats[7] << " sec" << " |" << std::endl;
      std::cout << "------------------------------|-----------------|"                                      << std::endl;
      std::cout << "Domain too small?             |" << std::setw(12) << (bool) stats[0]  << "     |"       << std::endl;
      std::cout << "Grid-step too large?          |" << std::setw(12) << (bool) stats[1]  << "     |"       << std::endl;
      std::cout << "------------------------------|-----------------|"                                      << std::endl;
      std::cout << "Discrete standard deviation   |" << std::setw(12) << (float) stats[2] << "     |"       << std::endl;
      std::cout << "Discrete mean value           |" << std::setw(12) << (float) stats[3] << "     |"       << std::endl;
      std::cout << "------------------------------|-----------------|"                                      << std::endl;
    }

    watch_start(&tictoc);
    S.io(nd, n, field, "fim_cart_whole_2d", &nwriters);
    watch_stop(&tictoc);

    if (world_rank == 0) {
      std::cout << "I/O time                      |" << std::setw(12) << (float) tictoc << " sec" << " |" << std::endl;
      std::cout << "*************************************************"                                    << std::endl;
    }

    // call destructor explicitly
    //S.~Initialize();
  }

  if (world_rank == 0){
    std::cout << ""                                                  << std::endl;
    std::cout << "*************************************************" << std::endl;
    std::cout << "***** FIM algorithm, 2D, non-cartesian mesh ******" << std::endl;
  }

  {

    Scarf3D::options.hurst = hurst;
    Scarf3D::options.alpha = alpha;
    Scarf3D::options.npoi  = npoi;
    Scarf3D::options.poi   = poi;
    Scarf3D::options.taper = taper;
    Scarf3D::options.mute  = mute;

    Scarf3D::Initialize<fim> S(nd, npts, x, y, dh, acf, cl, sigma);

    watch_start(&tictoc);
    S.execute(seed, field, stats);
    watch_stop(&tictoc);

    if (world_rank == 0){
      std::cout << ""                                                                                       << std::endl;
      std::cout << "Statistics for current simulation"                                                      << std::endl;
      std::cout << "*************************************************"                                      << std::endl;
      std::cout << "Elapsed time                  |" << std::setw(12) << (float) tictoc   << " sec" << " |" << std::endl;
      std::cout << "   + spectrum                 |" << std::setw(12) << (float) stats[4] << " sec" << " |" << std::endl;
      std::cout << "   + symmetry                 |" << std::setw(12) << (float) stats[5] << " sec" << " |" << std::endl;
      std::cout << "   + FFT                      |" << std::setw(12) << (float) stats[6] << " sec" << " |" << std::endl;
      std::cout << "   + interpolation            |" << std::setw(12) << (float) stats[7] << " sec" << " |" << std::endl;
      std::cout << "------------------------------|-----------------|"                                      << std::endl;
      std::cout << "Domain too small?             |" << std::setw(12) << (bool) stats[0]  << "     |"       << std::endl;
      std::cout << "Grid-step too large?          |" << std::setw(12) << (bool) stats[1]  << "     |"       << std::endl;
      std::cout << "------------------------------|-----------------|"                                      << std::endl;
      std::cout << "Discrete standard deviation   |" << std::setw(12) << (float) stats[2] << "     |"       << std::endl;
      std::cout << "Discrete mean value           |" << std::setw(12) << (float) stats[3] << "     |"       << std::endl;
      std::cout << "*************************************************"                                      << std::endl;
    }

    // I/O not possible for non-cartesian mesh because "fs" and "fe" are not defined

    // call destructor explicitly
    //S.~Initialize();
  }

  // ================================================================================================================================
  // --------------------------------------------------------------------------------------------------------------------------------
  // tests SRM algorithm
  // --------------------------------------------------------------------------------------------------------------------------------
  // ================================================================================================================================

#ifdef SPECTRAL

  if (world_rank == 0){
    std::cout << ""                                                  << std::endl;
    std::cout << "*************************************************" << std::endl;
    std::cout << "****** SRM algorithm, 2D, cartesian mesh *******" << std::endl;
  }

  {

    Scarf3D::options.hurst = hurst;
    Scarf3D::options.alpha = alpha;
    Scarf3D::options.npoi  = npoi;
    Scarf3D::options.poi   = poi;
    Scarf3D::options.taper = taper;
    Scarf3D::options.mute  = mute;

    Scarf3D::Initialize<srm> S(nd, fs, fe, ds, acf, cl, sigma);

    watch_start(&tictoc);
    S.execute(seed, field, stats);
    watch_stop(&tictoc);

    if (world_rank == 0){
      std::cout << ""                                                                               << std::endl;
      std::cout << "Statistics for current simulation"                                              << std::endl;
      std::cout << "*************************************************"                              << std::endl;
      std::cout << "Elapsed time                  |" << std::setw(12) << (float) tictoc   << " sec" << " |"<< std::endl;
      std::cout << "   + CPU (main loop)          |" << std::setw(12) << (float) stats[4] << " sec" << " |"<< std::endl;
      std::cout << "   + GPU (main loop)          |" << std::setw(12) << (float) stats[5] << " sec" << " |"<< std::endl;
      std::cout << "------------------------------|-----------------|"                              << std::endl;
      std::cout << "Domain too small?             |" << std::setw(12) << (bool) stats[0]            << " |"<< std::endl;
      std::cout << "Grid-step too large?          |" << std::setw(12) << (bool) stats[1]            << " |"<< std::endl;
      std::cout << "------------------------------|-----------------|"                              << std::endl;
      std::cout << "Discrete standard deviation   |" << std::setw(12) << (float) stats[2]           << " |"<< std::endl;
      std::cout << "Discrete mean value           |" << std::setw(12) << (float) stats[3]           << " |"<< std::endl;
      std::cout << "------------------------------|-----------------|"                              << std::endl;
    }


    watch_start(&tictoc);
    S.io(nd, n, field, "srm_cart_whole_2d", &nwriters);
    watch_stop(&tictoc);

    if (world_rank == 0) {
      std::cout << "I/O time                      |" << std::setw(12) << (float) tictoc << " sec" << " |" << std::endl;
      std::cout << "*************************************************"                                    << std::endl;
    }

    // call destructor explicitly
    // S.~Initialize();
    }

    if (world_rank == 0){
      std::cout << ""                                                  << std::endl;
      std::cout << "*************************************************" << std::endl;
      std::cout << "***** SRM algorithm, 2D, non-cartesian mesh ******" << std::endl;
    }

    {

      Scarf3D::options.hurst = hurst;
      Scarf3D::options.alpha = alpha;
      Scarf3D::options.npoi  = npoi;
      Scarf3D::options.poi   = poi;
      Scarf3D::options.taper = taper;
      Scarf3D::options.mute  = mute;

      Scarf3D::Initialize<srm> S(nd, npts, x, y, dh, acf, cl, sigma);

      watch_start(&tictoc);
      S.execute(seed, field, stats);
      watch_stop(&tictoc);

      if (world_rank == 0){
        std::cout << ""                                                                               << std::endl;
        std::cout << "Statistics for current simulation"                                              << std::endl;
        std::cout << "*************************************************"                              << std::endl;
        std::cout << "Elapsed time                  |" << std::setw(12) << (float) tictoc   << " sec" << " |"<< std::endl;
        std::cout << "   + CPU (main loop)          |" << std::setw(12) << (float) stats[4] << " sec" << " |"<< std::endl;
        std::cout << "   + GPU (main loop)          |" << std::setw(12) << (float) stats[5] << " sec" << " |"<< std::endl;
        std::cout << "------------------------------|-----------------|"                              << std::endl;
        std::cout << "Domain too small?             |" << std::setw(12) << (bool) stats[0]            << " |"<< std::endl;
        std::cout << "Grid-step too large?          |" << std::setw(12) << (bool) stats[1]            << " |"<< std::endl;
        std::cout << "------------------------------|-----------------|"                              << std::endl;
        std::cout << "Discrete standard deviation   |" << std::setw(12) << (float) stats[2]           << " |"<< std::endl;
        std::cout << "Discrete mean value           |" << std::setw(12) << (float) stats[3]           << " |"<< std::endl;
        std::cout << "*************************************************"                              << std::endl;
      }

      // call destructor explicitly
      //S.~Initialize();
    }

#endif

  // release resources
  delete[] x;
  delete[] y;
  delete[] field;
  delete[] poi;

  // ================================================================================================================================
  // --------------------------------------------------------------------------------------------------------------------------------
  // create sample cartesian mesh: 3D case
  // --------------------------------------------------------------------------------------------------------------------------------
  // ================================================================================================================================
  nd = 3;

  sample_mesh(&nd, &world_rank, &world_size, n, fs, fe);

  for (int i = 0; i < nd; i++){
    dims[i] = (fe[i] - fs[i] + 1);
  }

  x     = new real[dims[0]*dims[1]*dims[2]];
  y     = new real[dims[0]*dims[1]*dims[2]];
  z     = new real[dims[0]*dims[1]*dims[2]];
  field = new real[dims[0]*dims[1]*dims[2]];

  for (int k = 0; k < dims[2]; k++){
    for (int j = 0; j < dims[1]; j++){
      for (int i = 0; i < dims[0]; i++){
        c    = (k * dims[0] * dims[1]) + (j * dims[0]) + i;
        x[c] = (i + fs[0] - 1) * ds;
        y[c] = (j + fs[1] - 1) * ds;
        z[c] = (k + fs[2] - 1) * ds;
      }
    }
  }

  npts = dims[0] * dims[1] * dims[2];

  poi = new real[nd * npoi];

  poi[0] = 400 * ds;
  poi[1] = 250 * ds;
  poi[2] = 100 * ds;
  poi[3] = 200 * ds;
  poi[4] = 150 * ds;
  poi[5] = 50 * ds;

  // ================================================================================================================================
  // --------------------------------------------------------------------------------------------------------------------------------
  // tests FIM algorithm
  // --------------------------------------------------------------------------------------------------------------------------------
  // ================================================================================================================================

  if (world_rank == 0){
    std::cout << ""                                                  << std::endl;
    std::cout << "*************************************************" << std::endl;
    std::cout << "****** FIM algorithm, 3D, cartesian mesh *******" << std::endl;
  }

  {

    Scarf3D::options.hurst = hurst;
    /*
    Scarf3D::options.alpha = alpha;
    Scarf3D::options.beta  = beta;
    Scarf3D::options.npoi  = npoi;
    Scarf3D::options.poi   = poi;
    Scarf3D::options.taper = taper;
    Scarf3D::options.mute  = mute;
    */

    Scarf3D::Initialize<fim> S(nd, fs, fe, ds, acf, cl, sigma);

    watch_start(&tictoc);
    S.execute(seed, field, stats);
    watch_stop(&tictoc);

    if (world_rank == 0){
      std::cout << ""                                                                                       << std::endl;
      std::cout << "Statistics for current simulation"                                                      << std::endl;
      std::cout << "*************************************************"                                      << std::endl;
      std::cout << "Elapsed time                  |" << std::setw(12) << (float) tictoc   << " sec" << " |" << std::endl;
      std::cout << "   + spectrum                 |" << std::setw(12) << (float) stats[4] << " sec" << " |" << std::endl;
      std::cout << "   + symmetry                 |" << std::setw(12) << (float) stats[5] << " sec" << " |" << std::endl;
      std::cout << "   + FFT                      |" << std::setw(12) << (float) stats[6] << " sec" << " |" << std::endl;
      std::cout << "   + interpolation            |" << std::setw(12) << (float) stats[7] << " sec" << " |" << std::endl;
      std::cout << "------------------------------|-----------------|"                                      << std::endl;
      std::cout << "Domain too small?             |" << std::setw(12) << (bool) stats[0]  << "     |"       << std::endl;
      std::cout << "Grid-step too large?          |" << std::setw(12) << (bool) stats[1]  << "     |"       << std::endl;
      std::cout << "------------------------------|-----------------|"                                      << std::endl;
      std::cout << "Discrete standard deviation   |" << std::setw(12) << (float) stats[2] << "     |"       << std::endl;
      std::cout << "Discrete mean value           |" << std::setw(12) << (float) stats[3] << "     |"       << std::endl;
      std::cout << "------------------------------|-----------------|"                                      << std::endl;
    }

    // IO
    watch_start(&tictoc);
    S.io(n, "x", n[0]/2, field, "fim_cart_xslice");
    S.io(n, "y", n[1]/2, field, "fim_cart_yslice");
    S.io(n, "z", n[2]/2, field, "fim_cart_zslice");
    watch_stop(&tictoc);

    if (world_rank == 0) {
       std::cout << "I/O time (slices)             |" << std::setw(12) << (float) tictoc << " sec" << " |" << std::endl;
    }

    watch_start(&tictoc);
    S.io(nd, n, field, "fim_cart_whole_3d", &nwriters);
    watch_stop(&tictoc);

    if (world_rank == 0) {
      std::cout << "I/O time                      |" << std::setw(12) << (float) tictoc << " sec" << " |"<< std::endl;
      std::cout << "*************************************************"                                   << std::endl;
    }

    // call destructor explicitly
    //S.~Initialize();
  }

  if (world_rank == 0){
    std::cout << ""                                                  << std::endl;
    std::cout << "*************************************************" << std::endl;
    std::cout << "***** FIM algorithm, 3D, non-cartesian mesh ******" << std::endl;
  }

  {

    Scarf3D::options.hurst = hurst;
    /*
    Scarf3D::options.alpha = alpha;
    Scarf3D::options.beta  = beta;
    Scarf3D::options.npoi  = npoi;
    Scarf3D::options.poi   = poi;
    Scarf3D::options.taper = taper;
    Scarf3D::options.mute  = mute;
    */

    Scarf3D::Initialize<fim> S(nd, npts, x, y, z, dh, acf, cl, sigma);

    watch_start(&tictoc);
    S.execute(seed, field, stats);
    watch_stop(&tictoc);

    if (world_rank == 0){
      std::cout << ""                                                                                       << std::endl;
      std::cout << "Statistics for current simulation"                                                      << std::endl;
      std::cout << "*************************************************"                                      << std::endl;
      std::cout << "Elapsed time                  |" << std::setw(12) << (float) tictoc   << " sec" << " |" << std::endl;
      std::cout << "   + spectrum                 |" << std::setw(12) << (float) stats[4] << " sec" << " |" << std::endl;
      std::cout << "   + symmetry                 |" << std::setw(12) << (float) stats[5] << " sec" << " |" << std::endl;
      std::cout << "   + FFT                      |" << std::setw(12) << (float) stats[6] << " sec" << " |" << std::endl;
      std::cout << "   + interpolation            |" << std::setw(12) << (float) stats[7] << " sec" << " |" << std::endl;
      std::cout << "------------------------------|-----------------|"                                      << std::endl;
      std::cout << "Domain too small?             |" << std::setw(12) << (bool) stats[0]  << "     |"       << std::endl;
      std::cout << "Grid-step too large?          |" << std::setw(12) << (bool) stats[1]  << "     |"       << std::endl;
      std::cout << "------------------------------|-----------------|"                                      << std::endl;
      std::cout << "Discrete standard deviation   |" << std::setw(12) << (float) stats[2] << "     |"       << std::endl;
      std::cout << "Discrete mean value           |" << std::setw(12) << (float) stats[3] << "     |"       << std::endl;
      std::cout << "*************************************************"                                      << std::endl;
    }

    // call destructor explicitly
    //S.~Initialize();
  }


  // ================================================================================================================================
  // --------------------------------------------------------------------------------------------------------------------------------
  // tests SRM algorithm
  // --------------------------------------------------------------------------------------------------------------------------------
  // ================================================================================================================================

#ifdef SPECTRAL

  if (world_rank == 0){
    std::cout << ""                                                  << std::endl;
    std::cout << "*************************************************" << std::endl;
    std::cout << "****** SRM algorithm, 3D, cartesian mesh *******" << std::endl;
  }

  {

    Scarf3D::options.hurst = hurst;
    /*
    Scarf3D::options.alpha = alpha;
    Scarf3D::options.beta  = beta;
    Scarf3D::options.npoi  = npoi;
    Scarf3D::options.poi   = poi;
    Scarf3D::options.taper = taper;
    Scarf3D::options.mute  = mute;
    */

    Scarf3D::Initialize<srm> S(nd, fs, fe, ds, acf, cl, sigma);

    watch_start(&tictoc);
    S.execute(seed, field, stats);
    watch_stop(&tictoc);

    if (world_rank == 0){
      std::cout << ""                                                                               << std::endl;
      std::cout << "Statistics for current simulation"                                              << std::endl;
      std::cout << "*************************************************"                              << std::endl;
      std::cout << "Elapsed time                  |" << std::setw(12) << (float) tictoc   << " sec" << " |"<< std::endl;
      std::cout << "   + CPU (main loop)          |" << std::setw(12) << (float) stats[4] << " sec" << " |"<< std::endl;
      std::cout << "   + GPU (main loop)          |" << std::setw(12) << (float) stats[5] << " sec" << " |"<< std::endl;
      std::cout << "------------------------------|-----------------|"                              << std::endl;
      std::cout << "Domain too small?             |" << std::setw(12) << (bool) stats[0]            << " |"<< std::endl;
      std::cout << "Grid-step too large?          |" << std::setw(12) << (bool) stats[1]            << " |"<< std::endl;
      std::cout << "------------------------------|-----------------|"                              << std::endl;
      std::cout << "Discrete standard deviation   |" << std::setw(12) << (float) stats[2]           << " |"<< std::endl;
      std::cout << "Discrete mean value           |" << std::setw(12) << (float) stats[3]           << " |"<< std::endl;
      std::cout << "------------------------------|-----------------|"                              << std::endl;
    }

    // IO
    watch_start(&tictoc);
    S.io(n, "x", n[0]/2, field, "srm_cart_xslice");
    S.io(n, "y", n[1]/2, field, "srm_cart_yslice");
    S.io(n, "z", n[2]/2, field, "srm_cart_zslice");
    watch_stop(&tictoc);

    if (world_rank == 0) {
       std::cout << "I/O time (slices)             |" << std::setw(12) << (float) tictoc << " sec" << "  |" << std::endl;
    }

    watch_start(&tictoc);
    S.io(nd, n, field, "srm_cart_whole_3d", &nwriters);
    watch_stop(&tictoc);

    if (world_rank == 0) {
      std::cout << "I/O time                      |" << std::setw(12) << (float) tictoc << " sec" << " |"<< std::endl;
      std::cout << "*************************************************"                                   << std::endl;
    }

    // call destructor explicitly
    // S.~Initialize();
  }

  if (world_rank == 0){
    std::cout << ""                                                  << std::endl;
    std::cout << "*************************************************" << std::endl;
    std::cout << "***** SRM algorithm, 3D, non-cartesian mesh ******" << std::endl;
  }

  {

    Scarf3D::options.hurst = hurst;
    /*
    Scarf3D::options.alpha = alpha;
    Scarf3D::options.beta  = beta;
    Scarf3D::options.npoi  = npoi;
    Scarf3D::options.poi   = poi;
    Scarf3D::options.taper = taper;
    Scarf3D::options.mute  = mute;
    */

    Scarf3D::Initialize<srm> S(nd, npts, x, y, z, dh, acf, cl, sigma);

    watch_start(&tictoc);
    S.execute(seed, field, stats);
    watch_stop(&tictoc);

    if (world_rank == 0){
      std::cout << ""                                                                               << std::endl;
      std::cout << "Statistics for current simulation"                                              << std::endl;
      std::cout << "*************************************************"                              << std::endl;
      std::cout << "Elapsed time                  |" << std::setw(12) << (float) tictoc   << " sec" << " |"<< std::endl;
      std::cout << "   + CPU (main loop)          |" << std::setw(12) << (float) stats[4] << " sec" << " |"<< std::endl;
      std::cout << "   + GPU (main loop)          |" << std::setw(12) << (float) stats[5] << " sec" << " |"<< std::endl;
      std::cout << "------------------------------|-----------------|"                              << std::endl;
      std::cout << "Domain too small?             |" << std::setw(12) << (bool) stats[0]            << " |"<< std::endl;
      std::cout << "Grid-step too large?          |" << std::setw(12) << (bool) stats[1]            << " |"<< std::endl;
      std::cout << "------------------------------|-----------------|"                              << std::endl;
      std::cout << "Discrete standard deviation   |" << std::setw(12) << (float) stats[2]           << " |"<< std::endl;
      std::cout << "Discrete mean value           |" << std::setw(12) << (float) stats[3]           << " |"<< std::endl;
      std::cout << "*************************************************"                              << std::endl;
    }

    // call destructor explicitly
    //S.~Initialize();
  }

#endif

   // release resources
   delete[] x;
   delete[] y;
   delete[] z;
   delete[] field;
   delete[] poi;

   MPI_Barrier(MPI_COMM_WORLD);

   MPI_Finalize();

   return 0;

}
