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
   To a sample C driver program for the SCARF3D library

 Revisions:
     Date                    Description of change
     ====                    =====================
   04/05/20                  original version
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "scarf3d.h"


// FORTRAN subroutine prototype
extern void sample_mesh(const int* nd, const int* rank, const int* ntasks, const int n[], int fs[], int fe[]);

// timing functions
void watch_start(double* t){MPI_Barrier(MPI_COMM_WORLD); *t = MPI_Wtime();}
void watch_stop(double* t){MPI_Barrier(MPI_COMM_WORLD); *t = MPI_Wtime() - *t;}


int main(){

  // initialise MPI
  MPI_Init(NULL, NULL);

  // get number of tasks
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // get rank of each task
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (world_rank == 0){
    printf("\n");
    printf("This program will compute small 2D and 3D random fields\n");
    printf("using the FIM (and SRM - if selected at compile-time) algorithm(s).\n");
    printf("Test should complete in a few minutes. Output can be inspected\n");
    printf("with Matlab/Octave script 'scarf3d.m' and 'vizme.m'.\n");
    printf("\n");
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // ===========================================================================
  // ---------------------------------------------------------------------------
  // Set mandatory parameters

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

  // ===========================================================================
  // ---------------------------------------------------------------------------
  // Add extra parameters

  const real hurst = 0.25;

  // rotation angles
  const real alpha = 30;
  const real beta  = -10;

  // number of POI
  const int npoi = 2;

  const real taper = 2500;
  const real mute  = 500;

  // I/O writers (only for PFS)
  const int nwriters[1] = {2};

  double tictoc;
  int fs[3], fe[3], dims[3];
  int npts, nd;
  long c;
  real stats[8];
  real *x, *y, *z, *field, *poi;

  struct scarf_opt options;

  // ================================================================================================================================
  // --------------------------------------------------------------------------------------------------------------------------------
  // create sample cartesian mesh: 2D case
  // --------------------------------------------------------------------------------------------------------------------------------
  // ================================================================================================================================

  nd = 2;

  // domain decomposition
  sample_mesh(&nd, &world_rank, &world_size, n, fs, fe);

  for (int i = 0; i < nd; i++){
    dims[i] = (fe[i] - fs[i] + 1);
  }

  x     = (real*) malloc(dims[0]*dims[1] * sizeof(real));
  y     = (real*) malloc(dims[0]*dims[1] * sizeof(real));
  field = (real*) malloc(dims[0]*dims[1] * sizeof(real));

  for (int j = 0; j < dims[1]; j++){
    for (int i = 0; i < dims[0]; i++){
      c    = (j * dims[0]) + i;
      x[c] = (i + fs[0] - 1) * ds;
      y[c] = (j + fs[1] - 1) * ds;
    }
  }

  npts = dims[0] * dims[1];

  poi = (real*) malloc(nd * npoi * sizeof(real));

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
    printf("\n");
    printf("*************************************************\n");
    printf("****** FIM algorithm, 2D, cartesian mesh *******\n");
  }

  scarf_opt_init(&options);

  options.hurst = hurst;
  options.alpha = alpha;
  options.npoi  = npoi;
  options.poi   = poi;
  options.taper = taper;
  options.mute  = mute;

  scarf_cart_initialize(nd, fs, fe, ds, acf, cl, sigma, &options);

  watch_start(&tictoc);
  scarf_execute(seed, field, stats);
  watch_stop(&tictoc);

  if (world_rank == 0){
    printf("\n");
    printf("Statistics for current simulation\n");
    printf("*************************************************\n");
    printf("Elapsed time                  |%12.5f sec |\n", (float) tictoc);
    printf("   + spectrum                 |%12.5f sec |\n", (float) stats[4]);
    printf("   + symmetry                 |%12.5f sec |\n", (float) stats[5]);
    printf("   + FFT                      |%12.5f sec |\n", (float) stats[6]);
    printf("   + interpolation            |%12.5f sec |\n", (float) stats[7]);
    printf("------------------------------|-----------------|\n");
    printf("Domain too small?             |%12s     |\n", stats[0] ? "true" : "false");
    printf("Grid-step too large?          |%12s     |\n", stats[1] ? "true" : "false");
    printf("------------------------------|-----------------|\n");
    printf("Discrete standard deviation   |%12.5f     |\n", (float) stats[2]);
    printf("Discrete mean value           |%12.5f     |\n", (float) stats[3]);
    printf("------------------------------|-----------------|\n");
  }

  watch_start(&tictoc);
  scarf_io_one(nd, n, field, "fim_cart_whole_2d", nwriters);
  watch_stop(&tictoc);

  if (world_rank == 0) {
    printf("I/O time                      |%12.5f sec |\n", (float) tictoc);
    printf("*************************************************\n");
  }

  scarf_finalize();

  MPI_Barrier(MPI_COMM_WORLD);

  if (world_rank == 0){
    printf("\n");
    printf("*************************************************\n");
    printf("***** FIM algorithm, 2D, non-cartesian mesh ******\n");
  }

  scarf_opt_init(&options);

  options.hurst = hurst;
  options.alpha = alpha;
  options.npoi  = npoi;
  options.poi   = poi;
  options.taper = taper;
  options.mute  = mute;

  scarf_nocart_initialize(nd, npts, x, y, NULL, ds, acf, cl, sigma, &options);

  watch_start(&tictoc);
  scarf_execute(seed, field, stats);
  watch_stop(&tictoc);

  if (world_rank == 0){
    printf("\n");
    printf("Statistics for current simulation\n");
    printf("*************************************************\n");
    printf("Elapsed time                  |%12.5f sec |\n", (float) tictoc);
    printf("   + spectrum                 |%12.5f sec |\n", (float) stats[4]);
    printf("   + symmetry                 |%12.5f sec |\n", (float) stats[5]);
    printf("   + FFT                      |%12.5f sec |\n", (float) stats[6]);
    printf("   + interpolation            |%12.5f sec |\n", (float) stats[7]);
    printf("------------------------------|-----------------|\n");
    printf("Domain too small?             |%12s     |\n", stats[0] ? "true" : "false");
    printf("Grid-step too large?          |%12s     |\n", stats[1] ? "true" : "false");
    printf("------------------------------|-----------------|\n");
    printf("Discrete standard deviation   |%12.5f     |\n", (float) stats[2]);
    printf("Discrete mean value           |%12.5f     |\n", (float) stats[3]);
    printf("*************************************************\n");
  }

  scarf_finalize();

  // ================================================================================================================================
  // --------------------------------------------------------------------------------------------------------------------------------
  // tests SRM algorithm
  // --------------------------------------------------------------------------------------------------------------------------------
  // ================================================================================================================================

#ifdef SPECTRAL

  if (world_rank == 0){
    printf("\n");
    printf("*************************************************\n");
    printf("****** SRM algorithm, 2D, cartesian mesh *******\n");
  }

  scarf_opt_init(&options);

  options.hurst  = hurst;
  options.alpha  = alpha;
  options.npoi   = npoi;
  options.poi    = poi;
  options.taper  = taper;
  options.mute   = mute;
  options.method = 1;

  scarf_cart_initialize(nd, fs, fe, ds, acf, cl, sigma, &options);

  watch_start(&tictoc);
  scarf_execute(seed, field, stats);
  watch_stop(&tictoc);

  if (world_rank == 0){
    printf("\n");
    printf("Statistics for current simulation\n");
    printf("*************************************************\n");
    printf("Elapsed time                  |%12.5f sec |\n", (float) tictoc);
    printf("   + CPU (main loop)          |%12.5f sec |\n", (float) stats[4]);
    printf("   + GPU (main loop)          |%12.5f sec |\n", (float) stats[5]);
    printf("------------------------------|-----------------|\n");
    printf("Domain too small?             |%12s     |\n", stats[0] ? "true" : "false");
    printf("Grid-step too large?          |%12s     |\n", stats[1] ? "true" : "false");
    printf("------------------------------|-----------------|\n");
    printf("Discrete standard deviation   |%12.5f     |\n", (float) stats[2]);
    printf("Discrete mean value           |%12.5f     |\n", (float) stats[3]);
    printf("------------------------------|-----------------|\n");
  }

  watch_start(&tictoc);
  scarf_io_one(nd, n, field, "srm_cart_whole_2d", nwriters);
  watch_stop(&tictoc);

  if (world_rank == 0) {
    printf("I/O time                      |%12.5f sec |\n", (float) tictoc);
    printf("*************************************************\n");
  }

  scarf_finalize();

  if (world_rank == 0){
    printf("\n");
    printf("*************************************************\n");
    printf("****** SRM algorithm, 2D, cartesian mesh *******\n");
  }

  scarf_opt_init(&options);

  options.hurst  = hurst;
  options.alpha  = alpha;
  options.npoi   = npoi;
  options.poi    = poi;
  options.taper  = taper;
  options.mute   = mute;
  options.method = 1;

  scarf_nocart_initialize(nd, npts, x, y, NULL, ds, acf, cl, sigma, &options);

  watch_start(&tictoc);
  scarf_execute(seed, field, stats);
  watch_stop(&tictoc);

  if (world_rank == 0){
    printf("\n");
    printf("Statistics for current simulation\n");
    printf("*************************************************\n");
    printf("Elapsed time                  |%12.5f sec |\n", (float) tictoc);
    printf("   + CPU (main loop)          |%12.5f sec |\n", (float) stats[4]);
    printf("   + GPU (main loop)          |%12.5f sec |\n", (float) stats[5]);
    printf("------------------------------|-----------------|\n");
    printf("Domain too small?             |%12s     |\n", stats[0] ? "true" : "false");
    printf("Grid-step too large?          |%12s     |\n", stats[1] ? "true" : "false");
    printf("------------------------------|-----------------|\n");
    printf("Discrete standard deviation   |%12.5f     |\n", (float) stats[2]);
    printf("Discrete mean value           |%12.5f     |\n", (float) stats[3]);
    printf("*************************************************\n");
  }

  scarf_finalize();

#endif

  free(x);
  free(y);
  free(field);
  free(poi);

  // ================================================================================================================================
  // --------------------------------------------------------------------------------------------------------------------------------
  // create sample cartesian mesh: 3D case
  // --------------------------------------------------------------------------------------------------------------------------------
  // ================================================================================================================================

  nd = 3;

  // domain decomposition
  sample_mesh(&nd, &world_rank, &world_size, n, fs, fe);

  for (int i = 0; i < nd; i++){
    dims[i] = (fe[i] - fs[i] + 1);
  }

  x     = (real*) malloc(dims[0]*dims[1]*dims[2] * sizeof(real));
  y     = (real*) malloc(dims[0]*dims[1]*dims[2] * sizeof(real));
  z     = (real*) malloc(dims[0]*dims[1]*dims[2] * sizeof(real));
  field = (real*) malloc(dims[0]*dims[1]*dims[2] * sizeof(real));

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

  poi = (real*) malloc(nd * npoi * sizeof(real));

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
    printf("\n");
    printf("*************************************************\n");
    printf("****** FIM algorithm, 3D, cartesian mesh *******\n");
  }

  scarf_opt_init(&options);

  options.hurst  = hurst;
  /*
  options.alpha  = alpha;
  options.beta   = beta;
  options.npoi   = npoi;
  options.poi    = poi;
  options.taper  = taper;
  options.mute   = mute;
  */

  scarf_cart_initialize(nd, fs, fe, ds, acf, cl, sigma, &options);

  watch_start(&tictoc);
  scarf_execute(seed, field, stats);
  watch_stop(&tictoc);

  if (world_rank == 0){
    printf("\n");
    printf("Statistics for current simulation\n");
    printf("*************************************************\n");
    printf("Elapsed time                  |%12.5f sec |\n", (float) tictoc);
    printf("   + spectrum                 |%12.5f sec |\n", (float) stats[4]);
    printf("   + symmetry                 |%12.5f sec |\n", (float) stats[5]);
    printf("   + FFT                      |%12.5f sec |\n", (float) stats[6]);
    printf("   + interpolation            |%12.5f sec |\n", (float) stats[7]);
    printf("------------------------------|-----------------|\n");
    printf("Domain too small?             |%12s     |\n", stats[0] ? "true" : "false");
    printf("Grid-step too large?          |%12s     |\n", stats[1] ? "true" : "false");
    printf("------------------------------|-----------------|\n");
    printf("Discrete standard deviation   |%12.5f     |\n", (float) stats[2]);
    printf("Discrete mean value           |%12.5f     |\n", (float) stats[3]);
    printf("------------------------------|-----------------|\n");
  }

  watch_start(&tictoc);
  scarf_io_slice(n, "x", n[0]/2, field, "fim_cart_xslice");
  scarf_io_slice(n, "y", n[1]/2, field, "fim_cart_yslice");
  scarf_io_slice(n, "z", n[2]/2, field, "fim_cart_zslice");
  watch_stop(&tictoc);

  if (world_rank == 0) {
    printf("I/O time (slices)             |%12.5f sec |\n", (float) tictoc);
  }

  watch_start(&tictoc);
  scarf_io_one(nd, n, field, "fim_cart_whole_3d", nwriters);
  watch_stop(&tictoc);

  if (world_rank == 0) {
    printf("I/O time                      |%12.5f sec |\n", (float) tictoc);
    printf("*************************************************\n");
  }

  scarf_finalize();

  MPI_Barrier(MPI_COMM_WORLD);

  if (world_rank == 0){
    printf("\n");
    printf("*************************************************\n");
    printf("***** FIM algorithm, 3D, non-cartesian mesh ******\n");
  }

  scarf_opt_init(&options);

  options.hurst  = hurst;
  /*
  options.alpha  = alpha;
  options.beta   = beta;
  options.npoi   = npoi;
  options.poi    = poi;
  options.taper  = taper;
  options.mute   = mute;
  */

  scarf_nocart_initialize(nd, npts, x, y, z, ds, acf, cl, sigma, &options);

  watch_start(&tictoc);
  scarf_execute(seed, field, stats);
  watch_stop(&tictoc);

  if (world_rank == 0){
    printf("\n");
    printf("Statistics for current simulation\n");
    printf("*************************************************\n");
    printf("Elapsed time                  |%12.5f sec |\n", (float) tictoc);
    printf("   + spectrum                 |%12.5f sec |\n", (float) stats[4]);
    printf("   + symmetry                 |%12.5f sec |\n", (float) stats[5]);
    printf("   + FFT                      |%12.5f sec |\n", (float) stats[6]);
    printf("   + interpolation            |%12.5f sec |\n", (float) stats[7]);
    printf("------------------------------|-----------------|\n");
    printf("Domain too small?             |%12s     |\n", stats[0] ? "true" : "false");
    printf("Grid-step too large?          |%12s     |\n", stats[1] ? "true" : "false");
    printf("------------------------------|-----------------|\n");
    printf("Discrete standard deviation   |%12.5f     |\n", (float) stats[2]);
    printf("Discrete mean value           |%12.5f     |\n", (float) stats[3]);
    printf("*************************************************\n");
  }

  scarf_finalize();

  // ================================================================================================================================
  // --------------------------------------------------------------------------------------------------------------------------------
  // tests SRM algorithm
  // --------------------------------------------------------------------------------------------------------------------------------
  // ================================================================================================================================

#ifdef SPECTRAL

  if (world_rank == 0){
    printf("\n");
    printf("*************************************************\n");
    printf("****** SRM algorithm, 3D, cartesian mesh *******\n");
  }

  scarf_opt_init(&options);

  options.hurst  = hurst;
  /*
  options.alpha  = alpha;
  options.beta   = beta;
  options.npoi   = npoi;
  options.poi    = poi;
  options.taper  = taper;
  options.mute   = mute;
  options.method = 1;
  */

  scarf_cart_initialize(nd, fs, fe, ds, acf, cl, sigma, &options);

  watch_start(&tictoc);
  scarf_execute(seed, field, stats);
  watch_stop(&tictoc);

  if (world_rank == 0){
    printf("\n");
    printf("Statistics for current simulation\n");
    printf("*************************************************\n");
    printf("Elapsed time                  |%12.5f sec |\n", (float) tictoc);
    printf("   + CPU (main loop)          |%12.5f sec |\n", (float) stats[4]);
    printf("   + GPU (main loop)          |%12.5f sec |\n", (float) stats[5]);
    printf("------------------------------|-----------------|\n");
    printf("Domain too small?             |%12s     |\n", stats[0] ? "true" : "false");
    printf("Grid-step too large?          |%12s     |\n", stats[1] ? "true" : "false");
    printf("------------------------------|-----------------|\n");
    printf("Discrete standard deviation   |%12.5f     |\n", (float) stats[2]);
    printf("Discrete mean value           |%12.5f     |\n", (float) stats[3]);
    printf("------------------------------|-----------------|\n");
  }

  watch_start(&tictoc);
  scarf_io_slice(n, "x", n[0]/2, field, "srm_cart_xslice");
  scarf_io_slice(n, "y", n[1]/2, field, "srm_cart_yslice");
  scarf_io_slice(n, "z", n[2]/2, field, "srm_cart_zslice");
  watch_stop(&tictoc);

  if (world_rank == 0) {
    printf("I/O time (slices)             |%12.5f sec    |\n", (float) tictoc);
  }

  watch_start(&tictoc);
  scarf_io_one(nd, n, field, "srm_cart_whole_3d", nwriters);
  watch_stop(&tictoc);

  if (world_rank == 0) {
    printf("I/O time                      |%12.5f sec    |\n", (float) tictoc);
    printf("*************************************************\n");
  }

  scarf_finalize();

  if (world_rank == 0){
    printf("\n");
    printf("*************************************************\n");
    printf("****** SRM algorithm, 3D, non-cartesian mesh *******\n");
  }

  scarf_opt_init(&options);

  options.hurst  = hurst;
  /*
  options.alpha  = alpha;
  options.beta   = beta;
  options.npoi   = npoi;
  options.poi    = poi;
  options.taper  = taper;
  options.mute   = mute;
  options.method = 1;
  */

  scarf_nocart_initialize(nd, npts, x, y, z, ds, acf, cl, sigma, &options);

  watch_start(&tictoc);
  scarf_execute(seed, field, stats);
  watch_stop(&tictoc);

  if (world_rank == 0){
    printf("\n");
    printf("Statistics for current simulation\n");
    printf("*************************************************\n");
    printf("Elapsed time                  |%12.5f sec |\n", (float) tictoc);
    printf("   + CPU (main loop)          |%12.5f sec |\n", (float) stats[4]);
    printf("   + GPU (main loop)          |%12.5f sec |\n", (float) stats[5]);
    printf("------------------------------|-----------------|\n");
    printf("Domain too small?             |%12s     |\n", stats[0] ? "true" : "false");
    printf("Grid-step too large?          |%12s     |\n", stats[1] ? "true" : "false");
    printf("------------------------------|-----------------|\n");
    printf("Discrete standard deviation   |%12.5f     |\n", (float) stats[2]);
    printf("Discrete mean value           |%12.5f     |\n", (float) stats[3]);
    printf("*************************************************\n");
  }

  scarf_finalize();

#endif

  free(x);
  free(y);
  free(z);
  free(field);
  free(poi);

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();

  return 0;

}
