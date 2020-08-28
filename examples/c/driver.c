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

  MPI_Barrier(MPI_COMM_WORLD);

  // ===========================================================================
  // ---------------------------------------------------------------------------
  // Set mandatory parameters

  // number of grid points in model
  const int n[3] = {500, 450, 400};

  // grid-step
  const real ds = 50.;

  // autocorrelation function (0=von karman/exponential, 1=gaussian)
  const int acf = 0;

  // correlation length
  const real cl[3] = {2000., 500., 100.};

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

  const real taper = 5000;
  const real mute  = 1000;

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
  // create sample structured mesh: 2D case
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
    printf("****** FIM algorithm, 2D, structured mesh *******\n");
  }

  scarf_opt_init(&options);

  options.hurst = hurst;
  options.alpha = alpha;
  options.npoi  = npoi;
  options.poi   = poi;
  options.taper = taper;
  options.mute  = mute;

  scarf_struct_initialize(nd, fs, fe, ds, acf, cl, sigma, &options);

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
  scarf_io_one(nd, n, field, "fim_struct_whole_2d", nwriters);
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
    printf("***** FIM algorithm, 2D, unstructured mesh ******\n");
  }

  scarf_opt_init(&options);

  options.hurst = hurst;
  options.alpha = alpha;
  options.npoi  = npoi;
  options.poi   = poi;
  options.taper = taper;
  options.mute  = mute;

  scarf_unstruct_initialize(nd, npts, x, y, NULL, ds, acf, cl, sigma, &options);

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
    printf("****** SRM algorithm, 2D, structured mesh *******\n");
  }

  scarf_opt_init(&options);

  options.hurst  = hurst;
  options.alpha  = alpha;
  options.npoi   = npoi;
  options.poi    = poi;
  options.taper  = taper;
  options.mute   = mute;
  options.solver = 1;

  scarf_struct_initialize(nd, fs, fe, ds, acf, cl, sigma, &options);

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
  scarf_io_one(nd, n, field, "srm_struct_whole_2d", nwriters);
  watch_stop(&tictoc);

  if (world_rank == 0) {
    printf("I/O time                      |%12.5f sec |\n", (float) tictoc);
    printf("*************************************************\n");
  }

  scarf_finalize();

  if (world_rank == 0){
    printf("\n");
    printf("*************************************************\n");
    printf("****** SRM algorithm, 2D, structured mesh *******\n");
  }

  scarf_opt_init(&options);

  options.hurst  = hurst;
  options.alpha  = alpha;
  options.npoi   = npoi;
  options.poi    = poi;
  options.taper  = taper;
  options.mute   = mute;
  options.solver = 1;

  scarf_unstruct_initialize(nd, npts, x, y, NULL, ds, acf, cl, sigma, &options);

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
  // create sample structured mesh: 3D case
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
    printf("****** FIM algorithm, 3D, structured mesh *******\n");
  }

  scarf_opt_init(&options);

  options.hurst  = hurst;
  options.alpha  = alpha;
  options.beta   = beta;
  options.npoi   = npoi;
  options.poi    = poi;
  options.taper  = taper;
  options.mute   = mute;

  scarf_struct_initialize(nd, fs, fe, ds, acf, cl, sigma, &options);

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
  scarf_io_slice(n, "x", n[0]/2, field, "fim_struct_xslice");
  scarf_io_slice(n, "y", n[1]/2, field, "fim_struct_yslice");
  scarf_io_slice(n, "z", n[2]/2, field, "fim_struct_zslice");
  watch_stop(&tictoc);

  if (world_rank == 0) {
    printf("I/O time (slices)             |%12.5f sec |\n", (float) tictoc);
  }

  watch_start(&tictoc);
  scarf_io_one(nd, n, field, "fim_struct_whole_3d", nwriters);
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
    printf("***** FIM algorithm, 3D, unstructured mesh ******\n");
  }

  scarf_opt_init(&options);

  options.hurst  = hurst;
  options.alpha  = alpha;
  options.beta   = beta;
  options.npoi   = npoi;
  options.poi    = poi;
  options.taper  = taper;
  options.mute   = mute;

  scarf_unstruct_initialize(nd, npts, x, y, z, ds, acf, cl, sigma, &options);

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
    printf("****** SRM algorithm, 3D, structured mesh *******\n");
  }

  scarf_opt_init(&options);

  options.hurst  = hurst;
  options.alpha  = alpha;
  options.beta   = beta;
  options.npoi   = npoi;
  options.poi    = poi;
  options.taper  = taper;
  options.mute   = mute;
  options.solver = 1;

  scarf_struct_initialize(nd, fs, fe, ds, acf, cl, sigma, &options);

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
  scarf_io_slice(n, "x", n[0]/2, field, "srm_struct_xslice");
  scarf_io_slice(n, "y", n[1]/2, field, "srm_struct_yslice");
  scarf_io_slice(n, "z", n[2]/2, field, "srm_struct_zslice");
  watch_stop(&tictoc);

  if (world_rank == 0) {
    printf("I/O time (slices)             |%12.5f sec    |\n", (float) tictoc);
  }

  watch_start(&tictoc);
  scarf_io_one(nd, n, field, "srm_struct_whole_3d", nwriters);
  watch_stop(&tictoc);

  if (world_rank == 0) {
    printf("I/O time                      |%12.5f sec    |\n", (float) tictoc);
    printf("*************************************************\n");
  }

  scarf_finalize();

  if (world_rank == 0){
    printf("\n");
    printf("*************************************************\n");
    printf("****** SRM algorithm, 3D, structured mesh *******\n");
  }

  scarf_opt_init(&options);

  options.hurst  = hurst;
  options.alpha  = alpha;
  options.beta   = beta;
  options.npoi   = npoi;
  options.poi    = poi;
  options.taper  = taper;
  options.mute   = mute;
  options.solver = 1;

  scarf_unstruct_initialize(nd, npts, x, y, z, ds, acf, cl, sigma, &options);

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
