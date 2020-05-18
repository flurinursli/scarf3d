#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "scarf3d.h"


// FORTRAN subroutine prototype
extern void sample_mesh(const int* rank, const int* ntasks, const int n[], int fs[], int fe[]);

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
  const int n[3] = {500, 500, 500};

  // grid-step
  const fpp ds = 50.;

  // autocorrelation function (0=von karman/exponential, 1=gaussian)
  const int acf = 0;

  // correlation length
  const fpp cl[3] = {2000., 2000., 2000.};

  // standard deviation
  const fpp sigma = 0.05;

  // seed number
  const int seed = 1235;

  // ===========================================================================
  // ---------------------------------------------------------------------------
  // Add extra parameters

  const fpp hurst = 0.25;

  // ===========================================================================
  // ---------------------------------------------------------------------------
  // create sample structured mesh

  int fs[3], fe[3];

  // domain decomposition
  sample_mesh(&world_rank, &world_size, n, fs, fe);

  int dims[3];

  for (int i = 0; i < 3; i++){
    dims[i] = (fe[i] - fs[i] + 1);
  }

  fpp* x     = (fpp*) malloc(dims[0]*dims[1]*dims[2] * sizeof(fpp));
  fpp* y     = (fpp*) malloc(dims[0]*dims[1]*dims[2] * sizeof(fpp));
  fpp* z     = (fpp*) malloc(dims[0]*dims[1]*dims[2] * sizeof(fpp));
  fpp* field = (fpp*) malloc(dims[0]*dims[1]*dims[2] * sizeof(fpp));

  long c;

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

  const int npts = dims[0] * dims[1] * dims[2];

  fpp stats[8];

  double tictoc;

  // ===========================================================================
  // ---------------------------------------------------------------------------
  // FFT method test

  struct scarf_opt options;

  scarf_opt_init(&options);

  options.hurst = hurst;

  scarf_struct_initialize(fs, fe, ds, acf, cl, sigma, &options);

  watch_start(&tictoc);
  scarf_execute(seed, field, stats);
  watch_stop(&tictoc);

  if (world_rank == 0){
    printf("\n");
    printf("Structured Mesh Test completed in: %f\n", (float) tictoc);
    printf("Domain too small?                : %d\n", (int) stats[0]);
    printf("Grid-step too large?             : %d\n", (int) stats[1]);
    printf("Standard deviation               : %f\n", (float) stats[2]);
    printf("Mean value                       : %f\n", (float) stats[3]);
    printf("Timing for spectrum              : %f\n", (float) stats[4]);
    printf("Timing for symmetry              : %f\n", (float) stats[5]);
    printf("Timing for IFFT                  : %f\n", (float) stats[6]);
    printf("Timing for interpolation         : %f\n", (float) stats[7]);
  }

  // IO
  watch_start(&tictoc);
  scarf_io_slice(n, "x", n[0]/2, field, "fft_struct_xslice");
  scarf_io_slice(n, "y", n[1]/2, field, "fft_struct_yslice");
  scarf_io_slice(n, "z", n[2]/2, field, "fft_struct_zslice");
  watch_stop(&tictoc);

  if (world_rank == 0) {
    printf("Slice(s) written in              : %f\n", (float) tictoc);
  }

  int nwriters[1] = {3};

  watch_start(&tictoc);
  scarf_io_one(n, field, "fft_struct_whole", nwriters);
  watch_stop(&tictoc);

  if (world_rank == 0) {
    printf("Whole file written in            : %f\n", (float) tictoc);
  }

  scarf_finalize();

  MPI_Barrier(MPI_COMM_WORLD);

  scarf_unstruct_initialize(npts, x, y, z, ds, acf, cl, sigma, &options);

  watch_start(&tictoc);
  scarf_execute(seed, field, stats);
  watch_stop(&tictoc);

  if (world_rank == 0){
    printf("\n");
    printf("Unstructured Mesh Test completed in: %f\n", (float) tictoc);
    printf("Domain too small?                  : %d\n", (int) stats[0]);
    printf("Grid-step too large?               : %d\n", (int) stats[1]);
    printf("Standard deviation                 : %f\n", (float) stats[2]);
    printf("Mean value                         : %f\n", (float) stats[3]);
    printf("Timing for spectrum                : %f\n", (float) stats[4]);
    printf("Timing for symmetry                : %f\n", (float) stats[5]);
    printf("Timing for IFFT                    : %f\n", (float) stats[6]);
    printf("Timing for interpolation           : %f\n", (float) stats[7]);
  }

  scarf_finalize();

  // ===========================================================================
  // ---------------------------------------------------------------------------
  // SPEC method test

#ifdef SPECTRAL

  scarf_opt_init(&options);

  options.solver = 1;

  scarf_struct_initialize(fs, fe, ds, acf, cl, sigma, &options);

  watch_start(&tictoc);
  scarf_execute(seed, field, stats);
  watch_stop(&tictoc);

  if (world_rank == 0){
    printf("\n");
    printf("Structured Mesh Test completed in: %f\n", (float) tictoc);
    printf("Domain too small?                : %d\n", (int) stats[0]);
    printf("Grid-step too large?             : %d\n", (int) stats[1]);
    printf("Standard deviation               : %f\n", (float) stats[2]);
    printf("Mean value                       : %f\n", (float) stats[3]);
    printf("Timing for CPU execution         : %f\n", (float) stats[4]);
    printf("Timing for GPU execution         : %f\n", (float) stats[5]);
  }

  // IO
  watch_start(&tictoc);
  scarf_io_slice(n, "x", n[0]/2, field, "spec_struct_xslice");
  scarf_io_slice(n, "y", n[1]/2, field, "spec_struct_yslice");
  scarf_io_slice(n, "z", n[2]/2, field, "spec_struct_zslice");
  watch_stop(&tictoc);

  if (world_rank == 0) {
    printf("Slice(s) written in              : %f\n", (float) tictoc);
  }

  watch_start(&tictoc);
  scarf_io_one(n, field, "spec_struct_whole", nwriters);
  watch_stop(&tictoc);

  if (world_rank == 0) {
    printf("Whole file written in            : %f\n", (float) tictoc);
  }

  scarf_finalize();

  scarf_unstruct_initialize(npts, x, y, z, ds, acf, cl, sigma, &options);

  watch_start(&tictoc);
  scarf_execute(seed, field, stats);
  watch_stop(&tictoc);

  if (world_rank == 0){
    printf("\n");
    printf("Unstructured Mesh Test completed in: %f\n", (float) tictoc);
    printf("Domain too small?                  : %d\n", (int) stats[0]);
    printf("Grid-step too large?               : %d\n", (int) stats[1]);
    printf("Standard deviation                 : %f\n", (float) stats[2]);
    printf("Mean value                         : %f\n", (float) stats[3]);
    printf("Timing for CPU execution           : %f\n", (float) stats[4]);
    printf("Timing for GPU execution           : %f\n", (float) stats[5]);
  }

  scarf_finalize();

#endif

   free(x);
   free(y);
   free(z);
   free(field);

   MPI_Barrier(MPI_COMM_WORLD);

   MPI_Finalize();

  return 0;

}
