#include <mpi.h>
#include <iostream>
#include "scarf3d.h"

// FORTRAN subroutine prototype
extern "C"{
  extern void sample_mesh(const int* rank, const int* ntasks, const int n[], int fs[], int fe[]);
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

  MPI_Barrier(MPI_COMM_WORLD);

  // Set mandatory parameters

  // number of grid points in model
  const int n[3] = {501, 501, 501};

  // grid-step
  const fpp ds = 100.;


  // ===========================================================================
  // ---------------------------------------------------------------------------
  // create sample structured mesh
  // ---------------------------------------------------------------------------
  // ===========================================================================

  int fs[3], fe[3];

  // domain decomposition
  sample_mesh(&world_rank, &world_size, n, fs, fe);

  int dims[3];

  for (int i = 0; i < 3; i++){
    dims[i] = (fe[i] - fs[i] + 1);
  }

  fpp* x     = new fpp[dims[0]*dims[1]*dims[2]];
  fpp* y     = new fpp[dims[0]*dims[1]*dims[2]];
  fpp* z     = new fpp[dims[0]*dims[1]*dims[2]];
  fpp* field = new fpp[dims[0]*dims[1]*dims[2]];

  long c;

  for (int k = 0; k < dims[2]; k++){
    for (int j = 0; j < dims[1]; j++){
      for (int i = 0; i < dims[0]; i++){
        c    = (k * dims[0] * dims[1]) + (j * dims[0]) + i;
        x[c] = (i + fs[0] - 0.5) * ds;
        y[c] = (j + fs[1] - 0.5) * ds;
        z[c] = (k + fs[2] - 0.5) * ds;
      }
    }
  }

  const int npts = dims[0] * dims[1] * dims[2];

  fpp stats[8];

  // autocorrelation function (0=von karman/exponential, 1=gaussian)
  const int acf = 0;

  // correlation length
  const fpp cl[3] = {5000., 5000., 5000.};

  // standard deviation
  const fpp sigma = 0.05;

  // seed number
  const int seed = 1235;

  const fpp hurst[1] = {0.2};

  double tictoc;


  // ===========================================================================
  // ---------------------------------------------------------------------------
  // FFT method test
  // ---------------------------------------------------------------------------
  // ===========================================================================

  {

    Scarf3D::hurst = 0.1;

    Scarf3D::Initialize<fft> S(fs, fe, ds, acf, cl, sigma);

    watch_start(&tictoc);
    S.execute(seed, field, stats);
    watch_stop(&tictoc);

    if (world_rank == 0){
      std::cout << "" << std::endl;
      std::cout << "Structured Mesh Test completed in: " << (float) tictoc   << std::endl;
      std::cout << "Domain too small?                : " << (int) stats[0]   << std::endl;
      std::cout << "Grid-step too large?             : " << (int) stats[1]   << std::endl;
      std::cout << "Standard deviation               : " << (float) stats[2] << std::endl;
      std::cout << "Mean value                       : " << (float) stats[3] << std::endl;
      std::cout << "Timing for spectrum              : " << (float) stats[4] << std::endl;
      std::cout << "Timing for symmetry              : " << (float) stats[5] << std::endl;
      std::cout << "Timing for IFFT                  : " << (float) stats[6] << std::endl;
      std::cout << "Timing for interpolation         : " << (float) stats[7] << std::endl;
    }

    // IO
    watch_start(&tictoc);
    S.io(n, 1, 400, field, "fft_struct_xslice");
    S.io(n, 2, 250, field, "fft_struct_yslice");
    S.io(n, 3, 100, field, "fft_struct_zslice");
    watch_stop(&tictoc);

    if (world_rank == 0) {
      std::cout << "Slice(s) written in              : " << (float) tictoc << std::endl;
    }

    int nwriters[1] = {3};

    watch_start(&tictoc);
    S.io(n, field, "fft_struct_whole", nwriters);
    watch_stop(&tictoc);

    if (world_rank == 0) {
      std::cout << "Whole file written in            : " << (float) tictoc << std::endl;
    }

    // call destructor explicitly
    //S.~Initialize();
  }

  {

    Scarf3D::hurst = 0.1;

    Scarf3D::Initialize<fft> S(npts, x, y, z, ds, acf, cl, sigma);

    watch_start(&tictoc);
    S.execute(seed, field, stats);
    watch_stop(&tictoc);

    if (world_rank == 0){
      std::cout << "" << std::endl;
      std::cout << "Unstructured Mesh Test completed in: " << (float) tictoc   << std::endl;
      std::cout << "Domain too small?                  : " << (int) stats[0]   << std::endl;
      std::cout << "Grid-step too large?               : " << (int) stats[1]   << std::endl;
      std::cout << "Standard deviation                 : " << (float) stats[2] << std::endl;
      std::cout << "Mean value                         : " << (float) stats[3] << std::endl;
      std::cout << "Timing for spectrum                : " << (float) stats[4] << std::endl;
      std::cout << "Timing for symmetry                : " << (float) stats[5] << std::endl;
      std::cout << "Timing for IFFT                    : " << (float) stats[6] << std::endl;
      std::cout << "Timing for interpolation           : " << (float) stats[7] << std::endl;
    }

    // call destructor explicitly
    //S.~Initialize();
  }

  // ===========================================================================
  // ---------------------------------------------------------------------------
  // SPEC method test
  // ---------------------------------------------------------------------------
  // ===========================================================================

#ifdef SPECTRAL

  {

    Scarf3D::hurst = 0.1;

    Scarf3D::Initialize<spec> S(fs, fe, ds, acf, cl, sigma);

    watch_start(&tictoc);
    S.execute(seed, field, stats);
    watch_stop(&tictoc);

    if (world_rank == 0){
      std::cout << "" << std::endl;
      std::cout << "Structured Mesh Test completed in: " << (float) tictoc   << std::endl;
      std::cout << "Domain too small?                : " << (int) stats[0]   << std::endl;
      std::cout << "Grid-step too large?             : " << (int) stats[1]   << std::endl;
      std::cout << "Standard deviation               : " << (float) stats[2] << std::endl;
      std::cout << "Mean value                       : " << (float) stats[3] << std::endl;
      std::cout << "Timing for CPU execution         : " << (float) stats[4] << std::endl;
      std::cout << "Timing for GPU execution         : " << (float) stats[5] << std::endl;
    }

    // IO
    watch_start(&tictoc);
    S.io(n, 1, 400, field, "spec_struct_xslice");
    S.io(n, 2, 250, field, "spec_struct_yslice");
    S.io(n, 3, 100, field, "spec_struct_zslice");
    watch_stop(&tictoc);

    if (world_rank == 0) {
      std::cout << "Slice(s) written in              : " << (float) tictoc << std::endl;
    }

    int nwriters[1] = {3};

    watch_start(&tictoc);
    S.io(n, field, "spec_struct_whole", nwriters);
    watch_stop(&tictoc);

    if (world_rank == 0) {
      std::cout << "Whole file written in            : " << (float) tictoc << std::endl;
    }

    // call destructor explicitly
    // S.~Initialize();
  }

  {

    Scarf3D::hurst = 0.1;

    Scarf3D::Initialize<spec> S(npts, x, y, z, ds, acf, cl, sigma);

    watch_start(&tictoc);
    S.execute(seed, field, stats);
    watch_stop(&tictoc);

    if (world_rank == 0){
      std::cout << "" << std::endl;
      std::cout << "Unstructured Mesh Test completed in: " << (float) tictoc   << std::endl;
      std::cout << "Domain too small?                  : " << (int) stats[0]   << std::endl;
      std::cout << "Grid-step too large?               : " << (int) stats[1]   << std::endl;
      std::cout << "Standard deviation                 : " << (float) stats[2] << std::endl;
      std::cout << "Mean value                         : " << (float) stats[3] << std::endl;
      std::cout << "Timing for CPU execution           : " << (float) stats[4] << std::endl;
      std::cout << "Timing for GPU execution           : " << (float) stats[5] << std::endl;
    }

    // call destructor explicitly
    S.~Initialize();
  }

#endif

   delete[] x;
   delete[] y;
   delete[] z;
   delete[] field;

   MPI_Barrier(MPI_COMM_WORLD);

   MPI_Finalize();

   return 0;

}
