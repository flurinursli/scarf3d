#include <mpi.h>
#include "scarf3d.h"

#ifdef DOUBLE_PREC
  typedef double fpp;
  //using fpp = double;
#else
  typedef float fpp;
  //using fpp = float;
#endif

void watch_start(double* t){

  MPI_Barrier(MPI_COMM_WORLD);

  *t = MPI_Wtime();

}

void watch_stop(double* t){

  MPI_Barrier(MPI_COMM_WORLD);

  *t = MPI_Wtime() - *t;

}

int main(){

  // initialise MPI
  MPI_init(nullptr, nullptr);

  // get number of tasks
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // get rank of each task
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  MPI_Barrier(MPI_COMM_WORLD);

  // Set mandatory parameters

  // number of grid points in model
  const int n = {500, 500, 500};

  // grid-step
  const fpp ds = 100.;


  // ===========================================================================
  // ---------------------------------------------------------------------------
  // create sample structured mesh
  // ---------------------------------------------------------------------------
  // ===========================================================================

  int fs[3], fe[3];

  // domain decomposition
  create_mesh(world_rank, world_size, n, fs, fe);

  int dims[3];

  for (int i = 0; i < 3; i++){
    int dims[i] = (fe[0] - fs[0] + 1);
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

  fpp *stats = new fpp[8];

  // autocorrelation function (0=von karman/exponential, 1=gaussian)
  const int acf = 0;

  // correlation length
  const fpp cl = {500., 500., 500.};

  // standard deviation
  const fpp sigma = 0.05;

  // seed number
  const int seed = 1235;

  double tictoc;


  // ===========================================================================
  // ---------------------------------------------------------------------------
  // FFT method test
  // ---------------------------------------------------------------------------
  // ===========================================================================

  Scarf3D::Initialize<fft> S(fs, fe, ds, acf, cl, sigma, hurst = 0.5);

  watch_start(&tictoc);
  S.execute(seed, f, stats);
  watch_stop(&tictoc);

  if (world_rank == 0){
    std::cout << "" << std:endl;
    std::cout << "Structured Mesh Test completed in: " << (float) tictoc   << std::endl;
    std::cout << "Domain too small?                : " << (int) stats[1]   << std::endl;
    std::cout << "Grid-step too large?             : " << (int) stats[2]   << std::endl;
    std::cout << "Standard deviation               : " << (float) stats[3] << std::endl;
    std::cout << "Mean value                       : " << (float) stats[4] << std::endl;
    std::cout << "Timing for spectrum              : " << (float) stats[5] << std::endl;
    std::cout << "Timing for symmetry              : " << (float) stats[6] << std::endl;
    std::cout << "Timing for IFFT                  : " << (float) stats[7] << std::endl;
    std::cout << "Timing for interpolation         : " << (float) stats[8] << std::endl;
  }


  // IO

  // call destructor explicitly
  S.~Initialize();

  Scarf3D::Initialize<fft> S(x, y, z, dh, acf, cl, sigma, hurst = 0.5);

  watch_start(&tictoc);
  S.execute(seed, field, stats);
  watch_stop(&tictoc);

  if (world_rank == 0){
    std::cout << "" << std:endl;
    std::cout << "Unstructured Mesh Test completed in: " << (float) tictoc   << std::endl;
    std::cout << "Domain too small?                  : " << (int) stats[1]   << std::endl;
    std::cout << "Grid-step too large?               : " << (int) stats[2]   << std::endl;
    std::cout << "Standard deviation                 : " << (float) stats[3] << std::endl;
    std::cout << "Mean value                         : " << (float) stats[4] << std::endl;
    std::cout << "Timing for spectrum                : " << (float) stats[5] << std::endl;
    std::cout << "Timing for symmetry                : " << (float) stats[6] << std::endl;
    std::cout << "Timing for IFFT                    : " << (float) stats[7] << std::endl;
    std::cout << "Timing for interpolation           : " << (float) stats[8] << std::endl;
  }

  // IO


  // call destructor explicitly
  S.~Initialize();

  // ===========================================================================
  // ---------------------------------------------------------------------------
  // SPEC method test
  // ---------------------------------------------------------------------------
  // ===========================================================================

#ifdef SPECTRAL

  Scarf3D::Initialize<spec> S(fs, fe, ds, acf, cl, sigma, hurst = 0.5);

  watch_start(&tictoc);
  S.execute(seed, field, stats);
  watch_stop(&tictoc);

  if (world_rank == 0){
    std::cout << "" << std:endl;
    std::cout << "Structured Mesh Test completed in: " << (float) tictoc   << std::endl;
    std::cout << "Domain too small?                : " << (int) stats[1]   << std::endl;
    std::cout << "Grid-step too large?             : " << (int) stats[2]   << std::endl;
    std::cout << "Standard deviation               : " << (float) stats[3] << std::endl;
    std::cout << "Mean value                       : " << (float) stats[4] << std::endl;
    std::cout << "Timing for CPU execution         : " << (float) stats[5] << std::endl;
    std::cout << "Timing for GPU execution         : " << (float) stats[6] << std::endl;
  }


  // IO

  // call destructor explicitly
  S.~Initialize();

  Scarf3D::Initialize<spec> S(x, y, z, dh, acf, cl, sigma, hurst = 0.5);

  watch_start(&tictoc);
  S.execute(seed, field, stats);
  watch_stop(&tictoc);

  if (world_rank == 0){
    std::cout << "" << std:endl;
    std::cout << "Unstructured Mesh Test completed in: " << (float) tictoc   << std::endl;
    std::cout << "Domain too small?                  : " << (int) stats[1]   << std::endl;
    std::cout << "Grid-step too large?               : " << (int) stats[2]   << std::endl;
    std::cout << "Standard deviation                 : " << (float) stats[3] << std::endl;
    std::cout << "Mean value                         : " << (float) stats[4] << std::endl;
    std::cout << "Timing for CPU execution           : " << (float) stats[5] << std::endl;
    std::cout << "Timing for GPU execution           : " << (float) stats[6] << std::endl;
  }

  // call destructor explicitly
  S.~Initialize();

#endif

   delete[] field;
   delete[] stats;

   MPI_Finalize();

   return 0;

}
