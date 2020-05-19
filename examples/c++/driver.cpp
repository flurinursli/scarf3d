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
  const int n[3] = {500, 500, 500};

  // grid-step
  const real ds = 50.;

  // autocorrelation function (0=von karman/exponential, 1=gaussian)
  const int acf = 0;

  // correlation length
  const real cl[3] = {2000., 2000., 2000.};

  // standard deviation
  const real sigma = 0.05;

  // seed number
  const int seed = 1235;

  const real hurst = 0.25;

  real stats[8];

  // real x0[3], x1[3];
  //
  // for (int i = 0; i < 3; i++){
  //   x0[i] = 0;
  //   x1[i] = (n[i] - 1) * ds;
  // }
  //
  // real dh = ds;
  //
  // for (int i = 0; i < 3; i++) n[i] = n[i] / 4;
  // ds = ds * 4;

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

  real* x     = new real[dims[0]*dims[1]*dims[2]];
  real* y     = new real[dims[0]*dims[1]*dims[2]];
  real* z     = new real[dims[0]*dims[1]*dims[2]];
  real* field = new real[dims[0]*dims[1]*dims[2]];

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

  double tictoc;


  // ===========================================================================
  // ---------------------------------------------------------------------------
  // FFT method test
  // ---------------------------------------------------------------------------
  // ===========================================================================

  if (world_rank == 0){
    std::cout << ""                                                             << std::endl;
    std::cout << "************************************************************" << std::endl;
    std::cout << "************************ FIM method ************************" << std::endl;
  }

  {

    Scarf3D::options.hurst = hurst;
    //Scarf3D::options.dh    = dh;

    //Scarf3D::options.nc = x0;
    //Scarf3D::options.fc = x1;

    Scarf3D::Initialize<fft> S(fs, fe, ds, acf, cl, sigma);

    watch_start(&tictoc);
    S.execute(seed, field, stats);
    watch_stop(&tictoc);

    if (world_rank == 0){
      std::cout << ""                                                                << std::endl;
      std::cout << "Summary structured mesh"                                         << std::endl;
      std::cout << "  i)   test completed in       : " << (float) tictoc   << " sec" << std::endl;
      std::cout << "  ii)  domain too small?       : " << (int) stats[0]             << std::endl;
      std::cout << "  iii) grid-step too large?    : " << (int) stats[1]             << std::endl;
      std::cout << "  iv)  standard deviation      : " << (float) stats[2]           << std::endl;
      std::cout << "  v)   mean value              : " << (float) stats[3]           << std::endl;
      std::cout << "  vi)  timing for spectrum     : " << (float) stats[4] << " sec" << std::endl;
      std::cout << "  vii) timing for symmetry     : " << (float) stats[5] << " sec" << std::endl;
      std::cout << "  viii)timing for ifft         : " << (float) stats[6] << " sec" << std::endl;
      std::cout << "  ix)  timing for interpolation: " << (float) stats[7] << " sec" << std::endl;
    }

    // IO
    watch_start(&tictoc);
    S.io(n, "x", n[0]/2, field, "fft_struct_xslice");
    S.io(n, "y", n[1]/2, field, "fft_struct_yslice");
    S.io(n, "z", n[2]/2, field, "fft_struct_zslice");
    watch_stop(&tictoc);

    if (world_rank == 0) {
       std::cout << "  x)   slice(s) written in     : " << (float) tictoc << " sec" << std::endl;
    }

    int nwriters[1] = {3};

    watch_start(&tictoc);
    S.io(n, field, "fft_struct_whole", nwriters);
    watch_stop(&tictoc);

    if (world_rank == 0) {
      std::cout << "  xi)  whole file written in   : " << (float) tictoc << " sec" << std::endl;
    }

    // call destructor explicitly
    //S.~Initialize();
  }

  {

    Scarf3D::options.hurst = hurst;
    //Scarf3D::options.ds    = ds;

    //Scarf3D::options.nc = x0;
    //Scarf3D::options.fc = x1;

    Scarf3D::Initialize<fft> S(npts, x, y, z, ds, acf, cl, sigma);

    watch_start(&tictoc);
    S.execute(seed, field, stats);
    watch_stop(&tictoc);

    if (world_rank == 0){
      std::cout << ""                                                                << std::endl;
      std::cout << "Summary unstructured mesh"                                       << std::endl;
      std::cout << "  i)   test completed in       : " << (float) tictoc   << " sec" << std::endl;
      std::cout << "  ii)  domain too small?       : " << (int) stats[0]             << std::endl;
      std::cout << "  iii) grid-step too large?    : " << (int) stats[1]             << std::endl;
      std::cout << "  iv)  standard deviation      : " << (float) stats[2]           << std::endl;
      std::cout << "  v)   mean value              : " << (float) stats[3]           << std::endl;
      std::cout << "  vi)  timing for spectrum     : " << (float) stats[4] << " sec" << std::endl;
      std::cout << "  vii) timing for symmetry     : " << (float) stats[5] << " sec" << std::endl;
      std::cout << "  viii)timing for ifft         : " << (float) stats[6] << " sec" << std::endl;
      std::cout << "  ix)  timing for interpolation: " << (float) stats[7] << " sec" << std::endl;
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

  if (world_rank == 0){
    std::cout << ""                                                             << std::endl;
    std::cout << "************************************************************" << std::endl;
    std::cout << "************************ SRM method ************************" << std::endl;
  }

  {

    Scarf3D::options.hurst = hurst;

    Scarf3D::Initialize<spec> S(fs, fe, ds, acf, cl, sigma);

    watch_start(&tictoc);
    S.execute(seed, field, stats);
    watch_stop(&tictoc);

    if (world_rank == 0){
      std::cout << ""                                                             << std::endl;
      std::cout << "Summary structured mesh"                                      << std::endl;
      std::cout << "  i)   test completed in    : " << (float) tictoc   << " sec" << std::endl;
      std::cout << "  ii)  domain too small?    : " << (int) stats[0]             << std::endl;
      std::cout << "  iii) grid-step too large? : " << (int) stats[1]             << std::endl;
      std::cout << "  iv)  standard deviation   : " << (float) stats[2]           << std::endl;
      std::cout << "  v)   mean value           : " << (float) stats[3]           << std::endl;
      std::cout << "  vi)  CPU main loop        : " << (float) stats[4] << " sec" << std::endl;
      std::cout << "  vii) GPU main loop        : " << (float) stats[5] << " sec" << std::endl;
    }

    // IO
    watch_start(&tictoc);
    S.io(n, "x", n[0]/2, field, "spec_struct_xslice");
    S.io(n, "y", n[1]/2, field, "spec_struct_yslice");
    S.io(n, "z", n[2]/2, field, "spec_struct_zslice");
    watch_stop(&tictoc);

    if (world_rank == 0) {
      std::cout << "  viii)slice(s) written in  : " << (float) tictoc << " sec" << std::endl;
    }

    int nwriters[1] = {3};

    watch_start(&tictoc);
    S.io(n, field, "spec_struct_whole", nwriters);
    watch_stop(&tictoc);

    if (world_rank == 0) {
      std::cout << "  ix)  whole file written in: " << (float) tictoc << " sec" << std::endl;
    }

    // call destructor explicitly
    // S.~Initialize();
  }

  {

    Scarf3D::options.hurst = hurst;

    Scarf3D::Initialize<spec> S(npts, x, y, z, ds, acf, cl, sigma);

    watch_start(&tictoc);
    S.execute(seed, field, stats);
    watch_stop(&tictoc);

    if (world_rank == 0){
      std::cout << ""                                                            << std::endl;
      std::cout << "Summary unstructured mesh"                                   << std::endl;
      std::cout << "  i)   test completed in   : " << (float) tictoc   << " sec" << std::endl;
      std::cout << "  ii)  domain too small?   : " << (int) stats[0]             << std::endl;
      std::cout << "  iii) grid-step too large?: " << (int) stats[1]             << std::endl;
      std::cout << "  iv)  standard deviation  : " << (float) stats[2]           << std::endl;
      std::cout << "  v)   mean value          : " << (float) stats[3]           << std::endl;
      std::cout << "  vi)  CPU main loop       : " << (float) stats[4] << " sec" << std::endl;
      std::cout << "  vii) GPU main loop       : " << (float) stats[5] << " sec" << std::endl;
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
