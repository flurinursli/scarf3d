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

  MPI_Barrier(MPI_COMM_WORLD);

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

  // hurst exponent
  const real hurst = 0.25;

  const int npoi = 2;

  real poi[4] = {400*ds, 250*ds, 200*ds, 150*ds};

  // I/O writers (only for PFS)
  const int nwriters[1] = {2};

  // other variables used throughout the code
  double tictoc;
  int fs[3], fe[3], dims[3];
  int npts, nd;
  long c;
  real stats[8];
  real *x, *y, *z, *field;

  // set formatting
  std::cout << std::right << std::fixed << std::boolalpha;
  std::cout.precision(5);

  // ================================================================================================================================
  // --------------------------------------------------------------------------------------------------------------------------------
  // create sample structured mesh: 2D case
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

  // ================================================================================================================================
  // --------------------------------------------------------------------------------------------------------------------------------
  // tests FIM algorithm
  // --------------------------------------------------------------------------------------------------------------------------------
  // ================================================================================================================================

  if (world_rank == 0){
    std::cout << ""                                                  << std::endl;
    std::cout << "*************************************************" << std::endl;
    std::cout << "****** FIM algorithm, 2D, structured mesh *******" << std::endl;
  }

  {

    Scarf3D::options.hurst = hurst;
    Scarf3D::options.alpha = 30;
    Scarf3D::options.npoi  = npoi;
    Scarf3D::options.poi   = poi;

    Scarf3D::Initialize<fft> S(nd, fs, fe, ds, acf, cl, sigma);

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
    S.io(nd, n, field, "fim_struct_whole_2d", nwriters);
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
    std::cout << "***** FIM algorithm, 2D, unstructured mesh ******" << std::endl;
  }

  {

    Scarf3D::options.hurst = hurst;

    Scarf3D::Initialize<fft> S(nd, npts, x, y, ds, acf, cl, sigma);

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

    // I/O not possible for unstructured mesh because "fs" and "fe" are not defined

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
    std::cout << "****** SRM algorithm, 2D, structured mesh *******" << std::endl;
  }

  {

    Scarf3D::options.hurst = hurst;

    Scarf3D::options.alpha = 20;

    Scarf3D::Initialize<spec> S(nd, fs, fe, ds, acf, cl, sigma);

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
    S.io(nd, n, field, "srm_struct_whole_2d", nwriters);
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
      std::cout << "***** SRM algorithm, 2D, unstructured mesh ******" << std::endl;
    }

    {

      Scarf3D::options.hurst = hurst;
      Scarf3D::options.alpha = 50;

      Scarf3D::Initialize<spec> S(nd, npts, x, y, z, ds, acf, cl, sigma);

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

  // ================================================================================================================================
  // --------------------------------------------------------------------------------------------------------------------------------
  // create sample structured mesh: 3D case
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

  // ================================================================================================================================
  // --------------------------------------------------------------------------------------------------------------------------------
  // tests FIM algorithm
  // --------------------------------------------------------------------------------------------------------------------------------
  // ================================================================================================================================

  if (world_rank == 0){
    std::cout << ""                                                  << std::endl;
    std::cout << "*************************************************" << std::endl;
    std::cout << "****** FIM algorithm, 3D, structured mesh *******" << std::endl;
  }

  {

    Scarf3D::options.hurst = hurst;
    Scarf3D::options.alpha = 30;
    Scarf3D::options.beta  = 20;

    Scarf3D::Initialize<fft> S(nd, fs, fe, ds, acf, cl, sigma);

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
    S.io(n, "x", n[0]/2, field, "fim_struct_xslice");
    S.io(n, "y", n[1]/2, field, "fim_struct_yslice");
    S.io(n, "z", n[2]/2, field, "fim_struct_zslice");
    watch_stop(&tictoc);

    if (world_rank == 0) {
       std::cout << "I/O time (slices)             |" << std::setw(12) << (float) tictoc << " sec" << " |" << std::endl;
    }

    watch_start(&tictoc);
    S.io(nd, n, field, "fim_struct_whole_3d", nwriters);
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
    std::cout << "***** FIM algorithm, 3D, unstructured mesh ******" << std::endl;
  }

  {

    Scarf3D::options.hurst = hurst;

    Scarf3D::Initialize<fft> S(nd, npts, x, y, z, ds, acf, cl, sigma);

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
    std::cout << "****** SRM algorithm, 3D, structured mesh *******" << std::endl;
  }

  {

    Scarf3D::options.hurst = hurst;

    Scarf3D::options.alpha = 20;

    Scarf3D::Initialize<spec> S(nd, fs, fe, ds, acf, cl, sigma);

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
    S.io(n, "x", n[0]/2, field, "spec_struct_xslice");
    S.io(n, "y", n[1]/2, field, "spec_struct_yslice");
    S.io(n, "z", n[2]/2, field, "spec_struct_zslice");
    watch_stop(&tictoc);

    if (world_rank == 0) {
       std::cout << "I/O time (slices)             |" << std::setw(12) << (float) tictoc << " sec" << "  |" << std::endl;
    }

    watch_start(&tictoc);
    S.io(nd, n, field, "srm_struct_whole", nwriters);
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
    std::cout << "***** SRM algorithm, 3D, unstructured mesh ******" << std::endl;
  }

  {

    Scarf3D::options.hurst = hurst;

    Scarf3D::options.alpha = 50;
    Scarf3D::options.gamma = 10;

    Scarf3D::Initialize<spec> S(nd, npts, x, y, z, ds, acf, cl, sigma);

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

   MPI_Barrier(MPI_COMM_WORLD);

   MPI_Finalize();

   return 0;

}
