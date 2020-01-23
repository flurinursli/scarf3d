#include <iostream>
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>

extern "C"
{
#ifdef DOUBLE_PREC
double* prng(int seed, int *ls, int *le, int *npts);
void free_mem(double *p);
#else
  float* prng(int seed, int *ls, int *le, int *npts);
  void free_mem(float *p);
#endif
}

// -----------------------------------------------------------------------------

#ifdef DOUBLE_PREC
double* prng(int seed, int *ls, int *le, int *npts)
#else
float* prng(int seed, int *ls, int *le, int *npts)
#endif
{

  trng::yarn2 r;
  trng::uniform01_dist<> u;

  int n = (le[0] - ls[0] + 1) * (le[1] - ls[1] + 1) * (le[2] - ls[2] + 1);

  // allocate array for random numbers
#ifdef DOUBLE_PREC
  double* x = (double*) malloc(sizeof(double) * n);
#else
  float* x = (float*) malloc(sizeof(float) * n);
#endif

  unsigned long myseed = seed;

  // initialise with desired seed
  r.seed(myseed);

  // initialise counter
  int c = -1;
  int s;

  std::cout << "input " << seed << ' ' << ls[0] << ' ' << le[0] << std::endl;

  for (int k = ls[2]; k <= le[2]; ++k){

    // skip z-levels
    s = (k - 1) * npts[0] * npts[1];
    r.jump(s);

    for (int j = ls[1]; j <= le[1]; ++j){

      // skip points in the XY plane until "ls[0]"
      s = (j - 1) * npts[0] + ls[0] - 1;
      r.jump(s);

      for (int i = ls[0]; i <= le[0]; ++i){
        c = c + 1;
        x[c] = u(r);
      }

      // skip remaning points in the XY plane
      s = (npts[0] - le[0]) + (npts[1] - j) * npts[0];
      r.jump(s);

    }
  }

  return x;

}

// -----------------------------------------------------------------------------

#ifdef DOUBLE_PREC
void free_mem(double *p)
#else
void free_mem(float *p)
#endif
{
  free(p);
}
