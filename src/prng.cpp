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
   To provide an interface between the SCARF3D and TRNG4 libraries

 Revisions:
     Date                    Description of change
     ====                    =====================
   04/05/20                  original version
*/

#include <iostream>
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
#include <trng/normal_dist.hpp>

// g++ -c prng.cpp -Ipath-to-trng-include-folder -O3 -cpp -DDOUBLE_PREC

extern "C"
{
#ifdef DOUBLE_PREC
  double* prng(int seed, int* ls, int* le, int* npts);
  void srng(double* x);
  void normdist(double* x);
  void free_mem(double* p);
#else
  float* prng(int seed, int* ls, int* le, int* npts);
  void srng(float* x);
  void normdist(float* x);
  void free_mem(float* p);
#endif
void set_seed(int seed);
}

// declare object "r" as static because it gets initialised with seed number only once.
// this object is used in "srng" and "normdist".
static trng::yarn2 r;

// -----------------------------------------------------------------------------

// return a sequence of uniformly distributed random numbers in the interval (0, 1) based on
// the block splitting approach
#ifdef DOUBLE_PREC
double* prng(int seed, int* ls, int* le, int* npts)
#else
float* prng(int seed, int* ls, int* le, int* npts)
#endif
{

  trng::yarn2 r;
  trng::uniform01_dist<> u;

  unsigned long long n = (le[0] - ls[0] + 1) * (le[1] - ls[1] + 1) * (le[2] - ls[2] + 1);

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
  long long c = -1;
  unsigned long long s;

  // skip previous z-levels
  s = (ls[2] - 1) * npts[0];
  s = s * npts[1];
  r.jump(s);

  for (int k = ls[2]; k <= le[2]; ++k){

    // skip further for previous y-levels
    s = (ls[1] - 1) * npts[0];
    r.jump(s);

    for (int j = ls[1]; j <= le[1]; ++j){

      // skip first points along x
      s = ls[0] - 1;
      r.jump(s);

      for (int i = ls[0]; i <= le[0]; ++i){
        c = c + 1;
        x[c] = u(r);
      }

      // skip remaning points along x
      s = (npts[0] - le[0]);
      r.jump(s);

    }

    // skip remaning points in XY plane
    s = (npts[1] - le[1]) * npts[0];
    r.jump(s);

  }

  return x;

}

// -----------------------------------------------------------------------------

// deallocate memory allocated in function "prng"
#ifdef DOUBLE_PREC
void free_mem(double *p)
#else
void free_mem(float *p)
#endif
{
  free(p);
}

// -----------------------------------------------------------------------------

// set seed number for random generator
void set_seed(int seed){unsigned long myseed = seed; r.seed(myseed);}

// -----------------------------------------------------------------------------

// thi subroutine returns two random numbers uniformly distributed in the interval (0, 1)
#ifdef DOUBLE_PREC
void srng(double* x)
#else
void srng(float* x)
#endif
{

  trng::uniform01_dist<> u;

  for (int i = 0; i < 2; i++){
    x[i] = u(r);
  }

}

// -----------------------------------------------------------------------------

// thi subroutine returns two random numbers normally distributed around 0 with std.dev=1
#ifdef DOUBLE_PREC
void normdist(double* x)
#else
void normdist(float* x)
#endif
{

  trng::normal_dist<> u(0., 1.);

  for (int i = 0; i < 2; i++){
    x[i] = u(r);
  }

}
