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
  To provide an header file for the SCARF3D library in C/C++

Revisions:
    Date                    Description of change
    ====                    =====================
  04/05/20                  original version

*/

#ifndef SCARF3D_H_

#define SCARF3D_H_

#ifdef DOUBLE_PREC
typedef double real;
//using real = double;
#else
typedef float real;
//using real = float;
#endif

#ifdef __cplusplus
extern "C"{
#endif

enum algorithm {fim, srm};

struct scarf_opt{
  int method;
  real hurst;
  real dh;
  real ds;
  real * poi;
  int npoi;
  real mute;
  real taper;
  int rescale;
  int pad;
  real * nc;
  real * fc;
  real alpha;
  real beta;
  real gamma;
};

// declare C functions
void scarf_opt_init(struct scarf_opt * var);

void scarf_cart_initialize(const int nd, const int fs[], const int fe[], const real ds, const int acf, const real cl[], const real sigma, struct scarf_opt *var);

void scarf_nocart_initialize(const int nd, const int npts, const real* x, const real* y, const real* z, const real dh, const int acf, const real cl[], const real sigma, struct scarf_opt *var);

void scarf_execute(const int seed, real* field, real stats[]);

void scarf_finalize();

void scarf_io_one(const int nd, const int* npts, const real* field, const char fname[], const int* nwriters);

void scarf_io_slice(const int* npts, const char axis[], const int plane, const real* field, const char fname[]);

// define C++ class
#ifdef __cplusplus
}

namespace Scarf3D
{

  scarf_opt options {0, 0, 0, 0, nullptr, 0, 0, 0, 0, 0, nullptr, nullptr, 0, 0, 0};

  template <algorithm method> class Initialize
  {

    private:

    public:

      // constructor structured mesh
      Initialize(const int nd, const int fs[], const int fe[], const real ds, const int acf, const real cl[], const real sigma)
      {

        if (method == srm) options.method = 1;

        scarf_cart_initialize(nd, fs, fe, ds, acf, cl, sigma, &options);

      };

      // constructor unstructured mesh
      Initialize(const int nd, const int npts, const real* x, const real* y, const real* z, const real dh, const int acf, const real cl[], const real sigma)
      {

        if (method == srm) options.method = 1;

        scarf_nocart_initialize(nd, npts, x, y, z, dh, acf, cl, sigma, &options);

      };

      Initialize(const int nd, const int npts, const real* x, const real* y, const real dh, const int acf, const real cl[], const real sigma)
      {

        if (method == srm) options.method = 1;

        scarf_nocart_initialize(nd, npts, x, y, nullptr, dh, acf, cl, sigma, &options);

      };

      // destructor
      ~Initialize()
      {
        scarf_finalize();
        options = {0, 0, 0, 0, nullptr, 0, 0, 0, 0, 0, nullptr, nullptr, 0, 0, 0};
      };

      // simulate random field
      void execute(const int seed, real* field, real stats[])
      {
        scarf_execute(seed, field, stats);
      };

      // IO whole model
      void io(const int nd, const int npts[], const real* field, const char fname[], const int* nwriters = nullptr)
      {
        scarf_io_one(nd, npts, field, fname, nwriters);
      };

      // IO model slice
      void io(const int npts[], const char axis[], const int plane, const real* field, const char fname[])
      {
        scarf_io_slice(npts, axis, plane, field, fname);
      };

  };

}
#endif

#endif
