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
   To provide a C/C++ interface to the SCARF3D library

 Revisions:
     Date                    Description of change
     ====                    =====================
   04/05/20                  original version
*/

#include <stdlib.h>
#include "scarf3d.h"

// C functions are declared (prototyped) in "scarf3d.h"

// declare FORTRAN subroutines
extern void struct_initialize(const int* nd, const int fs[], const int fe[], const real* ds, const int* acf, const real cl[], const real* sigma,
                              const int* method, const real* hurst, const real* dh, real** poi, const int* npoi, const real* mute,
                              const real* taper, const int* rescale, const int* pad, const real nc[], const real fc[], const real* alpha,
                              const real* beta, const real* gamma);

extern void unstruct_initialize(const int* nd, const int* npts, const real* dh, const int* acf, const real cl[], const real* sigma, const real* x, const real* y, const real* z,
                                const int* method, const real* hurst, const real* ds, real** poi, const int* npoi, const real* mute, const real* taper,
                                const int* rescale, const int* pad, const real nc[], const real fc[], const real* alpha, const real* beta, const real* gamma);

extern void execute(const int* seed, real** field, real stats[]);

extern void finalize();

extern void io_one(const int* nd, const int* npts, const real** field, const char fname[], const int* nwriters);

extern void io_slice(const int* npts, const char axis[], const int* plane, const real** field, const char fname[]);

// end FORTRAN subroutines

// define C functions

void scarf_opt_init(struct scarf_opt * var){

  var -> method  = 0;
  var -> hurst   = 0;
  var -> dh      = 0;
  var -> ds      = 0;
  var -> poi     = NULL;
  var -> npoi    = 0;
  var -> mute    = 0;
  var -> taper   = 0;
  var -> rescale = 0;
  var -> pad     = 0;
  var -> nc      = NULL;
  var -> fc      = NULL;
  var -> alpha   = 0;
  var -> beta    = 0;
  var -> gamma   = 0;

}

// Note: many arguments are defined as pointers because they are optional on the FORTRAN side.
void scarf_cart_initialize(const int nd, const int fs[], const int fe[], const real ds, const int acf, const real cl[], const real sigma, struct scarf_opt *var){

   int *method = NULL, *npoi = NULL, *rescale = NULL, *pad = NULL;
   real *hurst = NULL, *dh = NULL, **poi = NULL, *taper = NULL, *mute = NULL, *nc = NULL, *fc = NULL;
   real *alpha = NULL, *beta = NULL, *gamma = NULL;

   if (var){
     if (var->method == 1) method = &var->method;
     if (var->hurst > 0) hurst = &var->hurst;
     if (var->dh > 0) dh = &var->dh;
     if (var->poi) poi = &var->poi;
     if (var->npoi > 0) npoi = &var->npoi;
     if (var->mute > 0) mute = &var->mute;
     if (var->taper > 0) taper = &var->taper;
     if (var->rescale == 1) rescale = &var->rescale;
     if (var->pad == 1) pad = &var->pad;
     if (var->fc && var->nc){fc = var -> fc; nc = var -> nc;}
     if (var->alpha != 0) alpha = &var->alpha;
     if (var->beta != 0) beta = &var->beta;
     if (var->gamma != 0) gamma = &var->gamma;
   }

   // call FORTRAN subroutine
   struct_initialize(&nd, fs, fe, &ds, &acf, cl, &sigma, method, hurst, dh, poi, npoi, mute, taper, rescale, pad, nc, fc, alpha, beta, gamma);

}

void scarf_nocart_initialize(const int nd, const int npts, const real* x, const real* y, const real* z, const real dh, const int acf, const real cl[], const real sigma, struct scarf_opt *var){

  int *method = NULL, *npoi = NULL, *rescale = NULL, *pad = NULL;
  real *hurst = NULL, *ds = NULL, **poi = NULL, *taper = NULL, *mute = NULL, *nc = NULL, *fc = NULL;
  real *alpha = NULL, *beta = NULL, *gamma = NULL;

  if (var){
    if (var->method == 1) method = &var->method;
    if (var->hurst > 0) hurst = &var->hurst;
    if (var->ds > 0) ds = &var->ds;
    if (var->poi) poi = &var->poi;
    if (var->npoi > 0) npoi = &var->npoi;
    if (var->mute > 0) mute = &var->mute;
    if (var->taper > 0) taper = &var->taper;
    if (var->rescale == 1) rescale = &var->rescale;
    if (var->pad == 1) pad = &var->pad;
    if (var->fc && var->nc){fc = var -> fc; nc = var -> nc;}
    if (var->alpha != 0) alpha = &var->alpha;
    if (var->beta != 0) beta = &var->beta;
    if (var->gamma != 0) gamma = &var->gamma;
  }

   // call FORTRAN subroutine
   //unstruct_initialize(&nd, &npts, x, y, z, &dh, &acf, cl, &sigma, method, hurst, ds, poi, npoi, mute, taper, rescale, pad, nc, fc, alpha, beta, gamma);
   unstruct_initialize(&nd, &npts, &dh, &acf, cl, &sigma, x, y, z, method, hurst, ds, poi, npoi, mute, taper, rescale, pad, nc, fc, alpha, beta, gamma);

}

void scarf_execute(const int seed, real* field, real stats[]){

   // call FORTRAN subroutine
   execute(&seed, &field, stats);

}

void scarf_finalize(){

   // call FORTRAN subroutine
   finalize();

}

// define C I/O functions

void scarf_io_one(const int nd, const int npts[], const real* field, const char fname[], const int* nwriters){

  io_one(&nd, npts, &field, fname, nwriters);

}

void scarf_io_slice(const int npts[], const char axis[], const int plane, const real* field, const char fname[]){

  io_slice(npts, axis, &plane, &field, fname);

}
