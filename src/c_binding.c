#include <stdlib.h>
#include "scarf3d.h"

// C functions are declared (prototyped) in "scarf3d.h"

// declare FORTRAN subroutines
extern void struct_initialize(const int fs[], const int fe[], const real* ds, const int* acf, const real cl[], const real* sigma,
                              const int* solver, const real* hurst, const real* dh, const real* poi, const int* npoi, const real* mute,
                              const real* taper, const int* rescale, const int* pad, const real nc[], const real fc[]);

extern void unstruct_initialize(const int* npts, const real* x, const real* y, const real* z, const real* dh, const int* acf, const real cl[], const real* sigma,
                                const int* solver, const real* hurst, const real* ds, const real* poi, const int* npoi, const real* mute, const real* taper,
                                const int* rescale, const int* pad, const real nc[], const real fc[]);

extern void execute(const int* seed, real** field, real stats[]);

extern void finalize();

extern void io_one(const int* npts, const real** field, const char fname[], const int* nwriters);

extern void io_slice(const int* npts, const char axis[], const int* plane, const real** field, const char fname[]);

// end FORTRAN subroutines

// define C functions

void scarf_opt_init(struct scarf_opt * var){

  var -> solver  = 0;
  var -> hurst   = 0;
  var -> dh      = 0;
  var -> ds      = 0;
  var -> poi     = NULL;
  var -> npoi    = 0;
  var -> mute    = 0;
  var -> taper   = 0;
  var -> rescale = 0;
  var -> pad     = 0;
  var -> nc = NULL;
  var -> fc = NULL;

}

// Note: many arguments are defined as pointers because they are optional on the FORTRAN side.
void scarf_struct_initialize(const int fs[], const int fe[], const real ds, const int acf, const real cl[], const real sigma, struct scarf_opt *var){

   int *solver = NULL, *npoi = NULL, *rescale = NULL, *pad = NULL;
   real *hurst = NULL, *dh = NULL, *poi = NULL, *taper = NULL, *mute = NULL, *nc = NULL, *fc = NULL;

   if (var){
     if (var->solver == 1) solver = &var->solver;
     if (var->hurst > 0) hurst = &var->hurst;
     if (var->dh > 0) dh = &var->dh;
     if (var->poi) poi = var->poi;
     if (var->npoi > 0) npoi = &var->npoi;
     if (var->mute > 0) mute = &var->mute;
     if (var->taper > 0) taper = &var->taper;
     if (var->rescale == 1) rescale = &var->rescale;
     if (var->pad == 1) pad = &var->pad;
     if (var->fc && var->nc){fc = var -> fc; nc = var -> nc;}
   }

   // call FORTRAN subroutine
   struct_initialize(fs, fe, &ds, &acf, cl, &sigma, solver, hurst, dh, poi, npoi, mute, taper, rescale, pad, nc, fc);

}

void scarf_unstruct_initialize(const int npts, const real* x, const real* y, const real* z, const real dh, const int acf, const real cl[], const real sigma, struct scarf_opt *var){

  int *solver = NULL, *npoi = NULL, *rescale = NULL, *pad = NULL;
  real *hurst = NULL, *ds = NULL, *poi = NULL, *taper = NULL, *mute = NULL, *nc = NULL, *fc = NULL;

  if (var){
    if (var->solver == 1) solver = &var->solver;
    if (var->hurst > 0) hurst = &var->hurst;
    if (var->ds > 0) ds = &var->ds;
    if (var->poi) poi = var->poi;
    if (var->npoi > 0) npoi = &var->npoi;
    if (var->mute > 0) mute = &var->mute;
    if (var->taper > 0) taper = &var->taper;
    if (var->rescale == 1) rescale = &var->rescale;
    if (var->pad == 1) pad = &var->pad;
    if (var->fc && var->nc){fc = var -> fc; nc = var -> nc;}
  }

   // call FORTRAN subroutine
   unstruct_initialize(&npts, x, y, z, &dh, &acf, cl, &sigma, solver, hurst, ds, poi, npoi, mute, taper, rescale, pad, nc, fc);

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

void scarf_io_one(const int npts[], const real* field, const char fname[], const int* nwriters){

  io_one(npts, &field, fname, nwriters);

}

void scarf_io_slice(const int npts[], const char axis[], const int plane, const real* field, const char fname[]){

  io_slice(npts, axis, &plane, &field, fname);

}
