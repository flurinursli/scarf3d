#include <stdlib.h>
#include "scarf3d.h"


// #ifdef DOUBLE_PREC
// typedef double fpp;
// //using fpp = double;
// #else
// typedef float fpp;
// //using fpp = float;
// #endif


// C functions are declared (prototyped) in "scarf3d.h"

// declare FORTRAN subroutines
extern void struct_initialize(const int fs[], const int fe[], const fpp* ds, const int* acf, const fpp cl[], const fpp* sigma,
                              const int* solver, const fpp* hurst, const fpp* dh, const fpp* poi, const int* npoi, const fpp* mute,
                              const fpp* taper, const int* rescale, const int* pad, const fpp nc[], const fpp fc[]);

extern void unstruct_initialize(const int* npts, const fpp* x, const fpp* y, const fpp* z, const fpp* dh, const int* acf, const fpp cl[], const fpp* sigma,
                                const int* solver, const fpp* hurst, const fpp* ds, const fpp* poi, const int* npoi, const fpp* mute, const fpp* taper,
                                const int* rescale, const int* pad, const fpp nc[], const fpp fc[]);

extern void execute(const int* seed, fpp** field, fpp stats[]);

extern void finalize();

extern void io_one(const int* npts, const fpp** field, const char fname[], const int* nwriters);

extern void io_slice(const int* npts, const char axis[], const int* plane, const fpp** field, const char fname[]);

// end FORTRAN subroutines

// define C functions

void scarf_opt_init(struct scarf_opt * var){

printf("hoi2\n");

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
  printf("hoi3\n");
  var -> nc[0]   = 0;
  var -> nc[1]   = 0;
  var -> nc[2]   = 0;
  var -> fc[0]   = 0;
  var -> fc[1]   = 0;
  var -> fc[2]   = 0;

printf("hoi4\n");

}



// Note: many arguments are defined as pointers because they are OPTIONAL on the FORTRAN side.
void scarf_struct_initialize(const int fs[], const int fe[], const fpp ds, const int acf, const fpp cl[], const fpp sigma, struct scarf_opt *var){

   int *solver = NULL, *npoi = NULL, *rescale = NULL, *pad = NULL;
   fpp *hurst = NULL, *dh = NULL, *poi = NULL, *taper = NULL, *mute = NULL, *nc = NULL, *fc = NULL;

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
     if ( (var->fc[0] - var->nc[0]) > 0 && (var->fc[1] - var->nc[1]) > 0 && (var->fc[2] - var->nc[2]) > 0){
       fc = var -> fc;
       nc = var -> nc;
     }
   }

   // call FORTRAN subroutine
   struct_initialize(fs, fe, &ds, &acf, cl, &sigma, solver, hurst, dh, poi, npoi, mute, taper, rescale, pad, nc, fc);

}

void scarf_unstruct_initialize(const int npts, const fpp* x, const fpp* y, const fpp* z, const fpp dh, const int acf, const fpp cl[], const fpp sigma, struct scarf_opt *var){

  int *solver = NULL, *npoi = NULL, *rescale = NULL, *pad = NULL;
  fpp *hurst = NULL, *ds = NULL, *poi = NULL, *taper = NULL, *mute = NULL, *nc = NULL, *fc = NULL;

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
    if ( (var->fc[0] - var->nc[0]) > 0 && (var->fc[1] - var->nc[1]) > 0 && (var->fc[2] - var->nc[2]) > 0){
      fc = var->fc;
      nc = var->nc;
    }
  }

   // call FORTRAN subroutine
   unstruct_initialize(&npts, x, y, z, &dh, &acf, cl, &sigma, solver, hurst, ds, poi, npoi, mute, taper, rescale, pad, nc, fc);

}

void scarf_execute(const int seed, fpp* field, fpp stats[]){

   // call FORTRAN subroutine
   execute(&seed, &field, stats);

}

void scarf_finalize(){

   // call FORTRAN subroutine
   finalize();

}

// define C I/O functions

void scarf_io_one(const int npts[], const fpp* field, const char fname[], const int* nwriters){

  io_one(npts, &field, fname, nwriters);

}

void scarf_io_slice(const int npts[], const char axis[], const int plane, const fpp* field, const char fname[]){

  io_slice(npts, axis, &plane, &field, fname);

}
