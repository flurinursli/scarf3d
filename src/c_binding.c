
// C functions are declared (prototyped) in "scarf3d.h"

// declare FORTRAN subroutines
extern void struct_initialize(const int fs[], const int fe[], const fpp* ds, const int* acf, const fpp cl[], const fpp* sigma,          &
                              const int* solver, const fpp* hurst, const fpp* dh, const fpp* poi, const int* npoi, const fpp* mute,     &
                              const fpp* taper, const int* rescale, const int* pad);

extern void unstruct_initialize(const int* npts, const fpp* x, const fpp* y, const fpp* z, const fpp* dh, const int* acf, const fpp cl[], const fpp* sigma,  &
                                const fpp* solver, const fpp* hurst, const fpp* poi, const int* npoi, const fpp* mute, const fpp* taper, const fpp* rescale, &
                                const fpp* pad);

extern void execute(const int* seed, fpp** field, fpp stats[]);

extern void finalize();
// end FORTRAN subroutines

// define C functions
// Note: many arguments are defined as pointers because that way we can simulate OPTIONAL arguments for FORTRAN
void scarf_struct_initialize(const int fs[], const int fe[], const fpp ds, const int acf, const fpp cl[], const fpp sigma,          &
                             const int* solver, const fpp* hurst, const fpp* dh, const fpp* poi, const int* npoi, const fpp* mute,  &
                             const fpp* taper, const int* rescale, const int* pad){

   // call FORTRAN subroutine
   struct_initialize(fs, fe, &ds, &acf, cl, &sigma, solver, hurst, dh, poi, npoi, mute, taper, rescale, pad);

}

void scarf_unstruct_initialize(const int npts[], const fpp* x, const fpp* y, const fpp* z, const fpp dh, const fpp cl[], const fpp sigma,  &
                               const int* solver, const fpp* hurst, const fpp* poi, const int* npoi, const fpp* mute, const fpp* taper,    &
                               const int* rescale, const int* pad){

   // call FORTRAN subroutine
   unstruct_initialize(&npts, x, y, z, &dh, cl, &sigma, solver, hurst, poi, npoi, mute, taper, rescale, pad);

}

void scarf_execute(const int seed, fpp* field, fpp stats[]){

   // call FORTRAN subroutine
   execute(&seed, &field, stats);

}

void scarf_finalize(){

   // call FORTRAN subroutine
   finalize();

}
