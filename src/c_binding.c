

void scarf_struct_initialize(const int fs[3], const int fe[3], const fpp ds, const int acf, const fpp cl[3], const fpp sigma,  &
                             const int* solver, const fpp* hurst, const fpp* dh, const fpp* poi, const int* npoi, const fpp* mute,  &
                             const fpp* taper, const int* rescale, const int* pad){

   struct_initialize(fs, fe, ds, acf, cl, sigma, solver, hurst, dh, poi, npoi, mute, taper, rescale, pad);

}

void scarf_unstruct_initialize(const int npts[3], const fpp* x, const fpp* y, const fpp* z, const fpp dh, const fpp cl[3], const fpp sigma,  &
                               const int* solver, const fpp* hurst, const fpp* poi, const int* npoi, const fpp* mute, const fpp* taper,      &
                               const int* rescale, const int* pad){

   unstruct_initialize(npts, x, y, z, dh, cl, sigma, solver, hurst, poi, npoi, mute, taper, rescale, pad);

}

void scarf_execute(const int seed, fpp* field, fpp* stats){

   execute(seed, field, stats);

}

void scarf_finalize(){

   finalize();

}
