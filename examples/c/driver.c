








// INITIALIZE_STRUCTURED(FS, FE, DS, ACF, CL, SIGMA, METHOD, HURST, DH, POI, MUTE, TAPER, RESCALE, PAD)
// INITIALIZE_UNSTRUCTURED(X, Y, Z, DH, ACF, CL, SIGMA, METHOD, HURST, POI, MUTE, TAPER, RESCALE, PAD)

//
scarf_obj S = scarf_struct_initialize(fs, fe, ds, acf, cl, sigma, 0, hurst, NULL, NULL, NULL, NULL, 1);
scarf_obj S = scarf_unstruct_initialize(x1, y1, z1, dh, acf, cl, sigma, 0, hurst, NULL, NULL, NULL, 1);

scarf_execute(S, seed, field);

scarf_finalize(&S);
