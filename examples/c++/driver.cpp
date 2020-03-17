
int main(){



  // INITIALIZE_STRUCTURED(FS, FE, DS, ACF, CL, SIGMA, METHOD, HURST, DH, POI, MUTE, TAPER, RESCALE, PAD)
  // INITIALIZE_UNSTRUCTURED(X, Y, Z, DH, ACF, CL, SIGMA, METHOD, HURST, POI, MUTE, TAPER, RESCALE, PAD)

  //

  Scarf3D::Initialize<spec> S(fs, fe, ds, acf, cl, sigma, hurst = hurst, pad = 1);
  Scarf3D::Initialize<fft> S(x, y, z, dh, acf, cl, sigma, hurst = hurst);

  S.execute(seed, field, stats);

  // S is destroyed automatically






}
