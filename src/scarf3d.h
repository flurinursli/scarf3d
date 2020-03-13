extern "C"{

#ifdef DOUBLE_PREC
  typedef double fpp;
#else
  typedef float fpp;
#endif

  typedef struct{
    const int acf, rescale, pad;
    const int[3] fs, fe;
    const fpp ds, dh;
    const fpp sigma, hurst;
    const fpp mute, taper;
    const fpp[3] cl;
    fpp* poi;
    fpp* stats;
  } scarf_struct;

  typedef struct{
    const int acf, rescale, pad;
    const fpp ds, dh;
    const fpp sigma, hurst;
    const fpp mute, taper;
    const fpp[3] cl;
    fpp* poi;
    fpp* stats;
    const fpp* x;
    const fpp* y;
    const fpp* z;
  } scarf_unstruct;


  scarf_struct scarf_struct_initialize(const int[3] fs, const int[3] fe, const fpp ds, const fpp[3] cl, const fpp sigma,
                                       const fpp hurst, const fpp dh, const fpp* poi, const fpp mute, const fpp taper,
                                       const int rescale, const int pad);

  void scarf_struct_execute(scarf_struct S, const int seed, fpp* field);

  scarf_unstruct scarf_unstruct_initialize(const fpp* x, const fpp* y, const fpp* z, const fpp dh, const fpp[3] cl, const fpp sigma,
                                           const fpp hurst, const fpp dh, const fpp* poi, const fpp mute, const fpp taper,
                                           const int rescale, const int pad);

  void scarf_unstruct_execute(scarf_unstruct S, const int seed, fpp* field);

}

namespace Scarf3D{

  template < int dims> class Initialize<dims, "fft">{
    private:
      scarf_struct c_obj;
    public:
      Initialize(const int[3] fs, const int[3] fe, const fpp ds, const ){
        c_obj = scarf_struct_initialize(fs, fe, ds, cl, sigma, hurst, dh, poi, mute, taper, rescale, pad);
      };

  };


}
