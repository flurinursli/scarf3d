extern "C"{

#ifdef DOUBLE_PREC
  typedef double fpp;
  //using fpp = double;
#else
  typedef float fpp;
  //using fpp = float;
#endif

  enum algorithm {fft, spec};

  typedef struct{
    const int acf, rescale, pad;
    const int[3] fs, fe;
    const fpp ds, dh;
    const fpp sigma, hurst;
    const fpp mute, taper;
    const fpp[3] cl;
    fpp* poi;
    fpp* stats;
    const fpp* x;
    const fpp* y;
    const fpp* z;
  } scarf_data;




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

  template <algorithm method> class Initialize{
    private:
      scarf_data scarf_obj;
    public:
      // structured mesh
      Initialize(const int[3] fs, const int[3] fe, const fpp ds, const fpp cl[3], const fpp sigma,  &
                 const fpp hurst = 0, const fpp dh = ds, const fpp* poi, const fpp mute = -1, const fpp taper = -1, const int rescale = 0, const int pad = 0){

                   scarf_obj = c_scarf3d_struct_initialize(fs, fe, ds, cl, sigma, hurst, dh, poi, mute, taper, rescale, pad);

                 };

      // unstructured mesh
      Initialize(const fpp* x, const fpp* y, const fpp* z, const fpp ds, const fpp cl[3], const fpp sigma,  &
                 const fpp hurst = 0, const fpp dh = ds, const fpp* poi, const fpp mute = -1, const fpp taper = -1, const int rescale = 0, const int pad = 0){};

      ~Initialize();

      void execute(scarf_data, const int seed, fpp* field);
   };






      {
        c_obj = scarf_struct_initialize(fs, fe, ds, cl, sigma, hurst, dh, poi, mute, taper, rescale, pad);
      };

  };


}
