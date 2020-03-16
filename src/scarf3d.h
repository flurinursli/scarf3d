extern "C"{

#ifdef DOUBLE_PREC
  typedef double fpp;
  //using fpp = double;
#else
  typedef float fpp;
  //using fpp = float;
#endif

  enum algorithm {fft, spec};

  void scarf_struct_initialize(const int fs[3], const int fe[3], const fpp ds, const int acf, const fpp cl[3], const fpp sigma, const int,   &
                               const fpp hurst, const fpp dh, const fpp* poi, const fpp mute, const fpp taper, const int rescale, const int pad);

  void scarf_unstruct_initialize(const fpp* x, const fpp* y, const fpp* z, const fpp dh, const fpp cl[3], const fpp sigma, const int,  &
                                 const fpp hurst, const fpp* poi, const fpp mute, const fpp taper, const int rescale, const int pad);

  void scarf_execute(const int seed, fpp* field, fpp stats[8]);

  void finalize();

}


namespace Scarf3D{

  template <algorithm method> class Initialize{
    private:
    public:
      // structured mesh
      Initialize(const int[3] fs, const int[3] fe, const fpp ds, const int acf, const fpp cl[3], const fpp sigma,  &
                 const fpp hurst = 0, const fpp dh = ds, const fpp* poi, const fpp mute = -1, const fpp taper = -1, const int rescale = 0, const int pad = 0){

         if method == fft{
            void scarf_struct_initialize(fs, fe, ds, acf, cl, sigma, 0, hurst, dh, poi, mute, taper, rescale, pad);
         }
         else{
            void scarf_struct_initialize(fs, fe, ds, acf, cl, sigma, 1, hurst, dh, poi, mute, taper, rescale, pad);
         }

      };

      // unstructured mesh
      Initialize(const fpp* x, const fpp* y, const fpp* z, const fpp dh, const fpp cl[3], const fpp sigma,  &
                 const fpp hurst = 0, const fpp* poi, const fpp mute = -1, const fpp taper = -1, const int rescale = 0, const int pad = 0){

                   if method == fft{
                     void scarf_unstruct_initialize(x, y, z, dh, acf, cl, sigma, 0, hurst, poi, mute, taper, rescale, pad);
                   }
                   else{
                     void scarf_unstruct_initialize(x, y, z, dh, acf, cl, sigma, 1, hurst, poi, mute, taper, rescale, pad);
                   }

                 };

      // destructor
      ~Initialize(){
        void scarf_finalize();
      };

      void execute(const int seed, fpp* field, fpp stats[8]){
        void scarf_execute(seed, field, stats);
      };

   };

}
