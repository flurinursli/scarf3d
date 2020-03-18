extern "C"{

#ifdef DOUBLE_PREC
  typedef double fpp;
  //using fpp = double;
#else
  typedef float fpp;
  //using fpp = float;
#endif

  enum algorithm {fft, spec};

  void scarf_struct_initialize(const int fs[3], const int fe[3], const fpp ds, const int acf, const fpp cl[3], const fpp sigma,  &
                               const int* solver, const fpp* hurst, const fpp* dh, const fpp* poi, const int* npoi, const fpp* mute,  &
                               const fpp* taper, const int* rescale, const int* pad);

  void scarf_unstruct_initialize(const int npts[3], const fpp* x, const fpp* y, const fpp* z, const fpp dh, const fpp cl[3], const fpp sigma,  &
                                 const int* solver, const fpp* hurst, const fpp* poi, const int* npoi, const fpp* mute, const fpp* taper,      &
                                 const int* rescale, const int* pad);

  void scarf_execute(const int seed, fpp* field, fpp* stats);

  void scarf_finalize();

}


// call as: Scarf3D::Initialize<fft> S(fs, fe, etc.);
// S.execute(seed, field, stats);

namespace Scarf3D{

  template <algorithm method> class Initialize{
    private:
    public:
      // structured mesh
      Initialize(const int fs[3], const int fe[3], const fpp ds, const int acf, const fpp cl[3], const fpp sigma,  &
                 const int* solver = nullptr, const fpp* hurst = nullptr, const fpp* dh = nullptr, const fpp* poi = nullptr, const int* npoi = nullptr,  &
                 const fpp* mute = nullptr, const fpp* taper = nullptr, const int* rescale = nullptr, const int* pad = nullptr){

         if method == spec{
           const int* solver = 1;
         }

         void scarf_struct_initialize(fs, fe, ds, acf, cl, sigma, solver, hurst, dh, poi, npoi, mute, taper, rescale, pad);

      };

      // unstructured mesh
      Initialize(const int npts[3], const fpp* x, const fpp* y, const fpp* z, const fpp dh, const fpp cl[3], const fpp sigma,  &
                 const int* solver = nullptr, const fpp* hurst = nullptr, const fpp* poi = nullptr, const int* npoi = nullptr,  &
                 const fpp* mute = nullptr, const fpp* taper = nullptr, const int* rescale = nullptr, const int* pad = nullptr){

         if method == spec{
            const int* solver = 1;
         }

         void scarf_unstruct_initialize(npts, x, y, z, dh, acf, cl, sigma, solver, hurst, poi, npoi, mute, taper, rescale, pad);

      };

      // destructor
      ~Initialize(){
        void scarf_finalize();
      };

      void execute(const int seed, fpp* field, fpp* stats){
        void scarf_execute(seed, field, stats);
      };

   };

}
