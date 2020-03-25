enum algorithm {fft, spec};

#ifdef DOUBLE_PREC
typedef double fpp;
//using fpp = double;
#else
typedef float fpp;
//using fpp = float;
#endif

#ifdef __cplusplus
extern "C"{
#endif

// declare FORTRAN subroutines
// extern void struct_initialize(const int fs[], const int fe[], const fpp* ds, const int* acf, const fpp cl[], const fpp* sigma,
//                               const int* solver, const fpp* hurst, const fpp* dh, const fpp* poi, const int* npoi, const fpp* mute,
//                               const fpp* taper, const int* rescale, const int* pad);
//
// extern void unstruct_initialize(const int* npts, const fpp* x, const fpp* y, const fpp* z, const fpp* dh, const int* acf, const fpp cl[], const fpp* sigma,
//                                 const int* solver, const fpp* hurst, const fpp* poi, const int* npoi, const fpp* mute, const fpp* taper, const int* rescale,
//                                 const int* pad);
//
// extern void execute(const int* seed, fpp** field, fpp stats[]);
//
// extern void finalize();
//
// extern void io_one(const int* npts, const fpp** field, const char fname[], const int* nwriters);
//
// extern void io_slice(const int* npts, const int* direction, const int* plane, const fpp** field, const char fname[]);

// declare C functions
void scarf_struct_initialize(const int fs[], const int fe[], const fpp ds, const int acf, const fpp cl[], const fpp sigma,
                             const int* solver, const fpp* hurst, const fpp* dh, const fpp* poi, const int* npoi, const fpp* mute,
                             const fpp* taper, const int* rescale, const int* pad);

void scarf_unstruct_initialize(const int npts, const fpp* x, const fpp* y, const fpp* z, const fpp dh, const int acf, const fpp cl[], const fpp sigma,
                               const int* solver, const fpp* hurst, const fpp* poi, const int* npoi, const fpp* mute, const fpp* taper,
                               const int* rescale, const int* pad);

void scarf_execute(const int seed, fpp* field, fpp stats[]);

void scarf_finalize();

void scarf_io_one(const int* npts, const fpp* field, const char fname[], const int* nwriters);

void scarf_io_slice(const int* npts, const int direction, const int plane, const fpp* field, const char fname[]);

#ifdef __cplusplus
}

namespace Scarf3D
{

  template <algorithm method> class Initialize
  {

    private:

      // variable used to define desired algorithm
      int flag = 0;

    public:

      // constructor structured mesh
      Initialize(const int fs[], const int fe[], const fpp ds, const int acf, const fpp cl[], const fpp sigma,
                 const fpp* hurst = nullptr, const fpp* dh = nullptr, const fpp* poi = nullptr, const int* npoi = nullptr,
                 const fpp* mute = nullptr, const fpp* taper = nullptr, const int* rescale = nullptr,const int* pad = nullptr)
      {

        if (method == spec) flag = 1;

        const int* solver = &flag;

        scarf_struct_initialize(fs, fe, ds, acf, cl, sigma, solver, hurst, dh, poi, npoi, mute, taper, rescale, pad);

      };

      // constructor unstructured mesh
      Initialize(const int npts, const fpp* x, const fpp* y, const fpp* z, const fpp dh, const int acf, const fpp cl[], const fpp sigma,
                 const fpp* hurst = nullptr, const fpp* poi = nullptr, const int* npoi = nullptr, const fpp* mute = nullptr,
                 const fpp* taper = nullptr, const int* rescale = nullptr, const int* pad = nullptr)
      {

        if (method == spec) flag = 1;

        const int* solver = &flag;

        scarf_unstruct_initialize(npts, x, y, z, dh, acf, cl, sigma, solver, hurst, poi, npoi, mute, taper, rescale, pad);

      };

      // destructor
      ~Initialize()
      {
        scarf_finalize();
      };

      // simulate random field
      void execute(const int seed, fpp* field, fpp stats[])
      {
        scarf_execute(seed, field, stats);
      };

      // IO whole model
      void io(const int npts[], const fpp* field, const char fname[], const int* nwriters = nullptr)
      {
        scarf_io_one(npts, field, fname, nwriters);
      };

      // IO model slice
      void io(const int npts[], const int direction, const int plane, const fpp* field, const char fname[])
      {
        scarf_io_slice(npts, direction, plane, field, fname);
      };

  };

}
#endif
