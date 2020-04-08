#ifndef SCARF3D_H_

#define SCARF3D_H_

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

enum algorithm {fft, spec};

struct scarf_opt{
  int solver;
  fpp hurst;
  fpp dh;
  fpp ds;
  fpp* poi;
  int npoi;
  fpp mute;
  fpp taper;
  int rescale;
  int pad;
  fpp* nc;
  fpp* fc;
};

// declare C functions
void scarf_opt_init(struct scarf_opt * var);

void scarf_struct_initialize(const int fs[], const int fe[], const fpp ds, const int acf, const fpp cl[], const fpp sigma, struct scarf_opt *var);
                             //const int* solver, const fpp* hurst, const fpp* dh, const fpp* poi, const int* npoi, const fpp* mute,
                             //const fpp* taper, const int* rescale, const int* pad);

void scarf_unstruct_initialize(const int npts, const fpp* x, const fpp* y, const fpp* z, const fpp dh, const int acf, const fpp cl[], const fpp sigma, struct scarf_opt *var);
                               //const int* solver, const fpp* hurst, const fpp* ds, const fpp* poi, const int* npoi, const fpp* mute, const fpp* taper,
                               //const int* rescale, const int* pad);

void scarf_execute(const int seed, fpp* field, fpp stats[]);

void scarf_finalize();

void scarf_io_one(const int* npts, const fpp* field, const char fname[], const int* nwriters);

void scarf_io_slice(const int* npts, const int direction, const int plane, const fpp* field, const char fname[]);

#ifdef __cplusplus
}

namespace Scarf3D
{

  scarf_opt options {0, 0, 0, 0, nullptr, 0, 0, 0, 0, 0, nullptr, nullptr};

  template <algorithm method> class Initialize
  {

    private:

    public:

      // constructor structured mesh
      Initialize(const int fs[], const int fe[], const fpp ds, const int acf, const fpp cl[], const fpp sigma)
      {

        if (method == spec) options.solver = 1;

        scarf_struct_initialize(fs, fe, ds, acf, cl, sigma, &options);

      };

      // constructor unstructured mesh
      Initialize(const int npts, const fpp* x, const fpp* y, const fpp* z, const fpp dh, const int acf, const fpp cl[], const fpp sigma)
      {

        if (method == spec) options.solver = 1;

        scarf_unstruct_initialize(npts, x, y, z, dh, acf, cl, sigma, &options);

      };

      // destructor
      ~Initialize()
      {
        scarf_finalize();
        options = {0, 0, 0, 0, nullptr, 0, 0, 0, 0, 0, nullptr, nullptr};
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

#endif
