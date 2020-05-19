/*
Purpose:
  To provide an header file for the SCARF3D library in C/C++

Revisions:
    Date                    Description of change
    ====                    =====================
  04/05/20                  original version

*/

#ifndef SCARF3D_H_

#define SCARF3D_H_

#ifdef DOUBLE_PREC
typedef double real;
//using real = double;
#else
typedef float real;
//using real = float;
#endif

#ifdef __cplusplus
extern "C"{
#endif

enum algorithm {fft, spec};

struct scarf_opt{
  int solver;
  real hurst;
  real dh;
  real ds;
  real* poi;
  int npoi;
  real mute;
  real taper;
  int rescale;
  int pad;
  real * nc;
  real * fc;
};

// declare C functions
void scarf_opt_init(struct scarf_opt * var);

void scarf_struct_initialize(const int fs[], const int fe[], const real ds, const int acf, const real cl[], const real sigma, struct scarf_opt *var);

void scarf_unstruct_initialize(const int npts, const real* x, const real* y, const real* z, const real dh, const int acf, const real cl[], const real sigma, struct scarf_opt *var);

void scarf_execute(const int seed, real* field, real stats[]);

void scarf_finalize();

void scarf_io_one(const int* npts, const real* field, const char fname[], const int* nwriters);

void scarf_io_slice(const int* npts, const char axis[], const int plane, const real* field, const char fname[]);

// define C++ class
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
      Initialize(const int fs[], const int fe[], const real ds, const int acf, const real cl[], const real sigma)
      {

        if (method == spec) options.solver = 1;

        scarf_struct_initialize(fs, fe, ds, acf, cl, sigma, &options);

      };

      // constructor unstructured mesh
      Initialize(const int npts, const real* x, const real* y, const real* z, const real dh, const int acf, const real cl[], const real sigma)
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
      void execute(const int seed, real* field, real stats[])
      {
        scarf_execute(seed, field, stats);
      };

      // IO whole model
      void io(const int npts[], const real* field, const char fname[], const int* nwriters = nullptr)
      {
        scarf_io_one(npts, field, fname, nwriters);
      };

      // IO model slice
      void io(const int npts[], const char axis[], const int plane, const real* field, const char fname[])
      {
        scarf_io_slice(npts, axis, plane, field, fname);
      };

  };

}
#endif

#endif
