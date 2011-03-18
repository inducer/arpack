/*
  ARPACK++ v1.0 8/1/1997
  c++ interface to ARPACK code.

  MODULE arpackf.h
  ARPACK FORTRAN routines.

  ARPACK Authors
     Richard Lehoucq
     Danny Sorensen
     Chao Yang
     Dept. of Computational & Applied Mathematics
     Rice University
     Houston, Texas

  Slightly modified for PyLinear by Andreas Kloeckner 2003.
*/

#ifndef ARPACKF_H
#define ARPACKF_H




#include <boost/scoped_array.hpp>  	
#include <complex>




typedef int integer;
typedef int logical;

#define ARPACK_F77NAME(x) x##_




extern "C"
{
// double precision symmetric routines.

  void ARPACK_F77NAME(dsaupd)(integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, double *tol, double *resid,
                       integer *ncv, double *V, integer *ldv,
                       integer *iparam, integer *ipntr, double *workd,
                       double *workl, integer *lworkl, integer *info);

  void ARPACK_F77NAME(dseupd)(logical *rvec, char *HowMny, logical *select,
                       double *d, double *Z, integer *ldz,
                       double *sigma, char *bmat, integer *n,
                       char *which, integer *nev, double *tol,
                       double *resid, integer *ncv, double *V,
                       integer *ldv, integer *iparam, integer *ipntr,
                       double *workd, double *workl,
                       integer *lworkl, integer *info);

// double precision nonsymmetric routines.

  void ARPACK_F77NAME(dnaupd)(integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, double *tol, double *resid,
                       integer *ncv, double *V, integer *ldv,
                       integer *iparam, integer *ipntr, double *workd,
                       double *workl, integer *lworkl, integer *info);

  void ARPACK_F77NAME(dneupd)(logical *rvec, char *HowMny, logical *select,
                       double *dr, double *di, double *Z,
                       integer *ldz, double *sigmar,
                       double *sigmai, double *workev,
                       char *bmat, integer *n, char *which,
                       integer *nev, double *tol, double *resid,
                       integer *ncv, double *V, integer *ldv,
                       integer *iparam, integer *ipntr,
                       double *workd, double *workl,
                       integer *lworkl, integer *info);

// single precision symmetric routines.

  void ARPACK_F77NAME(ssaupd)(integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, float *tol, float *resid,
                       integer *ncv, float *V, integer *ldv,
                       integer *iparam, integer *ipntr, float *workd,
                       float *workl, integer *lworkl, integer *info);

  void ARPACK_F77NAME(sseupd)(logical *rvec, char *HowMny, logical *select,
                       float *d, float *Z, integer *ldz,
                       float *sigma, char *bmat, integer *n,
                       char *which, integer *nev, float *tol,
                       float *resid, integer *ncv, float *V,
                       integer *ldv, integer *iparam, integer *ipntr,
                       float *workd, float *workl,
                       integer *lworkl, integer *info);

// single precision nonsymmetric routines.

  void ARPACK_F77NAME(snaupd)(integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, float *tol, float *resid,
                       integer *ncv, float *V, integer *ldv,
                       integer *iparam, integer *ipntr, float *workd,
                       float *workl, integer *lworkl, integer *info);

  void ARPACK_F77NAME(sneupd)(logical *rvec, char *HowMny, logical *select,
                       float *dr, float *di, float *Z,
                       integer *ldz, float *sigmar,
                       float *sigmai, float *workev, char *bmat,
                       integer *n, char *which, integer *nev,
                       float *tol, float *resid, integer *ncv,
                       float *V, integer *ldv, integer *iparam,
                       integer *ipntr, float *workd, float *workl,
                       integer *lworkl, integer *info);

// single precision complex routines.

  void ARPACK_F77NAME(cnaupd)(integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, float *tol, std::complex<float> *resid,
                       integer *ncv, std::complex<float> *V, integer *ldv,
                       integer *iparam, integer *ipntr, std::complex<float> *workd,
                       std::complex<float> *workl, integer *lworkl,
                       float *rwork, integer *info);

  void ARPACK_F77NAME(cneupd)(logical *rvec, char *HowMny, logical *select,
                       std::complex<float> *d, std::complex<float> *Z, integer *ldz,
                       std::complex<float> *sigma, std::complex<float> *workev,
                       char *bmat, integer *n, char *which, integer *nev,
                       float *tol, std::complex<float> *resid, integer *ncv,
                       std::complex<float> *V, integer *ldv, integer *iparam,
                       integer *ipntr, std::complex<float> *workd,
                       std::complex<float> *workl, integer *lworkl,
                       float *rwork, integer *info);

// double precision complex routines.

  void ARPACK_F77NAME(znaupd)(integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, double *tol, std::complex<double> *resid,
                       integer *ncv, std::complex<double> *V, integer *ldv,
                       integer *iparam, integer *ipntr, std::complex<double> *workd,
                       std::complex<double> *workl, integer *lworkl,
                       double *rwork, integer *info);

  void ARPACK_F77NAME(zneupd)(logical *rvec, char *HowMny, logical *select,
                       std::complex<double> *d, std::complex<double> *Z, integer *ldz,
                       std::complex<double> *sigma, std::complex<double> *workev,
                       char *bmat, integer *n, char *which, integer *nev,
                       double *tol, std::complex<double> *resid, integer *ncv,
                       std::complex<double> *V, integer *ldv, integer *iparam,
                       integer *ipntr, std::complex<double> *workd,
                       std::complex<double> *workl, integer *lworkl,
                       double *rwork, integer *info);

}



namespace boost { namespace numeric { namespace bindings {  namespace arpack { namespace detail {

#define DECLARE_NAUPD(TYPE, BASETYPE) \
  inline void naupd(integer *ido, char *bmat, integer *n, char *which, \
                       /*5*/integer *nev, BASETYPE *tol, TYPE *resid, \
                       /*9*/integer *ncv, TYPE *V, integer *ldv, \
                       /*12*/integer *iparam, integer *ipntr, TYPE *workd, \
                       /*15*/TYPE *workl, integer *lworkl, \
                       BASETYPE *rwork, integer *info) 

#define DEFINE_REAL_NAUPD(BASETYPE, LETTER) \
  DECLARE_NAUPD(BASETYPE, BASETYPE) \
  { \
    ARPACK_F77NAME(LETTER##naupd)(ido, bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, \
        lworkl, info); \
  }

#define DEFINE_COMPLEX_NAUPD(BASETYPE, LETTER) \
  DECLARE_NAUPD(std::complex<BASETYPE>, BASETYPE) \
  { \
    ARPACK_F77NAME(LETTER##naupd)(ido, bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, \
        lworkl, rwork, info); \
  }

#define DECLARE_NEUPD(TYPE, BASETYPE) \
  inline void neupd(logical *rvec, char *HowMny, logical *select, \
                       /*4*/std::complex<BASETYPE> *d, TYPE *Z, integer *ldz, \
                       /*7*/std::complex<BASETYPE> *sigma, TYPE *workev, \
                       /*9*/char *bmat, integer *n, char *which, integer *nev, \
                       /*13*/BASETYPE *tol, TYPE *resid, integer *ncv, \
                       TYPE *V, integer *ldv, integer *iparam, \
                       integer *ipntr, TYPE *workd, \
                       TYPE *workl, integer *lworkl, \
                       BASETYPE *rwork, integer *info) \

#define DEFINE_REAL_NEUPD(BASETYPE, LETTER) \
  DECLARE_NEUPD(BASETYPE, BASETYPE) \
  { \
    int nev_copy = *nev; /* BEWARE: neupd changes nev! */ \
    boost::scoped_array<BASETYPE> d_real(new BASETYPE[nev_copy + 1]); \
    boost::scoped_array<BASETYPE> d_imag(new BASETYPE[nev_copy + 1]); \
    for (int i = 0; i < nev_copy+1; i++) /* This silences valgrind. */ \
    { \
      d_real[i] = 0; \
      d_imag[i] = 0; \
    } \
    BASETYPE sigma_real = std::real(*sigma); \
    BASETYPE sigma_imag = std::imag(*sigma); \
    ARPACK_F77NAME(LETTER##neupd)(rvec, HowMny, select, d_real.get(), d_imag.get(), Z, ldz,  \
        &sigma_real, &sigma_imag, workev, bmat, n, which, nev, \
        tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, \
        info); \
    for (int i = 0; i < nev_copy+1; i++) \
      d[i] = std::complex<BASETYPE>(d_real[i], d_imag[i]); \
  }
#define DEFINE_COMPLEX_NEUPD(BASETYPE, LETTER) \
  DECLARE_NEUPD(std::complex<BASETYPE>, BASETYPE) \
  { \
    ARPACK_F77NAME(LETTER##neupd)(rvec, HowMny, select, d, Z, ldz, sigma, workev, bmat, n, which, nev, \
        tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, rwork,  \
        info); \
  }

  DEFINE_REAL_NAUPD(float, s);
  DEFINE_REAL_NAUPD(double, d);
  DEFINE_COMPLEX_NAUPD(float, c);
  DEFINE_COMPLEX_NAUPD(double, z);

  DEFINE_REAL_NEUPD(float, s);
  DEFINE_REAL_NEUPD(double, d);
  DEFINE_COMPLEX_NEUPD(float, c);
  DEFINE_COMPLEX_NEUPD(double, z);

}}}}}




#endif // ARPACKF_H
