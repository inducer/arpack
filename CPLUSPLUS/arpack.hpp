//
// Boost C++ ARPACK bindings.
//
// Copyright (c) 2004-2006
// Andreas Kloeckner
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee,
// provided that the above copyright notice appear in all copies and
// that both that copyright notice and this permission notice appear
// in supporting documentation.  The authors make no representations
// about the suitability of this software for any purpose.
// It is provided "as is" without express or implied warranty.
//

/* 
 * For ARPACK Copyright, License and Availability see 
 * detail/arpack_proto.hpp 
 */ 




#ifndef HEADER_SEEN_ARPACK_HPP
#define HEADER_SEEN_ARPACK_HPP




#include <vector>
#include <boost/numeric/bindings/traits/type_traits.hpp>  	
#include <boost/numeric/bindings/traits/vector_traits.hpp>  	
#include <boost/scoped_array.hpp>  	
#include <pyublasext/matrix_operator.hpp>
#include <pyublasext/arpack_proto.hpp>
#include <boost/lexical_cast.hpp>
//#include <boost/numeric/bindings/arpack/arpack_proto.hpp>




extern "C"
{
  void dlarnv_(int *idist, int *iseed, int *n, double *x);
}

namespace boost { namespace numeric { namespace bindings {  namespace arpack {

  enum which_eigenvalues { LARGEST_MAGNITUDE,
    SMALLEST_MAGNITUDE,
    LARGEST_REAL_PART,
    SMALLEST_REAL_PART,
    LARGEST_IMAGINARY_PART,
    SMALLEST_IMAGINARY_PART
  };

  namespace detail
  {
    const char *map_which_to_string(which_eigenvalues we)
    {
      switch (we)
      {
        case LARGEST_MAGNITUDE: 
          return "LM";
        case SMALLEST_MAGNITUDE: 
          return "SM";
        case LARGEST_REAL_PART: 
          return "LR";
        case SMALLEST_REAL_PART: 
          return "SR";
        case LARGEST_IMAGINARY_PART: 
          return "LI";
        case SMALLEST_IMAGINARY_PART: 
          return "SI";
        default: 
          throw std::runtime_error("arpack: invalid eigenvalue selector");
      }
    }
  }

  enum arpack_mode { 
    REGULAR_NON_GENERALIZED = 1,
    REGULAR_GENERALIZED = 2,
    SHIFT_AND_INVERT_GENERALIZED = 3
  };

  /* VectorType needs to be complex-valued. */
  template <typename VectorType>
  struct results
  {
      typedef typename boost::numeric::bindings::traits::type_traits<
        typename boost::numeric::bindings::traits::vector_traits<VectorType>::value_type
      >::real_type real_type;
      typedef std::complex<real_type> complex_type;

      typedef std::vector<complex_type> value_container;
      typedef std::vector<VectorType>
      vector_container;

      value_container m_ritz_values;
      vector_container m_ritz_vectors;
  };




  namespace detail
  {
    template <typename VectorType>
    void make_results(unsigned nconv, unsigned n, 
                     typename results<VectorType>::real_type *z, 
                     std::complex<typename results<VectorType>::real_type> *d, 
                     results<VectorType> &my_results)
    {
      // result generation for real types
      // slightly more complicated: take care of complex conjugate pairs

      typedef typename results<VectorType>::real_type real_type;
      typedef std::complex<real_type> complex_type;

      unsigned i = 0;

      my_results.m_ritz_values.clear();
      my_results.m_ritz_vectors.clear();

      while (i < nconv)
      {
        if (imag(d[i]) != 0)
        {
          // complex-conjugate pair
          if (i + 1 >= nconv)
            throw std::runtime_error("arpack: complex pair split up");

          my_results.m_ritz_values.push_back(d[i]);
          my_results.m_ritz_values.push_back(d[i+1]);

          VectorType ritz_vector(n);
          // un-conjugate ritz vector
          for (unsigned j = 0; j < n; j++)
            ritz_vector[j] = complex_type(z[i*n + j], z[(i+1)*n +j]);
          my_results.m_ritz_vectors.push_back(ritz_vector);

          VectorType ritz_vector_2(n);
          // conjugate ritz vector
          for (unsigned j = 0; j < n; j++)
            ritz_vector_2[j] = complex_type(z[i*n + j], -z[(i+1)*n +j]);
          my_results.m_ritz_vectors.push_back(ritz_vector_2);

          i += 2;
        }
        else
        {
          // real eigenvalue, single eigenvector
          my_results.m_ritz_values.push_back(d[i]);
          VectorType ritz_vector(n);
          for (unsigned j = 0; j < n; j++)
            ritz_vector[j] = z[i*n + j];
          my_results.m_ritz_vectors.push_back(ritz_vector);
          i++;
        }
      }
    }

    template <typename VectorType>
    void make_results(unsigned nconv, unsigned n, 
                     std::complex<typename results<VectorType>::real_type> *z, 
                     std::complex<typename results<VectorType>::real_type> *d, 
                     results<VectorType> &my_results)
    {
      // result generation for complex types
      my_results.m_ritz_values.clear();
      my_results.m_ritz_vectors.clear();

      // simple: just copy everything over.
      for (unsigned i = 0; i < nconv; i++)
      {
        my_results.m_ritz_values.push_back(d[i]);

        VectorType ritz_vector(n);
        for (unsigned j = 0; j < n; j++)
          ritz_vector[j] = z[i*n + j];
        my_results.m_ritz_vectors.push_back(ritz_vector);
      }
    }

    template <typename T>
    inline bool is_complex(const T &)
    {
      return false;
    }

    template <typename T2>
    inline bool is_complex(const std::complex<T2> &)
    {
      return true;
    }
  }




  template <typename MatrixOrOperator, 
           typename ResultsVectorType, 
           typename IterationVectorType>
  inline
  void perform_reverse_communication(
      const MatrixOrOperator &op, 
      const MatrixOrOperator *m,
      arpack_mode mode,
      std::complex<typename boost::numeric::bindings::traits::type_traits<
        typename MatrixOrOperator::operand_type::value_type
      >::real_type> spectral_shift,
      int number_of_eigenvalues,
      int number_of_arnoldi_vectors,
      results<ResultsVectorType> &res,
      const IterationVectorType &starting_vector,
      which_eigenvalues which_e = LARGEST_MAGNITUDE,
      typename boost::numeric::bindings::traits::type_traits<
        typename MatrixOrOperator::operand_type::value_type
      >::real_type tolerance = 1e-8,
      int max_iterations = 0
      )
  {
    typedef 
      typename MatrixOrOperator::operand_type::value_type
      value_type;
    typedef 
      typename boost::numeric::bindings::traits::type_traits<value_type>::real_type
      real_type;
    typedef
      std::complex<real_type>
      complex_type;

    int ido = 0;
    char bmat = m == 0 ? 'I' : 'G';
    int n = op.size1();

    if ((unsigned) n != op.size2())
      throw std::runtime_error("arpack: matrix sizes don't match.");
    if ((m != 0) && ((unsigned) n != m->size1() || (unsigned) n != m->size2()))
      throw std::runtime_error("arpack: matrix sizes don't match.");

    char *which = const_cast<char*>(detail::map_which_to_string(which_e));

    boost::scoped_array<value_type> residual(new value_type[n]);
    for (int i = 0; i < n; i++)
      residual[i] = starting_vector[i];

    boost::scoped_array<value_type> v(new value_type[number_of_arnoldi_vectors * n]);
    int ldv = n;

    int iparam[11];
    iparam[1-1] = 1; // exact shifts
    iparam[2-1] = 0; // unused
    iparam[3-1] = max_iterations != 0 ? max_iterations : 10000 * n;
    iparam[4-1] = 1; // block size
    iparam[5-1] = 0; // NCONV
    iparam[6-1] = 0; // IUPD, unused
    iparam[7-1] = mode;
    iparam[8-1] = 0; // NP, something to do with user-provided shifts
    iparam[9-1] = 0; // NUMOP
    iparam[10-1] = 0; // NUMOPB
    iparam[11-1] = 0; // NUMREO

    int info = 1; // we are specifying a previous residual

    int ipntr[14];

    boost::scoped_array<value_type> workd(new value_type[3*n]);
    int lworkl;
    lworkl = 3 * number_of_arnoldi_vectors * number_of_arnoldi_vectors 
      + 6 * number_of_arnoldi_vectors;

    boost::scoped_array<value_type> workl(new value_type[lworkl]);
    boost::scoped_array<double> rwork(new double[number_of_arnoldi_vectors]);

    do
    {
      detail::naupd(
          &ido,
          &bmat,
          &n,
          which,
          &number_of_eigenvalues,
          &tolerance,
          residual.get(),
          &number_of_arnoldi_vectors,
          v.get(),
          &ldv,
          iparam,
          ipntr,
          workd.get(),
          workl.get(),
          &lworkl,
          rwork.get(),
          &info
          );

      switch (info)
      {
        case 0:
          break;
        case 1:
          throw std::runtime_error("arpack, naupd: performed max. number of iterations (1)");
        case 3:
          throw std::runtime_error("arpack, naupd: could not apply shifts (3)");
        case -1:
          throw std::runtime_error("arpack, naupd: n not positive (-1)");
        case -2:
          throw std::runtime_error("arpack, naupd: nev not positive (-2)");
        case -3:
          throw std::runtime_error("arpack, naupd: ncv <= nev or ncv > n (-3)");
        case -4:
          throw std::runtime_error("arpack, naupd: max_iterations must be bigger than zero (-4)");
        case -5:
          throw std::runtime_error("arpack, naupd: invalid WHICH (-5)");
        case -6:
          throw std::runtime_error("arpack, naupd: invalid BMAT (-6)");
        case -7:
          throw std::runtime_error("arpack, naupd: work array too short (-7)");
        case -8:
          throw std::runtime_error("arpack, naupd: LAPACK error (-8)");
        case -9:
          throw std::runtime_error("arpack, naupd: starting vector is zero (-9)");
        case -10:
          throw std::runtime_error("arpack, naupd: invalid MODE (-10)");
        case -11:
          throw std::runtime_error("arpack, naupd: MODE and BMAT don't agree (-11)");
        case -12:
          throw std::runtime_error("arpack, naupd: ISHIFT invalid (-12)");
        case -9999:
          throw std::runtime_error("arpack, naupd: failed to build an Arnoldi factorization (-9999)");
        default:
          throw std::runtime_error("arpack, naupd: invalid naupd error code: "
                                   + boost::lexical_cast<std::string>(info));
      }

      if (ido == -1 || ido == 1 || ido == 2)
      {
        value_type *x = workd.get() + ipntr[1-1] - 1;
        value_type *y = workd.get() + ipntr[2-1] - 1;

        IterationVectorType operand;
        operand.resize(n);
        std::copy(x, x+n, operand.begin());

        IterationVectorType result;
        if (ido == 2) 
        {
          if (m == 0)
            throw std::runtime_error("arpack, rci: multiplication by m requested, but m not supplied");
          else
            result = prod(*m, operand);
        }
        else
          result = prod(op, operand);

        std::copy(result.begin(), result.end(), y);
      }
      else if (ido == 99) /*nothing*/;
      else
        throw std::runtime_error("arpack: reverse communication failure");
    }
    while (ido != 99);

    if (max_iterations != 0)
    {
      if (max_iterations <= iparam[2])
        throw std::runtime_error("arpack: hit iteration count limit");
      max_iterations -= iparam[2];
    }

    {
      // prepare for call to neupd
      int rvec = 1;
      char howmny = 'A';
      boost::scoped_array<int> select(new int[number_of_arnoldi_vectors]); // no-op
      boost::scoped_array<complex_type> d(new complex_type[number_of_eigenvalues+1]);

      unsigned z_size;
      unsigned workev_size;
      if (detail::is_complex(value_type()))
      {
        z_size = number_of_eigenvalues;
        workev_size = 2*number_of_arnoldi_vectors;
      }
      else
      {
        z_size = number_of_eigenvalues+1;
        workev_size = 3*number_of_arnoldi_vectors;
      }

      boost::scoped_array<value_type> z(new value_type[z_size*n]);
      int ldz = n;

      boost::scoped_array<value_type> workev(new value_type[workev_size]);

      detail::neupd(
          &rvec, &howmny, select.get(), d.get(), z.get(), &ldz, 
          &spectral_shift, workev.get(),
          // naupd parameters follow
          &bmat,
          &n,
          which,
          &number_of_eigenvalues,
          &tolerance,
          residual.get(),
          &number_of_arnoldi_vectors,
          v.get(),
          &ldv,
          iparam,
          ipntr,
          workd.get(),
          workl.get(),
          &lworkl,
          rwork.get(),
          &info
          );
      switch (info)
      {
        case 0:
          break;
        case 1:
          throw std::runtime_error("arpack, neupd: schur form could not be reordered (1)");
        case -1:
          throw std::runtime_error("arpack, neupd: n not positive (-1)");
        case -2:
          throw std::runtime_error("arpack, neupd: nev not positive (-2)");
        case -3:
          throw std::runtime_error("arpack, neupd: ncv <= nev or ncv > n (-3)");
        case -5:
          throw std::runtime_error("arpack, neupd: invalid WHICH (-5)");
        case -6:
          throw std::runtime_error("arpack, neupd: invalid BMAT (-6)");
        case -7:
          throw std::runtime_error("arpack, neupd: work array too short (-7)");
        case -8:
          throw std::runtime_error("arpack, neupd: LAPACK error (-8)");
        case -9:
          throw std::runtime_error("arpack, neupd: LAPACK _trevc failed (-9)");
        case -10:
          throw std::runtime_error("arpack, neupd: invalid MODE (-10)");
        case -11:
          throw std::runtime_error("arpack, neupd: MODE and BMAT don't agree (-11)");
        case -12:
          throw std::runtime_error("arpack, neupd: HOWNY = S invalid (-12)");
        case -13:
          throw std::runtime_error("arpack, neupd: HOWNY and RVEC don't agree (-13)");
        case -14:
          throw std::runtime_error("arpack, neupd: no eigenvalues found (-14)");
        default:
          throw std::runtime_error("arpack, neupd: invalid neupd error code: "
                                   + boost::lexical_cast<std::string>(info));
      }

      unsigned nconv = iparam[5-1];
      detail::make_results(nconv, n, z.get(), d.get(), res);
    }
  }
}}}} 




#endif
