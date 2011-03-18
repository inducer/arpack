//
// Copyright (c) 2004-2008 // Andreas Kloeckner
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee,
// provided that the above copyright notice appear in all copies and
// that both that copyright notice and this permission notice appear
// in supporting documentation.  The authors make no representations
// about the suitability of this software for any purpose.
// It is provided "as is" without express or implied warranty.
//




#ifndef HEADER_SEEN_MATRIX_OPERATOR_HPP
#define HEADER_SEEN_MATRIX_OPERATOR_HPP




#include <boost/type_traits/remove_reference.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/operation.hpp>




namespace pyublasext
{
  template <typename OperandType, typename ResultType = OperandType &>
  class matrix_operator
  {
    public:
      // types 
      typedef typename OperandType::size_type size_type;
      typedef OperandType operand_type;
      typedef typename boost::remove_reference<ResultType>::type result_type;
      typedef ResultType complete_result_type;

      // interface
      virtual ~matrix_operator() { }

      virtual unsigned size1() const = 0;
      virtual unsigned size2() const = 0;

      /** Before using apply, operand and result must have the correct size.
      */
      virtual void apply(const OperandType &operand, ResultType result) const
      {
        if (size2() != operand.size() || size1() != result.size())
          throw std::runtime_error("invalid vector sizes in matrix_operator::apply");
      }

      // matrix_expression compatibility
      const matrix_operator &operator()() const
      {
        return *this;
      }
  };




  template <typename OperandType, typename ResultType = OperandType &>
  class algorithm_matrix_operator : public matrix_operator<OperandType, ResultType>
  {
    protected:
      unsigned m_lastIterationCount;
      unsigned m_debugLevel;

    public:
      algorithm_matrix_operator()
        : m_lastIterationCount(0), m_debugLevel(0)
      { }

      unsigned get_debug_level() const
      {
        return m_debugLevel;
      }

      void set_debug_level(unsigned dl)
      {
        m_debugLevel = dl;
      }

      unsigned get_last_iteration_count() const 
      {
        return m_lastIterationCount;
      }
  };




  template <typename OperandType, typename ResultType = OperandType &>
  class iterative_solver_matrix_operator : public algorithm_matrix_operator<OperandType, ResultType>
  {
    protected:
      unsigned m_maxIterations;
      double m_tolerance;

    public:
      iterative_solver_matrix_operator(unsigned maxit = 0, double tol = 0)
        : m_maxIterations(maxit), m_tolerance(tol)
      { }

      unsigned get_max_iterations() const
      {
        return m_maxIterations;
      }

      void set_max_iterations(unsigned maxit)
      {
        m_maxIterations = maxit;
      }

      double get_tolerance() const
      {
        return m_tolerance;
      }

      void set_tolerance(double tol)
      {
        m_tolerance = tol;
      }
  };




  template <typename MatrixType, 
           typename OperandType, 
           typename ResultType = OperandType &,
           typename MatrixHolder = const MatrixType &>
  class ublas_matrix_operator : public matrix_operator<OperandType, ResultType>
  {
    private:
      MatrixHolder m_matrix;
      typedef matrix_operator<OperandType, ResultType> super;

    public:
      ublas_matrix_operator(const MatrixType &m)
        : m_matrix(m)
      { 
      }

      unsigned size1() const
      {
        return m_matrix.size1();
      }

      unsigned size2() const
      {
        return m_matrix.size2();
      }

      void apply(const OperandType &operand, ResultType result) const
      {
        super::apply(operand, result);

        using namespace boost::numeric::ublas;
        axpy_prod(m_matrix, operand, result, /* init */ true);
      }
  };




  template <typename OperandType, typename ResultType = OperandType &>
  class identity_matrix_operator : public matrix_operator<OperandType, ResultType>
  {
    private:
      typedef matrix_operator<OperandType, ResultType> super;
      unsigned m_size;

    public:
      identity_matrix_operator(unsigned size)
        : m_size(size)
      { }

      unsigned size2() const
      {
        return m_size;
      }
      unsigned size1() const
      {
        return m_size;
      }

      void apply(const OperandType &operand, ResultType result) const
      {
        super::apply(operand, result);
        result.assign(operand);
      }
  };




  template <typename OperandType, typename ResultType = OperandType &,
           typename IntermediateType = OperandType>
  class composite_matrix_operator : public matrix_operator<OperandType, ResultType>
  {
    private:
      typedef matrix_operator<OperandType, ResultType> super;
      const super       &m_outer, &m_inner;

    public:
      composite_matrix_operator(const super &outer, const super &inner)
        : m_outer(outer), m_inner(inner)
      { 
        if (m_inner.size1() != m_outer.size2())
          throw std::runtime_error("composite_matrix_operator: sizes do not match");
      }

      unsigned size1() const
      {
        return m_outer.size1();
      }

      unsigned size2() const
      {
        return m_inner.size2();
      }

      void apply(const OperandType &operand, ResultType result) const
      {
        super::apply(operand, result);

        OperandType temp(m_inner.size1());
        temp.clear();

        m_inner.apply(operand, temp);
        m_outer.apply(temp, result);
      }
  };




  template <typename OperandType, typename ResultType = OperandType &>
  class sum_of_matrix_operators : public matrix_operator<OperandType, ResultType>
  {
    private:
      typedef matrix_operator<OperandType, ResultType> super;

      const super       &m_op1, &m_op2;

    public:
      sum_of_matrix_operators (const super &op1, const super &op2)
        : m_op1(op1), m_op2(op2)
      { 
        if (m_op1.size1() != m_op2.size1())
          throw std::runtime_error("sum_of_matrix_operators: sizes do not match");
        if (m_op1.size2() != m_op2.size2())
          throw std::runtime_error("sum_of_matrix_operators: sizes do not match");
      }

      unsigned size1() const
      {
        return m_op1.size1();
      }

      unsigned size2() const
      {
        return m_op1.size2();
      }

      void apply(const OperandType &operand, ResultType result) const
      {
        super::apply(operand, result);

        ResultType temp(result);
        m_op1.apply(operand, temp);
        m_op2.apply(operand, result);

        result += temp;
      }
  };




  template <typename RealOp, 
           typename ResultType,
           typename OperandType = typename RealOp::operand_type >
  class complex_matrix_operator_adaptor : public matrix_operator<OperandType, ResultType>
  {
    private:
      typedef matrix_operator<OperandType, ResultType> super;

    public:
      typedef
        RealOp
        real_operator;

    private:
      const real_operator       &m_real, &m_imaginary;

    public:
      complex_matrix_operator_adaptor(
          const real_operator &real_part, 
          const real_operator &imaginary_part)
        : m_real(real_part), m_imaginary(imaginary_part)
      { 
        if (m_real.size1() != m_imaginary.size1())
          throw std::runtime_error("complex_matrix_operator_adaptor: sizes do not match");
        if (m_real.size2() != m_imaginary.size2())
          throw std::runtime_error("complex_matrix_operator_adaptor: sizes do not match");
      }

      unsigned size1() const
      { return m_real.size1(); }

      unsigned size2() const
      { return m_real.size2(); }

      void apply(const OperandType &operand, ResultType result) const
      {
        super::apply(operand, result);
        typedef typename real_operator::result_type::value_type real_t;
        typedef std::complex<real_t> complex_t;
        typename real_operator::operand_type 
          operand_real(real(operand)), operand_imag(imag(operand));
        typename real_operator::result_type 
          result_real_1(real(result)), result_imag_1(imag(result)),
          result_real_2(real(result)), result_imag_2(imag(result));

        m_real.apply(operand_real, result_real_1);
        m_imaginary.apply(operand_imag, result_real_2);
        result_real_2 *= -1;

        m_imaginary.apply(operand_real, result_imag_1);
        m_real.apply(operand_imag, result_imag_2);

        OperandType result_imag_12 = result_imag_1 + result_imag_2;

        result.assign(result_real_1 + result_real_2 + complex_t(0,1) * result_imag_12);
      }
  };




  template <typename OperandType, 
           typename ScalarType = typename OperandType::value_type,
           typename ResultType = OperandType &>
  class scalar_multiplication_matrix_operator : public matrix_operator<OperandType, ResultType>
  {
    private:
      typedef
        matrix_operator<OperandType, ResultType>
        super;

      ScalarType m_factor;
      unsigned m_size;

    public:
      scalar_multiplication_matrix_operator(const ScalarType &factor, unsigned size)
        : m_factor(factor), m_size(size)
      { 
      }

      unsigned size1() const
      {
        return m_size;
      }

      unsigned size2() const
      {
        return m_size;
      }

      void apply(const OperandType &operand, ResultType result) const
      {
        super::apply(operand, result);
        result.assign(m_factor * operand);
      }
  };




  // generic prod() interface ---------------------------------------------------
  template <typename OperandType, typename ResultType>
  ResultType prod(
      const matrix_operator<OperandType, ResultType> &mo,
      const OperandType &vec,
      bool clear_result = false)
  {
    ResultType result(mo.size1());
    if (clear_result)
      result.clear();

    mo.apply(vec, result);
    return result;
  }
}




#endif
