      COMPLEX FUNCTION AR_CLADIV( X, Y )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      COMPLEX            X, Y
*     ..
*
*  Purpose
*  =======
*
*  AR_CLADIV := X / Y, where X and Y are complex.  The computation of X / Y
*  will not overflow on an intermediary step unless the results
*  overflows.
*
*  Arguments
*  =========
*
*  X       (input) COMPLEX
*  Y       (input) COMPLEX
*          The complex scalars X and Y.
*
*  =====================================================================
*
*     .. Local Scalars ..
      REAL               ZI, ZR
*     ..
*     .. External Subroutines ..
      EXTERNAL           AR_SLADIV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          AIMAG, CMPLX, REAL
*     ..
*     .. Executable Statements ..
*
      CALL AR_SLADIV( REAL( X ), AIMAG( X ), REAL( Y ), AIMAG( Y ), ZR,
     $             ZI )
      AR_CLADIV = CMPLX( ZR, ZI )
*
      RETURN
*
*     End of AR_CLADIV
*
      END
