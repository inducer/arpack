      LOGICAL          FUNCTION AR_LSAMEN( N, CA, CB )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    CA, CB
      INTEGER            N
*     ..
*
*  Purpose
*  =======
*
*  AR_LSAMEN  tests if the first N letters of CA are the same as the
*  first N letters of CB, regardless of case.
*  AR_LSAMEN returns .TRUE. if CA and CB are equivalent except for case
*  and .FALSE. otherwise.  AR_LSAMEN also returns .FALSE. if LEN( CA )
*  or LEN( CB ) is less than N.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of characters in CA and CB to be compared.
*
*  CA      (input) CHARACTER*(*)
*  CB      (input) CHARACTER*(*)
*          CA and CB specify two character strings of length at least N.
*          Only the first N characters of each string will be accessed.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
*     ..
*     .. External Functions ..
      LOGICAL            AR_LSAME
      EXTERNAL           AR_LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LEN
*     ..
*     .. Executable Statements ..
*
      AR_LSAMEN = .FALSE.
      IF( LEN( CA ).LT.N .OR. LEN( CB ).LT.N )
     $   GO TO 20
*
*     Do for each character in the two strings.
*
      DO 10 I = 1, N
*
*        Test if the characters are equal using AR_LSAME.
*
         IF( .NOT.AR_LSAME( CA( I: I ), CB( I: I ) ) )
     $      GO TO 20
*
   10 CONTINUE
      AR_LSAMEN = .TRUE.
*
   20 CONTINUE
      RETURN
*
*     End of AR_LSAMEN
*
      END
