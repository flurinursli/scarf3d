MODULE SCARFLIB

  USE, INTRINSIC     :: ISO_FORTRAN_ENV
  USE, NON_INTRINSIC :: MPI

  IMPLICIT NONE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  PRIVATE

  PUBLIC :: SCARF3D_FFT !, SCARF3D_SPEC

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  ! SET PRECISION FOR FLOATING POINT NUMBERS

!#ifdef DOUBLE_PREC
!  INTEGER, PARAMETER, PUBLIC :: FPP       = REAL64
!  INTEGER, PARAMETER         :: REAL_TYPE = MPI_DOUBLE_PRECISION
!#else
  INTEGER, PARAMETER, PUBLIC :: FPP       = REAL32
  INTEGER, PARAMETER         :: REAL_TYPE = MPI_REAL
!#endif

  ! INTEGERS HAVE ALWAYS SAME "PRECISION"
  INTEGER, PARAMETER, PUBLIC :: IPP = INT32

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  ! INTERFACES
  INTERFACE
    MODULE SUBROUTINE SCARF3D_FFT(COMM, N, DH, ACF, CL, SIGMA, HURST, SEED, POI, MUTE, TAPER)
      INTEGER(IPP),                      INTENT(IN) :: COMM               !< MPI COMMUNICATOR
      INTEGER(IPP),     DIMENSION(3),    INTENT(IN) :: N                  !< MODEL SIZE: NUMBER OF POINTS ALONG X, Y, Z
      REAL(FPP),                         INTENT(IN) :: DH                 !< GRID-STEP
      CHARACTER(LEN=*),                  INTENT(IN) :: ACF                !< AUTOCORRELATION FUNCTION: "VK" OR "GAUSS"
      REAL(FPP),        DIMENSION(3),    INTENT(IN) :: CL                 !< CORRELATION LENGTH
      REAL(FPP),                         INTENT(IN) :: SIGMA              !< STANDARD DEVIATION
      REAL(FPP),                         INTENT(IN) :: HURST              !< HURST EXPONENT
      INTEGER(IPP),                      INTENT(IN) :: SEED               !< SEED NUMBER
      INTEGER(IPP),     DIMENSION(:,:),  INTENT(IN) :: POI                !< LOCATION OF POINT(S)-OF-INTEREST
      INTEGER(IPP),                      INTENT(IN) :: MUTE               !< NUMBER OF POINTS WHERE MUTING IS APPLIED
      INTEGER(IPP),                      INTENT(IN) :: TAPER              !< NUMBER OF POINTS WHERE TAPERING IS APPLIED
    END SUBROUTINE SCARF3D_FFT
  END INTERFACE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE SCARFLIB
