MODULE SCARF3D_C_BINDING

  USE, INTRINSIC     :: ISO_C_BINDING
  USE, NON_INTRINSIC :: SCARFLIB
  USE, NON_INTRINSIC :: SCARFLIB_COMMON, ONLY: C_FPP, FPP

  IMPLICIT NONE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS


    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

!scarf_struct_initialize(fs, fe, ds, acf, cl, sigma, 0, hurst, dh, poi, mute, taper, rescale, pad);

    SUBROUTINE SCARF_STRUCT_INITIALIZE(FS, FE, DS, ACF, CL, SIGMA, METHOD, HURST, DH, POI, NPOI, MUTE, TAPER, RESCALE, PAD)  &
      BIND(C, NAME = "scarf_struct_initialize")

      INTEGER(C_INT), DIMENSION(3),                     INTENT(IN) :: FS, FE
      REAL(C_FPP),                                      INTENT(IN) :: DS
      INTEGER(C_INT),                                   INTENT(IN) :: ACF
      REAL(C_FPP),    DIMENSION(3),                     INTENT(IN) :: CL
      REAL(C_FPP),                                      INTENT(IN) :: SIGMA
      INTEGER(C_INT),                                   INTENT(IN) :: METHOD
      REAL(C_FPP),                            OPTIONAL, INTENT(IN) :: DH
      TYPE(C_PTR),                            OPTIONAL, INTENT(IN) :: POI
      INTEGER(C_INT),                         OPTIONAL, INTENT(IN) :: NPOI
      REAL(C_FPP),                            OPTIONAL, INTENT(IN) :: MUTE, TAPER
      INTEGER(C_INT),                         OPTIONAL, INTENT(IN) :: RESCALE, PAD
      REAL(C_FPP),    DIMENSION(:,:), POINTER                      :: F_POI

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (PRESENT(POI)) THEN
        CALL C_F_POINTER(POI, F_POI, [3, NPOI])
      ELSE
        F_POI => NULL()
      ENDIF

      CALL SCARF_INITIALIZE(FS, FE, DS, ACF, CL, SIGMA, METHO, HURST, DH, F_POI, MUTE, TAPER, RESCALE, PAD)

    END SUBROUTINE SCARF_STRUCT_INITIALIZE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

! scarf_unstruct_initialize(npts, x, y, z, dh, acf, cl, sigma, 0, hurst, poi, npoi, mute, taper, rescale, pad);

    SUBROUTINE SCARF_UNSTRUCT_INITIALIZE(NPTS, X, Y, Z, DH, ACF, CL, SIGMA, METHOD, HURST, POI, NPOI, MUTE, TAPER, RESCALE, PAD)  &
      BIND(C, NAME = "scarf_unstruct_initialize")

      INTEGER(C_INT),                         INTENT(IN) :: NPTS
      TYPE(C_PTR),                            INTENT(IN) :: X, Y, Z
      REAL(C_FPP),                            INTENT(IN) :: DH
      INTEGER(C_INT),                         INTENT(IN) :: ACF
      REAL(C_FPP),    DIMENSION(3),           INTENT(IN) :: CL
      REAL(C_FPP),                            INTENT(IN) :: METHOD
      REAL(C_FPP),                  OPTIONAL, INTENT(IN) :: HURST
      TYPE(C_PTR),                  OPTIONAL, INTENT(IN) :: POI
      INTEGER(C_INT),               OPTIONAL, INTENT(IN) :: NPOI
      REAL(C_FPP),                  OPTIONAL, INTENT(IN) :: MUTE, TAPER
      INTEGER(C_INT),               OPTIONAL, INTENT(IN) :: RESCALE, PAD

      REAL(C_FPP),    DIMENSION(:),   POINTER            :: F_X, F_Y, F_Z
      REAL(C_FPP),    DIMENSION(:,:), POINTER            :: F_POI

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL C_F_POINTER(X, F_X, [NPTS])
      CALL C_F_POINTER(Y, F_Y, [NPTS])
      CALL C_F_POINTER(Z, F_Z, [NPTS])

      IF (PRESENT(POI)) THEN
        CALL C_F_POINTER(POI, F_POI, [3, NPOI])
      ELSE
        F_POI => NULL()
      ENDIF

      CALL SCARF_INITIALIZE(F_X, F_Y, F_Z, DH, ACF, CL, SIGMA, METHOD, HURST, F_POI, MUTE, TAPER, RESCALE, PAD)


    END SUBROUTINE SCARF_UNSTRUCT_INITIALIZE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE SCARF_EXECUTE_F2C


    END SUBROUTINE SCARF_EXECUTE_F2C

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE SCARF_FINALIZE_F2C


    END SUBROUTINE SCARF_FINALIZE_F2C

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


END MODULE SCARF3D_C_BINDING
