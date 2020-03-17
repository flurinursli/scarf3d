MODULE SCARF3D_C_BINDING

  USE, INTRINSIC     :: ISO_C_BINDING
  USE, NON_INTRINSIC :: SCARFLIB
  USE, NON_INTRINSIC :: SCARFLIB_COMMON, ONLY: C_FPP, FPP, IPP

  IMPLICIT NONE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTEGER(IPP),               PRIVATE :: STRUCTURED
  INTEGER(IPP),               PRIVATE :: N
  INTEGER(IPP), DIMENSION(3), PRIVATE :: DIMS

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE SCARF_STRUCT_INITIALIZE(FS, FE, DS, ACF, CL, SIGMA, METHOD, HURST, DH, POI, NPOI, MUTE, TAPER, RESCALE, PAD)  &
      BIND(C, NAME = "scarf_struct_initialize")

      INTEGER(C_INT), DIMENSION(3),           INTENT(IN) :: FS, FE
      REAL(C_FPP),                            INTENT(IN) :: DS
      INTEGER(C_INT),                         INTENT(IN) :: ACF
      REAL(C_FPP),    DIMENSION(3),           INTENT(IN) :: CL
      REAL(C_FPP),                            INTENT(IN) :: SIGMA
      TYPE(C_PTR),                            INTENT(IN) :: METHOD
      TYPE(C_PTR),                            INTENT(IN) :: HURST
      TYPE(C_PTR),                            INTENT(IN) :: DH
      TYPE(C_PTR),                            INTENT(IN) :: POI
      TYPE(C_PTR),                            INTENT(IN) :: NPOI
      TYPE(C_PTR),                            INTENT(IN) :: MUTE, TAPER
      TYPE(C_PTR),                            INTENT(IN) :: RESCALE, PAD
      INTEGER(C_INT),                 POINTER            :: F_METHOD, F_RESCALE, F_PAD, F_NPOI
      REAL(C_FPP),                    POINTER            :: F_HURST, F_DH, F_MUTE, F_TAPER
      REAL(C_FPP),    DIMENSION(:),   POINTER            :: F_X, F_Y, F_Z
      REAL(C_FPP),    DIMENSION(:,:), POINTER            :: F_POI

      !-----------------------------------------------------------------------------------------------------------------------------

      STRUCTURED = 1

      DIMS(:) = FE(:) - FS(:) + 1

      NULLIFY(F_METHOD, F_HURST, F_DH, F_POI, F_NPOI, F_MUTE, F_TAPER, F_RESCALE, F_PAD)

      IF (C_ASSOCIATED(METHOD))  CALL C_F_POINTER(METHOD, F_METHOD)
      IF (C_ASSOCIATED(HURST))   CALL C_F_POINTER(HURST, F_HURST)
      IF (C_ASSOCIATED(DH))      CALL C_F_POINTER(DH, F_DH)
      IF (C_ASSOCIATED(NPOI))    CALL C_F_POINTER(NPOI, F_NPOI)
      IF (C_ASSOCIATED(POI))     CALL C_F_POINTER(POI, F_POI, [3, F_NPOI])
      IF (C_ASSOCIATED(MUTE))    CALL C_F_POINTER(MUTE, F_MUTE)
      IF (C_ASSOCIATED(TAPER))   CALL C_F_POINTER(TAPER, F_TAPER)
      IF (C_ASSOCIATED(RESCALE)) CALL C_F_POINTER(RESCALE, F_RESCALE)
      IF (C_ASSOCIATED(PAD))     CALL C_F_POINTER(PAD, F_PAD)

      CALL SCARF_INITIALIZE(FS, FE, DS, ACF, CL, SIGMA, F_METHOD, F_HURST, F_DH, F_POI, F_MUTE, F_TAPER, F_RESCALE, F_PAD)

    END SUBROUTINE SCARF_STRUCT_INITIALIZE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE SCARF_UNSTRUCT_INITIALIZE(NPTS, X, Y, Z, DH, ACF, CL, SIGMA, METHOD, HURST, POI, NPOI, MUTE, TAPER, RESCALE, PAD)  &
      BIND(C, NAME = "scarf_unstruct_initialize")

      INTEGER(C_INT),                         INTENT(IN) :: NPTS
      TYPE(C_PTR),                            INTENT(IN) :: X, Y, Z
      REAL(C_FPP),                            INTENT(IN) :: DH
      INTEGER(C_INT),                         INTENT(IN) :: ACF
      REAL(C_FPP),    DIMENSION(3),           INTENT(IN) :: CL
      TYPE(C_PTR),                            INTENT(IN) :: METHOD
      TYPE(C_PTR),                            INTENT(IN) :: HURST
      TYPE(C_PTR),                            INTENT(IN) :: POI
      TYPE(C_PTR),                            INTENT(IN) :: NPOI
      TYPE(C_PTR),                            INTENT(IN) :: MUTE, TAPER
      TYPE(C_PTR),                            INTENT(IN) :: RESCALE, PAD
      INTEGER(C_INT),                 POINTER            :: F_METHOD, F_RESCALE, F_PAD, F_NPOI
      REAL(C_FPP),                    POINTER            :: F_HURST, F_MUTE, F_TAPER
      REAL(C_FPP),    DIMENSION(:),   POINTER            :: F_X, F_Y, F_Z
      REAL(C_FPP),    DIMENSION(:,:), POINTER            :: F_POI

      !-----------------------------------------------------------------------------------------------------------------------------

      STRUCTURED = 0

      N = NPTS

      CALL C_F_POINTER(X, F_X, [NPTS])
      CALL C_F_POINTER(Y, F_Y, [NPTS])
      CALL C_F_POINTER(Z, F_Z, [NPTS])

      NULLIFY(F_METHOD, F_HURST, F_POI, F_NPOI, F_MUTE, F_TAPER, F_RESCALE, F_PAD)

      IF (C_ASSOCIATED(METHOD))  CALL C_F_POINTER(METHOD, F_METHOD)
      IF (C_ASSOCIATED(HURST))   CALL C_F_POINTER(HURST, F_HURST)
      IF (C_ASSOCIATED(NPOI))    CALL C_F_POINTER(NPOI, F_NPOI)
      IF (C_ASSOCIATED(POI))     CALL C_F_POINTER(POI, F_POI, [3, F_NPOI])
      IF (C_ASSOCIATED(MUTE))    CALL C_F_POINTER(MUTE, F_MUTE)
      IF (C_ASSOCIATED(TAPER))   CALL C_F_POINTER(TAPER, F_TAPER)
      IF (C_ASSOCIATED(RESCALE)) CALL C_F_POINTER(RESCALE, F_RESCALE)
      IF (C_ASSOCIATED(PAD))     CALL C_F_POINTER(PAD, F_PAD)

      CALL SCARF_INITIALIZE(F_X, F_Y, F_Z, DH, ACF, CL, SIGMA, F_METHOD, F_HURST, F_POI, F_MUTE, F_TAPER, F_RESCALE, F_PAD)

    END SUBROUTINE SCARF_UNSTRUCT_INITIALIZE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE SCARF_C_EXECUTE(SEED, FIELD, STATS) BIND(C, NAME="scarf_execute")

      INTEGER(C_INT),                           INTENT(IN)    :: SEED
      TYPE(C_PTR),                              INTENT(INOUT) :: FIELD
      TYPE(C_PTR),                              INTENT(INOUT) :: STATS
      REAL(C_FPP),    DIMENSION(:),     POINTER               :: F_STATS
      REAL(C_FPP),    DIMENSION(:),     POINTER               :: F_UNSTRUCT
      REAL(C_FPP),    DIMENSION(:,:,:), POINTER               :: F_STRUCT

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL C_F_POINTER(STATS, F_STATS, [8])

      IF (STRUCTURED .EQ. 0) THEN
        CALL C_F_POINTER(FIELD, F_UNSTRUCT, [N])
        CALL SCARF_EXECUTE(SEED, F_UNSTRUCT, F_STATS)
      ELSE
        CALL C_F_POINTER(FIELD, F_STRUCT, [DIMS])
        CALL SCARF_EXECUTE(SEED, F_STRUCT, F_STATS)
      ENDIF

    END SUBROUTINE SCARF_C_EXECUTE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *



    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE SCARF_FINALIZE_F2C


    END SUBROUTINE SCARF_FINALIZE_F2C

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


END MODULE SCARF3D_C_BINDING
