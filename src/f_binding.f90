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

    SUBROUTINE STRUCT_INITIALIZE(FS, FE, DS, ACF, CL, SIGMA, METHOD, HURST, DH, POI, NP, MUTE, TAPER, RESCALE, PAD)  &
    BIND(C, NAME="struct_initialize")

      INTEGER(C_INT), DIMENSION(3),              INTENT(IN) :: FS, FE
      REAL(C_FPP),                               INTENT(IN) :: DS
      INTEGER(C_INT),                            INTENT(IN) :: ACF
      REAL(C_FPP),    DIMENSION(3),              INTENT(IN) :: CL
      REAL(C_FPP),                               INTENT(IN) :: SIGMA
      INTEGER(C_INT),                  OPTIONAL, INTENT(IN) :: METHOD
      REAL(C_FPP),                     OPTIONAL, INTENT(IN) :: HURST
      REAL(C_FPP),                     OPTIONAL, INTENT(IN) :: DH
      !REAL(C_FPP),    DIMENSION(3,NP), OPTIONAL, INTENT(IN) :: POI
      TYPE(C_PTR),                     OPTIONAL, INTENT(IN) :: POI
      INTEGER(C_INT),                  OPTIONAL, INTENT(IN) :: NP
      REAL(C_FPP),                     OPTIONAL, INTENT(IN) :: MUTE, TAPER
      INTEGER(C_INT),                  OPTIONAL, INTENT(IN) :: RESCALE, PAD

      ! INTEGER(C_INT),                 POINTER            :: F_METHOD, F_RESCALE, F_PAD, F_NPOI
      ! REAL(C_FPP),                    POINTER            :: F_HURST, F_DH, F_MUTE, F_TAPER
      ! REAL(C_FPP),    DIMENSION(:),   POINTER            :: F_X, F_Y, F_Z
      REAL(C_FPP),    DIMENSION(:,:), POINTER               :: PTR

      !-----------------------------------------------------------------------------------------------------------------------------

      STRUCTURED = 1

      DIMS(:) = FE(:) - FS(:) + 1

      ! NULLIFY(F_METHOD, F_HURST, F_DH, F_POI, F_NPOI, F_MUTE, F_TAPER, F_RESCALE, F_PAD)
      !
      ! IF (C_ASSOCIATED(METHOD))  CALL C_F_POINTER(METHOD, F_METHOD)
      ! IF (C_ASSOCIATED(HURST))   CALL C_F_POINTER(HURST, F_HURST)
      ! IF (C_ASSOCIATED(DH))      CALL C_F_POINTER(DH, F_DH)
      ! IF (C_ASSOCIATED(NPOI))    CALL C_F_POINTER(NPOI, F_NPOI)
      ! IF (C_ASSOCIATED(POI))     CALL C_F_POINTER(POI, F_POI, [3, F_NPOI])
      ! IF (C_ASSOCIATED(MUTE))    CALL C_F_POINTER(MUTE, F_MUTE)
      ! IF (C_ASSOCIATED(TAPER))   CALL C_F_POINTER(TAPER, F_TAPER)
      ! IF (C_ASSOCIATED(RESCALE)) CALL C_F_POINTER(RESCALE, F_RESCALE)
      ! IF (C_ASSOCIATED(PAD))     CALL C_F_POINTER(PAD, F_PAD)

      NULLIFY(PTR)

      IF (PRESENT(POI)) CALL C_F_POINTER(POI, PTR, [3, NP])

      CALL SCARF_INITIALIZE(FS, FE, DS, ACF, CL, SIGMA, METHOD, HURST, DH, PTR, MUTE, TAPER, RESCALE, PAD)

      NULLIFY(PTR)

    END SUBROUTINE STRUCT_INITIALIZE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE UNSTRUCT_INITIALIZE(NPTS, X, Y, Z, DH, ACF, CL, SIGMA, METHOD, HURST, POI, NP, MUTE, TAPER, RESCALE, PAD)  &
    BIND(C, NAME="unstruct_initialize")

      INTEGER(C_INT),                            INTENT(IN) :: NPTS
      REAL(C_FPP),    DIMENSION(NPTS),           INTENT(IN) :: X, Y, Z
      REAL(C_FPP),                               INTENT(IN) :: DH
      INTEGER(C_INT),                            INTENT(IN) :: ACF
      REAL(C_FPP),    DIMENSION(3),              INTENT(IN) :: CL
      REAL(C_FPP),                               INTENT(IN) :: SIGMA
      INTEGER(C_INT),                  OPTIONAL, INTENT(IN) :: METHOD
      REAL(C_FPP),                     OPTIONAL, INTENT(IN) :: HURST
      !REAL(C_FPP),    DIMENSION(3,NP), OPTIONAL, INTENT(IN) :: POI
      TYPE(C_PTR),                     OPTIONAL, INTENT(IN) :: POI
      INTEGER(C_INT),                  OPTIONAL, INTENT(IN) :: NP
      REAL(C_FPP),                     OPTIONAL, INTENT(IN) :: MUTE, TAPER
      INTEGER(C_INT),                  OPTIONAL, INTENT(IN) :: RESCALE, PAD
      ! INTEGER(C_INT),                 POINTER            :: F_METHOD, F_RESCALE, F_PAD, F_NPOI
      ! REAL(C_FPP),                    POINTER            :: F_HURST, F_MUTE, F_TAPER
      ! REAL(C_FPP),    DIMENSION(:),   POINTER            :: F_X, F_Y, F_Z
      REAL(C_FPP),    DIMENSION(:,:), POINTER            :: PTR

      !-----------------------------------------------------------------------------------------------------------------------------

      STRUCTURED = 0

      N = NPTS

      ! CALL C_F_POINTER(X, F_X, [NPTS])
      ! CALL C_F_POINTER(Y, F_Y, [NPTS])
      ! CALL C_F_POINTER(Z, F_Z, [NPTS])
      !
      ! NULLIFY(F_METHOD, F_HURST, F_POI, F_NPOI, F_MUTE, F_TAPER, F_RESCALE, F_PAD)
      !
      ! IF (C_ASSOCIATED(METHOD))  CALL C_F_POINTER(METHOD, F_METHOD)
      ! IF (C_ASSOCIATED(HURST))   CALL C_F_POINTER(HURST, F_HURST)
      ! IF (C_ASSOCIATED(NPOI))    CALL C_F_POINTER(NPOI, F_NPOI)
      ! IF (C_ASSOCIATED(POI))     CALL C_F_POINTER(POI, F_POI, [3, F_NPOI])
      ! IF (C_ASSOCIATED(MUTE))    CALL C_F_POINTER(MUTE, F_MUTE)
      ! IF (C_ASSOCIATED(TAPER))   CALL C_F_POINTER(TAPER, F_TAPER)
      ! IF (C_ASSOCIATED(RESCALE)) CALL C_F_POINTER(RESCALE, F_RESCALE)
      ! IF (C_ASSOCIATED(PAD))     CALL C_F_POINTER(PAD, F_PAD)

      IF (PRESENT(POI)) CALL C_F_POINTER(POI, PTR, [3, NP])

      CALL SCARF_INITIALIZE(X, Y, Z, DH, ACF, CL, SIGMA, METHOD, HURST, PTR, MUTE, TAPER, RESCALE, PAD)

      NULLIFY(PTR)

    END SUBROUTINE UNSTRUCT_INITIALIZE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE EXECUTE(SEED, FIELD, STATS) BIND(C, NAME="execute")

      INTEGER(C_INT),                           INTENT(IN)  :: SEED
      TYPE(C_PTR),                              INTENT(OUT) :: FIELD
      REAL(C_FPP),    DIMENSION(8),             INTENT(OUT) :: STATS
      REAL(C_FPP),    DIMENSION(:),     POINTER             :: F_UNSTRUCT
      REAL(C_FPP),    DIMENSION(:,:,:), POINTER             :: F_STRUCT

      !-----------------------------------------------------------------------------------------------------------------------------


print*, 'f_binding.f90 - execute ', SEED, STRUCTURED, DIMS, N

      IF (STRUCTURED .EQ. 0) THEN
        CALL C_F_POINTER(FIELD, F_UNSTRUCT, [N])
        CALL SCARF_EXECUTE(SEED, F_UNSTRUCT, STATS)
print*, 'f_binding.f90 - execute ', MINVAL(F_UNSTRUCT), MAXVAL(F_UNSTRUCT), STATS
        NULLIFY(F_UNSTRUCT)
      ELSE
        CALL C_F_POINTER(FIELD, F_STRUCT, [DIMS])
        CALL SCARF_EXECUTE(SEED, F_STRUCT, STATS)
print*, 'f_binding.f90 - execute ', MINVAL(F_STRUCT), MAXVAL(F_STRUCT), STATS
        NULLIFY(F_STRUCT)
      ENDIF

    END SUBROUTINE EXECUTE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE FINALIZE() BIND(C, NAME="finalize")

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL SCARF_FINALIZE()

    END SUBROUTINE FINALIZE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE IO_ONE(N, FIELD, FOUT, NWRITERS) BIND(C, NAME="io_one")

      INTEGER(C_INT),                DIMENSION(3),               INTENT(IN) :: N
      TYPE(C_PTR),                                               INTENT(IN) :: FIELD
      CHARACTER(C_CHAR),             DIMENSION(*),               INTENT(IN) :: FOUT
      INTEGER(C_INT),                                  OPTIONAL, INTENT(IN) :: NWRITERS
      CHARACTER(:),      ALLOCATABLE                                        :: FILENAME
      INTEGER(IPP)                                                          :: I
      REAL(C_FPP),                   DIMENSION(:,:,:), POINTER              :: PTR

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL C_F_POINTER(FIELD, PTR, [N])

      I = 1

      FILENAME = ""

      DO WHILE(FOUT(I) .NE. C_NULL_CHAR)
        FILENAME = FILENAME // FOUT(I)
        I        = I + 1
      ENDDO

      CALL SCARF_IO(N, PTR, FILENAME, NWRITERS)

      NULLIFY(PTR)

    END SUBROUTINE IO_ONE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE IO_SLICE(N, DIRECTION, PLANE, FIELD, FOUT) BIND(C, NAME="io_slice")

      INTEGER(C_INT),                DIMENSION(3),             INTENT(IN) :: N
      INTEGER(C_INT),                                          INTENT(IN) :: DIRECTION
      INTEGER(C_INT),                                          INTENT(IN) :: PLANE
      TYPE(C_PTR),                                             INTENT(IN) :: FIELD
      CHARACTER(C_CHAR),             DIMENSION(*),             INTENT(IN) :: FOUT
      CHARACTER(:),      ALLOCATABLE                                      :: FILENAME
      INTEGER(IPP)                                                        :: I
      REAL(C_FPP),                   DIMENSION(:,:,:), POINTER            :: PTR

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL C_F_POINTER(FIELD, PTR, [N])

      I = 1

      FILENAME = ""

      DO WHILE(FOUT(I) .NE. C_NULL_CHAR)
        FILENAME = FILENAME // FOUT(I)
        I        = I + 1
      ENDDO

      CALL SCARF_IO(N, DIRECTION, PLANE, PTR, FILENAME)

      NULLIFY(PTR)

    END SUBROUTINE IO_SLICE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE SCARF3D_C_BINDING
