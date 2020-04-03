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

    SUBROUTINE STRUCT_INITIALIZE(FS, FE, DS, ACF, CL, SIGMA, METHOD, HURST, DH, POI, NP, MUTE, TAPER, RESCALE, PAD, NC, FC)  &
    BIND(C, NAME="struct_initialize")

      INTEGER(C_INT), DIMENSION(3),              INTENT(IN) :: FS, FE
      REAL(C_FPP),                               INTENT(IN) :: DS
      INTEGER(C_INT),                            INTENT(IN) :: ACF
      REAL(C_FPP),    DIMENSION(3),              INTENT(IN) :: CL
      REAL(C_FPP),                               INTENT(IN) :: SIGMA
      INTEGER(C_INT),                  OPTIONAL, INTENT(IN) :: METHOD
      REAL(C_FPP),                     OPTIONAL, INTENT(IN) :: HURST
      REAL(C_FPP),                     OPTIONAL, INTENT(IN) :: DH
      TYPE(C_PTR),                     OPTIONAL, INTENT(IN) :: POI
      INTEGER(C_INT),                  OPTIONAL, INTENT(IN) :: NP
      REAL(C_FPP),                     OPTIONAL, INTENT(IN) :: MUTE, TAPER
      INTEGER(C_INT),                  OPTIONAL, INTENT(IN) :: RESCALE, PAD
      REAL(C_FPP),    DIMENSION(3),    OPTIONAL, INTENT(IN) :: NC, FC
      REAL(C_FPP),    DIMENSION(:,:),  POINTER              :: PTR
      REAL(C_FPP),    DIMENSION(0,0),  TARGET               :: VOID

      !-----------------------------------------------------------------------------------------------------------------------------

      STRUCTURED = 1

      DIMS(:) = FE(:) - FS(:) + 1

      IF (PRESENT(POI)) THEN
        CALL C_F_POINTER(POI, PTR, [3, NP])
      ELSE
        PTR => VOID
      ENDIF

      CALL SCARF_INITIALIZE(FS, FE, DS, ACF, CL, SIGMA, METHOD, HURST, DH, PTR, MUTE, TAPER, RESCALE, PAD, NC, FC)

      NULLIFY(PTR)

    END SUBROUTINE STRUCT_INITIALIZE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE UNSTRUCT_INITIALIZE(NPTS, X, Y, Z, DH, ACF, CL, SIGMA, METHOD, HURST, DS, POI, NP, MUTE, TAPER, RESCALE, PAD, NC,  &
                                   FC) BIND(C, NAME="unstruct_initialize")

      INTEGER(C_INT),                            INTENT(IN) :: NPTS
      REAL(C_FPP),    DIMENSION(NPTS),           INTENT(IN) :: X, Y, Z
      REAL(C_FPP),                               INTENT(IN) :: DH
      INTEGER(C_INT),                            INTENT(IN) :: ACF
      REAL(C_FPP),    DIMENSION(3),              INTENT(IN) :: CL
      REAL(C_FPP),                               INTENT(IN) :: SIGMA
      INTEGER(C_INT),                  OPTIONAL, INTENT(IN) :: METHOD
      REAL(C_FPP),                     OPTIONAL, INTENT(IN) :: HURST
      REAL(C_FPP),                     OPTIONAL, INTENT(IN) :: DS
      TYPE(C_PTR),                     OPTIONAL, INTENT(IN) :: POI
      INTEGER(C_INT),                  OPTIONAL, INTENT(IN) :: NP
      REAL(C_FPP),                     OPTIONAL, INTENT(IN) :: MUTE, TAPER
      INTEGER(C_INT),                  OPTIONAL, INTENT(IN) :: RESCALE, PAD
      REAL(C_FPP),    DIMENSION(3),    OPTIONAL, INTENT(IN) :: NC, FC
      REAL(C_FPP),    DIMENSION(:,:),  POINTER              :: PTR
      REAL(C_FPP),    DIMENSION(0,0),  TARGET               :: VOID

      !-----------------------------------------------------------------------------------------------------------------------------

      STRUCTURED = 0

      N = NPTS

      IF (PRESENT(POI)) THEN
        CALL C_F_POINTER(POI, PTR, [3, NP])
      ELSE
        PTR => VOID
      ENDIF

      CALL SCARF_INITIALIZE(X, Y, Z, DH, ACF, CL, SIGMA, METHOD, HURST, DS, PTR, MUTE, TAPER, RESCALE, PAD, NC, FC)

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

      IF (STRUCTURED .EQ. 0) THEN
        CALL C_F_POINTER(FIELD, F_UNSTRUCT, [N])
        CALL SCARF_EXECUTE(SEED, F_UNSTRUCT, STATS)
        NULLIFY(F_UNSTRUCT)
      ELSE
        CALL C_F_POINTER(FIELD, F_STRUCT, [DIMS])
        CALL SCARF_EXECUTE(SEED, F_STRUCT, STATS)
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

    SUBROUTINE IO_ONE(NPTS, FIELD, FOUT, NWRITERS) BIND(C, NAME="io_one")

      INTEGER(C_INT),                DIMENSION(3),               INTENT(IN) :: NPTS
      TYPE(C_PTR),                                               INTENT(IN) :: FIELD
      CHARACTER(C_CHAR),             DIMENSION(*),               INTENT(IN) :: FOUT
      INTEGER(C_INT),                                  OPTIONAL, INTENT(IN) :: NWRITERS
      CHARACTER(:),      ALLOCATABLE                                        :: FILENAME
      INTEGER(IPP)                                                          :: I
      REAL(C_FPP),                   DIMENSION(:,:,:), POINTER              :: PTR

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL C_F_POINTER(FIELD, PTR, [DIMS])

      I = 1

      FILENAME = ""

      DO WHILE(FOUT(I) .NE. C_NULL_CHAR)
        FILENAME = FILENAME // FOUT(I)
        I        = I + 1
      ENDDO

      CALL SCARF_IO(NPTS, PTR, FILENAME, NWRITERS)

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

      CALL C_F_POINTER(FIELD, PTR, [DIMS])

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
