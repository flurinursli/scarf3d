MODULE SCARFLIB

  ! THIS MODULE IS INTENDED TO BE USED AS WRAPPER FOR THE SCARF3D LIBRARY

  USE, NON_INTRINSIC :: MPI
  USE, NON_INTRINSIC :: SCARFLIB_COMMON
  USE, NON_INTRINSIC :: SCARFLIB_FFT
  USE, NON_INTRINSIC :: SCARFLIB_SPEC


  IMPLICIT NONE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! DECLARE ALL PRIVATE, INCLUDING ALSO PUBLIC VARIABLES AND SUBROUTINES ACCESSED FROM OTHER MODULES VIA HOST-ASSOCIATION
  PRIVATE

  ! THESE ARE THE ONLY SUBROUTINES/VARIABLES ACCESSIBLE BY CALLING PROGRAM
  PUBLIC :: SCARF_INITIALIZE, SCARF_EXECUTE, SCARF_FINALIZE, SCARF_IO, SCARF_DECOMP
  PUBLIC :: SCARF_OBJ

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTERFACE SCARF_INITIALIZE
    MODULE PROCEDURE INITIALIZE_STRUCTURED, INITIALIZE_UNSTRUCTURED
  END INTERFACE

  INTERFACE SCARF_EXECUTE
    MODULE PROCEDURE EXECUTE_STRUCTURED, EXECUTE_UNSTRUCTURED
  END INTERFACE

  INTERFACE SCARF_IO
    MODULE PROCEDURE SCARF_IO_SLICE, SCARF_IO_ONE
  END INTERFACE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  TYPE SCARF_OBJ

    INTEGER(IPP)                                       :: ACF
    INTEGER(IPP)                                       :: RESCALE
    INTEGER(IPP)                                       :: PAD
    INTEGER(IPP)                                       :: METHOD
    INTEGER(IPP),                       DIMENSION(3)   :: FS, FE
    REAL(FPP)                                          :: DS, DH
    REAL(FPP)                                          :: CL, SIGMA, HURST
    REAL(FPP)                                          :: MUTE, TAPER
    REAL(FPP),                 POINTER, DIMENSION(:)   :: X  => NULL(), Y => NULL(), Z => NULL()
    REAL(FPP),    ALLOCATABLE,          DIMENSION(:)   :: STATS
    REAL(FPP),    ALLOCATABLE,          DIMENSION(:,:) :: POI

  END TYPE SCARF_OBJ

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION INITIALIZE_STRUCTURED(FS, FE, DS, ACF, CL, SIGMA, METHOD, HURST, DH, POI, MUTE, TAPER, RESCALE, PAD) RESULT(OBJ)

      INTEGER(IPP),   DIMENSION(3),             INTENT(IN) :: FS, FE
      REAL(FPP),                                INTENT(IN) :: DS
      INTEGER(IPP),                             INTENT(IN) :: ACF
      REAL(FPP),                                INTENT(IN) :: CL, SIGMA
      INTEGER(IPP),                   OPTIONAL, INTENT(IN) :: METHOD
      REAL(FPP),                      OPTIONAL, INTENT(IN) :: HURST
      REAL(FPP),                      OPTIONAL, INTENT(IN) :: DH                    !< MAXIMUM GRID-STEP (CONTROLS MAX WAVENUMBER)
      REAL(FPP),      DIMENSION(:,:), OPTIONAL, INTENT(IN) :: POI
      REAL(FPP),                      OPTIONAL, INTENT(IN) :: MUTE, TAPER
      INTEGER(IPP),                   OPTIONAL, INTENT(IN) :: RESCALE, PAD
      TYPE(SCARF_OBJ)                                      :: OBJ

      !-----------------------------------------------------------------------------------------------------------------------------

      OBJ%FS  = FS
      OBJ%FE  = FE
      OBJ%DS  = DS

      OBJ%DH    = DH
      OBJ%ACF   = ACF
      OBJ%CL    = CL
      OBJ%SIGMA = SIGMA

      OBJ%METHOD = 0                      !< DEFAULT IS FFT METHOD

      OBJ%MUTE    = -999._FPP
      OBJ%TAPER   = -999._FPP
      OBJ%RESCALE = 0
      OBJ%PAD     = 0

      ALLOCATE(OBJ%POI(0,0))

      IF (.PRESENT(METHOD))  OBJ%METHOD  = METHOD
      IF (.PRESENT(HURST))   OBJ%HURST   = HURST
      IF (.PRESENT(DH))      OBJ%DH      = DH
      IF (.PRESENT(MUTE))    OBJ%MUTE    = MUTE
      IF (.PRESENT(TAPER))   OBJ%TAPER   = TAPER
      IF (.PRESENT(RESCALE)) OBJ%RESCALE = RESCALE
      IF (.PRESENT(PAD))     OBJ%PAD     = PAD
      IF (.PRESENT(POI))     OBJ%POI     = POI

      IF (OBJ%METHOD .EQ. 0) THEN
        ALLOCATE(OBJ%STATS(8))
      ELSEIF (OBJ%METHOD .EQ. 1) THEN
        ALLOCATE(OBJ%STATS(6))
      ENDIF

    END FUNCTION INITIALIZE_STRUCTURED

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION INITIALIZE_UNSTRUCTURED(X, Y, Z, DH, ACF, CL, SIGMA, METHOD, HURST, POI, MUTE, TAPER, RESCALE, PAD) RESULT(OBJ)

      REAL(FPP),      DIMENSION(:),   TARGET,           INTENT(IN) :: X, Y, Z
      REAL(FPP),                                        INTENT(IN) :: DH                 !< MAXIMUM GRID-STEP (CONTROLS MAX WAVENUMBER)
      INTEGER(IPP),                                     INTENT(IN) :: ACF
      REAL(FPP),                                        INTENT(IN) :: CL, SIGMA
      INTEGER(IPP),                           OPTIONAL, INTENT(IN) :: METHOD
      REAL(FPP),                              OPTIONAL, INTENT(IN) :: HURST
      REAL(FPP),      DIMENSION(:,:),         OPTIONAL, INTENT(IN) :: POI
      REAL(FPP),                              OPTIONAL, INTENT(IN) :: MUTE, TAPER
      INTEGER(IPP),                           OPTIONAL, INTENT(IN) :: RESCALE, PAD
      TYPE(SCARF_OBJ)                                              :: OBJ

      !-----------------------------------------------------------------------------------------------------------------------------

      OBJ%X     => X
      OBJ%Y     => Y
      OBJ%Z     => Z

      OBJ%DH    = DH
      OBJ%ACF   = ACF
      OBJ%CL    = CL
      OBJ%SIGMA = SIGMA

      OBJ%METHOD = 0                      !< DEFAULT IS FFT METHOD

      OBJ%MUTE    = -999._FPP
      OBJ%TAPER   = -999._FPP
      OBJ%RESCALE = 0
      OBJ%PAD     = 0

      ALLOCATE(OBJ%POI(0,0))

      IF (.PRESENT(METHOD))  OBJ%METHOD  = METHOD
      IF (.PRESENT(HURST))   OBJ%HURST   = HURST
      IF (.PRESENT(MUTE))    OBJ%MUTE    = MUTE
      IF (.PRESENT(TAPER))   OBJ%TAPER   = TAPER
      IF (.PRESENT(RESCALE)) OBJ%RESCALE = RESCALE
      IF (.PRESENT(PAD))     OBJ%PAD     = PAD
      IF (.PRESENT(POI))     OBJ%POI     = POI

      IF (OBJ%METHOD .EQ. 0) THEN
        ALLOCATE(OBJ%STATS(8))
      ELSEIF (OBJ%METHOD .EQ. 1) THEN
        ALLOCATE(OBJ%STATS(6))
      ENDIF

    END FUNCTION INITIALIZE_UNSTRUCTURED

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE EXECUTE_UNSTRUCTURED(OBJ, SEED, FIELD)

      TYPE(SCARF_OBJ),               INTENT(INOUT) :: OBJ
      INTEGER(IPP),                  INTENT(IN)    :: SEED
      REAL(FPP),       DIMENSION(:), INTENT(OUT)   :: FIELD

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (OBJ%METHOD .EQ. 0) THEN

        CALL SCARF3D_FFT(OBJ%X, OBJ%Y, OBJ%Z, OBJ%DH, OBJ%ACF, OBJ%CL, OBJ%SIGMA, OBJ%HURST, SEED, OBJ%POI, OBJ%MUTE, OBJ%TAPER,  &
                         OBJ%RESCALE, OBJ%PAD, FIELD, OBJ%STATS)

      ELSEIF (OBJ%METHOD .EQ. 1) THEN

        CALL SCARF3D_SPEC(OBJ%X, OBJ%Y, OBJ%Z, OBJ%DH, OBJ%ACF, OBJ%CL, OBJ%SIGMA, OBJ%HURST, SEED, OBJ%POI, OBJ%MUTE, OBJ%TAPER, &
                          FIELD, OBJ%STATS)

      ENDIF

    END SUBROUTINE EXECUTE_UNSTRUCTURED

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE EXECUTE_STRUCTURED(OBJ, SEED, FIELD)

      TYPE(SCARF_OBJ),                INTENT(INOUT) :: OBJ
      INTEGER(IPP),                   INTENT(IN)    :: SEED
      REAL(FPP),    DIMENSION(:,:,:), INTENT(OUT)   :: FIELD

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (OBJ%METHOD .EQ. 0) THEN

        CALL SCARF3D_FFT(OBJ%DS, OBJ%FS, OBJ%FE, OBJ%DH, OBJ%ACF, OBJ%CL, OBJ%SIGMA, OBJ%HURST, SEED, OBJ%POI, OBJ%MUTE,    &
                         OBJ%TAPER, OBJ%RESCALE, OBJ%PAD, FIELD, OBJ%STATS)

      ELSEIF (OBJ%METHOD .EQ. 1) THEN

        CALL SCARF3D_SPEC(OBJ%DS, OBJ%FS, OBJ%FE, OBJ%DH, OBJ%ACF, OBJ%CL, OBJ%SIGMA, OBJ%HURST, SEED, OBJ%POI, OBJ%MUTE,   &
                          OBJ%TAPER, FIELD, OBJ%STATS)

      ENDIF

    END SUBROUTINE EXECUTE_STRUCTURED

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE SCARF_FINALIZE(OBJ)

      TYPE(SCARF_OBJ), INTENT(IN) :: OBJ

      !-----------------------------------------------------------------------------------------------------------------------------

      DEALLOCATE(OBJ%POI)
      DEALLOCATE(OBJ%STATS)

      NULLIFY(OBJ%X, OBJ%Y, OBJ%Z)

    END SUBROUTINE SCARF_FINALIZE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE SCARF_IO_SLICE(DIRECTION, PLANE, FIELD, FILENAME)

      INTEGER(IPP),                       INTENT(IN) :: DIRECTION
      INTEGER(IPP),                       INTENT(IN) :: PLANE
      REAL(FPP),        DIMENSION(:,:,:), INTENT(IN) :: FIELD
      CHARACTER(LEN=*),                   INTENT(IN) :: FILENAME
      INTEGER(IPP)                                   :: RANK, NP, IERR

      !-----------------------------------------------------------------------------------------------------------------------------

      ! GET RANK NUMBER
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, RANK, IERR)

      ! GET NUMBER OF TASKS
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NP, IERR)

      ! THESE ARE GLOBAL VARIABLES STORED IN "SCARFLIB_COMMON"
      ALLOCATE(GS(3, NP), GE(3, NP))

      ! STORE GLOBAL INDICES
      GS(:, RANK) = OBJ%FS
      GE(:, RANK) = OBJ%FE

      ! SHARE INFO AMONGST ALL PROCESSES
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, GS, 3, MPI_INTEGER, MPI_COMM_WORLD, IERR)
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, GE, 3, MPI_INTEGER, MPI_COMM_WORLD, IERR)

      CALL WRITE_SLICE(DIRECTION, PLANE, FIELD, FILENAME)

      DEALLOCATE(GS, GE)

    END SUBROUTINE SCARF_IO_SLICE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE SCARF_IO_ONE(FIELD, FILENAME, NWRITERS)

      REAL(FPP),        DIMENSION(:,:,:), INTENT(IN) :: FIELD
      CHARACTER(LEN=*),                   INTENT(IN) :: FILENAME
      INTEGER(IPP),                       INTENT(IN) :: NWRITERS
      INTEGER(IPP)                                   :: RANK, NP, IERR

      !-----------------------------------------------------------------------------------------------------------------------------

      ! GET RANK NUMBER
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, RANK, IERR)

      ! GET NUMBER OF TASKS
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NP, IERR)

      ! THESE ARE GLOBAL VARIABLES STORED IN "SCARFLIB_COMMON"
      ALLOCATE(GS(3, NP), GE(3, NP))

      ! STORE GLOBAL INDICES
      GS(:, RANK) = OBJ%FS
      GE(:, RANK) = OBJ%FE

      ! SHARE INFO AMONGST ALL PROCESSES
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, GS, 3, MPI_INTEGER, MPI_COMM_WORLD, IERR)
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, GE, 3, MPI_INTEGER, MPI_COMM_WORLD, IERR)

      CALL WRITE_ONE(FIELD, FILENAME, NWRITERS)

      DEALLOCATE(GS, GE)

    END SUBROUTINE SCARF_IO_ONE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE SCARFLIB
