MODULE SCARFLIB

  USE, INTRINSIC     :: ISO_FORTRAN_ENV
  USE, INTRINSIC     :: ISO_C_BINDING
  USE, NON_INTRINSIC :: MPI

  IMPLICIT NONE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  !PRIVATE

  !PUBLIC :: SCARF3D_FFT, SCARF3D_SPEC

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  ! SET PRECISION FOR FLOATING POINT NUMBERS

!#ifdef DOUBLE_PREC
 INTEGER, PARAMETER :: FPP          = REAL64
 INTEGER, PARAMETER :: C_FPP        = C_DOUBLE
 INTEGER, PARAMETER :: C_CPP        = C_DOUBLE_COMPLEX
 INTEGER, PARAMETER :: REAL_TYPE    = MPI_DOUBLE_PRECISION
 INTEGER, PARAMETER :: COMPLEX_TYPE = MPI_DOUBLE_COMPLEX
!#else
  ! INTEGER, PARAMETER :: FPP          = REAL32
  ! INTEGER, PARAMETER :: C_FPP        = C_FLOAT
  ! INTEGER, PARAMETER :: C_CPP        = C_FLOAT_COMPLEX
  ! INTEGER, PARAMETER :: REAL_TYPE    = MPI_REAL
  ! INTEGER, PARAMETER :: COMPLEX_TYPE = MPI_COMPLEX
!#endif

  ! INTEGERS HAVE ALWAYS SAME "PRECISION"
  INTEGER, PARAMETER :: IPP = INT32

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  ! INTERFACES
  ! INTERFACE
  !   MODULE SUBROUTINE SCARF3D_FFT(COMM, N, DH, ACF, CL, SIGMA, HURST, SEED, POI, MUTE, TAPER)
  !     INTEGER(IPP),                      INTENT(IN) :: COMM               !< MPI COMMUNICATOR
  !     INTEGER(IPP),     DIMENSION(3),    INTENT(IN) :: N                  !< MODEL SIZE: NUMBER OF POINTS ALONG X, Y, Z
  !     REAL(FPP),                         INTENT(IN) :: DH                 !< GRID-STEP
  !     CHARACTER(LEN=*),                  INTENT(IN) :: ACF                !< AUTOCORRELATION FUNCTION: "VK" OR "GAUSS"
  !     REAL(FPP),        DIMENSION(3),    INTENT(IN) :: CL                 !< CORRELATION LENGTH
  !     REAL(FPP),                         INTENT(IN) :: SIGMA              !< STANDARD DEVIATION
  !     REAL(FPP),                         INTENT(IN) :: HURST              !< HURST EXPONENT
  !     INTEGER(IPP),                      INTENT(IN) :: SEED               !< SEED NUMBER
  !     INTEGER(IPP),     DIMENSION(:,:),  INTENT(IN) :: POI                !< LOCATION OF POINT(S)-OF-INTEREST
  !     INTEGER(IPP),                      INTENT(IN) :: MUTE               !< NUMBER OF POINTS WHERE MUTING IS APPLIED
  !     INTEGER(IPP),                      INTENT(IN) :: TAPER              !< NUMBER OF POINTS WHERE TAPERING IS APPLIED
  !   END SUBROUTINE SCARF3D_FFT
  !
  !   MODULE SUBROUTINE SCARF3D_SPEC(COMM, N, DH, ACF, CL, SIGMA, HURST, SEED, POI, MUTE, TAPER)
  !     INTEGER(IPP),                      INTENT(IN) :: COMM               !< MPI COMMUNICATOR
  !     INTEGER(IPP),     DIMENSION(3),    INTENT(IN) :: N                  !< MODEL SIZE: NUMBER OF POINTS ALONG X, Y, Z
  !     REAL(FPP),                         INTENT(IN) :: DH                 !< GRID-STEP
  !     CHARACTER(LEN=*),                  INTENT(IN) :: ACF                !< AUTOCORRELATION FUNCTION: "VK" OR "GAUSS"
  !     REAL(FPP),        DIMENSION(3),    INTENT(IN) :: CL                 !< CORRELATION LENGTH
  !     REAL(FPP),                         INTENT(IN) :: SIGMA              !< STANDARD DEVIATION
  !     REAL(FPP),                         INTENT(IN) :: HURST              !< HURST EXPONENT
  !     INTEGER(IPP),                      INTENT(IN) :: SEED               !< SEED NUMBER
  !     INTEGER(IPP),     DIMENSION(:,:),  INTENT(IN) :: POI                !< LOCATION OF POINT(S)-OF-INTEREST
  !     INTEGER(IPP),                      INTENT(IN) :: MUTE               !< NUMBER OF POINTS WHERE MUTING IS APPLIED
  !     INTEGER(IPP),                      INTENT(IN) :: TAPER              !< NUMBER OF POINTS WHERE TAPERING IS APPLIED
  !   END SUBROUTINE SCARF3D_SPEC
  ! END INTERFACE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  ! TOTAL NUMBER OF POINTS ALONG EACH AXIS (WHOLE MODEL)
  INTEGER(IPP),                DIMENSION(3)             :: N

  ! MPI STUFF
  INTEGER(IPP)                                          :: WORLD_RANK, WORLD_SIZE

  ! GLOBAL MODEL INDICES FOR ALL MPI PROCESSES
  INTEGER(IPP),   ALLOCATABLE, DIMENSION(:,:)           :: GS, GE

  ! VARIABLES
  REAL(FPP),                                  PARAMETER :: PI = 3.141592653589793_FPP

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE WATCH_START(TICTOC)

      ! START THE MPI-STOPWATCH

      REAL(FPP),   INTENT(OUT) :: TICTOC                            !< INITIAL TIME
      INTEGER(IPP)             :: IERR                              !< MPI STUFF

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

      TICTOC = MPI_WTIME()

    END SUBROUTINE WATCH_START

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE WATCH_STOP(TICTOC)

      ! STOP THE MPI-STOPWATCH AND RETURN ELAPSED TIME

      REAL(FPP),   INTENT(INOUT) :: TICTOC                          !< ELAPSED TIME
      INTEGER(IPP)               :: IERR                            !< MPI STUFF

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

      TICTOC = MPI_WTIME() - TICTOC

    END SUBROUTINE WATCH_STOP

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE SET_STREAM(SEED)

      ! INITIALSE THE RANDOM NUMBER GENERATOR

      INTEGER(IPP),                           INTENT(IN) :: SEED                      !< USER-DEFINED SEED NUMBER
      INTEGER(IPP)                                       :: SEED_SIZE                 !< STORE ARRAY SIZE
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)            :: TMP_SEED                  !< USED TO INITIALISE THE RANDOM GENERATOR

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL RANDOM_SEED(SIZE = SEED_SIZE)

      ALLOCATE(TMP_SEED(SEED_SIZE))

      TMP_SEED = SEED

      CALL RANDOM_SEED(PUT = TMP_SEED)

      DEALLOCATE(TMP_SEED)

    END SUBROUTINE SET_STREAM

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE IO_WRITE_ONE(V, FILENAME)

      ! CREATE A SINGLE FILE CONTAINING THE WHOLE RANDOM FIELD. EACH PROCESS WRITE ITS OWN PART BY MAKING USE OF THE SUBARRAY DATATYPE.

      REAL(FPP),                    DIMENSION(:,:,:), INTENT(IN) :: V                                  !< RANDOM FIELD
      CHARACTER(LEN=*),                               INTENT(IN) :: FILENAME                           !< NAME OF OUTPUT FILE
      INTEGER(IPP)                                               :: I                                  !< COUNTER
      INTEGER(IPP)                                               :: NEWTYPE, IERR, FH                  !< MPI STUFF
      INTEGER(KIND=MPI_OFFSET_KIND)                              :: FILESIZE, OFFSET                   !< MPI STUFF
      INTEGER(IPP),                 DIMENSION(3)                 :: SUBSIZES, STARTS                   !< ARRAYS FOR DERIVED DATATYPE

      !-----------------------------------------------------------------------------------------------------------------------------

      ! SET SUBSIZE AND STARTING INDICES FOR SUBARRAY DATATYPE
      DO I = 1, 3
        SUBSIZES(I) = SIZE(V, I)
        STARTS(I)   = GS(I, WORLD_RANK) - 1
      ENDDO

      ! CREATE SUBARRAY
      CALL MPI_TYPE_CREATE_SUBARRAY(3, N, SUBSIZES, STARTS, MPI_ORDER_FORTRAN, REAL_TYPE, NEWTYPE, IERR)

      CALL MPI_TYPE_COMMIT(NEWTYPE, IERR)

      CALL MPI_FILE_OPEN(MPI_COMM_WORLD, FILENAME, MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, FH, IERR)

      FILESIZE = 0_MPI_OFFSET_KIND

      ! THIS SUBROUTINE IS USED TO GUARANTEE OVERWRITING
      CALL MPI_FILE_SET_SIZE(FH, FILESIZE, IERR)

      OFFSET = 0_MPI_OFFSET_KIND

      ! SET FILE VIEW BASED ON SUBARRAY DATATYPE
      CALL MPI_FILE_SET_VIEW(FH, OFFSET, REAL_TYPE, NEWTYPE, 'NATIVE', MPI_INFO_NULL, IERR)

      CALL MPI_FILE_WRITE_ALL(FH, V, PRODUCT(SUBSIZES), REAL_TYPE, MPI_STATUS_IGNORE, IERR)

      CALL MPI_FILE_CLOSE(FH, IERR)

      ! RELEASE DERIVED DATATYPE
      CALL MPI_TYPE_FREE(NEWTYPE, IERR)

    END SUBROUTINE IO_WRITE_ONE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE IO_WRITE_SLICE(DIRECTION, PLANE, V, FILENAME)

      ! CREATE A SINGLE FILE CONTAINING A SLICE OF THE RANDOM FIELD. EACH PROCESS WRITE ITS OWN PART BY MAKING USE OF THE SUBARRAY
      ! DATATYPE. A SLICE CAN BE ORIENTED ALONG THE X-, Y- OR Z-AXIS.

      INTEGER(IPP),                                                INTENT(IN) :: DIRECTION                          !< DIRECTION OF DESIRED 2D SLICE
      INTEGER(IPP),                                                INTENT(IN) :: PLANE                              !< INDEX OF SLICE IN GLOBAL COORDINATES
      REAL(FPP),                                 DIMENSION(:,:,:), INTENT(IN) :: V                                  !< RANDOM FIELD
      CHARACTER(LEN=*),                                            INTENT(IN) :: FILENAME                           !< NAME OF OUTPUT FILE
      INTEGER(IPP)                                                            :: I, J, C                            !< COUNTERS
      INTEGER(IPP)                                                            :: NEWTYPE, IERR, FH, COMM, COLOR     !< MPI STUFF
      INTEGER(KIND=MPI_OFFSET_KIND)                                           :: FILESIZE, OFFSET                   !< MPI STUFF
      INTEGER(IPP),                              DIMENSION(2)                 :: SUBSIZES, STARTS, DIMS             !< ARRAYS FOR DERIVED DATATYPE
      LOGICAL                                                                 :: BOOL
      REAL(FPP),                    ALLOCATABLE, DIMENSION(:,:)               :: BUFFER                             !< STORE PART OF THE SLICE

      !-----------------------------------------------------------------------------------------------------------------------------

      ! CREAT A NEW COMMUNICATOR ONLY WITH PROCESSES CONTAINING DATA AT "PLANE". THESE WILL WRITE DATA TO DISK.
      COLOR = 0

      ! PREPARE SOME VARIABLES FOR WRITING TO DISK. "BOOL" DETERMINE IF A PROCESS BELONG TO THE LOCAL COMMUNICATOR. "C" IS USED FOR
      ! GLOBAL TO LOCAL INDEX MAPPING.
      ! SLICE ALONG X
      IF (DIRECTION .EQ. 1) THEN

        BOOL = (PLANE .GE. GS(1, WORLD_RANK)) .AND. (PLANE .LE. GE(1, WORLD_RANK))

        SUBSIZES = [SIZE(V, 2), SIZE(V, 3)]
        STARTS   = [GS(2, WORLD_RANK) - 1, GS(3, WORLD_RANK) - 1]

        DIMS = [N(2), N(3)]

        C = PLANE - GS(1, WORLD_RANK) + 1

      ! SLICE ALONG Y
      ELSEIF (DIRECTION .EQ. 2) THEN

        BOOL = (PLANE .GE. GS(2, WORLD_RANK)) .AND. (PLANE .LE. GE(2, WORLD_RANK))

        SUBSIZES = [SIZE(V, 1), SIZE(V, 3)]
        STARTS   = [GS(1, WORLD_RANK) - 1, GS(3, WORLD_RANK) - 1]

        DIMS = [N(1), N(3)]

        C = PLANE - GS(2, WORLD_RANK) + 1

      ! SLICE ALONG Z
      ELSEIF (DIRECTION .EQ. 3) THEN

        BOOL = (PLANE .GE. GS(3, WORLD_RANK)) .AND. (PLANE .LE. GE(3, WORLD_RANK))

        SUBSIZES = [SIZE(V, 1), SIZE(V, 2)]
        STARTS   = [GS(1, WORLD_RANK) - 1, GS(2, WORLD_RANK) - 1]

        DIMS = [N(1), N(2)]

        C = PLANE - GS(3, WORLD_RANK) + 1

      ENDIF

      ! UPDATE COLOR FOR CURRENT PROCESS IF IT CONTAINS PART OF THE SLICE
      IF (BOOL) COLOR = 1

      ! SPLIT COMMUNICATOR
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, COLOR, WORLD_RANK, COMM, IERR)

      ! ONLY PROCESSES IN THE NEW COMMUNICATOR ENTER SECTION BELOW
      IF (COLOR .EQ. 1) THEN

        ALLOCATE(BUFFER(SUBSIZES(1), SUBSIZES(2)))

        ! HERE BELOW WE COPY THE RIGHT DATA INTO "BUFFER"

        IF (DIRECTION .EQ. 1) THEN

          DO J = 1, SUBSIZES(2)
            DO I = 1, SUBSIZES(1)
              BUFFER(I, J) = V(C, I, J)
            ENDDO
          ENDDO

        ELSEIF (DIRECTION .EQ. 2) THEN

          DO J = 1, SUBSIZES(2)
            DO I = 1, SUBSIZES(1)
              BUFFER(I, J) = V(I, C, J)
            ENDDO
          ENDDO

        ELSEIF (DIRECTION .EQ. 3) THEN

          DO J = 1, SUBSIZES(2)
            DO I = 1, SUBSIZES(1)
              BUFFER(I, J) = V(I, J, C)
            ENDDO
          ENDDO

        ENDIF

        ! CREATE SUBARRAY
        CALL MPI_TYPE_CREATE_SUBARRAY(2, DIMS, SUBSIZES, STARTS, MPI_ORDER_FORTRAN, REAL_TYPE, NEWTYPE, IERR)

        CALL MPI_TYPE_COMMIT(NEWTYPE, IERR)

        CALL MPI_FILE_OPEN(COMM, FILENAME, MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, FH, IERR)

        FILESIZE = 0_MPI_OFFSET_KIND

        ! THIS SUBROUTINE IS USED TO GUARANTEE OVERWRITING
        CALL MPI_FILE_SET_SIZE(FH, FILESIZE, IERR)

        OFFSET = 0_MPI_OFFSET_KIND

        ! SET FILE VIEW BASED ON SUBARRAY DATATYPE
        CALL MPI_FILE_SET_VIEW(FH, OFFSET, REAL_TYPE, NEWTYPE, 'NATIVE', MPI_INFO_NULL, IERR)

        CALL MPI_FILE_WRITE_ALL(FH, BUFFER, PRODUCT(SUBSIZES), REAL_TYPE, MPI_STATUS_IGNORE, IERR)

        CALL MPI_FILE_CLOSE(FH, IERR)

        ! RELEASE DERIVED DATATYPE
        CALL MPI_TYPE_FREE(NEWTYPE, IERR)

        ! REALEASE LOCAL COMMUNICATOR
        CALL MPI_COMM_FREE(COMM, IERR)

        DEALLOCATE(BUFFER)

      ENDIF

      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

    END SUBROUTINE IO_WRITE_SLICE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE SCARFLIB
