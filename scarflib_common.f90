MODULE SCARFLIB_COMMON

  USE, INTRINSIC     :: ISO_FORTRAN_ENV
  USE, INTRINSIC     :: ISO_C_BINDING
  USE, NON_INTRINSIC :: MPI

  IMPLICIT NONE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  ! ALL VARIABLES/SUBROUTINE ARE ACCESSIBLE FROM OTHER MODULES VIA HOST-ASSOCIATION
  PUBLIC

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTERFACE VARIANCE
    MODULE PROCEDURE VARIANCE_1D, VARIANCE_3D
  END INTERFACE

  INTERFACE MEAN
    MODULE PROCEDURE MEAN_1D, MEAN_3D
  END INTERFACE

  INTERFACE TAPERING
    MODULE PROCEDURE TAPERING_UNSTRUCTURED, TAPERING_STRUCTURED
  END INTERFACE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  ! SET PRECISION

  ! INTEGERS HAVE ALWAYS SAME "PRECISION"
  INTEGER, PARAMETER :: IPP = INT32

#ifdef DOUBLE_PREC
  INTEGER(IPP), PARAMETER :: FPP          = REAL64
  INTEGER(IPP), PARAMETER :: C_FPP        = C_DOUBLE
  INTEGER(IPP), PARAMETER :: C_CPP        = C_DOUBLE_COMPLEX
  INTEGER(IPP), PARAMETER :: REAL_TYPE    = MPI_DOUBLE_PRECISION
  INTEGER(IPP), PARAMETER :: COMPLEX_TYPE = MPI_DOUBLE_COMPLEX
#else
  INTEGER(IPP), PARAMETER :: FPP          = REAL32
  INTEGER(IPP), PARAMETER :: C_FPP        = C_FLOAT
  INTEGER(IPP), PARAMETER :: C_CPP        = C_FLOAT_COMPLEX
  INTEGER(IPP), PARAMETER :: REAL_TYPE    = MPI_REAL
  INTEGER(IPP), PARAMETER :: COMPLEX_TYPE = MPI_COMPLEX
#endif

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  ! TOTAL NUMBER OF POINTS ALONG EACH AXIS (WHOLE MODEL)
  INTEGER(IPP),                DIMENSION(3)             :: NPTS

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

    SUBROUTINE WATCH_START(TICTOC, COMM)

      ! START THE MPI-STOPWATCH

      REAL(REAL64),           INTENT(OUT) :: TICTOC                            !< INITIAL TIME
      INTEGER(IPP), OPTIONAL, INTENT(IN)  :: COMM
      INTEGER(IPP)                        :: IERR                              !< MPI STUFF

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (.NOT.PRESENT(COMM)) THEN
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
      ELSE
        CALL MPI_BARRIER(COMM, IERR)
      ENDIF

      TICTOC = MPI_WTIME()

    END SUBROUTINE WATCH_START

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE WATCH_STOP(TICTOC, COMM)

      ! STOP THE MPI-STOPWATCH AND RETURN ELAPSED TIME

      REAL(REAL64),           INTENT(INOUT) :: TICTOC                          !< ELAPSED TIME
      INTEGER(IPP), OPTIONAL, INTENT(IN)    :: COMM
      INTEGER(IPP)                          :: IERR                            !< MPI STUFF

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (.NOT.PRESENT(COMM)) THEN
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
      ELSE
        CALL MPI_BARRIER(COMM, IERR)
      ENDIF

      TICTOC = MPI_WTIME() - TICTOC

    END SUBROUTINE WATCH_STOP

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    ! SUBROUTINE SET_STREAM(SEED)
    !
    !   ! INITIALSE THE RANDOM NUMBER GENERATOR
    !
    !   INTEGER(IPP),                           INTENT(IN) :: SEED                      !< USER-DEFINED SEED NUMBER
    !   INTEGER(IPP)                                       :: SEED_SIZE                 !< STORE ARRAY SIZE
    !   INTEGER(IPP), ALLOCATABLE, DIMENSION(:)            :: TMP_SEED                  !< USED TO INITIALISE THE RANDOM GENERATOR
    !
    !   !-----------------------------------------------------------------------------------------------------------------------------
    !
    !   CALL RANDOM_SEED(SIZE = SEED_SIZE)
    !
    !   ALLOCATE(TMP_SEED(SEED_SIZE))
    !
    !   TMP_SEED = SEED
    !
    !   CALL RANDOM_SEED(PUT = TMP_SEED)
    !
    !   DEALLOCATE(TMP_SEED)
    !
    ! END SUBROUTINE SET_STREAM

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE WRITE_SINGLE(V, FILENAME, NWRITERS)

      ! EACH PROCESS WRITES ITS OWN FILE WITH DATA. ONLY "NWRITERS" FILES ARE CREATED AT THE SAME TIME.

      REAL(FPP),                    DIMENSION(:,:,:),          INTENT(IN) :: V                                  !< RANDOM FIELD
      CHARACTER(LEN=*),                                        INTENT(IN) :: FILENAME                           !< NAME OF OUTPUT FILE
      INTEGER(IPP),                                            INTENT(IN) :: NWRITERS
      CHARACTER(:),     ALLOCATABLE                                       :: MYNAME
      INTEGER(IPP)                                                        :: I, N                          !< COUNTERS
      INTEGER(IPP)                                                        :: FH, IERR                          !< MPI STUFF
      INTEGER(IPP),                 DIMENSION(MPI_STATUS_SIZE)            :: STATUS                        !< MPI STUFF

      !-----------------------------------------------------------------------------------------------------------------------------

      DO I = 0, WORLD_SIZE - 1, NWRITERS

        N = I + NWRITERS

        ! ONLY PROCESSES IN THE RANGE [I, I + NWRITERS), I.E. MAX "NWRITERS", ENTER THIS BLOCK EVERY TIME
        IF ( (WORLD_RANK .GE. I) .AND. (WORLD_RANK .LT. N) ) THEN

          MYNAME = FILENAME // '_x=' // NUM2STR(GS(1, WORLD_RANK)) // '_' // NUM2STR(GE(1, WORLD_RANK)) //    &
                               '_y=' // NUM2STR(GS(2, WORLD_RANK)) // '_' // NUM2STR(GE(2, WORLD_RANK)) //    &
                               '_z=' // NUM2STR(GS(3, WORLD_RANK)) // '_' // NUM2STR(GE(3, WORLD_RANK))

          CALL MPI_FILE_OPEN(MPI_COMM_SELF, MYNAME, MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, FH, IERR)

          CALL MPI_FILE_WRITE(FH, V, SIZE(V), REAL_TYPE, STATUS, IERR)

          CALL MPI_FILE_CLOSE(FH, IERR)

        ENDIF

        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

      ENDDO

    END SUBROUTINE WRITE_SINGLE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE WRITE_ONE(V, FILENAME, NWRITERS)

      ! "NWRITERS" PROCESSES WRITE DATA OF ALL PROCESSES TO A SINGLE FILE

      REAL(FPP),                     DIMENSION(:,:,:),          INTENT(IN) :: V                                  !< RANDOM FIELD
      CHARACTER(LEN=*),                                         INTENT(IN) :: FILENAME                           !< NAME OF OUTPUT FILE
      INTEGER(IPP),                                             INTENT(IN) :: NWRITERS                           !< I/O PROCESSES
      CHARACTER(:),     ALLOCATABLE                                        :: DCHAR                              !< DUMMY
      INTEGER(IPP)                                                         :: I, J, N                            !< COUNTER
      INTEGER(IPP)                                                         :: NMAX, MAXTASKS                     !< COUNTER
      INTEGER(IPP)                                                         :: COLOR, NTASKS, NEWCOMM, WRTCOMM    !< MPI STUFF
      INTEGER(IPP)                                                         :: RANK, RCVRANK, REQUEST, INFO       !< MPI STUFF
      INTEGER(IPP)                                                         :: NEWTYPE, IERR, FH                  !< MPI STUFF
      INTEGER(IPP),                  DIMENSION(MPI_STATUS_SIZE)            :: STATUS                             !< MPI STUFF
      INTEGER(IPP),                  DIMENSION(3)                          :: SUBSIZES, STARTS                   !< ARRAYS FOR DERIVED DATATYPE
      INTEGER(IPP),     ALLOCATABLE, DIMENSION(:)                          :: RANK0, RANK1                       !< FIRST/LAST GLOBAL RANK
      INTEGER(IPP),     ALLOCATABLE, DIMENSION(:)                          :: NEW2WORLD                          !< MAP LOCAL TO GLOBAL RANKS
      REAL(FPP),        ALLOCATABLE, DIMENSION(:)                          :: BUFFER                             !< COLLECTED DATA

      !-----------------------------------------------------------------------------------------------------------------------------

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ORGANIZE PROCESSES IN "NWRITERS" COMMUNICATORS AND PREPARE A MAP OF THE LOCAL-TO-GLOBAL RANKS
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ALLOCATE(RANK0(NWRITERS), RANK1(NWRITERS))

      CALL SPLIT_TASK(WORLD_SIZE, NWRITERS, RANK0, RANK1)

      ! LOWEST RANK MUST BE ZERO
      RANK0 = RANK0 - 1
      RANK1 = RANK1 - 1

      DO I = 1, NWRITERS
        IF ( (WORLD_RANK .GE. RANK0(I)) .AND. (WORLD_RANK .LE. RANK1(I)) ) COLOR = I
      ENDDO

      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, COLOR, WORLD_RANK, NEWCOMM, IERR)

      CALL MPI_COMM_RANK(NEWCOMM, RANK, IERR)
      CALL MPI_COMM_SIZE(NEWCOMM, NTASKS, IERR)

      ALLOCATE(NEW2WORLD(0:NTASKS - 1))

      ! I-TH RANK IN "NEWCOMM" HAS ITS GLOBAL RANK CONTAINED IN "NEW2WORLD"
      CALL MAP_RANKS(NEWCOMM, MPI_COMM_WORLD, NEW2WORLD)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! PROCESSES WITH "RANK=0" IN EVERY "NEWCOMM" ARE WRITERS: GATHER THEM IN SEPARATE COMMUNICATOR "WRTCOMM"
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      COLOR = 0

      IF (RANK .EQ. 0) COLOR = 1

      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, COLOR, WORLD_RANK, WRTCOMM, IERR)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! DETERMINE MAXIMUM NUMBER OF TASKS AMONGST ALL "NEWCOMM" COMMUNICATORS, NEEDED TO FOR COLLECTIVE WRITE CALL BELOW
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      MAXTASKS = NTASKS

      CALL MPI_ALLREDUCE(MPI_IN_PLACE, MAXTASKS, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, IERR)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! DETERMINE MAXIMUM SIZE OF INPUT ARRAY "V" AND ALLOCATE "BUFFER" ACCORDINGLY. THIS WAY WE DON'T NEED TO REALLOCATE "BUFFER"
      ! EVERY TIME WE GET DATA FROM A PROCESS
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      N = SIZE(V)

      NMAX = N

      CALL MPI_ALLREDUCE(MPI_IN_PLACE, NMAX, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, IERR)

      ALLOCATE(BUFFER(NMAX))

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! EACH PROCESS MAKES NON-BLOCKING SEND. I/O PROCESSES COLLECT DATA AND WRITE THEM TO FILE AT THE RIGHT LOCATION. A WRITER MAY
      ! WRITE TWICE THE SAME DATA IF "WORLD_SIZE" NOT MULTIPLE OF "NWRITERS": THIS IS NECESSARY BECAUSE WE USE COLLECTIVE WRITE CALL
      ! THAT MUST BE CALLED BY ALL WRITERS AT THE SAME TIME.
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! EACH PROCESS SENDS ITS DATA TO WRITER INSIDE COMMUNICATOR "NEWCOMM"
      CALL MPI_ISEND(V, N, REAL_TYPE, 0, 1, NEWCOMM, REQUEST, IERR)

      ! ONLY WRITERS CAN WORK ON FILES
      IF (RANK .EQ. 0) THEN

        INFO = MPI_INFO_NULL

        ! SET "STRIPING_FACTOR" EQUAL TO NUMBER OF I/O PROCESSES. "STRIPING_UNIT" IS ARBITRAY SET TO 128MB (REASONABLE VALUE).
#ifdef PFS
        CALL MPI_INFO_CREATE(INFO, IERR)

        DCHAR = NUM2STR(NWRITERS)

        CALL MPI_INFO_SET(INFO, "striping_factor", DCHAR, IERR)

        DEALLOCATE(DCHAR)

        DCHAR = NUM2STR(128*1024*1024)

        CALL MPI_INFO_SET(INFO, "striping_unit", DCHAR, IERR)

        DEALLOCATE(DCHAR)
#endif

        CALL MPI_FILE_OPEN(WRTCOMM, FILENAME, MPI_MODE_CREATE + MPI_MODE_WRONLY, INFO, FH, IERR)

        ! LOOP UNTIL ALL DATA ARE RECEIVED
        DO J = 0, MAXTASKS - 1

          ! STOP RECEIVING IF WE HAVE ALREADY COLLECTED ALL DATA FOR "NEWCOMM"
          IF (J .LE. NTASKS - 1) THEN

            CALL MPI_RECV(BUFFER, SIZE(BUFFER), REAL_TYPE, MPI_ANY_SOURCE, MPI_ANY_TAG, NEWCOMM, STATUS, IERR)

            ! FIND WHO SENT DATA
            RCVRANK = STATUS(MPI_SOURCE)

            ! MAP ITS LOCAL RANK INTO GLOBAL ONE
            RCVRANK = NEW2WORLD(RCVRANK)

          ENDIF

          SUBSIZES(:) = GE(:, RCVRANK) - GS(:, RCVRANK) + 1
          STARTS(:)   = GS(:, RCVRANK) - 1

          CALL MPI_TYPE_CREATE_SUBARRAY(3, NPTS, SUBSIZES, STARTS, MPI_ORDER_FORTRAN, REAL_TYPE, NEWTYPE, IERR)

          CALL MPI_TYPE_COMMIT(NEWTYPE, IERR)

          ! VIEW IS CONTROLLED BY DERIVED DATATYPE
          CALL MPI_FILE_SET_VIEW(FH, 0_MPI_OFFSET_KIND, REAL_TYPE, NEWTYPE, 'NATIVE', MPI_INFO_NULL, IERR)

          CALL MPI_FILE_WRITE_ALL(FH, BUFFER, PRODUCT(SUBSIZES), REAL_TYPE, MPI_STATUS_IGNORE, IERR)

          CALL MPI_TYPE_FREE(NEWTYPE, IERR)

        ENDDO

        CALL MPI_FILE_CLOSE(FH, IERR)

#ifdef PFS
        CALL MPI_INFO_FREE(INFO, IERR)
#endif

      ENDIF

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! FREE RESOURCES
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      DEALLOCATE(RANK0, RANK1, NEW2WORLD, BUFFER)

      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

      CALL MPI_COMM_FREE(NEWCOMM, IERR)
      CALL MPI_COMM_FREE(WRTCOMM, IERR)

    END SUBROUTINE WRITE_ONE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION NUM2STR(V)

      ! CONVERT AN INTEGER INTO A STRING

      INTEGER(IPP),                       INTENT(IN) :: V
      CHARACTER(:),           ALLOCATABLE            :: NUM2STR
      CHARACTER(RANGE(V) + 2)                        :: DUM

      !-----------------------------------------------------------------------------------------------------------------------------

      WRITE(DUM, '(I0)') V

      NUM2STR = TRIM(DUM)

    END FUNCTION NUM2STR

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE SPLIT_TASK(NPTS, NTASKS, I0, I1)

      ! DISTRIBUTE "NPTS" POINTS AMONGST "NTASKS" PROCESSES AND, FOR EACH PROCESS, RETURN THE LOWEST AND HIGHEST INDEX.

      INTEGER(IPP),                        INTENT(IN)  :: NPTS                          !< NUMBER OF POINTS TO BE SPLIT
      INTEGER(IPP),                        INTENT(IN)  :: NTASKS                     !< NUMBER OF MPI PROCESSES
      INTEGER(IPP), DIMENSION(0:NTASKS-1), INTENT(OUT) :: I0, I1                     !< 1ST/LAST INDEX
      INTEGER(IPP)                                     :: P                          !< COUNTER

      !------------------------------------------------------------------------------------------------------------------------------

      DO P = 0, NTASKS - 1
        I0(P) = 1 + INT( REAL(NPTS, FPP) / REAL(NTASKS, FPP) * REAL(P, FPP) )
        I1(P) = INT( REAL(NPTS, FPP) / REAL(NTASKS, FPP) * REAL(P + 1, FPP) )
      ENDDO

    END SUBROUTINE SPLIT_TASK

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE MAP_RANKS(SRC_COMM, DEST_COMM, DEST_RANKS)

      ! MAP RANKS FROM LOCAL TO GLOBAL COMMUNICATOR

      INTEGER(IPP),                         INTENT(IN)    :: SRC_COMM                     !< COMMUNICATOR
      INTEGER(IPP),                         INTENT(IN)    :: DEST_COMM
      INTEGER(IPP), DIMENSION(:),           INTENT(INOUT) :: DEST_RANKS                   !< RANKS AFTER BEING MAPPED TO GLOBAL GROUP
      INTEGER(IPP)                                        :: I                        !< COUNTER
      INTEGER(IPP)                                        :: IERR, SRC_GROUP              !< MPI STUFF
      INTEGER(IPP)                                        :: DEST_GROUP
      INTEGER(IPP)                                        :: NTASKS
      INTEGER(IPP), DIMENSION(SIZE(DEST_RANKS))           :: SRC_RANKS

      !-----------------------------------------------------------------------------------------------------------------------------

      NTASKS = SIZE(SRC_RANKS)

      ! RANKS IN SOURCE COMMUNICATOR
      DO I = 1, NTASKS
        SRC_RANKS(I) = I - 1
      ENDDO

      ! GROUP ASSOCIATED TO SOURCE COMMUNICATOR
      CALL MPI_COMM_GROUP(SRC_COMM, SRC_GROUP, IERR)

      ! GROUP ASSOCIATED TO DESTINATION COMMUNICATOR
      CALL MPI_COMM_GROUP(DEST_COMM, DEST_GROUP, IERR)

      ! MAP RANKS IN "GROUP" INTO RANKS IN "MPI_GROUP_WORLD"
      CALL MPI_GROUP_TRANSLATE_RANKS(SRC_GROUP, NTASKS, SRC_RANKS, DEST_GROUP, DEST_RANKS, IERR)

      CALL MPI_GROUP_FREE(SRC_GROUP, IERR)
      CALL MPI_GROUP_FREE(DEST_GROUP, IERR)

    END SUBROUTINE MAP_RANKS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE WRITE_SLICE(DIRECTION, PLANE, V, FILENAME)

      ! CREATE A SINGLE FILE CONTAINING A SLICE OF THE RANDOM FIELD. EACH PROCESS WRITE ITS OWN PART BY MAKING USE OF THE SUBARRAY
      ! DATATYPE. A SLICE CAN BE ORIENTED ALONG THE X-, Y- OR Z-AXIS.

      INTEGER(IPP),                                    INTENT(IN) :: DIRECTION                          !< DIRECTION OF DESIRED 2D SLICE
      INTEGER(IPP),                                    INTENT(IN) :: PLANE                              !< INDEX OF SLICE IN GLOBAL COORDINATES
      REAL(FPP),                     DIMENSION(:,:,:), INTENT(IN) :: V                                  !< RANDOM FIELD
      CHARACTER(LEN=*),                                INTENT(IN) :: FILENAME                           !< NAME OF OUTPUT FILE
      INTEGER(IPP)                                                :: I, J, C                            !< COUNTERS
      INTEGER(IPP)                                                :: NEWTYPE, IERR, FH, COMM, COLOR     !< MPI STUFF
      INTEGER(IPP),                  DIMENSION(2)                 :: SUBSIZES, STARTS, DIMS             !< ARRAYS FOR DERIVED DATATYPE
      LOGICAL                                                     :: BOOL
      REAL(FPP),        ALLOCATABLE, DIMENSION(:,:)               :: BUFFER                             !< STORE PART OF THE SLICE

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

        DIMS = [NPTS(2), NPTS(3)]

        C = PLANE - GS(1, WORLD_RANK) + 1

      ! SLICE ALONG Y
      ELSEIF (DIRECTION .EQ. 2) THEN

        BOOL = (PLANE .GE. GS(2, WORLD_RANK)) .AND. (PLANE .LE. GE(2, WORLD_RANK))

        SUBSIZES = [SIZE(V, 1), SIZE(V, 3)]
        STARTS   = [GS(1, WORLD_RANK) - 1, GS(3, WORLD_RANK) - 1]

        DIMS = [NPTS(1), NPTS(3)]

        C = PLANE - GS(2, WORLD_RANK) + 1

      ! SLICE ALONG Z
      ELSEIF (DIRECTION .EQ. 3) THEN

        BOOL = (PLANE .GE. GS(3, WORLD_RANK)) .AND. (PLANE .LE. GE(3, WORLD_RANK))

        SUBSIZES = [SIZE(V, 1), SIZE(V, 2)]
        STARTS   = [GS(1, WORLD_RANK) - 1, GS(2, WORLD_RANK) - 1]

        DIMS = [NPTS(1), NPTS(2)]

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

        ! SET FILE VIEW BASED ON SUBARRAY DATATYPE
        CALL MPI_FILE_SET_VIEW(FH, 0_MPI_OFFSET_KIND, REAL_TYPE, NEWTYPE, 'NATIVE', MPI_INFO_NULL, IERR)

        CALL MPI_FILE_WRITE_ALL(FH, BUFFER, PRODUCT(SUBSIZES), REAL_TYPE, MPI_STATUS_IGNORE, IERR)

        CALL MPI_FILE_CLOSE(FH, IERR)

        ! RELEASE DERIVED DATATYPE
        CALL MPI_TYPE_FREE(NEWTYPE, IERR)

        ! REALEASE LOCAL COMMUNICATOR
        CALL MPI_COMM_FREE(COMM, IERR)

        DEALLOCATE(BUFFER)

      ENDIF

      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

    END SUBROUTINE WRITE_SLICE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(FPP) FUNCTION VARIANCE_1D(R)

      ! COMPUTE VARIANCE BASED ON THE COMPENSATED-SUMMATION VERSION OF THE TWO-PASS ALGORITHM. CALCULATIONS ARE ALWAYS IN DOUBLE
      ! PRECISION.

      REAL(FPP),   DIMENSION(:), INTENT(IN) :: R
      INTEGER(IPP)                          :: I
      INTEGER(IPP)                          :: N
      REAL(REAL64)                          :: MU, S1, S2, X, V

      !-----------------------------------------------------------------------------------------------------------------------------

      N = SIZE(R)

      MU = MEAN(R)

      S1 = 0._REAL64
      S2 = 0._REAL64

      DO I = 1, N
        X = REAL(R(I), REAL64) - MU
        S1 = S1 + X
        S2 = S2 + X**2
      ENDDO

      S1 = (S1**2) / REAL(N, REAL64)

      V = (S2 - S1) / REAL(N - 1, REAL64)

      VARIANCE_1D = REAL(V, FPP)

    END FUNCTION VARIANCE_1D

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(FPP) FUNCTION VARIANCE_3D(R)

      ! COMPUTE VARIANCE BASED ON THE COMPENSATED-SUMMATION VERSION OF THE TWO-PASS ALGORITHM. CALCULATIONS ARE ALWAYS IN DOUBLE
      ! PRECISION.

      REAL(FPP),   DIMENSION(:,:,:), INTENT(IN) :: R
      INTEGER(IPP)                              :: I, J, K
      INTEGER(IPP)                              :: NX, NY, NZ
      REAL(REAL64)                              :: MU, S1, S2, X, V

      !-----------------------------------------------------------------------------------------------------------------------------

      NX = SIZE(R, 1)
      NY = SIZE(R, 2)
      NZ = SIZE(R, 3)

      MU = MEAN(R)

      S1 = 0._REAL64
      S2 = 0._REAL64

      DO K = 1, NZ
        DO J = 1, NY
          DO I = 1, NX
            X = REAL(R(I, J, K), REAL64) - MU
            S1 = S1 + X
            S2 = S2 + X**2
          ENDDO
        ENDDO
      ENDDO

      S1 = (S1**2) / REAL(SIZE(R), REAL64)

      V = (S2 - S1) / REAL(SIZE(R) - 1, REAL64)

      VARIANCE_3D = REAL(V, FPP)

    END FUNCTION VARIANCE_3D

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(FPP) FUNCTION MEAN_1D(R)

      ! COMPUTE MEAN OF DATASET "R" AVOIDING FLOATING POINT INACCURACIES. CALCULATIONS ARE ALWAYS IN DOUBLE PRECISION.

      REAL(FPP),   DIMENSION(:), INTENT(IN) :: R
      INTEGER(IPP)                          :: I
      INTEGER(IPP)                          :: N
      REAL(REAL64)                          :: V, C

      !-------------------------------------------------------------------------------------------------------------------------------

      N = SIZE(R)

      V = 0._REAL64
      C = 1._REAL64

      DO I = 1, N
        V = V + (REAL(R(I), REAL64) - V) / C
        C = C + 1._REAL64
      ENDDO

      ! RETURN MEAN AT DESIRED PRECISION
      MEAN_1D = REAL(V, FPP)

    END FUNCTION MEAN_1D

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(FPP) FUNCTION MEAN_3D(R)

      ! COMPUTE MEAN OF DATASET "R" AVOIDING FLOATING POINT INACCURACIES. CALCULATIONS ARE ALWAYS IN DOUBLE PRECISION.

      REAL(FPP),   DIMENSION(:,:,:), INTENT(IN) :: R
      INTEGER(IPP)                              :: I, J, K
      INTEGER(IPP)                              :: NX, NY, NZ
      REAL(REAL64)                              :: V, C

      !-------------------------------------------------------------------------------------------------------------------------------

      NX = SIZE(R, 1)
      NY = SIZE(R, 2)
      NZ = SIZE(R, 3)

      V = 0._REAL64
      C = 1._REAL64

      DO K = 1, NZ
        DO J = 1, NY
          DO I = 1, NX
            V = V + (REAL(R(I, J, K), REAL64) - V) / C
            C = C + 1._REAL64
          ENDDO
        ENDDO
      ENDDO

      ! RETURN MEAN AT DESIRED PRECISION
      MEAN_3D = REAL(V, FPP)

    END FUNCTION MEAN_3D

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE PARALLEL_VARIANCE(VARSET, AVGSET, NSET, VAR, AVG)

      ! GIVEN A SEQUENCE OF VARIANCE "VARSET" AND MEAN "AVGSET" VALUES, EACH BASED ON "NSET" NUMBER OF POINTS, THIS SUBROUTINE RETURNS
      ! THE RESULTING TOTAL VARIANCE AND MEAN BASED ON THE ALGORITHM OF CHAN ET AL. (1979). CALCULATIONS ARE ALWAYS IN DOUBLE PRECISION.

      REAL(FPP),    DIMENSION(:), INTENT(IN)  :: VARSET                             !< SET OF VARIANCE VALUES
      REAL(FPP),    DIMENSION(:), INTENT(IN)  :: AVGSET                             !< SET OF MEAN VALUES
      INTEGER(IPP), DIMENSION(:), INTENT(IN)  :: NSET                               !< NUMBER OF POINTS FOR EACH SET
      REAL(FPP),                  INTENT(OUT) :: VAR                                !< RESULTING VARIANCE
      REAL(FPP),                  INTENT(OUT) :: AVG                                !< RESULTING AVERAGE
      INTEGER(IPP)                            :: I                                  !< COUNTER
      REAL(REAL64)                            :: SIG, MU, N, DELTA, M1, M2, M       !< LOCAL VARIABLES

      !-----------------------------------------------------------------------------------------------------------------------------

      ! VARIANCE, MEAN AND NUMBER OF POINTS OF FIRST SET
      SIG = REAL(VARSET(1), REAL64)
      MU  = REAL(AVGSET(1), REAL64)
      N   = REAL(NSET(1), REAL64)

      ! LOOP OVER SET OF VARIANCE/MEAN VALUES
      DO I = 2, SIZE(NSET)

        DELTA = REAL(AVGSET(I), REAL64) - MU

        M1 = SIG * (N - 1._REAL64)
        M2 = REAL(VARSET(I), REAL64) * (NSET(I) - 1._REAL64)

        M = M1 + M2 + DELTA**2 * N * REAL(NSET(I), REAL64) / (N + REAL(NSET(I), REAL64))

        ! RESULTING MEAN
        MU = (MU * N + REAL(AVGSET(I), REAL64) * REAL(NSET(I), REAL64)) / (N + REAL(NSET(I), REAL64))

        ! RESULTING NUMBER OF POINTS
        N = N + REAL(NSET(I), REAL64)

        ! RESULTING VARIANCE
        SIG = M / (N - 1._REAL64)

      ENDDO

      ! RETURN WITH DESIRED PRECISION
      VAR = REAL(SIG, FPP)
      AVG = REAL(MU, FPP)

    END SUBROUTINE PARALLEL_VARIANCE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE TAPERING_STRUCTURED (DH, FS, FE, FIELD, POI, MUTE, TAPER)

      ! MUTE AND/OR TAPER THE RANDOM FIELD AROUND A POINT-OF-INTEREST. MUTING OCCURS WITHIN A RADIUS OF "MUTE" POINTS; TAPERING IS
      ! ACHIEVED BY APPLYING A HANNING WINDOW WITHIN A RADIUS IN THE RANGE "MUTE + 1" AND "MUTE + TAPER" POINTS.

      REAL(FPP),                                                      INTENT(IN)    :: DH
      INTEGER(IPP), DIMENSION(3),                                     INTENT(IN)    :: FS, FE                !< START/END INDICES ALONG EACH DIRECTION
      REAL(FPP),    DIMENSION(FS(1):FE(1), FS(2):FE(2), FS(3):FE(3)), INTENT(INOUT) :: FIELD                !< RANDOM FIELD
      REAL(FPP),    DIMENSION(3),                                     INTENT(IN)    :: POI                   !< POINT-OF-INTEREST (TRASLATED)
      REAL(FPP),                                                      INTENT(IN)    :: MUTE, TAPER           !< MUTE/TAPER LENGTH (IN POINTS)
      INTEGER(IPP)                                                                  :: I, J, K               !< COUNTERS
      REAL(FPP)                                                                     :: D, DS                 !< VARIOUS DISTANCES
      REAL(FPP)                                                                     :: DX, DY, DZ            !< NODE-POI DISTANCE
      REAL(FPP)                                                                     :: T                     !< TAPER PARAMETER

      !-----------------------------------------------------------------------------------------------------------------------------

      ! TOTAL RADIUS OF TAPERED AND/OR MUTED VOLUME
      DS = MUTE + TAPER

      DO K = FS(3), FE(3)

        DZ = ((K - 1) * DH - POI(3))**2                                                !< DISTANCE ALONG Z FROM "POI"

        DO J = FS(2), FE(2)

          DY = ((J - 1) * DH - POI(2))**2                                              !< DISTANCE ALONG Y FROM "POI"

          DO I = FS(1), FE(1)

            DX = ((I - 1) * DH - POI(1))**2                                            !< DISTANCE ALONG X FROM "POI"

            D = SQRT(DX + DY + DZ)                                                     !< TOTAL DISTANCE

            ! MUTE IF DISTANCE NODE-POI IS BELOW "MUTE"
            IF (D .LE. MUTE) THEN

              FIELD(I, J, K) = 0._FPP

            ! TAPER IF DISTANCE NODE-POI IS BETWEEN "MUTE" AND "TAPER + MUTE"
            ELSEIF ( (D .GT. MUTE) .AND. (D .LE. DS) ) THEN

              T = (D - MUTE) / TAPER                                                   !< TAPER PARAMETER

              FIELD(I, J, K) = FIELD(I, J, K) * (0.5_FPP - 0.5_FPP * COS(T * PI))

            ENDIF

          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE TAPERING_STRUCTURED

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE TAPERING_UNSTRUCTURED(X, Y, Z, FIELD, POI, MUTE, TAPER)

      ! MUTE AND/OR TAPER THE RANDOM FIELD AROUND A POINT-OF-INTEREST. MUTING OCCURS WITHIN A RADIUS OF "MUTE" POINTS; TAPERING IS
      ! ACHIEVED BY APPLYING A HANNING WINDOW WITHIN A RADIUS IN THE RANGE "MUTE + 1" AND "MUTE + TAPER" POINTS.

      REAL(FPP),    DIMENSION(:), INTENT(IN)    :: X, Y, Z
      REAL(FPP),    DIMENSION(:), INTENT(INOUT) :: FIELD                 !< RANDOM FIELD
      REAL(FPP),    DIMENSION(3), INTENT(IN)    :: POI                   !< POINT-OF-INTEREST (NODES)
      REAL(FPP),                  INTENT(IN)    :: MUTE, TAPER           !< MUTE/TAPER LENGTH (IN POINTS)
      INTEGER(IPP)                              :: I                     !< COUNTERS
      REAL(FPP)                                 :: D, DS                 !< VARIOUS DISTANCES
      REAL(FPP)                                 :: T                     !< TAPER PARAMETER

      !-----------------------------------------------------------------------------------------------------------------------------

      ! TOTAL RADIUS OF TAPERED AND/OR MUTED VOLUME
      DS = MUTE + TAPER

      DO I = 1, SIZE(FIELD)

        ! DISTANCE BETWEEN GRID POINT AND "POI"
        D = SQRT( (X(I) - POI(1))**2 + (Y(I) - POI(2))**2 + (Z(I) - POI(3))**2 )

        ! MUTE IF DISTANCE NODE-POI IS BELOW "MUTE"
        IF (D .LE. MUTE) THEN

          FIELD(I) = 0._FPP

        ! TAPER IF DISTANCE NODE-POI IS BETWEEN "MUTE" AND "MUTE + TAPER"
        ELSEIF ( (D .GT. MUTE) .AND. (D .LE. DS) ) THEN

          T = (D - MUTE) / TAPER                                             !< TAPER PARAMETER

          FIELD(I) = FIELD(I) * (0.5_FPP - 0.5_FPP * COS(T * PI))

        ENDIF

      ENDDO

    END SUBROUTINE TAPERING_UNSTRUCTURED

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE BEST_CONFIG(DIMS)

      ! COMPUTE BEST GRID OF PROCESSES

      INTEGER(IPP),              DIMENSION(3), INTENT(OUT) :: DIMS
      !INTEGER(IPP),                            INTENT(OUT) :: NITER

      INTEGER(IPP)                                         :: I, J, K, L, C
      INTEGER(IPP)                                         :: N2, N3 !, N
      INTEGER(IPP)                                         :: A, B
      INTEGER(IPP),              DIMENSION(3)              :: V1, V2
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:,:)            :: FACT3, FACT2
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:,:)            :: LIST

      REAL(FPP) :: N, NITER

      !-----------------------------------------------------------------------------------------------------------------------------

      !NITER = HUGE(1)
      NITER = HUGE(1._FPP)

      C = 0

      ! FACTORISE NUMBER OF AVAILABLE PROCESSES
      CALL FACTORIZATION(WORLD_SIZE, FACT3)

      N3 = SIZE(FACT3, 2)

      ALLOCATE(LIST(3, N3 * 5))

      LIST = 0

      ! LOOP OVER FACTORISED PROCESSES
      DO L = 1, N3

        DO K = 1, 2

          IF (K .EQ. 1) THEN
            A = FACT3(1, L)
            B = FACT3(2, L)
          ELSE
            B = FACT3(1, L)
            A = FACT3(2, L)
          ENDIF

          CALL FACTORIZATION(B, FACT2)

          N2 = SIZE(FACT2, 2)

          DO J = 1, N2

            V1 = [A, FACT2(:, J)]

            ! SKIP TO NEXT PROCESSES GRID IF ALREADY ANALYSED
            IF (MATCH(V1, C, LIST) .EQV. .TRUE.) CYCLE

            C = C + 1

            ! ADD TO LIST
            LIST(:, C) = V1

            DO I = 0, 2

              V1 = CSHIFT(V1, 1)

              CALL TEST_CONFIG(V1, N)

              IF (N .LT. NITER) THEN
                DIMS  = V1
                NITER = N
              ENDIF

              V2 = [V1(1), V1(3), V1(2)]

              CALL TEST_CONFIG(V2, N)

              IF (N .LT. NITER) THEN
                DIMS  = V2
                NITER = N
              ENDIF

            ENDDO
            ! END PERMUTATIONS

          ENDDO
          ! END LOOP OVER FACTOR PAIRS FOR "A/B"

        ENDDO
        ! END LOOP OVER "A/B"

      ENDDO
      ! END LOOP OVER FACTOR PAIRS FOR "WORLD_SIZE"

      DEALLOCATE(FACT3, FACT2, LIST)

      !-----------------------------------------------------------------------------------------------------------------------------
      !-----------------------------------------------------------------------------------------------------------------------------

      CONTAINS

      LOGICAL FUNCTION MATCH(VEC, IMAX, LIST)

        INTEGER(IPP), DIMENSION(3),   INTENT(IN) :: VEC
        INTEGER(IPP),                 INTENT(IN) :: IMAX
        INTEGER(IPP), DIMENSION(:,:), INTENT(IN) :: LIST
        INTEGER(IPP)                             :: I

        !---------------------------------------------------------------------------------------------------------------------------

        DO I = 1, IMAX
          MATCH = ANY(V1(1) .EQ. LIST(:, I)) .AND. ANY(V1(2) .EQ. LIST(:, I)) .AND. ANY(V1(3) .EQ. LIST(:, I))
          IF (MATCH .EQV. .TRUE.) EXIT
        ENDDO

      END FUNCTION MATCH

    END SUBROUTINE BEST_CONFIG

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE SCARFLIB_COMMON
