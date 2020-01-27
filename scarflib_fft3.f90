MODULE SCARFLIB_FFT3

  ! ALL VARIABLES AND SUBMODULE PROCEDURES ARE GLOBAL WITHIN THE SUBMODULE, BUT LIMITED TO IT.

  ! mpif90 -O3 scarflib_fft2.f90 scarflib_spec.f90 scarflib.f90 driver.f90 *.f -I/home/walter/Backedup/Software/fftw-3.3.8/include
  ! -L/home/walter/Backedup/Software/fftw-3.3.8/lib -lfftw3 -fcheck=all -fbacktrace -fopenacc

  USE, NON_INTRINSIC ::  SCARFLIB_COMMON

  IMPLICIT NONE

  INCLUDE 'fftw3.f03'

  PRIVATE

  PUBLIC :: SCARF3D_FFT

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! PROCEDURE POINTERS
  PROCEDURE(VK_PSDF), POINTER :: FUN

  INTERFACE

    FUNCTION PRNG(SEED, LS, LE, NPTS) BIND(C, NAME="prng")
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT),                           VALUE :: SEED
      INTEGER(C_INT),  DIMENSION(3), INTENT(IN)       :: LS, LE
      INTEGER(C_INT),  DIMENSION(3), INTENT(IN)       :: NPTS
      TYPE(C_PTR)                                     :: PRNG
    END FUNCTION PRNG

    SUBROUTINE FREE_MEM(PTR) BIND(C, NAME="free_mem")
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      TYPE(C_PTR), VALUE, INTENT(IN) :: PTR
    END SUBROUTINE FREE_MEM

  END INTERFACE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! CARTESIAN GRID MUST BE 3-DIMENSIONAL
  INTEGER(IPP),                             PARAMETER :: NDIMS = 3

  ! MIN/MAX INDEX OF INPUT POINTS WHEN PROJECTED ONTO COMPUTATIONAL GRID
  INTEGER(IPP), ALLOCATABLE, DIMENSION(:,:)           :: PS, PE

  ! ALLOW GRID RE-ORDERING
  LOGICAL,                                  PARAMETER :: REORDER = .TRUE.

  ! SET CARTESIAN GRID PERIODICITY
  LOGICAL,                   DIMENSION(3),  PARAMETER :: ISPERIODIC = [.FALSE., .FALSE., .FALSE.]

  REAL(FPP),                                PARAMETER :: BIG = HUGE(1._FPP)

  ! ABSOLUTE POSITION OF MODEL'S FIRST POINT
  REAL(FPP),                 DIMENSION(3)             :: OFF_AXIS

  ! POINTER FOR FOURIER TRANSFORM
  COMPLEX(C_CPP),            DIMENSION(:),  POINTER   :: CDUM => NULL()

  ! POINTER FOR FOURIER TRANSFORM
  REAL(C_FPP),               DIMENSION(:),  POINTER   :: RDUM => NULL()

  ! POINTER TO FFTW PLAN
  TYPE(C_PTR)                                                   :: PC

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE SCARF3D_FFT(X, Y, Z, DH, ACF, CL, SIGMA, HURST, SEED, POI, MUTE, TAPER, FIELD, INFO)

      ! GLOBAL INDICES AND COMMUNICATOR SUBGROUPPING COULD BE HANDLED ELEGANTLY IF VIRTUAL TOPOLOGIES ARE USED IN CALLING PROGRAM.
      ! HOWEVER WE ASSUME THAT THESE ARE NOT USED AND THEREFORE WE ADOPT A SIMPLER APPROACH BASED ON COLLECTIVE CALLS.

      ! "SEED" AND NUMBER OF POINTS "NPTS" SHOULD BE THE SAME TO GUARANTEE SAME RANDOM FIELD

      REAL(FPP),                     DIMENSION(:),     INTENT(IN)  :: X, Y, Z              !< POSITION OF POINTS ALONG X, Y, Z
      REAL(FPP),                                       INTENT(IN)  :: DH                   !< GRID-STEP
      CHARACTER(LEN=*),                                INTENT(IN)  :: ACF                  !< AUTOCORRELATION FUNCTION: "VK" OR "GAUSS"
      REAL(FPP),                     DIMENSION(3),     INTENT(IN)  :: CL                   !< CORRELATION LENGTH
      REAL(FPP),                                       INTENT(IN)  :: SIGMA                !< STANDARD DEVIATION
      REAL(FPP),                                       INTENT(IN)  :: HURST                !< HURST EXPONENT
      INTEGER(IPP),                                    INTENT(IN)  :: SEED                 !< SEED NUMBER
      INTEGER(IPP),                  DIMENSION(:,:),   INTENT(IN)  :: POI                  !< LOCATION OF POINT(S)-OF-INTEREST
      INTEGER(IPP),                                    INTENT(IN)  :: MUTE                 !< NUMBER OF POINTS WHERE MUTING IS APPLIED
      INTEGER(IPP),                                    INTENT(IN)  :: TAPER                !< NUMBER OF POINTS WHERE TAPERING IS APPLIED
      REAL(FPP),                     DIMENSION(:),     INTENT(OUT) :: FIELD                !< VELOCITY FIELD TO BE PERTURBED
      REAL(FPP),                     DIMENSION(8),     INTENT(OUT) :: INFO                 !< ERRORS AND TIMING FOR PERFORMANCE ANALYSIS
      COMPLEX(FPP),     ALLOCATABLE, DIMENSION(:,:,:)              :: SPEC                 !< SPECTRUM/RANDOM FIELD
      INTEGER(IPP)                                                 :: IERR                 !< MPI STUFF
      INTEGER(IPP)                                                 :: I, J, K              !< COUNTERS
      INTEGER(IPP)                                                 :: N                    !< POINTS IN "BUFFER"
      INTEGER(IPP)                                                 :: ITER
      INTEGER(IPP)                                                 :: RANK, NTASKS
      INTEGER(IPP)                                                 :: CARTOPO
      INTEGER(IPP)                                                 :: COMM
      INTEGER(IPP)                                                 :: OFFSET               !< EXTRA POINTS ON EACH SIDE FOR FFT PERIODICITY
      INTEGER(IPP),                  DIMENSION(3)                  :: M                    !< POINTS FOR CALLING PROCESS
      INTEGER(IPP),                  DIMENSION(3)                  :: LS, LE               !< FIRST/LAST INDEX ALONG X, Y, Z
      INTEGER(IPP),                  DIMENSION(3)                  :: COORDS
      INTEGER(IPP),                  DIMENSION(3)                  :: DIMS
      INTEGER(IPP),     ALLOCATABLE, DIMENSION(:)                  :: COMPLETED
      INTEGER(IPP),     ALLOCATABLE, DIMENSION(:)                  :: ISBUSY
      INTEGER(IPP),     ALLOCATABLE, DIMENSION(:)                  :: NPOINTS
      INTEGER(IPP),     ALLOCATABLE, DIMENSION(:)                  :: BUDDIES, BUFFER
      INTEGER(IPP),     ALLOCATABLE, DIMENSION(:)                  :: I0, I1
      REAL(FPP)                                                    :: SCALING              !< SCALING FACTOR
      REAL(REAL64)                                                 :: TICTOC               !< TIMING
      REAL(FPP),                     DIMENSION(2)                  :: ET                   !< DUMMY FOR ELAPSED TIME
      REAL(FPP),                     DIMENSION(3)                  :: MIN_EXTENT           !< LOWER MODEL LIMITS (GLOBAL)
      REAL(FPP),                     DIMENSION(3)                  :: MAX_EXTENT           !< UPPER MODEL LIMITS (GLOBAL)
      REAL(FPP),        ALLOCATABLE, DIMENSION(:)                  :: VAR, MU              !< STATISTICS: VARIANCE AND AVERAGE
      REAL(FPP),        ALLOCATABLE, DIMENSION(:,:,:)              :: DELTA                !< RANDOM PERTURBATIONS ON CARTESIAN GRID

      INTEGER(IPP)                                                 :: FACTOR
      INTEGER(IPP),     ALLOCATABLE, DIMENSION(:)                  :: SENDCOUNTS, RECVCOUNTS

      !-----------------------------------------------------------------------------------------------------------------------------

      ! GET RANK NUMBER
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, WORLD_RANK, IERR)

      ! GET AVAILABLE MPI PROCESSES
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, WORLD_SIZE, IERR)

      ! MODEL LIMITS (PROCESS-WISE) ALONG EACH AXIS
      MIN_EXTENT = [MINVAL(X, DIM = 1), MINVAL(Y, DIM = 1), MINVAL(Z, DIM = 1)]
      MAX_EXTENT = [MAXVAL(X, DIM = 1), MAXVAL(Y, DIM = 1), MAXVAL(Z, DIM = 1)]

      ! (GLOBAL) MODEL LIMITS
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, MIN_EXTENT, 3, REAL_TYPE, MPI_MIN, MPI_COMM_WORLD, IERR)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, MAX_EXTENT, 3, REAL_TYPE, MPI_MAX, MPI_COMM_WORLD, IERR)

      !IF (WORLD_RANK == 0) PRINT*, 'MIN ', MIN_EXTENT, ' MAX ', MAX_EXTENT

      ! CYCLE OVER THE THREE MAIN DIRECTIONS TO DETERMINE THE NECESSARY MODEL SIZE (IN NUMBER OF POINTS) AND THE ABSOLUTE POSITION OF
      ! OF THE FIRST POINT TO COUNTERACT FFT PERIODICITY
      DO I = 1, 3

        ! MINIMUM NUMBER OF POINTS SUCH THAT RANDOM FIELD COVERS WHOLE DOMAIN
        NPTS(I) = NINT( (MAX_EXTENT(I) - MIN_EXTENT(I)) / DH) + 1

        ! POINTS FOR ONE CORRELATION LENGTH
        OFFSET = NINT(CL(I) / DH) + 1

        !OFFSET = 0

        ! EXTEND THE MODEL BY AT LEAST ONE CORRELATION LENGTH (ON EACH SIDE) TO COUNTERACT FFT PERIODICITY
        NPTS(I) = NPTS(I) + 2 * OFFSET

        ! MAKE SURE WE HAVE EVEN NUMBER OF POINTS
        IF (MOD(NPTS(I), 2) .NE. 0) NPTS(I) = NPTS(I) + 1

        ! ABSOLUTE POSITION OF FIRST POINT (IT COULD BE NEGATIVE)
        OFF_AXIS(I) = MIN_EXTENT(I) - OFFSET * DH
        !OFF_AXIS(I) = MIN_EXTENT(I) - (OFFSET - 0.5_FPP) * DH

      ENDDO

      INFO(:) = 0._FPP

      ! HERE WE CHECK IF THE MODEL IS LARGE ENOUGH TO CATCH THE LOWER PART OF THE SPECTRUM (LOW WAVENUMBERS)...
      IF (ANY(NPTS .LE. NINT(2._FPP * PI * CL / DH))) INFO(1) = 1._FPP

      ! ...AND HERE IF THE GRID-STEP IS SMALL ENOUGH TO CATCH THE UPPER PART OF THE SPECTRUM (HIGH WAVENUMBERS)
      IF (ANY(DH .GT. CL / 2._FPP)) INFO(2) = 1._FPP

      ! [ANOTHER WAY TO CHECK HOW CLOSE THE DISCRETE AND CONTINUOUS SPECTRA ARE, IS TO COMPARE THE RESPECTIVE STANDARD DEVIATIONS]

      ALLOCATE(PS(3, 0:WORLD_SIZE - 1), PE(3, 0:WORLD_SIZE - 1))

      ! MIN/MAX INDEX OF INPUT POINTS WHEN PROJECTED ONTO COMPUTATIONAL GRID
      PS(1, WORLD_RANK) = FLOOR((MINVAL(X, DIM = 1) - OFF_AXIS(1)) / DH) + 1
      PE(1, WORLD_RANK) = FLOOR((MAXVAL(X, DIM = 1) - OFF_AXIS(1)) / DH) + 1
      PS(2, WORLD_RANK) = FLOOR((MINVAL(Y, DIM = 1) - OFF_AXIS(2)) / DH) + 1
      PE(2, WORLD_RANK) = FLOOR((MAXVAL(Y, DIM = 1) - OFF_AXIS(2)) / DH) + 1
      PS(3, WORLD_RANK) = FLOOR((MINVAL(Z, DIM = 1) - OFF_AXIS(3)) / DH) + 1
      PE(3, WORLD_RANK) = FLOOR((MAXVAL(Z, DIM = 1) - OFF_AXIS(3)) / DH) + 1

      ! MAKE ALL PROCESSES AWARE OF THESE INDICES
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, PS, 3, MPI_INTEGER, MPI_COMM_WORLD, IERR)
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, PE, 3, MPI_INTEGER, MPI_COMM_WORLD, IERR)

      ALLOCATE(NPOINTS(0:WORLD_SIZE - 1))

      ! NUMBER OF INPUT POINTS
      NPOINTS(WORLD_RANK) = SIZE(X)

      ! SHARE THIS INFORMATION AMONGST PROCESSES
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, NPOINTS, 1, MPI_INTEGER, MPI_COMM_WORLD, IERR)


      ! PRINT*, WORLD_RANK, 'POINTS: ', NPOINTS(WORLD_RANK), ', PS/PE= ', PS(:,WORLD_RANK), PE(:,WORLD_RANK)


      ! HERE BELOW WE CREATE A REGULAR MESH AND A CARTESIAN TOPOLOGY. THE RANDOM FIELD WILL BE CALCULATED ON THIS MESH OF "NPTS" POINTS
      ! AND THEN INTERPOLATED AT THOSE POINTS (POSSIBLY IRREGULARLY DISTRIBUTED) OWNED BY EACH SINGLE PROCESS

      ! "GS"/"GE" STORE FIRST/LAST GLOBAL INDICES ALONG EACH DIRECTION FOR ALL PROCESS
      ALLOCATE(GS(3, 0:WORLD_SIZE-1), GE(3, 0:WORLD_SIZE-1))

      CALL WATCH_START(TICTOC)

      ! RETURN PROCESSORS GRID RESULTING MAXIMIZING THE NUMBER OF OVERLAPPING CALLS FOR INTERPOLATION (SEE BELOW)
      !CALL BEST_CONFIG(DIMS, ITER)
      CALL BEST_CONFIG(DIMS)

      CALL WATCH_STOP(TICTOC)

      IF (WORLD_RANK == 0) PRINT*, 'BEST DIMS ', DIMS, ' -- ', ITER, ' -- ', REAL(TICTOC, REAL32)

      ! CREATE TOPOLOGY
      CALL MPI_CART_CREATE(MPI_COMM_WORLD, NDIMS, DIMS, ISPERIODIC, REORDER, CARTOPO, IERR)

      ! RETURN PROCESS COORDINATES IN CURRENT TOPOLOGY
      CALL MPI_CART_COORDS(CARTOPO, WORLD_RANK, NDIMS, COORDS, IERR)

      ! RETURN FIRST/LAST-INDEX ("LS"/"LE") ALONG EACH DIRECTION FOR CALLING PROCESS. NOTE: FIRST POINT HAS ALWAYS INDEX EQUAL TO 1.
      CALL COORDS2INDEX(NPTS, DIMS, COORDS, LS, LE)

      ! DO I = 0, WORLD_SIZE - 1
      !   IF (WORLD_RANK .EQ. I) THEN
      !     PRINT*, 'FFT MESH: RANK', WORLD_RANK, '= ', LS(1), LE(1), ' -- ', LS(2), LE(2), ' -- ', LS(3), LE(3)
      !   ENDIF
      !   CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
      ! ENDDO
      !
      ! IF (WORLD_RANK == 0) PRINT*, 'NPTS: ', NPTS, ', OFF_AXIS= ', OFF_AXIS

      GS(:, WORLD_RANK) = LS
      GE(:, WORLD_RANK) = LE

      ! MAKE ALL PROCESSES AWARE OF GLOBAL INDICES ALONG EACH AXIS
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, GS, 3, MPI_INTEGER, MPI_COMM_WORLD, IERR)
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, GE, 3, MPI_INTEGER, MPI_COMM_WORLD, IERR)

      ! NUMBER OF POINTS FOR CALLING PROCESS
      DO I = 1, 3
        M(I) = LE(I) - LS(I) + 1
      ENDDO

      ! ALLOCATE MEMORY FOR SPECTRUM
      ALLOCATE(SPEC(M(1), M(2), M(3)))

      ! COMPUTE SPECTRUM AND APPLY HERMITIAN SYMMETRY
      CALL COMPUTE_SPECTRUM(LS, LE, DH, ACF, CL, SIGMA, HURST, SEED, SPEC, ET)

      INFO(3:4) = ET

      ! START TIMER
      CALL WATCH_START(TICTOC)

      ! TRANSFORM ALONG EACH DIRECTION
      CALL TRANSFORM_ALONG_Z(SPEC)
      CALL TRANSFORM_ALONG_Y(SPEC)
      CALL TRANSFORM_ALONG_X(SPEC)

      CALL WATCH_STOP(TICTOC)

      INFO(5) = TICTOC

      ! SCALING PARAMETER
      SCALING = 1._FPP / SQRT(REAL(NPTS(1), FPP) * REAL(NPTS(2), FPP) * REAL(NPTS(3), FPP) * DH**3)

      ALLOCATE(DELTA(0:M(1), 0:M(2), 0:M(3)))

      DELTA = 0._FPP

      ! SCALE IFFT
      DO K = 1, M(3)
        DO J = 1, M(2)
          DO I = 1, M(1)
            DELTA(I, J, K) = REAL(SPEC(I, J, K), FPP) * SCALING
          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE(SPEC)


      CALL WATCH_START(TICTOC)

      ! EXCHANGE HALO
      CALL EXCHANGE_HALO(CARTOPO, DELTA)

      ALLOCATE(SENDCOUNTS(0:WORLD_SIZE - 1), RECVCOUNTS(0:WORLD_SIZE - 1))

      CALL POINTS_DISTRIBUTION(X, Y, Z, SENDCOUNTS, RECVCOUNTS)

      FACTOR = AUTO_TUNING(M, SENDCOUNTS, RECVCOUNTS)





      CALL WATCH_STOP(TICTOC)

      INFO(6) = TICTOC

      ! RELEASE MEMORY
      DEALLOCATE(COMPLETED, ISBUSY)
      DEALLOCATE(DELTA)
      DEALLOCATE(GS, GE)

      CALL MPI_COMM_FREE(CARTOPO, IERR)

      !CALL IO_WRITE_ONE(REAL(SPEC, FPP), 'random_struct.bin')
      !CALL IO_WRITE_SLICE(1, 50, REAL(SPEC, FPP), 'random_struct_slice.bin')

      ALLOCATE(VAR(0:WORLD_SIZE - 1), MU(0:WORLD_SIZE - 1))

      ! COMPUTE VARIANCE AND MEAN OF RANDOM FIELD FOR EACH SINGLE PROCESS
      VAR(WORLD_RANK) = VARIANCE(FIELD)
      MU(WORLD_RANK)  = MEAN(FIELD)

      ! SHARE RESULTS
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, VAR, 1, REAL_TYPE, MPI_COMM_WORLD, IERR)
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, MU, 1, REAL_TYPE, MPI_COMM_WORLD, IERR)

      ! COMPUTE TOTAL VARIANCE AND MEAN
      CALL PARALLEL_VARIANCE(VAR, MU, NPOINTS, INFO(7), INFO(8))

      ! RETURN STANDARD DEVIATION
      INFO(7) = SQRT(INFO(7))

      ! SINCE WE HAVE SELECTED A SUBSET OF THE COMPUTED RANDOM FIELD (E.G. SEE "OFFSET"), MEAN MAY DEVIATE SOMEWHAT FROM ZERO. FOR
      ! THIS REASON, NOW REMOVE MEAN FROM OUTPUT RANDOM FIELD TO GET AN EFFECTIVELY ZERO-MEAN
      FIELD = FIELD - INFO(8)

      DEALLOCATE(VAR, MU, NPOINTS)

      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

    END SUBROUTINE SCARF3D_FFT

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE SHARE_AND_INTERPOLATE(DH, DELTA, X, Y, Z, FIELD)

      REAL(FPP),                          INTENT(IN)    :: DH                             !< FFT GRID-STEP
      REAL(FPP),    DIMENSION(:,:,:),     INTENT(IN)    :: DELTA                          !< RANDOM PERTURBATIONS ON FFT GRID
      REAL(FPP),    DIMENSION(:),         INTENT(IN)    :: X, Y, Z                        !< POINTS LOCATION
      REAL(FPP),    DIMENSION(:),         INTENT(INOUT) :: FIELD


      INTEGER(IPP)                                         :: I, P, N
      INTEGER(IPP)                                         :: NP, BLOCKWIDTH, NBLOCKS
      INTEGER(IPP),              DIMENSION(0:WORLD_SIZE-1) :: SCOUNTS, RCOUNTS
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)              :: P0, P1, MAP
      REAL(FPP),                 DIMENSION(6)              :: INFO
      REAL(FPP),    ALLOCATABLE, DIMENSION(:)              :: SV, RV
      REAL(FPP),    ALLOCATABLE, DIMENSION(:,:)            :: SXYZ, RXYZ

      !-----------------------------------------------------------------------------------------------------------------------------

      INFO = 0._FPP

      NP = SIZE(X)

      CALL MPI_ALLREDUCE(MPI_IN_PLACE, NP, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, IERR)

      ! FIND BEST NUMBER OF POINTS TO BE SENT/RECEIVED
      BLOCKWIDTH = 32 * 1

      NBLOCKS = NP / BLOCKWIDTH

      IF (MOD(NP, BLOCKWIDTH) .NE. 0) NBLOCKS = NBLOCKS + 1

      ALLOCATE(P0(NBLOCKS), P1(NBLOCKS))

      DO I = 1, NBLOCKS - 1
        P0(I) = (I - 1) * BLOCKWIDTH + 1
        P1(I) = I * BLOCKWIDTH
      ENDDO

      P0(I) = P1(I - 1) + 1
      P1(I) = NP


      DO P = 1, NBLOCKS

        N = P1(P) - P0(P) + 1

        ALLOCATE(SXYZ(3, N))

        CALL WATCH_START(TICTOC)
        CALL ORDER_POINTS(P0(P), P1(P), DH, X, Y, Z, SXYZ, SCOUNTS, MAP)
        CALL WATCH_STOP(TICTOC); INFO(1) = INFO(1) + TICTOC

        CALL WATCH_START(TICTOC)
        CALL MPI_ALLTOALL(SCOUNTS, 1, MPI_INTEGER, RCOUNTS, 1, MPI_INTEGER, MPI_COMM_WORLD, IERR)
        CALL WATCH_STOP(TICTOC); INFO(2) = INFO(2) + TICTOC

        RDISPLS(1) = 0

        DO I = 1, WORLD_SIZE - 1
          RDISPLS(I) = RDISPLS(I - 1) + RCOUNTS(I - 1)
        ENDDO

        N = SUM(RCOUNTS)

        ALLOCATE(RXYZ(3, N))

        CALL WATCH_START(TICTOC)
        CALL MPI_ALLTOALLV(SXYZ, SCOUNTS, SDISPLS, REAL_TYPE, RXYZ, RCOUNTS, RDISPLS, REAL_TYPE, MPI_COMM_WORLD, IERR)
        CALL WATCH_STOP(TICTOC); INFO(3) = INFO(3) + TICTOC

        ALLOCATE(SV(N))

        CALL WATCH_START(TICTOC)
        CALL INTERPOLATE(RXYZ, DELTA, SV)
        CALL WATCH_STOP(TICTOC); INFO(4) = INFO(4) + TICTOC

        ALLOCATE(RV(SIZE(MAP)))

        CALL WATCH_START(TICTOC)
        CALL MPI_ALLTOALLV(SV, RCOUNTS, RDISPLS, REAL_TYPE, RV, SCOUNTS, SDISPLS, REAL_TYPE, MPI_COMM_WORLD, IERR)
        CALL WATCH_STOP(TICTOC); INFO(5) = INFO(5) + TICTOC

        CALL WATCH_START(TICTOC)
        DO I = 1, SIZE(MAP)
          FIELD(MAP(I)) = RV(I)
        ENDDO
        CALL WATCH_STOP(TICTOC); INFO(6) = INFO(6) + TICTOC

        DEALLOCATE(RXYZ, SV, RV, SXYZ, MAP)

      ENDDO

      ! GET FIRST BATCH OF POINTS
      ! CALL ORDER_POINTS(P0(1), P1(1), X, Y, Z, SXYZ1, SCOUNTS1, SDISPLS1)
      !
      ! ! RECVCOUNTS / RDISPLS FOR FIRST BATCH
      ! CALL MPI_ALLTOALL(SCOUNTS1, 1, MPI_INTEGER, RCOUNTS1, 1, MPI_INTEGER, MPI_COMM_WORLD, IERR)
      ! CALL MPI_ALLTOALL(SDISPLS1, 1, MPI_INTEGER, RDISPLS1, 1, MPI_INTEGER, MPI_COMM_WORLD, IERR)
      !
      ! ! NON-BLOCKING FORWARD, FIRST BATCH
      ! CALL MPI_IALLTOALLV(SXYZ1, SCOUNTS1, SDISPLS1, REAL_TYPE, RXYZ1, RCOUNTS1, RDISPLS1, REAL_TYPE, MPI_COMM_WORLD, REQ1, IERR)
      !
      ! DO P = 2, NBLOCKS - 1
      !
      !   ! GET SECOND BATCH OF POINTS
      !   CALL ORDER_POINTS(P0(P), P1(P), X, Y, Z, SXYZ2, MAP2, SCOUNTS2, SDISPLS2)
      !
      !   ! NON-BLOCKING RECVCOUNTS / RDISPLS FOR SECOND BATCH
      !   CALL MPI_IALLTOALL(SCOUNTS2, 1, MPI_INTEGER, RCOUNTS2, 1, MPI_INTEGER, MPI_COMM_WORLD, REQ21, IERR)
      !   CALL MPI_IALLTOALL(SDISPLS2, 1, MPI_INTEGER, RDISPLS2, 1, MPI_INTEGER, MPI_COMM_WORLD, REQ22, IERR)
      !
      !   ! WAIT FOR NON-BLOCKING FORWARD, FIRST BATCH
      !   CALL MPI_WAIT(REQ1, STATUS, IERR)
      !
      !   ! INTERPOLATE FIRST BATCH
      !   CALL INTERPOLATE(RXYZ1, DELTA, SV1)
      !
      !   ! NON-BLOCKING BACKWARD, FIRST BATCH
      !   CALL MPI_IALLTOALLV(SV1, RCOUNTS1, RDISPLS1, REAL_TYPE, RV1, SCOUNTS1, SDISPLS1, REAL_TYPE, MPI_COMM_WORLD, REQ1, IERR)
      !
      !   ! WAIT FOR NON-BLOCKING RECVCOUNTS / RDISPLS, SECOND BATCH
      !   CALL MPI_WAIT(REQ21, STATUS, IERR)
      !   CALL MPI_WAIT(REQ22, STATUS, IERR)
      !
      !   ! NON-BLOCKING FORWARD, SECOND BATCH
      !   CALL MPI_IALLTOALLV(SXYZ2, SCOUNTS2, SDISPLS2, REAL_TYPE, RXYZ2, RCOUNTS2, RDISPLS2, REAL_TYPE, MPI_COMM_WORLD, REQ2, IERR)
      !
      !   ! WAIT FOR NON-BLOCKING BACKWARD, FIRST BATCH
      !   CALL MPI_WAIT(REQ1, STATUS, IERR)
      !
      !   ! COPY BACK FIRST BATCH
      !   CALL COPY_BACK(RV1, MAP1, FIELD)
      !
      !   ! WAIT FOR NON-BLOCKING FORWARD, SECOND BATCH
      !   CALL MPI_WAIT(REQ2, STATUS, IERR)
      !
      !   ! INTERPOLATE SECOND BATCH: THIS IS NOT OVERLAPPED BY ANY MESSAGING. EXPOSED!
      !   CALL INTERPOLATE(RXYZ2, DELTA, SV2)
      !
      !   ! NON-BLOCKING BACKWARD, SECOND BATCH
      !   CALL MPI_IALLTOALLV(SV2, RCOUNTS2, RDISPLS2, REAL_TYPE, RV2, SCOUNTS2, SDISPLS2, REAL_TYPE, MPI_COMM_WORLD, REQ2, IERR)
      !
      !   ! GET FIRST BATCH OF POINTS
      !   CALL ORDER_POINTS(P0(P), P1(P), X, Y, Z, SXYZ1, MAP1, SCOUNTS1, SDISPLS1)
      !
      !   ! NON-BLOCKING RECVCOUNTS / RDISPLS FOR FIRST BATCH
      !   CALL MPI_IALLTOALL(SCOUNTS1, 1, MPI_INTEGER, RCOUNTS1, 1, MPI_INTEGER, MPI_COMM_WORLD, REQ11, IERR)
      !   CALL MPI_IALLTOALL(SDISPLS1, 1, MPI_INTEGER, RDISPLS1, 1, MPI_INTEGER, MPI_COMM_WORLD, REQ12, IERR)
      !
      !   ! WAIT FOR NON-BLOCKING BACKWARD, SECOND BATCH
      !   CALL WAIT(REQ2, STATUS, IERR)
      !
      !   ! COPY BACK SECOND BATCH
      !   CALL COPY_BACK(RV2, MAP2, FIELD)
      !
      !   ! WAIT FOR NON-BLOCKING RECVCOUNTS / RDISPLS, FIRST BATCH
      !   CALL WAIT(REQ11, STATUS, IERR)
      !   CALL WAIT(REQ12, STATUS, IERR)
      !
      !   ! NON-BLOCKING FORWARD, FIRST BATCH
      !   CALL MPI_IALLTOALLV(SXYZ1, SCOUNTS1, SDISPLS1, REAL_TYPE, RXYZ1, RCOUNTS1, RDISPLS1, REAL_TYPE, MPI_COMM_WORLD, REQ1, IERR)
      !
      ! ENDDO

      ! WAIT FOR NON-BLOCKING FORWARD, FIRST BATCH
      CALL MPI_WAIT(REQ1, STATUS, IERR)

      ! INTERPOLATE FIRST BATCH
      CALL INTERPOLATE(RXYZ1, DELTA, SV1)

      ! BLOCKING BACKWARD, FIRST BATCH
      CALL MPI_ALLTOALLV(SV, RCOUNTS, RDISPLS, REAL_TYPE, RV, SCOUNTS, SDISPLS, REAL_TYPE, MPI_COMM_WORLD, IERR)

      ! COPY BACK FIRST BATCH
      CALL COPY_BACK(RV, MAP, FIELD)

    END SUBROUTINE SHARE_AND_INTERPOLATE


    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE ORDER_POINTS(P0, P1, DH, X, Y, Z, XYZ, SENDCOUNTS, MAP)

      INTEGER(IPP),                 INTENT(IN)  :: P0, P1
      REAL(FPP),                    INTENT(IN)  :: DH
      REAL(FPP),    DIMENSION(:),   INTENT(IN)  :: X, Y, Z
      REAL(FPP),    DIMENSION(:,:), INTENT(OUT) :: XYZ
      INTEGER(IPP), DIMENSION(:),   INTENT(OUT) :: SENDCOUNTS
      INTEGER(IPP), DIMENSION(:),   INTENT(OUT) :: MAP

      INTEGER(IPP)                              :: I, J, N
      LOGICAL                                   :: BOOL

      !-----------------------------------------------------------------------------------------------------------------------------

      ! FIND NUMBER OF POINTS THAT CALLING PROCESS MUST SEND TO THE OTHER PROCESSES
      DO I = 0, WORLD_SIZE - 1

        X0 = (GS(1, I) - 2) * DH - OFF_AXIS(1)
        X1 = (GE(1, I) - 1) * DH - OFF_AXIS(1)

        Y0 = (GS(2, I) - 2) * DH - OFF_AXIS(2)
        Y1 = (GE(2, I) - 1) * DH - OFF_AXIS(2)

        Z0 = (GS(3, I) - 2) * DH - OFF_AXIS(3)
        Z1 = (GE(3, I) - 1) * DH - OFF_AXIS(3)

        N = 0

        DO J = P0, P1

          BOOL = (X(J) .GE. X0) .AND. (X(J) .LT. X1) .AND.   &
                 (Y(J) .GE. Y0) .AND. (Y(J) .LT. Y1) .AND.   &
                 (Z(J) .GE. Z0) .AND. (Z(J) .LT. Z1)

          IF (BOOL .EQV. .TRUE.) THEN
            N         = N + 1
            XYZ(:, N) = []
            MAP(N)    = I
          ENDIF

        ENDDO

        SENDCOUNTS(I + 1) = N

      ENDDO

    END SUBROUTINE ORDER_POINTS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE COPY_BACK(BUFFER, MAP, FIELD)

      REAL(FPP),    DIMENSION(:), INTENT(IN)    :: BUFFER
      INTEGER(IPP), DIMENSION(:), INTENT(IN)    :: MAP
      REAL(FPP),    DIMENSION(:), INTENT(INOUT) :: FIELD

      !-----------------------------------------------------------------------------------------------------------------------------

      DO I = 1, SIZE(BUFFER)
        FIELD(MAP(I)) = BUFFER(I)
      ENDDO

    END SUBROUTINE COPY_BACK

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    ! SUBROUTINE SHARE_AND_INTERPOLATE(COMM, RANK, NTASKS, I0, I1, DH, DELTA, X, Y, Z, FIELD)
    !
    !   ! SHARE A SUBSET OF INPUT POINTS BELONGING TO MASTER PROCESS WITH SLAVE PROCESSES, INTERPOLATE FROM FFT GRID ONTO POINTS LOCATION,
    !   ! SEND EVERYTHING BACK TO MASTER WHICH WILL COPY ONLY THOSE POINTS THAT WERE INTERPOLATED.
    !
    !   INTEGER(IPP),                       INTENT(IN)    :: COMM                           !< COMMUNICATOR
    !   INTEGER(IPP),                       INTENT(IN)    :: RANK                           !< RANK OF CALLING PROCESS IN COMMUNICATOR
    !   INTEGER(IPP),                       INTENT(IN)    :: NTASKS                         !< SIZE OF COMMUNICATOR
    !   INTEGER(IPP),                       INTENT(IN)    :: I0, I1                         !< FIRST/LAST POINT INDEX
    !   REAL(FPP),                          INTENT(IN)    :: DH                             !< FFT GRID-STEP
    !   REAL(FPP),    DIMENSION(:,:,:),     INTENT(IN)    :: DELTA                          !< RANDOM PERTURBATIONS ON FFT GRID
    !   REAL(FPP),    DIMENSION(:),         INTENT(IN)    :: X, Y, Z                        !< POINTS LOCATION
    !   REAL(FPP),    DIMENSION(:),         INTENT(INOUT) :: FIELD                          !< INTERPOLATED RANDOM FIELD
    !   INTEGER(IPP)                                      :: I, N                           !< COUNTERS
    !   INTEGER(IPP)                                      :: IERR                           !< MPI STUFF
    !   INTEGER(IPP), DIMENSION(0:NTASKS-1)               :: RECVCOUNTS, DISPLS             !< MPI STUFF
    !   REAL(FPP),    DIMENSION(I0:I1)                    :: XP, YP, ZP, VP                 !< TEMPORARY ARRAYS
    !
    !   !-----------------------------------------------------------------------------------------------------------------------------
    !
    !   N = I1 - I0 + 1
    !
    !   ! COPY INPUT DATA INTO TEMPORARY ARRAYS OVERWRITTEN DURING BCASTING ("IF" IS MEANT ONLY TO AVOID THE POSSIBILITY THAT SLAVES
    !   ! HAVE MANY LESS INPUT POINTS THAN MASTER PROCESS)
    !   IF (RANK .EQ. 0) THEN
    !     DO I = I0, I1
    !       XP(I) = X(I)
    !       YP(I) = Y(I)
    !       ZP(I) = Z(I)
    !     ENDDO
    !   ENDIF
    !
    !   ! BROADCAST DATA TO SLAVES
    !   CALL MPI_BCAST(XP, N, REAL_TYPE, 0, COMM, IERR)
    !   CALL MPI_BCAST(YP, N, REAL_TYPE, 0, COMM, IERR)
    !   CALL MPI_BCAST(ZP, N, REAL_TYPE, 0, COMM, IERR)
    !
    !   ! MARK ALL POINTS WITH A VERY LARGE NUMBER. POINTS NOT INTERPOLATED WILL BEAR THIS VALUE.
    !   VP = BIG
    !
    !   CALL INTERPOLATE(DH, XP, YP, ZP, DELTA, VP)
    !
    !   ! SET INITIAL VALUE FOR COUNTER
    !   N = I0 - 1
    !
    !   ! RE-USE TEMPORARY ARRAYS BUT RENAME VARIABLES FOR CLARITY, KEEP ORIGINAL BOUNDS
    !   ASSOCIATE( VAL => XP, POS => YP )
    !
    !   ! MANUAL PACKING OF INTERPOLATED POINTS ONLY
    !   DO I = I0, I1
    !     IF (VP(I) .NE. BIG) THEN
    !       N      = N + 1
    !       VAL(N) = VP(I)
    !       POS(N) = REAL(I, FPP)
    !     ENDIF
    !   ENDDO
    !
    !   ! GET EFFECTIVE NUMBER OF INTERPOLATED POINTS
    !   N = N - I0 + 1
    !
    !   ! STORE NUMBER OF INTERPOLATED VALUES FOR EACH PROCESS
    !   RECVCOUNTS(RANK) = N
    !
    !   ! MAKE ALL PROCESSES AWARE
    !   CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, RECVCOUNTS, 1, MPI_INTEGER, COMM, IERR)
    !
    !   ! SET DISPLACEMENTS
    !   DISPLS(0) = 0
    !
    !   DO I = 1, NTASKS - 1
    !     DISPLS(I) = DISPLS(I - 1) + RECVCOUNTS(I - 1)
    !   ENDDO
    !
    !   ! COLLECT INTERPOLATED VALUES AND INDICES
    !   CALL MPI_GATHERV(POS, N, REAL_TYPE, POS, RECVCOUNTS, DISPLS, REAL_TYPE, 0, COMM, IERR)
    !   CALL MPI_GATHERV(VAL, N, REAL_TYPE, VAL, RECVCOUNTS, DISPLS, REAL_TYPE, 0, COMM, IERR)
    !
    !   ! ONLY MASTER PROCESS COPY INTERPOLATED VALUES (AT THEIR RIGHT LOCATION)
    !   IF (RANK .EQ. 0) THEN
    !     DO I = I0, I1
    !       FIELD(NINT(POS(I))) = VAL(I)
    !     ENDDO
    !   ENDIF
    !
    !   END ASSOCIATE
    !
    ! END SUBROUTINE SHARE_AND_INTERPOLATE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE FIND_BUDDIES(RANK, BUDDIES)

      ! MAKE A LIST OF PROCESSES WHOSE FFT GRID INCLUDE CURRENT SELECTION OF INPUT POINTS. NOTE THAT "GS" DO NOT INCLUDE HALOS! THESE
      ! CAN BE INCLUDED BY USING "GS - 1"

      INTEGER(IPP),                            INTENT(IN)  :: RANK                   !< GLOBAL PROCESS
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: BUDDIES                !< LIST OF PROCESSES
      INTEGER(IPP)                                         :: I, J                   !< COUNTERS
      LOGICAL,                   DIMENSION(3)              :: BOOL                   !< DETERMINE IF N-TH PROCESS MUST BE INCLUDED

      !-----------------------------------------------------------------------------------------------------------------------------

      ! FOR SURE WE NEED TO INCLUDE PROCESS OWNING POINTS
      BUDDIES = [RANK]

      ! LOOP OVER ALL PROCESSES
      DO I = 0, WORLD_SIZE - 1

        ! SKIP IF WE STUMBLED OVER CALLING PROCESS
        IF (I .EQ. RANK) CYCLE

        DO J = 1, 3
          BOOL(J) = (PS(J, RANK) .LT. GE(J, I)) .AND. (PS(J, RANK) .GE. (GS(J, I) - 1)) .OR.    &
                    (PE(J, RANK) .LT. GE(J, I)) .AND. (PE(J, RANK) .GE. (GS(J, I) - 1)) .OR.    &
                    (PE(J, RANK) .GE. GE(J, I)) .AND. (PS(J, RANK) .LT. (GS(J, I) - 1))
        ENDDO

        ! ADD I-TH PROCESS (USE AUTOMATIC REALLOCATION)
        IF (ALL(BOOL .EQV. .TRUE.)) BUDDIES = [BUDDIES, I]

      ENDDO

    END SUBROUTINE FIND_BUDDIES

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE MAKE_COMMUNICATOR(LIST, BUDDIES, COMM, RANK, NTASKS)

      ! BUILD A NEW COMMUNICATOR BASED ON PROCESSES CONTAINED IN "LIST". "BUDDIES", "COMM", "RANK" AND "NTASKS" ARE MODIFIED ONLY FOR
      ! THOSE PROCESSES IN THE NEW COMMUNICATOR. "LIST" IS INDENTICAL FOR ALL CALLING PROCESSES.

      INTEGER(IPP),              DIMENSION(:),             INTENT(IN)    :: LIST                         !< PROCESSES IN NEW COMMUNICATOR
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:),             INTENT(INOUT) :: BUDDIES                      !< PROCESSES IN NEW COMMUNICATOR
      INTEGER(IPP),                                        INTENT(INOUT) :: COMM                         !< NEW COMMUNICATOR
      INTEGER(IPP),                                        INTENT(INOUT) :: RANK                         !< RANK OF CALLING PROCESS IN NEW COMMUNICATOR
      INTEGER(IPP),                                        INTENT(INOUT) :: NTASKS                       !< NEW COMMUNICATOR SIZE
      INTEGER(IPP)                                                       :: I                            !< COUNTER
      INTEGER(IPP)                                                       :: IERR, KEY, DUM               !< MPI STUFF
      INTEGER(IPP),              DIMENSION(0:WORLD_SIZE-1)               :: COLOR                        !< MPI STUFF

      !-----------------------------------------------------------------------------------------------------------------------------

      COLOR(:) = MPI_UNDEFINED

      ! CHANGE COLOR ONLY TO THOSE PROCESSES IN "LIST"
      DO I = 1, SIZE(LIST)
        COLOR(LIST(I)) = 1
      ENDDO

      ! DEFAULT KEY FOR COMMUNICATOR-SPLITTING IS PROCESS' OLD RANK
      KEY = WORLD_RANK

      ! MAKE SURE THAT FIRST ELEMENT IN "BUDDY" GETS "RANK=0" WITHIN NEW COMMUNICATOR
      IF (WORLD_RANK .EQ. LIST(1)) KEY = -WORLD_RANK

      ! CREATE COMMUNICATOR SUBGROUP
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, COLOR(WORLD_RANK), KEY, DUM, IERR)

      ! UPDATE "BUDDIES", "COMM", "RANK" AND "NTASKS" ONLY FOR THOSE PROCESSES BEING PART OF THE NEWLY CREATED COMMUNICATOR
      IF (DUM .NE. MPI_COMM_NULL) THEN

        BUDDIES = LIST
        COMM    = DUM

        NTASKS = SIZE(BUDDIES)

        CALL MPI_COMM_RANK(COMM, RANK, IERR)
        !CALL MPI_COMM_SIZE(COMM, NTASKS, IERR)

      ENDIF

    END SUBROUTINE MAKE_COMMUNICATOR

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE TRANSFORM_ALONG_Z(SPECTRUM)

      ! COMPUTE MANY IN-PLACE COMPLEX-TO-COMPLEX IFFT ALONG Z-DIRECTION ON INPUT SPECTRUM. PROCESSES ARE GATHERED IN Z-ORIENTED PENCILS
      ! BASED ON THE INPUT GEOMETRY. INSIDE EACH PENCIL, THE ALGORITHM WORKS ON Y-SLICES (TO REQUIRE AS LESS MEMORY AS POSSIBLE) TO
      ! EXCHANGE DATA AND COMPUTE THE IFFT.

      COMPLEX(FPP),              DIMENSION(:,:,:), INTENT(INOUT) :: SPECTRUM                   !< INPUT FIELD
      COMPLEX(FPP), ALLOCATABLE, DIMENSION(:)                    :: SENDBUF, RECVBUF           !< SENDER/RECEIVER BUFFER
      INTEGER(IPP)                                               :: I, J, K, P, C              !< COUNTERS
      INTEGER(IPP)                                               :: RANK, NTASKS, IERR         !< MPI STUFF
      INTEGER(IPP)                                               :: PENCIL
      INTEGER(IPP)                                               :: NX, NY, NZ                 !< TOTAL NUMBER OF POINTS FOR PENCIL ALONG EACH AXIS
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: I0, I1                     !< LOCAL I-INDEX FOR PROCESSES BELONGING TO PENCIL
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: K0, K1                     !< LOCAL K-INDEX FOR PROCESSES BELONGING TO PENCIL
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: LX, LZ                     !< LOCAL NUMBER OF POINTS ALONG X AND Z AXIS
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: SENDCOUNTS, RECVCOUNTS     !< MPI STUFF
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: SDISPLS, RDISPLS           !< MPI STUFF
      TYPE(C_PTR)                                                :: P1

      !-------------------------------------------------------------------------------------------------------------------------------

      ! GROUP PROCESSES IN PENCILS ORIENTED ALONG Z-AXIS. "LZ" CONTAINS NUMBER OF POINTS FOR EACH PROCESS IN "PENCIL" ALONG Z
      CALL BUILD_PENCIL(2, PENCIL, RANK, NTASKS, LZ)

      ! THIS ARRAY WILL CONTAIN THE NUMBER OF POINTS LOCAL TO EACH PROCESS IN CURRENT PENCIL ALONG X
      ALLOCATE(LX(0:NTASKS-1))

      ! NUMBER OF POINTS ALONG X AND Y. THESE ARE INDENTICAL FOR ALL PROCESSES IN CURRENT PENCIL.
      NX = SIZE(SPECTRUM, 1)
      NY = SIZE(SPECTRUM, 2)

      ! TOTAL NUMBER OF POINTS ALONG Z-AXIS
      NZ = SUM(LZ)

      ALLOCATE(K0(0:NTASKS-1), K1(0:NTASKS-1))

      ! FIRST/LAST K-INDEX FOR EACH PROCESS
      DO P = 0, NTASKS - 1
        IF (P .EQ. 0) THEN
          K0(P) = 1
          K1(P) = LZ(P)
        ELSE
          K0(P) = K1(P - 1) + 1
          K1(P) = K1(P - 1) + LZ(P)
        ENDIF
      ENDDO

      ALLOCATE(I0(0:NTASKS-1), I1(0:NTASKS-1))

      ! DISTRIBUTE POINTS ALONG X-AXIS BETWEEN PROCESSES
      CALL MPI_SPLIT_TASK(NX, NTASKS, I0, I1)

      ! NUMBER OF POINTS ALONG X-AXIS FOR EACH PROCESS
      DO P = 0, NTASKS - 1
        LX(P) = I1(P) - I0(P) + 1
      ENDDO

      ALLOCATE(SENDCOUNTS(0:NTASKS-1), SDISPLS(0:NTASKS-1), RECVCOUNTS(0:NTASKS-1), RDISPLS(0:NTASKS-1))

      ! SET COUNTS/DISPLACEMENT ARRAYS
      DO P = 0, NTASKS - 1

        SENDCOUNTS(P) = LZ(RANK) * LX(P)              !< POINTS SEND FROM PROCESS "RANK" TO PROCESS "P"
        RECVCOUNTS(P) = LZ(P) * LX(RANK)              !< POINTS RECEIVED FROM PROCESS "RANK" BY PROCESS "P"

        IF (P .EQ. 0) THEN
          SDISPLS(P) = 0
          RDISPLS(P) = 0
        ELSE
          SDISPLS(P) = SDISPLS(P - 1) + SENDCOUNTS(P - 1)
          RDISPLS(P) = RDISPLS(P - 1) + RECVCOUNTS(P - 1)
        ENDIF

      ENDDO

      ALLOCATE(SENDBUF(NX * LZ(RANK)))
      ALLOCATE(RECVBUF(NZ * LX(RANK)))

      ! BUFFER "CDUM" WILL BE USED FOR TRANSFORM AND IT IS ALLOCATED USING FFTW MEMORY ALLOCATION UTILITY
#ifdef DOUBLE_PREC
      PC = FFTW_ALLOC_COMPLEX(INT(LX(RANK) * NZ, C_SIZE_T))
#else
      PC = FFTWF_ALLOC_COMPLEX(INT(LX(RANK) * NZ, C_SIZE_T))
#endif

      CALL C_F_POINTER(PC, CDUM, [LX(RANK) * NZ])

      ! PREPARE FFTW PLAN ALONG Z-AXIS ("LX(RANK)" TRANSFORMS, EACH HAVING "NZ" NUMBER OF POINTS)
#ifdef DOUBLE_PREC
      P1 = FFTW_PLAN_MANY_DFT(1, [NZ], LX(RANK), CDUM, [NZ], LX(RANK), 1, CDUM, [NZ], LX(RANK), 1, FFTW_BACKWARD, FFTW_ESTIMATE)
#else
      P1 = FFTWF_PLAN_MANY_DFT(1, [NZ], LX(RANK), CDUM, [NZ], LX(RANK), 1, CDUM, [NZ], LX(RANK), 1, FFTW_BACKWARD, FFTW_ESTIMATE)
#endif

      ! WORK Y-SLICES ONE BY ONE
      DO J = 1, NY

        C = 0

        ! REARRANGE ELEMENTS OF "SPECTRUM"
        DO P = 0, NTASKS - 1
          DO K = 1, LZ(RANK)
            DO I = I0(P), I1(P)
              C          = C + 1
              SENDBUF(C) = SPECTRUM(I, J, K)
            ENDDO
          ENDDO
        ENDDO

        ! EXCHANGE DATA
        CALL MPI_ALLTOALLV(SENDBUF, SENDCOUNTS, SDISPLS, COMPLEX_TYPE, CDUM, RECVCOUNTS, RDISPLS, COMPLEX_TYPE, PENCIL, IERR)

        ! IFFT
#ifdef DOUBLE_PREC
        CALL FFTW_EXECUTE_DFT(P1, CDUM, CDUM)
#else
        CALL FFTWF_EXECUTE_DFT(P1, CDUM, CDUM)
#endif

        C = 0

        ! REARRANGE ELEMENTS OF "CDUM" INTO "RECVBUF"
        DO P = 0, NTASKS - 1
          DO I = 1, LX(RANK)
            DO K = K0(P), K1(P)
              C          = C + 1
              RECVBUF(C) = CDUM(I + (K - 1)*LX(RANK))
            ENDDO
          ENDDO
        ENDDO

        ! SEND TRANSFORMED DATA BACK
        CALL MPI_ALLTOALLV(RECVBUF, RECVCOUNTS, RDISPLS, COMPLEX_TYPE, SENDBUF, SENDCOUNTS, SDISPLS, COMPLEX_TYPE, PENCIL, IERR)

        ! REARRANGE DATA INTO ORIGINAL LAYOUT
        DO K = 1, LZ(RANK)
          DO I = 1, NX
            SPECTRUM(I, J, K) = SENDBUF(K + (I - 1)*LZ(RANK))
          ENDDO
        ENDDO

      ENDDO

      ! DESTROY PLAN
#ifdef DOUBLE_PREC
      CALL FFTW_DESTROY_PLAN(P1)
#else
      CALL FFTWF_DESTROY_PLAN(P1)
#endif

      ! RELEASE MEMORY
      NULLIFY(CDUM)

#ifdef DOUBLE_PREC
      CALL FFTW_FREE(PC)
#else
      CALL FFTWF_FREE(PC)
#endif

      DEALLOCATE(SENDBUF, RECVBUF, SDISPLS, RDISPLS, SENDCOUNTS, RECVCOUNTS)
      DEALLOCATE(I0, I1, K0, K1, LX, LZ)

      ! RELEASE COMMUNICATOR
      CALL MPI_COMM_FREE(PENCIL, IERR)

    END SUBROUTINE TRANSFORM_ALONG_Z

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE TRANSFORM_ALONG_Y(SPECTRUM)

      ! COMPUTE MANY IN-PLACE COMPLEX-TO-COMPLEX IFFT ALONG Y-DIRECTION ON INPUT SPECTRUM. PROCESSES ARE GATHERED IN Y-ORIENTED PENCILS
      ! BASED ON THE INPUT GEOMETRY. INSIDE EACH PENCIL, THE ALGORITHM WORKS ON Z-SLICES (TO REQUIRE AS LESS MEMORY AS POSSIBLE) TO
      ! EXCHANGE DATA AND COMPUTE THE IFFT.

      COMPLEX(FPP),              DIMENSION(:,:,:), INTENT(INOUT) :: SPECTRUM                   !< INPUT FIELD
      COMPLEX(FPP), ALLOCATABLE, DIMENSION(:)                 :: SENDBUF, RECVBUF           !< SENDER/RECEIVER BUFFER
      INTEGER(IPP)                                            :: I, J, K, P, C              !< COUNTERS
      INTEGER(IPP)                                            :: RANK, NTASKS, IERR         !< MPI STUFF
      INTEGER(IPP)                                            :: PENCIL
      INTEGER(IPP)                                            :: NX, NY, NZ                 !< TOTAL NUMBER OF POINTS FOR PENCIL ALONG EACH AXIS
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                 :: I0, I1                     !< LOCAL I-INDEX FOR PROCESSES BELONGING TO PENCIL
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                 :: J0, J1                     !< LOCAL J-INDEX FOR PROCESSES BELONGING TO PENCIL
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                 :: LX, LY                     !< LOCAL NUMBER OF POINTS ALONG X AND Y AXIS
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                 :: SENDCOUNTS, RECVCOUNTS     !< MPI STUFF
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                 :: SDISPLS, RDISPLS           !< MPI STUFF
      TYPE(C_PTR)                                             :: P1

      !-------------------------------------------------------------------------------------------------------------------------------

      ! GROUP PROCESSES IN PENCILS ORIENTED ALONG Y-AXIS. "LY" CONTAINS NUMBER OF POINTS FOR EACH PROCESS IN "PENCIL" ALONG Y
      CALL BUILD_PENCIL(1, PENCIL, RANK, NTASKS, LY)

      ! THIS ARRAY WILL CONTAIN THE NUMBER OF POINTS LOCAL TO EACH PROCESS IN CURRENT PENCIL ALONG X
      ALLOCATE(LX(0:NTASKS-1))

      ! NUMBER OF POINTS ALONG X AND Z. THESE ARE INDENTICAL FOR ALL PROCESSES IN CURRENT PENCIL.
      NX = SIZE(SPECTRUM, 1)
      NZ = SIZE(SPECTRUM, 3)

      ! TOTAL NUMBER OF POINTS ALONG Y-AXIS
      NY = SUM(LY)

      ALLOCATE(J0(0:NTASKS-1), J1(0:NTASKS-1))

      ! FIRST/LAST K-INDEX FOR EACH PROCESS
      DO P = 0, NTASKS - 1
        IF (P .EQ. 0) THEN
          J0(P) = 1
          J1(P) = LY(P)
        ELSE
          J0(P) = J1(P - 1) + 1
          J1(P) = J1(P - 1) + LY(P)
        ENDIF
      ENDDO

      ALLOCATE(I0(0:NTASKS-1), I1(0:NTASKS-1))

      ! DISTRIBUTE POINTS ALONG X-AXIS BETWEEN PROCESSES
      CALL MPI_SPLIT_TASK(NX, NTASKS, I0, I1)

      ! NUMBER OF POINTS ALONG X-AXIS FOR EACH PROCESS
      DO P = 0, NTASKS - 1
        LX(P) = I1(P) - I0(P) + 1
      ENDDO

      ALLOCATE(SENDCOUNTS(0:NTASKS-1), SDISPLS(0:NTASKS-1), RECVCOUNTS(0:NTASKS-1), RDISPLS(0:NTASKS-1))

      ! SET COUNTS/DISPLACEMENT ARRAYS
      DO P = 0, NTASKS - 1

        SENDCOUNTS(P) = LY(RANK) * LX(P)              !< POINTS SEND FROM PROCESS "RANK" TO PROCESS "P"
        RECVCOUNTS(P) = LY(P) * LX(RANK)              !< POINTS RECEIVED FROM PROCESS "RANK" BY PROCESS "P"

        IF (P .EQ. 0) THEN
          SDISPLS(P) = 0
          RDISPLS(P) = 0
        ELSE
          SDISPLS(P) = SDISPLS(P - 1) + SENDCOUNTS(P - 1)
          RDISPLS(P) = RDISPLS(P - 1) + RECVCOUNTS(P - 1)
        ENDIF

      ENDDO

      ALLOCATE(SENDBUF(NX * LY(RANK)))
      ALLOCATE(RECVBUF(NY * LX(RANK)))

      ! BUFFER "CDUM" WILL BE USED FOR TRANSFORM AND IT IS ALLOCATED USING FFTW MEMORY ALLOCATION UTILITY
#ifdef DOUBLE_PREC
      PC = FFTW_ALLOC_COMPLEX(INT(LX(RANK) * NY, C_SIZE_T))
#else
      PC = FFTWF_ALLOC_COMPLEX(INT(LX(RANK) * NY, C_SIZE_T))
#endif

      CALL C_F_POINTER(PC, CDUM, [LX(RANK) * NY])

      ! PREPARE FFTW PLAN ALONG Y-AXIS ("LX(RANK)" TRANSFORMS, EACH HAVING "NY" NUMBER OF POINTS)
#ifdef DOUBLE_PREC
      P1 = FFTW_PLAN_MANY_DFT(1, [NY], LX(RANK), CDUM, [NY], LX(RANK), 1, CDUM, [NY], LX(RANK), 1, FFTW_BACKWARD, FFTW_ESTIMATE)
#else
      P1 = FFTWF_PLAN_MANY_DFT(1, [NY], LX(RANK), CDUM, [NY], LX(RANK), 1, CDUM, [NY], LX(RANK), 1, FFTW_BACKWARD, FFTW_ESTIMATE)
#endif

      ! WORK Z-SLICES ONE BY ONE
      DO K = 1, NZ

        C = 0

        ! REARRANGE ELEMENTS OF "SPECTRUM"
        DO P = 0, NTASKS - 1
          DO J = 1, LY(RANK)
            DO I = I0(P), I1(P)
              C          = C + 1
              SENDBUF(C) = SPECTRUM(I, J, K)
            ENDDO
          ENDDO
        ENDDO

        ! EXCHANGE DATA
        CALL MPI_ALLTOALLV(SENDBUF, SENDCOUNTS, SDISPLS, COMPLEX_TYPE, CDUM, RECVCOUNTS, RDISPLS, COMPLEX_TYPE, PENCIL, IERR)

        ! IFFT
#ifdef DOUBLE_PREC
        CALL FFTW_EXECUTE_DFT(P1, CDUM, CDUM)
#else
        CALL FFTWF_EXECUTE_DFT(P1, CDUM, CDUM)
#endif

        C = 0

        ! REARRANGE ELEMENTS OF "CDUM" INTO "RECVBUF"
        DO P = 0, NTASKS - 1
          DO I = 1, LX(RANK)
            DO J = J0(P), J1(P)
              C          = C + 1
              RECVBUF(C) = CDUM(I + (J - 1)*LX(RANK))
            ENDDO
          ENDDO
        ENDDO

        ! SEND TRANSFORMED DATA BACK
        CALL MPI_ALLTOALLV(RECVBUF, RECVCOUNTS, RDISPLS, COMPLEX_TYPE, SENDBUF, SENDCOUNTS, SDISPLS, COMPLEX_TYPE, PENCIL, IERR)

        ! REARRANGE DATA INTO ORIGINAL LAYOUT
        DO J = 1, LY(RANK)
          DO I = 1, NX
            SPECTRUM(I, J, K) = SENDBUF(J + (I - 1)*LY(RANK))
          ENDDO
        ENDDO

      ENDDO

      ! DESTROY PLAN
#ifdef DOUBLE_PREC
      CALL FFTW_DESTROY_PLAN(P1)
#else
      CALL FFTWF_DESTROY_PLAN(P1)
#endif

      ! RELEASE MEMORY
      NULLIFY(CDUM)

#ifdef DOUBLE_PREC
      CALL FFTW_FREE(PC)
#else
      CALL FFTWF_FREE(PC)
#endif

      DEALLOCATE(SENDBUF, RECVBUF, SDISPLS, RDISPLS, SENDCOUNTS, RECVCOUNTS)
      DEALLOCATE(I0, I1, J0, J1, LX, LY)

      ! RELEASE COMMUNICATOR
      CALL MPI_COMM_FREE(PENCIL, IERR)

    END SUBROUTINE TRANSFORM_ALONG_Y

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE TRANSFORM_ALONG_X(SPECTRUM)

      ! COMPUTE MANY IN-PLACE COMPLEX-TO-REAL IFFT ALONG X-DIRECTION ON INPUT SPECTRUM. PROCESSES ARE GATHERED IN X-ORIENTED PENCILS
      ! BASED ON THE INPUT GEOMETRY. INSIDE EACH PENCIL, THE ALGORITHM WORKS ON Z-SLICES (TO REQUIRE AS LESS MEMORY AS POSSIBLE) TO
      ! EXCHANGE DATA AND COMPUTE THE IFFT.

      COMPLEX(FPP),              DIMENSION(:,:,:), INTENT(INOUT) :: SPECTRUM                       !< INPUT FIELD
      COMPLEX(FPP), ALLOCATABLE, DIMENSION(:)                    :: SENDBUF                        !< SENDER BUFFER
      INTEGER(IPP)                                               :: I, J, K, P, C                  !< COUNTERS
      INTEGER(IPP)                                               :: RANK, NTASKS, IERR             !< MPI STUFF
      INTEGER(IPP)                                               :: NYQUIST, PENCIL
      INTEGER(IPP)                                               :: IMAX, P_IMAX
      INTEGER(IPP)                                               :: NX, NY, NZ                     !< TOTAL NUMBER OF POINTS FOR PENCIL ALONG EACH AXIS
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: OFFSET
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: I0, I1                         !< LOCAL I-INDEX FOR PROCESSES BELONGING TO PENCIL
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: J0, J1                         !< LOCAL J-INDEX FOR PROCESSES BELONGING TO PENCIL
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: LX, LY                         !< LOCAL NUMBER OF POINTS ALONG X AND Y AXIS
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: SENDCOUNTS, RECVCOUNTS         !< MPI STUFF
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: SDISPLS, RDISPLS               !< MPI STUFF
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: R_SENDCOUNTS, R_RECVCOUNTS     !< MPI STUFF
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: R_SDISPLS, R_RDISPLS           !< MPI STUFF
      REAL(FPP),    ALLOCATABLE, DIMENSION(:)                    :: R_SENDBUF, RECVBUF
      TYPE(C_PTR)                                                :: P1

      !-------------------------------------------------------------------------------------------------------------------------------

      ! INDEX OF NYQUIST FREQUENCY ALONG X-AXIS
      NYQUIST = NPTS(1) / 2 + 1

      ! GROUP PROCESSES IN PENCILS ORIENTED ALONG X-AXIS. "LX" CONTAINS NUMBER OF POINTS FOR EACH PROCESS IN "PENCIL" ALONG X
      CALL BUILD_PENCIL(0, PENCIL, RANK, NTASKS, LX)

      ! PARAMETER USED TO DETERMINE ELEMENTS TO BE SENT/RECEIVED
      ALLOCATE(OFFSET(0:NTASKS-1))

      OFFSET(RANK) = GS(1, WORLD_RANK) - 1

      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, OFFSET, 1, MPI_INTEGER, PENCIL, IERR)

      ! THIS ARRAY WILL CONTAIN THE NUMBER OF POINTS LOCAL TO EACH PROCESS IN CURRENT PENCIL ALONG Y
      ALLOCATE(LY(0:NTASKS-1))

      ! NUMBER OF POINTS ALONG Y AND Z. THESE ARE INDENTICAL FOR ALL PROCESSES IN CURRENT PENCIL.
      NY = SIZE(SPECTRUM, 2)
      NZ = SIZE(SPECTRUM, 3)

      ! TOTAL NUMBER OF POINTS ALONG X-AXIS
      NX = SUM(LX)

      ! MAX I-INDEX USED TO SEND DATA
      IMAX = MIN(LX(RANK) + OFFSET(RANK), NYQUIST) - OFFSET(RANK)

      IF (IMAX .LT. 1) IMAX = 0

      ALLOCATE(I0(0:NTASKS-1), I1(0:NTASKS-1))

      ! FIRST/LAST I-INDEX FOR EACH PROCESS
      DO P = 0, NTASKS - 1
        IF (P .EQ. 0) THEN
          I0(P) = 1
          I1(P) = LX(P)
        ELSE
          I0(P) = I1(P - 1) + 1
          I1(P) = I1(P - 1) + LX(P)
        ENDIF
      ENDDO

      ALLOCATE(J0(0:NTASKS-1), J1(0:NTASKS-1))

      ! DISTRIBUTE POINTS ALONG Y-AXIS BETWEEN PROCESSES
      CALL MPI_SPLIT_TASK(NY, NTASKS, J0, J1)

      ! NUMBER OF POINTS ALONG Y-AXIS FOR EACH PROCESS
      DO P = 0, NTASKS - 1
        LY(P) = J1(P) - J0(P) + 1
      ENDDO

      ALLOCATE(SENDCOUNTS(0:NTASKS-1), SDISPLS(0:NTASKS-1), RECVCOUNTS(0:NTASKS-1), RDISPLS(0:NTASKS-1))
      ALLOCATE(R_SENDCOUNTS(0:NTASKS-1), R_SDISPLS(0:NTASKS-1), R_RECVCOUNTS(0:NTASKS-1), R_RDISPLS(0:NTASKS-1))

      ! SET COUNTS/DISPLACEMENT ARRAYS
      DO P = 0, NTASKS - 1

        ! DEFINE "IMAX" FOR PROCESS "P"
        P_IMAX = MIN(LX(P) + OFFSET(P), NYQUIST) - OFFSET(P)

        IF (P_IMAX .LT. 1) P_IMAX = 0

        ! PROCESS COMPLETELY TO THE RIGHT OF "NYQUIST" WON'T SEND ANY ELEMENT ("SENDCOUNTS=0")
        SENDCOUNTS(P) = LY(P) * IMAX                                            !< POINTS SEND FROM PROCESS "RANK" TO PROCESS "P"

        ! A PROCESS WON'T RECEIVE ANYTHING FROM PROCESS "P" IF THE LATTER IS TO THE RIGHT OF "NYQUIST" ("RECVCOUNTS=0")
        RECVCOUNTS(P) = LY(RANK) * P_IMAX                                       !< POINTS RECEIVED FROM PROCESS "RANK" BY PROCESS "P"

        R_SENDCOUNTS(P) = LY(RANK) * LX(P)                                      !< POINTS SEND FROM PROCESS "RANK" TO PROCESS "P"
        R_RECVCOUNTS(P) = LY(P) * LX(RANK)                                      !< POINTS RECEIVED FROM PROCESS "RANK" BY PROCESS "P"

        IF (P .EQ. 0) THEN
          SDISPLS(P)   = 0
          RDISPLS(P)   = 0
          R_SDISPLS(P) = 0
          R_RDISPLS(P) = 0
        ELSE
          SDISPLS(P)   = SDISPLS(P - 1) + SENDCOUNTS(P - 1)
          RDISPLS(P)   = RDISPLS(P - 1) + RECVCOUNTS(P - 1)
          R_SDISPLS(P) = R_SDISPLS(P - 1) + R_SENDCOUNTS(P - 1)
          R_RDISPLS(P) = R_RDISPLS(P - 1) + R_RECVCOUNTS(P - 1)
        ENDIF

      ENDDO

      ALLOCATE(SENDBUF(NY * IMAX), R_SENDBUF(NX * LY(RANK)))
      ALLOCATE(RECVBUF(NY * LX(RANK)))

      ! ALLOCATE MEMORY USING FFTW ALLOCATION UTILITY
#ifdef DOUBLE_PREC
      PC = FFTW_ALLOC_COMPLEX(INT(LY(RANK) * NYQUIST, C_SIZE_T))
#else
      PC = FFTWF_ALLOC_COMPLEX(INT(LY(RANK) * NYQUIST, C_SIZE_T))
#endif

      ! FOR IN-PLACE TRANSFORMS, OUTPUT REAL ARRAY IS SLIGHTLY LONGER THAN ACTUAL PHYSICAL DIMENSION
      CALL C_F_POINTER(PC, CDUM, [LY(RANK) * NYQUIST])
      CALL C_F_POINTER(PC, RDUM, [LY(RANK) * NYQUIST * 2])

      ! PREPARE FFTW PLAN ALONG X-AXIS ("LY(RANK)" TRANSFORMS, EACH HAVING "NPTS(1)" NUMBER OF POINTS)
#ifdef DOUBLE_PREC
      P1 = FFTW_PLAN_MANY_DFT_C2R(1, [NPTS(1)], LY(RANK), CDUM, [NPTS(1)], LY(RANK), 1, RDUM, [NPTS(1)], LY(RANK), 1, FFTW_ESTIMATE)
#else
      P1 = FFTWF_PLAN_MANY_DFT_C2R(1, [NPTS(1)], LY(RANK), CDUM, [NPTS(1)], LY(RANK), 1, RDUM, [NPTS(1)], LY(RANK), 1,FFTW_ESTIMATE)
#endif

      ! WORK Z-SLICES ONE BY ONE
      DO K = 1, NZ

        C = 0

        ! REARRANGE ELEMENTS OF "SPECTRUM"
        DO P = 0, NTASKS - 1
          DO I = 1, IMAX
            DO J = J0(P), J1(P)
              C          = C + 1
              SENDBUF(C) = SPECTRUM(I, J, K)
            ENDDO
          ENDDO
        ENDDO

        ! EXCHANGE DATA
        CALL MPI_ALLTOALLV(SENDBUF, SENDCOUNTS, SDISPLS, COMPLEX_TYPE, CDUM, RECVCOUNTS, RDISPLS, COMPLEX_TYPE, PENCIL, IERR)

        ! IFFT
#ifdef DOUBLE_PREC
        CALL FFTW_EXECUTE_DFT_C2R(P1, CDUM, RDUM)
#else
        CALL FFTWF_EXECUTE_DFT_C2R(P1, CDUM, RDUM)
#endif

        C = 0

        ! REARRANGE ELEMENTS OF "RDUM" INTO "RECVBUF"
        DO P = 0, NTASKS - 1
          DO J = 1, LY(RANK)
            DO I = I0(P), I1(P)
              C            = C + 1
              R_SENDBUF(C) = RDUM(J + (I - 1)*LY(RANK))
            ENDDO
          ENDDO
        ENDDO

        ! SEND TRANSFORMED DATA BACK
        CALL MPI_ALLTOALLV(R_SENDBUF, R_SENDCOUNTS, R_SDISPLS, REAL_TYPE, RECVBUF, R_RECVCOUNTS, R_RDISPLS, REAL_TYPE, PENCIL, IERR)

        ! REARRANGE DATA INTO ORIGINAL LAYOUT
        DO J = 1, NY
          DO I = 1, LX(RANK)
            SPECTRUM(I, J, K) = CMPLX(RECVBUF(I + (J - 1)*LX(RANK)), 0._FPP, FPP)
          ENDDO
        ENDDO

      ENDDO

      ! DESTROY PLAN
#ifdef DOUBLE_PREC
      CALL FFTW_DESTROY_PLAN(P1)
#else
      CALL FFTWF_DESTROY_PLAN(P1)
#endif

      ! RELEASE MEMORY
      NULLIFY(CDUM, RDUM)

#ifdef DOUBLE_PREC
      CALL FFTW_FREE(PC)
#else
      CALL FFTWF_FREE(PC)
#endif

      DEALLOCATE(SENDBUF, RECVBUF, SDISPLS, RDISPLS, SENDCOUNTS, RECVCOUNTS)
      DEALLOCATE(R_SENDBUF, R_SDISPLS, R_RDISPLS, R_SENDCOUNTS, R_RECVCOUNTS)
      DEALLOCATE(I0, I1, J0, J1, LX, LY, OFFSET)

      ! RELEASE COMMUNICATOR
      CALL MPI_COMM_FREE(PENCIL, IERR)

    END SUBROUTINE TRANSFORM_ALONG_X

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE BUILD_PENCIL(DIR, NEWCOMM, RANK, NTASKS, N)

      ! GROUP PROCESSORS IN PENCILS ORIENTED ALONG A SPECIFIC DIRECTION.

      INTEGER(IPP),                                        INTENT(IN)  :: DIR                    !< STENCIL DIRECTION (0=X, 1=Y, 2=Z)
      INTEGER(IPP),                                        INTENT(OUT) :: NEWCOMM                !< HANDLE TO NEW COMMUNICATOR
      INTEGER(IPP),                                        INTENT(OUT) :: RANK                   !< RANK OF CALLING PROCESS IN NEW COMMUNICATOR
      INTEGER(IPP),                                        INTENT(OUT) :: NTASKS                 !< NUMBER OF PROCESSES IN NEW COMMUNICATOR
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:),             INTENT(OUT) :: N                      !< NUMBER OF POINTS PER PROCESS IN PENCIL DIRECTION
      INTEGER(IPP)                                                     :: I, IERR
      INTEGER(IPP),              DIMENSION(0:WORLD_SIZE-1)             :: COLOR
      LOGICAL,                   DIMENSION(2)                          :: BOOL

      !-----------------------------------------------------------------------------------------------------------------------------

      ! GROUP PROCESSES IN STENCILS
      DO I = 0, WORLD_SIZE - 1

        COLOR(I) = 0

        ! STENCIL ALONG X-AXIS
        IF (DIR .EQ. 0) THEN
          BOOL(1) = (GS(2, I) .EQ. GS(2, WORLD_RANK)) .AND. (GE(2, I) .EQ. GE(2, WORLD_RANK))
          BOOL(2) = (GS(3, I) .EQ. GS(3, WORLD_RANK)) .AND. (GE(3, I) .EQ. GE(3, WORLD_RANK))
        ! STENCIL ALONG Y-AXIS
        ELSEIF (DIR .EQ. 1) THEN
          BOOL(1) = (GS(1, I) .EQ. GS(1, WORLD_RANK)) .AND. (GE(1, I) .EQ. GE(1, WORLD_RANK))
          BOOL(2) = (GS(3, I) .EQ. GS(3, WORLD_RANK)) .AND. (GE(3, I) .EQ. GE(3, WORLD_RANK))
        ! STENCIL ALONG Z-AXIS
        ELSEIF (DIR .EQ. 2) THEN
          BOOL(1) = (GS(1, I) .EQ. GS(1, WORLD_RANK)) .AND. (GE(1, I) .EQ. GE(1, WORLD_RANK))
          BOOL(2) = (GS(2, I) .EQ. GS(2, WORLD_RANK)) .AND. (GE(2, I) .EQ. GE(2, WORLD_RANK))
        ENDIF

        IF (ALL(BOOL .EQV. .TRUE.)) COLOR(I) = I + 1

      ENDDO

      ! PROCESS BELONGING TO THE SAME STENCIL HAVE SAME COLOR
      COLOR(WORLD_RANK) = MAXVAL(COLOR, DIM = 1)

      ! CREATE COMMUNICATOR SUBGROUP
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, COLOR(WORLD_RANK), WORLD_RANK, NEWCOMM, IERR)

      ! PROCESS ID AND COMMUNICATOR SIZE
      CALL MPI_COMM_RANK(NEWCOMM, RANK, IERR)
      CALL MPI_COMM_SIZE(NEWCOMM, NTASKS, IERR)

      ALLOCATE(N(0:NTASKS - 1))

      ! NUMBER OF POINTS ALONG DIRECTION "DIR" FOR CALLING PROCESS
      N(RANK) = GE(DIR + 1, WORLD_RANK) - GS(DIR + 1, WORLD_RANK) + 1

      ! MAKE WHOLE COMMUNICATOR AWARE OF POINTS FOR EACH PROCESS
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, N, 1, MPI_INTEGER, NEWCOMM, IERR)

    END SUBROUTINE BUILD_PENCIL

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE COMPUTE_SPECTRUM(LS, LE, DH, ACF, CL, SIGMA, HURST, SEED, SPEC, TIME)

      ! COMPUTE THE SPECTRUM OF A RANDOM FIELD CHARACTERISED BY A SPECIFIC AUTOCORRELATION FUNCTION, CORRELATION LENGTH, STANDARD DEVIATION,
      ! AND HURST EXPONENT (NOT USED FOR GAUSSIAN FIELDS). BASED ON THE FOURIER INTEGRAL METHOD OF PARDO-IGUZQUIZA AND CHICA-OLMO.

      INTEGER(IPP),     DIMENSION(3),                                   INTENT(IN)  :: LS, LE                !< GLOBAL INDICES
      REAL(FPP),                                                        INTENT(IN)  :: DH                    !< GRID-STEP
      CHARACTER(LEN=*),                                                 INTENT(IN)  :: ACF                   !< AUTOCORRELATION FUNCTION
      REAL(FPP),        DIMENSION(3),                                   INTENT(IN)  :: CL                    !< CORRELATION LENGTH (CAN BE ANISOTROPIC)
      REAL(FPP),                                                        INTENT(IN)  :: SIGMA                 !< STANDARD DEVIATION (SIGMA%/100)
      REAL(FPP),                                                        INTENT(IN)  :: HURST                 !< HURST EXPONENT
      INTEGER(IPP),                                                     INTENT(IN)  :: SEED                  !< INITIAL SEED NUMBER
      COMPLEX(FPP),     DIMENSION(LS(1):LE(1),LS(2):LE(2),LS(3):LE(3)), INTENT(OUT) :: SPEC                  !< SPECTRUM
      REAL(FPP),        DIMENSION(2),                                   INTENT(OUT) :: TIME                  !< ELAPSED TIME
      INTEGER(IPP)                                                                  :: I, J, K               !< INDICES
      LOGICAL,          DIMENSION(3)                                                :: BOOL
      REAL(FPP)                                                                     :: BUTTER, NUM, AMP      !< USED TO COMPUTE SPECTRUM
      REAL(FPP)                                                                     :: KNYQ, KC, KR          !< USED TO COMPUTE SPECTRUM
      !REAL(FPP)                                                                     :: R
      REAL(REAL64)                                                                  :: TICTOC                !< USED FOR TIMING
      REAL(FPP),        DIMENSION(3)                                                :: DK                    !< RESOLUTION IN WAVENUMBER DOMAIN
      !REAL(FPP),        DIMENSION(NPTS(1))                                          :: HARVEST               !< RANDOM VALUES
      REAL(FPP),        DIMENSION(NPTS(1))                                          :: KX                    !< WAVENUMBER VECTOR (X)
      REAL(FPP),        DIMENSION(NPTS(2))                                          :: KY                    !< WAVENUMBER VECTOR (Y)
      REAL(FPP),        DIMENSION(NPTS(3))                                          :: KZ                    !< WAVENUMBER VECTOR (Z)
      REAL(C_FPP),      DIMENSION(:),                                   POINTER     :: R
      REAL(C_FPP),      DIMENSION(:,:,:),                               POINTER     :: HARVEST
      TYPE(C_PTR)                                                                   :: CPTR

      !-----------------------------------------------------------------------------------------------------------------------------

      ! RESOLUTION IN WAVENUMBER DOMAIN ALONG EACH DIRECTION
      DO I = 1, 3
        DK(I) = 2._FPP * PI / (REAL(NPTS(I), FPP) * DH)
      ENDDO

      ! NYQUIST WAVENUMBER
      KNYQ = PI / DH

      ! CORNER WAVENUMBER FOR SPECTRUM TAPERING TO AVOID ALIASING IS HALF NYQUIST WAVENUMBER
      KC = KNYQ / 2._FPP

      ! VECTORS GO FROM 0 TO NYQUIST AND THEN BACK AGAIN UNTIL DK
      KX = [[(I * DK(1), I = 0, NPTS(1)/2)], [(I * DK(1), I = NPTS(1)/2-1, 1, -1)]]
      KY = [[(J * DK(2), J = 0, NPTS(2)/2)], [(J * DK(2), J = NPTS(2)/2-1, 1, -1)]]
      KZ = [[(K * DK(3), K = 0, NPTS(3)/2)], [(K * DK(3), K = NPTS(3)/2-1, 1, -1)]]

      ! COMPUTE PART OF POWER SPECTRAL DENSITY OUTSIDE LOOP
      IF (ACF .EQ. 'VK') THEN
        NUM = 8._FPP * SQRT(PI**3) * GAMMA(HURST + 1.5_FPP) * SIGMA**2 * PRODUCT(CL) / GAMMA(HURST)
        FUN => VK_PSDF
      ELSEIF (ACF .EQ. 'GAUSS') THEN
        NUM = SIGMA**2 * PRODUCT(CL) * SQRT(PI**3)
        FUN => GS_PSDF
      ENDIF

      ! START TIMER
      CALL WATCH_START(TICTOC)

      ! EACH PROCESS GENERATE ITS SET OF RANDOM NUMBERS
      CPTR = PRNG(SEED, LS, LE, NPTS)

      CALL C_F_POINTER(CPTR, R, [(LE(1) - LS(1) + 1) * (LE(2) - LS(2) + 1) * (LE(3) - LS(3) + 1)])

      HARVEST(LS(1):LE(1), LS(2):LE(2), LS(3):LE(3)) => R

      ! COMPUTE SPECTRUM
      DO K = LS(3), LE(3)

        DO J = LS(2), LE(2)

          DO I = LS(1), LE(1)

            ! RADIAL WAVENUMBER
            !KR = SQRT(KX(I)**2 + KY(J)**2 + KZ(K)**2)

            ! FOURTH-ORDER LOW-PASS BUTTERWORTH FILTER
            !BUTTER = 1._FPP / SQRT(1._FPP + (KR / KC)**(2 * 4))

            ! NOW "KR" IS THE PRODUCT "K * CL"
            KR = (KX(I) * CL(1))**2 + (KY(J) * CL(2))**2 + (KZ(K) * CL(3))**2

            ! COMPLETE POWER SPECTRAL DENSITY AND GO TO AMPLITUDE SPECTRUM
            AMP = SQRT(NUM / FUN(KR, HURST))

            ! APPLY FILTER
            !AMP = AMP * BUTTER

            !R = HARVEST(I - LS(1) + 1, J - LS(2) + 1, K - LS(3) + 1) * 2._FPP * PI
            HARVEST(I, J, K) = HARVEST(I, J, K) * 2._FPP * PI

            ! COMBINE AMPLITUDE AND RANDOM PHASE
            SPEC(I, J, K) = CMPLX(COS(HARVEST(I, J, K)) * AMP, SIN(HARVEST(I, J, K)) * AMP, FPP)

          ENDDO
        ENDDO
      ENDDO

      NULLIFY(R, HARVEST)
      CALL FREE_MEM(CPTR)

      CALL WATCH_STOP(TICTOC)

      TIME(1) = TICTOC

      CALL WATCH_START(TICTOC)

      CALL ENFORCE_SYMMETRY(LS, LE, SPEC)

      CALL WATCH_STOP(TICTOC)

      TIME(2) = TICTOC

    END SUBROUTINE COMPUTE_SPECTRUM

    ! SUBROUTINE COMPUTE_SPECTRUM(LS, LE, DH, ACF, CL, SIGMA, HURST, SEED, SPEC, TIME)
    !
    !   ! COMPUTE THE SPECTRUM OF A RANDOM FIELD CHARACTERISED BY A SPECIFIC AUTOCORRELATION FUNCTION, CORRELATION LENGTH, STANDARD DEVIATION,
    !   ! AND HURST EXPONENT (NOT USED FOR GAUSSIAN FIELDS). BASED ON THE FOURIER INTEGRAL METHOD OF PARDO-IGUZQUIZA AND CHICA-OLMO.
    !
    !   INTEGER(IPP),     DIMENSION(3),                                   INTENT(IN)  :: LS, LE                !< GLOBAL INDICES
    !   REAL(FPP),                                                        INTENT(IN)  :: DH                    !< GRID-STEP
    !   CHARACTER(LEN=*),                                                 INTENT(IN)  :: ACF                   !< AUTOCORRELATION FUNCTION
    !   REAL(FPP),        DIMENSION(3),                                   INTENT(IN)  :: CL                    !< CORRELATION LENGTH (CAN BE ANISOTROPIC)
    !   REAL(FPP),                                                        INTENT(IN)  :: SIGMA                 !< STANDARD DEVIATION (SIGMA%/100)
    !   REAL(FPP),                                                        INTENT(IN)  :: HURST                 !< HURST EXPONENT
    !   INTEGER(IPP),                                                     INTENT(IN)  :: SEED                  !< INITIAL SEED NUMBER
    !   COMPLEX(FPP),     DIMENSION(LS(1):LE(1),LS(2):LE(2),LS(3):LE(3)), INTENT(OUT) :: SPEC                  !< SPECTRUM
    !   REAL(FPP),        DIMENSION(2),                                   INTENT(OUT) :: TIME                  !< ELAPSED TIME
    !   INTEGER(IPP)                                                                  :: I, J, K               !< INDICES
    !   LOGICAL,          DIMENSION(3)                                                :: BOOL
    !   REAL(FPP)                                                                     :: BUTTER, NUM, AMP      !< USED TO COMPUTE SPECTRUM
    !   REAL(FPP)                                                                     :: KNYQ, KC, KR          !< USED TO COMPUTE SPECTRUM
    !   REAL(REAL64)                                                                  :: TICTOC                !< USED FOR TIMING
    !   REAL(FPP),        DIMENSION(3)                                                :: DK                    !< RESOLUTION IN WAVENUMBER DOMAIN
    !   REAL(FPP),        DIMENSION(NPTS(1))                                          :: HARVEST               !< RANDOM VALUES
    !   REAL(FPP),        DIMENSION(NPTS(1))                                          :: KX                    !< WAVENUMBER VECTOR (X)
    !   REAL(FPP),        DIMENSION(NPTS(2))                                          :: KY                    !< WAVENUMBER VECTOR (Y)
    !   REAL(FPP),        DIMENSION(NPTS(3))                                          :: KZ                    !< WAVENUMBER VECTOR (Z)
    !
    !   !-----------------------------------------------------------------------------------------------------------------------------
    !
    !   ! RESOLUTION IN WAVENUMBER DOMAIN ALONG EACH DIRECTION
    !   DO I = 1, 3
    !     DK(I) = 2._FPP * PI / (REAL(NPTS(I), FPP) * DH)
    !   ENDDO
    !
    !   ! NYQUIST WAVENUMBER
    !   KNYQ = PI / DH
    !
    !   ! CORNER WAVENUMBER FOR SPECTRUM TAPERING TO AVOID ALIASING IS HALF NYQUIST WAVENUMBER
    !   KC = KNYQ / 2._FPP
    !
    !   ! VECTORS GO FROM 0 TO NYQUIST AND THEN BACK AGAIN UNTIL DK
    !   KX = [[(I * DK(1), I = 0, NPTS(1)/2)], [(I * DK(1), I = NPTS(1)/2-1, 1, -1)]]
    !   KY = [[(J * DK(2), J = 0, NPTS(2)/2)], [(J * DK(2), J = NPTS(2)/2-1, 1, -1)]]
    !   KZ = [[(K * DK(3), K = 0, NPTS(3)/2)], [(K * DK(3), K = NPTS(3)/2-1, 1, -1)]]
    !
    !   ! COMPUTE PART OF POWER SPECTRAL DENSITY OUTSIDE LOOP
    !   IF (ACF .EQ. 'VK') THEN
    !     NUM = 8._FPP * SQRT(PI**3) * GAMMA(HURST + 1.5_FPP) * SIGMA**2 * PRODUCT(CL) / GAMMA(HURST)
    !     FUN => VK_PSDF
    !   ELSEIF (ACF .EQ. 'GAUSS') THEN
    !     NUM = SIGMA**2 * PRODUCT(CL) * SQRT(PI**3)
    !     FUN => GS_PSDF
    !   ENDIF
    !
    !   ! INITIALISE RANDOM GENERATOR
    !   CALL SET_STREAM(SEED)
    !
    !   ! START TIMER
    !   CALL WATCH_START(TICTOC)
    !
    !   ! COMPUTE SPECTRUM
    !   DO K = 1, NPTS(3)
    !
    !     BOOL(3) = (K .GE. LS(3)) .AND. (K .LE. LE(3))
    !
    !     DO J = 1, NPTS(2)
    !
    !       BOOL(2) = (J .GE. LS(2)) .AND. (J .LE. LE(2))
    !
    !       CALL RANDOM_NUMBER(HARVEST)
    !
    !       IF ( BOOL(2) .AND. BOOL(3) ) THEN
    !
    !       !DO I = 1, NPTS(1)
    !       DO I = LS(1), LE(1)
    !
    !         !BOOL(1) = (I .GE. LS(1)) .AND. (I .LE. LE(1))
    !
    !         !IF (ALL(BOOL .EQV. .TRUE.)) THEN
    !
    !           ! RADIAL WAVENUMBER
    !           !KR = SQRT(KX(I)**2 + KY(J)**2 + KZ(K)**2)
    !
    !           ! FOURTH-ORDER LOW-PASS BUTTERWORTH FILTER
    !           !BUTTER = 1._FPP / SQRT(1._FPP + (KR / KC)**(2 * 4))
    !
    !           ! NOW "KR" IS THE PRODUCT "K * CL"
    !           KR = (KX(I) * CL(1))**2 + (KY(J) * CL(2))**2 + (KZ(K) * CL(3))**2
    !
    !           ! COMPLETE POWER SPECTRAL DENSITY AND GO TO AMPLITUDE SPECTRUM
    !           AMP = SQRT(NUM / FUN(KR, HURST))
    !
    !           ! APPLY FILTER
    !           !AMP = AMP * BUTTER
    !
    !           HARVEST(I) = HARVEST(I) * 2._FPP * PI
    !
    !           ! COMBINE AMPLITUDE AND RANDOM PHASE
    !           SPEC(I, J, K) = CMPLX(COS(HARVEST(I)) * AMP, SIN(HARVEST(I)) * AMP, FPP)
    !
    !         !ENDIF
    !
    !       ENDDO
    !
    !       ENDIF
    !
    !     ENDDO
    !   ENDDO
    !
    !   CALL WATCH_STOP(TICTOC)
    !
    !   TIME(1) = TICTOC
    !
    !   CALL WATCH_START(TICTOC)
    !
    !   CALL ENFORCE_SYMMETRY(LS, LE, SPEC)
    !
    !   CALL WATCH_STOP(TICTOC)
    !
    !   TIME(2) = TICTOC
    !
    ! END SUBROUTINE COMPUTE_SPECTRUM

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE MPI_SPLIT_TASK(NPTS, NTASKS, I0, I1)

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

    END SUBROUTINE MPI_SPLIT_TASK

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(FPP) FUNCTION VK_PSDF(M, HURST)

      ! FUNCTION USED TO GENERATE VON KARMAN (AND EXPONENTIAL) RANDOM FIELDS

      REAL(FPP), INTENT(IN) :: M                  !< PRODUCT (WAVENUMBER*CL)**2
      REAL(FPP), INTENT(IN) :: HURST              !< HURST EXPONENT

      !-----------------------------------------------------------------------------------------------------------------------------

      VK_PSDF = (1._FPP + M)**(HURST + 1.5_FPP)

    END FUNCTION VK_PSDF

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(FPP) FUNCTION GS_PSDF(M, DUMMY)

      ! FUNCTION USED TO GENERATE GAUSSIAN RANDOM FIELDS

      REAL(FPP), INTENT(IN) :: M                  !< PRODUCT (WAVENUMBER*CL)**2
      REAL(FPP), INTENT(IN) :: DUMMY              !< THIS VARIABLE IS NOT USED

      !-----------------------------------------------------------------------------------------------------------------------------

      GS_PSDF = EXP(0.25_FPP * M)

    END FUNCTION GS_PSDF

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE ENFORCE_SYMMETRY(FS, FE, SPEC)

      ! MAKE SURE THE COMPLEX SPECTRUM HAS HERMITIAN SYMMETRY. SYMMETRY CONDITIONS ARE ENFORCED AT "I=1" AND "I=NYQUIST" ON THE YZ-PLANE;
      ! SYMMETRY ALONG X-DIRECTION IS IMPLICIT AS WE MAKE USE OF A COMPLEX-TO-REAL IFFT. THE ALGORITHM RUNS OVER EACH PROCESS CUT BY
      ! THE YZ-PLANE(S) STEP BY STEP, IN ORDER TO REQUIRE AS LESS MEMORY AS POSSIBLE (THIS RESULTS IN A PERFORMANCE PERNALTY SINCE MORE
      ! COMMUNICATIONS ARE NEEDED). SEE FIGURE 4 AND 5 OF PARDO-IGUZQUIZA AND CHICA-OLMO FOR A BETTER UNDERSTANDING OF THE SYMMETRY
      ! CONDITIONS.

      INTEGER(IPP),              DIMENSION(3),                                   INTENT(IN)    :: FS, FE                       !< LOCAL INDICES
      COMPLEX(FPP),              DIMENSION(FS(1):FE(1),FS(2):FE(2),FS(3):FE(3)), INTENT(INOUT) :: SPEC                         !< SPECTRUM
      COMPLEX(FPP), ALLOCATABLE, DIMENSION(:)                                                  :: BUFFER                       !< BUFFER FOR MESSAGING
      INTEGER(IPP)                                                                             :: C, I, J, K, P                !< COUNTERS
      INTEGER(IPP)                                                                             :: I0, J0, J1, K0, K1           !< VARIOUS INDICES
      INTEGER(IPP)                                                                             :: NTASKS, RANK, COLOR, SLICE   !< MPI STUFF
      INTEGER(IPP)                                                                             :: COMM, IERR, KEY              !< MPI STUFF
      INTEGER(IPP),              DIMENSION(3)                                                  :: NYQUIST                      !< INDEX NYQUIST WAVENUMBER
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                                                  :: SLICE2GLOBAL, JMAP, KMAP     !< ARRAYS FOR MESSAGING
      LOGICAL                                                                                  :: BOOL

      !-----------------------------------------------------------------------------------------------------------------------------

      ! INDEX OF NYQUIST WAVENUMBER
      DO I = 1, 3
        NYQUIST(I) = NPTS(I) / 2 + 1
      ENDDO

      ! CONSIDER ONLY I-INDEX = 1 AND I-INDEX = NYQUIST(1)
      DO I = 1, NYQUIST(1), NYQUIST(1) - 1

        ! "TRUE" IF PROCESS CONTAINS POINT "ALPHA"
        BOOL = (FS(1) .EQ. 1) .AND. (FS(2) .EQ. 1) .AND. (FS(3) .EQ. 1)

        IF (BOOL) SPEC(1, 1, 1) = 0._FPP                                                         !< ALPHA MUST BE ZERO (ZERO-MEAN FIELD)

        ! "TRUE" IF PROCESS CONTAINS POINT "ALPHA"
        BOOL = (FS(1) .LE. I) .AND. (FE(1) .GE. I) .AND. (FS(2) .EQ. 1) .AND. (FS(3) .EQ. 1)

        IF (BOOL) SPEC(I, 1, 1) = REAL(SPEC(I, 1, 1), FPP)                                       !< ALPHA MUST BE REAL

        ! "TRUE" IF PROCESS CONTAINS POINT "PSI"
        BOOL = (FS(1) .LE. I) .AND. (FE(1) .GE. I) .AND. (FS(2) .LE. NYQUIST(2)) .AND. (FE(2) .GE. NYQUIST(2)) .AND.      &
               (FS(3) .LE. NYQUIST(3)) .AND. (FE(3) .GE. NYQUIST(3))

        IF (BOOL) SPEC(I, NYQUIST(2), NYQUIST(3)) = REAL(SPEC(I, NYQUIST(2), NYQUIST(3)), FPP)   !< PSI MUST BE REAL

        ! "TRUE" IF PROCESS CONTAINS POINT "XSI"
        BOOL = (FS(1) .LE. I) .AND. (FE(1) .GE. I) .AND. (FS(2) .EQ. 1) .AND. (FS(3) .LE. NYQUIST(3)) .AND. (FE(3) .GE. NYQUIST(3))

        IF (BOOL) SPEC(I, 1, NYQUIST(3)) = REAL(SPEC(I, 1, NYQUIST(3)), FPP)                     !< XSI MUST BE REAL

        ! "TRUE" IF PROCESS CONTAINS POINT "ETA"
        BOOL = (FS(1) .LE. I) .AND. (FE(1) .GE. I) .AND. (FS(2) .LE. NYQUIST(2)) .AND. (FE(2) .GE. NYQUIST(2)) .AND. (FS(3) .EQ. 1)

        IF (BOOL) SPEC(I, NYQUIST(2), 1) = REAL(SPEC(I, NYQUIST(2), 1), FPP)                     !< ETA MUST BE REAL

        ! "TRUE" ONLY FOR THOSE PROCESSES CONTAINING "I=1" OR "I=NYQUIST"
        BOOL = (FS(1) .LE. I) .AND. (FE(1) .GE. I)

        COLOR = 0

        IF (BOOL) COLOR = I

        ! FIRST SPLIT OF ORIGINAL COMMUNICATOR
        CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, COLOR, WORLD_RANK, SLICE, IERR)

        ! GLOBAL TO LOCAL MAPPING FOR I-INDEX
        I0 = I - FS(1) + 1

        ! ONLY PROCESSES BELONGING TO THE (YZ-)SLICE (I.E. WITH "COLOR=I") ENTER THIS IF STATEMENT
        IF (COLOR .EQ. I) THEN

          ! NUMBER OF PROCESSES IN THE NEW COMMUNICATOR
          CALL MPI_COMM_SIZE(SLICE, NTASKS, IERR)

          CALL MPI_COMM_RANK(SLICE, RANK, IERR)

          ALLOCATE(SLICE2GLOBAL(0:NTASKS-1))

          ! SLICE-RANK TO GLOBAL-RANK MAPPING
          SLICE2GLOBAL(RANK) = WORLD_RANK

          ! UPDATE GLOBAL-RANK MAPPING FOR ALL PROCESSES
          CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_INTEGER, SLICE2GLOBAL, 1, MPI_INTEGER, SLICE, IERR)

          ! LOOP OVER EACH PROCESS IN THE "SLICE" COMMUNICATOR: THOSE PROCESSES WHOSE DATA ARE NEEDED FOR COMPLEX CONJUGATION IN PROCESS
          ! "P" WILL BE PART OF A NEW COMMUNICATOR
          DO P = 0, NTASKS - 1

            ! ENTER STATEMENT IF PROCESS "P" HAS POINTS IN THE EPSILON*/THETA*/PHI* REGION
            IF ( (GE(2, SLICE2GLOBAL(P)) .GE. NYQUIST(2)) .AND. (GE(3, SLICE2GLOBAL(P)) .GE. NYQUIST(3)) ) THEN

              ! CORRESPONDING MIN/MAX INDICES IN THE EPSILON/THETA/PHI REGION
              J0 = NPTS(2) - GE(2, SLICE2GLOBAL(P)) + 2
              J1 = NPTS(2) - MAX(GS(2, SLICE2GLOBAL(P)), NYQUIST(2)) + 2

              K0 = NPTS(3) - GE(3, SLICE2GLOBAL(P)) + 2
              K1 = NPTS(3) - MAX(GS(3, SLICE2GLOBAL(P)), NYQUIST(3)) + 2

              ! "TRUE" IF CALLING PROCESS HAS DATA NEEDED BY PROCESS "P"
              BOOL = (( (FS(2) .LE. J0) .AND. (J0 .LE. FE(2)) ) .OR. ( (FS(2) .LE. J1) .AND. (J1 .LE. FE(2)) ) .OR.      &
                      ( (FS(2) .GE. J0) .AND. (FE(2) .LE. J1) ))                                               .AND.     &
                     (( (FS(3) .LE. K0) .AND. (K0 .LE. FE(3)) ) .OR. ( (FS(3) .LE. K1) .AND. (K1 .LE. FE(3)) ) .OR.      &
                      ( (FS(3) .GE. K0) .AND. (FE(3) .LE. K1) ))

              COLOR = 0

              ! CALLING PROCESS IS PART OF NEW COMMUNICATOR
              IF (BOOL) THEN
                COLOR = 1
                KEY   = RANK
              ENDIF

              ! ADD PROCESS "P" TO COMMUNICATOR AS IT WILL COLLECT DATA
              IF (P .EQ. RANK) THEN
                COLOR = 1
                KEY   = -1                                                  !< MAKE SURE IT GETS LOWEST (0) RANK VALUE
              ENDIF

              ! SPLIT "SLICE" COMMUNICATOR
              CALL MPI_COMM_SPLIT(SLICE, COLOR, KEY, COMM, IERR)

              ! ONLY PROCESSES IN THE NEW COMMUNICATOR EXECUTE STATEMENTS BELOW
              IF (COLOR .EQ. 1) THEN

                ! GATHER DATA FROM PROCESSES
                CALL GATHER_DATA(COMM, I0, SPEC, JMAP, KMAP, BUFFER)

                DO C = 1, SIZE(BUFFER)

                  ! "TRUE" IF CURRENT POINT IS INSIDE THE EPSILON/THETA/PHI REGION
                  BOOL = (JMAP(C) .GE. J0) .AND. (JMAP(C) .LE. J1) .AND. (KMAP(C) .GE. K0) .AND. (KMAP(C) .LE. K1)

                  ! INDICES WHERE COMPLEX-CONJUGATED VALUES MUST BE PLACED
                  J = NPTS(2) - JMAP(C) + 2
                  K = NPTS(3) - KMAP(C) + 2

                  ! UPDATE SPECTRUM (ONLY ROOT PROCESS IN COMMUNICATOR "COMM")
                  IF (BOOL .AND. (KEY .EQ. -1)) SPEC(I, J, K) = CONJG(BUFFER(C))

                ENDDO

              ENDIF

              CALL MPI_COMM_FREE(COMM, IERR)

            ENDIF
            ! END EPSILON*/THETA*/PHI* REGION

            ! ENTER STATEMENT IF PROCESS "P" HAS POINTS IN THE GAMMA* REGION
            IF ( (GS(2, SLICE2GLOBAL(P)) .LT. NYQUIST(2)) .AND. (GE(3, SLICE2GLOBAL(P)) .GT. NYQUIST(3)) ) THEN

              ! CORRESPONDING MIN/MAX INDICES IN THE GAMMA REGION
              J0 = NPTS(2) - MIN(GE(2, SLICE2GLOBAL(P)), NYQUIST(2) - 1) + 2
              J1 = NPTS(2) - GS(2, SLICE2GLOBAL(P)) + 2

              K0 = NPTS(3) - GE(3, SLICE2GLOBAL(P)) + 2
              K1 = NPTS(3) - MAX(GS(3, SLICE2GLOBAL(P)), NYQUIST(3) + 1) + 2

              ! "TRUE" IF CALLING PROCESS HAS DATA NEEDED BY PROCESS "P"
              BOOL = (( (FS(2) .LE. J0) .AND. (J0 .LE. FE(2)) ) .OR. ( (FS(2) .LE. J1) .AND. (J1 .LE. FE(2)) ) .OR.      &
                      ( (FS(2) .GE. J0) .AND. (FE(2) .LE. J1) ))                                               .AND.     &
                     (( (FS(3) .LE. K0) .AND. (K0 .LE. FE(3)) ) .OR. ( (FS(3) .LE. K1) .AND. (K1 .LE. FE(3)) ) .OR.      &
                      ( (FS(3) .GE. K0) .AND. (FE(3) .LE. K1) ))

              COLOR = 0

              ! CALLING PROCESS IS PART OF NEW COMMUNICATOR
              IF (BOOL) THEN
                COLOR = 1
                KEY   = RANK
              ENDIF

              ! ADD PROCESS "P" TO COMMUNICATOR AS IT WILL COLLECT DATA
              IF (P .EQ. RANK) THEN
                COLOR = 1
                KEY   = -1                                                  !< MAKE SURE IT GETS LOWEST (0) RANK VALUE
              ENDIF

              ! SPLIT "SLICE" COMMUNICATOR
              CALL MPI_COMM_SPLIT(SLICE, COLOR, KEY, COMM, IERR)

              ! ONLY PROCESSES IN THE NEW COMMUNICATOR EXECUTE STATEMENTS BELOW
              IF (COLOR .EQ. 1) THEN

                ! GATHER DATA FROM PROCESSES
                CALL GATHER_DATA(COMM, I0, SPEC, JMAP, KMAP, BUFFER)

                DO C = 1, SIZE(BUFFER)

                  ! "TRUE" IF CURRENT POINT IS INSIDE THE GAMMA REGION
                  BOOL = (JMAP(C) .GE. J0) .AND. (JMAP(C) .LE. J1) .AND. (KMAP(C) .GE. K0) .AND. (KMAP(C) .LE. K1)

                  ! INDICES WHERE COMPLEX-CONJUGATED VALUES MUST BE PLACED
                  J = NPTS(2) - JMAP(C) + 2
                  K = NPTS(3) - KMAP(C) + 2

                  ! UPDATE SPECTRUM (ONLY ROOT PROCESS IN COMMUNICATOR "COMM")
                  IF (BOOL .AND. (KEY .EQ. -1)) SPEC(I, J, K) = CONJG(BUFFER(C))

                ENDDO

              ENDIF

              CALL MPI_COMM_FREE(COMM, IERR)

            ENDIF
            ! END GAMMA* REGION

            ! ENTER STATEMENT IF PROCESS "P" HAS POINTS IN THE BETA* REGION
            IF ( (GE(2, SLICE2GLOBAL(P)) .GT. NYQUIST(2)) .AND. (GS(3, SLICE2GLOBAL(P)) .EQ. 1) ) THEN

              ! CORRESPONDING MIN/MAX INDICES IN THE BETA REGION
              J0 = NPTS(2) - GE(2, SLICE2GLOBAL(P)) + 2
              J1 = NPTS(2) - MAX(GS(2, SLICE2GLOBAL(P)), NYQUIST(2) + 1) + 2

              K0 = 1
              K1 = 1

              ! "TRUE" IF CALLING PROCESS HAS DATA NEEDED BY PROCESS "P"
              BOOL = (( (FS(2) .LE. J0) .AND. (J0 .LE. FE(2)) ) .OR. ( (FS(2) .LE. J1) .AND. (J1 .LE. FE(2)) ) .OR.      &
                      ( (FS(2) .GE. J0) .AND. (FE(2) .LE. J1) ))                                               .AND.     &
                     (( (FS(3) .LE. K0) .AND. (K0 .LE. FE(3)) ) .OR. ( (FS(3) .LE. K1) .AND. (K1 .LE. FE(3)) ) .OR.      &
                      ( (FS(3) .GE. K0) .AND. (FE(3) .LE. K1) ))

              COLOR = 0

              ! CALLING PROCESS IS PART OF NEW COMMUNICATOR
              IF (BOOL) THEN
                COLOR = 1
                KEY   = RANK
              ENDIF

              ! ADD PROCESS "P" TO COMMUNICATOR AS IT WILL COLLECT DATA
              IF (P .EQ. RANK) THEN
                COLOR = 1
                KEY   = -1                                                  !< MAKE SURE IT GETS LOWEST (0) RANK VALUE
              ENDIF

              ! SPLIT "SLICE" COMMUNICATOR
              CALL MPI_COMM_SPLIT(SLICE, COLOR, KEY, COMM, IERR)

              ! ONLY PROCESSES IN THE NEW COMMUNICATOR EXECUTE STATEMENTS BELOW
              IF (COLOR .EQ. 1) THEN

                ! GATHER DATA FROM PROCESSES
                CALL GATHER_DATA(COMM, I0, SPEC, JMAP, KMAP, BUFFER)

                DO C = 1, SIZE(BUFFER)

                  ! "TRUE" IF CURRENT POINT IS INSIDE THE BETA REGION
                  BOOL = (JMAP(C) .GE. J0) .AND. (JMAP(C) .LE. J1) .AND. (KMAP(C) .GE. K0) .AND. (KMAP(C) .LE. K1)

                  ! INDICES WHERE COMPLEX-CONJUGATED VALUES MUST BE PLACED
                  J = NPTS(2) - JMAP(C) + 2
                  K = 1

                  ! UPDATE SPECTRUM (ONLY ROOT PROCESS IN COMMUNICATOR "COMM")
                  IF (BOOL .AND. (KEY .EQ. -1)) SPEC(I, J, K) = CONJG(BUFFER(C))

                ENDDO

              ENDIF

              CALL MPI_COMM_FREE(COMM, IERR)

            ENDIF
            ! END BETA* REGION

            ! ENTER STATEMENT IF PROCESS "P" HAS POINTS IN THE DELTA* REGION
            IF ( (GS(2, SLICE2GLOBAL(P)) .EQ. 1) .AND. (GE(3, SLICE2GLOBAL(P)) .GT. NYQUIST(3)) ) THEN

              ! CORRESPONDING MIN/MAX INDICES IN THE DELTA REGION
              J0 = 1
              J1 = 1

              K0 = NPTS(3) - GE(3, SLICE2GLOBAL(P)) + 2
              K1 = NPTS(3) - MAX(GS(3, SLICE2GLOBAL(P)), NYQUIST(3) + 1) + 2

              ! "TRUE" IF CALLING PROCESS HAS DATA NEEDED BY PROCESS "P"
              BOOL = (( (FS(2) .LE. J0) .AND. (J0 .LE. FE(2)) ) .OR. ( (FS(2) .LE. J1) .AND. (J1 .LE. FE(2)) ) .OR.      &
                      ( (FS(2) .GE. J0) .AND. (FE(2) .LE. J1) ))                                               .AND.     &
                     (( (FS(3) .LE. K0) .AND. (K0 .LE. FE(3)) ) .OR. ( (FS(3) .LE. K1) .AND. (K1 .LE. FE(3)) ) .OR.      &
                      ( (FS(3) .GE. K0) .AND. (FE(3) .LE. K1) ))

              COLOR = 0

              ! CALLING PROCESS IS PART OF NEW COMMUNICATOR
              IF (BOOL) THEN
                COLOR = 1
                KEY   = RANK
              ENDIF

              ! ADD PROCESS "P" TO COMMUNICATOR AS IT WILL COLLECT DATA
              IF (P .EQ. RANK) THEN
                COLOR = 1
                KEY   = -1                                                  !< MAKE SURE IT GETS LOWEST (0) RANK VALUE
              ENDIF

              ! SPLIT "SLICE" COMMUNICATOR
              CALL MPI_COMM_SPLIT(SLICE, COLOR, KEY, COMM, IERR)

              ! ONLY PROCESSES IN THE NEW COMMUNICATOR EXECUTE STATEMENTS BELOW
              IF (COLOR .EQ. 1) THEN

                ! GATHER DATA FROM PROCESSES
                CALL GATHER_DATA(COMM, I0, SPEC, JMAP, KMAP, BUFFER)

                DO C = 1, SIZE(BUFFER)

                  ! "TRUE" IF CURRENT POINT IS INSIDE THE BETA REGION
                  BOOL = (JMAP(C) .GE. J0) .AND. (JMAP(C) .LE. J1) .AND. (KMAP(C) .GE. K0) .AND. (KMAP(C) .LE. K1)

                  ! INDICES WHERE COMPLEX-CONJUGATED VALUES MUST BE PLACED
                  J = 1
                  K = NPTS(3) - KMAP(C) + 2

                  ! UPDATE SPECTRUM (ONLY ROOT PROCESS IN COMMUNICATOR "COMM")
                  IF (BOOL .AND. (KEY .EQ. -1)) SPEC(I, J, K) = CONJG(BUFFER(C))

                ENDDO

              ENDIF

              CALL MPI_COMM_FREE(COMM, IERR)

            ENDIF
            ! END BETA* REGION


          ENDDO
          ! END LOOP OVER PROCESSES IN "SLICE" COMMUNICATOR

          DEALLOCATE(SLICE2GLOBAL)

        ENDIF
        ! END IF FOR COLOR = 1

        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

        CALL MPI_COMM_FREE(SLICE, IERR)


      ENDDO
      ! END LOOP OVER INDEX "I"

    END SUBROUTINE ENFORCE_SYMMETRY

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE GATHER_DATA(COMM, I0, SPEC, JMAP, KMAP, BUFFER)

      ! GATHER DATA FROM ALL PROCESSES IN A COMMUNICATOR IN ROOT PROCESS. GLOBAL INDICES OF GATHERED DATA ARE RETURNED AS WELL.

      INTEGER(IPP),                                INTENT(IN)  :: COMM                             !< MPI COMMUNICATOR
      INTEGER(IPP),                                INTENT(IN)  :: I0                               !< LOCAL INDEX
      COMPLEX(FPP),              DIMENSION(:,:,:), INTENT(IN)  :: SPEC                             !< SPECTRUM
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:),     INTENT(OUT) :: JMAP, KMAP                       !< GLOBAL INDICES OF COLLECTED DATA
      COMPLEX(FPP), ALLOCATABLE, DIMENSION(:),     INTENT(OUT) :: BUFFER                           !< COLLECTED DATA
      COMPLEX(FPP), ALLOCATABLE, DIMENSION(:,:)                :: DUM                              !< LOCAL BUFFER
      INTEGER(IPP)                                             :: C, J, K, P                       !< COUNTERS
      INTEGER(IPP)                                             :: NPTS
      INTEGER(IPP)                                             :: NTASKS, RANK, IERR               !< MPI STUFF
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                  :: COMM2GLOBAL, RECVCOUNTS, DISPLS  !< MPI STUFF

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL MPI_COMM_SIZE(COMM, NTASKS, IERR)

      CALL MPI_COMM_RANK(COMM, RANK, IERR)

      ALLOCATE(COMM2GLOBAL(0:NTASKS-1))

      ! COMM-RANK TO GLOBAL-RANK MAPPING
      COMM2GLOBAL(RANK) = WORLD_RANK

      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_INTEGER, COMM2GLOBAL, 1, MPI_INTEGER, COMM, IERR)

      ALLOCATE(RECVCOUNTS(0:NTASKS-1), DISPLS(0:NTASKS-1))

      RECVCOUNTS(RANK) = (GE(2, WORLD_RANK) - GS(2, WORLD_RANK) + 1) * (GE(3, WORLD_RANK) - GS(3, WORLD_RANK) + 1)

      ! NUMBER OF POINTS TO BE RECEIVED BY EACH PROCESS IN "COMM" COMMUNICATOR
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_INTEGER, RECVCOUNTS, 1, MPI_INTEGER, COMM, IERR)

      ! TOTAL NUMBER OF POINTS TO BE RECEIVED
      NPTS = SUM(RECVCOUNTS)

      ALLOCATE(BUFFER(NPTS), JMAP(NPTS), KMAP(NPTS))

      ! DETERMINE THEIR RELATIVE DISPLACEMENT
      DISPLS(0) = 0

      DO P = 1, NTASKS - 1
        DISPLS(P) = DISPLS(P - 1) + RECVCOUNTS(P - 1)
      ENDDO

      ALLOCATE(DUM(SIZE(SPEC, 2), SIZE(SPEC, 3)))

      ! PREPARE DATA TO BE TRANSFERRED
      DO K = 1, SIZE(SPEC, 3)
        DO J = 1, SIZE(SPEC, 2)
          DUM(J, K) = SPEC(I0, J, K)
        ENDDO
      ENDDO

      ! COLLECT DATA INTO "BUFFER"
      CALL MPI_GATHERV(DUM, SIZE(DUM), COMPLEX_TYPE, BUFFER, RECVCOUNTS, DISPLS, COMPLEX_TYPE, 0, COMM, IERR)

      DEALLOCATE(DUM, DISPLS, RECVCOUNTS)

      C = 0

      ! CREATE MAPS FOR "J" AND "K" INDICES FOR ALL POINTS IN "BUFFER"
      DO P = 0, NTASKS - 1
        DO K = GS(3, COMM2GLOBAL(P)), GE(3, COMM2GLOBAL(P))
          DO J = GS(2, COMM2GLOBAL(P)), GE(2, COMM2GLOBAL(P))
            C       = C + 1
            JMAP(C) = J
            KMAP(C) = K
          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE(COMM2GLOBAL)

    END SUBROUTINE GATHER_DATA

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE EXCHANGE_HALO(CARTOPO, BUFFER)

      ! EXCHANGE HALO BETWEEN NEIGHBORING PROCESSES. THE HALO IS REPRESENTED BY THE FIRST COLUMN/ROW (ALONG X/Y, RESPECTIVELY) OF
      ! "BUFFER"

      INTEGER(IPP),                   INTENT(IN)    :: CARTOPO                         !< HANDLE TO CARTESIAN TOPOLOGY
      REAL(FPP),    DIMENSION(:,:,:), INTENT(INOUT) :: BUFFER                          !< RANDOM FIELD
      INTEGER(IPP)                                  :: I
      INTEGER(IPP)                                  :: NX, NY, NZ                       !< POINTS ALONG X-/Y-DIRECTION
      INTEGER(IPP)                                  :: IERR, NEG, POS
      INTEGER(IPP)                                  :: FROM_NEG, TO_POS
      INTEGER(IPP), DIMENSION(3)                    :: SIZES, SUBSIZES, RSTARTS, SSTARTS

      !-----------------------------------------------------------------------------------------------------------------------------

      ! POINTS ALONG X AND Y WITHOUT EXTRA ROW/COLUMN
      NX = SIZE(BUFFER, 1) - 1
      NY = SIZE(BUFFER, 2) - 1
      NZ = SIZE(BUFFER, 3) - 1

      ! TOTAL SIZE OF "BUFFER" WITH EXTRA ROW/COLUMN FOR HALO EXCHANGE
      SIZES = [NX + 1, NY + 1, NZ + 1]

      ! EXCHANGE IN THE X-, Y-, Z-DIRECTION
      DO I = 0, 2

        ! DETERMINE NEIGHBORS FOR A POSITIVE UNITARY SHIFT ALONG CURRENT DIRECTION
        CALL MPI_CART_SHIFT(CARTOPO, I, 1, NEG, POS, IERR)

        IF (I .EQ. 0) THEN
          SUBSIZES = [1, NY, NZ]            !< SHAPE OF HALO
          RSTARTS  = [0, 1, 1]              !< INITIAL ADDRESS FOR DATA TO BE RECEIVED
          SSTARTS  = [NX, 1, 1]             !< INITIAL ADDRESS FOR DATA TO BE SENT
        ELSEIF (I .EQ. 1) THEN
          SUBSIZES = [NX + 1, 1, 1]         !< SHAPE OF HALO
          RSTARTS  = [0, 0, 1]              !< INITIAL ADDRESS FOR DATA TO BE RECEIVED
          SSTARTS  = [0, NY, 1]             !< INITIAL ADDRESS FOR DATA TO BE SENT
        ELSE
          SUBSIZES = [NX + 1, NY + 1, 1]    !< SHAPE OF HALO
          RSTARTS  = [0, 0, 0]              !< INITIAL ADDRESS FOR DATA TO BE RECEIVED
          SSTARTS  = [0, 0, NZ]             !< INITIAL ADDRESS FOR DATA TO BE SENT
        ENDIF

        ! DATA TO BE RECEIVED
        CALL MPI_TYPE_CREATE_SUBARRAY(3, SIZES, SUBSIZES, RSTARTS, MPI_ORDER_FORTRAN, REAL_TYPE, FROM_NEG, IERR)
        CALL MPI_TYPE_COMMIT(FROM_NEG, IERR)

        ! DATA TO BE SENT
        CALL MPI_TYPE_CREATE_SUBARRAY(3, SIZES, SUBSIZES, SSTARTS, MPI_ORDER_FORTRAN, REAL_TYPE, TO_POS, IERR)
        CALL MPI_TYPE_COMMIT(TO_POS, IERR)

        ! EXCHANGE HALO DATA WITH NEIGHBORS. SINCE WE OPERATE ON FIRST AND LAST COLUMNS, SEND "BUFFER" IS DISJOINT FROM RECV "BUFFER"
        CALL MPI_SENDRECV(BUFFER, 1, TO_POS, POS, 0, BUFFER, 1, FROM_NEG, NEG, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, IERR)

        ! REALEASE RESOURCES
        CALL MPI_TYPE_FREE(FROM_NEG, IERR)
        CALL MPI_TYPE_FREE(TO_POS, IERR)

      ENDDO

    END SUBROUTINE EXCHANGE_HALO

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE INTERPOLATE(DH, X, Y, Z, DELTA, FIELD)

      !

      REAL(FPP),                      INTENT(IN)  :: DH                               !< GRID-STEP
      REAL(FPP),    DIMENSION(:),     INTENT(IN)  :: X, Y, Z                          !< POINTS LOCATION
      REAL(FPP),    DIMENSION(:,:,:), INTENT(IN)  :: DELTA                            !< RANDOM FIELD
      REAL(FPP),    DIMENSION(:),     INTENT(OUT) :: FIELD                            !< INTERPOLATED VALUES
      INTEGER(IPP)                                :: P                                !< COUNTERS
      INTEGER(IPP)                                :: IS, JS, KS, I0, J0, K0           !< INDICES
      INTEGER(IPP), DIMENSION(3)                  :: N
      LOGICAL ,     DIMENSION(3)                  :: BOOL
      REAL(FPP)                                   :: I, J, K
      REAL(FPP)                                   :: PX, PY, PZ, IPX, IPY             !< INTERPOLATION STUFF
      REAL(FPP)                                   :: A, B                             !< INTERPOLATION STUFF
      REAL(FPP),    DIMENSION(2,2)                :: F                                !< INTERPOLATION STUFF

      !-----------------------------------------------------------------------------------------------------------------------------

      IS = GS(1, WORLD_RANK)
      JS = GS(2, WORLD_RANK)
      KS = GS(3, WORLD_RANK)

      N(1) = SIZE(DELTA, 1)
      N(2) = SIZE(DELTA, 2)
      N(3) = SIZE(DELTA, 3)

      ! LOOP OVER POINTS
      DO P = 1, SIZE(X)

        I = (X(P) - OFF_AXIS(1)) / DH + 3._FPP - IS
        J = (Y(P) - OFF_AXIS(2)) / DH + 3._FPP - JS
        K = (Z(P) - OFF_AXIS(3)) / DH + 3._FPP - KS

        I0 = FLOOR(I)
        J0 = FLOOR(J)
        K0 = FLOOR(K)

        BOOL(1) = (I0 .LT. 1) .OR. (I0 .GE. N(1))
        BOOL(2) = (J0 .LT. 1) .OR. (J0 .GE. N(2))
        BOOL(3) = (K0 .LT. 1) .OR. (K0 .GE. N(3))

        IF (ANY(BOOL .EQV. .TRUE.)) CYCLE

        F(:, 1) = DELTA(I0:I0 + 1, J0, K0)
        F(:, 2) = DELTA(I0:I0 + 1, J0 + 1, K0)

        PX = I - I0
        PY = J - J0

        IPX = (1._FPP - PX)
        IPY = (1._FPP - PY)

        ! BILINEAR INTERPOLATION AT LEVEL 1
        A = F(1, 1) * IPX * IPY + F(2, 1) * PX * IPY + F(1, 2) * IPX * PY + F(2, 2) * PX * PY

        F(:, 1) = DELTA(I0:I0 + 1, J0, K0 + 1)
        F(:, 2) = DELTA(I0:I0 + 1, J0 + 1, K0 + 1)

        ! BILINEAR INTERPOLATION AT LEVEL 2
        B = F(1, 1) * IPX * IPY + F(2, 1) * PX * IPY + F(1, 2) * IPX * PY + F(2, 2) * PX * PY

        PZ = K - K0

        ! LINEAR INTERPOLATED BETWEEN LEVEL 1 AND 2
        FIELD(P) = A * (1._FPP - PZ) + B * PZ

      ENDDO

    END SUBROUTINE INTERPOLATE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE COORDS2INDEX(NPTS, DIMS, COORDS, FS, FE)

      ! MAP PROCESS COORDINATES INTO FIRST/LAST INDICES

      INTEGER(IPP), DIMENSION(3), INTENT(IN)  :: NPTS                        !< POINTS ALONG EACH AXIS
      INTEGER(IPP), DIMENSION(3), INTENT(IN)  :: DIMS                        !< NUMBER OF PROCESSES ALONG EACH AXIS
      INTEGER(IPP), DIMENSION(3), INTENT(IN)  :: COORDS                      !< CALLING PROCESS COORDINATES
      INTEGER(IPP), DIMENSION(3), INTENT(OUT) :: FS, FE                      !< RESULTING FIRST/LAST INDICES
      INTEGER(IPP)                            :: I                           !< COUNTER

      !-----------------------------------------------------------------------------------------------------------------------------

      DO I = 1, 3
        FS(I) = 1 + INT( REAL(NPTS(I)) / REAL(DIMS(I)) * REAL(COORDS(I)) )
        FE(I) = INT( REAL(NPTS(I)) / REAL(DIMS(I)) * REAL(COORDS(I) + 1) )
      ENDDO

    END SUBROUTINE COORDS2INDEX

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

    ! SUBROUTINE TEST_CONFIG(DIMS, NITER)
    !
    !   INTEGER(IPP),              DIMENSION(3),             INTENT(IN)  :: DIMS
    !   INTEGER(IPP),                                        INTENT(OUT) :: NITER
    !   INTEGER(IPP)                                                     :: I                            !< COUNTER
    !   INTEGER(IPP)                                                     :: IERR, CARTOPO
    !   INTEGER(IPP),              DIMENSION(3)                          :: LS, LE
    !   INTEGER(IPP),              DIMENSION(3)                          :: COORDS
    !   INTEGER(IPP),              DIMENSION(0:WORLD_SIZE-1)             :: COMPLETED, ISBUSY
    !   INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                          :: BUFFER
    !
    !   !-----------------------------------------------------------------------------------------------------------------------------
    !
    !   ! CREATE TOPOLOGY
    !   CALL MPI_CART_CREATE(MPI_COMM_WORLD, NDIMS, DIMS, ISPERIODIC, REORDER, CARTOPO, IERR)
    !
    !   ! DO I = 0, WORLD_SIZE - 1
    !   !   CALL MPI_CART_COORDS(CARTOPO, I, NDIMS, COORDS, IERR)
    !   !   CALL COORDS2INDEX(NPTS, DIMS, COORDS, LS, LE)
    !   !   GS(:, I) = LS
    !   !   GE(:, I) = LE
    !   ! ENDDO
    !
    !   ! RETURN PROCESS COORDINATES IN CURRENT TOPOLOGY
    !   CALL MPI_CART_COORDS(CARTOPO, WORLD_RANK, NDIMS, COORDS, IERR)
    !
    !   ! RETURN FIRST/LAST-INDEX ("LS"/"LE") ALONG EACH DIRECTION FOR CALLING PROCESS. NOTE: FIRST POINT HAS ALWAYS INDEX EQUAL TO 1.
    !   CALL COORDS2INDEX(NPTS, DIMS, COORDS, LS, LE)
    !
    !   GS(:, WORLD_RANK) = LS
    !   GE(:, WORLD_RANK) = LE
    !
    !   ! MAKE ALL PROCESSES AWARE OF GLOBAL INDICES ALONG EACH AXIS
    !   CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, GS, 3, MPI_INTEGER, MPI_COMM_WORLD, IERR)
    !   CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, GE, 3, MPI_INTEGER, MPI_COMM_WORLD, IERR)
    !
    !   ! INITIALISE LIST WITH PROCESSES HAVING THEIR POINTS ASSIGNED
    !   COMPLETED(:) = 0
    !
    !   NITER = 0
    !
    !   ! CYCLE UNTIL ALL PROCESSES ARE COMPLETED
    !   DO WHILE(ANY(COMPLETED .EQ. 0))
    !
    !     ! RESET "BUSY PROCESS" FLAG
    !     ISBUSY(:) = 0
    !
    !     ! CREATE AS MANY COMMUNICATORS AS POSSIBLE
    !     DO I = 0, WORLD_SIZE - 1
    !
    !       IF (COMPLETED(I) .EQ. 0) THEN
    !
    !         ! PRODUCE A TENTATIVE LIST OF PROCESSES THAT MAY JOIN THE NEW COMMUNICATOR
    !         CALL FIND_BUDDIES(I, BUFFER)
    !
    !         ! BUILD NEW COMMUNICATOR ONLY IF ALL INVOLVED PROCESSES ARE NOT YET PART OF ANOTHER COMMUNICATOR
    !         IF (ALL(ISBUSY(BUFFER) .EQ. 0)) THEN
    !
    !           ! SET PROCESSES BELONGING TO NEW COMMUNICATOR AS BUSY
    !           ISBUSY(BUFFER) = 1
    !
    !           ! SET I-TH PROCESS AS COMPLETED
    !           COMPLETED(I) = 1
    !
    !         ENDIF
    !
    !         DEALLOCATE(BUFFER)
    !
    !       ENDIF
    !
    !     ENDDO
    !
    !     NITER = NITER + 1
    !
    !   ENDDO
    !
    !   CALL MPI_COMM_FREE(CARTOPO, IERR)
    !
    ! END SUBROUTINE TEST_CONFIG

    SUBROUTINE TEST_CONFIG(DIMS, MEASURE)

      INTEGER(IPP),              DIMENSION(3),             INTENT(IN)  :: DIMS
      REAL(FPP),                                           INTENT(OUT) :: MEASURE
      INTEGER(IPP)                                                     :: I                            !< COUNTER
      INTEGER(IPP)                                                     :: IERR, CARTOPO
      INTEGER(IPP),              DIMENSION(3)                          :: LS, LE
      INTEGER(IPP),              DIMENSION(3)                          :: COORDS
      INTEGER(IPP),              DIMENSION(0:WORLD_SIZE-1)             :: COMPLETED, ISBUSY
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                          :: BUFFER

      REAL(FPP), DIMENSION(3) :: SIDE

      !-----------------------------------------------------------------------------------------------------------------------------

      SIDE = REAL(NPTS, FPP) / REAL(DIMS, FPP)

      MEASURE = ABS(SIDE(1) - SIDE(2)) + ABS(SIDE(1) - SIDE(3)) + ABS(SIDE(2) - SIDE(3))


    END SUBROUTINE TEST_CONFIG


    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE FACTORIZATION(N, FACTORS)

      ! FACTORISE INTEGER "N" BASED ON TRIAL DIVISION METHOD

      INTEGER(IPP),                              INTENT(IN)  :: N
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: FACTORS
      INTEGER(IPP)                                           :: I, C, S
      INTEGER(IPP)                                           :: X
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:,:)              :: BUFFER

      !-----------------------------------------------------------------------------------------------------------------------------

      ! MAX POSSIBLE FACTOR
      S = FLOOR(SQRT(REAL(N, FPP)))

      ALLOCATE(BUFFER(2, S))

      BUFFER(:,:) = 0

      ! TEST FACTORS
      DO I = 1, S

        X = N / I

        IF (MOD(N, I) .EQ. 0) THEN
          BUFFER(1, I) = I                          !< ADD FACTOR ...
          BUFFER(2, I) = X                          !< ... AND ITS COMPANION
        ENDIF

      ENDDO

      ! ACTUAL FACTORS FOUND
      I = COUNT(BUFFER(1, :) .NE. 0)

      ALLOCATE(FACTORS(2, I))

      ! COPY FACTORS TO OUTPUT ARRAY
      C = 0
      DO I = 1, S
        IF (BUFFER(1, I) .NE. 0) THEN
          C = C + 1
          FACTORS(:, C) = BUFFER(:, I)
        ENDIF
      ENDDO

      DEALLOCATE(BUFFER)

    END SUBROUTINE FACTORIZATION

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE SCARFLIB_FFT3
