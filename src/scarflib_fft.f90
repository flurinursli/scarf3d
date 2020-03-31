MODULE SCARFLIB_FFT

  ! ALL VARIABLES AND SUBMODULE PROCEDURES ARE GLOBAL WITHIN THE SUBMODULE, BUT LIMITED TO IT.

  ! mpif90 -O3 scarflib_fft2.f90 scarflib_spec.f90 scarflib.f90 driver.f90 *.f -I/home/walter/Backedup/Software/fftw-3.3.8/include
  ! -L/home/walter/Backedup/Software/fftw-3.3.8/lib -lfftw3 -fcheck=all -fbacktrace -fopenacc

  USE, NON_INTRINSIC :: SCARFLIB_COMMON

  IMPLICIT NONE

  INCLUDE 'fftw3.f03'

  PRIVATE

  PUBLIC :: SCARF3D_FFT !, BEST_CONFIG

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTERFACE SCARF3D_FFT
    MODULE PROCEDURE SCARF3D_UNSTRUCTURED_FFT, SCARF3D_STRUCTURED_FFT
  END INTERFACE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! PROCEDURE POINTERS
  PROCEDURE(VK_PSDF), POINTER :: FUN

  ! INTERFACE TO C++ FUNCTIONS/SUBROUTINES
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
  !INTEGER(IPP), ALLOCATABLE, DIMENSION(:,:)           :: PS, PE

  ! MASK TO DETERMINE IF I-TH PROCESS HAS POINTS FALLING INSIDE DOMAIN OF CALLING PROCESS (0=NO, 1=YES)
  INTEGER(IPP), ALLOCATABLE, DIMENSION(:)             :: POINTS_WORLD2RANK

  ! MASK TO DETERMINE IF CALLING PROCESS HAS POINTS FALLING INSIDE DOMAIN OF I-TH PROCESS (0=NO, 1=YES)
  INTEGER(IPP), ALLOCATABLE, DIMENSION(:)             :: POINTS_RANK2WORLD

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
  TYPE(C_PTR)                                         :: PC

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE SCARF3D_UNSTRUCTURED_FFT(X, Y, Z, DH, ACF, CL, SIGMA, HURST, SEED, POI, MUTE, TAPER, RESCALE, PAD, FIELD, INFO)

      ! GLOBAL INDICES AND COMMUNICATOR SUBGROUPPING COULD BE HANDLED ELEGANTLY IF VIRTUAL TOPOLOGIES ARE USED IN CALLING PROGRAM.
      ! HOWEVER WE ASSUME THAT THESE ARE NOT USED AND THEREFORE WE ADOPT A SIMPLER APPROACH BASED ON COLLECTIVE CALLS.

      ! "SEED" AND NUMBER OF POINTS "NPTS" SHOULD BE THE SAME TO GUARANTEE SAME RANDOM FIELD

      REAL(FPP),                     DIMENSION(:),     INTENT(IN)  :: X, Y, Z              !< POSITION OF POINTS ALONG X, Y, Z
      REAL(FPP),                                       INTENT(IN)  :: DH                   !< GRID-STEP
      INTEGER(IPP),                                    INTENT(IN)  :: ACF                  !< AUTOCORRELATION FUNCTION: 0=VK, 1=GAUSS
      REAL(FPP),                     DIMENSION(3),     INTENT(IN)  :: CL                   !< CORRELATION LENGTH
      REAL(FPP),                                       INTENT(IN)  :: SIGMA                !< STANDARD DEVIATION
      REAL(FPP),                                       INTENT(IN)  :: HURST                !< HURST EXPONENT
      INTEGER(IPP),                                    INTENT(IN)  :: SEED                 !< SEED NUMBER
      REAL(FPP),                     DIMENSION(:,:),   INTENT(IN)  :: POI                  !< LOCATION OF POINT(S)-OF-INTEREST
      REAL(FPP),                                       INTENT(IN)  :: MUTE                 !< NUMBER OF POINTS WHERE MUTING IS APPLIED
      REAL(FPP),                                       INTENT(IN)  :: TAPER                !< NUMBER OF POINTS WHERE TAPERING IS APPLIED
      INTEGER(IPP),                                    INTENT(IN)  :: RESCALE              !< FLAG FOR RESCALING RANDOM FIELD TO DESIRED SIGMA
      INTEGER(IPP),                                    INTENT(IN)  :: PAD                  !< FLAG FOR HANDLE FFT PERIODICITY
      REAL(FPP),                     DIMENSION(:),     INTENT(OUT) :: FIELD                !< RANDOM FIELD AT X,Y,Z LOCATION
      REAL(FPP),                     DIMENSION(8),     INTENT(OUT) :: INFO                 !< ERRORS AND TIMING FOR PERFORMANCE ANALYSIS
      COMPLEX(FPP),     ALLOCATABLE, DIMENSION(:,:,:)              :: SPEC                 !< SPECTRUM/RANDOM FIELD
      INTEGER(IPP)                                                 :: IERR                 !< MPI STUFF
      INTEGER(IPP)                                                 :: I, J, K              !< COUNTERS
      INTEGER(IPP)                                                 :: CARTOPO
      INTEGER(IPP)                                                 :: OFFSET               !< EXTRA POINTS ON EACH SIDE FOR FFT PERIODICITY
      INTEGER(IPP),                  DIMENSION(3)                  :: M                    !< POINTS FOR CALLING PROCESS
      INTEGER(IPP),                  DIMENSION(3)                  :: LS, LE               !< FIRST/LAST INDEX ALONG X, Y, Z
      INTEGER(IPP),                  DIMENSION(3)                  :: COORDS
      INTEGER(IPP),                  DIMENSION(3)                  :: DIMS
      INTEGER(IPP),                  DIMENSION(3)                  :: PS, PE
      INTEGER(IPP),     ALLOCATABLE, DIMENSION(:)                  :: NPOINTS              !< NUMBER OF INPUT POINTS FOR EACH PROCESS
      LOGICAL,                       DIMENSION(3)                  :: BOOL
      REAL(FPP)                                                    :: SCALING              !< SCALING FACTOR
      REAL(REAL64)                                                 :: TICTOC               !< TIMING
      REAL(FPP),                     DIMENSION(2)                  :: ET                   !< DUMMY FOR ELAPSED TIME
      REAL(FPP),                     DIMENSION(3)                  :: MIN_EXTENT           !< LOWER MODEL LIMITS (GLOBAL)
      REAL(FPP),                     DIMENSION(3)                  :: MAX_EXTENT           !< UPPER MODEL LIMITS (GLOBAL)
      REAL(FPP),                     DIMENSION(3)                  :: PT
      REAL(FPP),        ALLOCATABLE, DIMENSION(:)                  :: VAR, MU              !< STATISTICS: VARIANCE AND AVERAGE
      REAL(FPP),        ALLOCATABLE, DIMENSION(:,:,:)              :: DELTA, BUFFER        !< RANDOM PERTURBATIONS ON COMPUTATIONAL GRID

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

      ! CYCLE OVER THE THREE MAIN DIRECTIONS TO DETERMINE THE NECESSARY MODEL SIZE (IN NUMBER OF POINTS) AND THE ABSOLUTE POSITION OF
      ! OF THE FIRST POINT TO COUNTERACT FFT PERIODICITY
      DO I = 1, 3

        ! MINIMUM NUMBER OF POINTS SUCH THAT RANDOM FIELD COVERS THE WHOLE DOMAIN. IN PRACTICE, WE DEFINE THE RANDOM FIELD ON A GRID
        ! SLIGHTLY LARGER THAN NECESSARY (HALF GRID-STEP) IN EACH DIRECTION, IN ORDER TO AVOID SIDE-EFFECTS DUE TO INTERPOLATION WHEN
        ! THE EXTRA EXTENSION TO HANDLE FFT PERIODICITY IS NOT DESIRED.
        NPTS(I) = NINT( (MAX_EXTENT(I) + DH * 0.5_FPP - MIN_EXTENT(I) + DH * 0.5_FPP) / DH) + 1

        ! POINTS FOR ONE CORRELATION LENGTH
        OFFSET = NINT(CL(I) / DH)

        ! DO NOT EXTEND MODEL UNLESS DESIRED
        IF (PAD .NE. 1) OFFSET = 0

        ! EXTEND THE MODEL BY AT LEAST ONE CORRELATION LENGTH (ON EACH SIDE) TO COUNTERACT FFT PERIODICITY
        NPTS(I) = NPTS(I) + 2 * OFFSET

        ! MAKE SURE WE HAVE EVEN NUMBER OF POINTS
        IF (MOD(NPTS(I), 2) .NE. 0) NPTS(I) = NPTS(I) + 1

        ! ABSOLUTE POSITION OF FIRST POINT (IT COULD BE NEGATIVE)
        ! OFF_AXIS(I) = MIN_EXTENT(I) - OFFSET * DH
        OFF_AXIS(I) = MIN_EXTENT(I) - (OFFSET + 0.5_FPP) * DH

      ENDDO

      INFO(:) = 0._FPP

      ! HERE WE CHECK IF THE MODEL IS LARGE ENOUGH TO CATCH THE LOWER PART OF THE SPECTRUM (LOW WAVENUMBERS)...
      IF (ANY(NPTS .LE. NINT(2._FPP * PI * CL / DH))) INFO(1) = 1._FPP

      ! ...AND HERE IF THE GRID-STEP IS SMALL ENOUGH TO CATCH THE UPPER PART OF THE SPECTRUM (HIGH WAVENUMBERS)
      IF (ANY(DH .GT. CL / 2._FPP)) INFO(2) = 1._FPP

      ! [ANOTHER WAY TO CHECK HOW CLOSE THE DISCRETE AND CONTINUOUS SPECTRA ARE, IS TO COMPARE THE RESPECTIVE STANDARD DEVIATIONS]

      ! HERE BELOW WE CREATE A REGULAR MESH AND A CARTESIAN TOPOLOGY. THE RANDOM FIELD WILL BE CALCULATED ON THIS MESH OF "NPTS" POINTS
      ! AND THEN INTERPOLATED AT THOSE POINTS (POSSIBLY IRREGULARLY DISTRIBUTED) OWNED BY EACH SINGLE PROCESS

      ! "GS"/"GE" STORE FIRST/LAST GLOBAL INDICES ALONG EACH DIRECTION FOR ALL PROCESS
      ALLOCATE(GS(3, 0:WORLD_SIZE-1), GE(3, 0:WORLD_SIZE-1))

      ! RETURN PROCESSORS GRID RESULTING MAXIMIZING THE NUMBER OF OVERLAPPING CALLS FOR INTERPOLATION (SEE BELOW)
      CALL BEST_CONFIG(DIMS)

      ! TEMPORARY MESSAGE
      !IF (WORLD_RANK == 0) PRINT*, 'BEST PROCESSORS GRID: ', DIMS

      ! CREATE TOPOLOGY
      CALL MPI_CART_CREATE(MPI_COMM_WORLD, NDIMS, DIMS, ISPERIODIC, REORDER, CARTOPO, IERR)

      ! RETURN PROCESS COORDINATES IN CURRENT TOPOLOGY
      CALL MPI_CART_COORDS(CARTOPO, WORLD_RANK, NDIMS, COORDS, IERR)

      ! RETURN FIRST/LAST-INDEX ("LS"/"LE") ALONG EACH DIRECTION FOR CALLING PROCESS. NOTE: FIRST POINT HAS ALWAYS INDEX EQUAL TO 1.
      CALL COORDS2INDEX(NPTS, DIMS, COORDS, LS, LE)

      GS(:, WORLD_RANK) = LS
      GE(:, WORLD_RANK) = LE

      ! MAKE ALL PROCESSES AWARE OF GLOBAL INDICES ALONG EACH AXIS
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, GS, 3, MPI_INTEGER, MPI_COMM_WORLD, IERR)
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, GE, 3, MPI_INTEGER, MPI_COMM_WORLD, IERR)

      ALLOCATE(POINTS_WORLD2RANK(0:WORLD_SIZE - 1), POINTS_RANK2WORLD(0:WORLD_SIZE - 1))

      ! MIN/MAX COMPUTATIONAL GRID INDICES CONTAINING INPUT POINTS
      PS(1) = FLOOR((MINVAL(X, DIM = 1) - OFF_AXIS(1)) / DH) + 1
      PE(1) = FLOOR((MAXVAL(X, DIM = 1) - OFF_AXIS(1)) / DH) + 1
      PS(2) = FLOOR((MINVAL(Y, DIM = 1) - OFF_AXIS(2)) / DH) + 1
      PE(2) = FLOOR((MAXVAL(Y, DIM = 1) - OFF_AXIS(2)) / DH) + 1
      PS(3) = FLOOR((MINVAL(Z, DIM = 1) - OFF_AXIS(3)) / DH) + 1
      PE(3) = FLOOR((MAXVAL(Z, DIM = 1) - OFF_AXIS(3)) / DH) + 1

      POINTS_RANK2WORLD = 0

      ! DETERMINE IF CALLING PROCESS HAS AT LEAST ONE POINT FALLING INSIDE DOMAIN OF "I-TH" PROCESS
      DO I = 0, WORLD_SIZE - 1
        DO J = 1, 3
          BOOL(J) = ( (PS(J) .LT. GE(J, I)) .AND. (PS(J) .GE. (GS(J, I) - 1)) ) .OR.     &
                    ( (PE(J) .LT. GE(J, I)) .AND. (PE(J) .GE. (GS(J, I) - 1)) ) .OR.     &
                    ( (PE(J) .GE. GE(J, I)) .AND. (PS(J) .LT. (GS(J, I) - 1)) )
        ENDDO
        IF (ALL(BOOL .EQV. .TRUE.)) POINTS_RANK2WORLD(I) = 1
      ENDDO

      ! DETERMINE IF PROCESS "I" HAS AT LEAST ONE POINT FALLING INSIDE DOMAIN OF CALLING PROCESS
      CALL MPI_ALLTOALL(POINTS_RANK2WORLD, 1, MPI_INTEGER, POINTS_WORLD2RANK, 1, MPI_INTEGER, MPI_COMM_WORLD, IERR)

      ! NUMBER OF POINTS FOR CALLING PROCESS
      DO I = 1, 3
        M(I) = LE(I) - LS(I) + 1
      ENDDO

      ! ALLOCATE MEMORY FOR SPECTRUM
      ALLOCATE(SPEC(M(1), M(2), M(3)))

      ! COMPUTE SPECTRUM AND APPLY HERMITIAN SYMMETRY
      CALL COMPUTE_SPECTRUM(DH, LS, LE, DH, ACF, CL, SIGMA, HURST, SEED, SPEC, ET)

#ifdef TIMING
      INFO(5:6) = ET

      ! START TIMER
      CALL WATCH_START(TICTOC)
#endif

      ! TRANSFORM ALONG EACH DIRECTION
      CALL TRANSFORM_ALONG_Z(SPEC)
      CALL TRANSFORM_ALONG_Y(SPEC)
      CALL TRANSFORM_ALONG_X(SPEC)

#ifdef TIMING
      CALL WATCH_STOP(TICTOC)

      INFO(7) = TICTOC
#endif

      ! SCALING PARAMETER
      SCALING = 1._FPP / SQRT(REAL(NPTS(1), FPP) * REAL(NPTS(2), FPP) * REAL(NPTS(3), FPP) * DH**3)

      ALLOCATE(BUFFER(M(1), M(2), M(3)))

      ! SCALE IFFT
      DO K = 1, M(3)
        DO J = 1, M(2)
          DO I = 1, M(1)
            BUFFER(I, J, K) = REAL(SPEC(I, J, K), FPP) * SCALING
          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE(SPEC)

      ALLOCATE(VAR(0:WORLD_SIZE - 1), MU(0:WORLD_SIZE - 1), NPOINTS(0:WORLD_SIZE - 1))

      ! COMPUTE VARIANCE AND MEAN OF RANDOM FIELD FOR EACH SINGLE PROCESS
      VAR(WORLD_RANK) = VARIANCE(BUFFER)
      MU(WORLD_RANK)  = MEAN(BUFFER)

      ! NUMBER OF GRID POINTS PER PROCESS
      NPOINTS(WORLD_RANK) = PRODUCT(M)

      ! SHARE RESULTS
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, VAR, 1, REAL_TYPE, MPI_COMM_WORLD, IERR)
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, MU, 1, REAL_TYPE, MPI_COMM_WORLD, IERR)
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, NPOINTS, 1, MPI_INTEGER, MPI_COMM_WORLD, IERR)

      ! COMPUTE TOTAL VARIANCE ("INFO(3)") AND MEAN ("INFO(4)")
      CALL PARALLEL_VARIANCE(VAR, MU, NPOINTS, INFO(3), INFO(4))

      ! RETURN STANDARD DEVIATION
      INFO(3) = SQRT(INFO(3))

      ! DO WE NEED TO RESCALE THE RANDOM FIELD TO DESIRED (I.E. CONTINUOUS) STANDARD DEVIATION?
      IF (RESCALE .EQ. 1) THEN

        SCALING = SIGMA / INFO(3)

        DO K = 1, M(3)
          DO J = 1, M(2)
            DO I = 1, M(1)
              BUFFER(I, J, K) = BUFFER(I, J, K) * SCALING - INFO(4)
            ENDDO
          ENDDO
        ENDDO

      ENDIF

      DEALLOCATE(VAR, MU, NPOINTS)

      ! APPLY TAPER/MUTE
      DO I = 1, SIZE(POI, 2)
        PT(:) = POI(:, I) - OFF_AXIS(:)
        CALL TAPERING(DH, LS, LE, BUFFER, PT, MUTE, TAPER)
      ENDDO

      ! COPY RANDOM FIELD TO ARRAY WITH HALO
      ALLOCATE(DELTA(0:M(1), 0:M(2), 0:M(3)))

      DELTA = 0._FPP

      ! SCALE IFFT
      DO K = 1, M(3)
        DO J = 1, M(2)
          DO I = 1, M(1)
            DELTA(I, J, K) = BUFFER(I, J, K)
          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE(BUFFER)

#ifdef TIMING
      CALL WATCH_START(TICTOC)
#endif

      ! EXCHANGE HALO
      CALL EXCHANGE_HALO(CARTOPO, DELTA)

      ! INTERPOLATE RANDOM FIELD VALUES AT DESIRED OUTPUT LOCATIONS
      CALL GRID_MAPPING(DH, DELTA, X, Y, Z, FIELD)

#ifdef TIMING
      CALL WATCH_STOP(TICTOC)

      INFO(8) = TICTOC
#endif

      ! RELEASE MEMORY
      DEALLOCATE(DELTA)
      DEALLOCATE(GS, GE)
      DEALLOCATE(POINTS_RANK2WORLD, POINTS_WORLD2RANK)

      CALL MPI_COMM_FREE(CARTOPO, IERR)

      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

    END SUBROUTINE SCARF3D_UNSTRUCTURED_FFT

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE SCARF3D_STRUCTURED_FFT(DS, FS, FE, DH, ACF, CL, SIGMA, HURST, SEED, POI, MUTE, TAPER, RESCALE, PAD, FIELD, INFO)

      ! GLOBAL INDICES AND COMMUNICATOR SUBGROUPPING COULD BE HANDLED ELEGANTLY IF VIRTUAL TOPOLOGIES ARE USED IN CALLING PROGRAM.
      ! HOWEVER WE ASSUME THAT THESE ARE NOT USED AND THEREFORE WE ADOPT A SIMPLER APPROACH BASED ON COLLECTIVE CALLS.

      ! "SEED" AND NUMBER OF POINTS "NPTS" SHOULD BE THE SAME TO GUARANTEE SAME RANDOM FIELD

      REAL(FPP),                                       INTENT(IN)  :: DS                   !< STRUCTURED GRID-STEP
      INTEGER(IPP),                  DIMENSION(3),     INTENT(IN)  :: FS, FE               !< FIRST/LAST STRUCTURED GRID INDICES
      REAL(FPP),                                       INTENT(IN)  :: DH                   !< GRID-STEP
      INTEGER(IPP),                                    INTENT(IN)  :: ACF                  !< AUTOCORRELATION FUNCTION: 0=VK, 1=GAUSS
      REAL(FPP),                     DIMENSION(3),     INTENT(IN)  :: CL                   !< CORRELATION LENGTH
      REAL(FPP),                                       INTENT(IN)  :: SIGMA                !< STANDARD DEVIATION
      REAL(FPP),                                       INTENT(IN)  :: HURST                !< HURST EXPONENT
      INTEGER(IPP),                                    INTENT(IN)  :: SEED                 !< SEED NUMBER
      REAL(FPP),                     DIMENSION(:,:),   INTENT(IN)  :: POI                  !< LOCATION OF POINT(S)-OF-INTEREST
      REAL(FPP),                                       INTENT(IN)  :: MUTE                 !< NUMBER OF POINTS WHERE MUTING IS APPLIED
      REAL(FPP),                                       INTENT(IN)  :: TAPER                !< NUMBER OF POINTS WHERE TAPERING IS APPLIED
      INTEGER(IPP),                                    INTENT(IN)  :: RESCALE              !< FLAG FOR RESCALING RANDOM FIELD TO DESIRED SIGMA
      INTEGER(IPP),                                    INTENT(IN)  :: PAD                  !< FLAG FOR HANDLE FFT PERIODICITY
      REAL(FPP),                     DIMENSION(:,:,:), INTENT(OUT) :: FIELD                !< RANDOM FIELD AT X,Y,Z LOCATION
      REAL(FPP),                     DIMENSION(8),     INTENT(OUT) :: INFO                 !< ERRORS AND TIMING FOR PERFORMANCE ANALYSIS
      COMPLEX(FPP),     ALLOCATABLE, DIMENSION(:,:,:)              :: SPEC                 !< SPECTRUM/RANDOM FIELD
      INTEGER(IPP)                                                 :: IERR                 !< MPI STUFF
      INTEGER(IPP)                                                 :: I, J, K              !< COUNTERS
      INTEGER(IPP)                                                 :: CARTOPO
      INTEGER(IPP)                                                 :: OFFSET               !< EXTRA POINTS ON EACH SIDE FOR FFT PERIODICITY
      INTEGER(IPP),                  DIMENSION(3)                  :: M                    !< POINTS FOR CALLING PROCESS
      INTEGER(IPP),                  DIMENSION(3)                  :: LS, LE               !< FIRST/LAST INDEX ALONG X, Y, Z
      INTEGER(IPP),                  DIMENSION(3)                  :: COORDS
      INTEGER(IPP),                  DIMENSION(3)                  :: DIMS
      INTEGER(IPP),                  DIMENSION(3)                  :: PS, PE
      INTEGER(IPP),     ALLOCATABLE, DIMENSION(:)                  :: NPOINTS              !< NUMBER OF INPUT POINTS FOR EACH PROCESS
      LOGICAL,                       DIMENSION(3)                  :: BOOL
      REAL(FPP)                                                    :: SCALING              !< SCALING FACTOR
      REAL(REAL64)                                                 :: TICTOC               !< TIMING
      REAL(FPP),                     DIMENSION(2)                  :: ET                   !< DUMMY FOR ELAPSED TIME
      REAL(FPP),                     DIMENSION(3)                  :: MIN_EXTENT           !< LOWER MODEL LIMITS (GLOBAL)
      REAL(FPP),                     DIMENSION(3)                  :: MAX_EXTENT           !< UPPER MODEL LIMITS (GLOBAL)
      REAL(FPP),                     DIMENSION(3)                  :: PT
      REAL(FPP),        ALLOCATABLE, DIMENSION(:)                  :: VAR, MU              !< STATISTICS: VARIANCE AND AVERAGE
      REAL(FPP),        ALLOCATABLE, DIMENSION(:,:,:)              :: DELTA, BUFFER        !< RANDOM PERTURBATIONS ON COMPUTATIONAL GRID

      !-----------------------------------------------------------------------------------------------------------------------------

      ! GET RANK NUMBER
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, WORLD_RANK, IERR)

      ! GET AVAILABLE MPI PROCESSES
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, WORLD_SIZE, IERR)

      ! MODEL LIMITS (PROCESS-WISE) ALONG EACH AXIS
      MIN_EXTENT = (FS - 1) * DS
      MAX_EXTENT = (FE - 1) * DS

      ! (GLOBAL) MODEL LIMITS
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, MIN_EXTENT, 3, REAL_TYPE, MPI_MIN, MPI_COMM_WORLD, IERR)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, MAX_EXTENT, 3, REAL_TYPE, MPI_MAX, MPI_COMM_WORLD, IERR)

      ! CYCLE OVER THE THREE MAIN DIRECTIONS TO DETERMINE THE NECESSARY MODEL SIZE (IN NUMBER OF POINTS) AND THE ABSOLUTE POSITION OF
      ! OF THE FIRST POINT TO COUNTERACT FFT PERIODICITY
      DO I = 1, 3

        ! MINIMUM NUMBER OF POINTS SUCH THAT RANDOM FIELD COVERS THE WHOLE DOMAIN. IN PRACTICE, WE DEFINE THE RANDOM FIELD ON A GRID
        ! SLIGHTLY LARGER THAN NECESSARY (HALF GRID-STEP) IN EACH DIRECTION, IN ORDER TO AVOID SIDE-EFFECTS DUE TO INTERPOLATION WHEN
        ! THE EXTRA EXTENSION TO HANDLE FFT PERIODICITY IS NOT DESIRED.
        NPTS(I) = NINT( (MAX_EXTENT(I) + DH * 0.5_FPP - MIN_EXTENT(I) + DH * 0.5_FPP) / DH) + 1

print*, 'npts ', npts(i), ' - ', MAX_EXTENT, ' ', MIN_EXTENT

        ! POINTS FOR ONE CORRELATION LENGTH
        OFFSET = NINT(CL(I) / DH)

        ! DO NOT EXTEND MODEL UNLESS DESIRED
        IF (PAD .NE. 1) OFFSET = 0

        ! EXTEND THE MODEL BY AT LEAST ONE CORRELATION LENGTH (ON EACH SIDE) TO COUNTERACT FFT PERIODICITY
        NPTS(I) = NPTS(I) + 2 * OFFSET

        ! MAKE SURE WE HAVE EVEN NUMBER OF POINTS
        IF (MOD(NPTS(I), 2) .NE. 0) NPTS(I) = NPTS(I) + 1

        ! ABSOLUTE POSITION OF FIRST POINT (IT COULD BE NEGATIVE)
        !OFF_AXIS(I) = MIN_EXTENT(I) - OFFSET * DH
        OFF_AXIS(I) = MIN_EXTENT(I) - (OFFSET + 0.5_FPP) * DH

      ENDDO

      INFO(:) = 0._FPP

      ! HERE WE CHECK IF THE MODEL IS LARGE ENOUGH TO CATCH THE LOWER PART OF THE SPECTRUM (LOW WAVENUMBERS)...
      IF (ANY(NPTS .LE. NINT(2._FPP * PI * CL / DH))) INFO(1) = 1._FPP

      ! ...AND HERE IF THE GRID-STEP IS SMALL ENOUGH TO CATCH THE UPPER PART OF THE SPECTRUM (HIGH WAVENUMBERS)
      IF (ANY(DH .GT. CL / 2._FPP)) INFO(2) = 1._FPP

      ! [ANOTHER WAY TO CHECK HOW CLOSE THE DISCRETE AND CONTINUOUS SPECTRA ARE, IS TO COMPARE THE RESPECTIVE STANDARD DEVIATIONS]

      ! HERE BELOW WE CREATE A REGULAR MESH AND A CARTESIAN TOPOLOGY. THE RANDOM FIELD WILL BE CALCULATED ON THIS MESH OF "NPTS" POINTS
      ! AND THEN INTERPOLATED AT THOSE POINTS (POSSIBLY IRREGULARLY DISTRIBUTED) OWNED BY EACH SINGLE PROCESS

      ! "GS"/"GE" STORE FIRST/LAST GLOBAL INDICES ALONG EACH DIRECTION FOR ALL PROCESS
      ALLOCATE(GS(3, 0:WORLD_SIZE-1), GE(3, 0:WORLD_SIZE-1))

      ! RETURN PROCESSORS GRID RESULTING MAXIMIZING THE NUMBER OF OVERLAPPING CALLS FOR INTERPOLATION (SEE BELOW)
      CALL BEST_CONFIG(DIMS)

      ! TEMPORARY MESSAGE
      IF (WORLD_RANK == 0) PRINT*, 'BEST PROCESSORS GRID: ', DIMS, ' -- ', NPTS

      ! CREATE TOPOLOGY
      CALL MPI_CART_CREATE(MPI_COMM_WORLD, NDIMS, DIMS, ISPERIODIC, REORDER, CARTOPO, IERR)

      ! RETURN PROCESS COORDINATES IN CURRENT TOPOLOGY
      CALL MPI_CART_COORDS(CARTOPO, WORLD_RANK, NDIMS, COORDS, IERR)

      ! RETURN FIRST/LAST-INDEX ("LS"/"LE") ALONG EACH DIRECTION FOR CALLING PROCESS. NOTE: FIRST POINT HAS ALWAYS INDEX EQUAL TO 1.
      CALL COORDS2INDEX(NPTS, DIMS, COORDS, LS, LE)

      GS(:, WORLD_RANK) = LS
      GE(:, WORLD_RANK) = LE

      ! MAKE ALL PROCESSES AWARE OF GLOBAL INDICES ALONG EACH AXIS
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, GS, 3, MPI_INTEGER, MPI_COMM_WORLD, IERR)
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, GE, 3, MPI_INTEGER, MPI_COMM_WORLD, IERR)

      ALLOCATE(POINTS_WORLD2RANK(0:WORLD_SIZE - 1), POINTS_RANK2WORLD(0:WORLD_SIZE - 1))

      ! MIN/MAX COMPUTATIONAL GRID INDICES CONTAINING INPUT POINTS
      PS = FLOOR(((FS - 1) * DS - OFF_AXIS) / DH) + 1
      PE = FLOOR(((FE - 1) * DS - OFF_AXIS) / DH) + 1

      POINTS_RANK2WORLD = 0

      ! DETERMINE IF CALLING PROCESS HAS AT LEAST ONE POINT FALLING INSIDE DOMAIN OF "I-TH" PROCESS
      DO I = 0, WORLD_SIZE - 1
        DO J = 1, 3
          BOOL(J) = ( (PS(J) .LT. GE(J, I)) .AND. (PS(J) .GE. (GS(J, I) - 1)) ) .OR.     &
                    ( (PE(J) .LT. GE(J, I)) .AND. (PE(J) .GE. (GS(J, I) - 1)) ) .OR.     &
                    ( (PE(J) .GE. GE(J, I)) .AND. (PS(J) .LT. (GS(J, I) - 1)) )
        ENDDO
        IF (ALL(BOOL .EQV. .TRUE.)) POINTS_RANK2WORLD(I) = 1
      ENDDO

      ! DETERMINE IF PROCESS "I" HAS AT LEAST ONE POINT FALLING INSIDE DOMAIN OF CALLING PROCESS
      CALL MPI_ALLTOALL(POINTS_RANK2WORLD, 1, MPI_INTEGER, POINTS_WORLD2RANK, 1, MPI_INTEGER, MPI_COMM_WORLD, IERR)

      ! NUMBER OF POINTS FOR CALLING PROCESS
      DO I = 1, 3
        M(I) = LE(I) - LS(I) + 1
      ENDDO

      ! ALLOCATE MEMORY FOR SPECTRUM
      ALLOCATE(SPEC(M(1), M(2), M(3)))

      ! COMPUTE SPECTRUM AND APPLY HERMITIAN SYMMETRY
      CALL COMPUTE_SPECTRUM(DS, LS, LE, DH, ACF, CL, SIGMA, HURST, SEED, SPEC, ET)

#ifdef TIMING
      INFO(5:6) = ET

      ! START TIMER
      CALL WATCH_START(TICTOC)
#endif

      ! TRANSFORM ALONG EACH DIRECTION
      CALL TRANSFORM_ALONG_Z(SPEC)
      CALL TRANSFORM_ALONG_Y(SPEC)
      CALL TRANSFORM_ALONG_X(SPEC)

#ifdef TIMING
      CALL WATCH_STOP(TICTOC)

      INFO(7) = TICTOC
#endif

      ! SCALING PARAMETER
      SCALING = 1._FPP / SQRT(REAL(NPTS(1), FPP) * REAL(NPTS(2), FPP) * REAL(NPTS(3), FPP) * DH**3)

      ALLOCATE(BUFFER(M(1), M(2), M(3)))

      ! SCALE IFFT
      DO K = 1, M(3)
        DO J = 1, M(2)
          DO I = 1, M(1)
            BUFFER(I, J, K) = REAL(SPEC(I, J, K), FPP) * SCALING
          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE(SPEC)

      ALLOCATE(VAR(0:WORLD_SIZE - 1), MU(0:WORLD_SIZE - 1), NPOINTS(0:WORLD_SIZE - 1))

      ! COMPUTE VARIANCE AND MEAN OF RANDOM FIELD FOR EACH SINGLE PROCESS
      VAR(WORLD_RANK) = VARIANCE(BUFFER)
      MU(WORLD_RANK)  = MEAN(BUFFER)

      ! NUMBER OF GRID POINTS PER PROCESS
      NPOINTS(WORLD_RANK) = PRODUCT(M)

      ! SHARE RESULTS
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, VAR, 1, REAL_TYPE, MPI_COMM_WORLD, IERR)
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, MU, 1, REAL_TYPE, MPI_COMM_WORLD, IERR)
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, NPOINTS, 1, MPI_INTEGER, MPI_COMM_WORLD, IERR)

      ! COMPUTE TOTAL VARIANCE ("INFO(3)") AND MEAN ("INFO(4)")
      CALL PARALLEL_VARIANCE(VAR, MU, NPOINTS, INFO(3), INFO(4))

      ! RETURN STANDARD DEVIATION
      INFO(3) = SQRT(INFO(3))

      ! DO WE NEED TO RESCALE THE RANDOM FIELD TO DESIRED (I.E. CONTINUOUS) STANDARD DEVIATION?
      IF (RESCALE .EQ. 1) THEN

        SCALING = SIGMA / INFO(3)

        DO K = 1, M(3)
          DO J = 1, M(2)
            DO I = 1, M(1)
              BUFFER(I, J, K) = BUFFER(I, J, K) * SCALING - INFO(4)
            ENDDO
          ENDDO
        ENDDO

      ENDIF

      DEALLOCATE(VAR, MU, NPOINTS)

      ! APPLY TAPER/MUTE
      DO I = 1, SIZE(POI, 2)
        PT(:) = POI(:, I) - OFF_AXIS(:)
        CALL TAPERING(DH, LS, LE, BUFFER, PT, MUTE, TAPER)
      ENDDO

      ! COPY RANDOM FIELD TO ARRAY WITH HALO
      ALLOCATE(DELTA(0:M(1), 0:M(2), 0:M(3)))

      DELTA = 0._FPP

      ! SCALE IFFT
      DO K = 1, M(3)
        DO J = 1, M(2)
          DO I = 1, M(1)
            DELTA(I, J, K) = BUFFER(I, J, K)
          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE(BUFFER)

#ifdef TIMING
      CALL WATCH_START(TICTOC)
#endif

      ! EXCHANGE HALO
      CALL EXCHANGE_HALO(CARTOPO, DELTA)

      ! INTERPOLATE RANDOM FIELD VALUES AT DESIRED OUTPUT LOCATIONS
      CALL GRID_MAPPING_STRUCT(DH, DELTA, DS, FS, FE, FIELD)

#ifdef TIMING
      CALL WATCH_STOP(TICTOC)

      INFO(8) = TICTOC
#endif

      ! RELEASE MEMORY
      DEALLOCATE(DELTA)
      DEALLOCATE(GS, GE)
      DEALLOCATE(POINTS_RANK2WORLD, POINTS_WORLD2RANK)

      CALL MPI_COMM_FREE(CARTOPO, IERR)

      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

    END SUBROUTINE SCARF3D_STRUCTURED_FFT

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE GRID_MAPPING(DH, DELTA, X, Y, Z, FIELD)

      REAL(FPP),                                       INTENT(IN)    :: DH                             !< FFT GRID-STEP
      REAL(FPP),                 DIMENSION(:,:,:),     INTENT(IN)    :: DELTA                          !< RANDOM PERTURBATIONS ON FFT GRID
      REAL(FPP),                 DIMENSION(:),         INTENT(IN)    :: X, Y, Z                        !< POINTS LOCATION
      REAL(FPP),                 DIMENSION(:),         INTENT(OUT)   :: FIELD
      INTEGER(IPP)                                                   :: I, J
      INTEGER(IPP)                                                   :: M, N, NP, NB
      INTEGER(IPP)                                                   :: MAXWIDTH, BLOCKWIDTH, NBLOCKS
      INTEGER(IPP)                                                   :: IERR, TRIP, REQ
      INTEGER(IPP),              DIMENSION(WORLD_SIZE)               :: SENDCOUNTS, RECVCOUNTS, SDISPLS, RDISPLS
      INTEGER(IPP),              DIMENSION(WORLD_SIZE)               :: SENDCOUNTS_N, RECVCOUNTS_N
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                        :: MAP, MAP_N
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                        :: P0, P1
      REAL(FPP),    ALLOCATABLE, DIMENSION(:)                        :: RV, SV
      REAL(FPP),    ALLOCATABLE, DIMENSION(:,:)                      :: SXYZ, RXYZ, SXYZ_N

      !-----------------------------------------------------------------------------------------------------------------------------

      FIELD = 0._FPP

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! IN ORDER NOT TO EXCEED THE PREVIOUS MEMORY PEAK, WE COMPUTE THE SIZE OF THE MINIMUM FFT GRID BLOCK "M", THE MAXIMUM NUMBER OF
      ! POINTS TO BE INTERPOLATED "NP" AMONGST ALL PROCESSES AND THE MAXIMUM NUMBER OF PROCESSES "NB" THAT COULD SEND POINTS TO ANOTHER
      ! ONE. WE USE THESE QUANTITIES TO PROVIDE AN UPPER BOUND "BLOCKWIDTH" TO THE NUMBER OF POINTS THAT CAN BE HANDLED BY ANY PROCESS
      ! AT ONCE. THIS WAY WE CAN ALLOCATE ARRAYS ONLY ONCE.
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      M = HUGE(1)

      ! MINIMUM COMPUTATIONAL DOMAIN SIZE
      DO I = 0, WORLD_SIZE - 1
        M = MIN(M, (GE(1, I) - GS(1, I) + 1) * (GE(2, I) - GS(2, I) + 1) * (GE(3, I) - GS(3, I) + 1))
      ENDDO

      ! NUMBER OF PROCESSES WITH AT LEAST ONE POINT FALLING INSIDE COMPUTATIONAL DOMAIN OF CALLING PROCESS
      NB = COUNT(POINTS_WORLD2RANK .EQ. 1)

      ! MAXIMUM NUMBER OF PROCESSES WITH AT LEAST ONE POINT FALLING INSIDE COMPUTATIONAL DOMAIN OF CALLING PROCESS (CONSERVATIVE ESTIMATE)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, NB, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, IERR)

      NP = SIZE(X)

      ! MAXIMUM NUMBER OF POINTS ANY PROCESS COULD SEND OR RECEIVE (CONSERVATIVE ESTIMATE)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, NP, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, IERR)

      ! MAX NUMBER OF POINTS THAT CAN BE HANDLED AT ONCE (*_N ARRAYS FOR NON-BLOCKING "ALLTOALLV" ARE INCLUDED)
      MAXWIDTH = (M * 2) / (9 + 4 * NB)

      BLOCKWIDTH = MAXWIDTH

      ! NUMBER OF BLOCKS REQUIRED TO COVER ALL POINTS
      NBLOCKS = NP / BLOCKWIDTH

      ! ADD ONE BLOCK MORE IF "NP" IS NOT MULTIPLE OF "BLOCKWIDTH"
      IF (MOD(NP, BLOCKWIDTH) .NE. 0) NBLOCKS = NBLOCKS + 1

      ALLOCATE(P0(NBLOCKS), P1(NBLOCKS))

      ! DEFINE FIRST/LAST POINT INDEX FOR EACH BLOCK
      DO I = 1, NBLOCKS
        P0(I) = (I - 1) * BLOCKWIDTH + 1
        P1(I) = I * BLOCKWIDTH
      ENDDO

      ! MAKE SURE THAT INDICES DO NOT EXCEED ACTUAL NUMBER OF POINTS IN CALLING PROCESS
      N = SIZE(X)

      ! LIMIT POINTS INDICES TO MAXIMUM NUMBER OF POINTS AVAILABLE TO CALLING PROCESS
      DO I = 1, NBLOCKS
        IF ( (P0(I) .LE. N) .AND. (P1(I) .GT. N) ) P1(I) = N
        IF ( (P0(I) .GT. N) .AND. (P1(I) .GT. N) ) P1(I) = P0(I) - 1
      ENDDO

      ! THE STRATEGY IS TO ALLOCATE ONLY ONCE THESE ARRAY, EVEN IF THEY MAY RESULT LARGER THAN NECESSARY
      ALLOCATE(SXYZ(3, BLOCKWIDTH), MAP(BLOCKWIDTH), RV(BLOCKWIDTH))
      ALLOCATE(RXYZ(3, BLOCKWIDTH * NB), SV(BLOCKWIDTH * NB))
      ALLOCATE(SXYZ_N(3, BLOCKWIDTH), MAP_N(BLOCKWIDTH))

      CALL MPI_TYPE_CONTIGUOUS(3, REAL_TYPE, TRIP, IERR)
      CALL MPI_TYPE_COMMIT(TRIP, IERR)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! FOR EACH BLOCK, WE BUILD AN ORDERED LIST OF POINTS TO BE SENT TO ALL OTHER PROCESSES AND THEIR NUMBER. THEN WE COMPUTE THE
      ! NUMBER OF POINTS TO BE RECEIVED. THESE POINTS ARE INTERPOLATED AND SENT BACK TO THE RIGHT OWNER. SIMILARLY, POINTS INTERPOLATED
      ! BY OTHER PROCESSES ARE COLLECTED AND STORED AT THE RIGHT LOCATION.
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! COLLECT DATA
      CALL ORDER_POINTS(P0(1), P1(1), DH, X, Y, Z, SXYZ, SENDCOUNTS, MAP)

      ! DETERMINE NUMBER OF POINTS TO BE RECEIVED FROM EACH PROCESS
      CALL MPI_ALLTOALL(SENDCOUNTS, 1, MPI_INTEGER, RECVCOUNTS, 1, MPI_INTEGER, MPI_COMM_WORLD, IERR)

      ! SET DISPLACEMENTS
      RDISPLS(1) = 0
      SDISPLS(1) = 0

      DO I = 2, WORLD_SIZE
        RDISPLS(I) = RDISPLS(I - 1) + RECVCOUNTS(I - 1)
        SDISPLS(I) = SDISPLS(I - 1) + SENDCOUNTS(I - 1)
      ENDDO

      DO J = 2, NBLOCKS

        ! FORWARD POINTS REFERRED TO PREVIOUS (J - 1) ITERATION
        CALL MPI_IALLTOALLV(SXYZ, SENDCOUNTS, SDISPLS, TRIP, RXYZ, RECVCOUNTS, RDISPLS, TRIP, MPI_COMM_WORLD, REQ, IERR)

        ! COLLECT DATA TO BE USED IN THE NEXT (J) ITERATION
        CALL ORDER_POINTS(P0(J), P1(J), DH, X, Y, Z, SXYZ_N, SENDCOUNTS_N, MAP_N)

        CALL MPI_ALLTOALL(SENDCOUNTS_N, 1, MPI_INTEGER, RECVCOUNTS_N, 1, MPI_INTEGER, MPI_COMM_WORLD, IERR)

        ! WAIT FOR DATA FROM PREVIOUS (J - 1) ITERATION
        CALL MPI_WAIT(REQ, MPI_STATUS_IGNORE, IERR)

        ! NUMBER OF POINTS RECEIVED: NEED TO INTEPOLATE ONLY THESE, NOT WHOLE "RXYZ" VECTOR
        N = SUM(RECVCOUNTS)

        ! INTERPOLATION
        CALL INTERPOLATE(N, RXYZ, DELTA, SV)

        ! COLLECT POINTS THAT HAVE BEEN INTERPOLATED BY OTHER PROCESSES
        CALL MPI_ALLTOALLV(SV, RECVCOUNTS, RDISPLS, REAL_TYPE, RV, SENDCOUNTS, SDISPLS, REAL_TYPE, MPI_COMM_WORLD, IERR)

        ! COPY INTERPOLATED DATA (IF ANY) TO RIGHT LOCATION
        DO I = 1, P1(J - 1) - P0(J - 1) + 1
          FIELD(MAP(I)) = RV(I)
        ENDDO

        ! ASSIGN DATA TO BE SENT AT NEXT ITERATION
        DO I = 1, BLOCKWIDTH
          SXYZ(1, I) = SXYZ_N(1, I)
          SXYZ(2, I) = SXYZ_N(2, I)
          SXYZ(3, I) = SXYZ_N(3, I)
          MAP(I)     = MAP_N(I)
        ENDDO

        RDISPLS(1) = 0
        SDISPLS(1) = 0

        SENDCOUNTS(1) = SENDCOUNTS_N(1)
        RECVCOUNTS(1) = RECVCOUNTS_N(1)

        DO I = 2, WORLD_SIZE
          SDISPLS(I)    = SDISPLS(I - 1) + SENDCOUNTS_N(I - 1)
          RDISPLS(I)    = RDISPLS(I - 1) + RECVCOUNTS_N(I - 1)
          SENDCOUNTS(I) = SENDCOUNTS_N(I)
          RECVCOUNTS(I) = RECVCOUNTS_N(I)
        ENDDO

      ENDDO

      ! FORWARD DATA
      CALL MPI_ALLTOALLV(SXYZ, SENDCOUNTS, SDISPLS, TRIP, RXYZ, RECVCOUNTS, RDISPLS, TRIP, MPI_COMM_WORLD, IERR)

      ! NUMBER OF POINTS RECEIVED: NEED TO INTEPOLATE ONLY THESE, NOT WHOLE "RXYZ" VECTOR
      N = SUM(RECVCOUNTS)

      ! INTERPOLATION
      CALL INTERPOLATE(N, RXYZ, DELTA, SV)

      ! BACKWARD DATA
      CALL MPI_ALLTOALLV(SV, RECVCOUNTS, RDISPLS, REAL_TYPE, RV, SENDCOUNTS, SDISPLS, REAL_TYPE, MPI_COMM_WORLD, IERR)

      ! COPY INTERPOLATED DATA TO RIGHT LOCATION. LIMIT THE MAXIMUM ITERATION INDEX IN CASE "NP" NOT MULTIPLE OF "BLOCKWIDTH"
      DO I = 1, P1(NBLOCKS) - P0(NBLOCKS) + 1
        FIELD(MAP(I)) = RV(I)
      ENDDO

      ! FREE RESOURCES
      CALL MPI_TYPE_FREE(TRIP, IERR)

      DEALLOCATE(SXYZ, RXYZ, MAP, RV, SV)
      DEALLOCATE(SXYZ_N, MAP_N)
      DEALLOCATE(P0, P1)

    END SUBROUTINE GRID_MAPPING

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE GRID_MAPPING_STRUCT(DH, DELTA, DS, FS, FE, FIELD)

      REAL(FPP),                                       INTENT(IN)    :: DH                             !< FFT GRID-STEP
      REAL(FPP),                 DIMENSION(:,:,:),     INTENT(IN)    :: DELTA                          !< RANDOM PERTURBATIONS ON FFT GRID
      REAL(FPP),                                       INTENT(IN)    :: DS
      INTEGER(IPP),              DIMENSION(3),         INTENT(IN)    :: FS, FE
      REAL(FPP),                 DIMENSION(:,:,:),     INTENT(OUT)   :: FIELD
      INTEGER(IPP)                                                   :: I, J, C
      INTEGER(IPP)                                                   :: M, N, NP, NB, NX, NY, NZ
      INTEGER(IPP)                                                   :: MAXWIDTH, BLOCKWIDTH, NBLOCKS
      INTEGER(IPP)                                                   :: IERR, TRIP, REQ
      INTEGER(IPP),              DIMENSION(3)                        :: DUM
      INTEGER(IPP),              DIMENSION(WORLD_SIZE)               :: SENDCOUNTS, RECVCOUNTS, SDISPLS, RDISPLS
      INTEGER(IPP),              DIMENSION(WORLD_SIZE)               :: SENDCOUNTS_N, RECVCOUNTS_N
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:,:)                      :: MAP, MAP_N
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:,:)                      :: P0, P1
      REAL(FPP),    ALLOCATABLE, DIMENSION(:)                        :: RV, SV
      REAL(FPP),    ALLOCATABLE, DIMENSION(:,:)                      :: SXYZ, RXYZ, SXYZ_N

      !-----------------------------------------------------------------------------------------------------------------------------

      FIELD = 0._FPP

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! IN ORDER NOT TO EXCEED THE PREVIOUS MEMORY PEAK, WE COMPUTE THE SIZE OF THE MINIMUM FFT GRID BLOCK "M", THE MAXIMUM NUMBER OF
      ! POINTS TO BE INTERPOLATED "NP" AMONGST ALL PROCESSES AND THE MAXIMUM NUMBER OF PROCESSES "NB" THAT COULD SEND POINTS TO ANOTHER
      ! ONE. WE USE THESE QUANTITIES TO PROVIDE AN UPPER BOUND "BLOCKWIDTH" TO THE NUMBER OF POINTS THAT CAN BE HANDLED BY ANY PROCESS
      ! AT ONCE. THIS WAY WE CAN ALLOCATE ARRAYS ONLY ONCE.
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      M = HUGE(1)

      ! MINIMUM COMPUTATIONAL DOMAIN SIZE
      DO I = 0, WORLD_SIZE - 1
        M = MIN(M, (GE(1, I) - GS(1, I) + 1) * (GE(2, I) - GS(2, I) + 1) * (GE(3, I) - GS(3, I) + 1))
      ENDDO

      ! NUMBER OF PROCESSES WITH AT LEAST ONE POINT FALLING INSIDE COMPUTATIONAL DOMAIN OF CALLING PROCESS
      NB = COUNT(POINTS_WORLD2RANK .EQ. 1)

      ! MAXIMUM NUMBER OF PROCESSES WITH AT LEAST ONE POINT FALLING INSIDE COMPUTATIONAL DOMAIN OF CALLING PROCESS (CONSERVATIVE ESTIMATE)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, NB, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, IERR)

      NX = FE(1) - FS(1) + 1
      NY = FE(2) - FS(2) + 1
      NZ = FE(3) - FS(3) + 1

      DUM = [NX, NY, NZ]

      ! MAXIMUM NUMBER OF POINTS ANY PROCESS COULD SEND OR RECEIVE (CONSERVATIVE ESTIMATE)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, DUM, 3, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, IERR)

      NX = DUM(1); NY = DUM(2); NZ = DUM(3)

      !NP = NX * NY * NZ

      ! MAX NUMBER OF POINTS THAT CAN BE HANDLED AT ONCE (*_N ARRAYS FOR NON-BLOCKING "ALLTOALLV" ARE INCLUDED)
      MAXWIDTH = (M * 2) / (9 + 4 * NB)

      ! NUMBER OF Z-PLANES THAT CAN BE PROCESSED AT ONCE
      BLOCKWIDTH = MAXWIDTH / (NX * NY)

      ! TEMPORARY WORKAROUND IF NO Z-PLANES CAN BE PREOCESSED: THIS MAY OCCURR IF INPUT MODEL DIMENSIONS ARE VERY SMALL!
      IF (BLOCKWIDTH .EQ. 0) THEN
        BLOCKWIDTH = 1
        MAXWIDTH   = NX * NY
      ENDIF

      ! NUMBER OF BLOCKS REQUIRED TO COVER HIGHEST NUMBER OF Z-PLANES
      NBLOCKS = NZ / BLOCKWIDTH

      IF (MOD(NZ, BLOCKWIDTH) .NE. 0) NBLOCKS = NBLOCKS + 1

      ALLOCATE(P0(3, NBLOCKS), P1(3, NBLOCKS))

      IF (BLOCKWIDTH .GE. 1) THEN

        DO I = 1, NBLOCKS
          P0(1, I) = 1
          P1(1, I) = FE(1) - FS(1) + 1
          P0(2, I) = 1
          P1(2, I) = FE(2) - FS(2) + 1
          P0(3, I) = (I - 1) * BLOCKWIDTH + 1
          P1(3, I) = I * BLOCKWIDTH
        ENDDO

        N = FE(3) - FS(3) + 1

        DO I = 1, NBLOCKS
          IF ( (P0(3, I) .LE. N) .AND. (P1(3, I) .GT. N) ) P1(3, I) = N
          IF ( (P0(3, I) .GT. N) .AND. (P1(3, I) .GT. N) ) P1(3, I) = P0(3, I) - 1
        ENDDO

      ELSE

        ! THIS PART MUST BE COMPLETED

        ! NUMBER OF Y-ROWS
        ! BLOCKWIDTH = MAXWIDTH / NX
        !
        ! IF (WORLD_RANK == 2) PRINT*, 'NUMBER Y ROWS ', BLOCKWIDTH
        !
        ! ! NUMBER OF BLOCKS OF Y-ROWS
        ! NBLOCKS = NY / BLOCKWIDTH
        !
        ! IF (MOD(NY, BLOCKWIDTH) .NE. 0) NBLOCKS = NBLOCKS + 1
        !
        ! ALLOCATE(P0(3, NBLOCKS*NZ), P1(3, NBLOCKS*NZ))
        !
        ! DO J = 1, NZ
        !   DO I = 1, NBLOCKS
        !     C        = (J - 1) * NBLOCKS + I
        !     P0(1, C) = 1
        !     P1(1, C) = FE(1) - FS(1) + 1
        !     P0(2, C) = (I - 1) * BLOCKWIDTH + 1
        !     P1(2, C) = I * BLOCKWIDTH
        !     P0(3, C) = J
        !     P1(3, C) = J
        !   ENDDO
        ! ENDDO
        !
        ! ! MAKE SURE WE DO NOT EXCEED INDICES
        !
        ! N = FE(2) - FS(2) + 1
        !
        ! DO J = 1, NZ
        !   DO I = 1, NBLOCKS
        !     C = (J - 1) * NBLOCKS + I
        !     IF ( (P0(2, C) .LE. N) .AND. (P1(2, C) .GT. N) ) P1(2, C) = N
        !     IF ( (P0(2, C) .GT. N) .AND. (P1(2, C) .GT. N) ) P1(2, C) = P0(2, C) - 1
        !   ENDDO
        ! ENDDO
        !
        ! N = FE(3) - FS(3) + 1
        !
        ! DO J = 1, NZ
        !   DO I = 1, NBLOCKS
        !     C = (J - 1) * NBLOCKS + I
        !     IF ( (P0(3, C) .GT. N) .AND. (P1(3, C) .GT. N) ) P1(3, C) = P0(3, C) - 1
        !   ENDDO
        ! ENDDO
        !
        ! BLOCKWIDTH = NX * BLOCKWIDTH
        !
        ! NBLOCKS = NBLOCKS * NZ
        !
        ! IF (WORLD_RANK == 2) PRINT*, 'NUMBER OF BLOCKWIDTH (POINTS) ', BLOCKWIDTH, 'BLOCKS ', NBLOCKS, NBLOCKS/NZ

      ENDIF

      ! THE STRATEGY IS TO ALLOCATE ONLY ONCE THESE ARRAY, EVEN IF THEY MAY RESULT LARGER THAN NECESSARY
      ALLOCATE(SXYZ(3, MAXWIDTH), MAP(3, MAXWIDTH), RV(MAXWIDTH))
      ALLOCATE(RXYZ(3, MAXWIDTH * NB), SV(MAXWIDTH * NB))
      ALLOCATE(SXYZ_N(3, MAXWIDTH), MAP_N(3, MAXWIDTH))

      CALL MPI_TYPE_CONTIGUOUS(3, REAL_TYPE, TRIP, IERR)
      CALL MPI_TYPE_COMMIT(TRIP, IERR)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! FOR EACH BLOCK, WE BUILD AN ORDERED LIST OF POINTS TO BE SENT TO ALL OTHER PROCESSES AND THEIR NUMBER. THEN WE COMPUTE THE
      ! NUMBER OF POINTS TO BE RECEIVED. THESE POINTS ARE INTERPOLATED AND SENT BACK TO THE RIGHT OWNER. SIMILARLY, POINTS INTERPOLATED
      ! BY OTHER PROCESSES ARE COLLECTED AND STORED AT THE RIGHT LOCATION.
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! COLLECT DATA
      CALL ORDER_POINTS_STRUCT(P0(:, 1), P1(:, 1), DH, FS, FE, DS, SXYZ, SENDCOUNTS, MAP)
      CALL MPI_ALLTOALL(SENDCOUNTS, 1, MPI_INTEGER, RECVCOUNTS, 1, MPI_INTEGER, MPI_COMM_WORLD, IERR)

      ! SET DISPLACEMENTS
      RDISPLS(1) = 0
      SDISPLS(1) = 0

      DO I = 2, WORLD_SIZE
        RDISPLS(I) = RDISPLS(I - 1) + RECVCOUNTS(I - 1)
        SDISPLS(I) = SDISPLS(I - 1) + SENDCOUNTS(I - 1)
      ENDDO

      DO J = 2, NBLOCKS

        ! FORWARD OLD DATA
        CALL MPI_IALLTOALLV(SXYZ, SENDCOUNTS, SDISPLS, TRIP, RXYZ, RECVCOUNTS, RDISPLS, TRIP, MPI_COMM_WORLD, REQ, IERR)

        ! COLLECT NEW DATA
        CALL ORDER_POINTS_STRUCT(P0(:, J), P1(:, J), DH, FS, FE, DS, SXYZ_N, SENDCOUNTS_N, MAP_N)
        CALL MPI_ALLTOALL(SENDCOUNTS_N, 1, MPI_INTEGER, RECVCOUNTS_N, 1, MPI_INTEGER, MPI_COMM_WORLD, IERR)

        ! WAIT FOR OLD DATA
        CALL MPI_WAIT(REQ, MPI_STATUS_IGNORE, IERR)

        ! NUMBER OF POINTS RECEIVED: NEED TO INTEPOLATE ONLY THESE, NOT WHOLE "RXYZ" VECTOR
        N = SUM(RECVCOUNTS)

        ! INTERPOLATION
        CALL INTERPOLATE(N, RXYZ, DELTA, SV)

        ! BACKWARD DATA
        CALL MPI_ALLTOALLV(SV, RECVCOUNTS, RDISPLS, REAL_TYPE, RV, SENDCOUNTS, SDISPLS, REAL_TYPE, MPI_COMM_WORLD, IERR)

        ! COPY INTERPOLATED DATA TO RIGHT LOCATION
        DO I = 1, PRODUCT(P1(:, J - 1) - P0(:, J - 1) + 1)
          FIELD(MAP(1, I), MAP(2, I), MAP(3, I)) = RV(I)
        ENDDO

        ! ASSIGN "NEW DATA"
        DO I = 1, PRODUCT(P1(:, J) - P0(:, J) + 1)
          SXYZ(1, I) = SXYZ_N(1, I)
          SXYZ(2, I) = SXYZ_N(2, I)
          SXYZ(3, I) = SXYZ_N(3, I)
          MAP(1, I)  = MAP_N(1, I)
          MAP(2, I)  = MAP_N(2, I)
          MAP(3, I)  = MAP_N(3, I)
        ENDDO

        RDISPLS(1) = 0
        SDISPLS(1) = 0

        ! ASSIGN "NEW DATA"
        SENDCOUNTS(1) = SENDCOUNTS_N(1)
        RECVCOUNTS(1) = RECVCOUNTS_N(1)

        DO I = 2, WORLD_SIZE
          SDISPLS(I)    = SDISPLS(I - 1) + SENDCOUNTS_N(I - 1)
          RDISPLS(I)    = RDISPLS(I - 1) + RECVCOUNTS_N(I - 1)
          SENDCOUNTS(I) = SENDCOUNTS_N(I)
          RECVCOUNTS(I) = RECVCOUNTS_N(I)
        ENDDO

      ENDDO

      ! FORWARD DATA
      CALL MPI_ALLTOALLV(SXYZ, SENDCOUNTS, SDISPLS, TRIP, RXYZ, RECVCOUNTS, RDISPLS, TRIP, MPI_COMM_WORLD, IERR)

      ! NUMBER OF POINTS RECEIVED: NEED TO INTEPOLATE ONLY THESE, NOT WHOLE "RXYZ" VECTOR
      N = SUM(RECVCOUNTS)

      ! INTERPOLATION
      CALL INTERPOLATE(N, RXYZ, DELTA, SV)

      ! BACKWARD DATA
      CALL MPI_ALLTOALLV(SV, RECVCOUNTS, RDISPLS, REAL_TYPE, RV, SENDCOUNTS, SDISPLS, REAL_TYPE, MPI_COMM_WORLD, IERR)

      ! COPY INTERPOLATED DATA TO RIGHT LOCATION. LIMIT THE MAXIMUM ITERATION INDEX IN CASE "NP" NOT MULTIPLE OF "BLOCKWIDTH"
      DO I = 1, PRODUCT(P1(:, NBLOCKS) - P0(:, NBLOCKS) + 1)
        FIELD(MAP(1, I), MAP(2, I), MAP(3, I)) = RV(I)
      ENDDO

      ! FREE RESOURCES
      CALL MPI_TYPE_FREE(TRIP, IERR)

      DEALLOCATE(SXYZ, RXYZ, MAP, RV, SV)
      DEALLOCATE(SXYZ_N, MAP_N)
      DEALLOCATE(P0, P1)

    END SUBROUTINE GRID_MAPPING_STRUCT

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE ORDER_POINTS(P0, P1, DH, X, Y, Z, XYZ, SENDCOUNTS, MAP)

      ! RE-ORDER POINTS TO BE SENT TO EACH PROCESS IN A LIST "XYZ" CONTAINING GRID INDICES FOR THAT PROCESS. THE NUMBER OF POINTS FOR
      ! EACH PROCESS IS STORED IN "SENDCOUNTS", WHILE THE ORIGINAL INDEX OF A POINT IS IN "MAP".

      INTEGER(IPP),                 INTENT(IN)  :: P0, P1
      REAL(FPP),                    INTENT(IN)  :: DH
      REAL(FPP),    DIMENSION(:),   INTENT(IN)  :: X, Y, Z
      REAL(FPP),    DIMENSION(:,:), INTENT(OUT) :: XYZ
      INTEGER(IPP), DIMENSION(:),   INTENT(OUT) :: SENDCOUNTS
      INTEGER(IPP), DIMENSION(:),   INTENT(OUT) :: MAP
      INTEGER(IPP)                              :: L, P, N, C
      INTEGER(IPP), DIMENSION(3)                :: CONST, M
      LOGICAL                                   :: BOOL
      REAL(FPP)                                 :: I, J, K

      !-----------------------------------------------------------------------------------------------------------------------------

      SENDCOUNTS = 0

      MAP = 0

      C = 0

      ! LOOP OVER PROCESSES
      DO L = 0, WORLD_SIZE - 1

        ! SKIP IF CALLING PROCESS HAS NO POINTS FALLING INSIDE CUBOID OF "L-TH" PROCESS
        IF (POINTS_RANK2WORLD(L) .EQ. 0) CYCLE

        DO P = 1, 3
          CONST(P) = GS(P, L)
          M(P)     = GE(P, L) - GS(P, L) + 2
        ENDDO

        N = C

        DO P = P0, P1

          ! TRANSFORM POINT COORDINATES INTO GRID INDICES FOR L-TH PROCESS
          I = (X(P) - OFF_AXIS(1)) / DH + 3._FPP - CONST(1)
          J = (Y(P) - OFF_AXIS(2)) / DH + 3._FPP - CONST(2)
          K = (Z(P) - OFF_AXIS(3)) / DH + 3._FPP - CONST(3)

          ! CHECK IF POINT IS WITHIN BLOCK OF L-TH PROCESS
          BOOL = (I .GE. 1) .AND. (I .LT. M(1)) .AND. (J .GE. 1) .AND. (J .LT. M(2)) .AND. (K .GE. 1) .AND. (K .LT. M(3))

          IF (BOOL .EQV. .TRUE.) THEN
            C         = C + 1
            XYZ(:, C) = [I, J, K]
            MAP(C)    = P
          ENDIF

        ENDDO

        ! NUMBER OF POINTS TO BE SENT TO L-TH PROCESS
        SENDCOUNTS(L + 1) = C - N

      ENDDO

    END SUBROUTINE ORDER_POINTS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE ORDER_POINTS_STRUCT(P0, P1, DH, FS, FE, DS, XYZ, SENDCOUNTS, MAP)

      INTEGER(IPP), DIMENSION(3),   INTENT(IN)  :: P0, P1
      REAL(FPP),                    INTENT(IN)  :: DH
      INTEGER(IPP), DIMENSION(3),   INTENT(IN)  :: FS, FE
      REAL(FPP),                    INTENT(IN)  :: DS
      REAL(FPP),    DIMENSION(:,:), INTENT(OUT) :: XYZ
      INTEGER(IPP), DIMENSION(:),   INTENT(OUT) :: SENDCOUNTS
      INTEGER(IPP), DIMENSION(:,:), INTENT(OUT) :: MAP
      INTEGER(IPP)                              :: L, P, N, C
      INTEGER(IPP)                              :: NX, NXY
      INTEGER(IPP)                              :: I0, J0, K0
      INTEGER(IPP), DIMENSION(3)                :: CONST, M
      LOGICAL                                   :: BOOL
      REAL(FPP)                                 :: I, J, K
      REAL(FPP),    DIMENSION(P0(3):P1(3))      :: Z
      REAL(FPP),    DIMENSION(P0(2):P1(2))      :: Y
      REAL(FPP),    DIMENSION(P0(1):P1(1))      :: X

      !-----------------------------------------------------------------------------------------------------------------------------

      DO K0 = P0(3), P1(3)
        Z(K0) = ((FS(3) - 1 + K0 - 1) * DS - OFF_AXIS(3)) / DH + 3._FPP
      ENDDO

      DO J0 = P0(2), P1(2)
        Y(J0) = ((FS(2) - 1 + J0 - 1) * DS - OFF_AXIS(2)) / DH + 3._FPP
      ENDDO

      DO I0 = P0(1), P1(1)
        X(I0) = ((FS(1) - 1 + I0 - 1) * DS - OFF_AXIS(1)) / DH + 3._FPP
      ENDDO

      SENDCOUNTS = 0

      MAP = 0

      C = 0

      ! LOOP OVER PROCESSES
      DO L = 0, WORLD_SIZE - 1

        ! SKIP IF CALLING PROCESS HAS NO POINTS FALLING INSIDE CUBOID OF "L-TH" PROCESS
        IF (POINTS_RANK2WORLD(L) .EQ. 0) CYCLE

        DO P = 1, 3
          CONST(P) = GS(P, L)
          M(P)     = GE(P, L) - GS(P, L) + 2
        ENDDO

        N = C

        DO K0 = P0(3), P1(3)

          K = Z(K0) - CONST(3)

          DO J0 = P0(2), P1(2)

            J = Y(J0) - CONST(2)

            DO I0 = P0(1), P1(1)

              I = X(I0) - CONST(1)

              BOOL = (I .GE. 1) .AND. (I .LT. M(1)) .AND. (J .GE. 1) .AND. (J .LT. M(2)) .AND. (K .GE. 1) .AND. (K .LT. M(3))

              IF (BOOL .EQV. .TRUE.) THEN
                C         = C + 1
                XYZ(:, C) = [I, J, K]
                MAP(:, C) = [I0, J0, K0]
              ENDIF

            ENDDO

          ENDDO

        ENDDO

        SENDCOUNTS(L + 1) = C - N

      ENDDO

    END SUBROUTINE ORDER_POINTS_STRUCT

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

    SUBROUTINE COMPUTE_SPECTRUM(DS, LS, LE, DH, ACF, CL, SIGMA, HURST, SEED, SPEC, TIME)

      ! COMPUTE THE SPECTRUM OF A RANDOM FIELD CHARACTERISED BY A SPECIFIC AUTOCORRELATION FUNCTION, CORRELATION LENGTH, STANDARD DEVIATION,
      ! AND HURST EXPONENT (NOT USED FOR GAUSSIAN FIELDS). BASED ON THE FOURIER INTEGRAL METHOD OF PARDO-IGUZQUIZA AND CHICA-OLMO.

      REAL(FPP),                                                        INTENT(IN)  :: DS
      INTEGER(IPP),     DIMENSION(3),                                   INTENT(IN)  :: LS, LE                !< GLOBAL INDICES
      REAL(FPP),                                                        INTENT(IN)  :: DH                    !< GRID-STEP
      INTEGER(IPP),                                                     INTENT(IN)  :: ACF                   !< AUTOCORRELATION FUNCTION
      REAL(FPP),        DIMENSION(3),                                   INTENT(IN)  :: CL                    !< CORRELATION LENGTH (CAN BE ANISOTROPIC)
      REAL(FPP),                                                        INTENT(IN)  :: SIGMA                 !< STANDARD DEVIATION (SIGMA%/100)
      REAL(FPP),                                                        INTENT(IN)  :: HURST                 !< HURST EXPONENT
      INTEGER(IPP),                                                     INTENT(IN)  :: SEED                  !< INITIAL SEED NUMBER
      COMPLEX(FPP),     DIMENSION(LS(1):LE(1),LS(2):LE(2),LS(3):LE(3)), INTENT(OUT) :: SPEC                  !< SPECTRUM
      REAL(FPP),        DIMENSION(2),                                   INTENT(OUT) :: TIME                  !< ELAPSED TIME
      INTEGER(IPP)                                                                  :: I, J, K               !< INDICES
      REAL(FPP)                                                                     :: BUTTER, NUM, AMP      !< USED TO COMPUTE SPECTRUM
      REAL(FPP)                                                                     :: KC, KR                !< USED TO COMPUTE SPECTRUM
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

      TIME(:) = 0._FPP

#ifdef TIMING
      ! START TIMER
      CALL WATCH_START(TICTOC)
#endif

      ! RESOLUTION IN WAVENUMBER DOMAIN ALONG EACH DIRECTION
      DO I = 1, 3
        DK(I) = 2._FPP * PI / (REAL(NPTS(I), FPP) * DH)
      ENDDO

      ! NYQUIST WAVENUMBER
      !KNYQ = PI / DH

      ! CORNER WAVENUMBER FOR FILTERING SPECTRUM
      KC = 2._FPP * PI / DS 

      ! VECTORS GO FROM 0 TO NYQUIST AND THEN BACK AGAIN UNTIL DK
      KX = [[(I * DK(1), I = 0, NPTS(1)/2)], [(I * DK(1), I = NPTS(1)/2-1, 1, -1)]]
      KY = [[(J * DK(2), J = 0, NPTS(2)/2)], [(J * DK(2), J = NPTS(2)/2-1, 1, -1)]]
      KZ = [[(K * DK(3), K = 0, NPTS(3)/2)], [(K * DK(3), K = NPTS(3)/2-1, 1, -1)]]

      ! COMPUTE PART OF POWER SPECTRAL DENSITY OUTSIDE LOOP
      IF (ACF .EQ. 0) THEN
        NUM = 8._FPP * SQRT(PI**3) * GAMMA(HURST + 1.5_FPP) * SIGMA**2 * PRODUCT(CL) / GAMMA(HURST)
        FUN => VK_PSDF
      ELSEIF (ACF .EQ. 1) THEN
        NUM = SIGMA**2 * PRODUCT(CL) * SQRT(PI**3)
        FUN => GS_PSDF
      ENDIF

      ! EACH PROCESS GENERATE ITS SET OF RANDOM NUMBERS
      CPTR = PRNG(SEED, LS, LE, NPTS)

      CALL C_F_POINTER(CPTR, R, [(LE(1) - LS(1) + 1) * (LE(2) - LS(2) + 1) * (LE(3) - LS(3) + 1)])

      HARVEST(LS(1):LE(1), LS(2):LE(2), LS(3):LE(3)) => R

      ! COMPUTE SPECTRUM
      DO K = LS(3), LE(3)

        DO J = LS(2), LE(2)

          DO I = LS(1), LE(1)

            ! RADIAL WAVENUMBER
            KR = SQRT(KX(I)**2 + KY(J)**2 + KZ(K)**2)

            ! FOURTH-ORDER LOW-PASS BUTTERWORTH FILTER
            BUTTER = 1._FPP / SQRT(1._FPP + (KR / KC)**(2 * 4))

            ! NOW "KR" IS THE PRODUCT "K * CL"
            KR = (KX(I) * CL(1))**2 + (KY(J) * CL(2))**2 + (KZ(K) * CL(3))**2

            ! COMPLETE POWER SPECTRAL DENSITY AND GO TO AMPLITUDE SPECTRUM
            AMP = SQRT(NUM / FUN(KR, HURST))

            ! APPLY FILTER
            AMP = AMP * BUTTER

            !R = HARVEST(I - LS(1) + 1, J - LS(2) + 1, K - LS(3) + 1) * 2._FPP * PI
            HARVEST(I, J, K) = HARVEST(I, J, K) * 2._FPP * PI

            ! COMBINE AMPLITUDE AND RANDOM PHASE
            SPEC(I, J, K) = CMPLX(COS(HARVEST(I, J, K)) * AMP, SIN(HARVEST(I, J, K)) * AMP, FPP)

          ENDDO
        ENDDO
      ENDDO

      NULLIFY(R, HARVEST)
      CALL FREE_MEM(CPTR)

#ifdef TIMING
      CALL WATCH_STOP(TICTOC)

      TIME(1) = TICTOC

      CALL WATCH_START(TICTOC)
#endif

      CALL ENFORCE_SYMMETRY(LS, LE, SPEC)

#ifdef TIMING
      CALL WATCH_STOP(TICTOC)

      TIME(2) = TICTOC
#endif

    END SUBROUTINE COMPUTE_SPECTRUM

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

      INTEGER(IPP),              DIMENSION(3),                                   INTENT(IN)    :: FS, FE
      COMPLEX(FPP),              DIMENSION(FS(1):FE(1),FS(2):FE(2),FS(3):FE(3)), INTENT(INOUT) :: SPEC
      INTEGER(IPP)                                                                             :: I
      INTEGER(IPP)                                                                             :: RANK, NTASKS
      INTEGER(IPP)                                                                             :: NEWCOMM, IERR, COLOR
      INTEGER(IPP),              DIMENSION(3)                                                  :: NYQUIST
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                                                  :: NEW2WORLD

      LOGICAL,                   DIMENSION(2)                                                  :: BOOL

      !-----------------------------------------------------------------------------------------------------------------------------

      ! INDEX OF NYQUIST WAVENUMBER
      DO I = 1, 3
        NYQUIST(I) = NPTS(I) / 2 + 1
      ENDDO

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! CREATE A NEW COMMUNICATOR CONTAINING ONLY THOSE PROCESSES CROSSING X=1 AND X=NPTS(1)/2 + 1
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      COLOR = 0

      BOOL(1) = FS(1) .EQ. 1
      BOOL(2) = (FS(1) .LE. NYQUIST(1)) .AND. (FE(1) .GE. NYQUIST(1))

      IF (BOOL(1)) COLOR = 1
      IF (BOOL(2)) COLOR = 2

      IF (ALL(BOOL)) COLOR = 3

      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, COLOR, WORLD_RANK, NEWCOMM, IERR)

      CALL MPI_COMM_RANK(NEWCOMM, RANK, IERR)
      CALL MPI_COMM_SIZE(NEWCOMM, NTASKS, IERR)

      ALLOCATE(NEW2WORLD(0:NTASKS - 1))

      ! I-TH RANK IN "NEWCOMM" HAS ITS GLOBAL RANK CONTAINED IN "NEW2WORLD"
      CALL MAP_RANKS(NEWCOMM, MPI_COMM_WORLD, NEW2WORLD)

      ! IF ONLY ONE OF CONDITIONS ABOVE IS TRUE, WE CAN WORK ON I=1 AND I=NYQUIST(1) AT THE SAME TIME...
      IF ( (COLOR .EQ. 1) .OR. (COLOR .EQ. 2) ) THEN

        IF (BOOL(1)) THEN
          I = 1
        ELSE
          I = NYQUIST(1)
        ENDIF

        CALL EXCHANGE_CONJG

      ! ... OTHERWISE WORK IN TWO STEPS
      ELSEIF (COLOR .EQ. 3) THEN

        I = 1

        CALL EXCHANGE_CONJG

        I = NYQUIST(1)

        CALL EXCHANGE_CONJG

      ENDIF

      DEALLOCATE(NEW2WORLD)

      CALL MPI_COMM_FREE(NEWCOMM, IERR)

      CONTAINS

      !-----------------------------------------------------------------------------------------------------------------------------

        SUBROUTINE EXCHANGE_CONJG

          COMPLEX(FPP), DIMENSION(FS(2):FE(2),FS(3):FE(3)) :: BUFFER, SENDBUF
          INTEGER(IPP)                                     :: J, K, P, L, C, CY, CZ
          INTEGER(IPP)                                     :: JS, JE, KS, KE
          INTEGER(IPP)                                     :: J0, J1, K0, K1
          INTEGER(IPP)                                     :: IERR
          INTEGER(IPP), DIMENSION(2)                       :: RSIZES, RSUBSIZES, RSTARTS
          INTEGER(IPP), DIMENSION(3)                       :: SSIZES, SSUBSIZES, SSTARTS
          INTEGER(IPP), DIMENSION(0:NTASKS-1)              :: SENDCOUNTS, RECVCOUNTS
          INTEGER(IPP), DIMENSION(0:NTASKS-1)              :: SDISPLS, RDISPLS
          INTEGER(IPP), DIMENSION(0:NTASKS-1)              :: SENDTYPES, RECVTYPES
          LOGICAL                                          :: BOOL

          !-------------------------------------------------------------------------------------------------------------------------

          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
          ! APPLY CONDITIONS OF FIGURE 4 FOR ALPHA, PSI, XSI AND ETA
          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

          BOOL = (I .EQ. 1) .AND. (FS(2) .EQ. 1) .AND. (FS(3) .EQ. 1)

          IF (BOOL) SPEC(I, 1, 1) = 0._FPP                                                         !< ALPHA MUST BE ZERO FOR ZERO-MEAN FIELD

          BOOL = (I .EQ. NYQUIST(1)) .AND. (FS(2) .EQ. 1) .AND. (FS(3) .EQ. 1)

          IF (BOOL) SPEC(I, 1, 1) = REAL(SPEC(I, 1, 1), FPP)                                       !< ALPHA MUST BE REAL

          BOOL = (FS(2) .LE. NYQUIST(2)) .AND. (FE(2) .GE. NYQUIST(2)) .AND. (FS(3) .LE. NYQUIST(3)) .AND. (FE(3) .GE. NYQUIST(3))

          IF (BOOL) SPEC(I, NYQUIST(2), NYQUIST(3)) = REAL(SPEC(I, NYQUIST(2), NYQUIST(3)), FPP)   !< PSI MUST BE REAL

          BOOL = (FS(2) .EQ. 1) .AND. (FS(3) .LE. NYQUIST(3)) .AND. (FE(3) .GE. NYQUIST(3))

          IF (BOOL) SPEC(I, 1, NYQUIST(3)) = REAL(SPEC(I, 1, NYQUIST(3)), FPP)                     !< XSI MUST BE REAL

          BOOL = (FS(2) .LE. NYQUIST(2)) .AND. (FE(2) .GE. NYQUIST(2)) .AND. (FS(3) .EQ. 1)

          IF (BOOL) SPEC(I, NYQUIST(2), 1) = REAL(SPEC(I, NYQUIST(2), 1), FPP)                     !< ETA MUST BE REAL

          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
          ! APPLY CONDITIONS OF FIGURE 4 AT REMAINING REGIONS. EXCHANGE DATA AMONGST PROCESSES. BETA- AND DELTA- REGIONS ARE TREATED
          ! SEPARATELY.
          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

          SENDCOUNTS = 0
          RECVCOUNTS = 0

          SDISPLS = 0
          RDISPLS = 0

          SSIZES = [FE(1) - FS(1) + 1, FE(2) - FS(2) + 1, FE(3) - FS(3) + 1]

          RSIZES = [FE(2) - FS(2) + 1, FE(3) - FS(3) + 1]

          ! COPY RELEVANT PART OF "SPEC". THIS IS NOT STRICTLY NECESSARY BUT MAY IMPROVE CACHE ACCESS.
          DO K = FS(3), FE(3)
            DO J = FS(2), FE(2)
              BUFFER(J, K) = SPEC(I, J, K)
            ENDDO
          ENDDO

          ! EXCHANGE DATA IN THE WHOLE REGION, EXCEPT BETA- AND DELTA-REGION
          DO P = 0, NTASKS - 1

            L = NEW2WORLD(P)                                                              !< GLOBAL INDEX OF P-TH PROCESS

            SENDTYPES(P) = COMPLEX_TYPE
            RECVTYPES(P) = COMPLEX_TYPE

            ! PROJECT POINTS OF PROCESS "P" INTO CONJUGATED FIELD.
            JS = NPTS(2) - GE(2, L) + 2
            JE = NPTS(2) - GS(2, L) + 2

            KS = NPTS(3) - GE(3, L) + 2
            KE = NPTS(3) - GS(3, L) + 2

            J0 = MAX(FS(2), JS)
            J1 = MIN(FE(2), JE)

            K0 = MAX(FS(3), KS)
            K1 = MIN(FE(3), KE)

            ! DETERMINE DATA TO BE SENT/RECEIVED TO/FROM PROCESS "P"
            IF ( (K0 .LE. K1) .AND. (J0 .LE. J1) ) THEN

              CY = J1 + J0
              CZ = K1 + K0

              DO K = K0, K1
                DO J = J0, J1
                  SENDBUF(CY - J, CZ - K) = BUFFER(J, K)
                ENDDO
              ENDDO

              RSUBSIZES = [J1 - J0 + 1, K1 - K0 + 1]
              RSTARTS   = [J0 - FS(2), K0 - FS(3)]

              CALL MPI_TYPE_CREATE_SUBARRAY(2, RSIZES, RSUBSIZES, RSTARTS, MPI_ORDER_FORTRAN, COMPLEX_TYPE, SENDTYPES(P), IERR)
              CALL MPI_TYPE_COMMIT(SENDTYPES(P), IERR)

              CALL MPI_TYPE_CREATE_SUBARRAY(2, RSIZES, RSUBSIZES, RSTARTS, MPI_ORDER_FORTRAN, COMPLEX_TYPE, RECVTYPES(P), IERR)
              CALL MPI_TYPE_COMMIT(RECVTYPES(P), IERR)

              SENDCOUNTS(P) = 1
              RECVCOUNTS(P) = 1

            ENDIF

          ENDDO

          CALL MPI_ALLTOALLW(SENDBUF, SENDCOUNTS, SDISPLS, SENDTYPES, BUFFER, RECVCOUNTS, RDISPLS, RECVTYPES, NEWCOMM, IERR)

          ! FREE RESOURCES
          DO P = 0, NTASKS - 1
            IF (SENDCOUNTS(P) .EQ. 1) THEN
              CALL MPI_TYPE_FREE(SENDTYPES(P), IERR)
              CALL MPI_TYPE_FREE(RECVTYPES(P), IERR)
            ENDIF
          ENDDO

          SENDCOUNTS = 0
          RECVCOUNTS = 0

          ! NOW EXCHANGE DATA IN THE BETA- AND DELTA-REGION
          DO P = 0, NTASKS - 1

            L = NEW2WORLD(P)                                                              !< GLOBAL INDEX OF P-TH PROCESS

            SENDTYPES(P) = COMPLEX_TYPE
            RECVTYPES(P) = COMPLEX_TYPE

            ! PROJECT POINTS OF PROCESS "P" INTO CONJUGATED FIELD.
            IF (GS(3, L) .EQ. 1) THEN

              JS = NPTS(2) - GE(2, L) + 2
              JE = NPTS(2) - GS(2, L) + 2

              KS = 1
              KE = 1

              J0 = MAX(FS(2), JS)
              J1 = MIN(FE(2), JE)

              K0 = MAX(FS(3), KS)
              K1 = MIN(FE(3), KE)

              ! DETERMINE DATA TO BE SENT/RECEIVED TO/FROM PROCESS "P"
              IF ( (K0 .LE. K1) .AND. (J0 .LE. J1) ) THEN

                RSUBSIZES = [J1 - J0 + 1, K1 - K0 + 1]
                RSTARTS   = [J0 - FS(2), K0 - FS(3)]

                CY = J1 + J0
                CZ = K1 + K0

                DO K = K0, K1
                  DO J = J0, J1
                    SENDBUF(CY - J, CZ - K) = BUFFER(J, K)
                  ENDDO
                ENDDO

                CALL MPI_TYPE_CREATE_SUBARRAY(2, RSIZES, RSUBSIZES, RSTARTS, MPI_ORDER_FORTRAN, COMPLEX_TYPE, SENDTYPES(P), IERR)
                CALL MPI_TYPE_COMMIT(SENDTYPES(P), IERR)

                CALL MPI_TYPE_CREATE_SUBARRAY(2, RSIZES, RSUBSIZES, RSTARTS, MPI_ORDER_FORTRAN, COMPLEX_TYPE, RECVTYPES(P), IERR)
                CALL MPI_TYPE_COMMIT(RECVTYPES(P), IERR)

                SENDCOUNTS(P) = 1
                RECVCOUNTS(P) = 1

              ENDIF

            ENDIF

            ! PROJECT POINTS OF PROCESS "P" INTO CONJUGATED FIELD.
            IF (GS(2, L) .EQ. 1) THEN

              JS = 1
              JE = 1

              KS = NPTS(3) - GE(3, L) + 2
              KE = NPTS(3) - GS(3, L) + 2

              J0 = MAX(FS(2), JS)
              J1 = MIN(FE(2), JE)

              K0 = MAX(FS(3), KS)
              K1 = MIN(FE(3), KE)

              ! DETERMINE DATA TO BE SENT/RECEIVED TO/FROM PROCESS "P"
              IF ( (K0 .LE. K1) .AND. (J0 .LE. J1) ) THEN

                RSUBSIZES = [J1 - J0 + 1, K1 - K0 + 1]
                RSTARTS   = [J0 - FS(2), K0 - FS(3)]

                CY = J1 + J0
                CZ = K1 + K0

                DO K = K0, K1
                  DO J = J0, J1
                    SENDBUF(CY - J, CZ - K) = BUFFER(J, K)
                  ENDDO
                ENDDO

                CALL MPI_TYPE_CREATE_SUBARRAY(2, RSIZES, RSUBSIZES, RSTARTS, MPI_ORDER_FORTRAN, COMPLEX_TYPE, SENDTYPES(P), IERR)
                CALL MPI_TYPE_COMMIT(SENDTYPES(P), IERR)

                CALL MPI_TYPE_CREATE_SUBARRAY(2, RSIZES, RSUBSIZES, RSTARTS, MPI_ORDER_FORTRAN, COMPLEX_TYPE, RECVTYPES(P), IERR)
                CALL MPI_TYPE_COMMIT(RECVTYPES(P), IERR)

                SENDCOUNTS(P) = 1
                RECVCOUNTS(P) = 1

              ENDIF

            ENDIF

          ENDDO

          CALL MPI_ALLTOALLW(SENDBUF, SENDCOUNTS, SDISPLS, SENDTYPES, BUFFER, RECVCOUNTS, RDISPLS, RECVTYPES, NEWCOMM, IERR)

          ! FREE RESOURCES
          DO P = 0, NTASKS - 1
            IF (SENDCOUNTS(P) .EQ. 1) THEN
              CALL MPI_TYPE_FREE(SENDTYPES(P), IERR)
              CALL MPI_TYPE_FREE(RECVTYPES(P), IERR)
            ENDIF
          ENDDO

          SENDCOUNTS = 0
          RECVCOUNTS = 0

          ! REPLACE BUFFER ON THE LHS WITH SPEC ON THE NEXT 3 BLOCKS
          ! TAKE COMPLEX CONJUGATE
          ! DELTA*, GAMMA*, PHI*, EPSILON*
          DO K = MAX(NYQUIST(3) + 1, FS(3)), FE(3)
            DO J = FS(2), FE(2)
              SPEC(I, J, K) = CONJG(BUFFER(J, K))
            ENDDO
          ENDDO

          ! BETA*
          DO K = MAX(1, FS(3)), 1
            DO J = MAX(NYQUIST(2) + 1, FS(2)), FE(2)
              SPEC(I, J, K) = CONJG(BUFFER(J, K))
            ENDDO
          ENDDO

          ! THETA*
          DO K = MAX(NYQUIST(3), FS(3)), MIN(NYQUIST(3), FE(3))
            DO J = MAX(NYQUIST(2) + 1, FS(2)), FE(2)
              SPEC(I, J, K) = CONJG(BUFFER(J, K))
            ENDDO
          ENDDO

        END SUBROUTINE EXCHANGE_CONJG

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
          SUBSIZES = [NX + 1, 1, NZ]        !< SHAPE OF HALO
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

    SUBROUTINE INTERPOLATE(PMAX, XYZ, DELTA, V)

      !

      INTEGER(IPP),                   INTENT(IN)  :: PMAX
      REAL(FPP),    DIMENSION(:,:),   INTENT(IN)  :: XYZ                          !< POINTS LOCATION
      REAL(FPP),    DIMENSION(:,:,:), INTENT(IN)  :: DELTA                            !< RANDOM FIELD
      REAL(FPP),    DIMENSION(:),     INTENT(OUT) :: V                            !< INTERPOLATED VALUES
      INTEGER(IPP)                                :: P                                !< COUNTERS
      INTEGER(IPP)                                :: I0, J0, K0                       !< INDICES
      REAL(FPP)                                   :: I, J, K
      REAL(FPP)                                   :: PX, PY, PZ, IPX, IPY             !< INTERPOLATION STUFF
      REAL(FPP)                                   :: A, B                             !< INTERPOLATION STUFF
      REAL(FPP),    DIMENSION(2,2)                :: F                                !< INTERPOLATION STUFF

      !-----------------------------------------------------------------------------------------------------------------------------

      ! LOOP OVER POINTS
      DO P = 1, PMAX

        I = XYZ(1, P)
        J = XYZ(2, P)
        K = XYZ(3, P)

        ! I0 = FLOOR(I)
        ! J0 = FLOOR(J)
        ! K0 = FLOOR(K)
        !
        ! F(:, 1) = DELTA(I0:I0 + 1, J0, K0)
        ! F(:, 2) = DELTA(I0:I0 + 1, J0 + 1, K0)
        !
        ! PX = I - I0
        ! PY = J - J0
        !
        ! IPX = (1._FPP - PX)
        ! IPY = (1._FPP - PY)
        !
        ! ! BILINEAR INTERPOLATION AT LEVEL 1
        ! A = F(1, 1) * IPX * IPY + F(2, 1) * PX * IPY + F(1, 2) * IPX * PY + F(2, 2) * PX * PY
        !
        ! F(:, 1) = DELTA(I0:I0 + 1, J0, K0 + 1)
        ! F(:, 2) = DELTA(I0:I0 + 1, J0 + 1, K0 + 1)
        !
        ! ! BILINEAR INTERPOLATION AT LEVEL 2
        ! B = F(1, 1) * IPX * IPY + F(2, 1) * PX * IPY + F(1, 2) * IPX * PY + F(2, 2) * PX * PY
        !
        ! PZ = K - K0
        !
        ! ! LINEAR INTERPOLATED BETWEEN LEVEL 1 AND 2
        ! V(P) = A * (1._FPP - PZ) + B * PZ

        I0 = NINT(I)
        J0 = NINT(J)
        K0 = NINT(K)

        V(P) = DELTA(I0, J0, K0)

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

      LIST(:,:) = 0

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

        MATCH = .FALSE.

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

      INTEGER(IPP), DIMENSION(3), INTENT(IN)  :: DIMS
      REAL(FPP),                  INTENT(OUT) :: MEASURE
      REAL(FPP),    DIMENSION(3)              :: SIDE

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

END MODULE SCARFLIB_FFT
