PROGRAM DRIVER

  ! driver program for scarf3d
  ! compile: mpif90 -O3 scarflib_fft2.f90 scarflib_spec.f90 scarflib.f90 driver.f90 *.f -I/home/walter/Backedup/Software/fftw-3.3.8/include
  !                 -L/home/walter/Backedup/Software/fftw-3.3.8/lib -lfftw3 -fcheck=all -fbacktrace -fopenacc
  !
  ! the product of "dims" below should be equal to the number of processes


  !USE :: MPI

  USE, NON_INTRINSIC :: SCARFLIB_COMMON
  USE, NON_INTRINSIC :: SCARFLIB_SPECTRAL
  USE, NON_INTRINSIC :: SCARFLIB_FFT3

  IMPLICIT NONE

  INTEGER(IPP)                                         :: I, J, K
  INTEGER(IPP)                                         :: IERR, RANK, NTASKS, TOPO, NDIMS
  INTEGER(IPP)                                         :: SEED, RESCALE, PAD, ACF
  INTEGER(IPP),              DIMENSION(3)              :: N, FS, FE, COORDS, DIMS
  LOGICAL                                              :: REORDER
  LOGICAL,                   DIMENSION(3)              :: ISPERIODIC
  REAL(FPP)                                            :: DH, DR, SIGMA, HURST
  REAL(FPP)                                            :: MUTE, TAPERING
  REAL(REAL64)                                         :: TICTOC
  REAL(FPP),                 DIMENSION(3)              :: CL
  REAL(FPP),                 DIMENSION(8)              :: INFO
  REAL(FPP),                 DIMENSION(3,2)            :: POI
  REAL(FPP),                 DIMENSION(:,:,:), POINTER :: X3, Y3, Z3, V3
  REAL(FPP),    ALLOCATABLE, DIMENSION(:),     TARGET  :: X1, Y1, Z1, V1

  !--------------------------------------------------------------------------------------------------------------------------------

  ! INITIALISE MPI
  CALL MPI_INIT(IERR)

  ! GET RANK NUMBER
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, RANK, IERR)

  ! GET NUMBER OF TASKS
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NTASKS, IERR)

  CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

  !-----------------------------------------------------------------------
  ! INPUT SECTION

  ! NUMBER OF POINTS FOR WHOLE MODEL
  !N = [2000*2, 3200*2, 1200*2]
  N = [600, 500, 400]

  ! GRID STEP FOR FFT GRID
  DH = 100._FPP

  ! GRID STEP FOR POINTS
  DR = 100._FPP

  ! AUTOCORRELATION (0=VON KARMAN/EXPONENTIAL, 1=GAUSSIAN)
  ACF = 0

  ! CORRELATION LENGTH
  CL = [5000._FPP, 5000._FPP, 5000._FPP]

  ! STANDARD DEVIATION (SIGMA%/100)
  SIGMA = 0.05_FPP

  ! HURST EXPONENT (NOT USED FOR GAUSSIAN ACF)
  HURST = 0.5_FPP

  ! SEED NUMBER
  SEED = 1235

  ! SET POSITION OF POINT-OF-INTEREST (MUTING/TAPERING)
  POI(:, 1) = [200000., 100000., 320000.]
  POI(:, 2) = [300000., 50000., 320000.]

  ! RADIUS FOR MUTING (AT POI)
  MUTE = 200.

  ! RADIUS FOR TAPERING (AT POI + MUTE)
  TAPERING = 5000.

  ! RESCALE TO DESIRED (CONTINUOUS) SIGMA
  RESCALE = 0

  ! EXPAND GRID TO HANDLE FFT PERIODICITY
  PAD = 1

  ! END INPUT SECTION
  !-----------------------------------------------------------------------

  ! STEPS BELOW:
  !
  ! A) CREATE SAMPLE STRUCTURED MESH, DESCRIBED BY CARTESIAN TOPOLOGY, REPRESENTING A PLAUSIBLE 3D DOMAIN AS IN, E.G., FINITE-DIFFERENCE
  ! B) COMPUTE THE RANDOM FIELD
  ! C) WRITE THE RANDOM FIELD TO DISK

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! STEP A)

  ! WE WORK IN 3D
  NDIMS = 3

  ! NUMBER OF PROCESSES ALONG EACH DIRECTION
  !DIMS  = [10*2, 16*2, 6*2]
  !DIMS = [2, 1, 2]

  WORLD_SIZE = NTASKS
  NPTS       = N

  CALL BEST_CONFIG(DIMS)

  IF (RANK == 0) PRINT*, 'CPU GRID IN DRIVER: ', DIMS

  ! ALLOW REORDERING
  REORDER = .TRUE.

  ! NO PERIODICITY
  ISPERIODIC = [.FALSE., .FALSE., .FALSE.]

  ! CREATE TOPOLOGY
  CALL MPI_CART_CREATE(MPI_COMM_WORLD, NDIMS, DIMS, ISPERIODIC, REORDER, TOPO, IERR)

  CALL MPI_CART_COORDS(TOPO, RANK, NDIMS, COORDS, IERR)

  ! SPLIT EACH AXIS AND ASSIGN POINTS TO CURRENT PROCESS
  CALL MPI_SPLIT_TASK(N, DIMS, COORDS, FS, FE)

  ! DO I = 0, NTASKS - 1
  !   IF (RANK .EQ. I) THEN
  !     PRINT*, 'INPUT MESH: ', RANK, FS(1), FE(1), ' -- ', FS(2), FE(2), ' -- ', FS(3), FE(3)
  !   ENDIF
  !   CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
  ! ENDDO

  ! PREPARE ARRAY FOR RANDOM FIELD
  ALLOCATE(V1((FE(1) - FS(1) + 1) * (FE(2) - FS(2) + 1) * (FE(3) - FS(3) + 1)))

  ! MAP 1D ARRAY INTO 3D ONE
  V3(FS(1):FE(1), FS(2):FE(2), FS(3):FE(3)) => V1

  ! PREPARE ARRAY WITH MODEL POINTS POSITION
  ALLOCATE(X1((FE(1) - FS(1) + 1) * (FE(2) - FS(2) + 1) * (FE(3) - FS(3) + 1)))
  ALLOCATE(Y1((FE(1) - FS(1) + 1) * (FE(2) - FS(2) + 1) * (FE(3) - FS(3) + 1)))
  ALLOCATE(Z1((FE(1) - FS(1) + 1) * (FE(2) - FS(2) + 1) * (FE(3) - FS(3) + 1)))

  X3(FS(1):FE(1), FS(2):FE(2), FS(3):FE(3)) => X1
  Y3(FS(1):FE(1), FS(2):FE(2), FS(3):FE(3)) => Y1
  Z3(FS(1):FE(1), FS(2):FE(2), FS(3):FE(3)) => Z1

  DO K = FS(3), FE(3)
    DO J = FS(2), FE(2)
      DO I = FS(1), FE(1)
        X3(I, J, K) = (I - 0.5) * DR
        Y3(I, J, K) = (J - 0.5) * DR
        Z3(I, J, K) = (K - 0.5) * DR
      ENDDO
    ENDDO
  ENDDO

  ! COMPUTE ZERO-MEAN RANDOM FIELD
  !CALL SCARF3D_FFT(MPI_COMM_WORLD, [NX, NY, NZ], 100._fpp, 'VK', [5000._fpp, 5000._fpp, 5000._fpp], 0.05_fpp, 0.1_fpp, 1235, POI, 4, 26)
  !CALL SCARF3D_FFT(MPI_COMM_WORLD, N, 100._fpp, 'VK', [5000._fpp, 5000._fpp, 5000._fpp], 0.05_fpp, 0.5_fpp, 1235, POI, 0, 0)

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! STEP B + C: GPU

  IF (RANK .EQ. 0) PRINT*, 'GPU'

  CALL WATCH_START(TICTOC)

  !CALL SCARF3D_SPEC(FS, FE, X, Y, Z, ACF, CL, SIGMA, HURST, SEED, POI, 0, 0, V, TICTOC)
  !CALL SCARF3D_SPEC(X1, Y1, Z1, DH, ACF, CL, SIGMA, HURST, SEED, POI, MUTE, TAPERING, V1, INFO(1))

  CALL WATCH_STOP(TICTOC)

  IF (RANK .EQ. 0) PRINT*, 'EXEC TIME GPU:', REAL(TICTOC, REAL32)

  ! ARRAYS WHERE FIRST AND LAST INDICES FOR EACH PROCESS ARE STORED. THESE WILL BE USED TO WRITE FIELD SLICES TO DISK.
  ! ALLOCATE(GS(3, 0:NTASKS-1), GE(3, 0:NTASKS-1))
  !
  ! ! STORE GLOBAL INDICES
  ! GS(:, RANK) = FS
  ! GE(:, RANK) = FE
  !
  ! ! MAKE ALL PROCESSES AWARE OF GLOBAL INDICES ALONG EACH AXIS
  ! CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, GS, 3, MPI_INTEGER, MPI_COMM_WORLD, IERR)
  ! CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, GE, 3, MPI_INTEGER, MPI_COMM_WORLD, IERR)

  ! WRITE A MODEL SLICE TO DISK
  !CALL IO_WRITE_SLICE(2, 1, V3, 'gpu_slice.bin')

  ! DEALLOCATE(GS, GE)

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! STEP B + C: FFT

  V1(:) = 0._FPP

  IF (RANK .EQ. 0) PRINT*, 'FFT'

  CALL WATCH_START(TICTOC)

  !CALL SCARF3D_FFT(X1, Y1, Z1, DH, ACF, CL, SIGMA, HURST, SEED, POI, MUTE, TAPERING, RESCALE, PAD, V1, INFO)
  CALL SCARF3D_FFT(DR, FS, FE, DH, ACF, CL, SIGMA, HURST, SEED, POI, MUTE, TAPERING, RESCALE, PAD, V3, INFO)

  CALL WATCH_STOP(TICTOC)

  IF (RANK .EQ. 0) PRINT*, 'EXEC TIME FFT:', REAL(TICTOC, REAL32)

  ! ARRAYS WHERE FIRST AND LAST INDICES FOR EACH PROCESS ARE STORED. THESE WILL BE USED TO WRITE FIELD SLICES TO DISK.
  ALLOCATE(GS(3, 0:NTASKS-1), GE(3, 0:NTASKS-1))

  ! STORE GLOBAL INDICES
  GS(:, RANK) = FS
  GE(:, RANK) = FE

  NPTS = N

  ! MAKE ALL PROCESSES AWARE OF GLOBAL INDICES ALONG EACH AXIS
  CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, GS, 3, MPI_INTEGER, MPI_COMM_WORLD, IERR)
  CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, GE, 3, MPI_INTEGER, MPI_COMM_WORLD, IERR)

  CALL WATCH_START(TICTOC)

  ! WRITE A MODEL SLICE TO DISK
  CALL SCARF3D_WRITE_SLICE(1, 100, V3, 'fft_xslice')
  CALL SCARF3D_WRITE_SLICE(2, 100, V3, 'fft_yslice')
  CALL SCARF3D_WRITE_SLICE(3, 100, V3, 'fft_zslice')

  CALL WATCH_STOP(TICTOC)

  IF (RANK .EQ. 0) PRINT*, 'EXEC SLICE IO:', REAL(TICTOC, REAL32)

  CALL WATCH_START(TICTOC)

  CALL SCARF3D_WRITE_ONE(V3, 'fft_whole', 3)
  !CALL SCARF3D_WRITE_SINGLE(V3, 'fft_single', 3)

  CALL WATCH_STOP(TICTOC)

  IF (RANK .EQ. 0) PRINT*, 'EXEC WHOLE IO:', REAL(TICTOC, REAL32)

  DEALLOCATE(GS, GE)

  ! DO I = 0, NTASKS - 1
  !   IF (RANK .EQ. I) PRINT*, RANK, MINVAL(V3), MAXVAL(V3)
  !   CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
  ! ENDDO

  ! SOME TIMING INFO
  IF (RANK .EQ. 0) THEN
    PRINT*, 'DOMAIN TOO SMALL?',   NINT(INFO(1))
    PRINT*, 'DH TOO LARGE?',       NINT(INFO(2))
    PRINT*, 'STAND. DEV.: ',       REAL(INFO(3), KIND=REAL32)
    PRINT*, 'MEAN: ',              REAL(INFO(4), KIND=REAL32)
    PRINT*, 'TIMING SPECTRUM: ',   REAL(INFO(5), KIND=REAL32)
    PRINT*, 'TIMING SYMMETRY: ',   REAL(INFO(6), KIND=REAL32)
    PRINT*, 'TIMING IFFT: ',       REAL(INFO(7), KIND=REAL32)
    PRINT*, 'TIMING INTERP: ',     REAL(INFO(8), KIND=REAL32)

  ENDIF

  NULLIFY(V3, X3, Y3, Z3)
  DEALLOCATE(V1, X1, Y1, Z1)

  CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

  CALL MPI_FINALIZE(IERR)

END PROGRAM DRIVER

! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
!===============================================================================================================================
! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

SUBROUTINE MPI_SPLIT_TASK(N, NTASKS, RANK, FS, FE)

  USE, NON_INTRINSIC :: SCARFLIB_COMMON

  IMPLICIT NONE

  INTEGER(IPP), DIMENSION(3), INTENT(IN)  :: N                          !< NUMBER OF POINTS TO BE SPLIT
  INTEGER(IPP), DIMENSION(3), INTENT(IN)  :: NTASKS                     !< NUMBER OF MPI PROCESSES
  INTEGER(IPP), DIMENSION(3), INTENT(IN)  :: RANK                       !< RANK OF CURRENT PREOCESS
  INTEGER(IPP), DIMENSION(3), INTENT(OUT) :: FS, FE                     !< 1ST/LAST INDICES
  INTEGER(IPP)                            :: I                          !< COUNTER

  !------------------------------------------------------------------------------------------------------------------------------

  DO I = 1, 3
    FS(I) = 1 + INT( REAL(N(I), FPP) / REAL(NTASKS(I), FPP) * REAL(RANK(I), FPP) )
    FE(I) = INT( REAL(N(I), FPP) / REAL(NTASKS(I), FPP) * REAL(RANK(I) + 1, FPP) )
  ENDDO

END SUBROUTINE MPI_SPLIT_TASK
