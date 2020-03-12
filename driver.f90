PROGRAM DRIVER

  ! driver program for scarf3d
  ! compile: mpif90 -O3 scarflib_fft2.f90 scarflib_spec.f90 scarflib.f90 driver.f90 *.f -I/home/walter/Backedup/Software/fftw-3.3.8/include
  !                 -L/home/walter/Backedup/Software/fftw-3.3.8/lib -lfftw3 -fcheck=all -fbacktrace -fopenacc
  !
  ! the product of "dims" below should be equal to the number of processes


  USE, INTRINSIC     :: ISO_FORTRAN_ENV
  USE, NON_INTRINSIC :: MPI
  USE, NON_INTRINSIC :: SCARFLIB
  USE, NON_INTRINSIC :: SCARFLIB_COMMON, ONLY: WATCH_START, WATCH_STOP
  USE, NON_INTRINSIC :: SCARFLIB_AUX

  IMPLICIT NONE

  INTEGER, PARAMETER :: IP = INT32
  !INTEGER, PARAMETER :: IP = INT64
  INTEGER(IP), PARAMETER :: FP = REAL32
  !INTEGER(IP), PARAMETER :: FP = REAL64

  INTEGER(IP)                               :: I, J, K
  INTEGER(IP)                               :: RANK, NTASKS, IERR
  INTEGER(IP)                               :: METHOD, ACF, RESCALE, PAD, SEED
  INTEGER(IP),    DIMENSION(3)              :: N, FS, FE
  REAL(FP)                                  :: DS, DH, SIGMA, HURST, MUTE, TAPER
  REAL(FP),       DIMENSION(3)              :: CL
  REAL(REAL64)                              :: TICTOC
  REAL(FP),       DIMENSION(3,2)            :: POI
  REAL(FP),       DIMENSION(:,:,:), POINTER :: X3, Y3, Z3, V3
  REAL(FP),       ALLOCATABLE, DIMENSION(:), TARGET :: X1, Y1, Z1, V1
  TYPE(SCARF_OBJ)                           :: PARAMS

  !--------------------------------------------------------------------------------------------------------------------------------

  ! INITIALISE MPI
  CALL MPI_INIT(IERR)

  ! GET RANK NUMBER
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, RANK, IERR)

  ! GET NUMBER OF TASKS
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NTASKS, IERR)

  CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

  !---------------------------------------------------------------------------------------------------------------------------------
  ! INPUT SECTION
  !---------------------------------------------------------------------------------------------------------------------------------

  ! FFT = 0, SPEC = 1
  METHOD = 0

  ! NUMBER OF POINTS FOR WHOLE MODEL
  !N = [2000*2, 3200*2, 1200*2]
  N = [500, 500, 500]

  ! GRID STEP
  DS = 100._FP

  DH = DS

  ! AUTOCORRELATION (0=VON KARMAN/EXPONENTIAL, 1=GAUSSIAN)
  ACF = 0

  ! CORRELATION LENGTH
  CL = [5000._FP, 5000._FP, 5000._FP]

  ! STANDARD DEVIATION (SIGMA%/100)
  SIGMA = 0.05_FP

  ! HURST EXPONENT (NOT USED FOR GAUSSIAN ACF)
  HURST = 0.5_FP

  ! SEED NUMBER
  SEED = 1235

  ! SET POSITION OF POINT-OF-INTEREST (MUTING/TAPERING), SAME UNITS AS "DH"
  POI(:, 1) = [400., 250., 100.] * DS
  POI(:, 2) = [350., 200.,  50.] * DS

  ! RADIUS FOR MUTING (AT POI), SAME UNITS AS "DH"
  MUTE = 1000. !1000.

  ! RADIUS FOR TAPERING (AT POI + MUTE), SAME UNITS AS "DH"
  TAPER = 5000. !5000.

  ! RESCALE TO DESIRED (CONTINUOUS) SIGMA
  RESCALE = 0

  ! EXPAND GRID TO HANDLE FFT PERIODICITY
  PAD = 0

  ! END INPUT SECTION
  !-----------------------------------------------------------------------

  ! STEPS BELOW:
  !
  ! A) CREATE SAMPLE STRUCTURED MESH, DESCRIBED BY CARTESIAN TOPOLOGY, REPRESENTING A PLAUSIBLE 3D DOMAIN AS IN, E.G., FINITE-DIFFERENCE
  ! B) COMPUTE THE RANDOM FIELD
  ! C) WRITE THE RANDOM FIELD TO DISK

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! STEP A)





  CALL SAMPLE_MESH(RANK, NTASKS, N, FS, FE, X1, Y1, Z1, V1)

  X3(FS(1):FE(1), FS(2):FE(2), FS(3):FE(3)) => X1
  Y3(FS(1):FE(1), FS(2):FE(2), FS(3):FE(3)) => Y1
  Z3(FS(1):FE(1), FS(2):FE(2), FS(3):FE(3)) => Z1
  V3(FS(1):FE(1), FS(2):FE(2), FS(3):FE(3)) => V1

  DO K = FS(3), FE(3)
    DO J = FS(2), FE(2)
      DO I = FS(1), FE(1)
        X3(I, J, K) = (I - 0.5) * DS
        Y3(I, J, K) = (J - 0.5) * DS
        Z3(I, J, K) = (K - 0.5) * DS
      ENDDO
    ENDDO
  ENDDO


  ! STRUCTURED MESH TEST
  PARAMS = SCARF_INITIALIZE(FS, FE, DS, ACF, CL, SIGMA, METHOD, HURST)

  ! PARAMS = SCARF_INITIALIZE(FS, FE, DS, ACF, CL, SIGMA, METHOD = 0, HURST, POI = POI, MUTE = MUTE, TAPER = TAPER, PAD = 1)
  ! PARAMS = SCARF_INITIALIZE(FS, FE, DS, ACF, CL, SIGMA, METHOD = 0, HURST, RESCALE = 1, PAD = 1)

  CALL WATCH_START(TICTOC)

  CALL SCARF_EXECUTE(PARAMS, SEED, V3)

  CALL WATCH_STOP(TICTOC)

  IF (RANK .EQ. 0) PRINT*, 'STRUCTURED MESH TEST COMPLETED IN: ', REAL(TICTOC, REAL32)

  IF ((RANK .EQ. 0) .AND. (METHOD .EQ. 0)) THEN
    PRINT*, 'DOMAIN TOO SMALL?', NINT(PARAMS%STATS(1))
    PRINT*, 'DS TOO LARGE?',     NINT(PARAMS%STATS(2))
    PRINT*, 'STAND. DEV.: ',     REAL(PARAMS%STATS(3), KIND=REAL32)
    PRINT*, 'MEAN: ',            REAL(PARAMS%STATS(4), KIND=REAL32)
    PRINT*, 'TIMING SPECTRUM: ', REAL(PARAMS%STATS(5), KIND=REAL32)
    PRINT*, 'TIMING SYMMETRY: ', REAL(PARAMS%STATS(6), KIND=REAL32)
    PRINT*, 'TIMING IFFT: ',     REAL(PARAMS%STATS(7), KIND=REAL32)
    PRINT*, 'TIMING INTERP: ',   REAL(PARAMS%STATS(8), KIND=REAL32)
  ENDIF

  IF ((RANK .EQ. 0) .AND. (METHOD .EQ. 1)) THEN
    PRINT*, 'DOMAIN TOO SMALL?', NINT(PARAMS%STATS(1))
    PRINT*, 'DS TOO LARGE?',     NINT(PARAMS%STATS(2))
    PRINT*, 'STAND. DEV.: ',     REAL(PARAMS%STATS(3), KIND=REAL32)
    PRINT*, 'MEAN: ',            REAL(PARAMS%STATS(4), KIND=REAL32)
    PRINT*, 'TIMING CPU: ',      REAL(PARAMS%STATS(5), KIND=REAL32)
    PRINT*, 'TIMING GPU: ',      REAL(PARAMS%STATS(6), KIND=REAL32)
  ENDIF

  CALL WATCH_START(TICTOC)

  CALL SCARF_IO(PARAMS, 1, 400, V3, 'fft_struct_xslice')
  CALL SCARF_IO(PARAMS, 2, 250, V3, 'fft_struct_yslice')
  CALL SCARF_IO(PARAMS, 3, 100, V3, 'fft_struct_zslice')

  CALL WATCH_STOP(TICTOC)

  IF (RANK .EQ. 0) PRINT*, 'SLICE(S) WRITTEN IN: ', REAL(TICTOC, REAL32)

  CALL WATCH_START(TICTOC)

  CALL SCARF_IO(PARAMS, V3, 'fft_struct_whole', 3)

  CALL WATCH_STOP(TICTOC)

  IF (RANK .EQ. 0) PRINT*, 'WHOLE FILE WRITTEN IN: ', REAL(TICTOC, REAL32)

  CALL SCARF_FINALIZE(PARAMS)

  ! UNSTRUCTURED MESH TES
  PARAMS = SCARF_INITIALIZE(X1, Y1, Z1, DH, ACF, CL, SIGMA, METHOD, HURST)
  ! PARAMS = SCARF_INITIALIZE(X1, Y1, Z1, DH, ACF, CL, SIGMA, METHOD, HURST, POI = POI, MUTE = MUTE, TAPER = TAPER, PAD = 1)
  ! PARAMS = SCARF_INITIALIZE(X1, Y1, Z1, DH, ACF, CL, SIGMA, METHOD, HURST, RESCALE = 1, PAD = 1)

  CALL WATCH_START(TICTOC)

  CALL SCARF_EXECUTE(PARAMS, SEED, V1)

  CALL WATCH_STOP(TICTOC)

  IF (RANK .EQ. 0) PRINT*, 'UNSTRUCTURED MESH TEST COMPLETED IN: ', REAL(TICTOC, REAL32)

  IF ((RANK .EQ. 0) .AND. (METHOD .EQ. 0)) THEN
    PRINT*, 'DOMAIN TOO SMALL?', NINT(PARAMS%STATS(1))
    PRINT*, 'DS TOO LARGE?',     NINT(PARAMS%STATS(2))
    PRINT*, 'STAND. DEV.: ',     REAL(PARAMS%STATS(3), KIND=REAL32)
    PRINT*, 'MEAN: ',            REAL(PARAMS%STATS(4), KIND=REAL32)
    PRINT*, 'TIMING SPECTRUM: ', REAL(PARAMS%STATS(5), KIND=REAL32)
    PRINT*, 'TIMING SYMMETRY: ', REAL(PARAMS%STATS(6), KIND=REAL32)
    PRINT*, 'TIMING IFFT: ',     REAL(PARAMS%STATS(7), KIND=REAL32)
    PRINT*, 'TIMING INTERP: ',   REAL(PARAMS%STATS(8), KIND=REAL32)
  ENDIF

  IF ((RANK .EQ. 0) .AND. (METHOD .EQ. 1)) THEN
    PRINT*, 'DOMAIN TOO SMALL?', NINT(PARAMS%STATS(1))
    PRINT*, 'DS TOO LARGE?',     NINT(PARAMS%STATS(2))
    PRINT*, 'STAND. DEV.: ',     REAL(PARAMS%STATS(3), KIND=REAL32)
    PRINT*, 'MEAN: ',            REAL(PARAMS%STATS(4), KIND=REAL32)
    PRINT*, 'TIMING CPU: ',      REAL(PARAMS%STATS(5), KIND=REAL32)
    PRINT*, 'TIMING GPU: ',      REAL(PARAMS%STATS(6), KIND=REAL32)
  ENDIF

  CALL WATCH_START(TICTOC)

  PARAMS%FS = FS
  PARAMS%FE = FE

  CALL SCARF_IO(PARAMS, 1, 400, V3, 'fft_unstruct_xslice')
  CALL SCARF_IO(PARAMS, 2, 250, V3, 'fft_unstruct_yslice')
  CALL SCARF_IO(PARAMS, 3, 100, V3, 'fft_unstruct_zslice')

  CALL WATCH_STOP(TICTOC)

  IF (RANK .EQ. 0) PRINT*, 'SLICE(S) WRITTEN IN: ', REAL(TICTOC, REAL32)

  CALL WATCH_START(TICTOC)

  CALL SCARF_IO(PARAMS, V3, 'fft_unstruct_whole', 3)

  CALL WATCH_STOP(TICTOC)

  IF (RANK .EQ. 0) PRINT*, 'WHOLE FILE WRITTEN IN: ', REAL(TICTOC, REAL32)

  CALL SCARF_FINALIZE(PARAMS)

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
