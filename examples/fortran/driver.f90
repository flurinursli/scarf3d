PROGRAM DRIVER

  ! driver program for scarf3d
  ! compile: mpif90 -O3 scarflib_fft2.f90 scarflib_spec.f90 scarflib.f90 driver.f90 *.f -I/home/walter/Backedup/Software/fftw-3.3.8/include
  !                 -L/home/walter/Backedup/Software/fftw-3.3.8/lib -lfftw3 -fcheck=all -fbacktrace -fopenacc
  !
  ! the product of "dims" below should be equal to the number of processes


  USE, INTRINSIC     :: ISO_FORTRAN_ENV
  USE, NON_INTRINSIC :: MPI
  USE, NON_INTRINSIC :: M_SCARFLIB
  USE, NON_INTRINSIC :: SCARFLIB_AUX

  IMPLICIT NONE

#ifdef DOUBLE_PREC
  INTEGER,     PARAMETER :: IP = INT32
  INTEGER(IP), PARAMETER :: FP = REAL64
#else
  INTEGER,     PARAMETER :: IP = INT32
  INTEGER(IP), PARAMETER :: FP = REAL32
#endif

  INTEGER(IP)                                         :: I, J, K
  INTEGER(IP)                                         :: RANK, NTASKS, IERR
  INTEGER(IP)                                         :: ACF, RESCALE, PAD, SEED
  INTEGER(IP),              DIMENSION(3)              :: N, FS, FE
  REAL(FP)                                            :: DS, DH, SIGMA, HURST, MUTE, TAPER
  REAL(FP),                 DIMENSION(3)              :: CL
  REAL(FP),                 DIMENSION(8)              :: STATS
  REAL(REAL64)                                        :: TICTOC
  REAL(FP),                 DIMENSION(3,2)            :: POI
  REAL(FP),                 DIMENSION(:,:,:), POINTER :: X3, Y3, Z3, V3
  REAL(FP),    ALLOCATABLE, DIMENSION(:),     TARGET  :: X1, Y1, Z1, V1

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

  ! NUMBER OF POINTS FOR WHOLE MODEL
  !N = [2000*2, 3200*2, 1200*2]
  N = [500, 500, 500]

  ! GRID STEP
  DS = 50._FP

  DH = 50._FP

  ! AUTOCORRELATION (0=VON KARMAN/EXPONENTIAL, 1=GAUSSIAN)
  ACF = 0

  ! CORRELATION LENGTH
  CL = [2000._FP, 2000._FP, 2000._FP]

  ! STANDARD DEVIATION (SIGMA%/100)
  SIGMA = 0.05_FP

  ! HURST EXPONENT (NOT USED FOR GAUSSIAN ACF)
  HURST = 0.25_FP

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
  !---------------------------------------------------------------------------------------------------------------------------------

  ! ================================================================================================================================-
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! CREATE SAMPLE STRUCTURED MESH
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! ================================================================================================================================-

  CALL SAMPLE_MESH(RANK, NTASKS, N, FS, FE)

  ! ALLOCATE MEMORY FOR MESH
  ALLOCATE(V1((FE(1) - FS(1) + 1) * (FE(2) - FS(2) + 1) * (FE(3) - FS(3) + 1)))
  ALLOCATE(X1((FE(1) - FS(1) + 1) * (FE(2) - FS(2) + 1) * (FE(3) - FS(3) + 1)))
  ALLOCATE(Y1((FE(1) - FS(1) + 1) * (FE(2) - FS(2) + 1) * (FE(3) - FS(3) + 1)))
  ALLOCATE(Z1((FE(1) - FS(1) + 1) * (FE(2) - FS(2) + 1) * (FE(3) - FS(3) + 1)))

  X3(FS(1):FE(1), FS(2):FE(2), FS(3):FE(3)) => X1
  Y3(FS(1):FE(1), FS(2):FE(2), FS(3):FE(3)) => Y1
  Z3(FS(1):FE(1), FS(2):FE(2), FS(3):FE(3)) => Z1
  V3(FS(1):FE(1), FS(2):FE(2), FS(3):FE(3)) => V1

  DO K = FS(3), FE(3)
    DO J = FS(2), FE(2)
      DO I = FS(1), FE(1)
        X3(I, J, K) = (I - 1) * DS
        Y3(I, J, K) = (J - 1) * DS
        Z3(I, J, K) = (K - 1) * DS
      ENDDO
    ENDDO
  ENDDO

  ! ================================================================================================================================-
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! FFT METHOD TESTS
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! ================================================================================================================================-

  IF (RANK .EQ. 0) PRINT*, ''
  IF (RANK .EQ. 0) PRINT*, '**********************************************'
  IF (RANK .EQ. 0) PRINT*, '***************** FFT METHOD *****************'

  ! STRUCTURED MESH TEST
  CALL SCARF_INITIALIZE(FS, FE, DS, ACF, CL, SIGMA, METHOD = 0, HURST = HURST, PAD = 0, DH = DH)
  !PARAMS = SCARF_INITIALIZE(FS, FE, DS, ACF, CL, SIGMA, METHOD = 0, HURST = HURST, PAD = 1)
  ! PARAMS = SCARF_INITIALIZE(FS, FE, DS, ACF, CL, SIGMA, METHOD = 0, HURST = HURST, POI = POI, MUTE = MUTE, TAPER = TAPER, PAD = 1)
  ! PARAMS = SCARF_INITIALIZE(FS, FE, DS, ACF, CL, SIGMA, METHOD = 0, HURST = HURST, RESCALE = 1, PAD = 1)

  CALL WATCH_START(TICTOC)

  CALL SCARF_EXECUTE(SEED, V3, STATS)

  CALL WATCH_STOP(TICTOC)

  CALL MPI_ALLREDUCE(MPI_IN_PLACE, STATS, 8, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, IERR)

  IF (RANK .EQ. 0) THEN
    PRINT*, ''
    PRINT*, 'STRUCTURED MESH TEST COMPLETED IN: ', REAL(TICTOC, REAL32)
    PRINT*, 'DOMAIN TOO SMALL?                : ', NINT(STATS(1))
    PRINT*, 'GRID-STEP TOO LARGE?             : ', NINT(STATS(2))
    PRINT*, 'STANDARD DEVIATION               : ', REAL(STATS(3), KIND=REAL32)
    PRINT*, 'MEAN VALUE                       : ', REAL(STATS(4), KIND=REAL32)
    PRINT*, 'TIMING FOR SPECTRUM              : ', REAL(STATS(5), KIND=REAL32)
    PRINT*, 'TIMING FOR SYMMETRY              : ', REAL(STATS(6), KIND=REAL32)
    PRINT*, 'TIMING FOR IFFT                  : ', REAL(STATS(7), KIND=REAL32)
    PRINT*, 'TIMING FOR INTERPOLATION         : ', REAL(STATS(8), KIND=REAL32)
  ENDIF

  CALL WATCH_START(TICTOC)

  CALL SCARF_IO(N, 'x', N(1)/2, V3, 'fft_struct_xslice')
  CALL SCARF_IO(N, 'y', N(2)/2, V3, 'fft_struct_yslice')
  CALL SCARF_IO(N, 'z', N(3)/2, V3, 'fft_struct_zslice')

  CALL WATCH_STOP(TICTOC)

  IF (RANK .EQ. 0) PRINT*, 'SLICE(S) WRITTEN IN              : ', REAL(TICTOC, REAL32)

  CALL WATCH_START(TICTOC)

  CALL SCARF_IO(N, V3, 'fft_struct_whole', 3)

  CALL WATCH_STOP(TICTOC)

  IF (RANK .EQ. 0) PRINT*, 'WHOLE FILE WRITTEN IN            : ', REAL(TICTOC, REAL32)

  CALL SCARF_FINALIZE()

  ! UNSTRUCTURED MESH TEST
  CALL SCARF_INITIALIZE(X1, Y1, Z1, DH, ACF, CL, SIGMA, METHOD = 0, HURST = HURST, PAD = 0, DS = DS)
  !PARAMS = SCARF_INITIALIZE(X1, Y1, Z1, DH, ACF, CL, SIGMA, METHOD = 0, HURST = HURST, PAD = 1)
  ! PARAMS = SCARF_INITIALIZE(X1, Y1, Z1, DH, ACF, CL, SIGMA, METHOD, HURST, POI = POI, MUTE = MUTE, TAPER = TAPER, PAD = 1)
  ! PARAMS = SCARF_INITIALIZE(X1, Y1, Z1, DH, ACF, CL, SIGMA, METHOD, HURST, RESCALE = 1, PAD = 1)

  CALL WATCH_START(TICTOC)

  CALL SCARF_EXECUTE(SEED, V1, STATS)

  CALL WATCH_STOP(TICTOC)

  CALL MPI_ALLREDUCE(MPI_IN_PLACE, STATS, 8, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, IERR)

  IF (RANK .EQ. 0) THEN
    PRINT*, ''
    PRINT*, 'UNSTRUCTURED MESH TEST COMPLETED IN: ', REAL(TICTOC, REAL32)
    PRINT*, 'DOMAIN TOO SMALL?                  : ', NINT(STATS(1))
    PRINT*, 'GRID-STEP TOO LARGE?               : ', NINT(STATS(2))
    PRINT*, 'STANDARD DEVIATION                 : ', REAL(STATS(3), KIND=REAL32)
    PRINT*, 'MEAN VALUE                         : ', REAL(STATS(4), KIND=REAL32)
    PRINT*, 'TIMING FOR SPECTRUM                : ', REAL(STATS(5), KIND=REAL32)
    PRINT*, 'TIMING FOR SYMMETRY                : ', REAL(STATS(6), KIND=REAL32)
    PRINT*, 'TIMING FOR IFFT                    : ', REAL(STATS(7), KIND=REAL32)
    PRINT*, 'TIMING FOR INTERPOLATION           : ', REAL(STATS(8), KIND=REAL32)
  ENDIF

  CALL WATCH_START(TICTOC)

  ! PARAMS%FS = FS
  ! PARAMS%FE = FE

  ! CALL SCARF_IO(N, 1, 400, V3, 'fft_unstruct_xslice')
  ! CALL SCARF_IO(N, 2, 250, V3, 'fft_unstruct_yslice')
  ! CALL SCARF_IO(N, 3, 100, V3, 'fft_unstruct_zslice')
  !
  ! CALL WATCH_STOP(TICTOC)
  !
  ! IF (RANK .EQ. 0) PRINT*, 'SLICE(S) WRITTEN IN                : ', REAL(TICTOC, REAL32)
  !
  ! CALL WATCH_START(TICTOC)
  !
  ! CALL SCARF_IO(N, V3, 'fft_unstruct_whole', 3)
  !
  ! CALL WATCH_STOP(TICTOC)
  !
  ! IF (RANK .EQ. 0) PRINT*, 'WHOLE FILE WRITTEN IN              : ', REAL(TICTOC, REAL32)

  CALL SCARF_FINALIZE()

  ! ================================================================================================================================-
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! SPECTRAL METHOD TESTS
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! ================================================================================================================================-

#ifdef SPECTRAL

  IF (RANK .EQ. 0) PRINT*, ''
  IF (RANK .EQ. 0) PRINT*, '**********************************************'
  IF (RANK .EQ. 0) PRINT*, '************** SPECTRAL METHOD ***************'

  ! STRUCTURED MESH TEST
  CALL SCARF_INITIALIZE(FS, FE, DS, ACF, CL, SIGMA, METHOD = 1, HURST = HURST, DH = DH)
  !PARAMS = SCARF_INITIALIZE(FS, FE, DS, ACF, CL, SIGMA, METHOD = 1, HURST = HURST)
  ! PARAMS = SCARF_INITIALIZE(FS, FE, DS, ACF, CL, SIGMA, METHOD = 1, HURST = HURST, POI = POI, MUTE = MUTE, TAPER = TAPER)
  ! PARAMS = SCARF_INITIALIZE(FS, FE, DS, ACF, CL, SIGMA, METHOD = 1, HURST = HURST)

  CALL WATCH_START(TICTOC)

  CALL SCARF_EXECUTE(SEED, V3, STATS)

  CALL WATCH_STOP(TICTOC)

  CALL MPI_ALLREDUCE(MPI_IN_PLACE, STATS, 8, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, IERR)

  IF (RANK .EQ. 0) THEN
    PRINT*, ''
    PRINT*, 'STRUCTURED MESH TEST COMPLETED IN: ', REAL(TICTOC, REAL32)
    PRINT*, 'DOMAIN TOO SMALL?                : ', NINT(STATS(1))
    PRINT*, 'GRID-STEP TOO LARGE?             : ', NINT(STATS(2))
    PRINT*, 'STANDARD DEVIATION               : ', REAL(STATS(3), KIND=REAL32)
    PRINT*, 'MEAN VALUE                       : ', REAL(STATS(4), KIND=REAL32)
    PRINT*, 'TIMING FOR CPU EXECUTION         : ', REAL(STATS(5), KIND=REAL32)
    PRINT*, 'TIMING FOR GPU EXECUTION         : ', REAL(STATS(6), KIND=REAL32)
  ENDIF

  CALL WATCH_START(TICTOC)

  CALL SCARF_IO(N, 'x', N(1)/2, V3, 'spec_struct_xslice')
  CALL SCARF_IO(N, 'y', N(2)/2, V3, 'spec_struct_yslice')
  CALL SCARF_IO(N, 'z', N(3)/2, V3, 'spec_struct_zslice')

  CALL WATCH_STOP(TICTOC)

  IF (RANK .EQ. 0) PRINT*, 'SLICE(S) WRITTEN IN              : ', REAL(TICTOC, REAL32)

  CALL WATCH_START(TICTOC)

  CALL SCARF_IO(N, V3, 'spec_struct_whole', 3)

  CALL WATCH_STOP(TICTOC)

  IF (RANK .EQ. 0) PRINT*, 'WHOLE FILE WRITTEN IN            : ', REAL(TICTOC, REAL32)

  CALL SCARF_FINALIZE()

  ! UNSTRUCTURED MESH TEST
  CALL SCARF_INITIALIZE(X1, Y1, Z1, DH, ACF, CL, SIGMA, METHOD = 1, HURST = HURST, DS = DS)
  !PARAMS = SCARF_INITIALIZE(X1, Y1, Z1, DH, ACF, CL, SIGMA, METHOD = 1, HURST = HURST)
  ! PARAMS = SCARF_INITIALIZE(X1, Y1, Z1, DH, ACF, CL, SIGMA, METHOD = 1, HURST, POI = POI, MUTE = MUTE, TAPER = TAPER)
  ! PARAMS = SCARF_INITIALIZE(X1, Y1, Z1, DH, ACF, CL, SIGMA, METHOD = 1, HURST)

  CALL WATCH_START(TICTOC)

  CALL SCARF_EXECUTE(SEED, V1, STATS)

  CALL WATCH_STOP(TICTOC)

  CALL MPI_ALLREDUCE(MPI_IN_PLACE, STATS, 8, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, IERR)

  IF (RANK .EQ. 0) THEN
    PRINT*, ''
    PRINT*, 'UNSTRUCTURED MESH TEST COMPLETED IN: ', REAL(TICTOC, REAL32)
    PRINT*, 'DOMAIN TOO SMALL?                  : ', NINT(STATS(1))
    PRINT*, 'GRID-STEP TOO LARGE?               : ', NINT(STATS(2))
    PRINT*, 'STANDARD DEVIATION                 : ', REAL(STATS(3), KIND=REAL32)
    PRINT*, 'MEAN VALUE                         : ', REAL(STATS(4), KIND=REAL32)
    PRINT*, 'TIMING FOR CPU EXECUTION           : ', REAL(STATS(5), KIND=REAL32)
    PRINT*, 'TIMING FOR GPU EXECUTION           : ', REAL(STATS(6), KIND=REAL32)
  ENDIF

  CALL WATCH_START(TICTOC)

  ! PARAMS%FS = FS
  ! PARAMS%FE = FE

  ! CALL SCARF_IO(N, 1, 400, V3, 'spec_unstruct_xslice')
  ! CALL SCARF_IO(N, 2, 250, V3, 'spec_unstruct_yslice')
  ! CALL SCARF_IO(N, 3, 100, V3, 'spec_unstruct_zslice')
  !
  ! CALL WATCH_STOP(TICTOC)
  !
  ! IF (RANK .EQ. 0) PRINT*, 'SLICE(S) WRITTEN IN                : ', REAL(TICTOC, REAL32)
  !
  ! CALL WATCH_START(TICTOC)
  !
  ! CALL SCARF_IO(N, V3, 'spec_unstruct_whole', 3)
  !
  ! CALL WATCH_STOP(TICTOC)
  !
  ! IF (RANK .EQ. 0) PRINT*, 'WHOLE FILE WRITTEN IN              : ', REAL(TICTOC, REAL32)

  CALL SCARF_FINALIZE()

#endif

  ! ================================================================================================================================-
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! ALL TESTS DONE, RELEASE MEMORY AND EXIT
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! ================================================================================================================================-

  NULLIFY(V3, X3, Y3, Z3)
  DEALLOCATE(V1, X1, Y1, Z1)

  CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

  IF (RANK .EQ. 0) THEN
    PRINT*, ''
    PRINT*, 'ALL TESTS COMPLETED'
  ENDIF

  CALL MPI_FINALIZE(IERR)

END PROGRAM DRIVER
