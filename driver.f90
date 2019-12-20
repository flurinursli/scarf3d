PROGRAM DRIVER

  ! driver program for scarf3d
  ! compile: mpif90 -O3 scarflib_fft2.f90 scarflib_spec.f90 scarflib.f90 driver.f90 *.f -I/home/walter/Backedup/Software/fftw-3.3.8/include
  !                 -L/home/walter/Backedup/Software/fftw-3.3.8/lib -lfftw3 -fcheck=all -fbacktrace -fopenacc
  !
  ! the product of "dims" below should be equal to the number of processes


  USE :: MPI

  USE, NON_INTRINSIC :: SCARFLIB
  !USE, NON_INTRINSIC :: SCARFLIB_FFT
  USE, NON_INTRINSIC :: SCARFLIB_SPECTRAL
  USE, NON_INTRINSIC :: SCARFLIB_FFT2

  IMPLICIT NONE

  CHARACTER(:), ALLOCATABLE                   :: ACF
  INTEGER(IPP)                                :: I, J, K
  INTEGER(IPP)                                :: IERR, RANK, NTASKS, TOPO, NDIMS
  INTEGER(IPP)                                :: SEED, MUTE, TAPERING
  INTEGER(IPP),              DIMENSION(3)     :: NPTS, FS, FE, COORDS, DIMS
  INTEGER(IPP),              DIMENSION(3,2)   :: POI
  LOGICAL                                     :: REORDER
  LOGICAL,                   DIMENSION(3)     :: ISPERIODIC
  REAL(FPP)                                           :: DH, SIGMA, HURST, TICTOC
  REAL(FPP),                 DIMENSION(3)             :: CL
  REAL(FPP),                 DIMENSION(3)             :: TOC
  REAL(FPP),                 DIMENSION(:,:,:),   POINTER :: XP, YP, ZP, VP
  REAL(FPP),    ALLOCATABLE, DIMENSION(:), TARGET  :: X, Y, Z
  REAL(FPP),    ALLOCATABLE, DIMENSION(:), TARGET  :: V

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
  NPTS = [200, 200, 120]

  ! GRID STEP
  DH = 100._FPP

  ! AUTOCORRELATION
  ACF = 'GAUSS'

  ! CORRELATION LENGTH
  CL = [1000._FPP, 1000._FPP, 1000._FPP]

  ! STANDARD DEVIATION (SIGMA%/100)
  SIGMA = 0.05_FPP

  ! HURST EXPONENT (NOT USED FOR GAUSSIAN ACF)
  HURST = 0.5_FPP

  ! SEED NUMBER
  SEED = 1235

  ! SET POSITION OF POINT-OF-INTEREST (MUTING/TAPERING)
  POI(:, 1) = [300, 250, 150]
  POI(:, 2) = [400, 250, 150]

  ! RADIUS FOR MUTING (AT POI)
  MUTE = 0

  ! RADIUS FOR TAPERING (AT POI + MUTE)
  TAPERING = 0

  ! END INPUT SECTION
  !-----------------------------------------------------------------------

  ! HERE BELOW WE FIRST CREATE A SAMPLE CARTESIAN TOPOLOGY THAT COULD REPRESENT HOW A 3D DOMAIN IS DECOMPOSED, E.G., IN FINITE-DIFFERENCE
  ! CODES. WE THEN COMPUTE THE RANDOM PERTURBATIONS AND ADD THEM TO AN HYPOTETHICAL VELOCITY FIELD. THE PERTURBED FIELD IS THEN WRITTEN
  ! TO DISK.

  ! WE WORK IN 3D
  NDIMS = 3

  ! NUMBER OF PROCESSES ALONG EACH DIRECTION
  DIMS  = [2, 1, 2]

  ! ALLOW REORDERING
  REORDER = .TRUE.

  ! NO PERIODICITY
  ISPERIODIC = [.FALSE., .FALSE., .FALSE.]

  ! CREATE TOPOLOGY
  CALL MPI_CART_CREATE(MPI_COMM_WORLD, NDIMS, DIMS, ISPERIODIC, REORDER, TOPO, IERR)

  CALL MPI_CART_COORDS(TOPO, RANK, NDIMS, COORDS, IERR)

  ! SPLIT EACH AXIS AND ASSIGN POINTS TO CURRENT PROCESS
  CALL MPI_SPLIT_TASK(NPTS, DIMS, COORDS, FS, FE)


  ! PREPARE ARRAY WITH SAMPLE VELOCITY FIELD
  ALLOCATE(V((FE(1) - FS(1) + 1) * (FE(2) - FS(2) + 1) * (FE(3) - FS(3) + 1)))

  VP(FS(1):FE(1), FS(2):FE(2), FS(3):FE(3)) => V

  VP = 1._FPP

  ! PREPARE ARRAY WITH POINTS POSITION
  ALLOCATE(X((FE(1) - FS(1) + 1) * (FE(2) - FS(2) + 1) * (FE(3) - FS(3) + 1)))
  ALLOCATE(Y((FE(1) - FS(1) + 1) * (FE(2) - FS(2) + 1) * (FE(3) - FS(3) + 1)))
  ALLOCATE(Z((FE(1) - FS(1) + 1) * (FE(2) - FS(2) + 1) * (FE(3) - FS(3) + 1)))

  XP(FS(1):FE(1), FS(2):FE(2), FS(3):FE(3)) => X
  YP(FS(1):FE(1), FS(2):FE(2), FS(3):FE(3)) => Y
  ZP(FS(1):FE(1), FS(2):FE(2), FS(3):FE(3)) => Z

  DO K = FS(3), FE(3)
    DO J = FS(2), FE(2)
      DO I = FS(1), FE(1)
        XP(I, J, K) = (I - 0.5) * DH
        YP(I, J, K) = (J - 0.5) * DH
        ZP(I, J, K) = (K - 0.5) * DH
      ENDDO
    ENDDO
  ENDDO

  ! COMPUTE ZERO-MEAN RANDOM FIELD
  !CALL SCARF3D_FFT(MPI_COMM_WORLD, [NX, NY, NZ], 100._fpp, 'VK', [5000._fpp, 5000._fpp, 5000._fpp], 0.05_fpp, 0.1_fpp, 1235, POI, 4, 26)
  !CALL SCARF3D_FFT(MPI_COMM_WORLD, NPTS, 100._fpp, 'VK', [5000._fpp, 5000._fpp, 5000._fpp], 0.05_fpp, 0.5_fpp, 1235, POI, 0, 0)

  !=========== GPU ============

  IF (RANK .EQ. 0) PRINT*, 'GPU'

  CALL SCARF3D_SPEC(FS, FE, X, Y, Z, ACF, CL, SIGMA, HURST, SEED, POI, 0, 0, V, TICTOC)

  ! REMOVE AVERAGE VALUE FROM PERTURBED FIELD
  V = V - 1._FPP

  IF (RANK .EQ. 0) PRINT*, MINVAL(V), MAXVAL(V)

  ! WRITE A MODEL SLICE TO DISK
  CALL IO_WRITE_SLICE(2, 1, VP, 'gpu_slice.bin')

  ! SOME TIMING INFO
  IF (RANK .EQ. 0) PRINT*, 'TIMING GPU: ', REAL(TICTOC, KIND=REAL32)

  NULLIFY(XP, YP, ZP)
  DEALLOCATE(X, Y, Z)

  DEALLOCATE(GS, GE)

  V = 1._FPP

  ! =========== FFT ==============

  IF (RANK .EQ. 0) PRINT*, 'FFT'

  CALL SCARF3D_STRUCTURED(FS, FE, DH, ACF, CL, SIGMA, HURST, SEED`, POI, 0, 0, VP, TOC)

  ! REMOVE AVERAGE VALUE FROM PERTURBED FIELD
  V = V - 1._FPP

  IF (RANK .EQ. 0) PRINT*, MINVAL(V), MAXVAL(V)

  ! WRITE A MODEL SLICE TO DISK
  CALL IO_WRITE_SLICE(2, 1, VP, 'fft_slice.bin')

  ! SOME TIMING INFO
  IF (RANK .EQ. 0) PRINT*, 'TIMING: ', REAL(TOC, KIND=REAL32)

  NULLIFY(VP)
  DEALLOCATE(V)

  CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

  CALL MPI_FINALIZE(IERR)

END PROGRAM DRIVER

! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
!===============================================================================================================================
! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

SUBROUTINE MPI_SPLIT_TASK(NPTS, NTASKS, RANK, FS, FE)

  USE, NON_INTRINSIC :: SCARFLIB

  IMPLICIT NONE

  INTEGER(IPP), DIMENSION(3), INTENT(IN)  :: NPTS                          !< NUMBER OF POINTS TO BE SPLIT
  INTEGER(IPP), DIMENSION(3), INTENT(IN)  :: NTASKS                     !< NUMBER OF MPI PROCESSES
  INTEGER(IPP), DIMENSION(3), INTENT(IN)  :: RANK                       !< RANK OF CURRENT PREOCESS
  INTEGER(IPP), DIMENSION(3), INTENT(OUT) :: FS, FE                     !< 1ST/LAST INDICES
  INTEGER(IPP)                            :: I                          !< COUNTER

  !------------------------------------------------------------------------------------------------------------------------------

  DO I = 1, 3
    FS(I) = 1 + INT( REAL(NPTS(I), FPP) / REAL(NTASKS(I), FPP) * REAL(RANK(I), FPP) )
    FE(I) = INT( REAL(NPTS(I), FPP) / REAL(NTASKS(I), FPP) * REAL(RANK(I) + 1, FPP) )
  ENDDO

END SUBROUTINE MPI_SPLIT_TASK
