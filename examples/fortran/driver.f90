PROGRAM driver

  USE, INTRINSIC     :: iso_fortran_env, stdout => output_unit
  USE, NON_INTRINSIC :: mpi
  USE, NON_INTRINSIC :: m_scarflib
  USE, NON_INTRINSIC :: m_scarflib_aux

  IMPLICIT none

#ifdef DOUBLE_PREC
  INTEGER, PARAMETER :: f_int = int32
  INTEGER, PARAMETER :: f_real = real64
#else
  INTEGER, PARAMETER :: f_int = int32
  INTEGER, PARAMETER :: f_real = real32
#endif
  INTEGER, PARAMETER :: f_sgle = real32
  INTEGER, PARAMETER :: f_dble = real64

  CHARACTER(40)                                          :: buffer
  INTEGER(f_int)                                         :: i, j, k
  INTEGER(f_int)                                         :: rank, ntasks, ierr
  INTEGER(f_int)                                         :: acf, rescale, pad, seed
  INTEGER(f_int),              DIMENSION(3)              :: n, fs, fe, samples
  REAL(f_real)                                           :: ds, dh, sigma, hurst, mute, taper
  REAL(f_real),                DIMENSION(3)              :: cl
  REAL(f_real),                DIMENSION(8)              :: stats
  REAL(f_dble)                                           :: tictoc
  REAL(f_real),                DIMENSION(3,2)            :: poi
  REAL(f_real),                DIMENSION(:,:,:), POINTER :: x3, y3, z3, v3
  REAL(f_real),                DIMENSION(:,:),   POINTER :: x2, y2, v2
  REAL(f_real),   ALLOCATABLE, DIMENSION(:),     TARGET  :: x1, y1, z1, v1

  !--------------------------------------------------------------------------------------------------------------------------------

  ! initialise mpi
  CALL mpi_init(ierr)

  ! get rank number
  CALL mpi_comm_rank(mpi_comm_world, rank, ierr)

  ! get number of tasks
  CALL mpi_comm_size(mpi_comm_world, ntasks, ierr)

  CALL mpi_barrier(mpi_comm_world, ierr)

  !---------------------------------------------------------------------------------------------------------------------------------
  ! input section
  !---------------------------------------------------------------------------------------------------------------------------------

  ! number of points for whole model
  n = [500, 450, 400]

  ! grid step
  ds = 50._f_real

  dh = 50._f_real

  ! autocorrelation (0=von karman/exponential, 1=gaussian)
  acf = 0

  ! correlation length
  cl = [2000._f_real, 500._f_real, 100._f_real]

  ! standard deviation (sigma%/100)
  sigma = 0.05_f_real

  ! hurst exponent (not used for gaussian acf)
  hurst = 0.25_f_real

  ! seed number
  seed = 1235

  ! set position of point-of-interest (muting/tapering), same units as "dh"
  poi(:, 1) = [400., 250., 100.] * ds
  poi(:, 2) = [200., 150.,  50.] * ds

  ! radius for muting (at poi), same units as "dh"
  mute = 1000.

  ! radius for tapering (at poi + mute), same units as "dh"
  taper = 5000. * 10

  ! rescale to desired (continuous) sigma
  rescale = 0

  ! expand grid to handle fft periodicity
  pad = 0

  ! end input section
  !---------------------------------------------------------------------------------------------------------------------------------

  ! ================================================================================================================================
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! create sample structured mesh: 2D case
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! ================================================================================================================================

  CALL sample_mesh(rank, ntasks, n(1:2), fs(1:2), fe(1:2))

  samples = [fe(1) - fs(1) + 1, fe(2) - fs(2) + 1, 1]

  ALLOCATE(v1(PRODUCT(samples)), x1(PRODUCT(samples)), y1(PRODUCT(samples)))

  x2(fs(1):fe(1), fs(2):fe(2)) => x1
  y2(fs(1):fe(1), fs(2):fe(2)) => y1
  v2(fs(1):fe(1), fs(2):fe(2)) => v1(1:PRODUCT(samples))

  DO j = fs(2), fe(2)
    DO i = fs(1), fe(1)
      x2(i, j) = (i - 1) * ds
      y2(i, j) = (j - 1) * ds
    ENDDO
  ENDDO

  ! ================================================================================================================================
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! test FIM algorithm
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! ================================================================================================================================

  IF (rank .eq. 0) THEN
    WRITE(stdout, *) ''
    WRITE(stdout, *) '*************************************************'
    WRITE(stdout, *) '****** FIM algorithm, 2D, structured mesh *******'
  ENDIF

  CALL scarf_initialize(fs(1:2), fe(1:2), ds, acf, cl(1:2), sigma, method = 0, hurst = hurst, alpha = 30._f_real,    &
                        poi = poi(1:2,:), taper = taper, mute = mute)

  CALL watch_start(tictoc)

  CALL scarf_execute(seed, v2, stats)

  CALL watch_stop(tictoc)

  CALL mpi_allreduce(mpi_in_place, stats, 8, mpi_real, mpi_max, mpi_comm_world, ierr)

  IF (rank .eq. 0) THEN
    WRITE(stdout, *)
    buffer = 'Statistics for current simulation'
    WRITE(stdout, '(X, A35)') ADJUSTL(buffer)
    WRITE(stdout, *) '*************************************************'
    buffer = 'Elapsed time'
    WRITE(stdout, '(X, A30, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', tictoc, ' sec', '|'
    buffer = '+ spectrum'
    WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(5), ' sec', '|'
    buffer = '+ symmetry'
    WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(6), ' sec', '|'
    buffer = '+ FFT'
    WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(7), ' sec', '|'
    buffer = '+ interpolation'
    WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(8), ' sec', '|'
    WRITE(stdout, *) '------------------------------|-----------------|'
    buffer = 'Domain too small?'
    WRITE(stdout, '(X, A30, A, L12, T50, A)') ADJUSTL(buffer), '|', NINT(stats(1)), '|'
    buffer = 'Grid-step too large?'
    WRITE(stdout, '(X, A30, A, L12, T50, A)') ADJUSTL(buffer), '|', NINT(stats(2)), '|'
    WRITE(stdout, *) '------------------------------|-----------------|'
    buffer = 'Discrete standard deviation'
    WRITE(stdout, '(X, A30, A, F12.5, T50, A)') ADJUSTL(buffer), '|', stats(3), '|'
    buffer = 'Discrete mean value'
    WRITE(stdout, '(X, A30, A, F12.5, T50, A)') ADJUSTL(buffer), '|', stats(4), '|'
    WRITE(stdout, *) '------------------------------|-----------------|'
  ENDIF

  CALL watch_start(tictoc)

  CALL scarf_io(n(1:2), v2, 'fim_struct_whole_2d', 3)

  CALL watch_stop(tictoc)

  IF (rank .eq. 0) THEN
    buffer = 'I/O time'
    WRITE(stdout, '(X, A30, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', tictoc, ' sec', '|'
    WRITE(stdout, *) '*************************************************'
    WRITE(stdout, *)
  ENDIF

  CALL scarf_finalize()

  IF (rank .eq. 0) THEN
    WRITE(stdout, *) ''
    WRITE(stdout, *) '*************************************************'
    WRITE(stdout, *) '***** FIM algorithm, 2D, unstructured mesh ******'
  ENDIF

  CALL scarf_initialize(dh, acf, cl(1:2), sigma, x1, y1, method = 0,  &
                        hurst = hurst, alpha = 30._f_real, poi = poi(1:2,:), taper = taper, mute = mute)

  CALL watch_start(tictoc)

  CALL scarf_execute(seed, v1, stats)

  CALL watch_stop(tictoc)

  CALL mpi_allreduce(mpi_in_place, stats, 8, mpi_real, mpi_max, mpi_comm_world, ierr)

  IF (rank .eq. 0) THEN
    WRITE(stdout, *)
    buffer = 'Statistics for current simulation'
    WRITE(stdout, '(X, A35)') ADJUSTL(buffer)
    WRITE(stdout, *) '*************************************************'
    buffer = 'Elapsed time'
    WRITE(stdout, '(X, A30, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', tictoc, ' sec', '|'
    buffer = '+ spectrum'
    WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(5), ' sec', '|'
    buffer = '+ symmetry'
    WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(6), ' sec', '|'
    buffer = '+ FFT'
    WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(7), ' sec', '|'
    buffer = '+ interpolation'
    WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(8), ' sec', '|'
    WRITE(stdout, *) '------------------------------|-----------------|'
    buffer = 'Domain too small?'
    WRITE(stdout, '(X, A30, A, L12, T50, A)') ADJUSTL(buffer), '|', NINT(stats(1)), '|'
    buffer = 'Grid-step too large?'
    WRITE(stdout, '(X, A30, A, L12, T50, A)') ADJUSTL(buffer), '|', NINT(stats(2)), '|'
    WRITE(stdout, *) '------------------------------|-----------------|'
    buffer = 'Discrete standard deviation'
    WRITE(stdout, '(X, A30, A, F12.5, T50, A)') ADJUSTL(buffer), '|', stats(3), '|'
    buffer = 'Discrete mean value'
    WRITE(stdout, '(X, A30, A, F12.5, T50, A)') ADJUSTL(buffer), '|', stats(4), '|'
    WRITE(stdout, *) '------------------------------|-----------------|'
  ENDIF

  CALL watch_start(tictoc)

  CALL scarf_io(n(1:2), v2, 'fim_unstruct_whole_2d', 3)

  CALL watch_stop(tictoc)

  IF (rank .eq. 0) THEN
    buffer = 'I/O time'
    WRITE(stdout, '(X, A30, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', tictoc, ' sec', '|'
    WRITE(stdout, *) '*************************************************'
    WRITE(stdout, *)
  ENDIF

  CALL scarf_finalize()

  ! ================================================================================================================================
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! tests SRM algorithm
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! ================================================================================================================================

#ifdef SPECTRAL

  IF (rank .eq. 0) THEN
    WRITE(stdout, *) ''
    WRITE(stdout, *) '*************************************************'
    WRITE(stdout, *) '****** SRM algorithm, 2D, structured mesh *******'
  ENDIF

  ! structured mesh test
  CALL scarf_initialize(fs(1:2), fe(1:2), ds, acf, cl(1:2), sigma, method = 1, hurst = hurst, alpha = -20._f_real,    &
                        poi = poi, taper = taper, mute = mute)

  CALL watch_start(tictoc)

  CALL scarf_execute(seed, v2, stats)

  CALL watch_stop(tictoc)

  CALL mpi_allreduce(mpi_in_place, stats, 8, mpi_real, mpi_max, mpi_comm_world, ierr)

  IF (rank .eq. 0) THEN
    WRITE(stdout, *)
    buffer = 'Statistics for current simulation'
    WRITE(stdout, '(X, A35)') ADJUSTL(buffer)
    WRITE(stdout, *) '*************************************************'
    buffer = 'Elapsed time'
    WRITE(stdout, '(X, A30, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', tictoc, ' sec', '|'
    buffer = '+ CPU (main loop)'
    WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(5), ' sec', '|'
    buffer = '+ GPU (main loop)'
    WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(6), ' sec', '|'
    WRITE(stdout, *) '------------------------------|-----------------|'
    buffer = 'Domain too small?'
    WRITE(stdout, '(X, A30, A, L12, T50, A)') ADJUSTL(buffer), '|', NINT(stats(1)), '|'
    buffer = 'Grid-step too large?'
    WRITE(stdout, '(X, A30, A, L12, T50, A)') ADJUSTL(buffer), '|', NINT(stats(2)), '|'
    WRITE(stdout, *) '------------------------------|-----------------|'
    buffer = 'Discrete standard deviation'
    WRITE(stdout, '(X, A30, A, F12.5, T50, A)') ADJUSTL(buffer), '|', stats(3), '|'
    buffer = 'Discrete mean value'
    WRITE(stdout, '(X, A30, A, F12.5, T50, A)') ADJUSTL(buffer), '|', stats(4), '|'
    WRITE(stdout, *) '------------------------------|-----------------|'
  ENDIF

  CALL watch_start(tictoc)

  CALL scarf_io(n(1:2), v2, 'srm_struct_whole_2d', 3)

  CALL watch_stop(tictoc)

  IF (rank .eq. 0) THEN
    buffer = 'I/O time'
    WRITE(stdout, '(X, A30, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', tictoc, ' sec', '|'
    WRITE(stdout, *) '*************************************************'
    WRITE(stdout, *)
  ENDIF

  CALL scarf_finalize()


  IF (rank .eq. 0) THEN
    WRITE(stdout, *)
    WRITE(stdout, *) '*************************************************'
    WRITE(stdout, *) '***** SRM algorithm, 2D, unstructured mesh ******'
  ENDIF

  ! unstructured mesh test
  CALL scarf_initialize(dh, acf, cl(1:2), sigma, x1, y1, method = 1,  &
                        hurst = hurst, alpha = -20._f_real, poi = poi(1:2,:), taper = taper, mute = mute)

  CALL watch_start(tictoc)

  CALL scarf_execute(seed, v1, stats)

  CALL watch_stop(tictoc)

  CALL mpi_allreduce(mpi_in_place, stats, 8, mpi_real, mpi_max, mpi_comm_world, ierr)

  IF (rank .eq. 0) THEN
    WRITE(stdout, *)
    buffer = 'Statistics for current simulation'
    WRITE(stdout, '(X, A35)') ADJUSTL(buffer)
    WRITE(stdout, *) '*************************************************'
    buffer = 'Elapsed time'
    WRITE(stdout, '(X, A30, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', tictoc, ' sec', '|'
    buffer = '+ CPU (main loop)'
    WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(5), ' sec', '|'
    buffer = '+ GPU (main loop)'
    WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(6), ' sec', '|'
    WRITE(stdout, *) '------------------------------|-----------------|'
    buffer = 'Domain too small?'
    WRITE(stdout, '(X, A30, A, L12, T50, A)') ADJUSTL(buffer), '|', NINT(stats(1)), '|'
    buffer = 'Grid-step too large?'
    WRITE(stdout, '(X, A30, A, L12, T50, A)') ADJUSTL(buffer), '|', NINT(stats(2)), '|'
    WRITE(stdout, *) '------------------------------|-----------------|'
    buffer = 'Discrete standard deviation'
    WRITE(stdout, '(X, A30, A, F12.5, T50, A)') ADJUSTL(buffer), '|', stats(3), '|'
    buffer = 'Discrete mean value'
    WRITE(stdout, '(X, A30, A, F12.5, T50, A)') ADJUSTL(buffer), '|', stats(4), '|'
    WRITE(stdout, *) '------------------------------|-----------------|'
  ENDIF

  CALL watch_start(tictoc)

  CALL scarf_io(n(1:2), v2, 'srm_unstruct_whole_2d', 3)

  CALL watch_stop(tictoc)

  IF (rank .eq. 0) THEN
    buffer = 'I/O time'
    WRITE(stdout, '(X, A30, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', tictoc, ' sec', '|'
    WRITE(stdout, *) '*************************************************'
    WRITE(stdout, *)
  ENDIF

  CALL scarf_finalize()

#endif

  ! release resources
  NULLIFY(v2, x2, y2)
  DEALLOCATE(v1, x1, y1)

! ================================================================================================================================
! --------------------------------------------------------------------------------------------------------------------------------
! create sample structured mesh: 3D case
! --------------------------------------------------------------------------------------------------------------------------------
! ================================================================================================================================

CALL sample_mesh(rank, ntasks, n, fs, fe)

samples = [fe(1) - fs(1) + 1, fe(2) - fs(2) + 1, fe(3) - fs(3) + 1]

ALLOCATE(v1(PRODUCT(samples)), x1(PRODUCT(samples)), y1(PRODUCT(samples)), z1(PRODUCT(samples)))

x3(fs(1):fe(1), fs(2):fe(2), fs(3):fe(3)) => x1
y3(fs(1):fe(1), fs(2):fe(2), fs(3):fe(3)) => y1
z3(fs(1):fe(1), fs(2):fe(2), fs(3):fe(3)) => z1
v3(fs(1):fe(1), fs(2):fe(2), fs(3):fe(3)) => v1

DO k = fs(3), fe(3)
  DO j = fs(2), fe(2)
    DO i = fs(1), fe(1)
      x3(i, j, k) = (i - 1) * ds
      y3(i, j, k) = (j - 1) * ds
      z3(i, j, k) = (k - 1) * ds
    ENDDO
  ENDDO
ENDDO

! ================================================================================================================================
! --------------------------------------------------------------------------------------------------------------------------------
! test FIM algorithm
! --------------------------------------------------------------------------------------------------------------------------------
! ================================================================================================================================

IF (rank .eq. 0) THEN
  WRITE(stdout, *) ''
  WRITE(stdout, *) '*************************************************'
  WRITE(stdout, *) '****** FIM algorithm, 3D, structured mesh *******'
ENDIF

CALL scarf_initialize(fs, fe, ds, acf, cl, sigma, method = 0, hurst = hurst, beta = 10._f_real, poi = poi, taper = taper,  &
                      mute = mute)

CALL watch_start(tictoc)

CALL scarf_execute(seed, v3, stats)

CALL watch_stop(tictoc)

CALL mpi_allreduce(mpi_in_place, stats, 8, mpi_real, mpi_max, mpi_comm_world, ierr)

IF (rank .eq. 0) THEN
  WRITE(stdout, *)
  buffer = 'Statistics for current simulation'
  WRITE(stdout, '(X, A35)') ADJUSTL(buffer)
  WRITE(stdout, *) '*************************************************'
  buffer = 'Elapsed time'
  WRITE(stdout, '(X, A30, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', tictoc, ' sec', '|'
  buffer = '+ spectrum'
  WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(5), ' sec', '|'
  buffer = '+ symmetry'
  WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(6), ' sec', '|'
  buffer = '+ FFT'
  WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(7), ' sec', '|'
  buffer = '+ interpolation'
  WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(8), ' sec', '|'
  WRITE(stdout, *) '------------------------------|-----------------|'
  buffer = 'Domain too small?'
  WRITE(stdout, '(X, A30, A, L12, T50, A)') ADJUSTL(buffer), '|', NINT(stats(1)), '|'
  buffer = 'Grid-step too large?'
  WRITE(stdout, '(X, A30, A, L12, T50, A)') ADJUSTL(buffer), '|', NINT(stats(2)), '|'
  WRITE(stdout, *) '------------------------------|-----------------|'
  buffer = 'Discrete standard deviation'
  WRITE(stdout, '(X, A30, A, F12.5, T50, A)') ADJUSTL(buffer), '|', stats(3), '|'
  buffer = 'Discrete mean value'
  WRITE(stdout, '(X, A30, A, F12.5, T50, A)') ADJUSTL(buffer), '|', stats(4), '|'
  WRITE(stdout, *) '------------------------------|-----------------|'
ENDIF

CALL watch_start(tictoc)

CALL scarf_io(n, v3, 'fim_struct_whole_3d', 3)

CALL watch_stop(tictoc)

IF (rank .eq. 0) THEN
  buffer = 'I/O time'
  WRITE(stdout, '(X, A30, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', tictoc, ' sec', '|'
  WRITE(stdout, *) '*************************************************'
  WRITE(stdout, *)
ENDIF

CALL scarf_finalize()

IF (rank .eq. 0) THEN
  WRITE(stdout, *) ''
  WRITE(stdout, *) '*************************************************'
  WRITE(stdout, *) '***** FIM algorithm, 3D, unstructured mesh ******'
ENDIF

CALL scarf_initialize(dh, acf, cl, sigma, x1, y1, z1, method = 0, hurst = hurst, beta = 10._f_real, poi = poi, taper = taper,  &
                      mute = mute)

CALL watch_start(tictoc)

CALL scarf_execute(seed, v1, stats)

CALL watch_stop(tictoc)

CALL mpi_allreduce(mpi_in_place, stats, 8, mpi_real, mpi_max, mpi_comm_world, ierr)

IF (rank .eq. 0) THEN
  WRITE(stdout, *)
  buffer = 'Statistics for current simulation'
  WRITE(stdout, '(X, A35)') ADJUSTL(buffer)
  WRITE(stdout, *) '*************************************************'
  buffer = 'Elapsed time'
  WRITE(stdout, '(X, A30, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', tictoc, ' sec', '|'
  buffer = '+ spectrum'
  WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(5), ' sec', '|'
  buffer = '+ symmetry'
  WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(6), ' sec', '|'
  buffer = '+ FFT'
  WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(7), ' sec', '|'
  buffer = '+ interpolation'
  WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(8), ' sec', '|'
  WRITE(stdout, *) '------------------------------|-----------------|'
  buffer = 'Domain too small?'
  WRITE(stdout, '(X, A30, A, L12, T50, A)') ADJUSTL(buffer), '|', NINT(stats(1)), '|'
  buffer = 'Grid-step too large?'
  WRITE(stdout, '(X, A30, A, L12, T50, A)') ADJUSTL(buffer), '|', NINT(stats(2)), '|'
  WRITE(stdout, *) '------------------------------|-----------------|'
  buffer = 'Discrete standard deviation'
  WRITE(stdout, '(X, A30, A, F12.5, T50, A)') ADJUSTL(buffer), '|', stats(3), '|'
  buffer = 'Discrete mean value'
  WRITE(stdout, '(X, A30, A, F12.5, T50, A)') ADJUSTL(buffer), '|', stats(4), '|'
  WRITE(stdout, *) '------------------------------|-----------------|'
ENDIF

CALL watch_start(tictoc)

CALL scarf_io(n, v3, 'fim_unstruct_whole_3d', 3)

CALL watch_stop(tictoc)

IF (rank .eq. 0) THEN
  buffer = 'I/O time'
  WRITE(stdout, '(X, A30, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', tictoc, ' sec', '|'
  WRITE(stdout, *) '*************************************************'
  WRITE(stdout, *)
ENDIF

CALL scarf_finalize()

! ================================================================================================================================-
! --------------------------------------------------------------------------------------------------------------------------------
! tests SRM algorithm
! --------------------------------------------------------------------------------------------------------------------------------
! ================================================================================================================================-

#ifdef SPECTRAL

IF (rank .eq. 0) THEN
  WRITE(stdout, *) ''
  WRITE(stdout, *) '*************************************************'
  WRITE(stdout, *) '****** SRM algorithm, 3D, structured mesh *******'
ENDIF

CALL scarf_initialize(fs, fe, ds, acf, cl, sigma, method = 1, hurst = hurst, beta = -20._f_real, poi=poi, taper=taper, mute=mute)

CALL watch_start(tictoc)

CALL scarf_execute(seed, v3, stats)

CALL watch_stop(tictoc)

CALL mpi_allreduce(mpi_in_place, stats, 8, mpi_real, mpi_max, mpi_comm_world, ierr)

IF (rank .eq. 0) THEN
  WRITE(stdout, *)
  buffer = 'Statistics for current simulation'
  WRITE(stdout, '(X, A35)') ADJUSTL(buffer)
  WRITE(stdout, *) '*************************************************'
  buffer = 'Elapsed time'
  WRITE(stdout, '(X, A30, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', tictoc, ' sec', '|'
  buffer = '+ CPU (main loop)'
  WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(5), ' sec', '|'
  buffer = '+ GPU (main loop)'
  WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(6), ' sec', '|'
  WRITE(stdout, *) '------------------------------|-----------------|'
  buffer = 'Domain too small?'
  WRITE(stdout, '(X, A30, A, L12, T50, A)') ADJUSTL(buffer), '|', NINT(stats(1)), '|'
  buffer = 'Grid-step too large?'
  WRITE(stdout, '(X, A30, A, L12, T50, A)') ADJUSTL(buffer), '|', NINT(stats(2)), '|'
  WRITE(stdout, *) '------------------------------|-----------------|'
  buffer = 'Discrete standard deviation'
  WRITE(stdout, '(X, A30, A, F12.5, T50, A)') ADJUSTL(buffer), '|', stats(3), '|'
  buffer = 'Discrete mean value'
  WRITE(stdout, '(X, A30, A, F12.5, T50, A)') ADJUSTL(buffer), '|', stats(4), '|'
  WRITE(stdout, *) '------------------------------|-----------------|'
ENDIF

CALL watch_start(tictoc)

CALL scarf_io(n, v3, 'srm_struct_whole_3d', 3)

CALL watch_stop(tictoc)

IF (rank .eq. 0) THEN
  buffer = 'I/O time'
  WRITE(stdout, '(X, A30, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', tictoc, ' sec', '|'
  WRITE(stdout, *) '*************************************************'
  WRITE(stdout, *)
ENDIF

CALL scarf_finalize()


IF (rank .eq. 0) THEN
  WRITE(stdout, *)
  WRITE(stdout, *) '*************************************************'
  WRITE(stdout, *) '***** SRM algorithm, 3D, unstructured mesh ******'
ENDIF

CALL scarf_initialize(dh, acf, cl, sigma, x1, y1, z1, method = 1, hurst = hurst, beta = -20._f_real, poi = poi, taper = taper,  &
                      mute = mute)

CALL watch_start(tictoc)

CALL scarf_execute(seed, v1, stats)

CALL watch_stop(tictoc)

CALL mpi_allreduce(mpi_in_place, stats, 8, mpi_real, mpi_max, mpi_comm_world, ierr)

IF (rank .eq. 0) THEN
  WRITE(stdout, *)
  buffer = 'Statistics for current simulation'
  WRITE(stdout, '(X, A35)') ADJUSTL(buffer)
  WRITE(stdout, *) '*************************************************'
  buffer = 'Elapsed time'
  WRITE(stdout, '(X, A30, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', tictoc, ' sec', '|'
  buffer = '+ CPU (main loop)'
  WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(5), ' sec', '|'
  buffer = '+ GPU (main loop)'
  WRITE(stdout, '(4X, A27, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', stats(6), ' sec', '|'
  WRITE(stdout, *) '------------------------------|-----------------|'
  buffer = 'Domain too small?'
  WRITE(stdout, '(X, A30, A, L12, T50, A)') ADJUSTL(buffer), '|', NINT(stats(1)), '|'
  buffer = 'Grid-step too large?'
  WRITE(stdout, '(X, A30, A, L12, T50, A)') ADJUSTL(buffer), '|', NINT(stats(2)), '|'
  WRITE(stdout, *) '------------------------------|-----------------|'
  buffer = 'Discrete standard deviation'
  WRITE(stdout, '(X, A30, A, F12.5, T50, A)') ADJUSTL(buffer), '|', stats(3), '|'
  buffer = 'Discrete mean value'
  WRITE(stdout, '(X, A30, A, F12.5, T50, A)') ADJUSTL(buffer), '|', stats(4), '|'
  WRITE(stdout, *) '------------------------------|-----------------|'
ENDIF

CALL watch_start(tictoc)

CALL scarf_io(n, v3, 'srm_unstruct_whole_3d', 3)

CALL watch_stop(tictoc)

IF (rank .eq. 0) THEN
  buffer = 'I/O time'
  WRITE(stdout, '(X, A30, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', tictoc, ' sec', '|'
  WRITE(stdout, *) '*************************************************'
  WRITE(stdout, *)
ENDIF

CALL scarf_finalize()

#endif

  ! ================================================================================================================================-
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! all tests done, release resources and exit
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! ================================================================================================================================-

  NULLIFY(v3, x3, y3, z3, v2, x2, y2)
  DEALLOCATE(v1, x1, y1, z1)

  CALL mpi_barrier(mpi_comm_world, ierr)

  IF (rank .eq. 0) THEN
    WRITE(stdout, *) ''
    WRITE(stdout, *) 'All tests completed'
  ENDIF

  CALL mpi_finalize(ierr)

END PROGRAM driver
