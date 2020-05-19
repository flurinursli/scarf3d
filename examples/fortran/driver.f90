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


  INTEGER(f_int)                                         :: i, j, k
  INTEGER(f_int)                                         :: rank, ntasks, ierr
  INTEGER(f_int)                                         :: acf, rescale, pad, seed
  INTEGER(f_int),              DIMENSION(3)              :: n, fs, fe
  REAL(f_real)                                           :: ds, dh, sigma, hurst, mute, taper
  REAL(f_real),                DIMENSION(3)              :: cl
  REAL(f_real),                DIMENSION(8)              :: stats
  REAL(f_dble)                                           :: tictoc
  REAL(f_real),                DIMENSION(3,2)            :: poi
  REAL(f_real),                DIMENSION(:,:,:), POINTER :: x3, y3, z3, v3
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
  !n = [2000*2, 3200*2, 1200*2]
  n = [500, 500, 500]

  ! grid step
  ds = 50._f_real

  dh = 50._f_real

  ! autocorrelation (0=von karman/exponential, 1=gaussian)
  acf = 0

  ! correlation length
  cl = [2000._f_real, 2000._f_real, 2000._f_real]

  ! standard deviation (sigma%/100)
  sigma = 0.05_f_real

  ! hurst exponent (not used for gaussian acf)
  hurst = 0.25_f_real

  ! seed number
  seed = 1235

  ! set position of point-of-interest (muting/tapering), same units as "dh"
  poi(:, 1) = [400., 250., 100.] * ds
  poi(:, 2) = [350., 200.,  50.] * ds

  ! radius for muting (at poi), same units as "dh"
  mute = 1000. !1000.

  ! radius for tapering (at poi + mute), same units as "dh"
  taper = 5000. !5000.

  ! rescale to desired (continuous) sigma
  rescale = 0

  ! expand grid to handle fft periodicity
  pad = 0

  ! end input section
  !---------------------------------------------------------------------------------------------------------------------------------

  ! ================================================================================================================================-
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! create sample structured mesh
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! ================================================================================================================================-

  CALL sample_mesh(rank, ntasks, n, fs, fe)

  ! allocate memory for mesh
  ALLOCATE(v1((fe(1) - fs(1) + 1) * (fe(2) - fs(2) + 1) * (fe(3) - fs(3) + 1)))
  ALLOCATE(x1((fe(1) - fs(1) + 1) * (fe(2) - fs(2) + 1) * (fe(3) - fs(3) + 1)))
  ALLOCATE(y1((fe(1) - fs(1) + 1) * (fe(2) - fs(2) + 1) * (fe(3) - fs(3) + 1)))
  ALLOCATE(z1((fe(1) - fs(1) + 1) * (fe(2) - fs(2) + 1) * (fe(3) - fs(3) + 1)))

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

  ! ================================================================================================================================-
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! fft method tests
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! ================================================================================================================================-

  IF (rank .eq. 0) WRITE(stdout, *) ''
  IF (rank .eq. 0) WRITE(stdout, *) '************************************************************'
  IF (rank .eq. 0) WRITE(stdout, *) '************************ FIM method ************************'

  ! structured mesh test
  CALL scarf_initialize(fs, fe, ds, acf, cl, sigma, method = 0, hurst = hurst)

  CALL watch_start(tictoc)

  CALL scarf_execute(seed, v3, stats)

  CALL watch_stop(tictoc)

  CALL mpi_allreduce(mpi_in_place, stats, 8, mpi_real, mpi_max, mpi_comm_world, ierr)

  IF (rank .eq. 0) THEN
    WRITE(stdout, *) ''
    WRITE(stdout, *) 'Summary structured mesh'
    WRITE(stdout, *) '  i)   test completed in       : ', REAL(tictoc, f_sgle),   ' sec'
    WRITE(stdout, *) '  ii)  domain too small?       : ', NINT(stats(1))
    WRITE(stdout, *) '  iii) grid-step too large?    : ', NINT(stats(2))
    WRITE(stdout, *) '  iv)  standard deviation      : ', REAL(stats(3), f_sgle)
    WRITE(stdout, *) '  vi)  mean value              : ', REAL(stats(4), f_sgle)
    WRITE(stdout, *) '  vii) timing for spectrum     : ', REAL(stats(5), f_sgle), ' sec'
    WRITE(stdout, *) '  viii)timing for symmetry     : ', REAL(stats(6), f_sgle), ' sec'
    WRITE(stdout, *) '  ix)  timing for ifft         : ', REAL(stats(7), f_sgle), ' sec'
    WRITE(stdout, *) '  x)   timing for interpolation: ', REAL(stats(8), f_sgle), ' sec'
  ENDIF

  CALL watch_start(tictoc)

  CALL scarf_io(n, 'x', n(1)/2, v3, 'fft_struct_xslice')
  CALL scarf_io(n, 'y', n(2)/2, v3, 'fft_struct_yslice')
  CALL scarf_io(n, 'z', n(3)/2, v3, 'fft_struct_zslice')

  CALL watch_stop(tictoc)
  IF (rank .eq. 0) WRITE(stdout, *) '  xi)  slice(s) written in     : ', REAL(tictoc, f_sgle), ' sec'

  CALL watch_start(tictoc)

  CALL scarf_io(n, v3, 'fft_struct_whole', 3)

  CALL watch_stop(tictoc)

  IF (rank .eq. 0) WRITE(stdout, *) '  xii) whole file written in   : ', REAL(tictoc, f_sgle), ' sec'

  CALL scarf_finalize()

  ! unstructured mesh test
  CALL scarf_initialize(x1, y1, z1, dh, acf, cl, sigma, method = 0, hurst = hurst)

  CALL watch_start(tictoc)

  CALL scarf_execute(seed, v1, stats)

  CALL watch_stop(tictoc)

  CALL mpi_allreduce(mpi_in_place, stats, 8, mpi_real, mpi_max, mpi_comm_world, ierr)

  IF (rank .eq. 0) THEN
    WRITE(stdout, *) ''
    WRITE(stdout, *) 'Summary unstructured mesh'
    WRITE(stdout, *) '  i)   test completed in       : ', REAL(tictoc, f_sgle),   ' sec'
    WRITE(stdout, *) '  ii)  domain too small?       : ', NINT(stats(1))
    WRITE(stdout, *) '  iii) grid-step too large?    : ', NINT(stats(2))
    WRITE(stdout, *) '  iv)  standard deviation      : ', REAL(stats(3), f_sgle)
    WRITE(stdout, *) '  v)   mean value              : ', REAL(stats(4), f_sgle)
    WRITE(stdout, *) '  vi)  timing for spectrum     : ', REAL(stats(5), f_sgle), ' sec'
    WRITE(stdout, *) '  vii) timing for symmetry     : ', REAL(stats(6), f_sgle), ' sec'
    WRITE(stdout, *) '  viii)timing for ifft         : ', REAL(stats(7), f_sgle), ' sec'
    WRITE(stdout, *) '  ix)  timing for interpolation: ', REAL(stats(8), f_sgle), ' sec'
  ENDIF

  CALL watch_start(tictoc)

  CALL scarf_finalize()

  ! ================================================================================================================================-
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! spectral method tests
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! ================================================================================================================================-

#ifdef spectral

  IF (rank .eq. 0) WRITE(stdout, *) ''
  IF (rank .eq. 0) WRITE(stdout, *) '************************************************************'
  IF (rank .eq. 0) WRITE(stdout, *) '************************ SRM method ************************'

  ! structured mesh test
  CALL scarf_initialize(fs, fe, ds, acf, cl, sigma, method = 1, hurst = hurst)

  CALL watch_start(tictoc)

  CALL scarf_execute(seed, v3, stats)

  CALL watch_stop(tictoc)

  CALL mpi_allreduce(mpi_in_place, stats, 8, mpi_real, mpi_max, mpi_comm_world, ierr)

  IF (rank .eq. 0) THEN
    WRITE(stdout, *) ''
    WRITE(stdout, *) 'Summary structured mesh'
    WRITE(stdout, *) '  i)   test completed in    : ', REAL(tictoc, f_sgle),   ' sec'
    WRITE(stdout, *) '  ii)  domain too small?    : ', NINT(stats(1))
    WRITE(stdout, *) '  iii) grid-step too large? : ', NINT(stats(2))
    WRITE(stdout, *) '  iv)  standard deviation   : ', REAL(stats(3), f_sgle)
    WRITE(stdout, *) '  vi)  mean value           : ', REAL(stats(4), f_sgle)
    WRITE(stdout, *) '  vii) CPU main loop        : ', REAL(stats(5), f_sgle), ' sec'
    WRITE(stdout, *) '  viii)GPU main loop        : ', REAL(stats(6), f_sgle), ' sec'
  ENDIF

  CALL watch_start(tictoc)

  CALL scarf_io(n, 'x', n(1)/2, v3, 'spec_struct_xslice')
  CALL scarf_io(n, 'y', n(2)/2, v3, 'spec_struct_yslice')
  CALL scarf_io(n, 'z', n(3)/2, v3, 'spec_struct_zslice')

  CALL watch_stop(tictoc)

  IF (rank .eq. 0) WRITE(stdout, *) '  ix)  slice(s) written in  : ', REAL(tictoc, f_sgle), ' sec'

  CALL watch_start(tictoc)

  CALL scarf_io(n, v3, 'spec_struct_whole', 3)

  CALL watch_stop(tictoc)

  IF (rank .eq. 0) WRITE(stdout, *) '  x)   whole file written in: ', REAL(tictoc, f_sgle), ' sec'

  CALL scarf_finalize()

  ! unstructured mesh test
  CALL scarf_initialize(x1, y1, z1, dh, acf, cl, sigma, method = 1, hurst = hurst)

  CALL watch_start(tictoc)

  CALL scarf_execute(seed, v1, stats)

  CALL watch_stop(tictoc)

  CALL mpi_allreduce(mpi_in_place, stats, 8, mpi_real, mpi_max, mpi_comm_world, ierr)

  IF (rank .eq. 0) THEN
    WRITE(stdout, *) ''
    WRITE(stdout, *) 'Summary unstructured mesh'
    WRITE(stdout, *) '  i)   test completed in   : ', REAL(tictoc, f_sgle),   ' sec'
    WRITE(stdout, *) '  ii)  domain too small?   : ', NINT(stats(1))
    WRITE(stdout, *) '  iii) grid-step too large?: ', NINT(stats(2))
    WRITE(stdout, *) '  iv)  standard deviation  : ', REAL(stats(3), f_sgle)
    WRITE(stdout, *) '  vi)  mean value          : ', REAL(stats(4), f_sgle)
    WRITE(stdout, *) '  vii) CPU main loop       : ', REAL(stats(5), f_sgle), ' sec'
    WRITE(stdout, *) '  viii)GPU main loop       : ', REAL(stats(6), f_sgle), ' sec'
  ENDIF

  CALL watch_start(tictoc)

  CALL scarf_finalize()

#endif

  ! ================================================================================================================================-
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! all tests done, release memory and exit
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! ================================================================================================================================-

  NULLIFY(v3, x3, y3, z3)
  DEALLOCATE(v1, x1, y1, z1)

  CALL mpi_barrier(mpi_comm_world, ierr)

  IF (rank .eq. 0) THEN
    WRITE(stdout, *) ''
    WRITE(stdout, *) 'All tests completed'
  ENDIF

  CALL mpi_finalize(ierr)

END PROGRAM driver
