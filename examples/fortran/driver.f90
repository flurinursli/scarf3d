PROGRAM driver

  ! Copyright (c) 2020, Eidgenoessische Technische Hochschule Zurich, ETHZ.
  !
  ! Written by:
  ! Walter Imperatori (walter.imperatori@sed.ethz.ch)
  !
  ! All rights reserved.
  !
  ! This file is part of SCARF3D, version: 2.4
  !
  ! SCARF3D is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU General Public License as published by
  ! the Free Software Foundation, either version 3 of the License, or
  ! (at your option) any later version.
  !
  ! SCARF3D is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ! GNU General Public License for more details.
  !
  ! Purpose:
  !   To provide a sample FORTRAN driver program for the SCARF3D library
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   04/05/20                  original version
  !

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
  REAL(f_real)                                           :: ds, dh, sigma, hurst, mute, taper, alpha, beta
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

  IF (rank .eq. 0) THEN
    WRITE(stdout, *) ''
    WRITE(stdout, *) 'This program will compute small 2D and 3D random fields'
    WRITE(stdout, *) 'using the FIM (and SRM - if selected at compile-time) algorithm(s).'
    WRITE(stdout, *) 'Test should complete in a few minutes. Output can be inspected'
    WRITE(stdout, *) 'with Matlab/Octave script "scarf3d.m" and "vizme.m".'
    WRITE(stdout, *) ''
  ENDIF

  !---------------------------------------------------------------------------------------------------------------------------------
  ! input section
  !---------------------------------------------------------------------------------------------------------------------------------

  ! number of points for whole model
  n = [500, 500, 500]

  ! grid step
  ds = 50._f_real

  dh = 50._f_real

  ! autocorrelation (0=von karman/exponential, 1=gaussian, 2=user-defined)
  acf = 0

  ! correlation length
  cl = [2000._f_real, 2000._f_real, 2000._f_real]

  ! standard deviation (sigma(%)/100)
  sigma = 0.05_f_real

  ! hurst exponent (not used for gaussian acf)
  hurst = 0.25_f_real

  ! seed number
  seed = 1235

  ! set position of point-of-interest (muting/tapering), same units as "dh"
  poi(:, 1) = [400., 250., 100.] * ds
  poi(:, 2) = [200., 150.,  50.] * ds

  ! radius for muting (at poi), same units as "dh"
  mute = 10 * ds

  ! radius for tapering (at poi + mute), same units as "dh"
  taper = 50 * ds

  ! rescale to desired (continuous) sigma
  rescale = 0

  ! expand grid to handle fft periodicity
  pad = 0

  ! rotation angles
  alpha = 30._f_real
  beta  = -10._f_real

  ! end input section
  !---------------------------------------------------------------------------------------------------------------------------------

  ! ================================================================================================================================
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! create sample cartesian mesh: 2D case
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
    WRITE(stdout, *) '****** FIM algorithm, 2D, cartesian mesh *******'
  ENDIF

  CALL scarf_initialize(fs(1:2), fe(1:2), ds, acf, cl(1:2), sigma, method = 0, hurst = hurst, alpha = alpha,     &
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

  CALL scarf_io(n(1:2), v2, 'fim_cart_whole_2d', 2)

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
    WRITE(stdout, *) '***** FIM algorithm, 2D, non-cartesian mesh ******'
  ENDIF

  CALL scarf_initialize(dh, acf, cl(1:2), sigma, x1, y1, method = 0,  &
                        hurst = hurst, alpha = alpha, poi = poi(1:2,:), taper = taper, mute = mute)

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

  CALL scarf_io(n(1:2), v2, 'fim_nocart_whole_2d', 2)

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
    WRITE(stdout, *) '****** SRM algorithm, 2D, cartesian mesh *******'
  ENDIF

  ! cartesian mesh test
  CALL scarf_initialize(fs(1:2), fe(1:2), ds, acf, cl(1:2), sigma, method = 1, hurst = hurst, alpha = alpha,    &
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

  CALL scarf_io(n(1:2), v2, 'srm_cart_whole_2d', 2)

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
    WRITE(stdout, *) '***** SRM algorithm, 2D, non-cartesian mesh ******'
  ENDIF

  ! non-cartesian mesh test
  CALL scarf_initialize(dh, acf, cl(1:2), sigma, x1, y1, method = 1,  &
                        hurst = hurst, alpha = alpha, poi = poi(1:2,:), taper = taper, mute = mute)

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

  CALL scarf_io(n(1:2), v2, 'srm_nocart_whole_2d', 2)

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
  ! create sample cartesian mesh: 3D case
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
    WRITE(stdout, *) '****** FIM algorithm, 3D, cartesian mesh *******'
  ENDIF

  CALL scarf_initialize(fs, fe, ds, acf, cl, sigma, method = 0, hurst = hurst)
  ! CALL scarf_initialize(fs, fe, ds, acf, cl, sigma, method = 0, hurst = hurst, alpha = alpha, beta = beta, poi = poi,   &
  !                       taper = taper, mute = mute)

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

  CALL scarf_io(n, "x", n(1)/2, v3, "fim_cart_xslice")
  CALL scarf_io(n, "y", n(2)/2, v3, "fim_cart_yslice")
  CALL scarf_io(n, "z", n(3)/2, v3, "fim_cart_zslice")

  CALL watch_stop(tictoc)

  IF (rank .eq. 0) THEN
    buffer = 'I/O time (slices)'
    WRITE(stdout, '(X, A30, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', tictoc, ' sec', '|'
  ENDIF

  CALL watch_start(tictoc)

  CALL scarf_io(n, v3, 'fim_cart_whole_3d', 2)

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
    WRITE(stdout, *) '***** FIM algorithm, 3D, non-cartesian mesh ******'
  ENDIF

  CALL scarf_initialize(dh, acf, cl, sigma, x1, y1, z1, method = 0, hurst = hurst)
  ! CALL scarf_initialize(dh, acf, cl, sigma, x1, y1, z1, method = 0, hurst = hurst, alpha = alpha, beta = beta, poi = poi,  &
  !                       taper = taper, mute = mute)

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

  CALL scarf_io(n, "x", n(1)/2, v3, "fim_nocart_xslice")
  CALL scarf_io(n, "y", n(2)/2, v3, "fim_nocart_yslice")
  CALL scarf_io(n, "z", n(3)/2, v3, "fim_nocart_zslice")

  CALL watch_stop(tictoc)

  IF (rank .eq. 0) THEN
    buffer = 'I/O time (slices)'
    WRITE(stdout, '(X, A30, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', tictoc, ' sec', '|'
  ENDIF

  CALL watch_start(tictoc)

  CALL scarf_io(n, v3, 'fim_nocart_whole_3d', 2)

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
    WRITE(stdout, *) '****** SRM algorithm, 3D, cartesian mesh *******'
  ENDIF

  CALL scarf_initialize(fs, fe, ds, acf, cl, sigma, method = 1, hurst = hurst)
  ! CALL scarf_initialize(fs, fe, ds, acf, cl, sigma, method = 1, hurst = hurst, alpha = alpha, beta = beta, poi=poi, taper=taper,  &
  !                       mute=mute)

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

  CALL scarf_io(n, "x", n(1)/2, v3, "srm_cart_xslice")
  CALL scarf_io(n, "y", n(2)/2, v3, "srm_cart_yslice")
  CALL scarf_io(n, "z", n(3)/2, v3, "srm_cart_zslice")

  CALL watch_stop(tictoc)

  IF (rank .eq. 0) THEN
    buffer = 'I/O time (slices)'
    WRITE(stdout, '(X, A30, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', tictoc, ' sec', '|'
  ENDIF

  CALL watch_start(tictoc)

  CALL scarf_io(n, v3, 'srm_cart_whole_3d', 2)

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
    WRITE(stdout, *) '***** SRM algorithm, 3D, non-cartesian mesh ******'
  ENDIF

  CALL scarf_initialize(dh, acf, cl, sigma, x1, y1, z1, method = 1, hurst = hurst)
  ! CALL scarf_initialize(dh, acf, cl, sigma, x1, y1, z1, method = 1, hurst = hurst, alpha = alpha, beta = beta, poi = poi,   &
  !                       taper = taper, mute = mute)

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

  CALL scarf_io(n, "x", n(1)/2, v3, "srm_nocart_xslice")
  CALL scarf_io(n, "y", n(2)/2, v3, "srm_nocart_yslice")
  CALL scarf_io(n, "z", n(3)/2, v3, "srm_nocart_zslice")

  CALL watch_stop(tictoc)

  IF (rank .eq. 0) THEN
    buffer = 'I/O time (slices)'
    WRITE(stdout, '(X, A30, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', tictoc, ' sec', '|'
  ENDIF

  CALL watch_start(tictoc)

  CALL scarf_io(n, v3, 'srm_nocart_whole_3d', 2)

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
