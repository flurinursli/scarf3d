PROGRAM main

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
  ! You should have received a copy of the GNU General Public License
  ! along with SCARF3D.  If not, see <https://www.gnu.org/licenses/>.
  !
  ! Purpose:
  !   To compute random fields on regular structured meshes based on the SCARF3D package by providing a standalone application. This
  !   application reads input parameters from a file provided by the user. The input file must follow a specific format, i.e.:
  !
  !   # field descriptor (1)
  !   numeric value(s)
  !   # field descriptor (2)
  !   numeric value(s)
  !
  !   If an error occurs, it is printed to screen and then the application stops.
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   04/05/20                  original version
  !   23/07/20                  return 2D/3D fields based on input
  !

  USE, INTRINSIC     :: iso_fortran_env, stdout => output_unit
  USE, NON_INTRINSIC :: mpi
  USE, NON_INTRINSIC :: m_scarflib
  USE, NON_INTRINSIC :: m_scarflib_aux, only: sample_mesh, watch_stop, watch_start

  IMPLICIT none

#ifdef DOUBLE_PREC
  INTEGER, PARAMETER :: f_int = int32
  INTEGER, PARAMETER :: f_real = real64
  INTEGER, PARAMETER :: real_type = mpi_double_precision
#else
  INTEGER, PARAMETER :: f_int = int32
  INTEGER, PARAMETER :: f_real = real32
  INTEGER, PARAMETER :: real_type = mpi_real
#endif
  INTEGER, PARAMETER :: f_sgle = real32
  INTEGER, PARAMETER :: f_dble = real64

  CHARACTER(1)                                  :: num2str
  CHARACTER(256)                                :: buffer
  CHARACTER(:),   ALLOCATABLE                   :: fmt
  INTEGER(f_int)                                :: ierr, rank, ntasks, ndim, npoi
  INTEGER(f_int)                                :: acf, seed, argc, rescale, pad, ios, algorithm
  INTEGER(f_int), ALLOCATABLE, DIMENSION(:)     :: n, fs, fe
  REAL(f_real)                                  :: ds, sigma, hurst, mute, taper, alpha, beta, gamma
  REAL(f_dble)                                  :: tictoc
  REAL(f_real),                DIMENSION(8)     :: stats
  REAL(f_real),   ALLOCATABLE, DIMENSION(:)     :: cl
  REAL(f_real),   ALLOCATABLE, DIMENSION(:,:)   :: poi, field2d
  REAL(f_real),   ALLOCATABLE, DIMENSION(:,:,:) :: field3d

  !---------------------------------------------------------------------------------------------------------------------------------

  ! initialise mpi
  CALL mpi_init(ierr)

  ! get rank number
  CALL mpi_comm_rank(mpi_comm_world, rank, ierr)

  ! get number of tasks
  CALL mpi_comm_size(mpi_comm_world, ntasks, ierr)

  IF (rank .eq. 0) THEN
    WRITE(stdout, '(A)') "This program comes with ABSOLUTELY NO WARRANTY; released under GPL."
    WRITE(stdout, '(A)') "This is free software, and you are welcome to redistribute it"
    WRITE(stdout, '(A)') "under certain conditions; see LICENSE.txt for more details."
    WRITE(stdout, '(A)') ""
  ENDIF

  CALL mpi_barrier(mpi_comm_world, ierr)

  ! set mandatory parameters to wrong values in order to trigger errors if not defined by user
  ndim  = 0
  ds    = 0._f_real
  hurst = -1._f_real
  acf   = -1
  sigma = -1._f_real
  seed  = 0

  ! set value for optional parameters
  mute      = 0._f_real
  taper     = 0._f_real
  rescale   = 0
  pad       = 0
  alpha     = 0._f_real
  beta      = 0._f_real
  gamma     = 0._f_real
  algorithm = 0

  ! ================================================================================================================================
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! retrieve input parameters
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! ================================================================================================================================

  ! process command line
  argc = command_argument_count()

  ! print error in case of no input file
  IF (argc .ne. 1) THEN
    IF (rank .eq. 0) WRITE(stdout, *) 'Please specify an input file!'
    CALL mpi_barrier(mpi_comm_world, ierr)
    CALL mpi_finalize(ierr)
    RETURN
  ENDIF

  ierr = 0
  ndim = 0

  ! only master reads input parameters
  IF (rank .eq. 0) THEN

    ! get file name
    CALL get_command_argument(1, buffer)

    OPEN(1, file = TRIM(buffer), status = 'unknown', form = 'formatted')

    ios = 0

    ! get number of dimension
    DO WHILE(ios .ne. -1)
      READ(1, '(A100)', iostat = ios) buffer
      IF (TRIM(buffer) .eq. '# dimension') READ(1, *) ndim
    ENDDO

    ! keep reading only if ndim = 2 or 3
    IF ( (ndim .eq. 2) .or. (ndim .eq. 3) ) THEN

      npoi = 0
      ios  = 0

      REWIND(1)

      ! get number of POI
      DO WHILE(ios .ne. -1)
        READ(1, '(A100)', iostat = ios) buffer
        IF (TRIM(buffer) .eq. '# poi') npoi = npoi + 1
      ENDDO

      ALLOCATE(n(ndim), cl(ndim), poi(ndim, npoi))

      n(:)     = 0
      cl(:)    = 0._f_real
      poi(:,:) = -1._f_real

      ios  = 0
      npoi = 0

      REWIND(1)

      ! start reading parameters
      DO WHILE(ios .ne. -1)

        READ(1, '(A100)', iostat = ios) buffer

        SELECT CASE(TRIM(buffer))
          CASE('# npts')
            READ(1, *) n
          CASE('# poi')
            npoi = npoi + 1
            READ(1, *) poi(:, npoi)
          CASE('# ds')
            READ(1, *) ds
          CASE('# cl')
            READ(1, *) cl
          CASE('# sigma')
            READ(1, *) sigma
          CASE('# acf')
            READ(1, *) acf
          CASE('# hurst')
            READ(1, *) hurst
          CASE('# seed')
            READ(1, *) seed
          CASE('# mute')
            READ(1, *) mute
          CASE('# taper')
            READ(1, *) taper
          CASE('# rescale')
            READ(1, *) rescale
          CASE('# pad')
            READ(1, *) pad                               !< pad model to avoid FFT wrap-around
          CASE('# alpha')
            READ(1, *) alpha
          CASE('# beta')
            READ(1, *) beta
          CASE('# gamma')
            READ(1, *) gamma
          CASE('# algorithm')
            READ(1, *) algorithm
        END SELECT

      ENDDO

      CLOSE(1)

    ENDIF

  ENDIF

  CALL mpi_barrier(mpi_comm_world, ierr)

  ! broadcast number of dimensions and POI
  CALL mpi_bcast(ndim, 1, mpi_integer, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(npoi, 1, mpi_integer, 0, mpi_comm_world, ierr)

  ! raise error if dimension not equal 2 or 3
  IF ( (ndim .ne. 2) .and. (ndim .ne. 3) ) THEN
    IF (rank .eq. 0) WRITE(stdout, *) 'Error: dimension must be equal to 2 or 3'
    CALL mpi_barrier(mpi_comm_world, ierr)
    CALL mpi_abort(mpi_comm_world, 1, ierr)
  ENDIF

  IF (.not.ALLOCATED(cl)) ALLOCATE(n(ndim), cl(ndim), poi(ndim, npoi))

  ! broadcast remaining parameters
  CALL mpi_bcast(n, ndim, mpi_integer, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(algorithm, 1, mpi_integer, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(acf, 1, mpi_integer, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(seed, 1, mpi_integer, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(rescale, 1, mpi_integer, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(mute, 1, real_type, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(taper, 1, real_type, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(poi, ndim*npoi, real_type, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(ds, 1, real_type, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(cl, ndim, real_type, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(sigma, 1, real_type, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(hurst, 1, real_type, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(alpha, 1, real_type, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(beta, 1, real_type, 0, mpi_comm_world, ierr)
  CALL mpi_bcast(gamma, 1, real_type, 0, mpi_comm_world, ierr)

  ! ================================================================================================================================
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! raise error if input parameters are not correct
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! ================================================================================================================================

  IF (ANY(n .le. 0)) THEN
    IF (rank .eq. 0) WRITE(stdout, *) 'Error: number of points must be larger than zero'
    CALL mpi_barrier(mpi_comm_world, ierr)
    CALL mpi_abort(mpi_comm_world, 1, ierr)
  ENDIF

  IF (ANY(poi .lt. 0._f_real)) THEN
    IF (rank .eq. 0) WRITE(stdout, *) 'Error: POI must be inside the model'
    CALL mpi_barrier(mpi_comm_world, ierr)
    CALL mpi_abort(mpi_comm_world, 1, ierr)
  ENDIF

  IF ( (acf .ne. 0) .and. (acf .ne. 1) ) THEN
    IF (rank .eq. 0) WRITE(stdout, *) 'Error: accepted ACF are Von Karman (0) and Gaussian (1)'
    CALL mpi_barrier(mpi_comm_world, ierr)
    CALL mpi_abort(mpi_comm_world, 1, ierr)
  ENDIF

  IF (mute .lt. 0._f_real) THEN
    IF (rank .eq. 0) WRITE(stdout, *) 'Error: mute parameter can be only positive'
    CALL mpi_barrier(mpi_comm_world, ierr)
    CALL mpi_abort(mpi_comm_world, 1, ierr)
  ENDIF

  IF (taper .lt. 0._f_real) THEN
    IF (rank .eq. 0) WRITE(stdout, *) 'Error: taper parameter can be only positive'
    CALL mpi_barrier(mpi_comm_world, ierr)
    CALL mpi_abort(mpi_comm_world, 1, ierr)
  ENDIF

  IF ( (rescale .ne. 0) .and. (rescale .ne. 1) ) THEN
    IF (rank .eq. 0) WRITE(stdout, *) 'Error: rescale parameter can be only off (0) or on (1)'
    CALL mpi_barrier(mpi_comm_world, ierr)
    CALL mpi_abort(mpi_comm_world, 1, ierr)
  ENDIF

  IF (sigma .le. 0._f_real) THEN
    IF (rank .eq. 0) WRITE(stdout, *) 'Error: standard deviation (sigma) must be larger than zero'
    CALL mpi_barrier(mpi_comm_world, ierr)
    CALL mpi_abort(mpi_comm_world, 1, ierr)
  ENDIF

  IF (ds .le. 0._f_real) THEN
    IF (rank .eq. 0) WRITE(stdout, *) 'Error: grid-step must be larger than zero'
    CALL mpi_barrier(mpi_comm_world, ierr)
    CALL mpi_abort(mpi_comm_world, 1, ierr)
  ENDIF

  IF ( (acf .eq. 0) .and. (hurst .lt. 0._f_real) .and. (hurst .gt. 1._f_real) ) THEN
    IF (rank .eq. 0) WRITE(stdout, *) 'Error: Hurst exponent must be in the range [0 1]'
    CALL mpi_barrier(mpi_comm_world, ierr)
    CALL mpi_abort(mpi_comm_world, 1, ierr)
  ENDIF

  IF (seed .eq. 0) THEN
    IF (rank .eq. 0) WRITE(stdout, *) 'Error: seed number not specified'
    CALL mpi_barrier(mpi_comm_world, ierr)
    CALL mpi_abort(mpi_comm_world, 1, ierr)
    STOP
  ENDIF

  ! ================================================================================================================================
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! make summary
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! ================================================================================================================================

  IF (rank .eq. (ntasks - 1)) THEN
    WRITE(stdout, *)
    WRITE(stdout, *) 'Summary of input parameters'
    WRITE(stdout, *) '****************************************************************'
    buffer = 'Number of dimensions'
    WRITE(stdout, '(X, A26, A, I12, T65, A)') ADJUSTL(buffer), '|', ndim,      '|'
    WRITE(stdout, *) '--------------------------|------------------------------------|'
    buffer = 'Algorithm'
    WRITE(stdout, '(X, A26, A, I12, T65, A)') ADJUSTL(buffer), '|', algorithm, '|'
    WRITE(stdout, *) '--------------------------|------------------------------------|'
    WRITE(num2str, '(I0)') ndim; fmt = '(X, A26, A, ' // num2str // 'I12, T65, A)'
    buffer = 'Number of points'
    WRITE(stdout, (fmt))                        ADJUSTL(buffer), '|', n, '|'
    WRITE(stdout, *) '--------------------------|------------------------------------|'
    buffer = 'Grid-step'
    WRITE(stdout, '(X, A26, A, F12.3, T65, A)') ADJUSTL(buffer), '|', ds,  '|'
    WRITE(stdout, *) '--------------------------|------------------------------------|'
    buffer = 'Auto-correlation function'
    WRITE(stdout, '(X, A26, A, I12, T65, A)')   ADJUSTL(buffer), '|', acf, '|'

    fmt = '(X, A26, A, ' // num2str // 'F12.3, T65, A)'
    buffer = 'Correlation length'
    WRITE(stdout, (fmt))                        ADJUSTL(buffer), '|', cl,    '|'
    buffer = 'Std. dev.'
    WRITE(stdout, '(X, A26, A, F12.3, T65, A)') ADJUSTL(buffer), '|', sigma, '|'

    IF (acf .eq. 0) THEN
      buffer = 'Hurst exponent'
      WRITE(stdout, '(X, A26, A, F12.3, T65, A)') ADJUSTL(buffer), '|', hurst, '|'
      WRITE(stdout, *) '--------------------------|------------------------------------|'
    ENDIF

    buffer = 'Seed number'
    WRITE(stdout, '(X, A26, A, I12, T65, A)') ADJUSTL(buffer), '|', seed, '|'
    WRITE(stdout, *) '--------------------------|------------------------------------|'

    IF (npoi .gt. 0) THEN
      buffer = 'Number of POIs'
      WRITE(stdout, '(X, A26, A, I12, T65, A)')   ADJUSTL(buffer), '|', npoi,  '|'
      buffer = 'Muting radius'
      WRITE(stdout, '(X, A26, A, F12.3, T65, A)') ADJUSTL(buffer), '|',  mute,  '|'
      buffer = 'Tapering radius'
      WRITE(stdout, '(X, A26, A, F12.3, T65, A)') ADJUSTL(buffer), '|',  taper, '|'
      WRITE(stdout, *) '--------------------------|------------------------------------|'
    ENDIF

    buffer = 'Alpha angle'
    WRITE(stdout, '(X, A26, A, F12.3, T65, A)') ADJUSTL(buffer), '|', alpha, '|'
    buffer = 'Beta angle'
    WRITE(stdout, '(X, A26, A, F12.3, T65, A)') ADJUSTL(buffer), '|', beta,  '|'
    buffer = 'Gamma angle'
    WRITE(stdout, '(X, A26, A, F12.3, T65, A)') ADJUSTL(buffer), '|', gamma, '|'

    WRITE(stdout, *) '****************************************************************'
    WRITE(stdout, *)

    FLUSH(stdout)

  ENDIF

  CALL mpi_barrier(mpi_comm_world, ierr)

  ! ================================================================================================================================
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! create structured mesh
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! ================================================================================================================================

  ALLOCATE(fs(ndim), fe(ndim))

  CALL sample_mesh(rank, ntasks, n, fs, fe)

  IF (ndim .eq. 2) THEN
    ALLOCATE(field2d(fe(1) - fs(1) + 1, fe(2) - fs(2) + 1))
  ELSE
    ALLOCATE(field3d(fe(1) - fs(1) + 1, fe(2) - fs(2) + 1, fe(3) - fs(3) + 1))
  ENDIF

  ! ================================================================================================================================
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! compute random field
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! ================================================================================================================================

  CALL scarf_initialize(fs, fe, ds, acf, cl, sigma, method = algorithm, hurst = hurst, alpha = alpha, beta = beta, gamma = gamma, &
                        poi = poi, taper = taper, mute = mute)

  CALL watch_start(tictoc)

  IF (ndim .eq. 2) THEN
    CALL scarf_execute(seed, field2d, stats)
  ELSE
    CALL scarf_execute(seed, field3d, stats)
  ENDIF

  CALL watch_stop(tictoc)

  CALL mpi_allreduce(mpi_in_place, stats, 8, mpi_real, mpi_max, mpi_comm_world, ierr)

  IF (rank .eq. 0) THEN
    WRITE(stdout, *) ''
    buffer = 'Statistics for current simulation'
    WRITE(stdout, '(X, A35)') ADJUSTL(buffer)
    WRITE(stdout, *) '*************************************************'
    buffer = 'Elapsed time'
    WRITE(stdout, '(X, A30, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', tictoc, ' sec', '|'
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

  ! ================================================================================================================================
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! write random field to disk
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! ================================================================================================================================

  CALL watch_start(tictoc)

  IF (ndim .eq. 2) THEN
    CALL scarf_io(n, field2d, 'random_field.bin', 3)
  ELSE
    CALL scarf_io(n, field3d, 'random_field.bin', 3)
  ENDIF

  CALL watch_stop(tictoc)

  IF (rank .eq. 0) THEN
    buffer = 'I/O time'
    WRITE(stdout, '(X, A30, A, F12.5, A, T50, A)') ADJUSTL(buffer), '|', tictoc, ' sec', '|'
    WRITE(stdout, *) '*************************************************'
    WRITE(stdout, *)
  ENDIF

  ! ================================================================================================================================
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! cleanup
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! ================================================================================================================================

  CALL scarf_finalize()

  CALL mpi_barrier(mpi_comm_world, ierr)

  IF (ndim .eq. 2) THEN
    DEALLOCATE(field2d)
  ELSE
    DEALLOCATE(field3d)
  ENDIF

  DEALLOCATE(cl, poi, n, fs, fe)

  IF (rank .eq. 0) THEN
    WRITE(stdout, *) ''
    WRITE(stdout, *) 'Program completed successfully'
  ENDIF

  CALL mpi_finalize(ierr)


END PROGRAM main
