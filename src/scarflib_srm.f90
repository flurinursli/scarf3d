MODULE m_scarflib_srm

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
  !   To compute a random field according to the Spectral Representation Method (SRM). Subprograms rely on the TRNG (version 4) library
  !   to generate random numbers.
  !   Subroutines 'scarf3d_unstructured_srm' and 'scarf3d_structured_srm' are virtually indentical, with only few differences sparse
  !   at different points in the code (this inhibits somehow the use of 'include' statements to handle a single code)
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   04/05/20                  original version
  !   11/05/20                  updated macro for double-precision
  !   21/07/20                  corr.len. refers to arbitrary axis
  !   23/07/20                  return 2D/3D fields based on input
  !

  USE, NON_INTRINSIC :: m_scarflib_common
  USE, NON_INTRINSIC :: m_psdf

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  PRIVATE

  PUBLIC :: scarf3d_srm

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTERFACE scarf3d_srm
    MODULE PROCEDURE scarf3d_unstructured_srm, scarf3d_structured_srm
  END INTERFACE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  ! interface to c++ functions/subroutines
  INTERFACE

    SUBROUTINE srng(rand) BIND(c, name="srng")
      USE, INTRINSIC     :: iso_c_binding
      USE, NON_INTRINSIC :: m_scarflib_common
      IMPLICIT none
      REAL(c_real), DIMENSION(2), INTENT(OUT) :: rand
    END SUBROUTINE srng

    SUBROUTINE normdist(rand) BIND(c, name="normdist")
      USE, INTRINSIC     :: iso_c_binding
      USE, NON_INTRINSIC :: m_scarflib_common
      IMPLICIT none
      REAL(c_real), DIMENSION(2), INTENT(OUT) :: rand
    END SUBROUTINE normdist

    SUBROUTINE set_seed(seed) BIND(c, name="set_seed")
      USE, INTRINSIC :: iso_c_binding
      IMPLICIT none
      INTEGER(c_int), value :: seed
    END SUBROUTINE set_seed

  END INTERFACE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! number of dimensions (2 or 3)
  INTEGER(f_int)            :: ndim

  ! number of harmonics in spectral summation
  INTEGER(f_int), PARAMETER :: nharm = 10000

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! pointer to PSD function for integration
  PROCEDURE(fn_vk), POINTER :: fun

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  !=================================================================================================================================
  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  SUBROUTINE scarf3d_unstructured_srm(nc, fc, ds, x, y, z, dh, acf, cl, sigma, hurst, seed, poi, mute, taper, rescale, alpha,  &
                                      beta, gamma, field, info)

    ! Purpose:
    ! To compute a random field characterised by a selected ACF on unstructured meshes. For a given ACF, the actual distribution of
    ! heterogeneity may be influenced by the size of the external mesh (determined by the near-corner 'nc' and far-corner 'fc') and
    ! the minimum desired grid-step 'dh'. Harmonics contribution outside the corrisponding interval [kmin kmax] will be rejected if
    ! minimum and actual external grid-steps do not coincide.
    !
    ! Revisions:
    !     Date                    Description of change
    !     ====                    =====================
    !   04/05/20                  original version
    !   21/07/20                  rotation angles
    !   23/07/20                  2D fields
    !   08/09/20                  sigma normalisation
    !

    REAL(f_real),                DIMENSION(3),   INTENT(IN)  :: nc, fc       !< min/max extent of external grid, control min wavenumber
    REAL(f_real),                                INTENT(IN)  :: ds           !< current grid-step of external grid, control max resolvable wavenumber
    REAL(f_real),                DIMENSION(:),   INTENT(IN)  :: x, y, z      !< position of grid points along x, y, z
    REAL(f_real),                                INTENT(IN)  :: dh           !< minimum grid-step of external grid, control max wavenumber
    INTEGER(f_int),                              INTENT(IN)  :: acf          !< autocorrelation function: 0=vk, 1=gauss
    REAL(f_real),                DIMENSION(3),   INTENT(IN)  :: cl           !< correlation length
    REAL(f_real),                                INTENT(IN)  :: sigma        !< standard deviation
    REAL(f_real),                                INTENT(IN)  :: hurst        !< hurst exponent
    INTEGER(f_int),                              INTENT(IN)  :: seed         !< seed number
    REAL(f_real),                DIMENSION(:,:), INTENT(IN)  :: poi          !< location of point(s)-of-interest where muting/tapering is applied
    REAL(f_real),                                INTENT(IN)  :: mute         !< radius for muting
    REAL(f_real),                                INTENT(IN)  :: taper        !< radius for tapering
    INTEGER(f_int),                              INTENT(IN)  :: rescale      !< flag for rescaling random field to desired sigma
    REAL(f_real),                                INTENT(IN)  :: alpha        !< angle about z-axis
    REAL(f_real),                                INTENT(IN)  :: beta         !< angle about y-axis
    REAL(f_real),                                INTENT(IN)  :: gamma        !< angle about x-axis
    REAL(f_real),                DIMENSION(:),   INTENT(OUT) :: field        !< random field at x,y,z location
    REAL(f_real),                DIMENSION(6),   INTENT(OUT) :: info         !< errors flags and TIMING for performance analysis
    CHARACTER(30)                                            :: string
    INTEGER(f_int)                                           :: i, l
    INTEGER(f_int)                                           :: npts
    INTEGER(f_int)                                           :: ierr
    INTEGER(c_int)                                           :: c_seed
    INTEGER(f_int), ALLOCATABLE, DIMENSION(:)                :: npoints
    REAL(f_real)                                             :: scaling
    REAL(f_real)                                             :: kmax, kmin, kc
    REAL(f_real)                                             :: k, d, const
    REAL(f_real)                                             :: phi, theta
    REAL(f_real)                                             :: v1, v2, v3
    REAL(f_real)                                             :: a, b, arg
    REAL(f_real)                                             :: calpha, salpha
    REAL(f_real)                                             :: cbeta, sbeta
    REAL(f_real)                                             :: cgamma, sgamma
    REAL(f_real)                                             :: hurst_factor
    REAL(f_dble)                                             :: tictoc
    REAL(c_real),                DIMENSION(2)                :: r
    REAL(f_real),                DIMENSION(3)                :: span, pt
    REAL(f_real),                DIMENSION(3)                :: min_extent, max_extent
    REAL(f_real),                DIMENSION(3)                :: bar, obar
    REAL(f_real),   ALLOCATABLE, DIMENSION(:)                :: mu, var
    REAL(f_real),                DIMENSION(3,3)              :: matrix

    !-------------------------------------------------------------------------------------------------------------------------------

    CALL mpi_comm_SIZE(mpi_comm_world, world_size, ierr)
    CALL mpi_comm_rank(mpi_comm_world, world_rank, ierr)

#ifdef DEBUG
    IF (world_rank .eq. 0) THEN
      WRITE(output_unit, *) 'DEBUG MODE: input parameters for "scarf3d_unstructured_srm"'
      WRITE(output_unit, *) '****************************************************************'
      string = 'near corner (NC)'
      WRITE(output_unit, '(X, A26, A, 3F12.3, T65, A)') ADJUSTL(string), '|', nc, '|'
      string = 'far corner (FC)'
      WRITE(output_unit, '(X, A26, A, 3F12.3, T65, A)') ADJUSTL(string), '|', fc, '|'
      string = '(maximum) grid-step (DS)'
      WRITE(output_unit, '(X, A26, A, F12.3, T65, A)') ADJUSTL(string), '|', ds, '|'
      string = 'number of points (X-Y-Z)'
      WRITE(output_unit, '(X, A26, A, 3I12, T65, A)') ADJUSTL(string), '|', SIZE(x), SIZE(y), SIZE(z), '|'
      string = 'internal grid-step (DH)'
      WRITE(output_unit, '(X, A26, A, F12.3, T65, A)') ADJUSTL(string), '|', dh, '|'
      string = 'autocorr. fun. (ACF)'
      WRITE(output_unit, '(X, A26, A, I12, T65, A)') ADJUSTL(string), '|', acf, '|'
      string = 'correlation length (CL)'
      WRITE(output_unit, '(X, A26, A, 3F12.3, T65, A)') ADJUSTL(string), '|', cl, '|'
      string = 'std. dev. (SIGMA)'
      WRITE(output_unit, '(X, A26, A, F12.3, T65, A)') ADJUSTL(string), '|', sigma, '|'
      string = 'hurst exp. (HURST)'
      WRITE(output_unit, '(X, A26, A, F12.3, T65, A)') ADJUSTL(string), '|', hurst, '|'
      string = 'seed number (SEED)'
      WRITE(output_unit, '(X, A26, A, I12, T65, A)') ADJUSTL(string), '|', seed, '|'
      string = 'number of POIs'
      WRITE(output_unit, '(X, A26, A, I12, T65, A)') ADJUSTL(string), '|', SIZE(poi, 2), '|'
      string = 'muting radius (MUTE)'
      WRITE(output_unit, '(X, A26, A, F12.3, T65, A)') ADJUSTL(string), '|', mute, '|'
      string = 'tapering radius (TAPER)'
      WRITE(output_unit, '(X, A26, A, F12.3, T65, A)') ADJUSTL(string), '|', taper, '|'
      string = 'rescaling opt. (RESCALE)'
      WRITE(output_unit, '(X, A26, A, L12, T65, A)') ADJUSTL(string), '|', rescale, '|'
      string = 'alpha angle (ALPHA)'
      WRITE(output_unit, '(X, A26, A, F12.3, T65, A)') ADJUSTL(string), '|', alpha, '|'
      string = 'beta angle (BETA)'
      WRITE(output_unit, '(X, A26, A, F12.3, T65, A)') ADJUSTL(string), '|', beta, '|'
      string = 'gamma angle (GAMMA)'
      WRITE(output_unit, '(X, A26, A, F12.3, T65, A)') ADJUSTL(string), '|', gamma, '|'
      WRITE(output_unit, *) '****************************************************************'
    ENDIF

    CALL mpi_barrier(mpi_comm_world, ierr)
#endif

    info(:) = 0._f_real

    ! discriminate between 2D and 3D case
    IF (cl(3) .gt. 0._f_real) THEN
      ndim = 3
    ELSE
      ndim = 2
    ENDIF

    ! initialise random generator
    c_seed = seed
    CALL set_seed(c_seed)

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    ! set min/max radial wavenumber such that "v1", "v2", "v3" below are in the range ["kmin" "kmax"] when the trigonometric terms
    ! equal 1. This is achieved by taking the largest "kmin*cl" and the smallest "kmax*cl". "cl" is considered in the calculations
    ! because the radial avenumber "k" is divided by "cl" (see Eq.8 of Raess et al., 2019)
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    ! model limits (for each process) along each axis
    min_extent = nc
    max_extent = fc

    ! (global) model limits
    CALL mpi_allreduce(mpi_in_place, min_extent, 3, real_type, mpi_min, mpi_comm_world, ierr)
    CALL mpi_allreduce(mpi_in_place, max_extent, 3, real_type, mpi_max, mpi_comm_world, ierr)

    ! (global) effective grid size
    span = max_extent - min_extent

    kmin = 2._f_real * pi * maxval(cl(1:ndim) / span(1:ndim))

    kmax = pi / dh * minval(cl(1:ndim))

    ! corner wavenumber for filtering spectrum is controlled by external mesh grid-step
    kc = pi / ds * minval(cl(1:ndim))

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    ! check if the following conditions for a correct wavenumber representation are violated:
    ! 1) kmin < 1/cl
    ! 2) kmax > 1/cl, or ds <= cl / 2 (frenje & juhlin, 2000)
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    IF (ANY(2._f_real * pi / span(1:ndim) .ge. 1._f_real / cl(1:ndim))) info(1) = 1._f_real

    IF (ANY(ds .gt. cl(1:ndim) / 2._f_real)) info(2) = 1._f_real

    ! set scaling factor
    scaling = sigma / SQRT(REAL(nharm, f_real))

    ! number of points where random field must be calculated
    npts = SIZE(x)

    ! initialise random field
    field(:) = 0._f_real

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    ! compute random field directly at each grid point. USE uniform distribution in the range [0 1] as majorant function.
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    ! angle about z-axis
    calpha = cos(alpha * pi / 180._f_real)
    salpha = sin(alpha * pi / 180._f_real)

    ! angle about y-axis
    cbeta = cos(beta * pi / 180._f_real)
    sbeta = sin(beta * pi / 180._f_real)

    ! angle about x-axis
    cgamma = cos(gamma * pi / 180._f_real)
    sgamma = sin(gamma * pi / 180._f_real)

    matrix(:, 1) = [calpha*cbeta                       , salpha*cbeta                       , -sbeta      ]
    matrix(:, 2) = [calpha*sbeta*sgamma - salpha*cgamma, salpha*sbeta*sgamma + calpha*cgamma, cbeta*sgamma]
    matrix(:, 3) = [calpha*sbeta*cgamma + salpha*sgamma, salpha*sbeta*cgamma - calpha*sgamma, cbeta*cgamma]

    ! baricenter in original reference frame
    bar = (max_extent - min_extent) / 2._f_real

    ! baricenter after rotation
    obar = MATMUL(matrix, bar)

    ! translation vector to rotate around baricenter in original reference frame
    bar = bar - obar

    ! define power term for vk psdf, "nu" is a variable from module "m_psdf"
    nu = hurst + ndim / 2._f_real

    ! set PSD function
    SELECT CASE (acf)
      CASE(0)
        const = 1._f_real
        fun => fn_vk
      CASE(1)
        IF (ndim .eq. 3) const = 1._f_real / SQRT(pi)                         !< scaling factor to be below unity
        fun => fn_gs
      CASE(2)
        const = 1._f_real
        fun => fn_ud
    END SELECT

    ! update scaling factor, assuming that discrete sigma is always equal to continuous sigma if random field is not normalised.
    ! Note that sigma_discr = sigma_cont * sqrt(intg(kmin,kmax) / intg(0, inf))
    scaling = scaling * SQRT(intg(kmin, kmax)/intg())

    ! loop over harmonics
    ! openacc: "field" is copied in&out, all the others are only copied in. "arg" is created locally on the accelerator
    !$acc data copy(field) copyin(x, y, z, scaling, npts, a, b, v1, v2, v3, bar, matrix)
    DO l = 1, nharm

#ifdef TIMING
      CALL watch_start(tictoc, mpi_comm_self)
#endif

      DO

        ! first set of random numbers
        CALL srng(r)

        ! random variable "k`" must span the range [kmin, kmax]
        k = r(1) * (kmax - kmin) + kmin

        ! evaluate original pdf (require input in double precision)
        d = srm_fn(real(k, f_dble)) * const

        ! take "k" and exit if inequality r < d is verified
        IF (r(2) .lt. d) EXIT

      ENDDO

      ! second set of random numbers
      CALL srng(r)

      ! neglect contribution of current harmonic if wavenumber "u" was larger than corner wavenumber "kc" (always false is "ds" = "dh")
      IF (k .gt. kc) CYCLE

      ! compute azimuth and polar angles
      phi = r(1) * 2._f_real * pi

      IF (ndim .eq. 3) THEN                                     !< spherical coordinates
        theta = ACOS(1._f_real - 2._f_real * r(2))
        v1    = k * SIN(phi) * SIN(theta) / cl(1)
        v2    = k * COS(phi) * SIN(theta) / cl(2)
        v3    = k * COS(theta)            / cl(3)
      ELSE                                                   !< polar coordinates
        v1 = k * COS(phi) / cl(1)
        v2 = k * SIN(phi) / cl(2)
        v3 = 0._f_real
      ENDIF

      ! compute harmonics coefficient "a" and "b"
      ! CALL normdist(r)
      !
      ! a = r(1)
      ! b = r(2)

      ! box-mueller transform
      CALL srng(r)
      a = SQRT(-2._f_real * log(r(1))) * COS(2._f_real * pi * r(2))

      CALL srng(r)
      b = SQRT(-2._f_real * log(r(1))) * COS(2._f_real * pi * r(2))

#ifdef TIMING
      CALL watch_stop(tictoc, mpi_comm_self)

      info(5) = info(5) + REAL(tictoc, f_real)

      CALL watch_start(tictoc, mpi_comm_self)
#endif

      !$acc wait

      ! update selected variables in gpu memory
      !$acc update device(a, b, v1, v2, v3)

      ! !$acc parallel loop private(arg) async
      ! DO i = 1, npts
      !   arg = v1 * x(i) + v2 * y(i) + v3 * z(i)
      !   field(i) = field(i) + a * SIN(arg) + b * COS(arg)
      ! ENDDO
      ! !$acc end parallel loop

      !$acc kernels async
      !$acc loop independent private(arg)
      DO i = 1, npts
        !arg = v1 * x(i) + v2 * y(i) + v3 * z(i)
        arg = v1 * ((matrix(1,1) * x(i) + matrix(1,2) * y(i) + matrix(1,3) * z(i)) + bar(1)) +   &
              v2 * ((matrix(2,1) * x(i) + matrix(2,2) * y(i) + matrix(2,3) * z(i)) + bar(2)) +   &
              v3 * ((matrix(3,1) * x(i) + matrix(3,2) * y(i) + matrix(3,3) * z(i)) + bar(3))
        field(i) = field(i) + a * SIN(arg) + b * COS(arg)
      ENDDO
      !$acc end kernels


#ifdef TIMING
      CALL watch_stop(tictoc, mpi_comm_self)

      info(6) = info(6) + REAL(tictoc, f_real)
#endif

    ENDDO
    ! END loop over harmonics

#ifdef TIMING
    CALL watch_start(tictoc, mpi_comm_self)
#endif

    !$acc wait

    ! normalise random field
    ! !$acc parallel loop
    ! DO i = 1, npts
    !   field(i) = field(i) * scaling
    ! ENDDO
    ! !$acc END parallel loop

    !$acc kernels
    !$acc loop independent
    DO i = 1, npts
      field(i) = field(i) * scaling
    ENDDO
    !$acc end kernels

    ! copy "field" back to host and free memory on device
    !$acc end data

#ifdef TIMING
    CALL watch_stop(tictoc, mpi_comm_self)

    info(6) = info(6) + REAL(tictoc, f_real)
#endif

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
    ! compute mean and variance of the whole random field. THEN rescale, if desired, discrete std. dev. to its continuous counterpart
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

    ALLOCATE(var(0:world_size - 1), mu(0:world_size - 1), npoints(0:world_size - 1))

    ! variance and mean for each single process
    var(world_rank) = variance(field)
    mu(world_rank)  = mean(field)

    ! number of grid points per process
    npoints(world_rank) = npts

    ! share results
    CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, var, 1, real_type, mpi_comm_world, ierr)
    CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, mu, 1, real_type, mpi_comm_world, ierr)
    CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, npoints, 1, mpi_integer, mpi_comm_world, ierr)

    ! return total variance in "info(3)" and mean in "info(4)"
    CALL parallel_variance(var, mu, npoints, info(3), info(4))

    ! return standard deviation
    info(3) = SQRT(info(3))

    IF (rescale .eq. 1) THEN

      scaling = sigma / info(3)

      DO i = 1, npts
        field(i) = field(i) * scaling - info(4)
      ENDDO

      info(3) = sigma
      info(4) = 0._f_real

    ENDIF

    DEALLOCATE(var, mu, npoints)

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    ! taper/mute random field at desired locations.
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    pt(:) = 0._f_real

    DO i = 1, SIZE(poi, 2)
      pt(1:ndim) = poi(1:ndim, i)
      CALL tapering(x, y, z, field, pt, mute, taper)
    ENDDO

    CALL mpi_barrier(mpi_comm_world, ierr)

  END SUBROUTINE scarf3d_unstructured_srm

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  !=================================================================================================================================
  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  SUBROUTINE scarf3d_structured_srm(nc, fc, ds, fs, fe, dh, acf, cl, sigma, hurst, seed, poi, mute, taper, rescale, alpha, beta, &
                                    gamma, field, info)

    ! Purpose:
    ! To compute a random field characterised by a selected ACF on unstructured meshes. For a given ACF, the actual distribution of
    ! heterogeneity may be influenced by the size of the external mesh (determined by the near-corner 'nc' and far-corner 'fc') and
    ! the minimum desired grid-step 'dh'. Harmonics contribution outside the corrisponding interval [kmin kmax] will be rejected if
    ! minimum and actual external grid-steps do not coincide.
    !
    ! Revisions:
    !     Date                    Description of change
    !     ====                    =====================
    !   04/05/20                  original version
    !   21/07/20                  rotation angles
    !   23/07/20                  2D fields
    !   08/09/20                  sigma normalisation
    !

    REAL(f_real),                DIMENSION(3),     INTENT(IN)  :: nc, fc          !< min/max extent of external grid, control min wavenumber
    REAL(f_real),                                  INTENT(IN)  :: ds              !< current grid-step of external grid, control max resolvable wavenumber
    INTEGER(f_int),              DIMENSION(3),     INTENT(IN)  :: fs, fe          !< first/last external structured grid index
    REAL(f_real),                                  INTENT(IN)  :: dh              !< minimum grid-step of external grid, control max wavenumber
    INTEGER(f_int),                                INTENT(IN)  :: acf             !< autocorrelation function: 0=vk, 1=gauss
    REAL(f_real),                DIMENSION(3),     INTENT(IN)  :: cl              !< correlation length
    REAL(f_real),                                  INTENT(IN)  :: sigma           !< standard deviation
    REAL(f_real),                                  INTENT(IN)  :: hurst           !< hurst exponent
    INTEGER(f_int),                                INTENT(IN)  :: seed            !< seed number
    REAL(f_real),                DIMENSION(:,:),   INTENT(IN)  :: poi             !< location of point(s)-of-interest where muting/tapering is applied
    REAL(f_real),                                  INTENT(IN)  :: mute            !< radius for mutin
    REAL(f_real),                                  INTENT(IN)  :: taper           !< radius for tapering
    INTEGER(f_int),                                INTENT(IN)  :: rescale         !< flag for rescaling random field to desired sigma
    REAL(f_real),                                  INTENT(IN)  :: alpha           !< angle about z-axis
    REAL(f_real),                                  INTENT(IN)  :: beta            !< angle about y-axis
    REAL(f_real),                                  INTENT(IN)  :: gamma           !< angle about x-axis
    REAL(f_real),                DIMENSION(:,:,:), INTENT(OUT) :: field           !< random field
    REAL(f_real),                DIMENSION(6),     INTENT(OUT) :: info            !< errors flags and TIMING for performance analysis
    CHARACTER(30)                                              :: string
    INTEGER(f_int)                                             :: i, j, k, l
    INTEGER(f_int)                                             :: ierr
    INTEGER(c_int)                                             :: c_seed
    INTEGER(f_int),              DIMENSION(3)                  :: npts
    INTEGER(f_int), ALLOCATABLE, DIMENSION(:)                  :: npoints
    REAL(f_real)                                               :: scaling
    REAL(f_real)                                               :: kmax, kmin, kc
    REAL(f_real)                                               :: u, d, const
    REAL(f_real)                                               :: phi, theta
    REAL(f_real)                                               :: v1, v2, v3
    REAL(f_real)                                               :: a, b, arg
    REAL(f_real)                                               :: x, y, z
    REAL(f_real)                                               :: calpha, salpha
    REAL(f_real)                                               :: cbeta, sbeta
    REAL(f_real)                                               :: cgamma, sgamma
    REAL(f_real)                                               :: hurst_factor
    REAL(f_dble)                                               :: tictoc
    REAL(c_real),                DIMENSION(2)                  :: r
    REAL(f_real),                DIMENSION(3)                  :: span, pt
    REAL(f_real),                DIMENSION(3)                  :: min_extent, max_extent
    REAL(f_real),                DIMENSION(3)                  :: bar, obar
    REAL(f_real),   ALLOCATABLE, DIMENSION(:)                  :: mu, var
    REAL(f_real),                DIMENSION(3,3)                :: matrix

    !-------------------------------------------------------------------------------------------------------------------------------

    CALL mpi_comm_size(mpi_comm_world, world_size, ierr)
    CALL mpi_comm_rank(mpi_comm_world, world_rank, ierr)

#ifdef DEBUG
    IF (world_rank .eq. 0) THEN
      WRITE(output_unit, *) 'DEBUG MODE: input parameters for "scarf3d_structured_srm"'
      WRITE(output_unit, *) '****************************************************************'
      string = 'near corner (NC)'
      WRITE(output_unit, '(X, A26, A, 3F12.3, T65, A)') ADJUSTL(string), '|', nc, '|'
      string = 'far corner (FC)'
      WRITE(output_unit, '(X, A26, A, 3F12.3, T65, A)') ADJUSTL(string), '|', fc, '|'
      string = 'grid-step (DS)'
      WRITE(output_unit, '(X, A26, A, F12.3, T65, A)') ADJUSTL(string), '|', ds, '|'
      string = 'first sample (FS)'
      WRITE(output_unit, '(X, A26, A, 3I12, T65, A)') ADJUSTL(string), '|', fs, '|'
      string = 'last sample (FE)'
      WRITE(output_unit, '(X, A26, A, 3I12, T65, A)') ADJUSTL(string), '|', fe, '|'
      string = 'internal grid-step (DH)'
      WRITE(output_unit, '(X, A26, A, F12.3, T65, A)') ADJUSTL(string), '|', dh, '|'
      string = 'autocorr. fun. (ACF)'
      WRITE(output_unit, '(X, A26, A, I12, T65, A)') ADJUSTL(string), '|', acf, '|'
      string = 'correlation length (CL)'
      WRITE(output_unit, '(X, A26, A, 3F12.3, T65, A)') ADJUSTL(string), '|', cl, '|'
      string = 'std. dev. (SIGMA)'
      WRITE(output_unit, '(X, A26, A, F12.3, T65, A)') ADJUSTL(string), '|', sigma, '|'
      string = 'hurst exp. (HURST)'
      WRITE(output_unit, '(X, A26, A, F12.3, T65, A)') ADJUSTL(string), '|', hurst, '|'
      string = 'seed number (SEED)'
      WRITE(output_unit, '(X, A26, A, I12, T65, A)') ADJUSTL(string), '|', seed, '|'
      string = 'number of POIs'
      WRITE(output_unit, '(X, A26, A, I12, T65, A)') ADJUSTL(string), '|', SIZE(poi, 2), '|'
      string = 'muting radius (MUTE)'
      WRITE(output_unit, '(X, A26, A, F12.3, T65, A)') ADJUSTL(string), '|', mute, '|'
      string = 'tapering radius (TAPER)'
      WRITE(output_unit, '(X, A26, A, F12.3, T65, A)') ADJUSTL(string), '|', taper, '|'
      string = 'rescaling opt. (RESCALE)'
      WRITE(output_unit, '(X, A26, A, L12, T65, A)') ADJUSTL(string), '|', rescale, '|'
      string = 'alpha angle (ALPHA)'
      WRITE(output_unit, '(X, A26, A, F12.3, T65, A)') ADJUSTL(string), '|', alpha, '|'
      string = 'beta angle (BETA)'
      WRITE(output_unit, '(X, A26, A, F12.3, T65, A)') ADJUSTL(string), '|', beta, '|'
      string = 'gamma angle (GAMMA)'
      WRITE(output_unit, '(X, A26, A, F12.3, T65, A)') ADJUSTL(string), '|', gamma, '|'
      WRITE(output_unit, *) '****************************************************************'
    ENDIF

    CALL mpi_barrier(mpi_comm_world, ierr)
#endif

    info(:) = 0._f_real

    ! discriminate between 2D and 3D case, "nd" (alias for "n") is a variable from module "m_psdf"
    IF (cl(3) .gt. 0._f_real) THEN
      ndim = 3
    ELSE
      ndim = 2
    ENDIF

    ! initialise random generator
    c_seed = seed
    CALL set_seed(c_seed)

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    ! set min/max radial wavenumber such that "v1", "v2", "v3" below are in the range ["kmin" "kmax"] when the trigonometric terms
    ! equal 1. This is achieved by taking the largest "kmin*cl" and the smallest "kmax*cl". "cl" is considered in the calculations
    ! because the radial avenumber "k" is divided by "cl" (see Eq.8 of Raess et al., 2019)
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    ! model limits (for each process) along each axis
    min_extent = nc
    max_extent = fc

    ! (global) model limits
    CALL mpi_allreduce(mpi_in_place, min_extent, 3, real_type, mpi_min, mpi_comm_world, ierr)
    CALL mpi_allreduce(mpi_in_place, max_extent, 3, real_type, mpi_max, mpi_comm_world, ierr)

    ! global effective grid size
    span = max_extent - min_extent

    kmin = 2._f_real * pi * maxval(cl(1:ndim) / span(1:ndim))

    kmax = pi / dh * minval(cl(1:ndim))

    ! corner wavenumber for filtering spectrum is controlled by external mesh grid-step
    kc = pi / ds * minval(cl(1:ndim))

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    ! check if the following conditions for a correct wavenumber representation are violated:
    ! 1) kmin < 1/cl
    ! 2) kmax > 1/cl, or ds <= cl / 2 (frenje & juhlin, 2000)
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    IF (ANY(2._f_real * pi / span(1:ndim) .ge. 1._f_real / cl(1:ndim))) info(1) = 1._f_real

    IF (ANY(ds .gt. cl(1:ndim) / 2._f_real)) info(2) = 1._f_real

    ! set scaling factor
    scaling = sigma / SQRT(REAL(nharm, f_real))

    ! number of points where random field must be calculated
    npts(:) = fe(:) - fs(:) + 1

    ! initialise random field
    field(:,:,:) = 0._f_real

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    ! compute random field directly at each grid point. Use uniform distribution in the range [0 1] as majorant function.
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    ! angle about z-axis
    calpha = cos(alpha * pi / 180._f_real)
    salpha = sin(alpha * pi / 180._f_real)

    ! angle about y-axis
    cbeta = cos(beta * pi / 180._f_real)
    sbeta = sin(beta * pi / 180._f_real)

    ! angle about x-axis
    cgamma = cos(gamma * pi / 180._f_real)
    sgamma = sin(gamma * pi / 180._f_real)

    matrix(:, 1) = [calpha*cbeta                       , salpha*cbeta                       , -sbeta      ]
    matrix(:, 2) = [calpha*sbeta*sgamma - salpha*cgamma, salpha*sbeta*sgamma + calpha*cgamma, cbeta*sgamma]
    matrix(:, 3) = [calpha*sbeta*cgamma + salpha*sgamma, salpha*sbeta*cgamma - calpha*sgamma, cbeta*cgamma]

    ! baricenter in original reference frame
    bar = (max_extent - min_extent) / 2._f_real

    ! baricenter after rotation
    obar = MATMUL(matrix, bar)

    ! translation vector to rotate around baricenter in original reference frame
    bar = bar - obar

    ! define power term for vk psdf, "nu" is a variable from module "m_psdf"
    nu = hurst + ndim / 2._f_real

    ! set PSD function
    SELECT CASE (acf)
      CASE(0)
        const = 1._f_real
        fun => fn_vk
      CASE(1)
        IF (ndim .eq. 3) const = 1._f_real / SQRT(pi)                         !< scaling factor to be below unity
        fun => fn_gs
      CASE(2)
        const = 1._f_real
        fun => fn_ud
    END SELECT

    ! update scaling factor, assuming that discrete sigma is always equal to continuous sigma if random field is not normalised.
    ! Note that sigma_discr = sigma_cont * sqrt(intg(kmin,kmax) / intg(0, inf))
    scaling = scaling * SQRT(intg(kmin, kmax)/intg())

    ! loop over harmonics
    ! openacc: "field" is copied in&out, all the others are only copied in. "arg" is created locally on the accelerator
    !$acc data copy(field) copyin(ds, scaling, npts, fs, a, b, v1, v2, v3, matrix, bar)
    DO l = 1, nharm

#ifdef TIMING
      CALL watch_start(tictoc, mpi_comm_self)
#endif

      DO

        ! first set of random numbers
        CALL srng(r)

        ! random variable "k`" must span the range [kmin, kmax]
        u = r(1) * (kmax - kmin) + kmin

        ! evaluate original pdf (require input in double precision)
        d = srm_fn(real(u, f_dble)) * const

        ! take "k" and exit if inequality r < d is verified
        IF (r(2) .lt. d) EXIT

      ENDDO

      ! second set of random numbers
      CALL srng(r)

      ! neglect contribution of current harmonic if wavenumber "u" was larger than corner wavenumber "kc" (always false is "ds" = "dh")
      IF (u .gt. kc) CYCLE

      ! compute azimuth and polar angles
      phi = r(1) * 2._f_real * pi

      IF (ndim .eq. 3) THEN                                             !< spherical coordinates
        theta = ACOS(1._f_real - 2._f_real * r(2))
        v1    = u * SIN(phi) * SIN(theta) / cl(1)
        v2    = u * COS(phi) * SIN(theta) / cl(2)
        v3    = u * COS(theta)            / cl(3)
      ELSE                                                              !< polar coordinates
        v1 = u * COS(phi) / cl(1)
        v2 = u * SIN(phi) / cl(2)
        v3 = 0._f_real
      ENDIF

      ! compute harmonics coefficient "a" and "b"
      ! CALL normdist(r)
      !
      ! a = r(1)
      ! b = r(2)

      ! box-mueller transform
      CALL srng(r)
      a = SQRT(-2._f_real * log(r(1))) * COS(2._f_real * pi * r(2))

      CALL srng(r)
      b = SQRT(-2._f_real * log(r(1))) * COS(2._f_real * pi * r(2))


#ifdef TIMING
      CALL watch_stop(tictoc, mpi_comm_self)

      info(5) = info(5) + REAL(tictoc, f_real)

      CALL watch_start(tictoc, mpi_comm_self)
#endif

      !$acc wait

      ! update selected variables in gpu memory
      !$acc update device(a, b, v1, v2, v3)

      !$acc parallel loop private(x, y, z, arg) collapse(3) async
      DO k = 1, npts(3)
        DO j = 1, npts(2)
          DO i = 1, npts(1)
            x              = (i + fs(1) - 2) * ds
            y              = (j + fs(2) - 2) * ds
            z              = (k + fs(3) - 2) * ds
            ! arg            = v1 * x + v2 * y + v3 * z
            arg            = v1 * ((matrix(1,1)*x + matrix(1,2)*y + matrix(1,3)*z) + bar(1)) +   &
                             v2 * ((matrix(2,1)*x + matrix(2,2)*y + matrix(2,3)*z) + bar(2)) +   &
                             v3 * ((matrix(3,1)*x + matrix(3,2)*y + matrix(3,3)*z) + bar(3))
            field(i, j, k) = field(i, j, k) + a * SIN(arg) + b * COS(arg)
          ENDDO
        ENDDO
      ENDDO
      !$acc end parallel loop

      ! !$acc kernels async
      ! !$acc loop independent private(x, y, z, arg)
      ! DO k = 1, npts(3)
      !  z = (k + fs(3) - 2) * ds
      !  DO j = 1, npts(2)
      !    y = (j + fs(2) - 2) * ds
      !    DO i = 1, npts(1)
      !      x              = (i + fs(1) - 2) * ds
      !      arg            = v1 * x + v2 * y + v3 * z
      !      field(i, j, k) = field(i, j, k) + a * SIN(arg) + b * COS(arg)
      !    ENDDO
      !  ENDDO
      ! ENDDO
      ! !$acc end kernels


#ifdef TIMING
      CALL watch_stop(tictoc, mpi_comm_self)

      info(6) = info(6) + REAL(tictoc, f_real)
#endif

    ENDDO
    ! end loop over harmonics

#ifdef TIMING
    CALL watch_start(tictoc, mpi_comm_self)
#endif

    !$acc wait

    ! normalise random field
    !$acc parallel loop collapse(3)
    DO k = 1, npts(3)
      DO j = 1, npts(2)
        DO i = 1, npts(1)
          field(i, j, k) = field(i, j, k) * scaling
        ENDDO
      ENDDO
    ENDDO
    !$acc end parallel loop

    ! !$acc kernels
    ! !$acc loop independent
    ! DO k = 1, npts(3)
    !  DO j = 1, npts(2)
    !    DO i = 1, npts(1)
    !      field(i, j, k) = field(i, j, k) * scaling
    !    ENDDO
    !  ENDDO
    ! ENDDO
    ! !$acc end kernels

    ! copy "field" back to host and free memory on device
    !$acc end data

#ifdef TIMING
    CALL watch_stop(tictoc, mpi_comm_self)

    info(6) = info(6) + REAL(tictoc, f_real)
#endif

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    ! compute mean and variance of the whole random field. THEN rescale, if desired, discrete std. dev. to its continuous counterpart
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    ALLOCATE(var(0:world_size - 1), mu(0:world_size - 1), npoints(0:world_size - 1))

    ! compute variance and mean of random field for each single process
    var(world_rank) = variance(field)
    mu(world_rank)  = mean(field)

    ! number of grid points per process
    npoints(world_rank) = npts(1) * npts(2) * npts(3)

    ! share results
    CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, var, 1, real_type, mpi_comm_world, ierr)
    CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, mu, 1, real_type, mpi_comm_world, ierr)
    CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, npoints, 1, mpi_integer, mpi_comm_world, ierr)

    ! return total variance in "info(3)" and mean in "info(4)"
    CALL parallel_variance(var, mu, npoints, info(3), info(4))

    ! return standard deviation
    info(3) = SQRT(info(3))

    IF (rescale .eq. 1) THEN

      scaling = sigma / info(3)

      DO k = 1, npts(3)
        DO j = 1, npts(2)
          DO i = 1, npts(1)
            field(i, j, k) = field(i, j, k) * scaling - info(4)
          ENDDO
        ENDDO
      ENDDO

      info(3) = sigma
      info(4) = 0._f_real

    ENDIF

    DEALLOCATE(var, mu, npoints)

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    ! taper/mute random field at desired locations.
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    pt(:) = 0._f_real

    DO i = 1, SIZE(poi, 2)
      pt(1:ndim) = poi(1:ndim, i)
      CALL tapering(ds, fs, fe, field, pt, mute, taper)
    ENDDO

    CALL mpi_barrier(mpi_comm_world, ierr)

  END SUBROUTINE scarf3d_structured_srm

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  !=================================================================================================================================
  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  REAL(f_real) FUNCTION intg(k0, k1)

    ! Purpose:
    ! To evaluate a definite integrate of the PSD function. The actual integration algorithm is chosen depending on whether the
    ! integral is improper (optional args not specified) or not. In the former case, integration endpoints are 0 and +inf.
    !
    ! Revisions:
    !     Date                    Description of change
    !     ====                    =====================
    !   08/09/20                  original version
    !

    REAL(f_real),                    OPTIONAL, INTENT(IN) :: k0, k1
    REAL(f_dble)                                          :: kmin, kmax, solv
    INTEGER(f_int)                                        :: neval, ierr, last
    INTEGER(f_int),                  PARAMETER            :: limlst = 50, limit = 500, leniw = 2 * limit + limlst
    INTEGER(f_int),                  PARAMETER            :: maxp1 = 21, lenw = leniw * 2 + maxp1 * 25
    INTEGER(f_int), DIMENSION(leniw)                      :: iwork
    REAL(f_dble)                                          :: abserr
    REAL(f_dble),                    PARAMETER            :: epsabs = 1.e-08_f_dble, epsrel = 1.e-05_f_dble
    REAL(f_dble),   DIMENSION(lenw)                       :: work

    !-------------------------------------------------------------------------------------------------------------------------------

    IF (PRESENT(k0) .and. PRESENT(k1)) THEN
      kmin = REAL(k0, f_dble)
      kmax = REAL(k1, f_dble)
      CALL dqng(srm_fn, kmin, kmax, epsabs, epsrel, solv, abserr, neval, ierr)
    ELSE
      CALL dqagi(srm_fn, 0._f_dble, 1, epsabs, epsrel, solv, abserr, neval, ierr, limit, lenw, last, iwork, work)
    ENDIF

    intg = REAL(solv, f_real)

  END FUNCTION intg

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  !=================================================================================================================================
  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  REAL(f_dble) FUNCTION srm_fn(x)

    ! Purpose:
    ! To compute 2D and 3D Gaussian PSDF for the SRM. Input and result must be in double precision as requested by numerical
    ! integration routines.
    !
    ! Revisions:
    !     Date                    Description of change
    !     ====                    =====================
    !   04/05/20                  original version
    !

    REAL(f_dble), INTENT(IN) :: x

    !-----------------------------------------------------------------------------------------------------------------------------

    srm_fn = x * fun(REAL(x, f_real))

    IF (ndim .eq. 3) srm_fn = srm_fn * x

  END FUNCTION srm_fn

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  !=================================================================================================================================
  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

END MODULE m_scarflib_srm
