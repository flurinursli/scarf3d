MODULE m_scarflib_srm

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
  !

  USE, NON_INTRINSIC :: m_scarflib_common

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

  ! number of harmonics in spectral summation
  INTEGER(f_int), PARAMETER :: nharm = 10000

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  !=================================================================================================================================
  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  SUBROUTINE scarf3d_unstructured_srm(nc, fc, ds, x, y, z, dh, acf, cl, sigma, hurst, seed, poi, mute, taper, rescale, field, info)

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
    REAL(f_real),                DIMENSION(:),   INTENT(OUT) :: field        !< random field at x,y,z location
    REAL(f_real),                DIMENSION(6),   INTENT(OUT) :: info         !< errors flags and TIMING for performance analysis
    INTEGER(f_int)                                           :: i, l
    INTEGER(f_int)                                           :: npts
    INTEGER(f_int)                                           :: ierr
    INTEGER(c_int)                                           :: c_seed
    INTEGER(f_int), ALLOCATABLE, DIMENSION(:)                :: npoints
    REAL(f_real)                                             :: scaling
    REAL(f_real)                                             :: kmax, kmin, kc
    REAL(f_real)                                             :: k, d
    REAL(f_real)                                             :: phi, theta
    REAL(f_real)                                             :: v1, v2, v3
    REAL(f_real)                                             :: a, b, arg
    REAL(f_dble)                                             :: tictoc
    REAL(c_real),                DIMENSION(2)                :: r
    REAL(f_real),                DIMENSION(3)                :: span
    REAL(f_real),                DIMENSION(3)                :: min_extent, max_extent
    REAL(f_real),   ALLOCATABLE, DIMENSION(:)                :: mu, var

    !-------------------------------------------------------------------------------------------------------------------------------

    CALL mpi_comm_SIZE(mpi_comm_world, world_size, ierr)
    CALL mpi_comm_rank(mpi_comm_world, world_rank, ierr)

    IF (world_rank == 0) THEN
      print*, 'input params @ unstructured-spec'
      print*, 'nc ', nc
      print*, 'fc ', fc
      print*, 'ds ', ds
      print*, 'x ', SIZE(x)
      print*, 'z ', SIZE(z)
      print*, 'dh ', dh
      print*, 'acf ', acf
      print*, 'cl ', cl
      print*, 'sigma ', sigma
      print*, 'hurst ', hurst
      print*, 'seed ', seed
      print*, 'poi ', SIZE(poi)
      print*, 'mute ', mute
      print*, 'taper ', taper
    ENDIF
    CALL mpi_barrier(mpi_comm_world, ierr)

    info(:) = 0._f_real

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

    kmin = 2._f_real * pi * maxval(cl / span)

    kmax = pi / dh * minval(cl)

    ! corner wavenumber for filtering spectrum is controlled by external mesh grid-step
    kc = pi / ds * minval(cl)

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    ! check if the following conditions for a correct wavenumber representation are violated:
    ! 1) kmin < 1/cl
    ! 2) kmax > 1/cl, or ds <= cl / 2 (frenje & juhlin, 2000)
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    IF (ANY(2._f_real * pi / span .ge. 1._f_real / cl)) info(1) = 1._f_real

    IF (ANY(ds .gt. cl / 2._f_real)) info(2) = 1._f_real

    ! set scaling factor
    scaling = sigma / SQRT(REAL(nharm, f_real))

    ! number of points where random field must be calculated
    npts = SIZE(x)

    ! initialise random field
    field(:) = 0._f_real

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    ! compute random field directly at each grid point. USE uniform distribution in the range [0 1] as majorant function.
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    ! loop over harmonics
    ! openacc: "field" is copied in&out, all the others are only copied in. "arg" is created locally on the accelerator
    !$acc data copy(field) copyin(x, y, z, scaling, npts, a, b, v1, v2, v3)
    DO l = 1, nharm

#ifdef TIMING
      CALL watch_start(tictoc, mpi_comm_self)
#endif

      DO

        ! first set of random numbers
        CALL srng(r)

        ! random variable "k`" must span the range [kmin, kmax]
        k = r(1) * (kmax - kmin) + kmin

        ! evaluate original pdf
        IF (acf .eq. 0) THEN
          d = k**2 / (1._f_real + k**2)**(1.5_f_real + hurst)
        ELSE
          d = k**2 * EXP(-0.25_f_real * k**2) / SQRT(pi)
        ENDIF

        ! take "k" and exit if inequality r < d is verified
        IF (r(2) .lt. d) EXIT

      ENDDO

      ! second set of random numbers
      CALL srng(r)

      ! neglect contribution of current harmonic if wavenumber "u" was larger than corner wavenumber "kc" (always false is "ds" = "dh")
      IF (k .gt. kc) CYCLE

      ! compute azimuth and polar angles
      phi   = r(1) * 2._f_real * pi
      theta = ACOS(1._f_real - 2._f_real * r(2))

      v1 = k * SIN(phi) * SIN(theta) / cl(1)
      v2 = k * COS(phi) * SIN(theta) / cl(2)
      v3 = k * COS(theta)            / cl(3)

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
        arg = v1 * x(i) + v2 * y(i) + v3 * z(i)
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
    !$acc END kernels

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

    DO i = 1, SIZE(poi, 2)
      CALL tapering(x, y, z, field, poi(:, i), mute, taper)
    ENDDO

    CALL mpi_barrier(mpi_comm_world, ierr)

  END SUBROUTINE scarf3d_unstructured_srm

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  !=================================================================================================================================
  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  SUBROUTINE scarf3d_structured_srm(nc, fc, ds, fs, fe, dh, acf, cl, sigma, hurst, seed, poi, mute, taper, rescale, field, info)

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
    REAL(f_real),                DIMENSION(:,:,:), INTENT(OUT) :: field           !< random field
    REAL(f_real),                DIMENSION(6),     INTENT(OUT) :: info            !< errors flags and TIMING for performance analysis
    INTEGER(f_int)                                             :: i, j, k, l
    INTEGER(f_int)                                             :: ierr
    INTEGER(c_int)                                             :: c_seed
    INTEGER(f_int),              DIMENSION(3)                  :: npts
    INTEGER(f_int), ALLOCATABLE, DIMENSION(:)                  :: npoints
    REAL(f_real)                                               :: scaling
    REAL(f_real)                                               :: kmax, kmin, kc
    REAL(f_real)                                               :: u, d
    REAL(f_real)                                               :: phi, theta
    REAL(f_real)                                               :: v1, v2, v3
    REAL(f_real)                                               :: a, b, arg
    REAL(f_real)                                               :: x, y, z
    REAL(f_real)                                               :: calpha, salpha, cbeta, sbeta, cgamma, sgamma
    REAL(f_dble)                                               :: tictoc
    REAL(c_real),                DIMENSION(2)                  :: r
    REAL(f_real),                DIMENSION(3)                  :: span
    REAL(f_real),                DIMENSION(3)                  :: min_extent, max_extent
    REAL(f_real),                DIMENSION(3)                  :: bar, obar
    REAL(f_real),   ALLOCATABLE, DIMENSION(:)                  :: mu, var
    REAL(f_real),                DIMENSION(3,3)                :: matrix

    !-------------------------------------------------------------------------------------------------------------------------------

    CALL mpi_comm_SIZE(mpi_comm_world, world_size, ierr)
    CALL mpi_comm_rank(mpi_comm_world, world_rank, ierr)

    IF (world_rank == 0) THEN
      print*, 'input params @ structured-spec'
      print*, 'nc ', nc
      print*, 'fc ', fc
      print*, 'ds ', ds
      print*, 'fs ', fs
      print*, 'fe ', fe
      print*, 'dh ', dh
      print*, 'acf ', acf
      print*, 'cl ', cl
      print*, 'sigma ', sigma
      print*, 'hurst ', hurst
      print*, 'seed ', seed
      print*, 'poi ', SIZE(poi)
      print*, 'mute ', mute
      print*, 'taper ', taper
    ENDIF
    CALL mpi_barrier(mpi_comm_world, ierr)


    info(:) = 0._f_real

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

    kmin = 2._f_real * pi * maxval(cl / span)

    kmax = pi / dh * minval(cl)

    ! corner wavenumber for filtering spectrum is controlled by external mesh grid-step
    kc = pi / ds * minval(cl)

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    ! check if the following conditions for a correct wavenumber representation are violated:
    ! 1) kmin < 1/cl
    ! 2) kmax > 1/cl, or ds <= cl / 2 (frenje & juhlin, 2000)
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    IF (ANY(2._f_real * pi / span .ge. 1._f_real / cl)) info(1) = 1._f_real

    IF (ANY(ds .gt. cl / 2._f_real)) info(2) = 1._f_real

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
    calpha = cos(20._f_real * pi / 180._f_real)
    salpha = sin(20._f_real * pi / 180._f_real)

    ! angle about y-axis
    cbeta = cos(0._f_real * pi / 180._f_real)
    sbeta = sin(0._f_real * pi / 180._f_real)

    ! angle about x-axis
    cgamma = cos(0._f_real * pi / 180._f_real)
    sgamma = sin(0._f_real * pi / 180._f_real)

    matrix(:, 1) = [calpha*cbeta, salpha*cbeta, -sbeta]
    matrix(:, 2) = [calpha*sbeta*sgamma - salpha*cgamma, salpha*sbeta*sgamma + calpha*cgamma, cbeta*sgamma]
    matrix(:, 3) = [calpha*sbeta*cgamma + salpha*sgamma, salpha*sbeta*cgamma - calpha*sgamma, cbeta*cgamma]

    ! baricenter in original reference frame
    bar = (max_extent - min_extent) / 2._f_real

    ! baricenter after rotation
    obar = MATMUL(matrix, bar)

    ! translation vector to rotate around baricenter in original reference frame
    bar = bar - obar


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

        ! evaluate original pdf
        IF (acf .eq. 0) THEN
          d = u**2 / (1._f_real + u**2)**(1.5_f_real + hurst)
        ELSE
          d = u**2 * EXP(-0.25_f_real * u**2) / SQRT(pi)
        ENDIF

        ! take "k" and exit if inequality r < d is verified
        IF (r(2) .lt. d) EXIT

      ENDDO

      ! second set of random numbers
      CALL srng(r)

      ! neglect contribution of current harmonic if wavenumber "u" was larger than corner wavenumber "kc" (always false is "ds" = "dh")
      IF (u .gt. kc) CYCLE

      ! compute azimuth and polar angles
      phi   = r(1) * 2._f_real * pi
      theta = ACOS(1._f_real - 2._f_real * r(2))

      v1 = u * SIN(phi) * SIN(theta) / cl(1)
      v2 = u * COS(phi) * SIN(theta) / cl(2)
      v3 = u * COS(theta)            / cl(3)

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
            !arg            = v1 * (x * calpha + z * salpha) + v2 * y + v3 * (-x * salpha + z * calpha)
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

    DO i = 1, SIZE(poi, 2)
      CALL tapering(ds, fs, fe, field, poi(:, i), mute, taper)
    ENDDO

    CALL mpi_barrier(mpi_comm_world, ierr)

  END SUBROUTINE scarf3d_structured_srm

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  !=================================================================================================================================
  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --


END MODULE m_scarflib_srm
