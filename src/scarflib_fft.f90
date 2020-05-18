MODULE m_scarflib_fim

  ! Purpose:
  !   To compute a random field according to the Fourier Integral Method (FIM). Subprograms rely on the FFTW (version 3) library to
  !   compute the fast Fourier transform and on the TRNG (version 4) library to generate random numbers.
  !   Subroutines 'scarf3d_unstructured_fim' and 'scarf3d_structured_fim' are virtually indentical, with only few differences sparse
  !   at different points in the code (this inhibits somehow the use of 'include' statements to handle a single code)
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   04/05/20                  original version
  !   11/05/20                  updated macro for double-precision
  !

  ! mpif90 -o3 scarflib_fft2.f90 scarflib_spec.f90 scarflib.f90 driver.f90 *.f -i/home/walter/backedup/software/FFTW-3.3.8/INCLUDE
  ! -l/home/walter/backedup/software/FFTW-3.3.8/lib -lfftw3 -fcheck=ALL -fbacktrace -fopenacc

  ! import mpi library and intrinsic modules, data types and subprograms
  USE, NON_INTRINSIC :: m_scarflib_common

  IMPLICIT none

  ! include FFTW interface
  INCLUDE 'fftw3.f03'

  PRIVATE

  ! make only generic interface accessible
  PUBLIC :: scarf3d_fim

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTERFACE scarf3d_fim
    MODULE PROCEDURE scarf3d_unstructured_fim, scarf3d_structured_fim
  END INTERFACE

  INTERFACE grid_mapping
    MODULE PROCEDURE grid_mapping_unstruct, grid_mapping_struct
  END INTERFACE

  INTERFACE order_points
    MODULE PROCEDURE order_points_unstruct, order_points_struct
  END INTERFACE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! define a procedure pointer to compute spectrum according to selected ACF
  PROCEDURE(vk_psdf), POINTER :: fun

  ! interface to C++ functions/subroutines
  INTERFACE

    FUNCTION prng(seed, ls, le, npts) BIND(c, name="prng")
      USE, INTRINSIC :: iso_c_binding
      IMPLICIT none
      INTEGER(c_int),                          VALUE :: seed
      INTEGER(c_int), DIMENSION(3), INTENT(IN)       :: ls, le
      INTEGER(c_int), DIMENSION(3), INTENT(IN)       :: npts
      TYPE(c_ptr)                                    :: prng
    END FUNCTION prng

    SUBROUTINE free_mem(ptr) BIND(c, name="free_mem")
      USE, INTRINSIC :: iso_c_binding
      IMPLICIT none
      TYPE(c_ptr), VALUE, INTENT(IN) :: ptr
    END SUBROUTINE free_mem

  END INTERFACE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

  ! work with 3D cartesian grids
  INTEGER(f_int),                            PARAMETER  :: ndims = 3

  ! mask to determine if i-th process has points falling inside domain of calling process (0=no, 1=yes)
  INTEGER(f_int), ALLOCATABLE, DIMENSION(:)             :: points_world2rank

  ! mask to determine if calling process has points falling inside domain of i-th process (0=no, 1=yes)
  INTEGER(f_int), ALLOCATABLE, DIMENSION(:)             :: points_rank2world

  ! allow grid re-ordering
  LOGICAL,                                   PARAMETER  :: reorder = .true.

  ! set cartesian grid periodicity
  LOGICAL,                     DIMENSION(3), PARAMETER  :: isperiodic = [.false., .false., .false.]

  ! set a very large number
  REAL(f_real),                              PARAMETER  :: big = HUGE(1._f_real)

  ! absolute position of model's first point
  REAL(f_real),                DIMENSION(3)             :: off_axis

  ! pointer for fourier transform
  COMPLEX(c_cplx),             DIMENSION(:), POINTER    :: cdum => NULL()

  ! pointer for fourier transform
  REAL(c_real),                DIMENSION(:), POINTER    :: rdum => NULL()

  ! pointer to FFTW plan
  TYPE(c_ptr)                                           :: pc

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE scarf3d_unstructured_fim(nc, fc, ds, x, y, z, dh, acf, cl, sigma, hurst, seed, poi, mute, taper, rescale, pad, &
                                        field, info)

      ! Purpose:
      ! To compute a random field characterised by a selected ACF on unstructured meshes. For a given ACF, the actual distribution of
      ! heterogeneity depends on number points of the internal (FFT) grid and the seed number. The former may be influenced by optional
      ! padding (to avoid wrap-around correlation effects). The internal grid is determined by the near-corner 'nc', far-corner 'fc'
      ! and the grid-step 'dh'. The spectrum will be filtered according to 'ds' (grid-step external grid) if the latter is not equal
      ! to 'dh'.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_real),                 DIMENSION(3),    INTENT(IN)  :: nc, fc       !< min/max extent of input grid (m or km)
      REAL(f_real),                                  INTENT(IN)  :: ds           !< grid-step external grid, control max resolvable wavenumber
      REAL(f_real),                 DIMENSION(:),    INTENT(IN)  :: x, y, z      !< position of points along x, y, z
      REAL(f_real),                                  INTENT(IN)  :: dh           !< grid-step of internal (FFT) grid
      INTEGER(f_int),                                INTENT(IN)  :: acf          !< autocorrelation function: 0=vk, 1=gauss
      REAL(f_real),                 DIMENSION(3),    INTENT(IN)  :: cl           !< correlation length
      REAL(f_real),                                  INTENT(IN)  :: sigma        !< standard deviation
      REAL(f_real),                                  INTENT(IN)  :: hurst        !< hurst exponent
      INTEGER(f_int),                                INTENT(IN)  :: seed         !< seed number
      REAL(f_real),                 DIMENSION(:,:),  INTENT(IN)  :: poi          !< location of point(s)-of-interest where muting/tapering is applied
      REAL(f_real),                                  INTENT(IN)  :: mute         !< radius for muting
      REAL(f_real),                                  INTENT(IN)  :: taper        !< radius for tapering
      INTEGER(f_int),                                INTENT(IN)  :: rescale      !< flag for rescaling random field to desired sigma
      INTEGER(f_int),                                INTENT(IN)  :: pad          !< flag for handle fft periodicity (grid padding)
      REAL(f_real),                 DIMENSION(:),    INTENT(OUT) :: field        !< random field at x,y,z location
      REAL(f_real),                 DIMENSION(8),    INTENT(OUT) :: info         !< errors flags and timing for performance analysis
      COMPLEX(f_real), ALLOCATABLE, DIMENSION(:,:,:)             :: spec
      INTEGER(f_int)                                             :: ierr
      INTEGER(f_int)                                             :: i, j, k
      INTEGER(f_int)                                             :: cartopo
      INTEGER(f_int)                                             :: offset
      INTEGER(f_int),               DIMENSION(3)                 :: m
      INTEGER(f_int),               DIMENSION(3)                 :: ls, le
      INTEGER(f_int),               DIMENSION(3)                 :: coords
      INTEGER(f_int),               DIMENSION(3)                 :: dims
      INTEGER(f_int),               DIMENSION(3)                 :: ps, pe
      INTEGER(f_int),  ALLOCATABLE, DIMENSION(:)                 :: npoints
      LOGICAL,                      DIMENSION(3)                 :: bool
      REAL(f_real)                                               :: scaling
      REAL(f_dble)                                               :: tictoc
      REAL(f_real),                 DIMENSION(2)                 :: et
      REAL(f_real),                 DIMENSION(3)                 :: min_extent
      REAL(f_real),                 DIMENSION(3)                 :: max_extent
      REAL(f_real),                 DIMENSION(3)                 :: pt
      REAL(f_real),    ALLOCATABLE, DIMENSION(:)                 :: var, mu
      REAL(f_real),    ALLOCATABLE, DIMENSION(:,:,:)             :: delta, buffer

      !-----------------------------------------------------------------------------------------------------------------------------

      ! get rank number
      CALL mpi_comm_rank(mpi_comm_world, world_rank, ierr)

      ! get available mpi processes
      CALL mpi_comm_size(mpi_comm_world, world_size, ierr)

      ! initialise variables
      info(:) = 0._f_real

      IF (world_rank == 0) THEN
        print*, 'input params @ unstructured-fft'
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
        print*, 'rescale ', rescale
        print*, 'pad ', pad
      ENDIF
      CALL mpi_barrier(mpi_comm_world, ierr)


      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! determine internal (FFT) grid size
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! model limits (for each process) along each axis
      min_extent = nc
      max_extent = fc

      ! (global) model limits
      CALL mpi_allreduce(mpi_in_place, min_extent, 3, real_type, mpi_min, mpi_comm_world, ierr)
      CALL mpi_allreduce(mpi_in_place, max_extent, 3, real_type, mpi_max, mpi_comm_world, ierr)

      ! cycle over the three main directions to determine the necessary model size (in number of points) and the absolute position
      ! of the first point to counteract fft periodicity
      DO i = 1, 3

        ! minimum number of points such that random field covers the whole domain. In practice, we define the random field on a grid
        ! slightly larger than necessary (half grid-step) in each direction to avoid side-effects due to interpolation (useful when
        ! the extra extension to handle fft periodicity is not desired)
        npts(i) = NINT( (max_extent(i) + dh * 0.5_f_real - min_extent(i) + dh * 0.5_f_real) / dh) + 1

        IF (world_rank == 0) print*, 'npts ', npts(i), ' - ', min_extent, ' ', max_extent

        ! padding depends on correlation length
        offset = NINT(cl(i) / dh)

        ! do not pad model unless desired
        IF (pad .ne. 1) offset = 0

        ! new model size including padding
        npts(i) = npts(i) + 2 * offset

        ! make sure we have even number of points (allow to handle Nyquist frequency explicitly)
        IF (MOD(npts(i), 2) .ne. 0) npts(i) = npts(i) + 1

        ! absolute position of first point (it could be negative)
        off_axis(i) = min_extent(i) - (offset + 0.5_f_real) * dh

      ENDDO

      ! check if model is large enough to catch the lower part of the spectrum (low wavenumbers)...
      IF (ANY(npts .le. NINT(2._f_real * pi * cl / dh))) info(1) = 1._f_real

      ! ...and if internal grid-step is small enough to catch the upper part of the spectrum (high wavenumbers)
      IF (ANY(dh .gt. cl / 2._f_real)) info(2) = 1._f_real

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! create regular mesh (for FFT) and associated cartesian topology: random field will be computed on this mesh and then interpolated
      ! on the external one
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! "gs"/"ge" are specified in "m_scarflib_common" and store first/last global index along each direction for all processes
      ALLOCATE(gs(3, 0:world_size-1), ge(3, 0:world_size-1))

      ! return processors grid
      CALL best_config(dims)

      ! create topology
      CALL mpi_cart_create(mpi_comm_world, ndims, dims, isperiodic, reorder, cartopo, ierr)

      ! return process coordinates in current topology
      CALL mpi_cart_coords(cartopo, world_rank, ndims, coords, ierr)

      ! return first/last index ("ls"/"le") along each direction for the calling process. Note: first point has "ls = 1".
      CALL coords2index(npts, dims, coords, ls, le)

      gs(:, world_rank) = ls
      ge(:, world_rank) = le

      ! make all processes aware of global indices
      CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, gs, 3, mpi_integer, mpi_comm_world, ierr)
      CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, ge, 3, mpi_integer, mpi_comm_world, ierr)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! determine if calling process has any external grid point belonging to other processes and which processes have external grid
      ! points belonging to calling process. These info are needed to define communicator groups when interpolating data.
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ALLOCATE(points_world2rank(0:world_size - 1), points_rank2world(0:world_size - 1))

      ! min/max internal grid index containing input points
      ps(1) = FLOOR((MINVAL(x, dim = 1) - off_axis(1)) / dh) + 1
      pe(1) = FLOOR((MAXVAL(x, dim = 1) - off_axis(1)) / dh) + 1
      ps(2) = FLOOR((MINVAL(y, dim = 1) - off_axis(2)) / dh) + 1
      pe(2) = FLOOR((MAXVAL(y, dim = 1) - off_axis(2)) / dh) + 1
      ps(3) = FLOOR((MINVAL(z, dim = 1) - off_axis(3)) / dh) + 1
      pe(3) = FLOOR((MAXVAL(z, dim = 1) - off_axis(3)) / dh) + 1

      points_rank2world = 0

      ! determine if calling process has at least one point falling inside domain of "i-th" process
      DO i = 0, world_size - 1
        DO j = 1, 3
          bool(j) = ( (ps(j) .lt. ge(j, i)) .and. (ps(j) .ge. (gs(j, i) - 1)) ) .or.     &
                    ( (pe(j) .lt. ge(j, i)) .and. (pe(j) .ge. (gs(j, i) - 1)) ) .or.     &
                    ( (pe(j) .ge. ge(j, i)) .and. (ps(j) .lt. (gs(j, i) - 1)) )
        ENDDO
        IF (ALL(bool .eqv. .true.)) points_rank2world(i) = 1
      ENDDO

      ! determine if process "i" has at least one point falling inside domain of calling process
      CALL mpi_alltoall(points_rank2world, 1, mpi_integer, points_world2rank, 1, mpi_integer, mpi_comm_world, ierr)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! compute power spectrum and inverse FFT
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! number of points for calling process
      DO i = 1, 3
        m(i) = le(i) - ls(i) + 1
      ENDDO

      ! allocate resources to hold spectrum
      ALLOCATE(spec(m(1), m(2), m(3)))

      ! compute spectrum and apply hermitian symmetry
      CALL compute_spectrum(ds, ls, le, dh, acf, cl, sigma, hurst, seed, spec, et)

#ifdef TIMING
      info(5:6) = et

      ! start timer
      CALL watch_start(tictoc)
#endif

      ! (inverse) transform along each direction
      CALL transform_along_z(spec)
      CALL transform_along_y(spec)
      CALL transform_along_x(spec)

#ifdef TIMING
      CALL watch_stop(tictoc)

      info(7) = tictoc
#endif

      ! define scaling factor
      scaling = 1._f_real / SQRT(REAL(npts(1), f_real) * REAL(npts(2), f_real) * REAL(npts(3), f_real) * dh**3)

      ALLOCATE(buffer(m(1), m(2), m(3)))

      ! scale field
      DO k = 1, m(3)
        DO j = 1, m(2)
          DO i = 1, m(1)
            buffer(i, j, k) = REAL(spec(i, j, k), f_real) * scaling
          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE(spec)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! compute mean and variance of the whole random field. Then rescale, if desired, discrete std. dev. to its continuous counterpart
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ALLOCATE(var(0:world_size - 1), mu(0:world_size - 1), npoints(0:world_size - 1))

      ! variance and mean for each single process
      var(world_rank) = variance(buffer)
      mu(world_rank)  = mean(buffer)

      ! number of grid points per process
      npoints(world_rank) = PRODUCT(m)

      ! share results
      CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, var, 1, real_type, mpi_comm_world, ierr)
      CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, mu, 1, real_type, mpi_comm_world, ierr)
      CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, npoints, 1, mpi_integer, mpi_comm_world, ierr)

      ! return total variance in "info(3)" and mean in "info(4)"
      CALL parallel_variance(var, mu, npoints, info(3), info(4))

      ! we want standard deviation
      info(3) = SQRT(info(3))

      IF (rescale .eq. 1) THEN

        scaling = sigma / info(3)

        DO k = 1, m(3)
          DO j = 1, m(2)
            DO i = 1, m(1)
              buffer(i, j, k) = buffer(i, j, k) * scaling - info(4)
            ENDDO
          ENDDO
        ENDDO

        info(3) = sigma
        info(4) = 0._f_real

      ENDIF

      DEALLOCATE(var, mu, npoints)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! taper/mute random field at desired locations.
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      DO i = 1, SIZE(poi, 2)
        pt(:) = poi(:, i) - off_axis(:)                         !< update "poi" according to coordinates of internal grid
        CALL tapering(dh, ls, le, buffer, pt, mute, taper)
      ENDDO

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! each process exchange data with neighbors to define halos to be used during interpolation
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! copy random field to array with halo
      ALLOCATE(delta(0:m(1), 0:m(2), 0:m(3)))

      delta = 0._f_real

      ! copy field
      DO k = 1, m(3)
        DO j = 1, m(2)
          DO i = 1, m(1)
            delta(i, j, k) = buffer(i, j, k)
          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE(buffer)

#ifdef TIMING
      CALL watch_start(tictoc)
#endif

      ! exchange halo
      CALL exchange_halo(cartopo, delta)

      ! interpolate random field values at desired output locations
      CALL grid_mapping(dh, delta, x, y, z, field)

#ifdef TIMING
      CALL watch_stop(tictoc)

      info(8) = tictoc
#endif

      ! release memory
      DEALLOCATE(delta)
      DEALLOCATE(gs, ge)
      DEALLOCATE(points_rank2world, points_world2rank)

      CALL mpi_comm_free(cartopo, ierr)

      CALL mpi_barrier(mpi_comm_world, ierr)

    END SUBROUTINE scarf3d_unstructured_fim

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE scarf3d_structured_fim(nc, fc, ds, fs, fe, dh, acf, cl, sigma, hurst, seed, poi, mute, taper, rescale, pad, field, &
                                      info)

      ! Purpose:
      ! To compute a random field characterised by a selected ACF on unstructured meshes. For a given ACF, the actual distribution of
      ! heterogeneity depends on number points of the internal (FFT) grid and the seed number. The former may be influenced by optional
      ! padding (to avoid wrap-around correlation effects). The internal grid is determined by the near-corner 'nc', far-corner 'fc'
      ! and the grid-step 'dh'. The spectrum will be filtered according to 'ds' (grid-step external grid) if the latter is not equal
      ! to 'dh'.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_real),                 DIMENSION(3),     INTENT(IN)  :: nc, fc          !< min/max extent of input grid (m or km)
      REAL(f_real),                                   INTENT(IN)  :: ds              !< grid-step external grid, control max resolvable wavenumber
      INTEGER(f_int),               DIMENSION(3),     INTENT(IN)  :: fs, fe          !< first/last external structured grid index
      REAL(f_real),                                   INTENT(IN)  :: dh              !< grid-step of internal (FFT) grid
      INTEGER(f_int),                                 INTENT(IN)  :: acf             !< autocorrelation function: 0=vk, 1=gauss
      REAL(f_real),                 DIMENSION(3),     INTENT(IN)  :: cl              !< correlation length
      REAL(f_real),                                   INTENT(IN)  :: sigma           !< standard deviation
      REAL(f_real),                                   INTENT(IN)  :: hurst           !< hurst exponent
      INTEGER(f_int),                                 INTENT(IN)  :: seed            !< seed number
      REAL(f_real),                 DIMENSION(:,:),   INTENT(IN)  :: poi             !< location of point(s)-of-interest where muting/tapering is applied
      REAL(f_real),                                   INTENT(IN)  :: mute            !< radius for mutin
      REAL(f_real),                                   INTENT(IN)  :: taper           !< radius for tapering
      INTEGER(f_int),                                 INTENT(IN)  :: rescale         !< flag for rescaling random field to desired sigma
      INTEGER(f_int),                                 INTENT(IN)  :: pad             !< flag for handle fft periodicity (grid padding)
      REAL(f_real),                 DIMENSION(:,:,:), INTENT(OUT) :: field           !< random field
      REAL(f_real),                 DIMENSION(8),     INTENT(OUT) :: info            !< errors flags and timing for performance analysis
      COMPLEX(f_real), ALLOCATABLE, DIMENSION(:,:,:)              :: spec
      INTEGER(f_int)                                              :: ierr
      INTEGER(f_int)                                              :: i, j, k
      INTEGER(f_int)                                              :: cartopo
      INTEGER(f_int)                                              :: offset
      INTEGER(f_int),               DIMENSION(3)                  :: m
      INTEGER(f_int),               DIMENSION(3)                  :: ls, le
      INTEGER(f_int),               DIMENSION(3)                  :: coords
      INTEGER(f_int),               DIMENSION(3)                  :: dims
      INTEGER(f_int),               DIMENSION(3)                  :: ps, pe
      INTEGER(f_int),  ALLOCATABLE, DIMENSION(:)                  :: npoints
      LOGICAL,                      DIMENSION(3)                  :: bool
      REAL(f_real)                                                :: scaling
      REAL(f_dble)                                                :: tictoc
      REAL(f_real),                 DIMENSION(2)                  :: et
      REAL(f_real),                 DIMENSION(3)                  :: min_extent
      REAL(f_real),                 DIMENSION(3)                  :: max_extent
      REAL(f_real),                 DIMENSION(3)                  :: pt
      REAL(f_real),    ALLOCATABLE, DIMENSION(:)                  :: var, mu
      REAL(f_real),    ALLOCATABLE, DIMENSION(:,:,:)              :: delta, buffer

      !-----------------------------------------------------------------------------------------------------------------------------

      ! get rank number
      CALL mpi_comm_rank(mpi_comm_world, world_rank, ierr)

      ! get available mpi processes
      CALL mpi_comm_size(mpi_comm_world, world_size, ierr)

      ! initialise variables
      info(:) = 0._f_real

      IF (world_rank == 0) THEN
        print*, 'input params @ structured-fft'
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
        print*, 'poi ', SIZE(poi), SIZE(poi, 1), SIZE(poi, 2)
        print*, 'mute ', mute
        print*, 'taper ', taper
        print*, 'rescale ', rescale
        print*, 'pad ', pad
      ENDIF
      CALL mpi_barrier(mpi_comm_world, ierr)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! determine internal (FFT) grid size
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! model limits (for each process) along each axis
      min_extent = nc
      max_extent = fc

      ! (global) model limits
      CALL mpi_allreduce(mpi_in_place, min_extent, 3, real_type, mpi_min, mpi_comm_world, ierr)
      CALL mpi_allreduce(mpi_in_place, max_extent, 3, real_type, mpi_max, mpi_comm_world, ierr)

      ! CYCLE over the three main directions to determine the necessary model size (in number of points) and the absolute position
      ! of the first point to counteract fft periodicity
      DO i = 1, 3

        ! minimum number of points such that random field covers the whole domain. In practice, we define the random field on a grid
        ! slightly larger than necessary (half grid-step) in each direction to avoid side-effects due to interpolation (useful when
        ! the extra extension to handle fft periodicity is not desired)
        npts(i) = NINT( (max_extent(i) + dh * 0.5_f_real - min_extent(i) + dh * 0.5_f_real) / dh) + 1

        IF (world_rank == 0) print*, 'component', i, ' npts ', npts(i), ' - ', min_extent, ' ', max_extent

        ! padding depends on correlation length
        offset = NINT(cl(i) / dh)

        ! do not pad model unless desired
        IF (pad .ne. 1) offset = 0

        ! new model size including padding
        npts(i) = npts(i) + 2 * offset

        ! make sure we have even number of points (allow to handle Nyquist frequency explicitly)
        IF (MOD(npts(i), 2) .ne. 0) npts(i) = npts(i) + 1

        ! absolute position of first point (it could be negative)
        off_axis(i) = min_extent(i) - (offset + 0.5_f_real) * dh

      ENDDO

      ! check if model is large enough to catch the lower part of the spectrum (low wavenumbers)...
      IF (ANY(npts .le. NINT(2._f_real * pi * cl / dh))) info(1) = 1._f_real

      ! ...and if internal grid-step is small enough to catch the upper part of the spectrum (high wavenumbers)
      IF (ANY(dh .gt. cl / 2._f_real)) info(2) = 1._f_real

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! create regular mesh (for FFT) and associated cartesian topology: random field will be computed on this mesh and then interpolated
      ! on the external one
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! "gs"/"ge" are specified in "m_scarflib_common" and store first/last global index along each direction for all processes
      ALLOCATE(gs(3, 0:world_size-1), ge(3, 0:world_size-1))

      ! return processors grid
      CALL best_config(dims)

      ! create topology
      CALL mpi_cart_create(mpi_comm_world, ndims, dims, isperiodic, reorder, cartopo, ierr)

      ! return process coordinates in current topology
      CALL mpi_cart_coords(cartopo, world_rank, ndims, coords, ierr)

      ! return first/last index ("ls"/"le") along each direction for the calling process. Note: first point has "ls = 1".
      CALL coords2index(npts, dims, coords, ls, le)

      gs(:, world_rank) = ls
      ge(:, world_rank) = le

      ! make all processes aware of global indices
      CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, gs, 3, mpi_integer, mpi_comm_world, ierr)
      CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, ge, 3, mpi_integer, mpi_comm_world, ierr)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! determine if calling process has any external grid point belonging to other processes and which processes have external grid
      ! points belonging to calling process. These info are needed to define communicator groups when interpolating data.
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ALLOCATE(points_world2rank(0:world_size - 1), points_rank2world(0:world_size - 1))

      ! min/max internal grid index containing input points
      ps = FLOOR(((fs - 1) * ds - off_axis) / dh) + 1
      pe = FLOOR(((fe - 1) * ds - off_axis) / dh) + 1

      points_rank2world = 0

      ! determine if calling process has at least one point falling inside domain of "i-th" process
      DO i = 0, world_size - 1
        DO j = 1, 3
          bool(j) = ( (ps(j) .lt. ge(j, i)) .and. (ps(j) .ge. (gs(j, i) - 1)) ) .or.     &
                    ( (pe(j) .lt. ge(j, i)) .and. (pe(j) .ge. (gs(j, i) - 1)) ) .or.     &
                    ( (pe(j) .ge. ge(j, i)) .and. (ps(j) .lt. (gs(j, i) - 1)) )
        ENDDO
        IF (ALL(bool .eqv. .true.)) points_rank2world(i) = 1
      ENDDO

      ! determine IF process "i" has at least one point falling inside domain of calling process
      CALL mpi_alltoall(points_rank2world, 1, mpi_integer, points_world2rank, 1, mpi_integer, mpi_comm_world, ierr)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! compute power spectrum and inverse FFT
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! number of points for calling process
      DO i = 1, 3
        m(i) = le(i) - ls(i) + 1
      ENDDO

      ! allocate resources to hold spectrum
      ALLOCATE(spec(m(1), m(2), m(3)))

      ! compute spectrum and apply hermitian symmetry
      CALL compute_spectrum(ds, ls, le, dh, acf, cl, sigma, hurst, seed, spec, et)

#ifdef TIMING
      info(5:6) = et

      ! start timer
      CALL watch_start(tictoc)
#endif

      ! (inverse) transform along each direction
      CALL transform_along_z(spec)
      CALL transform_along_y(spec)
      CALL transform_along_x(spec)

#ifdef TIMING
      CALL watch_stop(tictoc)

      info(7) = tictoc
#endif

      ! define scaling factor
      scaling = 1._f_real / SQRT(REAL(npts(1), f_real) * REAL(npts(2), f_real) * REAL(npts(3), f_real) * dh**3)

      ALLOCATE(buffer(m(1), m(2), m(3)))

      ! scale field
      DO k = 1, m(3)
        DO j = 1, m(2)
          DO i = 1, m(1)
            buffer(i, j, k) = REAL(spec(i, j, k), f_real) * scaling
          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE(spec)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! compute mean and variance of the whole random field. Then rescale, if desired, discrete std. dev. to its continuous counterpart
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ALLOCATE(var(0:world_size - 1), mu(0:world_size - 1), npoints(0:world_size - 1))

      ! variance and mean for each single process
      var(world_rank) = variance(buffer)
      mu(world_rank)  = mean(buffer)

      ! number of grid points per process
      npoints(world_rank) = PRODUCT(m)

      ! share results
      CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, var, 1, real_type, mpi_comm_world, ierr)
      CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, mu, 1, real_type, mpi_comm_world, ierr)
      CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, npoints, 1, mpi_integer, mpi_comm_world, ierr)

      ! return total variance in "info(3)" and mean in "info(4)"
      CALL parallel_variance(var, mu, npoints, info(3), info(4))

      ! we want standard deviation
      info(3) = SQRT(info(3))

      IF (rescale .eq. 1) THEN

        scaling = sigma / info(3)

        DO k = 1, m(3)
          DO j = 1, m(2)
            DO i = 1, m(1)
              buffer(i, j, k) = buffer(i, j, k) * scaling - info(4)
            ENDDO
          ENDDO
        ENDDO

        info(3) = sigma
        info(4) = 0._f_real

      ENDIF

      DEALLOCATE(var, mu, npoints)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! taper/mute random field at desired locations.
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      DO i = 1, SIZE(poi, 2)
        pt(:) = poi(:, i) - off_axis(:)                           !< update "poi" according to coordinates of internal grid
        CALL tapering(dh, ls, le, buffer, pt, mute, taper)
      ENDDO

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! each process exchange data with neighbors to define halos to be used during interpolation
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! copy random field to array with halo
      ALLOCATE(delta(0:m(1), 0:m(2), 0:m(3)))

      delta = 0._f_real

      ! copy field
      DO k = 1, m(3)
        DO j = 1, m(2)
          DO i = 1, m(1)
            delta(i, j, k) = buffer(i, j, k)
          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE(buffer)

#ifdef TIMING
      CALL watch_start(tictoc)
#endif

      ! exchange halo
      CALL exchange_halo(cartopo, delta)

      ! interpolate random field values at desired output locations
      CALL grid_mapping(dh, delta, ds, fs, fe, field)

#ifdef TIMING
      CALL watch_stop(tictoc)

      info(8) = tictoc
#endif

      ! release memory
      DEALLOCATE(delta)
      DEALLOCATE(gs, ge)
      DEALLOCATE(points_rank2world, points_world2rank)

      CALL mpi_comm_free(cartopo, ierr)

      CALL mpi_barrier(mpi_comm_world, ierr)

    END SUBROUTINE scarf3d_structured_fim

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE grid_mapping_unstruct(dh, delta, x, y, z, field)

      ! Purpose:
      ! To map (interpolate) the random field from the internal (FFT) grid to an external unstructured one. Take advantage of asynchronous
      ! MPI calls to overlap communitions and calculations as much as possible
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_real),                                      INTENT(IN)  :: dh                             !< fft grid-step
      REAL(f_real),                DIMENSION(:,:,:),     INTENT(IN)  :: delta                          !< random perturbations on fft grid
      REAL(f_real),                DIMENSION(:),         INTENT(IN)  :: x, y, z                        !< points location
      REAL(f_real),                DIMENSION(:),         INTENT(OUT) :: field                          !< random field
      INTEGER(f_int)                                                 :: i, j
      INTEGER(f_int)                                                 :: m, n, np, nb
      INTEGER(f_int)                                                 :: maxwidth, blockwidth, nblocks
      INTEGER(f_int)                                                 :: ierr, trip, req
      INTEGER(f_int),              DIMENSION(world_size)             :: sendcounts, recvcounts, sdispls, rdispls
      INTEGER(f_int),              DIMENSION(world_size)             :: sendcounts_n, recvcounts_n
      INTEGER(f_int), ALLOCATABLE, DIMENSION(:)                      :: map, map_n
      INTEGER(f_int), ALLOCATABLE, DIMENSION(:)                      :: p0, p1
      REAL(f_real),   ALLOCATABLE, DIMENSION(:)                      :: rv, sv
      REAL(f_real),   ALLOCATABLE, DIMENSION(:,:)                    :: sxyz, rxyz, sxyz_n

      !-----------------------------------------------------------------------------------------------------------------------------

      field = 0._f_real

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! in order not to exceed the previous memory peak, we compute the size of the minimum FFT grid block "m", the maximum number of
      ! points to be interpolated "np" amongst all processes and the maximum number of processes "nb" that could send points to another
      ! one. We use these quantities to provide an upper bound "blockwidth" to the number of points that can be handled by any process
      ! at once. As a secondary benefit, this allow us to allocate arrays only once.
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      m = HUGE(1)

      ! minimum computational domain size
      DO i = 0, world_size - 1
        m = MIN(m, (ge(1, i) - gs(1, i) + 1) * (ge(2, i) - gs(2, i) + 1) * (ge(3, i) - gs(3, i) + 1))
      ENDDO

      ! number of processes with at least one point falling inside computational domain of calling process
      nb = COUNT(points_world2rank .eq. 1)

      ! maximum number of processes with at least one point falling inside computational domain of calling process (conservative estimate)
      CALL mpi_allreduce(mpi_in_place, nb, 1, mpi_integer, mpi_max, mpi_comm_world, ierr)

      np = SIZE(x)

      ! maximum number of points aby process could send or receive (conservative estimate)
      CALL mpi_allreduce(mpi_in_place, np, 1, mpi_integer, mpi_max, mpi_comm_world, ierr)

      ! max number of points that can be handled at once (*_n arrays for non-blocking "alltoallv" are included)
      maxwidth = (m * 2) / (9 + 4 * nb)

      blockwidth = maxwidth

      ! number of blocks required to cover ALL points
      nblocks = np / blockwidth

      ! add one more block if "np" is not multiple of "blockwidth"
      IF (MOD(np, blockwidth) .ne. 0) nblocks = nblocks + 1

      ALLOCATE(p0(nblocks), p1(nblocks))

      ! define first/last point index for each block
      DO i = 1, nblocks
        p0(i) = (i - 1) * blockwidth + 1
        p1(i) = i * blockwidth
      ENDDO

      ! make sure that indices DO not exceed actual number of points in calling process
      n = SIZE(x)

      ! limit points indices to maximum number of points available to calling process
      DO i = 1, nblocks
        IF ( (p0(i) .le. n) .and. (p1(i) .gt. n) ) p1(i) = n
        IF ( (p0(i) .gt. n) .and. (p1(i) .gt. n) ) p1(i) = p0(i) - 1
      ENDDO

      ! the strategy is to allocate only once these array, even if they may result larger than necessary
      ALLOCATE(sxyz(3, blockwidth), map(blockwidth), rv(blockwidth))
      ALLOCATE(rxyz(3, blockwidth * nb), sv(blockwidth * nb))
      ALLOCATE(sxyz_n(3, blockwidth), map_n(blockwidth))

      CALL mpi_type_contiguous(3, real_type, trip, ierr)
      CALL mpi_type_commit(trip, ierr)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! for each block, we build an ordered list of points to be sent to all other processes and their number. Then we compute the
      ! number of points to be received. These points are interpolated and sent back to the right owner. Similarly, points interpolated
      ! by other processes are collected and stored at the right location.
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! collect data
      CALL order_points(p0(1), p1(1), dh, x, y, z, sxyz, sendcounts, map)

      ! determine number of points to be received from each process
      CALL mpi_alltoall(sendcounts, 1, mpi_integer, recvcounts, 1, mpi_integer, mpi_comm_world, ierr)

      ! set displacements
      rdispls(1) = 0
      sdispls(1) = 0

      DO i = 2, world_size
        rdispls(i) = rdispls(i - 1) + recvcounts(i - 1)
        sdispls(i) = sdispls(i - 1) + sendcounts(i - 1)
      ENDDO

      DO j = 2, nblocks

        ! forward points referred to previous (j - 1) iteration
        CALL mpi_ialltoallv(sxyz, sendcounts, sdispls, trip, rxyz, recvcounts, rdispls, trip, mpi_comm_world, req, ierr)

        ! collect data to be used in the next (j) iteration
        CALL order_points(p0(j), p1(j), dh, x, y, z, sxyz_n, sendcounts_n, map_n)

        CALL mpi_alltoall(sendcounts_n, 1, mpi_integer, recvcounts_n, 1, mpi_integer, mpi_comm_world, ierr)

        ! wait for data from previous (j - 1) iteration
        CALL mpi_wait(req, mpi_status_ignore, ierr)

        ! number of points received: need to intepolate only these, not whole "rxyz" vector
        n = SUM(recvcounts)

        ! interpolation
        CALL interpolate(n, rxyz, delta, sv)

        ! collect points that have been interpolated by other processes
        CALL mpi_alltoallv(sv, recvcounts, rdispls, real_type, rv, sendcounts, sdispls, real_type, mpi_comm_world, ierr)

        ! copy interpolated data (if any) to right location
        DO i = 1, p1(j - 1) - p0(j - 1) + 1
          field(map(i)) = rv(i)
        ENDDO

        ! assign data to be sent at next iteration
        DO i = 1, blockwidth
          sxyz(1, i) = sxyz_n(1, i)
          sxyz(2, i) = sxyz_n(2, i)
          sxyz(3, i) = sxyz_n(3, i)
          map(i)     = map_n(i)
        ENDDO

        rdispls(1) = 0
        sdispls(1) = 0

        sendcounts(1) = sendcounts_n(1)
        recvcounts(1) = recvcounts_n(1)

        DO i = 2, world_size
          sdispls(i)    = sdispls(i - 1) + sendcounts_n(i - 1)
          rdispls(i)    = rdispls(i - 1) + recvcounts_n(i - 1)
          sendcounts(i) = sendcounts_n(i)
          recvcounts(i) = recvcounts_n(i)
        ENDDO

      ENDDO

      ! forward data
      CALL mpi_alltoallv(sxyz, sendcounts, sdispls, trip, rxyz, recvcounts, rdispls, trip, mpi_comm_world, ierr)

      ! number of points received: need to intepolate only these, not whole "rxyz" vector
      n = SUM(recvcounts)

      ! interpolation
      CALL interpolate(n, rxyz, delta, sv)

      ! backward data
      CALL mpi_alltoallv(sv, recvcounts, rdispls, real_type, rv, sendcounts, sdispls, real_type, mpi_comm_world, ierr)

      ! copy interpolated data to right location. limit the maximum iteration index in case "np" not multiple of "blockwidth"
      DO i = 1, p1(nblocks) - p0(nblocks) + 1
        field(map(i)) = rv(i)
      ENDDO

      ! free resources
      CALL mpi_type_free(trip, ierr)

      DEALLOCATE(sxyz, rxyz, map, rv, sv)
      DEALLOCATE(sxyz_n, map_n)
      DEALLOCATE(p0, p1)

    END SUBROUTINE grid_mapping_unstruct

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE grid_mapping_struct(dh, delta, ds, fs, fe, field)

      ! Purpose:
      ! To map (interpolate) the random field from the internal (FFT) grid to an external structured one. Take advantage of asynchronous
      ! MPI calls to overlap communitions and calculations as much as possible
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_real),                                      INTENT(IN)  :: dh                             !< internal grid grid-step
      REAL(f_real),                DIMENSION(:,:,:),     INTENT(IN)  :: delta                          !< random field on internal grid
      REAL(f_real),                                      INTENT(IN)  :: ds                             !< external grid grid-step
      INTEGER(f_int),              DIMENSION(3),         INTENT(IN)  :: fs, fe                         !< min/max external grid index
      REAL(f_real),                DIMENSION(:,:,:),     INTENT(OUT) :: field
      INTEGER(f_int)                                                 :: i, j, c
      INTEGER(f_int)                                                 :: m, n, np, nb, nx, ny, nz
      INTEGER(f_int)                                                 :: maxwidth, blockwidth, nblocks
      INTEGER(f_int)                                                 :: ierr, trip, req
      INTEGER(f_int),              DIMENSION(3)                      :: dum
      INTEGER(f_int),              DIMENSION(world_size)             :: sendcounts, recvcounts, sdispls, rdispls
      INTEGER(f_int),              DIMENSION(world_size)             :: sendcounts_n, recvcounts_n
      INTEGER(f_int), ALLOCATABLE, DIMENSION(:,:)                    :: map, map_n
      INTEGER(f_int), ALLOCATABLE, DIMENSION(:,:)                    :: p0, p1
      REAL(f_real),   ALLOCATABLE, DIMENSION(:)                      :: rv, sv
      REAL(f_real),   ALLOCATABLE, DIMENSION(:,:)                    :: sxyz, rxyz, sxyz_n

      !-----------------------------------------------------------------------------------------------------------------------------

      field = 0._f_real

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! in order not to exceed the previous memory peak, we compute the size of the minimum FFT grid block "m", the maximum number of
      ! points to be interpolated "np" amongst all processes and the maximum number of processes "nb" that could send points to another
      ! one. We use these quantities to provide an upper bound "blockwidth" to the number of points that can be handled by any process
      ! at once. As a secondary benefit, this allow us to allocate arrays only once.
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      m = HUGE(1)

      ! minimum computational domain SIZE
      DO i = 0, world_size - 1
        m = MIN(m, (ge(1, i) - gs(1, i) + 1) * (ge(2, i) - gs(2, i) + 1) * (ge(3, i) - gs(3, i) + 1))
      ENDDO

      ! number of processes with at least one point falling inside computational domain of calling process
      nb = COUNT(points_world2rank .eq. 1)

      ! maximum number of processes with at least one point falling inside computational domain of calling process (conservative estimate)
      CALL mpi_allreduce(mpi_in_place, nb, 1, mpi_integer, mpi_max, mpi_comm_world, ierr)

      nx = fe(1) - fs(1) + 1
      ny = fe(2) - fs(2) + 1
      nz = fe(3) - fs(3) + 1

      dum = [nx, ny, nz]

      ! maximum number of points ANY process could send or receive (conservative estimate)
      CALL mpi_allreduce(mpi_in_place, dum, 3, mpi_integer, mpi_max, mpi_comm_world, ierr)

      nx = dum(1); ny = dum(2); nz = dum(3)

      ! max number of points that can be handled at once (*_n arrays for non-blocking "alltoallv" are included)
      maxwidth = (m * 2) / (9 + 4 * nb)

      ! number of z-planes that can be processed at once
      blockwidth = maxwidth / (nx * ny)

      ! temporary workaround if no z-planes can be preocessed: this may occurr if input model dimensions are very small!
      IF (blockwidth .eq. 0) THEN
        blockwidth = 1
        maxwidth   = nx * ny
      ENDIF

      ! number of blocks required to cover highest number of z-planes
      nblocks = nz / blockwidth

      IF (MOD(nz, blockwidth) .ne. 0) nblocks = nblocks + 1

      ALLOCATE(p0(3, nblocks), p1(3, nblocks))

      ! define first/last point index for each block
      DO i = 1, nblocks
        p0(1, i) = 1
        p1(1, i) = fe(1) - fs(1) + 1
        p0(2, i) = 1
        p1(2, i) = fe(2) - fs(2) + 1
        p0(3, i) = (i - 1) * blockwidth + 1
        p1(3, i) = i * blockwidth
      ENDDO

      ! make sure that indices do not exceed actual number of z-planes in calling process
      n = fe(3) - fs(3) + 1

      DO i = 1, nblocks
        IF ( (p0(3, i) .le. n) .and. (p1(3, i) .gt. n) ) p1(3, i) = n
        IF ( (p0(3, i) .gt. n) .and. (p1(3, i) .gt. n) ) p1(3, i) = p0(3, i) - 1
      ENDDO

      ! the strategy is to allocate only once these array, even if they may result larger than necessary
      ALLOCATE(sxyz(3, maxwidth), map(3, maxwidth), rv(maxwidth))
      ALLOCATE(rxyz(3, maxwidth * nb), sv(maxwidth * nb))
      ALLOCATE(sxyz_n(3, maxwidth), map_n(3, maxwidth))

      CALL mpi_type_contiguous(3, real_type, trip, ierr)
      CALL mpi_type_commit(trip, ierr)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! for each block, we build an ordered list of points to be sent to all other processes and their number. Then we compute the
      ! number of points to be received. These points are interpolated and sent back to the right owner. Similarly, points interpolated
      ! by other processes are collected and stored at the right location.
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! collect data
      CALL order_points(p0(:, 1), p1(:, 1), dh, fs, fe, ds, sxyz, sendcounts, map)

      ! determine number of points to be received from each process
      CALL mpi_alltoall(sendcounts, 1, mpi_integer, recvcounts, 1, mpi_integer, mpi_comm_world, ierr)

      ! set displacements
      rdispls(1) = 0
      sdispls(1) = 0

      DO i = 2, world_size
        rdispls(i) = rdispls(i - 1) + recvcounts(i - 1)
        sdispls(i) = sdispls(i - 1) + sendcounts(i - 1)
      ENDDO

      DO j = 2, nblocks

        ! forward old data
        CALL mpi_ialltoallv(sxyz, sendcounts, sdispls, trip, rxyz, recvcounts, rdispls, trip, mpi_comm_world, req, ierr)

        ! collect new data
        CALL order_points_struct(p0(:, j), p1(:, j), dh, fs, fe, ds, sxyz_n, sendcounts_n, map_n)
        CALL mpi_alltoall(sendcounts_n, 1, mpi_integer, recvcounts_n, 1, mpi_integer, mpi_comm_world, ierr)

        ! wait for old data
        CALL mpi_wait(req, mpi_status_ignore, ierr)

        ! number of points received: need to intepolate only these, not whole "rxyz" vector
        n = SUM(recvcounts)

        ! interpolation
        CALL interpolate(n, rxyz, delta, sv)

        ! backward data
        CALL mpi_alltoallv(sv, recvcounts, rdispls, real_type, rv, sendcounts, sdispls, real_type, mpi_comm_world, ierr)

        ! copy interpolated data to right location
        DO i = 1, PRODUCT(p1(:, j - 1) - p0(:, j - 1) + 1)
          field(map(1, i), map(2, i), map(3, i)) = rv(i)
        ENDDO

        ! assign data to be sent at next iteration
        DO i = 1, PRODUCT(p1(:, j) - p0(:, j) + 1)
          sxyz(1, i) = sxyz_n(1, i)
          sxyz(2, i) = sxyz_n(2, i)
          sxyz(3, i) = sxyz_n(3, i)
          map(1, i)  = map_n(1, i)
          map(2, i)  = map_n(2, i)
          map(3, i)  = map_n(3, i)
        ENDDO

        rdispls(1) = 0
        sdispls(1) = 0

        sendcounts(1) = sendcounts_n(1)
        recvcounts(1) = recvcounts_n(1)

        DO i = 2, world_size
          sdispls(i)    = sdispls(i - 1) + sendcounts_n(i - 1)
          rdispls(i)    = rdispls(i - 1) + recvcounts_n(i - 1)
          sendcounts(i) = sendcounts_n(i)
          recvcounts(i) = recvcounts_n(i)
        ENDDO

      ENDDO

      ! forward data
      CALL mpi_alltoallv(sxyz, sendcounts, sdispls, trip, rxyz, recvcounts, rdispls, trip, mpi_comm_world, ierr)

      ! number of points received: need to intepolate only these, not whole "rxyz" vector
      n = SUM(recvcounts)

      ! interpolation
      CALL interpolate(n, rxyz, delta, sv)

      ! backward data
      CALL mpi_alltoallv(sv, recvcounts, rdispls, real_type, rv, sendcounts, sdispls, real_type, mpi_comm_world, ierr)

      ! copy interpolated data to right location. Limit the maximum iteration index in case "np" not multiple of "blockwidth"
      DO i = 1, PRODUCT(p1(:, nblocks) - p0(:, nblocks) + 1)
        field(map(1, i), map(2, i), map(3, i)) = rv(i)
      ENDDO

      ! free resources
      CALL mpi_type_free(trip, ierr)

      DEALLOCATE(sxyz, rxyz, map, rv, sv)
      DEALLOCATE(sxyz_n, map_n)
      DEALLOCATE(p0, p1)

    END SUBROUTINE grid_mapping_struct

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE order_points_unstruct(p0, p1, dh, x, y, z, xyz, sendcounts, map)

      ! Purpose:
      ! To reorder external grid nodes to be sent to each process in rank-increasing order, storing indices of the nearest internal
      ! (FFT) grid nodes, to return the number of points to be sent to each process and to store the original index of the external
      ! grid nodes.
      ! Array 'xyz' stores floating point numbers instead of integers to speed-up interpolation algorithm if bilinear interpolation
      ! is used.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(f_int),                 INTENT(IN)  :: p0, p1             !< first/last index of points to be sent
      REAL(f_real),                   INTENT(IN)  :: dh                 !< grid step of internal (FFT) grid
      REAL(f_real),   DIMENSION(:),   INTENT(IN)  :: x, y, z            !< vector of points (external grid)
      REAL(f_real),   DIMENSION(:,:), INTENT(OUT) :: xyz                !< list with nearest internal grid nodes
      INTEGER(f_int), DIMENSION(:),   INTENT(OUT) :: sendcounts         !< number of points to be sent to each process
      INTEGER(f_int), DIMENSION(:),   INTENT(OUT) :: map                !< position of a point inside vector
      INTEGER(f_int)                              :: l, p, n, c
      INTEGER(f_int), DIMENSION(3)                :: const, m
      LOGICAL                                     :: bool
      REAL(f_real)                                :: i, j, k

      !-----------------------------------------------------------------------------------------------------------------------------

      sendcounts(:) = 0
      map(:)        = 0

      c = 0

      ! loop over processes
      DO l = 0, world_size - 1

        ! skip rest of the loop if calling process has no points falling inside domain of "l-th" process
        IF (points_rank2world(l) .eq. 0) CYCLE

        DO p = 1, 3
          const(p) = gs(p, l)
          m(p)     = ge(p, l) - gs(p, l) + 2
        ENDDO

        n = c

        DO p = p0, p1

          ! transform point coordinates into grid indices for l-th process
          i = (x(p) - off_axis(1)) / dh + 3._f_real - const(1)
          j = (y(p) - off_axis(2)) / dh + 3._f_real - const(2)
          k = (z(p) - off_axis(3)) / dh + 3._f_real - const(3)

          ! check if point is within block of l-th process
          bool = (i .ge. 1) .and. (i .lt. m(1)) .and. (j .ge. 1) .and. (j .lt. m(2)) .and. (k .ge. 1) .and. (k .lt. m(3))

          IF (bool .eqv. .true.) THEN
            c         = c + 1
            xyz(:, c) = [i, j, k]
            map(c)    = p
          ENDIF

        ENDDO

        ! number of points to be sent to l-th process
        sendcounts(l + 1) = c - n

      ENDDO

    END SUBROUTINE order_points_unstruct

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE order_points_struct(p0, p1, dh, fs, fe, ds, xyz, sendcounts, map)

      ! Purpose:
      ! To reorder external grid nodes to be sent to each process in rank-increasing order, storing indices of the nearest internal
      ! (FFT) grid nodes, to return the number of points to be sent to each process and to store the original index of the external
      ! grid nodes.
      ! Array 'xyz' stores floating point numbers instead of integers to speed-up interpolation algorithm if bilinear interpolation
      ! is used.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(f_int), DIMENSION(3),          INTENT(IN)  :: p0, p1              !< first/last index (internal FFT grid) of points to be sent
      REAL(f_real),                          INTENT(IN)  :: dh                  !< grid step of internal (FFT) grid
      INTEGER(f_int), DIMENSION(3),          INTENT(IN)  :: fs, fe              !< first/last index of external grid for calling process
      REAL(f_real),                          INTENT(IN)  :: ds                  !< grid-step external grid
      REAL(f_real),   DIMENSION(:,:),        INTENT(OUT) :: xyz                 !< list with nearest internal grid nodes
      INTEGER(f_int), DIMENSION(:),          INTENT(OUT) :: sendcounts          !< number of points to be sent to each process
      INTEGER(f_int), DIMENSION(:,:),        INTENT(OUT) :: map                 !< position of a point inside vector
      INTEGER(f_int)                                     :: l, p, n, c
      INTEGER(f_int)                                     :: nx, nxy
      INTEGER(f_int)                                     :: i0, j0, k0
      INTEGER(f_int), DIMENSION(3)                       :: const, m
      LOGICAL                                            :: bool
      REAL(f_real)                                       :: i, j, k
      REAL(f_real),   DIMENSION(p0(3):p1(3))             :: z
      REAL(f_real),   DIMENSION(p0(2):p1(2))             :: y
      REAL(f_real),   DIMENSION(p0(1):p1(1))             :: x

      !-----------------------------------------------------------------------------------------------------------------------------

      ! transform point coordinates into grid indices
      DO k0 = p0(3), p1(3)
        z(k0) = ((fs(3) - 1 + k0 - 1) * ds - off_axis(3)) / dh + 3._f_real
      ENDDO

      DO j0 = p0(2), p1(2)
        y(j0) = ((fs(2) - 1 + j0 - 1) * ds - off_axis(2)) / dh + 3._f_real
      ENDDO

      DO i0 = p0(1), p1(1)
        x(i0) = ((fs(1) - 1 + i0 - 1) * ds - off_axis(1)) / dh + 3._f_real
      ENDDO

      sendcounts(:) = 0
      map(:,:)      = 0

      c = 0

      ! loop over processes
      DO l = 0, world_size - 1

        ! skip rest of the loop if calling process has no points falling inside domain of "l-th" process
        IF (points_rank2world(l) .eq. 0) CYCLE

        DO p = 1, 3
          const(p) = gs(p, l)
          m(p)     = ge(p, l) - gs(p, l) + 2
        ENDDO

        n = c

        DO k0 = p0(3), p1(3)

          k = z(k0) - const(3)

          DO j0 = p0(2), p1(2)

            j = y(j0) - const(2)

            DO i0 = p0(1), p1(1)

              i = x(i0) - const(1)

              ! check if point is within block of l-th process
              bool = (i .ge. 1) .and. (i .lt. m(1)) .and. (j .ge. 1) .and. (j .lt. m(2)) .and. (k .ge. 1) .and. (k .lt. m(3))

              IF (bool .eqv. .true.) THEN
                c         = c + 1
                xyz(:, c) = [i, j, k]
                map(:, c) = [i0, j0, k0]
              ENDIF

            ENDDO

          ENDDO

        ENDDO

        ! number of points to be sent to l-th process
        sendcounts(l + 1) = c - n

      ENDDO

    END SUBROUTINE order_points_struct

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE transform_along_z(spectrum)

      ! Purpose:
      ! To compute many in-place complex-to-complex inverse FFT of 3D spectrum array along the z-axis. Processes are gathered in pencils
      ! oriented along the z-axis. Inside each pencil, inverse FFT and data exchange occurs on slices perpendicular to the y-axis in
      ! order to save memory.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      COMPLEX(f_real),              DIMENSION(:,:,:), INTENT(INOUT) :: spectrum                   !< input/output spectrum
      COMPLEX(f_real), ALLOCATABLE, DIMENSION(:)                    :: sendbuf, recvbuf
      INTEGER(f_int)                                                :: i, j, k, p, c
      INTEGER(f_int)                                                :: rank, ntasks, ierr
      INTEGER(f_int)                                                :: pencil
      INTEGER(f_int)                                                :: nx, ny, nz
      INTEGER(f_int), ALLOCATABLE, DIMENSION(:)                     :: i0, i1
      INTEGER(f_int), ALLOCATABLE, DIMENSION(:)                     :: k0, k1
      INTEGER(f_int), ALLOCATABLE, DIMENSION(:)                     :: lx, lz
      INTEGER(f_int), ALLOCATABLE, DIMENSION(:)                     :: sendcounts, recvcounts
      INTEGER(f_int), ALLOCATABLE, DIMENSION(:)                     :: sdispls, rdispls
      TYPE(c_ptr)                                                   :: p1

      !-----------------------------------------------------------------------------------------------------------------------------

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! set up pencils oriented along z-axis, prepare vectors for data exchange
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! group processes in pencils oriented along z-axis. "lz" contains the number of points for each process in the pencil
      CALL build_pencil(2, pencil, rank, ntasks, lz)

      ! this array will contain the number of points local to each process in current pencil along x-axis
      ALLOCATE(lx(0:ntasks-1))

      ! number of points along x and y. These are indentical for all processes inside same pencil
      nx = SIZE(spectrum, 1)
      ny = SIZE(spectrum, 2)

      ! total number of points along z-axis
      nz = SUM(lz)

      ALLOCATE(k0(0:ntasks-1), k1(0:ntasks-1))

      ! first/last k-index for each process
      DO p = 0, ntasks - 1
        IF (p .eq. 0) THEN
          k0(p) = 1
          k1(p) = lz(p)
        ELSE
          k0(p) = k1(p - 1) + 1
          k1(p) = k1(p - 1) + lz(p)
        ENDIF
      ENDDO

      ALLOCATE(i0(0:ntasks-1), i1(0:ntasks-1))

      ! distribute points along x-axis between processes
      CALL split_task(nx, ntasks, i0, i1)

      ! number of points along x-axis for each process
      DO p = 0, ntasks - 1
        lx(p) = i1(p) - i0(p) + 1
      ENDDO

      ALLOCATE(sendcounts(0:ntasks-1), sdispls(0:ntasks-1), recvcounts(0:ntasks-1), rdispls(0:ntasks-1))

      ! set counts/displacement arrays
      DO p = 0, ntasks - 1

        sendcounts(p) = lz(rank) * lx(p)              !< points send from process "rank" to process "p"
        recvcounts(p) = lz(p) * lx(rank)              !< points received from process "rank" by process "p"

        IF (p .eq. 0) THEN
          sdispls(p) = 0
          rdispls(p) = 0
        ELSE
          sdispls(p) = sdispls(p - 1) + sendcounts(p - 1)
          rdispls(p) = rdispls(p - 1) + recvcounts(p - 1)
        ENDIF

      ENDDO

      ALLOCATE(sendbuf(nx * lz(rank)))
      ALLOCATE(recvbuf(nz * lx(rank)))

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! allocate memory and prepare plan for FFTW library
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! buffer "cdum" will be used for transform and it is allocated using FFTW memory allocation utility
#ifdef DOUBLE_PREC
      pc = fftw_alloc_complex(INT(lx(rank) * nz, c_size_t))
#else
      pc = fftwf_alloc_complex(INT(lx(rank) * nz, c_size_t))
#endif

      CALL C_F_POINTER(pc, cdum, [lx(rank) * nz])

      ! prepare FFTW plan along z-axis ("lx(rank)" transforms, each having "nz" number of points)
#ifdef DOUBLE_PREC
      p1 = fftw_plan_many_dft(1, [nz], lx(rank), cdum, [nz], lx(rank), 1, cdum, [nz], lx(rank), 1, fftw_backward, fftw_estimate)
#else
      p1 = fftwf_plan_many_dft(1, [nz], lx(rank), cdum, [nz], lx(rank), 1, cdum, [nz], lx(rank), 1, fftw_backward, fftw_estimate)
#endif

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! collect data, apply transform and send data back
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! work y-slices one by one
      DO j = 1, ny

        c = 0

        ! rearrange elements of "spectrum"
        DO p = 0, ntasks - 1
          DO k = 1, lz(rank)
            DO i = i0(p), i1(p)
              c          = c + 1
              sendbuf(c) = spectrum(i, j, k)
            ENDDO
          ENDDO
        ENDDO

        ! exchange data
        CALL mpi_alltoallv(sendbuf, sendcounts, sdispls, cplx_type, cdum, recvcounts, rdispls, cplx_type, pencil, ierr)

        ! ifft
#ifdef DOUBLE_PREC
        CALL fftw_execute_dft(p1, cdum, cdum)
#else
        CALL fftwf_execute_dft(p1, cdum, cdum)
#endif

        c = 0

        ! rearrange elements of "cdum" into "recvbuf"
        DO p = 0, ntasks - 1
          DO i = 1, lx(rank)
            DO k = k0(p), k1(p)
              c          = c + 1
              recvbuf(c) = cdum(i + (k - 1)*lx(rank))
            ENDDO
          ENDDO
        ENDDO

        ! send transformed data back
        CALL mpi_alltoallv(recvbuf, recvcounts, rdispls, cplx_type, sendbuf, sendcounts, sdispls, cplx_type, pencil, ierr)

        ! rearrange data into original layout
        DO k = 1, lz(rank)
          DO i = 1, nx
            spectrum(i, j, k) = sendbuf(k + (i - 1)*lz(rank))
          ENDDO
        ENDDO

      ENDDO

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! release resources
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! destroy FFTW plan
#ifdef DOUBLE_PREC
      CALL fftw_destroy_plan(p1)
#else
      CALL fftwf_destroy_plan(p1)
#endif

      ! release memory
      NULLIFY(cdum)

#ifdef DOUBLE_PREC
      CALL fftw_free(pc)
#else
      CALL fftwf_free(pc)
#endif

      DEALLOCATE(sendbuf, recvbuf, sdispls, rdispls, sendcounts, recvcounts)
      DEALLOCATE(i0, i1, k0, k1, lx, lz)

      ! release communicator
      CALL mpi_comm_free(pencil, ierr)

    END SUBROUTINE transform_along_z

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE transform_along_y(spectrum)

      ! Purpose:
      ! To compute many in-place complex-to-complex inverse FFT of 3D spectrum array along the y-axis. Processes are gathered in pencils
      ! oriented along the y-axis. Inside each pencil, inverse FFT and data exchange occurs on slices perpendicular to the z-axis in
      ! order to save memory.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      COMPLEX(f_real),              DIMENSION(:,:,:), INTENT(INOUT) :: spectrum                   !< input/output spectrum
      COMPLEX(f_real), ALLOCATABLE, DIMENSION(:)                    :: sendbuf, recvbuf
      INTEGER(f_int)                                                :: i, j, k, p, c
      INTEGER(f_int)                                                :: rank, ntasks, ierr
      INTEGER(f_int)                                                :: pencil
      INTEGER(f_int)                                                :: nx, ny, nz
      INTEGER(f_int),  ALLOCATABLE, DIMENSION(:)                    :: i0, i1
      INTEGER(f_int),  ALLOCATABLE, DIMENSION(:)                    :: j0, j1
      INTEGER(f_int),  ALLOCATABLE, DIMENSION(:)                    :: lx, ly
      INTEGER(f_int),  ALLOCATABLE, DIMENSION(:)                    :: sendcounts, recvcounts
      INTEGER(f_int),  ALLOCATABLE, DIMENSION(:)                    :: sdispls, rdispls
      TYPE(c_ptr)                                                   :: p1

      !-------------------------------------------------------------------------------------------------------------------------------

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! set up pencils oriented along y-axis, prepare vectors for data exchange
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! group processes in pencils oriented along y-axis. "ly" contains the number of points for each process in the pencil
      CALL build_pencil(1, pencil, rank, ntasks, ly)

      ! this array will contain the number of points local to each process in current pencil along x-axis
      ALLOCATE(lx(0:ntasks-1))

      ! number of points along x and z. These are indentical for all processes inside same pencil
      nx = SIZE(spectrum, 1)
      nz = SIZE(spectrum, 3)

      ! total number of points along y-axis
      ny = SUM(ly)

      ALLOCATE(j0(0:ntasks-1), j1(0:ntasks-1))

      ! first/last j-index for each process
      DO p = 0, ntasks - 1
        IF (p .eq. 0) THEN
          j0(p) = 1
          j1(p) = ly(p)
        ELSE
          j0(p) = j1(p - 1) + 1
          j1(p) = j1(p - 1) + ly(p)
        ENDIF
      ENDDO

      ALLOCATE(i0(0:ntasks-1), i1(0:ntasks-1))

      ! distribute points along x-axis between processes
      CALL split_task(nx, ntasks, i0, i1)

      ! number of points along x-axis for each process
      DO p = 0, ntasks - 1
        lx(p) = i1(p) - i0(p) + 1
      ENDDO

      ALLOCATE(sendcounts(0:ntasks-1), sdispls(0:ntasks-1), recvcounts(0:ntasks-1), rdispls(0:ntasks-1))

      ! set counts/displacement arrays
      DO p = 0, ntasks - 1

        sendcounts(p) = ly(rank) * lx(p)              !< points send from process "rank" to process "p"
        recvcounts(p) = ly(p) * lx(rank)              !< points received from process "rank" by process "p"

        IF (p .eq. 0) THEN
          sdispls(p) = 0
          rdispls(p) = 0
        ELSE
          sdispls(p) = sdispls(p - 1) + sendcounts(p - 1)
          rdispls(p) = rdispls(p - 1) + recvcounts(p - 1)
        ENDIF

      ENDDO

      ALLOCATE(sendbuf(nx * ly(rank)))
      ALLOCATE(recvbuf(ny * lx(rank)))

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! allocate memory and prepare plan for FFTW library
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! buffer "cdum" will be used for transform and it is allocated using FFTW memory allocation utility
#ifdef DOUBLE_PREC
      pc = fftw_alloc_complex(INT(lx(rank) * ny, c_size_t))
#else
      pc = fftwf_alloc_complex(INT(lx(rank) * ny, c_size_t))
#endif

      CALL C_F_POINTER(pc, cdum, [lx(rank) * ny])

      ! prepare FFTW plan along y-axis ("lx(rank)" transforms, each having "ny" number of points)
#ifdef DOUBLE_PREC
      p1 = fftw_plan_many_dft(1, [ny], lx(rank), cdum, [ny], lx(rank), 1, cdum, [ny], lx(rank), 1, fftw_backward, fftw_estimate)
#else
      p1 = fftwf_plan_many_dft(1, [ny], lx(rank), cdum, [ny], lx(rank), 1, cdum, [ny], lx(rank), 1, fftw_backward, fftw_estimate)
#endif

      ! work z-slices one by one
      DO k = 1, nz

        c = 0

        ! rearrange elements of "spectrum"
        DO p = 0, ntasks - 1
          DO j = 1, ly(rank)
            DO i = i0(p), i1(p)
              c          = c + 1
              sendbuf(c) = spectrum(i, j, k)
            ENDDO
          ENDDO
        ENDDO

        ! exchange data
        CALL mpi_alltoallv(sendbuf, sendcounts, sdispls, cplx_type, cdum, recvcounts, rdispls, cplx_type, pencil, ierr)

        ! ifft
#ifdef DOUBLE_PREC
        CALL fftw_execute_dft(p1, cdum, cdum)
#else
        CALL fftwf_execute_dft(p1, cdum, cdum)
#endif

        c = 0

        ! rearrange elements of "cdum" into "recvbuf"
        DO p = 0, ntasks - 1
          DO i = 1, lx(rank)
            DO j = j0(p), j1(p)
              c          = c + 1
              recvbuf(c) = cdum(i + (j - 1)*lx(rank))
            ENDDO
          ENDDO
        ENDDO

        ! send transformed data back
        CALL mpi_alltoallv(recvbuf, recvcounts, rdispls, cplx_type, sendbuf, sendcounts, sdispls, cplx_type, pencil, ierr)

        ! rearrange data into original layout
        DO j = 1, ly(rank)
          DO i = 1, nx
            spectrum(i, j, k) = sendbuf(j + (i - 1)*ly(rank))
          ENDDO
        ENDDO

      ENDDO

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! release resources
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! destroy FFTW plan
#ifdef DOUBLE_PREC
      CALL fftw_destroy_plan(p1)
#else
      CALL fftwf_destroy_plan(p1)
#endif

      ! release memory
      NULLIFY(cdum)

#ifdef DOUBLE_PREC
      CALL fftw_free(pc)
#else
      CALL fftwf_free(pc)
#endif

      DEALLOCATE(sendbuf, recvbuf, sdispls, rdispls, sendcounts, recvcounts)
      DEALLOCATE(i0, i1, j0, j1, lx, ly)

      ! release communicator
      CALL mpi_comm_free(pencil, ierr)

    END SUBROUTINE transform_along_y

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE transform_along_x(spectrum)

      ! Purpose:
      ! To compute many in-place complex-to-complex inverse FFT of 3D spectrum array along the x-axis. Processes are gathered in pencils
      ! oriented along the x-axis. Inside each pencil, inverse FFT and data exchange occurs on slices perpendicular to the z-axis in
      ! order to save memory.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      COMPLEX(f_real),              DIMENSION(:,:,:), INTENT(INOUT) :: spectrum                       !< input/output spectrum
      COMPLEX(f_real), ALLOCATABLE, DIMENSION(:)                    :: sendbuf
      INTEGER(f_int)                                                :: i, j, k, p, c
      INTEGER(f_int)                                                :: rank, ntasks, ierr
      INTEGER(f_int)                                                :: nyquist, pencil
      INTEGER(f_int)                                                :: imax, p_imax
      INTEGER(f_int)                                                :: nx, ny, nz
      INTEGER(f_int),  ALLOCATABLE, DIMENSION(:)                    :: offset
      INTEGER(f_int),  ALLOCATABLE, DIMENSION(:)                    :: i0, i1
      INTEGER(f_int),  ALLOCATABLE, DIMENSION(:)                    :: j0, j1
      INTEGER(f_int),  ALLOCATABLE, DIMENSION(:)                    :: lx, ly
      INTEGER(f_int),  ALLOCATABLE, DIMENSION(:)                    :: sendcounts, recvcounts
      INTEGER(f_int),  ALLOCATABLE, DIMENSION(:)                    :: sdispls, rdispls
      INTEGER(f_int),  ALLOCATABLE, DIMENSION(:)                    :: r_sendcounts, r_recvcounts
      INTEGER(f_int),  ALLOCATABLE, DIMENSION(:)                    :: r_sdispls, r_rdispls
      REAL(f_real),    ALLOCATABLE, DIMENSION(:)                    :: r_sendbuf, recvbuf
      TYPE(c_ptr)                                                   :: p1

      !-------------------------------------------------------------------------------------------------------------------------------

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! set up pencils oriented along x-axis, prepare vectors for data exchange
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! index of nyquist frequency along x-axis
      nyquist = npts(1) / 2 + 1

      ! group processes in pencils oriented along x-axis. "lx" contains number of points for each process in the pencil
      CALL build_pencil(0, pencil, rank, ntasks, lx)

      ! parameter used to determine elements to be sent/received
      ALLOCATE(offset(0:ntasks-1))

      offset(rank) = gs(1, world_rank) - 1

      CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, offset, 1, mpi_integer, pencil, ierr)

      ! this array will contain the number of points local to each process in current pencil along y-axis
      ALLOCATE(ly(0:ntasks-1))

      ! number of points along y and z. these are indentical for all processes inside same pencil
      ny = SIZE(spectrum, 2)
      nz = SIZE(spectrum, 3)

      ! total number of points along x-axis
      nx = SUM(lx)

      ! max i-index used to send data
      imax = MIN(lx(rank) + offset(rank), nyquist) - offset(rank)

      IF (imax .lt. 1) imax = 0

      ALLOCATE(i0(0:ntasks-1), i1(0:ntasks-1))

      ! first/last i-index for each process
      DO p = 0, ntasks - 1
        IF (p .eq. 0) THEN
          i0(p) = 1
          i1(p) = lx(p)
        ELSE
          i0(p) = i1(p - 1) + 1
          i1(p) = i1(p - 1) + lx(p)
        ENDIF
      ENDDO

      ALLOCATE(j0(0:ntasks-1), j1(0:ntasks-1))

      ! distribute points along y-axis between processes
      CALL split_task(ny, ntasks, j0, j1)

      ! number of points along y-axis for each process
      DO p = 0, ntasks - 1
        ly(p) = j1(p) - j0(p) + 1
      ENDDO

      ALLOCATE(sendcounts(0:ntasks-1), sdispls(0:ntasks-1), recvcounts(0:ntasks-1), rdispls(0:ntasks-1))
      ALLOCATE(r_sendcounts(0:ntasks-1), r_sdispls(0:ntasks-1), r_recvcounts(0:ntasks-1), r_rdispls(0:ntasks-1))

      ! set counts/displacement arrays
      DO p = 0, ntasks - 1

        ! define "imax" for process "p"
        p_imax = MIN(lx(p) + offset(p), nyquist) - offset(p)

        IF (p_imax .lt. 1) p_imax = 0

        ! process completely to the right of "nyquist" won't send any element ("sendcounts=0")
        sendcounts(p) = ly(p) * imax                                            !< points send from process "rank" to process "p"

        ! a process won't receive anything from process "p" if the latter is to the right of "nyquist" ("recvcounts=0")
        recvcounts(p) = ly(rank) * p_imax                                       !< points received from process "rank" by process "p"

        r_sendcounts(p) = ly(rank) * lx(p)                                      !< points send from process "rank" to process "p"
        r_recvcounts(p) = ly(p) * lx(rank)                                      !< points received from process "rank" by process "p"

        IF (p .eq. 0) THEN
          sdispls(p)   = 0
          rdispls(p)   = 0
          r_sdispls(p) = 0
          r_rdispls(p) = 0
        ELSE
          sdispls(p)   = sdispls(p - 1) + sendcounts(p - 1)
          rdispls(p)   = rdispls(p - 1) + recvcounts(p - 1)
          r_sdispls(p) = r_sdispls(p - 1) + r_sendcounts(p - 1)
          r_rdispls(p) = r_rdispls(p - 1) + r_recvcounts(p - 1)
        ENDIF

      ENDDO

      ALLOCATE(sendbuf(ny * imax), r_sendbuf(nx * ly(rank)))
      ALLOCATE(recvbuf(ny * lx(rank)))

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! allocate memory and prepare plan for FFTW library (complex-to-real transforms)
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

#ifdef DOUBLE_PREC
      pc = fftw_alloc_complex(int(ly(rank) * nyquist, c_size_t))
#else
      pc = fftwf_alloc_complex(int(ly(rank) * nyquist, c_size_t))
#endif

      ! for in-place transforms, output real array is slightly longer than actual physical dimensions
      CALL C_F_POINTER(pc, cdum, [ly(rank) * nyquist])
      CALL C_F_POINTER(pc, rdum, [ly(rank) * nyquist * 2])

      ! prepare FFTW plan along x-axis ("ly(rank)" transforms, each having "npts(1)" number of points)
#ifdef DOUBLE_PREC
      p1 = fftw_plan_many_dft_c2r(1, [npts(1)], ly(rank), cdum, [npts(1)], ly(rank), 1, rdum, [npts(1)], ly(rank), 1, fftw_estimate)
#else
      p1 = fftwf_plan_many_dft_c2r(1, [npts(1)], ly(rank), cdum, [npts(1)], ly(rank), 1, rdum, [npts(1)], ly(rank), 1,fftw_estimate)
#endif

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! collect data, apply transform and send data back
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! work z-slices one by one
      DO k = 1, nz

        c = 0

        ! rearrange elements of "spectrum"
        DO p = 0, ntasks - 1
          DO i = 1, imax
            DO j = j0(p), j1(p)
              c          = c + 1
              sendbuf(c) = spectrum(i, j, k)
            ENDDO
          ENDDO
        ENDDO

        ! exchange data
        CALL mpi_alltoallv(sendbuf, sendcounts, sdispls, cplx_type, cdum, recvcounts, rdispls, cplx_type, pencil, ierr)

        ! ifft
#ifdef DOUBLE_PREC
        CALL fftw_execute_dft_c2r(p1, cdum, rdum)
#else
        CALL fftwf_execute_dft_c2r(p1, cdum, rdum)
#endif

        c = 0

        ! rearrange elements of "rdum" into "recvbuf"
        DO p = 0, ntasks - 1
          DO j = 1, ly(rank)
            DO i = i0(p), i1(p)
              c            = c + 1
              r_sendbuf(c) = rdum(j + (i - 1)*ly(rank))
            ENDDO
          ENDDO
        ENDDO

        ! send transformed data back
        CALL mpi_alltoallv(r_sendbuf, r_sendcounts, r_sdispls, real_type, recvbuf, r_recvcounts, r_rdispls, real_type, pencil, ierr)

        ! rearrange data into original layout
        DO j = 1, ny
          DO i = 1, lx(rank)
            spectrum(i, j, k) = cmplx(recvbuf(i + (j - 1)*lx(rank)), 0._f_real, f_real)
          ENDDO
        ENDDO

      ENDDO

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! release resources
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! destroy FFTW plan
#ifdef DOUBLE_PREC
      CALL fftw_destroy_plan(p1)
#else
      CALL fftwf_destroy_plan(p1)
#endif

      ! release memory
      NULLIFY(cdum, rdum)

#ifdef DOUBLE_PREC
      CALL fftw_free(pc)
#else
      CALL fftwf_free(pc)
#endif

      DEALLOCATE(sendbuf, recvbuf, sdispls, rdispls, sendcounts, recvcounts)
      DEALLOCATE(r_sendbuf, r_sdispls, r_rdispls, r_sendcounts, r_recvcounts)
      DEALLOCATE(i0, i1, j0, j1, lx, ly, offset)

      ! release communicator
      CALL mpi_comm_free(pencil, ierr)

    END SUBROUTINE transform_along_x

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE build_pencil(dir, newcomm, rank, ntasks, n)

      ! Purpose:
      ! To group processes falling inside the same pencil oriented along a specific direction.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(f_int),                                        INTENT(IN)  :: dir            !< pencil direction (0=x, 1=y, 2=z)
      INTEGER(f_int),                                        INTENT(OUT) :: newcomm        !< handle to new communicator
      INTEGER(f_int),                                        INTENT(OUT) :: rank           !< rank of calling process in new communicator
      INTEGER(f_int),                                        INTENT(OUT) :: ntasks         !< number of processes in new communicator
      INTEGER(f_int), ALLOCATABLE, DIMENSION(:),             INTENT(OUT) :: n              !< number of points per process in pencil direction
      INTEGER(f_int)                                                     :: i, ierr
      INTEGER(f_int),              DIMENSION(0:world_size-1)             :: color
      LOGICAL,                     DIMENSION(2)                          :: bool

      !-----------------------------------------------------------------------------------------------------------------------------

      ! group processes into pencils
      DO i = 0, world_size - 1

        color(i) = 0

        ! pencil oriented along x-axis
        IF (dir .eq. 0) THEN
          bool(1) = (gs(2, i) .eq. gs(2, world_rank)) .and. (ge(2, i) .eq. ge(2, world_rank))
          bool(2) = (gs(3, i) .eq. gs(3, world_rank)) .and. (ge(3, i) .eq. ge(3, world_rank))
        ! pencil oriented along y-axis
        ELSEIF (dir .eq. 1) THEN
          bool(1) = (gs(1, i) .eq. gs(1, world_rank)) .and. (ge(1, i) .eq. ge(1, world_rank))
          bool(2) = (gs(3, i) .eq. gs(3, world_rank)) .and. (ge(3, i) .eq. ge(3, world_rank))
        ! pencil oriented along z-axis
        ELSEIF (dir .eq. 2) THEN
          bool(1) = (gs(1, i) .eq. gs(1, world_rank)) .and. (ge(1, i) .eq. ge(1, world_rank))
          bool(2) = (gs(2, i) .eq. gs(2, world_rank)) .and. (ge(2, i) .eq. ge(2, world_rank))
        ENDIF

        IF (ALL(bool .eqv. .true.)) color(i) = i + 1

      ENDDO

      ! process belonging to the same pencil have same color
      color(world_rank) = MAXVAL(color, dim = 1)

      ! create communicator subgroup
      CALL mpi_comm_split(mpi_comm_world, color(world_rank), world_rank, newcomm, ierr)

      ! process id and communicator size
      CALL mpi_comm_rank(newcomm, rank, ierr)
      CALL mpi_comm_size(newcomm, ntasks, ierr)

      ALLOCATE(n(0:ntasks - 1))

      ! number of points along pencil direction for calling process
      n(rank) = ge(dir + 1, world_rank) - gs(dir + 1, world_rank) + 1

      ! make whole communicator aware of points for each process
      CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, n, 1, mpi_integer, newcomm, ierr)

    END SUBROUTINE build_pencil

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE compute_spectrum(ds, ls, le, dh, acf, cl, sigma, hurst, seed, spec, time)

      ! Purpose:
      ! To compute the power spectral density function (PSDF) corresponding to a specific autocorrelation function (ACF). Spectrum
      ! is tapered if grid-step of external grid is not equal to grid-step of internal (FFT) grid.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_real),                                                    INTENT(IN)  :: ds
      INTEGER(f_int),  DIMENSION(3),                                   INTENT(IN)  :: ls, le                !< global indices
      REAL(f_real),                                                    INTENT(IN)  :: dh                    !< grid-step
      INTEGER(f_int),                                                  INTENT(IN)  :: acf                   !< autocorrelation FUNCTION
      REAL(f_real),    DIMENSION(3),                                   INTENT(IN)  :: cl                    !< correlation length (can be anisotropic)
      REAL(f_real),                                                    INTENT(IN)  :: sigma                 !< standard deviation (sigma%/100)
      REAL(f_real),                                                    INTENT(IN)  :: hurst                 !< hurst exponent
      INTEGER(f_int),                                                  INTENT(IN)  :: seed                  !< initial seed number
      COMPLEX(f_real), DIMENSION(ls(1):le(1),ls(2):le(2),ls(3):le(3)), INTENT(OUT) :: spec                  !< spectrum
      REAL(f_real),    DIMENSION(2),                                   INTENT(OUT) :: time                  !< elapsed time
      INTEGER(f_int)                                                               :: i, j, k
      LOGICAL                                                                      :: filter
      REAL(f_real)                                                                 :: butter, num, amp
      REAL(f_real)                                                                 :: kc, kr
      REAL(f_dble)                                                                 :: tictoc
      REAL(f_real),    DIMENSION(3)                                                :: dk
      REAL(f_real),    DIMENSION(npts(1))                                          :: kx
      REAL(f_real),    DIMENSION(npts(2))                                          :: ky
      REAL(f_real),    DIMENSION(npts(3))                                          :: kz
      REAL(c_real),    DIMENSION(:),                                   POINTER     :: r
      REAL(c_real),    DIMENSION(:,:,:),                               POINTER     :: harvest
      TYPE(c_ptr)                                                                  :: cptr

      !-----------------------------------------------------------------------------------------------------------------------------

      time(:) = 0._f_real

#ifdef TIMING
      ! start timer
      CALL watch_start(tictoc)
#endif

      ! resolution in wavenumber domain along each direction
      DO i = 1, 3
        dk(i) = 2._f_real * pi / (REAL(npts(i), f_real) * dh)
      ENDDO

      ! nyquist wavenumber
      !knyq = pi / dh

      ! corner wavenumber for filtering spectrum is controlled by external mesh grid-step
      ! kc = pi / (ds) * SQRT(3._f_real)
      kc = pi / ds

      filter = .false.

      ! spectrum is filtered only if grid-step of internal mesh is not equal to grid-step of external mesh
      IF (ds .ne. dh) filter = .true.

      ! vectors go from 0 to nyquist and THEN back again until dk
      kx = [[(i * dk(1), i = 0, npts(1)/2)], [(i * dk(1), i = npts(1)/2-1, 1, -1)]]
      ky = [[(j * dk(2), j = 0, npts(2)/2)], [(j * dk(2), j = npts(2)/2-1, 1, -1)]]
      kz = [[(k * dk(3), k = 0, npts(3)/2)], [(k * dk(3), k = npts(3)/2-1, 1, -1)]]

      ! compute part of power spectral density outside loop
      IF (acf .eq. 0) THEN
        num = 8._f_real * SQRT(pi**3) * GAMMA(hurst + 1.5_f_real) * sigma**2 * PRODUCT(cl) / GAMMA(hurst)
        fun => vk_psdf
      ELSEIF (acf .eq. 1) THEN
        num = sigma**2 * PRODUCT(cl) * SQRT(pi**3)
        fun => gs_psdf
      ENDIF

      ! each process generate its set of random numbers
      cptr = prng(seed, ls, le, npts)

      CALL C_F_POINTER(cptr, r, [(le(1) - ls(1) + 1) * (le(2) - ls(2) + 1) * (le(3) - ls(3) + 1)])

      harvest(ls(1):le(1), ls(2):le(2), ls(3):le(3)) => r

      ! compute spectrum
      DO k = ls(3), le(3)

        DO j = ls(2), le(2)

          DO i = ls(1), le(1)

            ! radial wavenumber
            kr = SQRT(kx(i)**2 + ky(j)**2 + kz(k)**2)

            ! fourth-order low-pass butterworth filter
            ! butter = 1._f_real / SQRT(1._f_real + (kr / kc)**(2 * 4))

            butter = 1._f_real

            IF (filter .and. (kr .gt. kc)) butter = 0._f_real

            ! now "kr" is the product "k * cl"
            kr = (kx(i) * cl(1))**2 + (ky(j) * cl(2))**2 + (kz(k) * cl(3))**2

            ! complete power spectral density and go to amplitude spectrum
            amp = SQRT(num / fun(kr, hurst))

            ! apply filter
            amp = amp * butter

            harvest(i, j, k) = harvest(i, j, k) * 2._f_real * pi

            ! combine amplitude and random phase
            spec(i, j, k) = cmplx(cos(harvest(i, j, k)) * amp, sin(harvest(i, j, k)) * amp, f_real)

          ENDDO
        ENDDO
      ENDDO

      NULLIFY(r, harvest)
      CALL free_mem(cptr)

#ifdef TIMING
      CALL watch_stop(tictoc)

      time(1) = tictoc

      CALL watch_start(tictoc)
#endif

      ! apply symmetry conditions
      CALL enforce_symmetry(ls, le, spec)

#ifdef TIMING
      CALL watch_stop(tictoc)

      time(2) = tictoc
#endif

    END SUBROUTINE compute_spectrum

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(f_real) FUNCTION vk_psdf(m, hurst)

      ! Purpose:
      ! To compute the denominator of Von Karman (or exponential, if hurst = 0.5) PSDFs.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_real), INTENT(IN) :: m                  !< quantity (wavenumber*cl)**2
      REAL(f_real), INTENT(IN) :: hurst              !< hurst exponent

      !-----------------------------------------------------------------------------------------------------------------------------

      vk_psdf = (1._f_real + m)**(hurst + 1.5_f_real)

    END FUNCTION vk_psdf

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(f_real) FUNCTION gs_psdf(m, void)

      ! Purpose:
      ! To compute the denominator of Gaussian PSDFs. Second argument is not used but is kept to allow generic procedure pointers.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_real), INTENT(IN) :: m                 !< quantity (wavenumber*cl)**2
      REAL(f_real), INTENT(IN) :: void              !< this variable is not used

      !-----------------------------------------------------------------------------------------------------------------------------

      gs_psdf = EXP(0.25_f_real * m)

    END FUNCTION gs_psdf

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE enforce_symmetry(fs, fe, spec)

      ! Purpose:
      ! To apply Hermitian symmetry on the PSDF.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(f_int),              DIMENSION(3),                                   INTENT(IN)    :: fs, fe        !< grid indices
      COMPLEX(f_real),             DIMENSION(fs(1):fe(1),fs(2):fe(2),fs(3):fe(3)), INTENT(INOUT) :: spec          !< spectrum
      INTEGER(f_int)                                                                             :: i
      INTEGER(f_int)                                                                             :: rank, ntasks
      INTEGER(f_int)                                                                             :: newcomm, ierr, color
      INTEGER(f_int),              DIMENSION(3)                                                  :: nyquist
      INTEGER(f_int), ALLOCATABLE, DIMENSION(:)                                                  :: new2world
      LOGICAL,                     DIMENSION(2)                                                  :: bool

      !-----------------------------------------------------------------------------------------------------------------------------

      ! index of nyquist wavenumber
      DO i = 1, 3
        nyquist(i) = npts(i) / 2 + 1
      ENDDO

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! create a new communicator containing only those processes crossing x-axis at 1 and npts(1)/2 + 1
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      color = 0

      bool(1) = fs(1) .eq. 1
      bool(2) = (fs(1) .le. nyquist(1)) .and. (fe(1) .ge. nyquist(1))

      IF (bool(1)) color = 1
      IF (bool(2)) color = 2

      IF (ALL(bool)) color = 3

      CALL mpi_comm_split(mpi_comm_world, color, world_rank, newcomm, ierr)

      CALL mpi_comm_rank(newcomm, rank, ierr)
      CALL mpi_comm_size(newcomm, ntasks, ierr)

      ALLOCATE(new2world(0:ntasks - 1))

      ! i-th rank in "newcomm" has its global rank contained in "new2world"
      CALL map_ranks(newcomm, mpi_comm_world, new2world)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! exchange data and enforce symmetry
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! if only one of conditions above is true, we can work on i=1 and i=nyquist(1) at the same time...
      IF ( (color .eq. 1) .or. (color .eq. 2) ) THEN

        IF (bool(1)) THEN
          i = 1
        ELSE
          i = nyquist(1)
        ENDIF

        CALL exchange_conjg

      ! ... otherwise work in two steps
      ELSEIF (color .eq. 3) THEN

        i = 1

        CALL exchange_conjg

        i = nyquist(1)

        CALL exchange_conjg

      ENDIF

      DEALLOCATE(new2world)

      CALL mpi_comm_free(newcomm, ierr)

      CONTAINS

      !-----------------------------------------------------------------------------------------------------------------------------

        SUBROUTINE exchange_conjg

          COMPLEX(f_real), DIMENSION(fs(2):fe(2),fs(3):fe(3)) :: buffer, sendbuf
          INTEGER(f_int)                                      :: j, k, p, l, c, cy, cz
          INTEGER(f_int)                                      :: js, je, ks, ke
          INTEGER(f_int)                                      :: j0, j1, k0, k1
          INTEGER(f_int)                                      :: ierr
          INTEGER(f_int), DIMENSION(2)                        :: rsizes, rsubsizes, rstarts
          INTEGER(f_int), DIMENSION(3)                        :: ssizes, ssubsizes, sstarts
          INTEGER(f_int), DIMENSION(0:ntasks-1)               :: sendcounts, recvcounts
          INTEGER(f_int), DIMENSION(0:ntasks-1)               :: sdispls, rdispls
          INTEGER(f_int), DIMENSION(0:ntasks-1)               :: sendtypes, recvtypes
          LOGICAL                                             :: bool

          !-------------------------------------------------------------------------------------------------------------------------

          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
          ! apply conditions of figure 4 for alpha, psi, xsi and eta
          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

          bool = (i .eq. 1) .and. (fs(2) .eq. 1) .and. (fs(3) .eq. 1)

          IF (bool) spec(i, 1, 1) = 0._f_real                                                         !< alpha must be zero for zero-mean field

          bool = (i .eq. nyquist(1)) .and. (fs(2) .eq. 1) .and. (fs(3) .eq. 1)

          IF (bool) spec(i, 1, 1) = REAL(spec(i, 1, 1), f_real)                                       !< alpha must be REAL

          bool = (fs(2) .le. nyquist(2)) .and. (fe(2) .ge. nyquist(2)) .and. (fs(3) .le. nyquist(3)) .and. (fe(3) .ge. nyquist(3))

          IF (bool) spec(i, nyquist(2), nyquist(3)) = REAL(spec(i, nyquist(2), nyquist(3)), f_real)   !< psi must be REAL

          bool = (fs(2) .eq. 1) .and. (fs(3) .le. nyquist(3)) .and. (fe(3) .ge. nyquist(3))

          IF (bool) spec(i, 1, nyquist(3)) = REAL(spec(i, 1, nyquist(3)), f_real)                     !< xsi must be REAL

          bool = (fs(2) .le. nyquist(2)) .and. (fe(2) .ge. nyquist(2)) .and. (fs(3) .eq. 1)

          IF (bool) spec(i, nyquist(2), 1) = REAL(spec(i, nyquist(2), 1), f_real)                     !< eta must be REAL

          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
          ! apply conditions of figure 4 at remaining regions. exchange data amongst processes. beta- and delta- regions are treated
          ! separately.
          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

          sendcounts = 0
          recvcounts = 0

          sdispls = 0
          rdispls = 0

          ssizes = [fe(1) - fs(1) + 1, fe(2) - fs(2) + 1, fe(3) - fs(3) + 1]

          rsizes = [fe(2) - fs(2) + 1, fe(3) - fs(3) + 1]

          ! copy relevant part of "spec". This is not strictly necessary but may improve cache access.
          DO k = fs(3), fe(3)
            DO j = fs(2), fe(2)
              buffer(j, k) = spec(i, j, k)
            ENDDO
          ENDDO

          ! exchange data in the whole region, except beta- and delta-region
          DO p = 0, ntasks - 1

            l = new2world(p)                                                              !< global index of p-th process

            sendtypes(p) = cplx_type
            recvtypes(p) = cplx_type

            ! project points of process "p" into conjugated field.
            js = npts(2) - ge(2, l) + 2
            je = npts(2) - gs(2, l) + 2

            ks = npts(3) - ge(3, l) + 2
            ke = npts(3) - gs(3, l) + 2

            j0 = MAX(fs(2), js)
            j1 = MIN(fe(2), je)

            k0 = MAX(fs(3), ks)
            k1 = MIN(fe(3), ke)

            ! determine data to be sent/received to/from process "p"
            IF ( (k0 .le. k1) .and. (j0 .le. j1) ) THEN

              cy = j1 + j0
              cz = k1 + k0

              DO k = k0, k1
                DO j = j0, j1
                  sendbuf(cy - j, cz - k) = buffer(j, k)
                ENDDO
              ENDDO

              rsubsizes = [j1 - j0 + 1, k1 - k0 + 1]
              rstarts   = [j0 - fs(2), k0 - fs(3)]

              CALL mpi_type_create_subarray(2, rsizes, rsubsizes, rstarts, mpi_order_fortran, cplx_type, sendtypes(p), ierr)
              CALL mpi_type_commit(sendtypes(p), ierr)

              CALL mpi_type_create_subarray(2, rsizes, rsubsizes, rstarts, mpi_order_fortran, cplx_type, recvtypes(p), ierr)
              CALL mpi_type_commit(recvtypes(p), ierr)

              sendcounts(p) = 1
              recvcounts(p) = 1

            ENDIF

          ENDDO

          CALL mpi_alltoallw(sendbuf, sendcounts, sdispls, sendtypes, buffer, recvcounts, rdispls, recvtypes, newcomm, ierr)

          ! free resources
          DO p = 0, ntasks - 1
            IF (sendcounts(p) .eq. 1) THEN
              CALL mpi_type_free(sendtypes(p), ierr)
              CALL mpi_type_free(recvtypes(p), ierr)
            ENDIF
          ENDDO

          sendcounts = 0
          recvcounts = 0

          ! now exchange data in the beta- and delta-region
          DO p = 0, ntasks - 1

            l = new2world(p)                                                              !< global index of p-th process

            sendtypes(p) = cplx_type
            recvtypes(p) = cplx_type

            ! project points of process "p" into conjugated field.
            IF (gs(3, l) .eq. 1) THEN

              js = npts(2) - ge(2, l) + 2
              je = npts(2) - gs(2, l) + 2

              ks = 1
              ke = 1

              j0 = MAX(fs(2), js)
              j1 = MIN(fe(2), je)

              k0 = MAX(fs(3), ks)
              k1 = MIN(fe(3), ke)

              ! determine data to be sent/received to/from process "p"
              IF ( (k0 .le. k1) .and. (j0 .le. j1) ) THEN

                rsubsizes = [j1 - j0 + 1, k1 - k0 + 1]
                rstarts   = [j0 - fs(2), k0 - fs(3)]

                cy = j1 + j0
                cz = k1 + k0

                DO k = k0, k1
                  DO j = j0, j1
                    sendbuf(cy - j, cz - k) = buffer(j, k)
                  ENDDO
                ENDDO

                CALL mpi_type_create_subarray(2, rsizes, rsubsizes, rstarts, mpi_order_fortran, cplx_type, sendtypes(p), ierr)
                CALL mpi_type_commit(sendtypes(p), ierr)

                CALL mpi_type_create_subarray(2, rsizes, rsubsizes, rstarts, mpi_order_fortran, cplx_type, recvtypes(p), ierr)
                CALL mpi_type_commit(recvtypes(p), ierr)

                sendcounts(p) = 1
                recvcounts(p) = 1

              ENDIF

            ENDIF

            ! project points of process "p" into conjugated field.
            IF (gs(2, l) .eq. 1) THEN

              js = 1
              je = 1

              ks = npts(3) - ge(3, l) + 2
              ke = npts(3) - gs(3, l) + 2

              j0 = MAX(fs(2), js)
              j1 = MIN(fe(2), je)

              k0 = MAX(fs(3), ks)
              k1 = MIN(fe(3), ke)

              ! determine data to be sent/received to/from process "p"
              IF ( (k0 .le. k1) .and. (j0 .le. j1) ) THEN

                rsubsizes = [j1 - j0 + 1, k1 - k0 + 1]
                rstarts   = [j0 - fs(2), k0 - fs(3)]

                cy = j1 + j0
                cz = k1 + k0

                DO k = k0, k1
                  DO j = j0, j1
                    sendbuf(cy - j, cz - k) = buffer(j, k)
                  ENDDO
                ENDDO

                CALL mpi_type_create_subarray(2, rsizes, rsubsizes, rstarts, mpi_order_fortran, cplx_type, sendtypes(p), ierr)
                CALL mpi_type_commit(sendtypes(p), ierr)

                CALL mpi_type_create_subarray(2, rsizes, rsubsizes, rstarts, mpi_order_fortran, cplx_type, recvtypes(p), ierr)
                CALL mpi_type_commit(recvtypes(p), ierr)

                sendcounts(p) = 1
                recvcounts(p) = 1

              ENDIF

            ENDIF

          ENDDO

          CALL mpi_alltoallw(sendbuf, sendcounts, sdispls, sendtypes, buffer, recvcounts, rdispls, recvtypes, newcomm, ierr)

          ! free resources
          DO p = 0, ntasks - 1
            IF (sendcounts(p) .eq. 1) THEN
              CALL mpi_type_free(sendtypes(p), ierr)
              CALL mpi_type_free(recvtypes(p), ierr)
            ENDIF
          ENDDO

          sendcounts = 0
          recvcounts = 0

          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
          ! update spectrum with conjugate values
          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

          ! delta*, gamma*, phi*, epsilon*
          DO k = MAX(nyquist(3) + 1, fs(3)), fe(3)
            DO j = fs(2), fe(2)
              spec(i, j, k) = CONJG(buffer(j, k))
            ENDDO
          ENDDO

          ! beta*
          DO k = MAX(1, fs(3)), 1
            DO j = MAX(nyquist(2) + 1, fs(2)), fe(2)
              spec(i, j, k) = CONJG(buffer(j, k))
            ENDDO
          ENDDO

          ! theta*
          DO k = MAX(nyquist(3), fs(3)), MIN(nyquist(3), fe(3))
            DO j = MAX(nyquist(2) + 1, fs(2)), fe(2)
              spec(i, j, k) = CONJG(buffer(j, k))
            ENDDO
          ENDDO

        END SUBROUTINE exchange_conjg

    END SUBROUTINE enforce_symmetry

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE exchange_halo(cartopo, buffer)

      ! Purpose:
      ! To exchange points inside halos between neighboring processes. An halo is represented by the first column/row (along x/y,
      ! respectively) of input/output array.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(f_int),                   INTENT(IN)    :: cartopo                         !< handle to cartesian topology
      REAL(f_real),   DIMENSION(:,:,:), INTENT(INOUT) :: buffer                         !< random field with halo
      INTEGER(f_int)                                  :: i
      INTEGER(f_int)                                  :: nx, ny, nz
      INTEGER(f_int)                                  :: ierr, neg, pos
      INTEGER(f_int)                                  :: from_neg, to_pos
      INTEGER(f_int), DIMENSION(3)                    :: sizes, subsizes, rstarts, sstarts

      !-----------------------------------------------------------------------------------------------------------------------------

      ! points along x and y without extra row/column
      nx = SIZE(buffer, 1) - 1
      ny = SIZE(buffer, 2) - 1
      nz = SIZE(buffer, 3) - 1

      ! total SIZE of "buffer" with extra row/column for halo exchange
      sizes = [nx + 1, ny + 1, nz + 1]

      ! exchange in the x-, y-, z-direction
      DO i = 0, 2

        ! determine neighbors for a positive unitary shift along current direction
        CALL mpi_cart_shift(cartopo, i, 1, neg, pos, ierr)

        IF (i .eq. 0) THEN
          subsizes = [1, ny, nz]            !< shape of halo
          rstarts  = [0, 1, 1]              !< initial address for data to be received
          sstarts  = [nx, 1, 1]             !< initial address for data to be sent
        ELSEIF (i .eq. 1) THEN
          subsizes = [nx + 1, 1, nz]        !< shape of halo
          rstarts  = [0, 0, 1]              !< initial address for data to be received
          sstarts  = [0, ny, 1]             !< initial address for data to be sent
        ELSE
          subsizes = [nx + 1, ny + 1, 1]    !< shape of halo
          rstarts  = [0, 0, 0]              !< initial address for data to be received
          sstarts  = [0, 0, nz]             !< initial address for data to be sent
        ENDIF

        ! data to be received
        CALL mpi_type_create_subarray(3, sizes, subsizes, rstarts, mpi_order_fortran, real_type, from_neg, ierr)
        CALL mpi_type_commit(from_neg, ierr)

        ! data to be sent
        CALL mpi_type_create_subarray(3, sizes, subsizes, sstarts, mpi_order_fortran, real_type, to_pos, ierr)
        CALL mpi_type_commit(to_pos, ierr)

        ! exchange halo data with neighbors. since we operate on first and last columns, send "buffer" is disjoint from recv "buffer"
        CALL mpi_sendrecv(buffer, 1, to_pos, pos, 0, buffer, 1, from_neg, neg, 0, mpi_comm_world, mpi_status_ignore, ierr)

        ! realease resources
        CALL mpi_type_free(from_neg, ierr)
        CALL mpi_type_free(to_pos, ierr)

      ENDDO

    END SUBROUTINE exchange_halo

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE interpolate(pmax, xyz, delta, v)

      ! Purpose:
      ! To interpolate values from internal (FFT) grid into external one. Possible interpolation schemes are nearest neighbor and
      ! bilinear. Note that bilinear interpolation has a smoothing effect stronger than nearest neighbor.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version (nearest neighbor interp.)
      !

      INTEGER(f_int),                   INTENT(IN)  :: pmax                         !< highest index to consider
      REAL(f_real),   DIMENSION(:,:),   INTENT(IN)  :: xyz                          !< external mesh nodes projected onto internal mesh
      REAL(f_real),   DIMENSION(:,:,:), INTENT(IN)  :: delta                        !< input field
      REAL(f_real),   DIMENSION(:),     INTENT(OUT) :: v                            !< interpolated values
      INTEGER(f_int)                                :: p
      INTEGER(f_int)                                :: i0, j0, k0
      REAL(f_real)                                  :: i, j, k
      REAL(f_real)                                  :: px, py, pz, ipx, ipy
      REAL(f_real)                                  :: a, b
      REAL(f_real),   DIMENSION(2,2)                :: f

      !-----------------------------------------------------------------------------------------------------------------------------

      ! loop over points
      DO p = 1, pmax

        i = xyz(1, p)
        j = xyz(2, p)
        k = xyz(3, p)

        ! nearest neighbor
        ! {
        i0 = NINT(i)
        j0 = NINT(j)
        k0 = NINT(k)

        v(p) = delta(i0, j0, k0)
        ! }

        ! bilinear
        ! {
        ! i0 = FLOOR(i)
        ! j0 = FLOOR(j)
        ! k0 = FLOOR(k)
        !
        ! f(:, 1) = delta(i0:i0 + 1, j0, k0)
        ! f(:, 2) = delta(i0:i0 + 1, j0 + 1, k0)
        !
        ! px = i - i0
        ! py = j - j0
        !
        ! ipx = (1._f_real - px)
        ! ipy = (1._f_real - py)
        !
        ! ! bilinear interpolation at level 1
        ! a = f(1, 1) * ipx * ipy + f(2, 1) * px * ipy + f(1, 2) * ipx * py + f(2, 2) * px * py
        !
        ! f(:, 1) = delta(i0:i0 + 1, j0, k0 + 1)
        ! f(:, 2) = delta(i0:i0 + 1, j0 + 1, k0 + 1)
        !
        ! ! bilinear interpolation at level 2
        ! b = f(1, 1) * ipx * ipy + f(2, 1) * px * ipy + f(1, 2) * ipx * py + f(2, 2) * px * py
        !
        ! pz = k - k0
        !
        ! ! linear interpolated between level 1 and 2
        ! v(p) = a * (1._f_real - pz) + b * pz
        ! }

      ENDDO

    END SUBROUTINE interpolate

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE coords2index(npts, dims, coords, fs, fe)

      ! Purpose:
      ! To return first/last grid index for the calling process in a given topology.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(f_int), DIMENSION(3), INTENT(IN)  :: npts                        !< total number of grid points along each axis
      INTEGER(f_int), DIMENSION(3), INTENT(IN)  :: dims                        !< number of processes along each axis
      INTEGER(f_int), DIMENSION(3), INTENT(IN)  :: coords                      !< calling process coordinates
      INTEGER(f_int), DIMENSION(3), INTENT(OUT) :: fs, fe                      !< resulting first/last index
      INTEGER(f_int)                            :: i

      !-----------------------------------------------------------------------------------------------------------------------------

      DO i = 1, 3
        fs(i) = 1 + int( REAL(npts(i)) / REAL(dims(i)) * REAL(coords(i)) )
        fe(i) = int( REAL(npts(i)) / REAL(dims(i)) * REAL(coords(i) + 1) )
      ENDDO

    END SUBROUTINE coords2index

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE best_config(dims)

      ! Purpose:
      ! To return an optimal grid of process according to some specific cost function.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(f_int),              DIMENSION(3), INTENT(OUT) :: dims
      INTEGER(f_int)                                         :: i, j, k, l, c
      INTEGER(f_int)                                         :: n2, n3
      INTEGER(f_int)                                         :: a, b
      INTEGER(f_int),              DIMENSION(3)              :: v1, v2
      INTEGER(f_int), ALLOCATABLE, DIMENSION(:,:)            :: fact3, fact2
      INTEGER(f_int), ALLOCATABLE, DIMENSION(:,:)            :: list
      REAL(f_real)                                           :: val, minimum

      !-----------------------------------------------------------------------------------------------------------------------------

      minimum = HUGE(1._f_real)

      c = 0

      ! factorise the number of available processes (return pairs)
      CALL factorization(world_size, fact3)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! factorise each pair member (thus resulting in a triplet of numbers), evaluate cost function and eventually store triplet
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! number of pairs found
      n3 = SIZE(fact3, 2)

      ! allocate a list where to store all the factorized triplets. Factor of 10 seems to be large enough for 200k processes
      ALLOCATE(list(3, n3 * 10))

      list(:,:) = 0

      ! loop over factorized processes
      DO l = 1, n3

        ! loop over pair members
        DO k = 1, 2

          IF (k .eq. 1) THEN
            a = fact3(1, l)
            b = fact3(2, l)
          ELSE
            b = fact3(1, l)
            a = fact3(2, l)
          ENDIF

          ! factorization of a number returns a new pair
          CALL factorization(b, fact2)

          n2 = SIZE(fact2, 2)

          ! loop over new pairs
          DO j = 1, n2

            ! build candidate triplet
            v1 = [a, fact2(:, j)]

            ! skip to next pair if current triplet already analysed ("v1" is already in "list")
            IF (match(v1, c, list) .eqv. .true.) CYCLE

            c = c + 1

            ! triplet is new: add it to the list
            list(:, c) = v1

            ! evaluate cost function for all three possible arrangements (permutations) of the triplet
            DO i = 0, 2

              v1 = CSHIFT(v1, 1)

              ! evaluate cost function
              CALL cost_fun(v1, val)

              IF (val .lt. minimum) THEN
                dims    = v1
                minimum = val
              ENDIF

              v2 = [v1(1), v1(3), v1(2)]

              CALL cost_fun(v2, val)

              IF (val .lt. minimum) THEN
                dims    = v2
                minimum = val
              ENDIF

            ENDDO
            ! end permutations

          ENDDO
          ! end loop over factor pairs for "a/b"

        ENDDO
        ! end loop over "a/b"

      ENDDO
      ! end loop over factor pairs for "world_size"

      DEALLOCATE(fact3, fact2, list)

      !-----------------------------------------------------------------------------------------------------------------------------
      !-----------------------------------------------------------------------------------------------------------------------------

      CONTAINS

      LOGICAL FUNCTION match(vec, imax, list)

        INTEGER(f_int), DIMENSION(3),   INTENT(IN) :: vec
        INTEGER(f_int),                 INTENT(IN) :: imax
        INTEGER(f_int), DIMENSION(:,:), INTENT(IN) :: list
        INTEGER(f_int)                             :: i

        !---------------------------------------------------------------------------------------------------------------------------

        match = .false.

        DO i = 1, imax
          match = ANY(v1(1) .eq. list(:, i)) .and. ANY(v1(2) .eq. list(:, i)) .and. ANY(v1(3) .eq. list(:, i))
          IF (match .eqv. .true.) exit
        ENDDO

      END FUNCTION match

    END SUBROUTINE best_config

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    ! SUBROUTINE cost_fun(dims, niter)
    !
    !   INTEGER(f_int),              DIMENSION(3),             INTENT(IN)  :: dims
    !   INTEGER(f_int),                                        INTENT(OUT) :: niter
    !   INTEGER(f_int)                                                     :: i                            !< counter
    !   INTEGER(f_int)                                                     :: ierr, cartopo
    !   INTEGER(f_int),              DIMENSION(3)                          :: ls, le
    !   INTEGER(f_int),              DIMENSION(3)                          :: coords
    !   INTEGER(f_int),              DIMENSION(0:world_size-1)             :: completed, isbusy
    !   INTEGER(f_int), ALLOCATABLE, DIMENSION(:)                          :: buffer
    !
    !   !-----------------------------------------------------------------------------------------------------------------------------
    !
    !   ! create topology
    !   CALL mpi_cart_create(mpi_comm_world, ndims, dims, isperiodic, reorder, cartopo, ierr)
    !
    !   ! DO i = 0, world_size - 1
    !   !   CALL mpi_cart_coords(cartopo, i, ndims, coords, ierr)
    !   !   CALL coords2index(npts, dims, coords, ls, le)
    !   !   gs(:, i) = ls
    !   !   ge(:, i) = le
    !   ! ENDDO
    !
    !   ! return process coordinates in current topology
    !   CALL mpi_cart_coords(cartopo, world_rank, ndims, coords, ierr)
    !
    !   ! return first/last-index ("ls"/"le") along each direction for calling process. note: first point has always index equal to 1.
    !   CALL coords2index(npts, dims, coords, ls, le)
    !
    !   gs(:, world_rank) = ls
    !   ge(:, world_rank) = le
    !
    !   ! make ALL processes aware of global indices along each axis
    !   CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, gs, 3, mpi_integer, mpi_comm_world, ierr)
    !   CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, ge, 3, mpi_integer, mpi_comm_world, ierr)
    !
    !   ! initialise list with processes having their points assigned
    !   completed(:) = 0
    !
    !   niter = 0
    !
    !   ! CYCLE until ALL processes are completed
    !   DO while(ANY(completed .eq. 0))
    !
    !     ! reset "busy process" flag
    !     isbusy(:) = 0
    !
    !     ! create as many communicators as possible
    !     DO i = 0, world_size - 1
    !
    !       IF (completed(i) .eq. 0) THEN
    !
    !         ! produce a tentative list of processes that may join the new communicator
    !         CALL find_buddies(i, buffer)
    !
    !         ! build new communicator only IF ALL involved processes are not yet part of another communicator
    !         IF (ALL(isbusy(buffer) .eq. 0)) THEN
    !
    !           ! set processes belonging to new communicator as busy
    !           isbusy(buffer) = 1
    !
    !           ! set i-th process as completed
    !           completed(i) = 1
    !
    !         ENDIF
    !
    !         DEALLOCATE(buffer)
    !
    !       ENDIF
    !
    !     ENDDO
    !
    !     niter = niter + 1
    !
    !   ENDDO
    !
    !   CALL mpi_comm_free(cartopo, ierr)
    !
    ! END SUBROUTINE cost_fun

    SUBROUTINE cost_fun(dims, measure)

      ! Purpose:
      ! To return a cost estimate based on the size of a grid block assigned to each process, where the cost is represented by differences
      ! along each side (the goal is to having blocks in shape of cubes).
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(f_int), DIMENSION(3), INTENT(IN)  :: dims           !< processes in each direction
      REAL(f_real),                 INTENT(OUT) :: measure        !< cost
      REAL(f_real),   DIMENSION(3)              :: side           !< grid nodes in a block along each direction

      !-----------------------------------------------------------------------------------------------------------------------------

      side = REAL(npts, f_real) / REAL(dims, f_real)

      measure = abs(side(1) - side(2)) + abs(side(1) - side(3)) + abs(side(2) - side(3))

    END SUBROUTINE cost_fun

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE factorization(n, factors)

      ! Purpose:
      ! To compute the factorization of an integer number based on the trial division method. For example, for "n=12" the output looks
      ! like "[1 12; 2 6; 3 4]"
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(f_int),                              INTENT(IN)  :: n
      INTEGER(f_int), ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: factors
      INTEGER(f_int)                                           :: i, c, s
      INTEGER(f_int)                                           :: x
      INTEGER(f_int), ALLOCATABLE, DIMENSION(:,:)              :: buffer

      !-----------------------------------------------------------------------------------------------------------------------------

      ! max possible number of factors
      s = FLOOR(SQRT(REAL(n, f_real)))

      ALLOCATE(buffer(2, s))

      buffer(:,:) = 0

      ! test factors
      DO i = 1, s

        x = n / i

        IF (MOD(n, i) .eq. 0) THEN
          buffer(1, i) = i                          !< add factor ...
          buffer(2, i) = x                          !< ... and its companion
        ENDIF

      ENDDO

      ! actual factors found
      i = COUNT(buffer(1, :) .ne. 0)

      ALLOCATE(factors(2, i))

      ! copy factors to output array
      c = 0
      DO i = 1, s
        IF (buffer(1, i) .ne. 0) THEN
          c = c + 1
          factors(:, c) = buffer(:, i)
        ENDIF
      ENDDO

      DEALLOCATE(buffer)

    END SUBROUTINE factorization

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE m_scarflib_fim
