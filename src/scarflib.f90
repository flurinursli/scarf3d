MODULE m_scarflib

  ! Purpose:
  !   To provide an interface for subroutines at the core of the SCARF3D library
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   04/05/20                  original version
  !

  !USE, NON_INTRINSIC :: mpi                     !< this is included in common???
  USE, NON_INTRINSIC :: m_scarflib_common
  USE, NON_INTRINSIC :: m_scarflib_fim
  USE, NON_INTRINSIC :: m_scarflib_srm

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! declare all private, including also public variables and subroutines accessed from other modules via host-association
  PRIVATE

  ! these are the only bits and pieces accessible by an external calling program
  PUBLIC :: scarf_initialize, scarf_execute, scarf_finalize, scarf_io

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTERFACE scarf_initialize
    MODULE PROCEDURE initialize_structured, initialize_unstructured
  END INTERFACE

  INTERFACE scarf_execute
    MODULE PROCEDURE execute_structured, execute_unstructured
  END INTERFACE

  INTERFACE scarf_io
    MODULE PROCEDURE scarf_io_slice, scarf_io_one
  END INTERFACE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! derived datatype containing input parameters and output statistics
  TYPE scarf_obj

    INTEGER(f_int)                                       :: acf
    INTEGER(f_int)                                       :: rescale
    INTEGER(f_int)                                       :: pad
    INTEGER(f_int)                                       :: method
    INTEGER(f_int),                       DIMENSION(3)   :: fs, fe
    REAL(f_real)                                         :: ds, dh
    REAL(f_real)                                         :: sigma, hurst
    REAL(f_real)                                         :: mute, taper
    REAL(f_real),                         DIMENSION(3)   :: cl
    REAL(f_real),                         DIMENSION(3)   :: nc, fc
    REAL(f_real),                POINTER, DIMENSION(:)   :: x => NULL(), y => NULL(), z => NULL()
    REAL(f_real),   ALLOCATABLE,          DIMENSION(:)   :: stats
    REAL(f_real),   ALLOCATABLE,          DIMENSION(:,:) :: poi

  END TYPE scarf_obj

  TYPE(scarf_obj) :: obj

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE initialize_structured(fs, fe, ds, acf, cl, sigma, method, hurst, dh, poi, mute, taper, rescale, pad, nc, fc)

      ! Purpose:
      !   To setup most of the parameters needed to compute random fields: structured mesh version. Parameters concerning the external
      !   grid refer to calling process view.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(f_int), DIMENSION(3),             INTENT(IN) :: fs, fe         !< first/last index of external grid
      REAL(f_real),                             INTENT(IN) :: ds             !< grid-step external grid (controls max resolvable wavenumber)
      INTEGER(f_int),                           INTENT(IN) :: acf            !< autocorrelation function (0=von karman, 1=gaussian)
      REAL(f_real),   DIMENSION(3),             INTENT(IN) :: cl             !< correlation length
      REAL(f_real),                             INTENT(IN) :: sigma          !< standard deviation
      INTEGER(f_int),                 OPTIONAL, INTENT(IN) :: method         !< algorithm of choice (empty=0=fim, 1=srm)
      REAL(f_real),                   OPTIONAL, INTENT(IN) :: hurst          !< hurst exponent (needed only for acf=0)
      REAL(f_real),                   OPTIONAL, INTENT(IN) :: dh             !< min grid-step external grid (controls max absolute wavenumber)
      REAL(f_real),   DIMENSION(:,:), OPTIONAL, INTENT(IN) :: poi            !< points where taper/muting should be applied
      REAL(f_real),                   OPTIONAL, INTENT(IN) :: mute, taper    !< radius for taper/muting
      INTEGER(f_int),                 OPTIONAL, INTENT(IN) :: rescale        !< rescale discrete std.dev. to continuous one
      INTEGER(f_int),                 OPTIONAL, INTENT(IN) :: pad            !< pad internal grid (fim only, empty=0=no, 1=yes)
      REAL(f_real),   DIMENSION(3),   OPTIONAL, INTENT(IN) :: nc, fc         !< min/max extent external grid (empty = use fs,fe)

      !-----------------------------------------------------------------------------------------------------------------------------

      obj%fs  = fs
      obj%fe  = fe
      obj%ds  = ds

      ! by default minimum grid-step of external grid equals current grid-step: this implies Kc = Kmax
      obj%dh    = ds
      obj%acf   = acf
      obj%cl    = cl
      obj%sigma = sigma

      obj%method = 0                       !< default to FIM algorithm

      obj%hurst = 0._f_real

      obj%mute    = -999._f_real
      obj%taper   = -999._f_real
      obj%rescale = 0
      obj%pad     = 0

      obj%nc      = (fs - 1) * ds          !< by default, "nc" and "fc" are given by "fs" and "fe"
      obj%fc      = (fe - 1) * ds

      ALLOCATE(obj%poi(0,0))               !< by default, there are no points where tapering/muting should be applied

      ! override default values
      IF (PRESENT(method))  obj%method  = method
      IF (PRESENT(hurst))   obj%hurst   = hurst
      IF (PRESENT(dh))      obj%dh      = dh
      IF (PRESENT(mute))    obj%mute    = mute
      IF (PRESENT(taper))   obj%taper   = taper
      IF (PRESENT(rescale)) obj%rescale = rescale
      IF (PRESENT(pad))     obj%pad     = pad
      IF (PRESENT(poi))     obj%poi     = poi
      IF (PRESENT(nc))      obj%nc      = nc
      IF (PRESENT(fc))      obj%fc      = fc

      IF (obj%method .eq. 0) THEN
        ALLOCATE(obj%stats(8))
      ELSEIF (obj%method .eq. 1) THEN
        ALLOCATE(obj%stats(6))
      ENDIF

    END SUBROUTINE initialize_structured

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE initialize_unstructured(x, y, z, dh, acf, cl, sigma, method, hurst, ds, poi, mute, taper, rescale, pad, nc, fc)

      ! Purpose:
      !   To setup most of the parameters needed to compute random fields: unstructured mesh version. Parameters concerning the external
      !   grid refer to calling process view.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_real),   DIMENSION(:),   TARGET,           INTENT(IN) :: x, y, z            !< external grid nodes
      REAL(f_real),                                     INTENT(IN) :: dh                 !< min grid-step external grid (controls max absolute wavenumber)
      INTEGER(f_int),                                   INTENT(IN) :: acf                !< autocorrelation function (0=von karman, 1=gaussian)
      REAL(f_real),   DIMENSION(3),                     INTENT(IN) :: cl                 !< correlation length
      REAL(f_real),                                     INTENT(IN) :: sigma              !< continuous tandard deviation
      INTEGER(f_int),                         OPTIONAL, INTENT(IN) :: method             !< algorithm of choice (empty=0=fim, 1=srm)
      REAL(f_real),                           OPTIONAL, INTENT(IN) :: hurst              !< hurst exponent (needed only for acf=0)
      REAL(f_real),                           OPTIONAL, INTENT(IN) :: ds                 !< grid-step external grid (controls max resolvable wavenumber)
      REAL(f_real),   DIMENSION(:,:),         OPTIONAL, INTENT(IN) :: poi                !< points where taper/muting should be applied
      REAL(f_real),                           OPTIONAL, INTENT(IN) :: mute, taper        !< radius for taper/muting
      INTEGER(f_int),                         OPTIONAL, INTENT(IN) :: rescale            !< rescale discrete std.dev. to continuous one
      INTEGER(f_int),                         OPTIONAL, INTENT(IN) :: pad                !< pad internal grid (FIM only, empty=0=no, 1=yes)
      REAL(f_real),   DIMENSION(3),           OPTIONAL, INTENT(IN) :: nc, fc             !< near corner, far corner

      !-----------------------------------------------------------------------------------------------------------------------------

      obj%x     => x
      obj%y     => y
      obj%z     => z

      obj%dh    = dh
      obj%acf   = acf
      obj%cl    = cl
      obj%sigma = sigma

      ! by default current grid-step of external grid equals min grid-step: this implies Kc = Kmax
      obj%ds     = dh

      obj%method = 0                        !< default to FIM algorithm

      obj%hurst = 0._f_real

      obj%mute    = -999._f_real
      obj%taper   = -999._f_real
      obj%rescale = 0
      obj%pad     = 0

      ! by default, "nc" and "fc" are determined by input grid nodes
      obj%nc      = [MINVAL(x, dim = 1), MINVAL(y, dim = 1), MINVAL(z, dim = 1)]
      obj%fc      = [MAXVAL(x, dim = 1), MAXVAL(y, dim = 1), MAXVAL(z, dim = 1)]

      ALLOCATE(obj%poi(0,0))                 !< by default, there are no points where tapering/muting should be applied

      ! override default values
      IF (PRESENT(method))  obj%method  = method
      IF (PRESENT(hurst))   obj%hurst   = hurst
      IF (PRESENT(ds))      obj%ds      = ds
      IF (PRESENT(mute))    obj%mute    = mute
      IF (PRESENT(taper))   obj%taper   = taper
      IF (PRESENT(rescale)) obj%rescale = rescale
      IF (PRESENT(pad))     obj%pad     = pad
      IF (PRESENT(poi))     obj%poi     = poi
      IF (PRESENT(nc))      obj%nc      = nc
      IF (PRESENT(fc))      obj%fc      = fc

      IF (obj%method .eq. 0) THEN
        ALLOCATE(obj%stats(8))
      ELSEIF (obj%method .eq. 1) THEN
        ALLOCATE(obj%stats(6))
      ENDIF

    END SUBROUTINE initialize_unstructured

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE execute_unstructured(seed, field, stats)

      ! Purpose:
      !   To launch random field calculations for unstructured meshes
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(f_int),               INTENT(IN)  :: seed         !< seed number
      REAL(f_real),   DIMENSION(:), INTENT(OUT) :: field        !< random field
      REAL(f_real),   DIMENSION(8), INTENT(OUT) :: stats        !< vector with accuracy (1:2), std.dev. (3), mean(4) and timing flags (5:6/8)

      !-----------------------------------------------------------------------------------------------------------------------------

      stats(:) = 0._f_real

      IF (obj%method .eq. 0) THEN

        CALL scarf3d_fim(obj%nc, obj%fc, obj%ds, obj%x, obj%y, obj%z, obj%dh, obj%acf, obj%cl, obj%sigma, obj%hurst, seed,   &
                         obj%poi, obj%mute, obj%taper, obj%rescale, obj%pad, field, obj%stats)

        stats = obj%stats

      ELSEIF (obj%method .eq. 1) THEN

        CALL scarf3d_srm(obj%nc, obj%fc, obj%ds, obj%x, obj%y, obj%z, obj%dh, obj%acf, obj%cl, obj%sigma, obj%hurst, seed,  &
                          obj%poi, obj%mute, obj%taper, obj%rescale, field, obj%stats)

        stats(1:6) = obj%stats(:)

      ENDIF

    END SUBROUTINE execute_unstructured

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE execute_structured(seed, field, stats)

      ! Purpose:
      !   To launch random field calculations for structured meshes
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(f_int),                   INTENT(IN)  :: seed         !< seed number
      REAL(f_real),   DIMENSION(:,:,:), INTENT(OUT) :: field        !< random field
      REAL(f_real),   DIMENSION(8),     INTENT(OUT) :: stats        !< vector with accuracy (1:2), std.dev. (3), mean(4) and timing flags (5:6/8)

      !-----------------------------------------------------------------------------------------------------------------------------

      stats(:) = 0._f_real

      IF (obj%method .eq. 0) THEN

        CALL scarf3d_fim(obj%nc, obj%fc, obj%ds, obj%fs, obj%fe, obj%dh, obj%acf, obj%cl, obj%sigma, obj%hurst, seed, obj%poi, &
                         obj%mute, obj%taper, obj%rescale, obj%pad, field, obj%stats)

        stats = obj%stats

      ELSEIF (obj%method .eq. 1) THEN

        CALL scarf3d_srm(obj%nc, obj%fc, obj%ds, obj%fs, obj%fe, obj%dh, obj%acf, obj%cl, obj%sigma, obj%hurst, seed, obj%poi, &
                          obj%mute, obj%taper, obj%rescale, field, obj%stats)

        stats(1:6) = obj%stats(:)

      ENDIF

    END SUBROUTINE execute_structured

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE scarf_finalize()

      ! Purpose:
      !   To release resources allocated by the library.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (ALLOCATED(obj%poi))   DEALLOCATE(obj%poi)
      IF (ALLOCATED(obj%stats)) DEALLOCATE(obj%stats)

      NULLIFY(obj%x, obj%y, obj%z)

    END SUBROUTINE scarf_finalize

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE scarf_io_slice(n, axis, plane, field, filename)

      ! Purpose:
      !   To write to disk slices of random field perpendicular to cartesian axes. Works for structured meshes only.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(f_int),   DIMENSION(3),     INTENT(IN) :: n                  !< grid points in the global mesh
      CHARACTER(LEN=1),                   INTENT(IN) :: axis               !< axis perpendicular to the slice
      INTEGER(f_int),                     INTENT(IN) :: plane              !< slice location in terms of grid nodes
      REAL(f_real),     DIMENSION(:,:,:), INTENT(IN) :: field              !< random field
      CHARACTER(LEN=*),                   INTENT(IN) :: filename           !< output file name
      INTEGER(f_int)                                 :: rank, np, ierr, cut

      !-----------------------------------------------------------------------------------------------------------------------------

      ! get rank number
      CALL mpi_comm_rank(mpi_comm_world, rank, ierr)

      ! get number of tasks
      CALL mpi_comm_size(mpi_comm_world, np, ierr)

      ! these are global variables stored in "m_scarflib_common" and used during slicing
      ALLOCATE(gs(3, 0:np - 1), ge(3, 0:np - 1))

      ! store process-dependent grid nodes into global variable
      gs(:, rank) = obj%fs
      ge(:, rank) = obj%fe

      ! store total (process-independent) grid nodes into gloabl variable (defined in  "m_scarflib_common")
      npts = n

      ! share info amongst all processes
      CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, gs, 3, mpi_integer, mpi_comm_world, ierr)
      CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, ge, 3, mpi_integer, mpi_comm_world, ierr)

      ! IF (ierr .ne. 0) THEN
      !   stderr%code('error code', ierr)
      !   stderr%tree('m_scarflib -> scarf_io_slice -> mpi_allgather')
      ! ENDIF

      SELECT CASE(axis)
      CASE('x')
        cut = 1
      CASE('y')
        cut = 2
      CASE('z')
        cut = 3
      CASE DEFAULT
        cut = 0
      END SELECT

      !IF (plane .gt. npts(cut)) stderr%logic('plane', plane)
      !IF (cut .eq. 0)           stderr%arg('axis', axis)

      !IF (stderr%new) stderr%tree('m_scarflib -> scarf_io_slice')

      CALL write_slice(cut, plane, field, filename)

      DEALLOCATE(gs, ge)

    END SUBROUTINE scarf_io_slice

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE scarf_io_one(n, field, filename, nwriters)

      ! Purpose:
      !   To write to disk whole random field in a single file. The number of MPI processes with I/O permit can be defined in the
      !   arguments list. Works for structured meshes only.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(f_int),   DIMENSION(3),               INTENT(IN) :: n                     !< grid points in the global mesh
      REAL(f_real),     DIMENSION(:,:,:),           INTENT(IN) :: field                 !< random field
      CHARACTER(len=*),                             INTENT(IN) :: filename              !< output file name
      INTEGER(f_int),                     OPTIONAL, INTENT(IN) :: nwriters              !< number of concurrent I/O processes (empty = 1)
      INTEGER(f_int)                                           :: rank, np, nw, ierr

      !-----------------------------------------------------------------------------------------------------------------------------

      ! get rank number
      CALL mpi_comm_rank(mpi_comm_world, rank, ierr)

      ! get number of tasks
      CALL mpi_comm_size(mpi_comm_world, np, ierr)

      ! these are global variables stored in "scarflib_common"
      ALLOCATE(gs(3, 0:np - 1), ge(3, 0:np - 1))

      ! store global indices
      gs(:, rank) = obj%fs
      ge(:, rank) = obj%fe

      npts = n

      ! share info amongst all processes
      CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, gs, 3, mpi_integer, mpi_comm_world, ierr)
      CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, ge, 3, mpi_integer, mpi_comm_world, ierr)

      ! default case is to USE one writer only
      nw = 1

      ! we cannot have more writers than available processes
      IF (PRESENT(nwriters)) nw = MIN(np, nwriters)

      CALL write_one(field, filename, nw)

      DEALLOCATE(gs, ge)

    END SUBROUTINE scarf_io_one

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE m_scarflib
