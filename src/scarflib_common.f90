MODULE m_scarflib_common

  ! Purpose:
  !   To define constants and subprograms common within the scarf3d library.
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   04/05/20                  original version
  !   11/05/20                  updated macro for double-precision
  !

  USE, INTRINSIC     :: iso_fortran_env
  USE, INTRINSIC     :: iso_c_binding
  USE, NON_INTRINSIC :: mpi

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! make modules, constants, shared variables and subprograms accessible to other modules via host-association
  PUBLIC

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTERFACE variance
    MODULE PROCEDURE variance_1d, variance_3d
  END INTERFACE

  INTERFACE mean
    MODULE PROCEDURE mean_1d, mean_3d
  END INTERFACE

  INTERFACE tapering
    MODULE PROCEDURE tapering_structured, tapering_unstructured
  END INTERFACE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! define type for integers
  INTEGER,        PARAMETER :: f_int = int32

  ! define precision of floating point numbers
#ifdef DOUBLE_PREC
  INTEGER(f_int), PARAMETER :: f_real    = real64
  INTEGER(f_int), PARAMETER :: f_dble    = real64
  INTEGER(f_int), PARAMETER :: c_real    = c_double
  INTEGER(f_int), PARAMETER :: c_cplx    = c_double_complex
  INTEGER(f_int), PARAMETER :: real_type = mpi_double_precision
  INTEGER(f_int), PARAMETER :: cplx_type = mpi_double_complex
#else
  INTEGER(f_int), PARAMETER :: f_real    = real32
  INTEGER(f_int), PARAMETER :: f_dble    = real64
  INTEGER(f_int), PARAMETER :: c_real    = c_float
  INTEGER(f_int), PARAMETER :: c_cplx    = c_float_complex
  INTEGER(f_int), PARAMETER :: real_type = mpi_real
  INTEGER(f_int), PARAMETER :: cplx_type = mpi_complex
#endif

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  ! total number of points along each axis (whole model)
  INTEGER(f_int),              DIMENSION(3)             :: npts

  ! mpi stuff
  INTEGER(f_int)                                        :: world_rank, world_size

  ! global model indices for all mpi processes
  INTEGER(f_int), ALLOCATABLE, DIMENSION(:,:)           :: gs, ge

  ! variables
  REAL(f_real),                               PARAMETER :: pi = 3.141592653589793_f_real

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE watch_start(tictoc, comm)

      ! Purpose:
      !   To start the MPI stopwatch. Timing is in double-precision. If specific communicator handle not given, mpi_comm_world is
      !   used.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_dble),             INTENT(OUT) :: tictoc                            !< initial time
      INTEGER(f_int), OPTIONAL, INTENT(IN)  :: comm                              !< communicator handle
      INTEGER(f_int)                        :: ierr

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (.NOT.PRESENT(comm)) THEN
        CALL mpi_barrier(mpi_comm_world, ierr)
      ELSE
        CALL mpi_barrier(comm, ierr)
      ENDIF

      tictoc = mpi_wtime()

    END SUBROUTINE watch_start

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE watch_stop(tictoc, comm)

      ! Purpose:
      !   To stop the MPI stopwatch and return elapsed time. Timing is in double-precision.  If specific communicator handle not given,
      !   mpi_comm_world is used.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_dble),             INTENT(INOUT) :: tictoc                          !< elapsed time
      INTEGER(f_int), OPTIONAL, INTENT(IN)    :: comm                            !< communicator handle
      INTEGER(f_int)                          :: ierr

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (.NOT.PRESENT(comm)) THEN
        CALL mpi_barrier(mpi_comm_world, ierr)
      ELSE
        CALL mpi_barrier(comm, ierr)
      ENDIF

      tictoc = mpi_wtime() - tictoc

    END SUBROUTINE watch_stop

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    ! SUBROUTINE set_stream(seed)
    !
    !   ! initialse the random number generator
    !
    !   INTEGER(f_int),                           INTENT(IN) :: seed                      !< USEr-defined seed number
    !   INTEGER(f_int)                                       :: seed_size                 !< store array SIZE
    !   INTEGER(f_int), ALLOCATABLE, DIMENSION(:)            :: tmp_seed                  !< USEd to initialise the random generator
    !
    !   !-----------------------------------------------------------------------------------------------------------------------------
    !
    !   CALL random_seed(SIZE = seed_size)
    !
    !   ALLOCATE(tmp_seed(seed_size))
    !
    !   tmp_seed = seed
    !
    !   CALL random_seed(put = tmp_seed)
    !
    !   DEALLOCATE(tmp_seed)
    !
    ! END SUBROUTINE set_stream

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE write_single(v, filename, nwriters)

      ! Purpose:
      !   To allow each MPI process writes its own file with data. Only 'nwriters' files at most are written at the same time. This
      !   subroutine works for structured meshes only
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_real),               DIMENSION(:,:,:),          INTENT(IN) :: v                !< random field
      CHARACTER(LEN=*),                                         INTENT(IN) :: filename         !< base-name of output file
      INTEGER(f_int),                                        INTENT(IN) :: nwriters         !< number of concurrent files written
      CHARACTER(:),      ALLOCATABLE                                       :: myname
      INTEGER(f_int)                                                    :: i, n
      INTEGER(f_int)                                                    :: fh, ierr
      INTEGER(f_int),             DIMENSION(mpi_status_size)            :: status

      !-----------------------------------------------------------------------------------------------------------------------------

      DO i = 0, world_size - 1, nwriters

        n = i + nwriters

        ! only processes in the range [i, i + nwriters), i.e. max "nwriters", enter this block every time
        IF ( (world_rank .ge. i) .and. (world_rank .lt. n) ) THEN

          myname = filename // '_x=' // num2str(gs(1, world_rank)) // '_' // num2str(ge(1, world_rank)) //    &
                               '_y=' // num2str(gs(2, world_rank)) // '_' // num2str(ge(2, world_rank)) //    &
                               '_z=' // num2str(gs(3, world_rank)) // '_' // num2str(ge(3, world_rank))

          CALL mpi_file_open(mpi_comm_self, myname, mpi_mode_create + mpi_mode_wronly, mpi_info_null, fh, ierr)

          CALL mpi_file_write(fh, v, SIZE(v), real_type, status, ierr)

          CALL mpi_file_close(fh, ierr)

        ENDIF

        CALL mpi_barrier(mpi_comm_world, ierr)

      ENDDO

    END SUBROUTINE write_single

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE write_one(v, filename, nwriters)

      ! Purpose:
      !   To write data to a single file. Only 'nwriters' MPI process are allowed to write to disk after gathering data from all other
      !   processess. This subroutine works for structured meshes only.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_real),                 DIMENSION(:,:,:),           INTENT(IN) :: v                                  !< random field
      CHARACTER(LEN=*),                                         INTENT(IN) :: filename                           !< name of output file
      INTEGER(f_int),                                           INTENT(IN) :: nwriters                           !< number of I/O processes
      CHARACTER(:),     ALLOCATABLE                                        :: dchar
      INTEGER(f_int)                                                       :: i, j, n
      INTEGER(f_int)                                                       :: nmax, maxtasks
      INTEGER(f_int)                                                       :: color, ntasks, newcomm, wrtcomm
      INTEGER(f_int)                                                       :: rank, rcvrank, request, info
      INTEGER(f_int)                                                       :: newtype, ierr, fh
      INTEGER(f_int),                DIMENSION(mpi_status_size)            :: status
      INTEGER(f_int),                DIMENSION(3)                          :: subsizes, starts
      INTEGER(f_int),   ALLOCATABLE, DIMENSION(:)                          :: rank0, rank1
      INTEGER(f_int),   ALLOCATABLE, DIMENSION(:)                          :: new2world
      REAL(f_real),     ALLOCATABLE, DIMENSION(:)                          :: buffer

      !-----------------------------------------------------------------------------------------------------------------------------

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! organize processes in 'nwriters' communicators and prepare a map of the local-to-global ranks
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ALLOCATE(rank0(nwriters), rank1(nwriters))

      ! distribute 'world_size' elements between 'nwriters'
      CALL split_task(world_size, nwriters, rank0, rank1)

      ! lowest rank must be zero
      rank0 = rank0 - 1
      rank1 = rank1 - 1

      DO i = 1, nwriters
        IF ( (world_rank .ge. rank0(i)) .and. (world_rank .le. rank1(i)) ) color = i
      ENDDO

      CALL mpi_comm_split(mpi_comm_world, color, world_rank, newcomm, ierr)

      CALL mpi_comm_rank(newcomm, rank, ierr)
      CALL mpi_comm_size(newcomm, ntasks, ierr)

      ALLOCATE(new2world(0:ntasks - 1))

      ! i-th rank in "newcomm" has its global rank contained in "new2world"
      CALL map_ranks(newcomm, mpi_comm_world, new2world)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! processes with "rank=0" in every "newcomm" are writers: gather them in separate communicator "wrtcomm"
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      color = 0

      IF (rank .eq. 0) color = 1

      CALL mpi_comm_split(mpi_comm_world, color, world_rank, wrtcomm, ierr)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! determine maximum number of tasks amongst all "newcomm" communicators, needed to for collective WRITE CALL below
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      maxtasks = ntasks

      CALL mpi_allreduce(mpi_in_place, maxtasks, 1, mpi_integer, mpi_max, mpi_comm_world, ierr)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! determine maximum SIZE of input array "v" and ALLOCATE "buffer" accordingly. this way we don't need to reallocate "buffer"
      ! every time we get data from a process
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      n = SIZE(v)

      nmax = n

      CALL mpi_allreduce(mpi_in_place, nmax, 1, mpi_integer, mpi_max, mpi_comm_world, ierr)

      ALLOCATE(buffer(nmax))

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! each process makes non-blocking send. i/o processes collect data and WRITE them to file at the right location. a writer may
      ! WRITE twice the same data IF "world_size" not multiple of "nwriters": this is necessary becaUSE we USE collective WRITE CALL
      ! that must be called by all writers at the same time.
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! each process of communicator "newcomm" sends its data to the writer of the communicator
      CALL mpi_isend(v, n, real_type, 0, 1, newcomm, request, ierr)

      ! only writers can work on files
      IF (rank .eq. 0) THEN

        info = mpi_info_null

        ! set "striping_factor" equal to number of i/o processes. "striping_unit" is arbitray set to 128mb (reasonable value).
#ifdef PFS
        CALL mpi_info_create(info, ierr)

        dchar = num2str(nwriters)

        CALL mpi_info_set(info, "striping_factor", dchar, ierr)

        DEALLOCATE(dchar)

        dchar = num2str(128*1024*1024)

        CALL mpi_info_set(info, "striping_unit", dchar, ierr)

        DEALLOCATE(dchar)
#endif

        CALL mpi_file_open(wrtcomm, filename, mpi_mode_create + mpi_mode_wronly, info, fh, ierr)

        ! loop until all data are received
        DO j = 0, maxtasks - 1

          ! stop receiving if we have already collected all data for "newcomm"
          IF (j .le. ntasks - 1) THEN

            ! receive data
            CALL mpi_recv(buffer, SIZE(buffer), real_type, mpi_any_source, mpi_any_tag, newcomm, status, ierr)

            ! find who sent data
            rcvrank = status(mpi_source)

            ! map its local rank into global one
            rcvrank = new2world(rcvrank)

          ENDIF

          subsizes(:) = ge(:, rcvrank) - gs(:, rcvrank) + 1
          !starts(:)   = gs(:, rcvrank) - 1
          starts(:)   = gs(:, rcvrank) - MINVAL(gs, dim = 2)

          CALL mpi_type_create_subarray(3, npts, subsizes, starts, mpi_order_fortran, real_type, newtype, ierr)

          CALL mpi_type_commit(newtype, ierr)

          ! view is controlled by derived datatype
          CALL mpi_file_set_view(fh, 0_mpi_offset_kind, real_type, newtype, 'native', mpi_info_null, ierr)

          CALL mpi_file_write_all(fh, buffer, PRODUCT(subsizes), real_type, mpi_status_ignore, ierr)

          CALL mpi_type_free(newtype, ierr)

        ENDDO

        CALL mpi_file_close(fh, ierr)

#ifdef PFS
        CALL mpi_info_free(info, ierr)
#endif

      ENDIF

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! free resources
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      DEALLOCATE(rank0, rank1, new2world, buffer)

      CALL mpi_barrier(mpi_comm_world, ierr)

      CALL mpi_comm_free(newcomm, ierr)
      CALL mpi_comm_free(wrtcomm, ierr)

    END SUBROUTINE write_one

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION num2str(v)

      ! Purpose:
      !   To convert an integer into a string.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(f_int),                     INTENT(IN) :: v               !< integer to be converted
      CHARACTER(:),           ALLOCATABLE            :: num2str
      CHARACTER(RANGE(v) + 2)                        :: dum

      !-----------------------------------------------------------------------------------------------------------------------------

      WRITE(dum, '(i0)') v

      num2str = TRIM(dum)

    END FUNCTION num2str

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE split_task(npts, ntasks, i0, i1)

      ! Purpose:
      !   To evenly distribute elements of vector amongst processes, returning first and last index for each process.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(f_int),                        INTENT(IN)  :: npts                       !< number of elements to be distributed
      INTEGER(f_int),                        INTENT(IN)  :: ntasks                     !< number of MPI processes
      INTEGER(f_int), DIMENSION(0:ntasks-1), INTENT(OUT) :: i0, i1                     !< 1st/last index
      INTEGER(f_int)                                     :: p

      !------------------------------------------------------------------------------------------------------------------------------

      DO p = 0, ntasks - 1
        i0(p) = 1 + INT( REAL(npts, f_real) / REAL(ntasks, f_real) * REAL(p, f_real) )
        i1(p) = INT( REAL(npts, f_real) / REAL(ntasks, f_real) * REAL(p + 1, f_real) )
      ENDDO

    END SUBROUTINE split_task

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE map_ranks(src_comm, dest_comm, dest_ranks)

      ! Purpose:
      !   To map ranks from one communicator to another.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(f_int),                             INTENT(IN)    :: src_comm                     !< source communicator handle
      INTEGER(f_int),                             INTENT(IN)    :: dest_comm                    !< destination communicator handle
      INTEGER(f_int), DIMENSION(:),               INTENT(INOUT) :: dest_ranks                   !< ranks after mapping
      INTEGER(f_int)                                            :: i
      INTEGER(f_int)                                            :: ierr, src_group
      INTEGER(f_int)                                            :: dest_group
      INTEGER(f_int)                                            :: ntasks
      INTEGER(f_int), DIMENSION(SIZE(dest_ranks))               :: src_ranks

      !-----------------------------------------------------------------------------------------------------------------------------

      ntasks = SIZE(src_ranks)

      ! ranks in source communicator
      DO i = 1, ntasks
        src_ranks(i) = i - 1
      ENDDO

      ! group associated to source communicator
      CALL mpi_comm_group(src_comm, src_group, ierr)

      ! group associated to destination communicator
      CALL mpi_comm_group(dest_comm, dest_group, ierr)

      ! map ranks in "group" into ranks in "mpi_group_world"
      CALL mpi_group_translate_ranks(src_group, ntasks, src_ranks, dest_group, dest_ranks, ierr)

      CALL mpi_group_free(src_group, ierr)
      CALL mpi_group_free(dest_group, ierr)

    END SUBROUTINE map_ranks

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE write_slice(axis, plane, v, filename)

      ! Purpose:
      !   To write to disk a single file containing a slice of the random field. Each MPI process writes its own part based on the
      !   subarray datatype. A slice can be oriented along the x-, y- or z-axis. This subroutine works for structured meshes only.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !   11/05/20                  handle case when min(gs) .ne. 1
      !

      INTEGER(f_int),                                  INTENT(IN) :: axis                               !< axis (1=x,2=y,3=z) cut by the slice
      INTEGER(f_int),                                  INTENT(IN) :: plane                              !< slice index in global coordinates
      REAL(f_real),                  DIMENSION(:,:,:), INTENT(IN) :: v                                  !< random field
      CHARACTER(LEN=*),                                INTENT(IN) :: filename                           !< name  of output file
      INTEGER(f_int)                                              :: i, j, c
      INTEGER(f_int)                                              :: newtype, ierr, fh, comm, color
      INTEGER(f_int),                DIMENSION(2)                 :: dir
      INTEGER(f_int),                DIMENSION(2)                 :: subsizes, starts, dims
      LOGICAL                                                     :: bool
      REAL(f_real),     ALLOCATABLE, DIMENSION(:,:)               :: buffer

      !-----------------------------------------------------------------------------------------------------------------------------

      SELECT CASE(axis)
      CASE(1)
        dir  = [2, 3]
      CASE (2)
        dir  = [1, 3]
      CASE (3)
        dir  = [1, 2]
      END SELECT

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! create a new communicator only with processes having data on the slice
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      color = 0

      bool = (plane .ge. gs(axis, world_rank)) .and. (plane .le. ge(axis, world_rank))

      DO i = 1, 2
        subsizes(i) = SIZE(v, dir(i))
        starts(i)   = gs(dir(i), world_rank) - MINVAL(gs(dir(i), :))
        dims(i)     = npts(dir(i))
      ENDDO

      ! global-to-local index mapping
      c = plane - gs(axis, world_rank) + 1

      ! update color of calling process only if it contains part of the slice
      IF (bool) color = 1

      ! split communicator
      CALL mpi_comm_split(mpi_comm_world, color, world_rank, comm, ierr)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! copy data into a contiguous buffer and write to disk
      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

      ! only processes in the new communicator enter section below
      IF (color .eq. 1) THEN

        ALLOCATE(buffer(subsizes(1), subsizes(2)))

        IF (axis .eq. 1) THEN

          DO j = 1, subsizes(2)
            DO i = 1, subsizes(1)
              buffer(i, j) = v(c, i, j)
            ENDDO
          ENDDO

        ELSEIF (axis .eq. 2) THEN

          DO j = 1, subsizes(2)
            DO i = 1, subsizes(1)
              buffer(i, j) = v(i, c, j)
            ENDDO
          ENDDO

        ELSEIF (axis .eq. 3) THEN

          DO j = 1, subsizes(2)
            DO i = 1, subsizes(1)
              buffer(i, j) = v(i, j, c)
            ENDDO
          ENDDO

        ENDIF

        ! create subarray
        CALL mpi_type_create_subarray(2, dims, subsizes, starts, mpi_order_fortran, real_type, newtype, ierr)

        CALL mpi_type_commit(newtype, ierr)

        CALL mpi_file_open(comm, filename, mpi_mode_create + mpi_mode_wronly, mpi_info_null, fh, ierr)

        ! set file view based on subarray datatype
        CALL mpi_file_set_view(fh, 0_mpi_offset_kind, real_type, newtype, 'native', mpi_info_null, ierr)

        CALL mpi_file_write_all(fh, buffer, PRODUCT(subsizes), real_type, mpi_status_ignore, ierr)

        CALL mpi_file_close(fh, ierr)

        ! release derived datatype
        CALL mpi_type_free(newtype, ierr)

        ! realease local communicator
        CALL mpi_comm_free(comm, ierr)

        DEALLOCATE(buffer)

      ENDIF

      CALL mpi_barrier(mpi_comm_world, ierr)

    END SUBROUTINE write_slice

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(f_real) FUNCTION variance_1d(r)

      ! Purpose:
      !   To compute variance according to the compensated-summation version of the two-pass algorithm. Calculations in double-precision.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_real),  DIMENSION(:), INTENT(IN) :: r                      !< vector of data
      INTEGER(f_int)                          :: i, n
      REAL(f_dble)                            :: mu, s1, s2, x, v

      !-----------------------------------------------------------------------------------------------------------------------------

      n = SIZE(r)

      mu = mean(r)

      s1 = 0._f_dble
      s2 = 0._f_dble

      DO i = 1, n
        x = REAL(r(i), f_dble) - mu
        s1 = s1 + x
        s2 = s2 + x**2
      ENDDO

      s1 = (s1**2) / REAL(n, f_dble)

      v = (s2 - s1) / REAL(n - 1, f_dble)

      variance_1d = REAL(v, f_real)

    END FUNCTION variance_1d

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(f_real) FUNCTION variance_3d(r)

      ! Purpose:
      !   To compute variance according to the compensated-summation version of the two-pass algorithm. Calculations in double-precision.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_real),  DIMENSION(:,:,:), INTENT(IN) :: r                    !< array of data
      INTEGER(f_int)                              :: i, j, k
      INTEGER(f_int)                              :: nx, ny, nz
      REAL(f_dble)                                :: mu, s1, s2, x, v

      !-----------------------------------------------------------------------------------------------------------------------------

      nx = SIZE(r, 1)
      ny = SIZE(r, 2)
      nz = SIZE(r, 3)

      mu = mean(r)

      s1 = 0._f_dble
      s2 = 0._f_dble

      DO k = 1, nz
        DO j = 1, ny
          DO i = 1, nx
            x = REAL(r(i, j, k), f_dble) - mu
            s1 = s1 + x
            s2 = s2 + x**2
          ENDDO
        ENDDO
      ENDDO

      s1 = (s1**2) / REAL(SIZE(r), f_dble)

      v = (s2 - s1) / REAL(SIZE(r) - 1, f_dble)

      variance_3d = REAL(v, f_real)

    END FUNCTION variance_3d

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(f_real) FUNCTION mean_1d(r)

      ! Purpose:
      !   To compute average. Calculations in double-precision to reduce risk of cancellation.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_real),  DIMENSION(:), INTENT(IN) :: r               !< vector of data
      INTEGER(f_int)                          :: i
      INTEGER(f_int)                          :: n
      REAL(f_dble)                            :: v, c

      !-----------------------------------------------------------------------------------------------------------------------------

      n = SIZE(r)

      v = 0._f_dble
      c = 1._f_dble

      DO i = 1, n
        v = v + (REAL(r(i), f_dble) - v) / c
        c = c + 1._f_dble
      ENDDO

      ! return mean at desired precision
      mean_1d = REAL(v, f_real)

    END FUNCTION mean_1d

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(f_real) FUNCTION mean_3d(r)

      ! Purpose:
      !   To compute average. Calculations in double-precision to reduce risk of cancellation.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_real),  DIMENSION(:,:,:), INTENT(IN) :: r                       !< array of data
      INTEGER(f_int)                              :: i, j, k
      INTEGER(f_int)                              :: nx, ny, nz
      REAL(f_dble)                                :: v, c

      !-----------------------------------------------------------------------------------------------------------------------------

      nx = SIZE(r, 1)
      ny = SIZE(r, 2)
      nz = SIZE(r, 3)

      v = 0._f_dble
      c = 1._f_dble

      DO k = 1, nz
        DO j = 1, ny
          DO i = 1, nx
            v = v + (REAL(r(i, j, k), f_dble) - v) / c
            c = c + 1._f_dble
          ENDDO
        ENDDO
      ENDDO

      ! return mean with desired precision
      mean_3d = REAL(v, f_real)

    END FUNCTION mean_3d

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE parallel_variance(varset, avgset, nset, var, avg)

      ! Purpose:
      !   To compute variance and mean of a sequence of variance and mean values, each computed from a set of "nset" elements, based
      !   on the algorithm of Chan et al. (1979). Internal calculations in double-precision.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_real),   DIMENSION(:), INTENT(IN)  :: varset                             !< set of variance values
      REAL(f_real),   DIMENSION(:), INTENT(IN)  :: avgset                             !< set of mean values
      INTEGER(f_int), DIMENSION(:), INTENT(IN)  :: nset                               !< number of points for each set
      REAL(f_real),                 INTENT(OUT) :: var                                !< resulting variance
      REAL(f_real),                 INTENT(OUT) :: avg                                !< resulting average
      INTEGER(f_int)                            :: i
      REAL(f_dble)                              :: sig, mu, n, delta, m1, m2, m

      !-----------------------------------------------------------------------------------------------------------------------------

      ! variance, mean and number of points of first set
      sig = REAL(varset(1), f_dble)
      mu  = REAL(avgset(1), f_dble)
      n   = REAL(nset(1), f_dble)

      ! loop over set of variance/mean values
      DO i = 2, SIZE(nset)

        delta = REAL(avgset(i), f_dble) - mu

        m1 = sig * (n - 1._f_dble)
        m2 = REAL(varset(i), f_dble) * (nset(i) - 1._f_dble)

        m = m1 + m2 + delta**2 * n * REAL(nset(i), f_dble) / (n + REAL(nset(i), f_dble))

        ! resulting mean
        mu = (mu * n + REAL(avgset(i), f_dble) * REAL(nset(i), f_dble)) / (n + REAL(nset(i), f_dble))

        ! resulting number of points
        n = n + REAL(nset(i), f_dble)

        ! resulting variance
        sig = m / (n - 1._f_dble)

      ENDDO

      ! return with desired precision
      var = REAL(sig, f_real)
      avg = REAL(mu, f_real)

    END SUBROUTINE parallel_variance

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE tapering_structured (dh, fs, fe, field, poi, mute, taper)

      ! Purpose:
      !   To mute and/or taper the random field in a volume around a specific point. The random field at those grid points whose
      !   distance from the point is <= 'mute' are set to zero, while a cosine taper (ranging from 0 to 1) is applied if the distance
      !   is in the range ('mute' 'mute + taper']
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_real),                                                     INTENT(IN)    :: dh             !< grid-step
      INTEGER(f_int), DIMENSION(3),                                     INTENT(IN)    :: fs, fe         !< first/last index along each direction
      REAL(f_real),   DIMENSION(fs(1):fe(1), fs(2):fe(2), fs(3):fe(3)), INTENT(INOUT) :: field          !< random field
      REAL(f_real),   DIMENSION(3),                                     INTENT(IN)    :: poi            !< point-of-interest
      REAL(f_real),                                                     INTENT(IN)    :: mute, taper    !< radius for muting/tapering
      INTEGER(f_int)                                                                  :: i, j, k
      REAL(f_real)                                                                    :: d, ds
      REAL(f_real)                                                                    :: dx, dy, dz
      REAL(f_real)                                                                    :: t

      !-----------------------------------------------------------------------------------------------------------------------------

      ! total radius of tapered and/or muted volume
      ds = mute + taper

      DO k = fs(3), fe(3)

        dz = ((k - 1) * dh - poi(3))**2                                                !< distance along z from "poi"

        DO j = fs(2), fe(2)

          dy = ((j - 1) * dh - poi(2))**2                                              !< distance along y from "poi"

          DO i = fs(1), fe(1)

            dx = ((i - 1) * dh - poi(1))**2                                            !< distance along x from "poi"

            d = SQRT(dx + dy + dz)                                                     !< total distance

            ! mute if distance node-poi is below "mute"
            IF (d .le. mute) THEN

              field(i, j, k) = 0._f_real

            ! taper if distance node-poi is between "mute" and "taper + mute"
            ELSEIF ( (d .gt. mute) .and. (d .le. ds) ) THEN

              t = (d - mute) / taper                                                   !< taper parameter

              field(i, j, k) = field(i, j, k) * (0.5_f_real - 0.5_f_real * COS(t * pi))

            ENDIF

          ENDDO
        ENDDO
      ENDDO

    END SUBROUTINE tapering_structured

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE tapering_unstructured(x, y, z, field, poi, mute, taper)

      ! Purpose:
      !   To mute and/or taper the random field in a volume around a specific point. The random field at those grid points whose
      !   distance from the point is <= 'mute' are set to zero, while a cosine taper (ranging from 0 to 1) is applied if the distance
      !   is in the range ('mute' 'mute + taper']
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_real),  DIMENSION(:), INTENT(IN)    :: x, y, z
      REAL(f_real),  DIMENSION(:), INTENT(INOUT) :: field                 !< random field
      REAL(f_real),  DIMENSION(3), INTENT(IN)    :: poi                   !< point-of-interest (nodes)
      REAL(f_real),                INTENT(IN)    :: mute, taper           !< radius for muting/tapering
      INTEGER(f_int)                             :: i
      REAL(f_real)                               :: d, ds
      REAL(f_real)                               :: t

      !-----------------------------------------------------------------------------------------------------------------------------

      ! total radius of tapered and/or muted volume
      ds = mute + taper

      DO i = 1, SIZE(field)

        ! distance between grid point and "poi"
        d = SQRT( (x(i) - poi(1))**2 + (y(i) - poi(2))**2 + (z(i) - poi(3))**2 )

        ! mute IF distance node-poi is below "mute"
        IF (d .le. mute) THEN

          field(i) = 0._f_real

        ! taper IF distance node-poi is between "mute" and "mute + taper"
        ELSEIF ( (d .gt. mute) .and. (d .le. ds) ) THEN

          t = (d - mute) / taper                                             !< taper parameter

          field(i) = field(i) * (0.5_f_real - 0.5_f_real * COS(t * pi))

        ENDIF

      ENDDO

    END SUBROUTINE tapering_unstructured

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE m_scarflib_common
