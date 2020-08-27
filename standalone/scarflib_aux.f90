MODULE m_scarflib_aux

  ! Purpose:
  !   To provide routines necessary to set-up a sample driver program.
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   04/05/20                  original version
  !

  USE, INTRINSIC     :: iso_fortran_env, only: stdout => output_unit
  USE, NON_INTRINSIC :: m_scarflib_common

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: sample_mesh
  PUBLIC :: watch_start, watch_stop                      !< these subroutines are defined in "m_scarflib_common"

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE sample_mesh(rank, ntasks, n, fs, fe)

      ! Purpose:
      !   To create a sample structured mesh to be used in the driver program and to distribute it amonst available processes.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(f_int),               INTENT(IN)  :: rank, ntasks
      INTEGER(f_int), DIMENSION(:), INTENT(IN)  :: n
      INTEGER(f_int), DIMENSION(:), INTENT(OUT) :: fs, fe
      INTEGER(f_int)                            :: topo, ierr
      INTEGER(f_int), DIMENSION(3)              :: dims, coords, fs3, fe3

      !-----------------------------------------------------------------------------------------------------------------------------

      ! set global variables in "scarflib_common" and needed by "best_config"
      world_size = ntasks

      IF (SIZE(n) .eq. 3) THEN
        npts = n
      ELSE
        npts = [n(1), n(2), 1]
      ENDIF

      CALL best_config(dims)

      ! IF (rank .eq. 0) WRITE(stdout, *) 'CPU-grid in driver program: ', dims

      ! create topology
      CALL mpi_cart_create(mpi_comm_world, 3, dims, [.false., .false., .false.], .true., topo, ierr)

      CALL mpi_cart_coords(topo, rank, 3, coords, ierr)

      ! split each axis and assign points to current process
      CALL mpi_split_task(n, dims, coords, fs3, fe3)

      IF (SIZE(fs) .eq. 3) THEN
        fs = fs3
        fe = fe3
      ELSE
        fs = fs3(1:2)
        fe = fe3(1:2)
      ENDIF

    END SUBROUTINE sample_mesh

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE mpi_split_task(n, ntasks, rank, fs, fe)

      ! Purpose:
      !   To evenly distribute elements of vector amongst processes, returning first and last index for each process.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(f_int), DIMENSION(3), INTENT(IN)  :: n                          !< number of points to be split
      INTEGER(f_int), DIMENSION(3), INTENT(IN)  :: ntasks                     !< number of mpi processes
      INTEGER(f_int), DIMENSION(3), INTENT(IN)  :: rank                       !< rank of current preocess
      INTEGER(f_int), DIMENSION(3), INTENT(OUT) :: fs, fe                     !< 1st/last indices
      INTEGER(f_int)                            :: i                          !< counter

      !------------------------------------------------------------------------------------------------------------------------------

      DO i = 1, 3
        fs(i) = 1 + INT( REAL(n(i), f_real) / REAL(ntasks(i), f_real) * REAL(rank(i), f_real) )
        fe(i) = INT( REAL(n(i), f_real) / REAL(ntasks(i), f_real) * REAL(rank(i) + 1, f_real) )
      ENDDO

    END SUBROUTINE mpi_split_task

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
            IF (match() .eqv. .true.) CYCLE

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

      LOGICAL FUNCTION match()

        INTEGER(f_int) :: i

        !---------------------------------------------------------------------------------------------------------------------------

        match = .false.

        DO i = 1, c
          match = ANY(v1(1) .eq. list(:, i)) .and. ANY(v1(2) .eq. list(:, i)) .and. ANY(v1(3) .eq. list(:, i))
          IF (match .eqv. .true.) EXIT
        ENDDO

      END FUNCTION match

    END SUBROUTINE best_config

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

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

      measure = ABS(side(1) - side(2)) + ABS(side(1) - side(3)) + ABS(side(2) - side(3))

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


END MODULE m_scarflib_aux
