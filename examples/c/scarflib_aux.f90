MODULE SCARFLIB_AUX

  USE, INTRINSIC     :: ISO_C_BINDING
  !USE, NON_INTRINSIC :: MPI
  USE, NON_INTRINSIC :: M_SCARFLIB_COMMON

  IMPLICIT NONE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: SAMPLE_MESH

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE SAMPLE_MESH(RANK, NTASKS, N, FS, FE) BIND(C, NAME="sample_mesh")

      INTEGER(C_INT),               INTENT(IN)  :: RANK, NTASKS
      INTEGER(C_INT), DIMENSION(3), INTENT(IN)  :: N
      INTEGER(C_INT), DIMENSION(3), INTENT(OUT) :: FS, FE
      INTEGER(f_int)                              :: TOPO, IERR
      INTEGER(f_int),   DIMENSION(3)              :: DIMS, COORDS

      !-----------------------------------------------------------------------------------------------------------------------------

      ! SET GLOBAL VARIABLES IN "SCARFLIB_COMMON" AND NEEDED BY "BEST_CONFIG"
      WORLD_SIZE = NTASKS
      NPTS       = N

      CALL BEST_CONFIG(DIMS)

      IF (RANK .EQ. 0) PRINT*, 'CPU GRID IN DRIVER PROGRAM: ', DIMS

      ! CREATE TOPOLOGY
      CALL MPI_CART_CREATE(MPI_COMM_WORLD, 3, DIMS, [.FALSE., .FALSE., .FALSE.], .TRUE., TOPO, IERR)

      CALL MPI_CART_COORDS(TOPO, RANK, 3, COORDS, IERR)

      ! SPLIT EACH AXIS AND ASSIGN POINTS TO CURRENT PROCESS
      CALL MPI_SPLIT_TASK(N, DIMS, COORDS, FS, FE)

    END SUBROUTINE SAMPLE_MESH

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE MPI_SPLIT_TASK(N, NTASKS, RANK, FS, FE)

      INTEGER(f_int), DIMENSION(3), INTENT(IN)  :: N                          !< NUMBER OF POINTS TO BE SPLIT
      INTEGER(f_int), DIMENSION(3), INTENT(IN)  :: NTASKS                     !< NUMBER OF MPI PROCESSES
      INTEGER(f_int), DIMENSION(3), INTENT(IN)  :: RANK                       !< RANK OF CURRENT PREOCESS
      INTEGER(f_int), DIMENSION(3), INTENT(OUT) :: FS, FE                     !< 1ST/LAST INDICES
      INTEGER(f_int)                            :: I                          !< COUNTER

      !------------------------------------------------------------------------------------------------------------------------------

      DO I = 1, 3
        FS(I) = 1 + INT( REAL(N(I), f_real) / REAL(NTASKS(I), f_real) * REAL(RANK(I), f_real) )
        FE(I) = INT( REAL(N(I), f_real) / REAL(NTASKS(I), f_real) * REAL(RANK(I) + 1, f_real) )
      ENDDO

    END SUBROUTINE MPI_SPLIT_TASK

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE BEST_CONFIG(DIMS)

      ! COMPUTE BEST GRID OF PROCESSES

      INTEGER(f_int),              DIMENSION(3), INTENT(OUT) :: DIMS
      !INTEGER(f_int),                            INTENT(OUT) :: NITER

      INTEGER(f_int)                                         :: I, J, K, L, C
      INTEGER(f_int)                                         :: N2, N3 !, N
      INTEGER(f_int)                                         :: A, B
      INTEGER(f_int),              DIMENSION(3)              :: V1, V2
      INTEGER(f_int), ALLOCATABLE, DIMENSION(:,:)            :: FACT3, FACT2
      INTEGER(f_int), ALLOCATABLE, DIMENSION(:,:)            :: LIST

      REAL(f_real) :: N, NITER

      !-----------------------------------------------------------------------------------------------------------------------------

      !NITER = HUGE(1)
      NITER = HUGE(1._f_real)

      C = 0

      ! FACTORISE NUMBER OF AVAILABLE PROCESSES
      CALL FACTORIZATION(WORLD_SIZE, FACT3)

      N3 = SIZE(FACT3, 2)

      ALLOCATE(LIST(3, N3 * 10))

      LIST(:,:) = 0

      ! LOOP OVER FACTORISED PROCESSES
      DO L = 1, N3

        DO K = 1, 2

          IF (K .EQ. 1) THEN
            A = FACT3(1, L)
            B = FACT3(2, L)
          ELSE
            B = FACT3(1, L)
            A = FACT3(2, L)
          ENDIF

          CALL FACTORIZATION(B, FACT2)

          N2 = SIZE(FACT2, 2)

          DO J = 1, N2

            V1 = [A, FACT2(:, J)]

            ! SKIP TO NEXT PROCESSES GRID IF ALREADY ANALYSED
            IF (MATCH(V1, C, LIST) .EQV. .TRUE.) CYCLE

            C = C + 1

            ! ADD TO LIST
            LIST(:, C) = V1

            DO I = 0, 2

              V1 = CSHIFT(V1, 1)

              CALL TEST_CONFIG(V1, N)

              IF (N .LT. NITER) THEN
                DIMS  = V1
                NITER = N
              ENDIF

              V2 = [V1(1), V1(3), V1(2)]

              CALL TEST_CONFIG(V2, N)

              IF (N .LT. NITER) THEN
                DIMS  = V2
                NITER = N
              ENDIF

            ENDDO
            ! END PERMUTATIONS

          ENDDO
          ! END LOOP OVER FACTOR PAIRS FOR "A/B"

        ENDDO
        ! END LOOP OVER "A/B"

      ENDDO
      ! END LOOP OVER FACTOR PAIRS FOR "WORLD_SIZE"

      DEALLOCATE(FACT3, FACT2, LIST)

      !-----------------------------------------------------------------------------------------------------------------------------
      !-----------------------------------------------------------------------------------------------------------------------------

      CONTAINS

      LOGICAL FUNCTION MATCH(VEC, IMAX, LIST)

        INTEGER(f_int), DIMENSION(3),   INTENT(IN) :: VEC
        INTEGER(f_int),                 INTENT(IN) :: IMAX
        INTEGER(f_int), DIMENSION(:,:), INTENT(IN) :: LIST
        INTEGER(f_int)                             :: I

        !---------------------------------------------------------------------------------------------------------------------------

        MATCH = .FALSE.

        DO I = 1, IMAX
          MATCH = ANY(V1(1) .EQ. LIST(:, I)) .AND. ANY(V1(2) .EQ. LIST(:, I)) .AND. ANY(V1(3) .EQ. LIST(:, I))
          IF (MATCH .EQV. .TRUE.) EXIT
        ENDDO

      END FUNCTION MATCH

    END SUBROUTINE BEST_CONFIG

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE TEST_CONFIG(DIMS, MEASURE)

      INTEGER(f_int), DIMENSION(3), INTENT(IN)  :: DIMS
      REAL(f_real),                  INTENT(OUT) :: MEASURE
      REAL(f_real),    DIMENSION(3)              :: SIDE

      !-----------------------------------------------------------------------------------------------------------------------------

      SIDE = REAL(NPTS, f_real) / REAL(DIMS, f_real)

      MEASURE = ABS(SIDE(1) - SIDE(2)) + ABS(SIDE(1) - SIDE(3)) + ABS(SIDE(2) - SIDE(3))

    END SUBROUTINE TEST_CONFIG

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE FACTORIZATION(N, FACTORS)

      ! FACTORISE INTEGER "N" BASED ON TRIAL DIVISION METHOD

      INTEGER(f_int),                              INTENT(IN)  :: N
      INTEGER(f_int), ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: FACTORS
      INTEGER(f_int)                                           :: I, C, S
      INTEGER(f_int)                                           :: X
      INTEGER(f_int), ALLOCATABLE, DIMENSION(:,:)              :: BUFFER

      !-----------------------------------------------------------------------------------------------------------------------------

      ! MAX POSSIBLE FACTOR
      S = FLOOR(SQRT(REAL(N, f_real)))

      ALLOCATE(BUFFER(2, S))

      BUFFER(:,:) = 0

      ! TEST FACTORS
      DO I = 1, S

        X = N / I

        IF (MOD(N, I) .EQ. 0) THEN
          BUFFER(1, I) = I                          !< ADD FACTOR ...
          BUFFER(2, I) = X                          !< ... AND ITS COMPANION
        ENDIF

      ENDDO

      ! ACTUAL FACTORS FOUND
      I = COUNT(BUFFER(1, :) .NE. 0)

      ALLOCATE(FACTORS(2, I))

      ! COPY FACTORS TO OUTPUT ARRAY
      C = 0
      DO I = 1, S
        IF (BUFFER(1, I) .NE. 0) THEN
          C = C + 1
          FACTORS(:, C) = BUFFER(:, I)
        ENDIF
      ENDDO

      DEALLOCATE(BUFFER)

    END SUBROUTINE FACTORIZATION

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


END MODULE SCARFLIB_AUX
