MODULE SCARFLIB_AUX

  USE, NON_INTRINSIC :: MPI
  USE, NON_INTRINSIC :: SCARFLIB_COMMON
  USE, NON_INTRINSIC :: SCARFLIB_FFT

  IMPLICIT NONE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: SAMPLE_MESH

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE SAMPLE_MESH(RANK, NTASKS, N, FS, FE, X, Y, Z, V)

      INTEGER(IPP),                            INTENT(IN)  :: RANK, NTASKS
      INTEGER(IPP),              DIMENSION(3), INTENT(IN)  :: N
      INTEGER(IPP),              DIMENSION(3), INTENT(OUT) :: FS, FE
      REAL(FPP),    ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: X, Y, Z, V
      INTEGER(IPP)                                         :: TOPO, IERR
      INTEGER(IPP), DIMENSION(3)                           :: DIMS, COORDS

      !-----------------------------------------------------------------------------------------------------------------------------

      ! SET GLOBAL VARIABLES IN "SCARFLIB_COMMON" AND NEEDED BY "BEST_CONFIG"
      WORLD_SIZE = NTASKS
      NPTS       = N

      CALL BEST_CONFIG(DIMS)

      PRINT*, 'CPU GRID IN DRIVER: ', DIMS

      ! CREATE TOPOLOGY
      CALL MPI_CART_CREATE(MPI_COMM_WORLD, 3, DIMS, [.FALSE., .FALSE., .FALSE.], .TRUE., TOPO, IERR)

      CALL MPI_CART_COORDS(TOPO, RANK, 3, COORDS, IERR)

      ! SPLIT EACH AXIS AND ASSIGN POINTS TO CURRENT PROCESS
      CALL MPI_SPLIT_TASK(N, DIMS, COORDS, FS, FE)

      ! PREPARE ARRAY FOR RANDOM FIELD
      ALLOCATE(V((FE(1) - FS(1) + 1) * (FE(2) - FS(2) + 1) * (FE(3) - FS(3) + 1)))

      ! PREPARE ARRAY WITH MODEL POINTS POSITION
      ALLOCATE(X((FE(1) - FS(1) + 1) * (FE(2) - FS(2) + 1) * (FE(3) - FS(3) + 1)))
      ALLOCATE(Y((FE(1) - FS(1) + 1) * (FE(2) - FS(2) + 1) * (FE(3) - FS(3) + 1)))
      ALLOCATE(Z((FE(1) - FS(1) + 1) * (FE(2) - FS(2) + 1) * (FE(3) - FS(3) + 1)))

    END SUBROUTINE SAMPLE_MESH

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE MPI_SPLIT_TASK(N, NTASKS, RANK, FS, FE)

      INTEGER(IPP), DIMENSION(3), INTENT(IN)  :: N                          !< NUMBER OF POINTS TO BE SPLIT
      INTEGER(IPP), DIMENSION(3), INTENT(IN)  :: NTASKS                     !< NUMBER OF MPI PROCESSES
      INTEGER(IPP), DIMENSION(3), INTENT(IN)  :: RANK                       !< RANK OF CURRENT PREOCESS
      INTEGER(IPP), DIMENSION(3), INTENT(OUT) :: FS, FE                     !< 1ST/LAST INDICES
      INTEGER(IPP)                            :: I                          !< COUNTER

      !------------------------------------------------------------------------------------------------------------------------------

      DO I = 1, 3
        FS(I) = 1 + INT( REAL(N(I), FPP) / REAL(NTASKS(I), FPP) * REAL(RANK(I), FPP) )
        FE(I) = INT( REAL(N(I), FPP) / REAL(NTASKS(I), FPP) * REAL(RANK(I) + 1, FPP) )
      ENDDO

    END SUBROUTINE MPI_SPLIT_TASK

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE SCARFLIB_AUX
