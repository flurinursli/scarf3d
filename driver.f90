PROGRAM DRIVER

  USE :: MPI
  USE :: OMP_LIB

  USE, NON_INTRINSIC :: PRECISIONS
  USE, NON_INTRINSIC :: SCARFLIB

  IMPLICIT NONE

  INTEGER(ISP)               :: I, I0, I1, J0, J1, K0, K1
  INTEGER(ISP)               :: NX, NY, NZ

  ! MPI STUFF
  INTEGER(ISP)               :: MCW, IERR, NTASKS, RANK, TOPOZ, NDIMS
  INTEGER(ISP)               :: DIR, DISPL, SRC, DEST
  INTEGER(ISP)               :: NONEIGHBOR
  INTEGER(ISP), DIMENSION(3) :: DIMS, COORDS

  LOGICAL                    :: REORDER
  LOGICAL, DIMENSION(3)      :: ISPERIODIC


  !--------------------------------------------------------------------------------------------------------------------------------

  ! INITIALISE MPI
  CALL MPI_INIT(IERR)

  ! DUPLICATE COMMUNICATOR (I.E. MPI_COMM_WORLD INTO MCW)
  CALL MPI_COMM_DUP(MPI_COMM_WORLD, MCW, IERR)

  ! GET RANK NUMBER
  CALL MPI_COMM_RANK(MCW, RANK, IERR)

  ! GET NUMBER OF TASKS
  CALL MPI_COMM_SIZE(MCW, NTASKS, IERR)

  ! SET BARRIER
  CALL MPI_BARRIER(MCW, IERR)


  NX = 11
  NY = 20
  NZ = 5

  CALL SCARF3D(MCW, NX, NY, NZ, 10)

  ! WE WORKING IN 3D
  ! NDIMS = 3
  !
  ! ! NUMBER OF BLOCKS ALONG EACH DIRECTION
  ! DIMS  = [2, 1, 2]
  !
  ! ! ALLOW REORDERING
  ! REORDER = .TRUE.
  !
  ! ! NO PERIODICITY
  ! ISPERIODIC = [.FALSE., .FALSE., .FALSE.]
  !
  ! ! CREATE TOPOLOGY
  ! CALL MPI_CART_CREATE(MCW, NDIMS, DIMS, ISPERIODIC, REORDER, TOPOZ, IERR)

  ! DO I = 0, NTASKS - 1
  !
  !   ! PRINT COORDINATES OF CALLING PROCESS
  !   IF (RANK .EQ. I) THEN
  !     CALL MPI_CART_GET(TOPOZ, NDIMS, DIMS, ISPERIODIC, COORDS, IERR)
  !
  !     PRINT*, 'PROC ', RANK, ': ', COORDS
  !
  !     CALL MPI_SPLIT_TASK(NX, DIMS(1), COORDS(1), I0, I1)
  !     CALL MPI_SPLIT_TASK(NY, DIMS(2), COORDS(2), J0, J1)
  !     CALL MPI_SPLIT_TASK(NZ, DIMS(3), COORDS(3), K0, K1)
  !
  !     PRINT*, 'INDICES: ', 'X= ', I0, I1, ' Y= ', J0, J1, ' Z= ', K0, K1
  !
  !     DIR   = 0
  !     DISPL = 1
  !
  !     CALL MPI_CART_SHIFT(TOPOZ, DIR, DISPL, SRC, DEST, IERR)
  !
  !     PRINT*, 'X-DIR: ', 'LEFT= ', SRC, ' MID= ', RANK, ' RIGHT= ', DEST
  !
  !     DIR   = 1
  !     DISPL = 1
  !
  !     CALL MPI_CART_SHIFT(TOPOZ, DIR, DISPL, SRC, DEST, IERR)
  !
  !     PRINT*, 'Y-DIR: ', 'LEFT= ', SRC, ' MID= ', RANK, ' RIGHT= ', DEST
  !
  !     DIR   = 2
  !     DISPL = 1
  !
  !     CALL MPI_CART_SHIFT(TOPOZ, DIR, DISPL, SRC, DEST, IERR)
  !
  !     PRINT*, 'Z-DIR: ', 'LEFT= ', SRC, ' MID= ', RANK, ' RIGHT= ', DEST
  !
  !   ENDIF
  !
  ! ENDDO



  CALL MPI_BARRIER(MCW, IERR)

  CALL MPI_FINALIZE(IERR)

END PROGRAM DRIVER
