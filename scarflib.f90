MODULE SCARFLIB

  USE                :: MPI
  USE                :: OMP_LIB

  USE, NON_INTRINSIC :: PRECISIONS

  IMPLICIT NONE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


  PRIVATE

  PUBLIC :: SCARF3D

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  ! INTERFACE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  ! VARIABLES
  REAL(FPP), PARAMETER :: BIT2GB = 1._FPP / 8._FPP / 1024._FPP / 1024._FPP / 1024._FPP


  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE SCARF3D(MCW, NX, NY, NZ, MAXMEM)

      USE DECOMP_2D
      USE DECOMP_2D_FFT

      INTEGER(ISP), INTENT(IN) :: MCW                              !< MPI COMMUNICATOR
      INTEGER(ISP), INTENT(IN) :: NX, NY, NZ                       !< NUMBER OF POINTS ALONG X, Y, Z
      INTEGER(ISP), INTENT(IN) :: MAXMEM                           !< MAXIMUM AVAILABLE MEMORY PER NODE

      INTEGER(ISP)               :: RANK, NTASKS, IERR
      INTEGER(ISP)               :: XTOPO, YTOPO, ZTOPO
      INTEGER(ISP)               :: I0, I1, J0, J1, K0, K1
      INTEGER(ISP)               :: I, J, JSTAR, K
      INTEGER(ISP)               :: P_ROW, P_COL
      INTEGER(ISP)               :: SENDER, RECEIVER
      INTEGER(ISP), DIMENSION(3) :: NYQUIST
      INTEGER(ISP), DIMENSION(3) :: FFT_START, FFT_END, FFT_SIZE
      INTEGER(ISP), DIMENSION(3) :: COORDS
      INTEGER(ISP), DIMENSION(NY) :: LIST

      COMPLEX(RDP),              DIMENSION(NZ)    :: STRIPE
      COMPLEX(RDP), ALLOCATABLE, DIMENSION(:,:,:) :: SPEC

      REAL(RDP)                                   :: SCALING

      REAL(RDP),    ALLOCATABLE, DIMENSION(:,:,:) :: FIELD

      !-----------------------------------------------------------------------------------------------------------------------------

      ! GET RANK NUMBER
      CALL MPI_COMM_RANK(MCW, RANK, IERR)

      ! GET AVAILABLE MPI PROCESSES
      CALL MPI_COMM_SIZE(MCW, NTASKS, IERR)

      CALL MPI_BARRIER(MCW, IERR)

      ! INDEX OF NYQUIST FREQUENCY
      NYQUIST = [NX, NY, NZ] / 2 + 1


      P_ROW = 0
      P_COL = 0

      ! DECOMPOSE DOMAIN
      CALL DECOMP_2D_INIT(NX, NY, NZ, P_ROW, P_COL)

      ! INITIALISE FFT ENGINE
      CALL DECOMP_2D_FFT_INIT

      ! GET COMPLEX ARRAY SIZE FOR EACH MPI PROCESS
      CALL DECOMP_2D_FFT_GET_SIZE(FFT_START, FFT_END, FFT_SIZE)

      PRINT*, RANK, ' X ', FFT_START(1), FFT_END(1), ' -- ', FFT_START(2), FFT_END(2), ' -- ', FFT_START(3), FFT_END(3)

!call MPI_ALLGATHER(i_size,1,MPI_INTEGER,counts,1,MPI_INTEGER,MCW,err)

      ! ALLOCATE ARRAY FOR SPECTRUM (STENCIL ALONG Z-AXIS)
      ALLOCATE(SPEC(FFT_START(1):FFT_END(1), FFT_START(2):FFT_END(2), FFT_START(3):FFT_END(3)))

      ! COMPUTE SPECTRUM
      DO K = FFT_START(3), FFT_END(3)
        DO J = FFT_START(2), FFT_END(2)
          DO I = FFT_START(1), FFT_END(1)
            !SPEC(I, J, K) =
          ENDDO
        ENDDO
      ENDDO

      IF ( (I .EQ. 1) .AND. (J .EQ. 1) .AND. (K .EQ. 1) ) SPEC(I, J, K) = 0.               !< ZERO MEAN

      ! MAKE VALUES AT ALL NYQUIST REAL

      ! IF ( (I .EQ. 1) .AND. (J .EQ. 1) ) THEN
      !
      !   DO K = 2, N(3)/2
      !      SPEC(I, J, N(3) - K + 2) = CONJG(SPEC(I, J, K))
      !   ENDDO

      ! AT THIS POINT, "SPEC" AT I = 1 AND I = NYQUIST(1) IS GOOD BETWEEN J=1,NYQUIST(2) AND K=1,NYQUIST(3) (SEE FIG.4 PARDO)

      LIST = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]

      DO I = FFT_START(1), FFT_END(1)

        IF ( (I .EQ. 1) .OR. (I .EQ. NYQUIST(1)) ) THEN

          DO J = FFT_START(2), FFT_END(2)

            IF ( (J .EQ. 1) .OR. (J .EQ. NYQUIST(2)) ) THEN

              DO K = 2, NZ/2
                SPEC(I, J, NZ - K + 2) = CONJG(SPEC(I, J, K))                 !< BETA* = CONJG(BETA), THETA* = CONJG(THETA)
              ENDDO

              SPEC(I, J, 1)          = REAL(SPEC(I, J, 1), RDP)               !< ETA MUST BE REAL (ALPHA IS ALREADY ZERO)
              SPEC(I, J, NYQUIST(3)) = REAL(SPEC(I, J, NYQUIST(3)), RDP)      !< XSI AND PSI MUST BE REAL

            ELSEIF ( ( J .GE. 2) .AND. (J .LT. NYQUIST(2)) ) THEN

              ! "JSTAR" MUST CONTAIN THE CONJG OF "J"
              JSTAR = NY - J + 2

              ! PROCESS WITH INDEX "J" MUST SEND WHOLE STRIPE (DELTA, EPSILON, PHI, GAMMA) TO PROCESS "JSTAR"
              !STRIPE = SPEC(I, J, :)
              STRIPE = [CMPLX(RANK, 1), CMPLX(RANK, 2), CMPLX(RANK, 3), CMPLX(RANK, 4), CMPLX(RANK, 5)]

              SENDER   = LIST(J)
              RECEIVER = LIST(JSTAR)

              !PRINT*, 'SR ',  RANK, SENDER, RECEIVER

              IF (SENDER .EQ. RANK)   CALL MPI_SEND(STRIPE, NZ, MPI_DOUBLE_COMPLEX, RECEIVER, 0, MCW, IERR)
              IF (RECEIVER .EQ. RANK) CALL MPI_RECV(STRIPE, NZ, MPI_DOUBLE_COMPLEX, SENDER, 0, MCW, MPI_STATUS_IGNORE, IERR)

              ! AT THIS POINT STRIPE CONTAINS VALUES FROM "SENDER"

              PRINT*, RANK, STRIPE

              SPEC(I, JSTAR, 1) = CONJG(STRIPE(1))

              DO K = 2, NZ/2 + 1
                SPEC(I, JSTAR, NZ - K + 2) = CONJG(STRIPE(K))
              ENDDO

  ! DOWN HERE: 30.10

            ENDIF

          ENDDO

        ENDIF

      ENDDO

              ! do k=2,N(3)/2
              !    do j=2,N(2)/2
              !       SPEC(N(2)-j+2,N(3)-k+2) = conjg(SPEC(j,k))



      ! ALLOCATE OUTPUT ARRAY (STENCIL ALONG X-AXIS)
      ALLOCATE(FIELD(XSTART(1):XEND(1), XSTART(2):XEND(2), XSTART(3):XEND(3)))

      ! ONCE SPECTRUM IS READY, CALL COMPLEX TO REAL TRANSFORM
      CALL DECOMP_2D_FFT_3D(SPEC, FIELD)

      DO K = XSTART(3), XEND(3)
        DO J = XSTART(2), XEND(2)
          DO I = XSTART(1), XEND(1)
            FIELD(I, J, K) = FIELD(I, J, K) / SCALING
          ENDDO
        ENDDO
      ENDDO

      DO I = 0, NTASKS-1
        IF (RANK .EQ. I)   &
          PRINT*, 'RANK= ', RANK , '| X= ', ZSTART(1), ZEND(1), '| Y= ', ZSTART(2), ZEND(2), '| Z= ', ZSTART(3), ZEND(3)
          !PRINT*, 'RANK= ', RANK , '| X= ', YSTART(1), YEND(1), '| Y= ', YSTART(2), YEND(2), '| Z= ', YSTART(3), YEND(3)
      ENDDO

      !IF (RANK .EQ. 0) PRINT*, STORAGE_SIZE(CMPLX(0., 0.), KIND=C_DOUBLE_COMPLEX) * 1000000000 * BIT2GB

      CALL DECOMP_2D_FFT_FINALIZE

      CALL DECOMP_2D_FINALIZE

    END SUBROUTINE SCARF3D

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *




    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE MPI_SPLIT_TASK(N, TASKS, RANK, I0, I1)

      INTEGER(ISP), INTENT(IN)  :: N
      INTEGER(ISP), INTENT(IN)  :: TASKS, RANK
      INTEGER(ISP), INTENT(OUT) :: I0, I1

      !------------------------------------------------------------------------------------------------------------------------------

      I0 = 1 + INT( REAL(N) / REAL(TASKS) * REAL(RANK) )
      I1 = INT( REAL(N) / REAL(TASKS) * REAL(RANK + 1) )

    END SUBROUTINE MPI_SPLIT_TASK

! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
!===============================================================================================================================
! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *






END MODULE SCARFLIB
