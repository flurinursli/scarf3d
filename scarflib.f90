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

  ! PROCEDURE POINTERS
  PROCEDURE(VK_PSDF), POINTER :: FUN

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  ! VARIABLES
  REAL(RDP), PARAMETER :: PI = 3.141592653589793_RDP

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE SCARF3D(MCW, N, DH, ACF, CL, SIGMA, HURST, SEED, POI, MUTE, TAPER)

      USE DECOMP_2D
      USE DECOMP_2D_FFT
      USE DECOMP_2D_IO

      INTEGER(ISP),                                   INTENT(IN) :: MCW                                         !< MPI COMMUNICATOR
      INTEGER(ISP),                  DIMENSION(3),    INTENT(IN) :: N                                           !< NUMBER OF POINTS ALONG X, Y, Z
      REAL(RDP),                                      INTENT(IN) :: DH
      CHARACTER(LEN=*),                               INTENT(IN) :: ACF
      REAL(RDP),                     DIMENSION(3),    INTENT(IN) :: CL
      REAL(RDP),                                      INTENT(IN) :: SIGMA
      REAL(RDP),                                      INTENT(IN) :: HURST
      INTEGER(ISP),                                   INTENT(IN) :: SEED
      INTEGER(ISP),                  DIMENSION(:,:),  INTENT(IN) :: POI                                         !< LOCATION OF POINT(S)-OF-INTEREST
      INTEGER(ISP),                                   INTENT(IN) :: MUTE                                        !< NUMBER OF POINTS WHERE MUTING IS APPLIED
      INTEGER(ISP),                                   INTENT(IN) :: TAPER                                       !< NUMBER OF POINTS WHERE TAPERING IS APPLIED
      COMPLEX(RDP),     ALLOCATABLE, DIMENSION(:,:,:)            :: SPEC
      INTEGER(ISP)                                               :: RANK, NTASKS, IERR
      INTEGER(ISP)                                               :: I, J, K
      INTEGER(ISP)                                               :: P_ROW, P_COL
      INTEGER(ISP),                  DIMENSION(3)                :: NYQUIST
      INTEGER(ISP),                  DIMENSION(3)                :: FFT_START, FFT_END, FFT_SIZE
      INTEGER(ISP),     ALLOCATABLE, DIMENSION(:)                :: NPTS
      REAL(RDP)                                                  :: KNYQ, KC, KR
      REAL(RDP)                                                  :: BUTTER, NUM, AMP
      REAL(RDP)                                                  :: SCALING, MU, VAR
      REAL(RDP)                                                  :: TIC, TAC
      REAL(RDP),                     DIMENSION(3)                :: DK
      REAL(RDP),        ALLOCATABLE, DIMENSION(:)                :: KX, KY, KZ
      REAL(RDP),        ALLOCATABLE, DIMENSION(:)                :: AVERAGE, STD_DEV
      REAL(RDP),        ALLOCATABLE, DIMENSION(:,:,:)            :: FIELD, R, XCORR

      !-----------------------------------------------------------------------------------------------------------------------------

      ! GET RANK NUMBER
      CALL MPI_COMM_RANK(MCW, RANK, IERR)

      ! GET AVAILABLE MPI PROCESSES
      CALL MPI_COMM_SIZE(MCW, NTASKS, IERR)

      CALL MPI_BARRIER(MCW, IERR)

      ! INDEX OF NYQUIST FREQUENCY
      NYQUIST = N / 2 + 1

      P_ROW = 0 !1
      P_COL = 0  !3

      TIC = MPI_WTIME()

      ! DECOMPOSE DOMAIN
      CALL DECOMP_2D_INIT(N(1), N(2), N(3), P_ROW, P_COL)

      ! IF (RANK .EQ. 0) PRINT*, 'PENCILS ALONG Z-AXIS'
      !
      ! CALL MPI_BARRIER(MCW, IERR)
      !
      ! DO I = 0, NTASKS-1
      !   IF (RANK .EQ. I)   &
      !     PRINT*, 'RANK= ', RANK , '| X= ', ZSTART(1), ZEND(1), '| Y= ', ZSTART(2), ZEND(2), '| Z= ', ZSTART(3), ZEND(3)
      ! ENDDO
      !
      ! CALL MPI_BARRIER(MCW, IERR)

      ! INITIALISE FFT ENGINE
      CALL DECOMP_2D_FFT_INIT

      ! GET COMPLEX ARRAY SIZE FOR EACH MPI PROCESS
      CALL DECOMP_2D_FFT_GET_SIZE(FFT_START, FFT_END, FFT_SIZE)

      CALL MPI_BARRIER(MCW, IERR)

      TAC = MPI_WTIME()

      IF (RANK .EQ. 0) PRINT*, 'INITIALISATION TIME: ', REAL(TAC - TIC, RSP)


    !
    !   IF (RANK .EQ. 0) PRINT*, 'FFT PENCILS'
    !
    !
    !   DO I = 0, NTASKS-1
    !     IF (RANK .EQ. I)   &
    ! PRINT*, 'RANK= ', RANK, '| X= ', FFT_START(1), FFT_END(1), '| Y= ', FFT_START(2), FFT_END(2), '| Z= ', FFT_START(3), FFT_END(3)
    !   ENDDO

      TIC = MPI_WTIME()

      ALLOCATE(R(FFT_START(1):FFT_END(1), FFT_START(2):FFT_END(2), FFT_START(3):FFT_END(3)))

      ! GENERATE A SET OF RANDOM NUMBERS IN THE RANGE [0, 2*PI)
      CALL RANDOM_SEQUENCE(SEED, FFT_START, FFT_END, N, R)

      CALL MPI_BARRIER(MCW, IERR)

      TAC = MPI_WTIME()

      IF (RANK .EQ. 0) PRINT*, 'RANDOM SEQUENCE TIMING: ', REAL(TAC - TIC, RSP)

      ! ALLOCATE ARRAY FOR SPECTRUM (STENCILS ORIENTED ALONG Z-AXIS)
      ALLOCATE(SPEC(FFT_START(1):FFT_END(1), FFT_START(2):FFT_END(2), FFT_START(3):FFT_END(3)))

      ! RESOLUTION IN WAVENUMBER DOMAIN ALONG EACH DIRECTION
      DO I = 1, 3
        DK(I) = 2._RDP * PI / (REAL(N(I), RDP) * DH)
      ENDDO

      ! NYQUIST WAVENUMBER
      KNYQ = PI / DH

      ! CORNER WAVENUMBER FOR SPECTRUM TAPERING TO AVOID ALIASING IS HALF NYQUIST WAVENUMBER
      KC = KNYQ / 2._RDP

      ! WAVENUMBER VECTORS
      ALLOCATE(KX(N(1)), KY(N(2)), KZ(N(3)))

      ! VECTORS GO FROM 0 TO NYQUIST AND THEN BACK AGAIN UNTIL DK
      KX = [[(I * DK(1), I = 0, N(1)/2)], [(I * DK(1), I = N(1)/2-1, 1, -1)]]
      KY = [[(J * DK(2), J = 0, N(2)/2)], [(J * DK(2), J = N(2)/2-1, 1, -1)]]
      KZ = [[(K * DK(3), K = 0, N(3)/2)], [(K * DK(3), K = N(3)/2-1, 1, -1)]]

      ! COMPUTE PART OF POWER SPECTRAL DENSITY OUTSIDE LOOP
      IF (ACF .EQ. 'VK') THEN
        NUM = 8._RDP * SQRT(PI**3) * GAMMA(HURST + 1.5_RDP) * SIGMA**2 * PRODUCT(CL) / GAMMA(HURST)
        FUN => VK_PSDF
      ELSEIF (ACF .EQ. 'GAUSS') THEN
        NUM = SIGMA**2 * PRODUCT(CL) * SQRT(PI**3)
        FUN => GS_PSDF
      ENDIF

      ! SCALING PARAMETER
      SCALING = 1._RDP / SQRT(PRODUCT(N) * DH**3)

      ! COMPUTE SPECTRUM
      DO K = FFT_START(3), FFT_END(3)
        DO J = FFT_START(2), FFT_END(2)
          DO I = FFT_START(1), FFT_END(1)

            ! RADIAL WAVENUMBER
            !KR = SQRT(KX(I)**2 + KY(J)**2 + KZ(K)**2)

            ! FOURTH-ORDER LOW-PASS BUTTERWORTH FILTER
            !BUTTER = 1._RDP / SQRT(1._RDP + (KR / KC)**(2 * 4))

            ! NOW "KR" IS THE PRODUCT "K * CL"
            KR = (KX(I) * CL(1))**2 + (KY(J) * CL(2))**2 + (KZ(K) * CL(3))**2

            ! COMPLETE POWER SPECTRAL DENSITY AND GO TO AMPLITUDE SPECTRUM
            AMP = SQRT(NUM / FUN(KR, HURST))

            ! APPLY FILTER
            !AMP = AMP * BUTTER

            ! COMBINE AMPLITUDE AND RANDOM PHASE
            SPEC(I, J, K) = CMPLX(AMP * COS(R(I, J, K)), AMP * SIN(R(I, J, K)))

          ENDDO
        ENDDO
      ENDDO

      CALL MPI_BARRIER(MCW, IERR)

      TIC = MPI_WTIME()

      IF (RANK .EQ. 0) PRINT*, 'SPECTRUM TIMING: ', REAL(TIC - TAC, RSP)

      ! ENFORCE SYMMETRY CONDITIONS AT Y-Z PLANES FOR "FFT_START(1)=1" AND "FFT_END(1)=NYQUIST(1)"
      CALL ENFORCE_SYMMETRY(MCW, FFT_START, FFT_END, NYQUIST, SPEC)

      CALL MPI_BARRIER(MCW, IERR)

      TAC = MPI_WTIME()

      IF (RANK .EQ. 0) PRINT*, 'SYMMETRY TIMING: ', REAL(TAC - TIC, RSP)

      ! ALLOCATE OUTPUT ARRAY (STENCILS ORIENTED ALONG X-AXIS)
      ALLOCATE(FIELD(XSTART(1):XEND(1), XSTART(2):XEND(2), XSTART(3):XEND(3)))

      ! ONCE SPECTRUM IS READY, CALL COMPLEX TO REAL TRANSFORM
      CALL DECOMP_2D_FFT_3D(SPEC, FIELD)

      CALL MPI_BARRIER(MCW, IERR)

      TIC = MPI_WTIME()

      IF (RANK .EQ. 0) PRINT*, 'IFFT TIMING: ', REAL(TIC - TAC, RSP)

      CALL DECOMP_2D_FFT_FINALIZE

      ! APPLY SCALING FACTOR
      DO K = XSTART(3), XEND(3)
        DO J = XSTART(2), XEND(2)
          DO I = XSTART(1), XEND(1)
            FIELD(I, J, K) = FIELD(I, J, K) * SCALING
          ENDDO
        ENDDO
      ENDDO

      ALLOCATE(AVERAGE(NTASKS), STD_DEV(NTASKS), NPTS(NTASKS))

      VAR = VARIANCE(FIELD)
      MU  = MEAN(FIELD)

      CALL MPI_ALLGATHER(VAR, 1, MPI_DOUBLE_PRECISION, STD_DEV, 1, MPI_DOUBLE_PRECISION, MCW, IERR)
      CALL MPI_ALLGATHER(MU, 1, MPI_DOUBLE_PRECISION, AVERAGE, 1, MPI_DOUBLE_PRECISION, MCW, IERR)
      CALL MPI_ALLGATHER(PRODUCT(FFT_SIZE), 1, MPI_INTEGER, NPTS, 1, MPI_INTEGER, MCW, IERR)

      DO I = 1, NTASKS - 1
        CALL PARALLEL_VARIANCE(STD_DEV(I:I + 1), AVERAGE(I:I + 1), NPTS(I:I + 1))
      ENDDO

      IF (RANK .EQ. 0) PRINT*, 'MEAN, STD.DEV: ', AVERAGE(NTASKS), SQRT(STD_DEV(NTASKS))


      CALL MPI_BARRIER(MCW, IERR)

      TIC = MPI_WTIME()

      ! TAPER / MUTE AROUND POINT-OF-INTEREST
      DO I = 1, SIZE(POI, 2)
        CALL TAPERING(XSTART, XEND, FIELD, POI(:, I), MUTE, TAPER)
      ENDDO

      CALL MPI_BARRIER(MCW, IERR)

      TAC = MPI_WTIME()

      IF (RANK .EQ. 0) PRINT*, 'TAPERING TIMING: ', REAL(TAC - TIC, RSP)


!       IF (RANK .EQ. 0) PRINT*, 'AUTOCORRELATION'
!
!       ALLOCATE(XCORR(XSTART(1):XEND(1), XSTART(2):XEND(2), XSTART(3):XEND(3)))
!
!       ! AUTOCORRELATION ALONG X-AXIS
!       CALL AUTOCORRELATION(1, FIELD, XCORR)
!
!       CALL DECOMP_2D_WRITE_ONE(1, XCORR, 'xcorr.bin')


      CALL MPI_BARRIER(MCW, IERR)

      TIC = MPI_WTIME()

      ! WRITE TO DISK
      CALL DECOMP_2D_WRITE_ONE(1, FIELD, 'random.bin')

      CALL DECOMP_2D_WRITE_PLANE(1, FIELD, 1, N(1)/2, 'slice.bin')

      CALL MPI_BARRIER(MCW, IERR)

      TAC = MPI_WTIME()

      IF (RANK .EQ. 0) PRINT*, 'IO TIMING: ', REAL(TAC - TIC, RSP)

      ! CLEAN UP MEMORY
      DEALLOCATE(SPEC, FIELD)

      CALL DECOMP_2D_FINALIZE

    END SUBROUTINE SCARF3D

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE ENFORCE_SYMMETRY(COMM, FS, FE, NYQUIST, SPEC)

      INTEGER(ISP),                                                              INTENT(IN)    :: COMM                              !< MPI COMMUNICATOR
      INTEGER(ISP),              DIMENSION(3),                                   INTENT(IN)    :: FS, FE
      INTEGER(ISP),              DIMENSION(3),                                   INTENT(IN)    :: NYQUIST
      COMPLEX(RDP),              DIMENSION(FS(1):FE(1),FS(2):FE(2),FS(3):FE(3)), INTENT(INOUT) :: SPEC
      COMPLEX(RDP), ALLOCATABLE, DIMENSION(:)                                                  :: STRIPE
      INTEGER(ISP)                                                                             :: I, J, JSTAR, K                    !< COUNTERS
      INTEGER(ISP)                                                                             :: J0, J1, STEP                      !< COUNTERS
      INTEGER(ISP)                                                                             :: RANK, SENDER, RECEIVER, NTASKS    !< MPI STUFF
      INTEGER(ISP)                                                                             :: IERR
      INTEGER(ISP),              DIMENSION(3)                                                  :: N                                 !< POINTS ALONG EACH DIMENSION
      INTEGER(ISP), ALLOCATABLE, DIMENSION(:)                                                  :: LIST
      INTEGER(ISP), ALLOCATABLE, DIMENSION(:,:)                                                :: LBOUND, UBOUND

      !-----------------------------------------------------------------------------------------------------------------------------

      ! POINTS ALONG EACH DIMENSION
      DO I = 1, 3
        N(I) = 2 * (NYQUIST(I) - 1)
      ENDDO

      ! GET RANK NUMBER
      CALL MPI_COMM_RANK(COMM, RANK, IERR)

      ! GET AVAILABLE MPI PROCESSES
      CALL MPI_COMM_SIZE(COMM, NTASKS, IERR)

      ALLOCATE(LBOUND(3, 0:NTASKS-1), UBOUND(3, 0:NTASKS-1))

      ! GATHER START/END FFT INDEX ALONG EACH DIMENSION FOR ALL MPI PROCESSES
      CALL MPI_ALLGATHER(FS, 3, MPI_INTEGER, LBOUND, 3, MPI_INTEGER, COMM, IERR)
      CALL MPI_ALLGATHER(FE, 3, MPI_INTEGER, UBOUND, 3, MPI_INTEGER, COMM, IERR)

      ! PRINT*, RANK, ' ++ ', UBOUND(2,:)

      ! ALLOCATE ARRAYS WITH AS MANY POINTS AS NY/NZ
      ALLOCATE(LIST(N(2)), STRIPE(N(3)))

      DO I = FS(1), FE(1)

        ! ENFORCE SYMMETRY AT X = 1 AND X = N(1)/2 + 1
        IF ( (I .EQ. 1) .OR. (I .EQ. NYQUIST(1)) ) THEN

          DO J = 0, NTASKS - 1
            IF ( (LBOUND(1, J) .EQ. I) .OR. (UBOUND(1, J) .EQ. I) ) THEN
              LIST(LBOUND(2, J):UBOUND(2, J)) = J
            ENDIF
          ENDDO

          ! PRINT*, RANK, '-- ', LIST

          ! ASSIGN LOOP INDICES FOR SENDING MPI TASKS
          J0   = FS(2)
          J1   = FE(2)
          STEP = 1

          ! SWAP INDICES TO CREATE RECEIVING MPI TASKS MATCHING THE SENDING ONES
          IF (FS(2) .GE. NYQUIST(2)) THEN
            J0   = FE(2)
            J1   = FS(2)
            STEP = -1
          ENDIF

          ! "J" IS VERTICAL AXIS IN FIGURE 4
          DO J = J0, J1, STEP

            IF ( (J .EQ. 1) .OR. (J .EQ. NYQUIST(2)) ) THEN

              ! "K" IN HORIZONTAL AXIS IN FIGURE 4
              DO K = 2, NYQUIST(3) - 1
                SPEC(I, J, N(3) - K + 2) = CONJG(SPEC(I, J, K))               !< BETA* = CONJG(BETA), THETA* = CONJG(THETA)
              ENDDO

              SPEC(I, J, 1)          = REAL(SPEC(I, J, 1), RDP)               !< ETA MUST BE REAL
              SPEC(I, J, NYQUIST(3)) = REAL(SPEC(I, J, NYQUIST(3)), RDP)      !< XSI AND PSI MUST BE REAL

              IF ( (I .EQ. 1) .AND. (J .EQ. 1) ) SPEC(I, J, 1) = 0._RDP       !< ALPHA MUST BE ZERO

            ELSE

              ! "JSTAR" MUST CONTAIN THE CONJG OF "J"
              JSTAR = N(2) - J + 2

              IF (J .LT. NYQUIST(2)) THEN
                SENDER   = LIST(J)
                RECEIVER = LIST(JSTAR)
              ELSE
                SENDER   = LIST(JSTAR)
                RECEIVER = LIST(J)
              ENDIF

              ! MPI PROCESS SENDING DATA DIFFERS MPI PROCESS RECEIVING DATA
              IF (SENDER .NE. RECEIVER) THEN

                ! TAKE WHOLE STRIPE ALONG "K"-AXIS
                STRIPE = SPEC(I, J, :)
                !STRIPE = [CMPLX(J, 1), CMPLX(J, 2), CMPLX(J, 3), CMPLX(J, 4), CMPLX(J, 5), CMPLX(J, 6)]

                !print*, rank, '--', sender, receiver, '--', i, '++', j, jstar

                IF (SENDER .EQ. RANK)   CALL MPI_SEND(STRIPE, N(3), MPI_DOUBLE_COMPLEX, RECEIVER, JSTAR, COMM, IERR)
                IF (RECEIVER .EQ. RANK) CALL MPI_RECV(STRIPE, N(3), MPI_DOUBLE_COMPLEX, SENDER, J, COMM, MPI_STATUS_IGNORE, IERR)

                JSTAR = J

              ! MPI PROCESSES ARE THE SAME
              ELSE

                ! NO NEED TO COMPUTE COMPLEX CONJUGATE IF WE ARE ABOVE NYQUIST INDEX (BECAUSE IT HAS BEEN COMPUTE ALREADY)
                IF (J .GT. NYQUIST(2)) CYCLE

                ! TAKE WHOLE STRIPE ALONG "K"-AXIS
                STRIPE = SPEC(I, J, :)
                !STRIPE = [CMPLX(J, 1), CMPLX(J, 2), CMPLX(J, 3), CMPLX(J, 4), CMPLX(J, 5), CMPLX(J, 6)]

              ENDIF

              ! UPDATE STRIPE ALONG "K"-AXIS WITH COMPLEX CONJUGATE
              IF (RANK .EQ. RECEIVER) THEN

                !PRINT*, RANK, SENDER, RECEIVER, JSTAR, nint(real(STRIPE(1))), ' - ', nint(aimag(stripe))

                SPEC(I, JSTAR, 1) = CONJG(STRIPE(1))                        !< DELTA* = CONJG(DELTA)

                DO K = 2, N(3)
                  SPEC(I, JSTAR, K) = CONJG(STRIPE(N(3) - K + 2))
                ENDDO

                ! IF (I .EQ. NYQUIST(1)) &
                ! PRINT*, RANK, SENDER, RECEIVER, JSTAR, ' - ', NINT(REAL(SPEC(I, JSTAR, 1))), ' - ', NINT(AIMAG(SPEC(I, JSTAR, :)))

              ENDIF

            ENDIF

          ENDDO

        ENDIF

      ENDDO

      ! REALEASE MEMORY
      DEALLOCATE(LIST, STRIPE, LBOUND, UBOUND)

    END SUBROUTINE ENFORCE_SYMMETRY

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE RANDOM_SEQUENCE(SEED, FS, FE, N, R)

      ! GENERATE A SEQUENCE OF RANDOM NUMBERS IN THE RANGE [0 2*PI)

      INTEGER(ISP),                                                              INTENT(IN)  :: SEED                 !< INITIAL SEED
      INTEGER(ISP),              DIMENSION(3),                                   INTENT(IN)  :: FS, FE               !< START/END INDICES ALONG EACH DIRECTION
      INTEGER(ISP),              DIMENSION(3),                                   INTENT(IN)  :: N                    !< SIZE OF WHOLE RANDOM FIELD (POINTS)
      REAL(RDP),                 DIMENSION(FS(1):FE(1),FS(2):FE(2),FS(3):FE(3)), INTENT(OUT) :: R                    !< RANDOM SEQUENCE
      INTEGER(ISP)                                                                           :: SEED_SIZE
      INTEGER(ISP)                                                                           :: I, J                 !< COUNTERS
      INTEGER(ISP), ALLOCATABLE, DIMENSION(:)                                                :: TMP_SEED
      REAL(RDP),                 DIMENSION(FS(3):FE(3))                                      :: HARVEST

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL RANDOM_SEED(SIZE = SEED_SIZE)

      ALLOCATE(TMP_SEED(SEED_SIZE))

      DO I = FS(1), FE(1)
        DO J = FS(2), FE(2)

          ! EACH LONG-SIDE OF A PENCIL ORIENTED ALONG Z-AXIS HAS ITS OWN UNIQUE SEED NUMBER. THIS GUARANTEES IDENTICAL RESULTS EVEN
          ! IF THE TOTAL NUMBER OF MPI TASKS OR THE NUMBER OF TASKS ALONG X AND Y CHANGES.
          TMP_SEED = SEED + J + (I - 1) * N(1)

          CALL RANDOM_SEED(PUT = TMP_SEED)

          CALL RANDOM_NUMBER(HARVEST)

          R(I, J, :) = HARVEST * 2._RDP * PI                               !< SEQUENCE IN THE RANGE [0, 2*PI)

        ENDDO
      ENDDO

      DEALLOCATE(TMP_SEED)

    END SUBROUTINE RANDOM_SEQUENCE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(RDP) FUNCTION VK_PSDF(M, HURST)

      REAL(RDP), INTENT(IN) :: M
      REAL(RDP), INTENT(IN) :: HURST

      !-----------------------------------------------------------------------------------------------------------------------------

      VK_PSDF = (1._RDP + M)**(HURST + 1.5_RDP)

    END FUNCTION VK_PSDF

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(RDP) FUNCTION GS_PSDF(M, HURST)

      REAL(RDP), INTENT(IN) :: M
      REAL(RDP), INTENT(IN) :: HURST

      !-----------------------------------------------------------------------------------------------------------------------------

      GS_PSDF = EXP(0.25_RDP * M)

    END FUNCTION GS_PSDF

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION COMPENSATED_SUM(R) RESULT(V)

      ! SUM ELEMENTS OF ARRAY BASED ON KAHAN-BABUSKA ALGORITHM

      REAL(RDP),   DIMENSION(:,:,:), INTENT(IN) :: R
      INTEGER(ISP)                              :: I, J, K
      REAL(RDP)                                 :: V, C, T

      !-------------------------------------------------------------------------------------------------------------------------------

      V = 0._RDP
      C = 0._RDP

      DO K = 1, SIZE(R, 3)
        DO J = 1, SIZE(R, 2)
          DO I = 1, SIZE(R, 1)

            T = V + R(I, J, K)

            IF (ABS(V) .GE. ABS(R(I, J, K))) THEN
              C = C + (V - T) + R(I, J, K)
            ELSE
              C = C + (R(I, J, K) - T) + V
            ENDIF

            V = T

          ENDDO
        ENDDO
      ENDDO

      V = V + C

    END FUNCTION COMPENSATED_SUM

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION MEAN(R) RESULT(V)

      REAL(RDP),   DIMENSION(:,:,:), INTENT(IN) :: R
      INTEGER(ISP)                              :: I, J, K
      REAL(RDP)                                 :: V, C

      !-------------------------------------------------------------------------------------------------------------------------------

      !V = COMPENSATED_SUM(R) / REAL(SIZE(R), RDP)

      V = 0._RDP
      C = 1._RDP

      DO K = 1, SIZE(R, 3)
        DO J = 1, SIZE(R, 2)
          DO I = 1, SIZE(R, 1)
            V = V + (R(I, J, K) - V) / C
            C = C + 1._RDP
          ENDDO
        ENDDO
      ENDDO

    END FUNCTION MEAN

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION VARIANCE(R) RESULT(V)

      ! COMPUTE VARIANCE BASED ON THE COMPENSATED-SUMMATION VERSION OF THE TWO-PASS ALGORITHM

      REAL(RDP),   DIMENSION(:,:,:), INTENT(IN) :: R
      INTEGER(ISP)                              :: I, J, K
      REAL(RDP)                                 :: MU, S1, S2, X, V

      !-----------------------------------------------------------------------------------------------------------------------------

      MU = MEAN(R)

      S1 = 0._RDP
      S2 = 0._RDP

      DO K = 1, SIZE(R, 3)
        DO J = 1, SIZE(R, 2)
          DO I = 1, SIZE(R, 1)
            X = R(I, J, K) - MU
            S1 = S1 + X
            S2 = S2 + X**2
          ENDDO
        ENDDO
      ENDDO

      S1 = (S1**2) / REAL(SIZE(R), RDP)

      V = (S2 - S1) / REAL(SIZE(R) - 1, RDP)

    END FUNCTION VARIANCE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE PARALLEL_VARIANCE(VAR, MU, N)

      REAL(RDP),    DIMENSION(2), INTENT(INOUT) :: VAR
      REAL(RDP),    DIMENSION(2), INTENT(INOUT) :: MU
      INTEGER(ISP), DIMENSION(2), INTENT(INOUT) :: N
      REAL(RDP)                                 :: DELTA, M2
      REAL(RDP),    DIMENSION(2)                :: M

      !-----------------------------------------------------------------------------------------------------------------------------

      DELTA = MU(2) - MU(1)

      M(1) = VAR(1) * REAL(N(1) - 1, RDP)
      M(2) = VAR(2) * REAL(N(2) - 1, RDP)

      M2 = SUM(M) + DELTA**2 * REAL(PRODUCT(N), RDP) / REAL(SUM(N), RDP)

      M2 = M2 / REAL(SUM(N) - 1, RDP)

      ! OUTPUT MEAN
      !MU(2) = MU(1) + DELTA * REAL(N(2), RDP) / REAL(SUM(N), RDP)
      MU(2) = (N(1) * MU(1) + N(2) * MU(2)) / REAL(SUM(N), RDP)          !< THIS SHOULD BE MORE STABLE WHEN N(1)~N(2)

      ! OUTPUT VARIANCE
      VAR(2) = M2

      ! UPDATE TOTAL NUMBER OF POINTS
      N(2) = SUM(N)

    END SUBROUTINE PARALLEL_VARIANCE

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

    SUBROUTINE AUTOCORRELATION(AXIS, FIELD, CORR)

      USE, INTRINSIC                                                   :: ISO_C_BINDING

      INCLUDE 'fftw3.f03'

      INTEGER(ISP),                                        INTENT(IN)  :: AXIS
      REAL(RDP),                 DIMENSION(:,:,:),         INTENT(IN)  :: FIELD
      REAL(RDP),                 DIMENSION(:,:,:),         INTENT(OUT) :: CORR
      INTEGER(ISP)                                                     :: I, J, K, M
      INTEGER(ISP)                                                     :: NX, NY, NZ, NPTS
      COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:),     POINTER             :: TU => NULL()
      REAL(C_DOUBLE),            DIMENSION(:),     POINTER             :: U => NULL()
      TYPE(C_PTR)                                                      :: P1, P2, P

      !-----------------------------------------------------------------------------------------------------------------------------

      NX = SIZE(FIELD, 1)
      NY = SIZE(FIELD, 2)
      NZ = SIZE(FIELD, 3)

      IF (AXIS .EQ. 1) NPTS = NX                                      !< AUTOCORRELATION ALONG X-DIRECTION
      IF (AXIS .EQ. 2) NPTS = NY                                      !< AUTOCORRELATION ALONG Y-DIRECTION
      IF (AXIS .EQ. 3) NPTS = NZ                                      !< AUTOCORRELATION ALONG Z-DIRECTION

      ! NUMBER OF POINTS IN CORRELATION (INCL. NEGATIVE AND POSITIVE LAGS)
      NPTS = 2 * NPTS - 1

      M = NPTS / 2 + 1

      ! ALLOCATE CONTIGUOUS MEMORY
      P = FFTW_ALLOC_COMPLEX(INT(M, C_SIZE_T))

      CALL C_F_POINTER(P, U, [2*M])
      CALL C_F_POINTER(P, TU, [M])

      P1 = FFTW_PLAN_DFT_R2C_1D(NPTS, U, TU, FFTW_MEASURE)
      P2 = FFTW_PLAN_DFT_C2R_1D(NPTS, TU, U, FFTW_MEASURE)

      IF (AXIS .EQ. 1) THEN

        DO K = 1, NZ
          DO J = 1, NY

            U(:) = 0._RDP

            ! COPY INPUT VALUES INTO CONTIGUOUS MEMORY
            DO I = 1, NX
              U(I) = FIELD(I, J, K)
            ENDDO

            CALL FFTW_EXECUTE_DFT_R2C(P1, U, TU)

            DO I = 1, M
              TU(I) = TU(I) * CONJG(TU(I))
            ENDDO

            CALL FFTW_EXECUTE_DFT_C2R(P2, TU, U)

            ! SCALING IFFT
            DO I = 1, NPTS
              U(I) = U(I) / REAL(NPTS, RDP)
            ENDDO

            ! COMPUTE BIASED AUTOCORRELATION, RETURN ONLY ZERO AND POSITIVE LAGS
            DO I = 1, M
              CORR(I, J, K) = U(I) / REAL(M, RDP)
            ENDDO

          ENDDO
        ENDDO

      ELSEIF (AXIS .EQ. 2) THEN

        DO K = 1, NZ
          DO I = 1, NX

            U(:) = 0._RDP

            ! COPY INPUT VALUES INTO CONTIGUOUS MEMORY
            DO J = 1, NY
              U(J) = FIELD(I, J, K)
            ENDDO

            CALL FFTW_EXECUTE_DFT_R2C(P1, U, TU)

            DO J = 1, M
              TU(J) = TU(J) * CONJG(TU(J))
            ENDDO

            CALL FFTW_EXECUTE_DFT_C2R(P2, TU, U)

            ! SCALING IFFT
            DO J = 1, NPTS
              U(J) = U(J) / REAL(NPTS, RDP)
            ENDDO

            ! COMPUTE BIASED AUTOCORRELATION, RETURN ONLY ZERO AND POSITIVE LAGS
            DO J = 1, M
              CORR(I, J, K) = U(J) / REAL(M, RDP)
            ENDDO

          ENDDO
        ENDDO

      ELSEIF (AXIS .EQ. 3) THEN

        DO J = 1, NY
          DO I = 1, NX

            U(:) = 0._RDP

            ! COPY INPUT VALUES INTO CONTIGUOUS MEMORY
            DO K = 1, NZ
              U(K) = FIELD(I, J, K)
            ENDDO

            CALL FFTW_EXECUTE_DFT_R2C(P1, U, TU)

            DO K = 1, M
              TU(K) = TU(K) * CONJG(TU(K))
            ENDDO

            CALL FFTW_EXECUTE_DFT_C2R(P2, TU, U)

            ! SCALING IFFT
            DO K = 1, NPTS
              U(K) = U(K) / REAL(NPTS, RDP)
            ENDDO

            ! COMPUTE BIASED AUTOCORRELATION, RETURN ONLY ZERO AND POSITIVE LAGS
            DO K = 1, M
              CORR(I, J, K) = U(K) / REAL(M, RDP)
            ENDDO

          ENDDO
        ENDDO

      ENDIF

      ! REALEASE MEMORY
      CALL FFTW_DESTROY_PLAN(P1)
      CALL FFTW_DESTROY_PLAN(P2)

      NULLIFY(U, TU)

      CALL FFTW_FREE(P)

    END SUBROUTINE AUTOCORRELATION

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE TAPERING (FS, FE, X, POI, MUTE, TAPER)

      ! MUTE AND/OR TAPER THE RANDOM FIELD AROUND A POINT-OF-INTEREST. MUTING OCCURS WITHIN A RADIUS OF "MUTE" POINTS; TAPERING IS
      ! ACHIEVED BY APPLYING A HANNING WINDOW WITHIN A RADIUS IN THE RANGE "MUTE + 1" AND "MUTE + TAPER" POINTS.

      INTEGER(ISP), DIMENSION(3),                                     INTENT(IN)    :: FS, FE                !< START/END INDICES ALONG EACH DIRECTION
      REAL(RDP),    DIMENSION(FS(1):FE(1), FS(2):FE(2), FS(3):FE(3)), INTENT(INOUT) :: X                     !< RANDOM FIELD
      INTEGER(ISP), DIMENSION(3),                                     INTENT(IN)    :: POI                   !< POINT-OF-INTEREST (NODES)
      INTEGER(ISP),                                                   INTENT(IN)    :: MUTE, TAPER           !< MUTE/TAPER LENGTH (IN POINTS)
      INTEGER(ISP)                                                                  :: I, J, K               !< COUNTERS
      INTEGER(ISP)                                                                  :: D, DS, DM, DT         !< VARIOUS DISTANCES
      INTEGER(ISP)                                                                  :: DX, DY, DZ            !< NODE-POI DISTANCE
      REAL(RDP)                                                                     :: T                     !< TAPER ARGUMENT

      !-----------------------------------------------------------------------------------------------------------------------------

      DS = (MUTE + TAPER)**2
      DM = MUTE**2
      DT = TAPER**2

      DO K = FS(3), FE(3)

        DZ = (K - POI(3))**2                                                    !< DISTANCE ALONG Z

        DO J = FS(2), FE(2)

          DY = (J - POI(2))**2                                                  !< DISTANCE ALONG Y

          DO I = FS(1), FE(1)

            DX = (I - POI(1))**2                                                !< DISTANCE ALONG X

            D = DX + DY + DZ                                                    !< TOTAL DISTANCE

            IF (D .LE. DM) THEN

              X(I, J, K) = 0._RDP                                               !< MUTING

            ELSEIF ( (D .GT. DM) .AND. (D .LE. DS) ) THEN

              T = REAL(D - DM, RDP) / REAL(DT, RDP)

              X(I, J, K) = X(I, J, K) * (1._RDP - COS(T * PI / 2._RDP))         !< TAPERING

            ENDIF

          ENDDO
        ENDDO
      ENDDO


    END SUBROUTINE TAPERING

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *






END MODULE SCARFLIB
