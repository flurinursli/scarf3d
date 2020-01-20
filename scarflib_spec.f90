MODULE SCARFLIB_SPECTRAL
  ! ALL VARIABLES AND SUBMODULE PROCEDURES ARE GLOBAL WITHIN THE SUBMODULE, BUT LIMITED TO IT.

!    USE, INTRINSIC     :: OMP_LIB

  USE, NON_INTRINSIC :: SCARFLIB_COMMON

  IMPLICIT NONE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  PRIVATE

  PUBLIC :: SCARF3D_SPEC

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! PROCEDURE POINTERS
  ! PROCEDURE(VK_PSDF), POINTER :: FUN

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! NUMBER OF HARMONICS IN SPECTRAL SUMMATION
  INTEGER(IPP), PARAMETER :: NHARM = 5000*4

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  !=================================================================================================================================
  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  !SUBROUTINE SCARF3D_SPEC(LS, LE, X, Y, Z, ACF, CL, SIGMA, HURST, SEED, POI, MUTE, TAPER, FIELD, TIME)
  SUBROUTINE SCARF3D_SPEC(X, Y, Z, DH, ACF, CL, SIGMA, HURST, SEED, POI, MUTE, TAPER, FIELD, INFO)

    ! GLOBAL INDICES AND COMMUNICATOR SUBGROUPPING COULD BE HANDLED ELEGANTLY IF VIRTUAL TOPOLOGIES ARE USED IN CALLING PROGRAM.
    ! HOWEVER WE ASSUME THAT THESE ARE NOT USED AND THEREFORE WE ADOPT A SIMPLER APPROACH BASED ON COLLECTIVE CALLS.

    REAL(FPP),                     DIMENSION(:),     INTENT(IN)  :: X, Y, Z              !< POSITION OF POINTS ALONG X, Y, Z
    REAL(FPP),                                       INTENT(IN)  :: DH                   !< GRID-STEP
    CHARACTER(LEN=*),                                INTENT(IN)  :: ACF                  !< AUTOCORRELATION FUNCTION: "VK" OR "GAUSS"
    REAL(FPP),                     DIMENSION(3),     INTENT(IN)  :: CL                   !< CORRELATION LENGTH
    REAL(FPP),                                       INTENT(IN)  :: SIGMA                !< STANDARD DEVIATION
    REAL(FPP),                                       INTENT(IN)  :: HURST                !< HURST EXPONENT
    INTEGER(IPP),                                    INTENT(IN)  :: SEED                 !< SEED NUMBER
    INTEGER(IPP),                  DIMENSION(:,:),   INTENT(IN)  :: POI                  !< LOCATION OF POINT(S)-OF-INTEREST
    INTEGER(IPP),                                    INTENT(IN)  :: MUTE                 !< NUMBER OF POINTS WHERE MUTING IS APPLIED
    INTEGER(IPP),                                    INTENT(IN)  :: TAPER                !< NUMBER OF POINTS WHERE TAPERING IS APPLIED
    REAL(FPP),                     DIMENSION(:),     INTENT(OUT) :: FIELD                !< VELOCITY FIELD TO BE PERTURBED
    REAL(FPP),                     DIMENSION(1),     INTENT(OUT) :: INFO                 !< ERRORS AND TIMING FOR PERFORMANCE ANALYSIS

    ! INTEGER(IPP),     DIMENSION(3),     INTENT(IN)    :: LS, LE         !< FIRST/LAST INDEX ALONG X, Y, Z
    ! REAL(FPP),        DIMENSION(:),     INTENT(IN)    :: X, Y, Z        !< POSITION OF POINTS ALONG X, Y, Z`
    ! CHARACTER(LEN=*),                   INTENT(IN)    :: ACF            !< AUTOCORRELATION FUNCTION: "VK" OR "GAUSS"
    ! REAL(FPP),        DIMENSION(3),     INTENT(IN)    :: CL             !< CORRELATION LENGTH
    ! REAL(FPP),                          INTENT(IN)    :: SIGMA          !< STANDARD DEVIATION
    ! REAL(FPP),                          INTENT(IN)    :: HURST          !< HURST EXPONENT
    ! INTEGER(IPP),                       INTENT(IN)    :: SEED           !< SEED NUMBER
    ! INTEGER(IPP),     DIMENSION(:,:),   INTENT(IN)    :: POI            !< LOCATION OF POINT(S)-OF-INTEREST
    ! INTEGER(IPP),                       INTENT(IN)    :: MUTE           !< NUMBER OF POINTS WHERE MUTING IS APPLIED
    ! INTEGER(IPP),                       INTENT(IN)    :: TAPER          !< NUMBER OF POINTS WHERE TAPERING IS APPLIED
    ! REAL(FPP),        DIMENSION(:),     INTENT(INOUT) :: FIELD          !< VELOCITY FIELD TO BE PERTURBED
    ! REAL(FPP),                          INTENT(OUT)   :: TIME           !< TIMING FOR PERFORMANCE ANALYSIS

    INTEGER(IPP)                                      :: I, L
    INTEGER(IPP)                                      :: NPTS
    INTEGER(IPP)                                      :: IERR
    REAL(FPP)                                         :: SCALING
    REAL(FPP)                                         :: KMAX, UMAX
    REAL(FPP)                                         :: K, D, U        !< USED TO COMPUTE THE COVARIANCE FUNCTION
    REAL(FPP)                                         :: PHI, THETA
    REAL(FPP)                                         :: V1, V2, V3
    REAL(FPP)                                         :: A, B, ARG
    REAL(FPP)                                         :: TICTOC
    REAL(FPP),        DIMENSION(2)                    :: R              !< RANDOM NUMBER

    !-------------------------------------------------------------------------------------------------------------------------------

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, WORLD_SIZE, IERR)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, WORLD_RANK, IERR)

    ! ALLOCATE(GS(3, 0:WORLD_SIZE-1), GE(3, 0:WORLD_SIZE-1))
    !
    ! ! STORE GLOBAL INDICES
    ! GS(:, WORLD_RANK) = LS
    ! GE(:, WORLD_RANK) = LE
    !
    ! ! MAKE ALL PROCESSES AWARE OF GLOBAL INDICES ALONG EACH AXIS
    ! CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, GS, 3, MPI_INTEGER, MPI_COMM_WORLD, IERR)
    ! CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, GE, 3, MPI_INTEGER, MPI_COMM_WORLD, IERR)
    !
    ! ! FIND TOTAL NUMBER OF POINTS ALONG EACH AXIS
    ! N = MAXVAL(GE, DIM = 2)

    ! INITIALISE RANDOM NUMBERS GENERATOR
    CALL SET_STREAM(SEED)

    ! SET NYQUIST WAVENUMBER BASED ON MINIMUM GRID-STEP
    KMAX = PI / DH * SQRT(CL(1)**2 + CL(2)**2 + CL(3)**2)

    ! SET SCALING FACTOR
    SCALING = SIGMA / SQRT(REAL(NHARM, FPP))

    ! NUMBER OF POINTS WHERE RANDOM FIELD MUST BE CALCULATED
    NPTS = SIZE(X)

    ! INITIALISE RANDOM FIELD
    FIELD(:) = 0._FPP

    ! SELECT PARAMETERS ACCORDING TO AUTOCORRELATION FUNCTION
    IF (ACF .EQ. 'VK') THEN

      IF (HURST .LT. 0.25_FPP) THEN
        UMAX = LOG(KMAX + 1._FPP)
      ELSEIF ( (HURST .GE. 0.25_FPP) .AND. (HURST .LT. 0.5_FPP) ) THEN
        UMAX = 2._FPP * (1._FPP - 1._FPP / SQRT(KMAX + 1._FPP))
      ELSEIF (HURST .GE. 0.5_FPP) THEN
        UMAX = ATAN(KMAX)
      ENDIF

    ELSEIF (ACF .EQ. 'GAUSS') THEN

      UMAX = SQRT(PI) * (ERF(KMAX * 0.5_FPP - 1._FPP) - ERF(-1._FPP))

    ENDIF

    ! START TIMER
    CALL WATCH_START(TICTOC)

    ! LOOP OVER HARMONICS
    ! OPENACC: "FIELD" IS COPIED IN&OUT, ALL THE OTHERS ARE ONLY COPIED IN. "ARG" IS CREATED LOCALLY ON THE ACCELERATOR
    !$ACC DATA COPY(FIELD) COPYIN(X, Y, Z, SCALING, NPTS, A, B, V1, V2, V3) CREATE(ARG)
    DO L = 1, NHARM

      DO

        ! FIRST SET OF RANDOM NUMBERS
        CALL RANDOM_NUMBER(R)

        ! RANDOM VARIABLE "U" MUST SPAN THE RANGE [0, UMAX]
        U = R(1) * UMAX

        ! OBTAIN "K" FROM MAJORANT CDF
        K = CDF2K(ACF, HURST, U)

        ! EVALUATE ORIGINAL PDF
        IF (ACF .EQ. 'VK') THEN
          D = K**2 / (1._FPP + K**2)**(1.5_FPP + HURST)
        ELSE
          D = K**2 * EXP(-0.25_FPP * K**2)
        ENDIF

        ! DIVIDE ORIGINAL PDF BY MAJORANT PDF
        D = D / PDF(ACF, HURST, K)

        ! TAKE "K" AND EXIT IF INEQUALITY R < F/G IS VERIFIED
        IF (R(2) .LT. D) EXIT

      ENDDO

      ! SECOND SET OF RANDOM NUMBERS
      CALL RANDOM_NUMBER(R)

      ! COMPUTE AZIMUTH AND POLAR ANGLES
      PHI   = R(1) * 2._FPP * PI
      THETA = ACOS(1._FPP - 2._FPP * R(2))

      V1 = K * SIN(PHI) * SIN(THETA) / CL(1)
      V2 = K * COS(PHI) * SIN(THETA) / CL(2)
      V3 = K * COS(THETA)            / CL(3)

      ! COMPUTE HARMONICS COEFFICIENT "A" AND "B"
      CALL RANDOM_NUMBER(R)
      A = BOX_MUELLER(R)

      CALL RANDOM_NUMBER(R)
      B = BOX_MUELLER(R)

      !$ACC WAIT

      ! UPDATE SELECTED VARIABLES IN GPU MEMORY
      !$ACC UPDATE DEVICE(A, B, V1, V2, V3)

      !$ACC KERNELS ASYNC
      !$ACC LOOP INDEPENDENT
      DO I = 1, NPTS
        ARG = V1 * X(I) + V2 * Y(I) + V3 * Z(I)
        FIELD(I) = FIELD(I) + A * SIN(ARG) + B * COS(ARG)
      ENDDO
      !$ACC END KERNELS

    ENDDO
    ! END LOOP OVER HARMONICS

    !$ACC WAIT

    ! NORMALISE RANDOM FIELD
    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT
    DO I = 1, NPTS
      FIELD(I) = FIELD(I) * SCALING
    ENDDO
    !$ACC END KERNELS

    ! COPY "FIELD" BACK TO HOST AND FREE MEMORY ON DEVICE
    !$ACC END DATA

    CALL WATCH_STOP(TICTOC)

    INFO = TICTOC

    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

  END SUBROUTINE SCARF3D_SPEC

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  !=================================================================================================================================
  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  REAL(FPP) FUNCTION BOX_MUELLER(RAND)

    ! GIVEN TWO UNIFORMLY DISTRIBUTED RANDOM NUMBERS, USE BOX-MUELLER ALGORITHM TO DRAW A NORMALLY DISTRIBUTED NUMBER IN THE RANGE
    ! (-INF INF)

    REAL(FPP), DIMENSION(2), INTENT(IN) :: RAND

    !-------------------------------------------------------------------------------------------------------------------------------

    BOX_MUELLER = SQRT(-2._FPP * LOG(RAND(1))) * COS(2._FPP * PI * RAND(2))

  END FUNCTION BOX_MUELLER

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  !=================================================================================================================================
  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  REAL(FPP) FUNCTION PDF(ACF, HURST, K)

    ! RETURN THE VALUE OF THE MAJORANT PDF AT WAVENUMBER "K". THE MAJORANT IS SELECTED BASED ON THE AUTOCORRELATION FUNCTION AND THE
    ! HURST EXPONENT (EXCEPT WHEN "ACF" IS GAUSSIAN).

    CHARACTER(LEN=*), INTENT(IN) :: ACF                        !< CORRELATION FUNCTION
    REAL(FPP),        INTENT(IN) :: HURST                      !< HURST EXPONENT
    REAL(FPP),        INTENT(IN) :: K                          !< WAVENUMBER

    !-------------------------------------------------------------------------------------------------------------------------------

    IF (ACF .EQ. 'VK') THEN

      IF (HURST .LT. 0.25_FPP) THEN
        PDF = 1.2_FPP / (1._FPP + K)
      ELSEIF ( (HURST .GE. 0.25_FPP) .AND. (HURST .LT. 0.5_FPP) ) THEN
        PDF = 1.3_FPP / (1._FPP + K)**1.5_FPP
      ELSE
        PDF = 1._FPP / (1._FPP + K**2)
      ENDIF

    ELSEIF (ACF .EQ. 'GAUSS') THEN

      PDF = 1.5_FPP * EXP(-0.25_FPP * (K - 2._FPP)**2)

    ENDIF

  END FUNCTION PDF

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  !=================================================================================================================================
  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  REAL(FPP) FUNCTION CDF2K(ACF, HURST, U)

    ! RETURN THE VALUE OF THE MAJORANT PDF AT WAVENUMBER "K". THE MAJORANT IS SELECTED BASED ON THE AUTOCORRELATION FUNCTION AND THE
    ! HURST EXPONENT (EXCEPT WHEN "ACF" IS GAUSSIAN).

    CHARACTER(LEN=*), INTENT(IN) :: ACF                        !< CORRELATION FUNCTION
    REAL(FPP),        INTENT(IN) :: HURST                      !< HURST EXPONENT
    REAL(FPP),        INTENT(IN) :: U                          !< RANDOM NUMBER IN THE RANGE [0, 1]
    REAL(FPP)                    :: DERFI

    !-------------------------------------------------------------------------------------------------------------------------------

    IF (ACF .EQ. 'VK') THEN

      IF (HURST .LT. 0.25_FPP) THEN
        CDF2K = EXP(U) - 1._FPP
      ELSEIF ( (HURST .GE. 0.25_FPP) .AND. (HURST .LT. 0.5_FPP) ) THEN
        CDF2K = 1._FPP / (1._FPP - U / 2._FPP)**2 - 1._FPP
      ELSE
        CDF2K = TAN(U)
      ENDIF

    ELSEIF (ACF .EQ. 'GAUSS') THEN

      CDF2K = 2._FPP * (1._FPP + DERFI(ERF(-1._FPP) + U / SQRT(PI)))

    ENDIF


  END FUNCTION CDF2K

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  !=================================================================================================================================
  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! FUNCTION CUMSUM(X) RESULT(V)
  !
  !   REAL(FPP),   DIMENSION(:),      INTENT(IN) :: X                        !< INPUT SEQUENCE
  !   REAL(FPP),   DIMENSION(SIZE(X))            :: V                        !< RESULT
  !   INTEGER(IPP)                               :: I, N                     !< COUNTERS
  !
  !   !-------------------------------------------------------------------------------------------------------------------------------
  !
  !   N = SIZE(X)
  !
  !   V(1) = X(1)
  !
  !   DO I = 2, N
  !     V(I) = V(I - 1) + X(I)
  !   ENDDO
  !
  ! END FUNCTION CUMSUM

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  !=================================================================================================================================
  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! FUNCTION VARIANCE(R) RESULT(V)
  !
  !   ! COMPUTE VARIANCE BASED ON THE COMPENSATED-SUMMATION VERSION OF THE TWO-PASS ALGORITHM
  !
  !   REAL(FPP),   DIMENSION(:,:,:), INTENT(IN) :: R
  !   INTEGER(IPP)                              :: I, J, K
  !   REAL(FPP)                                 :: MU, S1, S2, X, V
  !
  !   !-----------------------------------------------------------------------------------------------------------------------------
  !
  !   MU = MEAN(R)
  !
  !   S1 = 0._FPP
  !   S2 = 0._FPP
  !
  !   DO K = 1, SIZE(R, 3)
  !     DO J = 1, SIZE(R, 2)
  !       DO I = 1, SIZE(R, 1)
  !         X = R(I, J, K) - MU
  !         S1 = S1 + X
  !         S2 = S2 + X**2
  !       ENDDO
  !     ENDDO
  !   ENDDO
  !
  !   S1 = (S1**2) / REAL(SIZE(R), FPP)
  !
  !   V = (S2 - S1) / REAL(SIZE(R) - 1, FPP)
  !
  ! END FUNCTION VARIANCE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
  !===============================================================================================================================
  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
  ! 
  ! FUNCTION MEAN(R) RESULT(V)
  !
  !   REAL(FPP),   DIMENSION(:,:,:), INTENT(IN) :: R
  !   INTEGER(IPP)                              :: I, J, K
  !   REAL(FPP)                                 :: V, C
  !
  !   !-------------------------------------------------------------------------------------------------------------------------------
  !
  !   !V = COMPENSATED_SUM(R) / REAL(SIZE(R), FPP)
  !
  !   V = 0._FPP
  !   C = 1._FPP
  !
  !   DO K = 1, SIZE(R, 3)
  !     DO J = 1, SIZE(R, 2)
  !       DO I = 1, SIZE(R, 1)
  !         V = V + (R(I, J, K) - V) / C
  !         C = C + 1._FPP
  !       ENDDO
  !     ENDDO
  !   ENDDO
  !
  ! END FUNCTION MEAN

END MODULE SCARFLIB_SPECTRAL
