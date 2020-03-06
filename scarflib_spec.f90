MODULE SCARFLIB_SPEC
  ! ALL VARIABLES AND SUBMODULE PROCEDURES ARE GLOBAL WITHIN THE SUBMODULE, BUT LIMITED TO IT.

!    USE, INTRINSIC     :: OMP_LIB

  USE, NON_INTRINSIC :: SCARFLIB_COMMON

  IMPLICIT NONE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  PRIVATE

  PUBLIC :: SCARF3D_SPEC

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  ! INTERFACE TO C++ FUNCTIONS/SUBROUTINES
  INTERFACE

    SUBROUTINE SRNG(RAND) BIND(C, NAME="srng")
      USE, INTRINSIC     :: ISO_C_BINDING
      USE, NON_INTRINSIC :: SCARFLIB_COMMON
      IMPLICIT NONE
      REAL(C_FPP), DIMENSION(2), INTENT(OUT) :: RAND
    END SUBROUTINE SRNG

    SUBROUTINE SET_SEED(SEED) BIND(C, NAME="set_seed")
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      INTEGER(C_INT), VALUE :: SEED
    END SUBROUTINE SET_SEED

  END INTERFACE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! NUMBER OF HARMONICS IN SPECTRAL SUMMATION
  INTEGER(IPP), PARAMETER :: NHARM = 5000*4 / 2

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  !=================================================================================================================================
  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  !SUBROUTINE SCARF3D_SPEC(LS, LE, X, Y, Z, ACF, CL, SIGMA, HURST, SEED, POI, MUTE, TAPER, FIELD, TIME)
  SUBROUTINE SCARF3D_SPEC(X, Y, Z, DH, ACF, CL, SIGMA, HURST, SEED, POI, MUTE, TAPER, FIELD, INFO)

    ! GLOBAL INDICES AND COMMUNICATOR SUBGROUPPING COULD BE HANDLED ELEGANTLY IF VIRTUAL TOPOLOGIES ARE USED IN CALLING PROGRAM.
    ! HOWEVER WE ASSUME THAT THESE ARE NOT USED AND THEREFORE WE ADOPT A SIMPLER APPROACH BASED ON COLLECTIVE CALLS.

    REAL(FPP),     DIMENSION(:),     INTENT(IN)  :: X, Y, Z              !< POSITION OF POINTS ALONG X, Y, Z
    REAL(FPP),                       INTENT(IN)  :: DH                   !< GRID-STEP
    INTEGER(IPP),                    INTENT(IN)  :: ACF                  !< AUTOCORRELATION FUNCTION: "VK" OR "GAUSS"
    REAL(FPP),     DIMENSION(3),     INTENT(IN)  :: CL                   !< CORRELATION LENGTH
    REAL(FPP),                       INTENT(IN)  :: SIGMA                !< STANDARD DEVIATION
    REAL(FPP),                       INTENT(IN)  :: HURST                !< HURST EXPONENT
    INTEGER(IPP),                    INTENT(IN)  :: SEED                 !< SEED NUMBER
    REAL(FPP),     DIMENSION(:,:),   INTENT(IN)  :: POI                  !< LOCATION OF POINT(S)-OF-INTEREST
    REAL(FPP),                       INTENT(IN)  :: MUTE                 !< NUMBER OF POINTS WHERE MUTING IS APPLIED
    REAL(FPP),                       INTENT(IN)  :: TAPER                !< NUMBER OF POINTS WHERE TAPERING IS APPLIED
    REAL(FPP),     DIMENSION(:),     INTENT(OUT) :: FIELD                !< VELOCITY FIELD TO BE PERTURBED
    REAL(FPP),     DIMENSION(2),     INTENT(OUT) :: INFO                 !< ERRORS AND TIMING FOR PERFORMANCE ANALYSIS
    INTEGER(IPP)                                 :: I, L
    INTEGER(IPP)                                 :: NPTS
    INTEGER(IPP)                                 :: IERR
    INTEGER(C_INT)                               :: C_SEED
    REAL(FPP)                                    :: SCALING
    REAL(FPP)                                    :: KMAX, UMAX
    REAL(FPP)                                    :: K, D, U        !< USED TO COMPUTE THE COVARIANCE FUNCTION
    REAL(FPP)                                    :: PHI, THETA
    REAL(FPP)                                    :: V1, V2, V3
    REAL(FPP)                                    :: A, B, ARG
    REAL(REAL64)                                 :: TICTOC
    REAL(C_FPP),    DIMENSION(2)                 :: R

    !-------------------------------------------------------------------------------------------------------------------------------

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, WORLD_SIZE, IERR)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, WORLD_RANK, IERR)

    INFO = 0._FPP

    ! INITIALISE RANDOM NUMBERS GENERATOR
    !CALL SET_STREAM(SEED)

    C_SEED = SEED
    CALL SET_SEED(C_SEED)

    ! SET NYQUIST WAVENUMBER BASED ON MINIMUM GRID-STEP
    KMAX = PI / DH * SQRT(CL(1)**2 + CL(2)**2 + CL(3)**2)


    IF (WORLD_RANK == 0) PRINT*, 'KMAX ', KMAX

    ! SET SCALING FACTOR
    SCALING = SIGMA / SQRT(REAL(NHARM, FPP))

    ! NUMBER OF POINTS WHERE RANDOM FIELD MUST BE CALCULATED
    NPTS = SIZE(X)

    ! INITIALISE RANDOM FIELD
    FIELD(:) = 0._FPP

    ! SELECT PARAMETERS ACCORDING TO AUTOCORRELATION FUNCTION
    IF (ACF .EQ. 0) THEN

      IF (HURST .LT. 0.25_FPP) THEN
        UMAX = LOG(KMAX + 1._FPP)
        !UMAX = KMAX
      ELSEIF ( (HURST .GE. 0.25_FPP) .AND. (HURST .LT. 0.5_FPP) ) THEN
        UMAX = 2._FPP * (1._FPP - 1._FPP / SQRT(KMAX + 1._FPP))
      ELSEIF (HURST .GE. 0.5_FPP) THEN
        UMAX = ATAN(KMAX)
      ENDIF

    ELSEIF (ACF .EQ. 1) THEN

      UMAX = SQRT(PI) * (ERF(KMAX * 0.5_FPP - 1._FPP) - ERF(-1._FPP))

    ENDIF

    UMAX = KMAX

    ! LOOP OVER HARMONICS
    ! OPENACC: "FIELD" IS COPIED IN&OUT, ALL THE OTHERS ARE ONLY COPIED IN. "ARG" IS CREATED LOCALLY ON THE ACCELERATOR
    !$ACC DATA COPY(FIELD) COPYIN(X, Y, Z, SCALING, NPTS, A, B, V1, V2, V3) CREATE(ARG)
    DO L = 1, NHARM

      IF ((MOD(L, 100) == 0) .AND. (WORLD_RANK == 0)) PRINT*, WORLD_RANK, L

      CALL WATCH_START(TICTOC, MPI_COMM_SELF)

      DO

        ! FIRST SET OF RANDOM NUMBERS
        !CALL RANDOM_NUMBER(R)
        CALL SRNG(R)

        ! RANDOM VARIABLE "U" MUST SPAN THE RANGE [0, UMAX]
        U = R(1) * UMAX

        ! OBTAIN "K" FROM MAJORANT CDF
        K = CDF2K(ACF, HURST, U)

        ! EVALUATE ORIGINAL PDF
        IF (ACF .EQ. 0) THEN
          D = K**2 / (1._FPP + K**2)**(1.5_FPP + HURST)
        ELSE
          D = K**2 * EXP(-0.25_FPP * K**2) / SQRT(PI)
        ENDIF

        ! DIVIDE ORIGINAL PDF BY MAJORANT PDF
        D = D / PDF(ACF, HURST, K)

        ! TAKE "K" AND EXIT IF INEQUALITY R < F/G IS VERIFIED
        IF (R(2) .LT. D) EXIT

      ENDDO

      ! SECOND SET OF RANDOM NUMBERS
      !CALL RANDOM_NUMBER(R)
      CALL SRNG(R)

      ! COMPUTE AZIMUTH AND POLAR ANGLES
      PHI   = R(1) * 2._FPP * PI
      THETA = ACOS(1._FPP - 2._FPP * R(2))

      V1 = K * SIN(PHI) * SIN(THETA) / CL(1)
      V2 = K * COS(PHI) * SIN(THETA) / CL(2)
      V3 = K * COS(THETA)            / CL(3)

      ! COMPUTE HARMONICS COEFFICIENT "A" AND "B"
      !CALL RANDOM_NUMBER(R)
      CALL SRNG(R)
      A = BOX_MUELLER(R)

      !CALL RANDOM_NUMBER(R)
      CALL SRNG(R)
      B = BOX_MUELLER(R)

      CALL WATCH_STOP(TICTOC, MPI_COMM_SELF)

      INFO(1) = INFO(1) + TICTOC

      CALL WATCH_START(TICTOC, MPI_COMM_SELF)

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

      CALL WATCH_STOP(TICTOC, MPI_COMM_SELF)

      INFO(2) = INFO(2) + TICTOC

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

    INTEGER(IPP), INTENT(IN) :: ACF                        !< CORRELATION FUNCTION
    REAL(FPP),    INTENT(IN) :: HURST                      !< HURST EXPONENT
    REAL(FPP),    INTENT(IN) :: K                          !< WAVENUMBER

    !-------------------------------------------------------------------------------------------------------------------------------

    ! IF (ACF .EQ. 0) THEN
    !
    !   IF (HURST .LT. 0.25_FPP) THEN
    !     PDF = 1.2_FPP / (1._FPP + K)
    !   ELSEIF ( (HURST .GE. 0.25_FPP) .AND. (HURST .LT. 0.5_FPP) ) THEN
    !     PDF = 1.3_FPP / (1._FPP + K)**1.5_FPP
    !   ELSE
    !     PDF = 1._FPP / (1._FPP + K**2)
    !   ENDIF
    !
    ! ELSEIF (ACF .EQ. 1) THEN
    !
    !   PDF = 1.5_FPP * EXP(-0.25_FPP * (K - 2._FPP)**2)
    !
    ! ENDIF

    PDF = 1._FPP

  END FUNCTION PDF

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  !=================================================================================================================================
  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  REAL(FPP) FUNCTION CDF2K(ACF, HURST, U)

    ! RETURN THE VALUE OF THE MAJORANT PDF AT WAVENUMBER "K". THE MAJORANT IS SELECTED BASED ON THE AUTOCORRELATION FUNCTION AND THE
    ! HURST EXPONENT (EXCEPT WHEN "ACF" IS GAUSSIAN).

    INTEGER(IPP), INTENT(IN) :: ACF                        !< CORRELATION FUNCTION
    REAL(FPP),    INTENT(IN) :: HURST                      !< HURST EXPONENT
    REAL(FPP),    INTENT(IN) :: U                          !< RANDOM NUMBER IN THE RANGE [0, 1]
    REAL(FPP)                :: DERFI

    !-------------------------------------------------------------------------------------------------------------------------------

    ! IF (ACF .EQ. 0) THEN
    !
    !   IF (HURST .LT. 0.25_FPP) THEN
    !     CDF2K = EXP(U) - 1._FPP
    !     !CDF2K = U
    !   ELSEIF ( (HURST .GE. 0.25_FPP) .AND. (HURST .LT. 0.5_FPP) ) THEN
    !     CDF2K = 1._FPP / (1._FPP - U / 2._FPP)**2 - 1._FPP
    !   ELSE
    !     CDF2K = TAN(U)
    !   ENDIF
    !
    ! ELSEIF (ACF .EQ. 1) THEN
    !
    !   CDF2K = 2._FPP * (1._FPP + DERFI(ERF(-1._FPP) + U / SQRT(PI)))
    !
    ! ENDIF

    CDF2K = U

  END FUNCTION CDF2K

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  !=================================================================================================================================
  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --


END MODULE SCARFLIB_SPEC
