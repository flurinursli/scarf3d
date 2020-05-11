MODULE SCARFLIB_SPEC

  USE, NON_INTRINSIC :: SCARFLIB_COMMON

  IMPLICIT NONE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  PRIVATE

  PUBLIC :: SCARF3D_SPEC

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTERFACE SCARF3D_SPEC
    MODULE PROCEDURE SCARF3D_UNSTRUCTURED_SPEC, SCARF3D_STRUCTURED_SPEC
  END INTERFACE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  ! INTERFACE TO C++ FUNCTIONS/SUBROUTINES
  INTERFACE

    SUBROUTINE SRNG(RAND) BIND(C, NAME="srng")
      USE, INTRINSIC     :: ISO_C_BINDING
      USE, NON_INTRINSIC :: SCARFLIB_COMMON
      IMPLICIT NONE
      REAL(C_FPP), DIMENSION(2), INTENT(OUT) :: RAND
    END SUBROUTINE SRNG

    SUBROUTINE NORMDIST(RAND) BIND(C, NAME="normdist")
      USE, INTRINSIC     :: ISO_C_BINDING
      USE, NON_INTRINSIC :: SCARFLIB_COMMON
      IMPLICIT NONE
      REAL(C_FPP), DIMENSION(2), INTENT(OUT) :: RAND
    END SUBROUTINE NORMDIST

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

  SUBROUTINE scarf3d_unstructured_spec(nc, fc, ds, x, y, z, dh, acf, cl, sigma, hurst, seed, poi, mute, taper, field, stats)

    ! Purpose:
    !   To compute a random field on structured/unstructured meshes based on the spectral representation method.
    !
    ! Revisions:
    !     Date                    Description of change
    !     ====                    =====================
    !   11/05/20                  Original version
    !

    REAL(FPP),                   DIMENSION(3),   INTENT(IN)  :: NC, FC
    REAL(FPP),                                   INTENT(IN)  :: DS
    REAL(FPP),                   DIMENSION(:),   INTENT(IN)  :: X, Y, Z              !< POSITION OF GRID POINTS ALONG X, Y, Z
    REAL(FPP),                                   INTENT(IN)  :: DH                   !< MIN DISTANCE BETWEEN GRID POINTS
    INTEGER(IPP),                                INTENT(IN)  :: ACF                  !< AUTOCORRELATION FUNCTION: "VK" OR "GAUSS"
    REAL(FPP),                   DIMENSION(3),   INTENT(IN)  :: CL                   !< CORRELATION LENGTH
    REAL(FPP),                                   INTENT(IN)  :: SIGMA                !< STANDARD DEVIATION
    REAL(FPP),                                   INTENT(IN)  :: HURST                !< HURST EXPONENT
    INTEGER(IPP),                                INTENT(IN)  :: SEED                 !< SEED NUMBER
    REAL(FPP),                   DIMENSION(:,:), INTENT(IN)  :: POI                  !< LOCATION OF POINT(S)-OF-INTEREST
    REAL(FPP),                                   INTENT(IN)  :: MUTE                 !< NUMBER OF POINTS WHERE MUTING IS APPLIED
    REAL(FPP),                                   INTENT(IN)  :: TAPER                !< NUMBER OF POINTS WHERE TAPERING IS APPLIED
    REAL(FPP),                   DIMENSION(:),   INTENT(OUT) :: FIELD                !< VELOCITY FIELD TO BE PERTURBED
    REAL(FPP),                   DIMENSION(6),   INTENT(OUT) :: STATS                 !< ERRORS AND TIMING FOR PERFORMANCE ANALYSIS
    INTEGER(IPP)                                             :: I, L
    INTEGER(IPP)                                             :: NPTS
    INTEGER(IPP)                                             :: IERR
    INTEGER(C_INT)                                           :: C_SEED
    INTEGER(IPP),  ALLOCATABLE,  DIMENSION(:)                :: NPOINTS
    REAL(FPP)                                                :: SCALING
    REAL(FPP)                                                :: KMAX, KMIN, KC
    REAL(FPP)                                                :: K, D                 !< USED TO COMPUTE THE COVARIANCE FUNCTION
    REAL(FPP)                                                :: PHI, THETA
    REAL(FPP)                                                :: V1, V2, V3
    REAL(FPP)                                                :: A, B, ARG
    REAL(REAL64)                                             :: TICTOC
    REAL(C_FPP),                 DIMENSION(2)                :: R
    REAL(FPP),                   DIMENSION(3)                :: SPAN
    REAL(FPP),                   DIMENSION(3)                :: MIN_EXTENT, MAX_EXTENT
    REAL(FPP),      ALLOCATABLE, DIMENSION(:)                :: MU, VAR

    !-------------------------------------------------------------------------------------------------------------------------------

    CALL mpi_comm_size(mpi_comm_world, world_size, ierr)
    CALL mpi_comm_rank(mpi_comm_world, world_rank, ierr)

    IF (world_rank == 0) THEN
      PRINT*, 'INPUT PARAMS @ UNSTRUCTURED-SPEC'
      PRINT*, 'NC ', NC
      PRINT*, 'FC ', FC
      PRINT*, 'DS ', DS
      PRINT*, 'X ', SIZE(X)
      PRINT*, 'Z ', SIZE(Z)
      PRINT*, 'DH ', DH
      PRINT*, 'ACF ', ACF
      PRINT*, 'CL ', CL
      PRINT*, 'SIGMA ', SIGMA
      PRINT*, 'HURST ', HURST
      PRINT*, 'SEED ', SEED
      PRINT*, 'POI ', SIZE(POI)
      PRINT*, 'MUTE ', MUTE
      PRINT*, 'TAPER ', TAPER
    ENDIF
    CALL mpi_barrier(mpi_comm_world, ierr)

    STATS = 0._FPP

    ! INITIALISE RANDOM GENERATOR
    c_seed = seed
    CALL set_seed(c_seed)

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    ! SET MIN/MAX RADIAL WAVENUMBER SUCH THAT "V1", "V2", "V3" BELOW ARE IN THE RANGE ["KMIN" "KMAX"] WHEN THE TRIGONOMETRIC TERMS
    ! EQUAL UNITY. THIS IS ACHIEVED BY TAKING THE LARGEST "KMIN*CL" AND THE SMALLEST "KMAX*CL". "CL" IS CONSIDERED IN THE CALCULATIONS
    ! BECAUSE THE RADIAL AVENUMBER "K" IS DIVIDED BY "CL" (SEE .EQ. 8 OF RAESS ET AL., 2019)
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    MIN_EXTENT = NC
    MAX_EXTENT = FC

    ! (GLOBAL) MODEL LIMITS
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, MIN_EXTENT, 3, REAL_TYPE, MPI_MIN, MPI_COMM_WORLD, IERR)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, MAX_EXTENT, 3, REAL_TYPE, MPI_MAX, MPI_COMM_WORLD, IERR)

    ! GLOBAL EFFECTIVE GRID SIZE
    SPAN = MAX_EXTENT - MIN_EXTENT

    KMIN = 2._fpp * PI * MAXVAL(CL / SPAN)

    KMAX = PI / DH * MINVAL(CL)

    ! CORNER WAVENUMBER FOR FILTERING SPECTRUM IS CONTROLLED BY MESH GRID-STEP
    KC = PI / DS * MINVAL(CL) !* SQRT(3._FPP)

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    ! CHECK IF THE FOLLOWING CONDITIONS FOR A CORRECT WAVENUMBER REPRESENTATION ARE VIOLATED:
    ! 1) KMIN < 1/CL
    ! 2) KMAX > 1/CL, OR DS <= CL / 2 (FRENJE & JUHLIN, 2000)
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    IF (ANY(2._FPP * PI / SPAN .GE. 1._FPP / CL)) STATS(1) = 1._FPP

    IF (ANY(DS .GT. CL / 2._FPP)) STATS(2) = 1._FPP

    !IF (WORLD_RANK == 0) PRINT*, 'KMAX ', KMIN, KMAX

    ! SET SCALING FACTOR
    SCALING = SIGMA / SQRT(REAL(NHARM, FPP))

    ! NUMBER OF POINTS WHERE RANDOM FIELD MUST BE CALCULATED
    NPTS = SIZE(X)

    ! INITIALISE RANDOM FIELD
    FIELD(:) = 0._FPP

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    ! COMPUTE RANDOM FIELD DIRECTLY AT EACH GRID POINT. USE UNIFORM DISTRIBUTION IN THE RANGE [0 1] AS MAJORANT FUNCTION.
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    ! LOOP OVER HARMONICS
    ! OPENACC: "FIELD" IS COPIED IN&OUT, ALL THE OTHERS ARE ONLY COPIED IN. "ARG" IS CREATED LOCALLY ON THE ACCELERATOR
    !$ACC DATA COPY(FIELD) COPYIN(X, Y, Z, SCALING, NPTS, A, B, V1, V2, V3)
    DO L = 1, NHARM

      !IF ((MOD(L, 100) == 0) .AND. (WORLD_RANK == 0)) PRINT*, WORLD_RANK, L

#ifdef TIMING
      CALL WATCH_START(TICTOC, MPI_COMM_SELF)
#endif

      DO

        ! FIRST SET OF RANDOM NUMBERS
        CALL SRNG(R)

        ! RANDOM VARIABLE "K`" MUST SPAN THE RANGE [KMIN, KMAX]
        K = R(1) * (KMAX - KMIN) + KMIN

        ! EVALUATE ORIGINAL PDF
        IF (ACF .EQ. 0) THEN
          D = K**2 / (1._FPP + K**2)**(1.5_FPP + HURST)
        ELSE
          D = K**2 * EXP(-0.25_FPP * K**2) / SQRT(PI)
        ENDIF

        ! TAKE "K" AND EXIT IF INEQUALITY R < D IS VERIFIED
        IF (R(2) .LT. D) EXIT

      ENDDO

      ! SECOND SET OF RANDOM NUMBERS
      CALL SRNG(R)

      ! NEGLECT CONTRIBUTION OF CURRENT HARMONIC IF WAVENUMBER "U" WAS LARGER THAN CORNER WAVENUMBER "KC"
      IF (K .GT. KC) CYCLE

      ! COMPUTE AZIMUTH AND POLAR ANGLES
      PHI   = R(1) * 2._FPP * PI
      THETA = ACOS(1._FPP - 2._FPP * R(2))

      V1 = K * SIN(PHI) * SIN(THETA) / CL(1)
      V2 = K * COS(PHI) * SIN(THETA) / CL(2)
      V3 = K * COS(THETA)            / CL(3)

      ! COMPUTE HARMONICS COEFFICIENT "A" AND "B"
      ! CALL NORMDIST(R)
      !
      ! A = R(1)
      ! B = R(2)

      ! BOX-MUELLER TRANSFORM
      CALL SRNG(R)
      A = SQRT(-2._FPP * LOG(R(1))) * COS(2._FPP * PI * R(2))

      CALL SRNG(R)
      B = SQRT(-2._FPP * LOG(R(1))) * COS(2._FPP * PI * R(2))


#ifdef TIMING
      CALL WATCH_STOP(TICTOC, MPI_COMM_SELF)

      STATS(5) = STATS(5) + TICTOC

      CALL WATCH_START(TICTOC, MPI_COMM_SELF)
#endif

      !$ACC WAIT

      ! UPDATE SELECTED VARIABLES IN GPU MEMORY
      !$ACC UPDATE DEVICE(A, B, V1, V2, V3)

      ! !$ACC PARALLEL LOOP PRIVATE(ARG) ASYNC
      ! DO I = 1, NPTS
      !   ARG = V1 * X(I) + V2 * Y(I) + V3 * Z(I)
      !   FIELD(I) = FIELD(I) + A * SIN(ARG) + B * COS(ARG)
      ! ENDDO
      ! !$ACC END PARALLEL LOOP

      !$ACC KERNELS ASYNC
      !$ACC LOOP INDEPENDENT PRIVATE(ARG)
      DO I = 1, NPTS
        ARG = V1 * X(I) + V2 * Y(I) + V3 * Z(I)
        FIELD(I) = FIELD(I) + A * SIN(ARG) + B * COS(ARG)
      ENDDO
      !$ACC END KERNELS


#ifdef TIMING
      CALL WATCH_STOP(TICTOC, MPI_COMM_SELF)

      STATS(6) = STATS(6) + TICTOC
#endif

    ENDDO
    ! END LOOP OVER HARMONICS

#ifdef TIMING
    CALL WATCH_START(TICTOC, MPI_COMM_SELF)
#endif

    !$ACC WAIT

    ! NORMALISE RANDOM FIELD
    ! !$ACC PARALLEL LOOP
    ! DO I = 1, NPTS
    !   FIELD(I) = FIELD(I) * SCALING
    ! ENDDO
    ! !$ACC END PARALLEL LOOP

    !$ACC KERNELS
    !$ACC LOOP INDEPENDENT
    DO I = 1, NPTS
      FIELD(I) = FIELD(I) * SCALING
    ENDDO
    !$ACC END KERNELS

    ! COPY "FIELD" BACK TO HOST AND FREE MEMORY ON DEVICE
    !$ACC END DATA

#ifdef TIMING
    CALL WATCH_STOP(TICTOC, MPI_COMM_SELF)

    STATS(6) = STATS(6) + TICTOC
#endif

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    ! COMPUTE VARIANCE AND MEAN OF RANDOM FIELD. IF DESIRED, REMOVE MEAN AND SET ACTUAL STANDARD DEVIATION TO DESIRED VALUE.
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    ALLOCATE(VAR(0:WORLD_SIZE - 1), MU(0:WORLD_SIZE - 1), NPOINTS(0:WORLD_SIZE - 1))

    ! COMPUTE VARIANCE AND MEAN OF RANDOM FIELD FOR EACH SINGLE PROCESS
    VAR(WORLD_RANK) = VARIANCE(FIELD)
    MU(WORLD_RANK)  = MEAN(FIELD)

    ! NUMBER OF GRID POINTS PER PROCESS
    NPOINTS(WORLD_RANK) = NPTS

    ! SHARE RESULTS
    CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, VAR, 1, REAL_TYPE, MPI_COMM_WORLD, IERR)
    CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, MU, 1, REAL_TYPE, MPI_COMM_WORLD, IERR)
    CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, NPOINTS, 1, MPI_INTEGER, MPI_COMM_WORLD, IERR)

    ! COMPUTE TOTAL VARIANCE ("STATS(3)") AND MEAN ("STATS(4)")
    CALL PARALLEL_VARIANCE(VAR, MU, NPOINTS, STATS(3), STATS(4))

    ! RETURN STANDARD DEVIATION
    STATS(3) = SQRT(STATS(3))

    ! DO WE NEED TO RESCALE THE RANDOM FIELD TO DESIRED (I.E. CONTINUOUS) STANDARD DEVIATION?
    ! IF (RESCALE .EQ. 1) THEN
    !
    !   SCALING = SIGMA / STATS(3)
    !
    !   DO I = 1, NPTS
    !     FIELD(I) = FIELD(I) * SCALING - STATS(4)
    !   ENDDO
    !
    !   STATS(3) = SIGMA
    !   STATS(4) = 0._FPP
    !
    ! ENDIF

    DEALLOCATE(VAR, MU, NPOINTS)

    ! APPLY TAPER/MUTE
    DO I = 1, SIZE(POI, 2)
      CALL TAPERING(X, Y, Z, FIELD, POI(:, I), MUTE, TAPER)
    ENDDO

    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

  END SUBROUTINE SCARF3D_UNSTRUCTURED_SPEC

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  !=================================================================================================================================
  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  SUBROUTINE SCARF3D_STRUCTURED_SPEC(NC, FC, DS, FS, FE, DH, ACF, CL, SIGMA, HURST, SEED, POI, MUTE, TAPER, FIELD, STATS)

    ! GLOBAL INDICES AND COMMUNICATOR SUBGROUPPING COULD BE HANDLED ELEGANTLY IF VIRTUAL TOPOLOGIES ARE USED IN CALLING PROGRAM.
    ! HOWEVER WE ASSUME THAT THESE ARE NOT USED AND THEREFORE WE ADOPT A SIMPLER APPROACH BASED ON COLLECTIVE CALLS.

    REAL(FPP),                   DIMENSION(3),     INTENT(IN)  :: NC, FC
    REAL(FPP),                                     INTENT(IN)  :: DS                   !< STRUCTURED GRID-STEP
    INTEGER(IPP),                DIMENSION(3),     INTENT(IN)  :: FS, FE               !< FIRST/LAST STRUCTURED GRID INDICES
    REAL(FPP),                                     INTENT(IN)  :: DH                   !< MIN GRID-STEP
    INTEGER(IPP),                                  INTENT(IN)  :: ACF                  !< AUTOCORRELATION FUNCTION: "VK" OR "GAUSS"
    REAL(FPP),                   DIMENSION(3),     INTENT(IN)  :: CL                   !< CORRELATION LENGTH
    REAL(FPP),                                     INTENT(IN)  :: SIGMA                !< STANDARD DEVIATION
    REAL(FPP),                                     INTENT(IN)  :: HURST                !< HURST EXPONENT
    INTEGER(IPP),                                  INTENT(IN)  :: SEED                 !< SEED NUMBER
    REAL(FPP),                   DIMENSION(:,:),   INTENT(IN)  :: POI                  !< LOCATION OF POINT(S)-OF-INTEREST
    REAL(FPP),                                     INTENT(IN)  :: MUTE                 !< NUMBER OF POINTS WHERE MUTING IS APPLIED
    REAL(FPP),                                     INTENT(IN)  :: TAPER                !< NUMBER OF POINTS WHERE TAPERING IS APPLIED
    REAL(FPP),                   DIMENSION(:,:,:), INTENT(OUT) :: FIELD                !< VELOCITY FIELD TO BE PERTURBED
    REAL(FPP),                   DIMENSION(6),     INTENT(OUT) :: STATS                !< ERRORS AND TIMING FOR PERFORMANCE ANALYSIS
    INTEGER(IPP)                                               :: I, J, K, L
    INTEGER(IPP)                                               :: IERR
    INTEGER(C_INT)                                             :: C_SEED
    INTEGER(IPP),                DIMENSION(3)                  :: NPTS
    INTEGER(IPP),  ALLOCATABLE,  DIMENSION(:)                  :: NPOINTS
    REAL(FPP)                                                  :: SCALING
    REAL(FPP)                                                  :: KMAX, KMIN, KC
    REAL(FPP)                                                  :: U, D                 !< USED TO COMPUTE THE COVARIANCE FUNCTION
    REAL(FPP)                                                  :: PHI, THETA
    REAL(FPP)                                                  :: V1, V2, V3
    REAL(FPP)                                                  :: A, B, ARG
    REAL(FPP)                                                  :: X, Y, Z
    REAL(REAL64)                                               :: TICTOC
    REAL(C_FPP),                 DIMENSION(2)                  :: R
    REAL(FPP),                   DIMENSION(3)                  :: SPAN
    REAL(FPP),                   DIMENSION(3)                  :: MIN_EXTENT, MAX_EXTENT
    REAL(FPP),      ALLOCATABLE, DIMENSION(:)                  :: MU, VAR

    !-------------------------------------------------------------------------------------------------------------------------------

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, WORLD_SIZE, IERR)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, WORLD_RANK, IERR)

    if (world_rank == 0) then
      print*, 'input params @ STRUCTURED-SPEC'
      print*, 'NC ', nc
      print*, 'FC ', fc
      print*, 'DS ', ds
      print*, 'FS ', fs
      print*, 'FE ', fe
      print*, 'DH ', dh
      print*, 'ACF ', acf
      print*, 'CL ', cl
      print*, 'SIGMA ', sigma
      print*, 'HURST ', hurst
      print*, 'SEED ', seed
      print*, 'POI ', size(poi)
      print*, 'MUTE ', mute
      print*, 'TAPER ', taper
    endif
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)


    STATS = 0._FPP

    ! INITIALISE RANDOM GENERATOR
    C_SEED = SEED
    CALL SET_SEED(C_SEED)

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    ! SET MIN/MAX RADIAL WAVENUMBER SUCH THAT "V1", "V2", "V3" BELOW ARE IN THE RANGE ["KMIN" "KMAX"] WHEN THE TRIGONOMETRIC TERMS
    ! EQUAL UNITY. THIS IS ACHIEVED BY TAKING THE LARGEST "KMIN*CL" AND THE SMALLEST "KMAX*CL". "CL" IS CONSIDERED IN THE CALCULATIONS
    ! BECAUSE THE RADIAL AVENUMBER "K" IS DIVIDED BY "CL" (SEE .EQ. 8 OF RAESS ET AL., 2019)
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    MIN_EXTENT = NC
    MAX_EXTENT = FC

    ! (GLOBAL) MODEL LIMITS
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, MIN_EXTENT, 3, REAL_TYPE, MPI_MIN, MPI_COMM_WORLD, IERR)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, MAX_EXTENT, 3, REAL_TYPE, MPI_MAX, MPI_COMM_WORLD, IERR)

    ! GLOBAL EFFECTIVE GRID SIZE
    SPAN = MAX_EXTENT - MIN_EXTENT

    KMIN = 2._FPP * PI * MAXVAL(CL / SPAN)

    KMAX = PI / DH * MINVAL(CL)

    ! CORNER WAVENUMBER FOR FILTERING SPECTRUM IS CONTROLLED BY MESH GRID-STEP
    KC = PI / DS * MINVAL(CL) !* SQRT(3._FPP)

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    ! CHECK IF THE FOLLOWING CONDITIONS FOR A CORRECT WAVENUMBER REPRESENTATION ARE VIOLATED:
    ! 1) KMIN < 1/CL
    ! 2) KMAX > 1/CL, OR DS <= CL / 2 (FRENJE & JUHLIN, 2000)
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    IF (ANY(2._FPP * PI / SPAN .GE. 1._FPP / CL)) STATS(1) = 1._FPP

    IF (ANY(DS .GT. CL / 2._FPP)) STATS(2) = 1._FPP

    !IF (WORLD_RANK == 0) PRINT*, 'KMAX ', KMIN, KMAX

    ! SET SCALING FACTOR
    SCALING = SIGMA / SQRT(REAL(NHARM, FPP))

    ! NUMBER OF POINTS WHERE RANDOM FIELD MUST BE CALCULATED
    NPTS(:) = FE(:) - FS(:) + 1

    ! INITIALISE RANDOM FIELD
    FIELD(:,:,:) = 0._FPP

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    ! COMPUTE RANDOM FIELD DIRECTLY AT EACH GRID POINT. USE UNIFORM DISTRIBUTION IN THE RANGE [0 1] AS MAJORANT FUNCTION.
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    ! LOOP OVER HARMONICS
    ! OPENACC: "FIELD" IS COPIED IN&OUT, ALL THE OTHERS ARE ONLY COPIED IN. "ARG" IS CREATED LOCALLY ON THE ACCELERATOR
    !$ACC DATA COPY(FIELD) COPYIN(DS, SCALING, NPTS, FS, A, B, V1, V2, V3)
    DO L = 1, NHARM

      !IF ((MOD(L, 100) == 0) .AND. (WORLD_RANK == 0)) PRINT*, WORLD_RANK, L

#ifdef TIMING
      CALL WATCH_START(TICTOC, MPI_COMM_SELF)
#endif

      DO

        ! FIRST SET OF RANDOM NUMBERS
        CALL SRNG(R)

        ! RANDOM VARIABLE "K`" MUST SPAN THE RANGE [KMIN, KMAX]
        U = R(1) * (KMAX - KMIN) + KMIN

        ! EVALUATE ORIGINAL PDF
        IF (ACF .EQ. 0) THEN
          D = U**2 / (1._FPP + U**2)**(1.5_FPP + HURST)
        ELSE
          D = U**2 * EXP(-0.25_FPP * U**2) / SQRT(PI)
        ENDIF

        ! TAKE "K" AND EXIT IF INEQUALITY R < D IS VERIFIED
        IF (R(2) .LT. D) EXIT

      ENDDO

      ! SECOND SET OF RANDOM NUMBERS
      CALL SRNG(R)

      ! NEGLECT CONTRIBUTION OF CURRENT HARMONIC IF WAVENUMBER "U" WAS LARGER THAN CORNER WAVENUMBER "KC"
      IF (U .GT. KC) CYCLE

      ! COMPUTE AZIMUTH AND POLAR ANGLES
      PHI   = R(1) * 2._FPP * PI
      THETA = ACOS(1._FPP - 2._FPP * R(2))

      V1 = U * SIN(PHI) * SIN(THETA) / CL(1)
      V2 = U * COS(PHI) * SIN(THETA) / CL(2)
      V3 = U * COS(THETA)            / CL(3)

      ! COMPUTE HARMONICS COEFFICIENT "A" AND "B"
      ! CALL NORMDIST(R)
      !
      ! A = R(1)
      ! B = R(2)

      ! BOX-MUELLER TRANSFORM
      CALL SRNG(R)
      A = SQRT(-2._FPP * LOG(R(1))) * COS(2._FPP * PI * R(2))

      CALL SRNG(R)
      B = SQRT(-2._FPP * LOG(R(1))) * COS(2._FPP * PI * R(2))


#ifdef TIMING
      CALL WATCH_STOP(TICTOC, MPI_COMM_SELF)

      STATS(5) = STATS(5) + TICTOC

      CALL WATCH_START(TICTOC, MPI_COMM_SELF)
#endif

      !$ACC WAIT

      ! UPDATE SELECTED VARIABLES IN GPU MEMORY
      !$ACC UPDATE DEVICE(A, B, V1, V2, V3)

      !$ACC PARALLEL LOOP PRIVATE(X, Y, Z, ARG) COLLAPSE(3) ASYNC
      DO K = 1, NPTS(3)
        DO J = 1, NPTS(2)
          DO I = 1, NPTS(1)
            X              = (I + FS(1) - 2) * DS
            Y              = (J + FS(2) - 2) * DS
            Z              = (K + FS(3) - 2) * DS
            ARG            = V1 * X + V2 * Y + V3 * Z
            FIELD(I, J, K) = FIELD(I, J, K) + A * SIN(ARG) + B * COS(ARG)
          ENDDO
        ENDDO
      ENDDO
      !$ACC END PARALLEL LOOP

      ! !$ACC KERNELS ASYNC
      ! !$ACC LOOP INDEPENDENT PRIVATE(X, Y, Z, ARG)
      ! DO K = 1, NPTS(3)
      !  Z = (K + FS(3) - 2) * DS
      !  DO J = 1, NPTS(2)
      !    Y = (J + FS(2) - 2) * DS
      !    DO I = 1, NPTS(1)
      !      X              = (I + FS(1) - 2) * DS
      !      ARG            = V1 * X + V2 * Y + V3 * Z
      !      FIELD(I, J, K) = FIELD(I, J, K) + A * SIN(ARG) + B * COS(ARG)
      !    ENDDO
      !  ENDDO
      ! ENDDO
      ! !$ACC END KERNELS


#ifdef TIMING
      CALL WATCH_STOP(TICTOC, MPI_COMM_SELF)

      STATS(6) = STATS(6) + TICTOC
#endif

    ENDDO
    ! END LOOP OVER HARMONICS

#ifdef TIMING
    CALL WATCH_START(TICTOC, MPI_COMM_SELF)
#endif

    !$ACC WAIT

    ! NORMALISE RANDOM FIELD
    !$ACC PARALLEL LOOP COLLAPSE(3)
    DO K = 1, NPTS(3)
      DO J = 1, NPTS(2)
        DO I = 1, NPTS(1)
          FIELD(I, J, K) = FIELD(I, J, K) * SCALING
        ENDDO
      ENDDO
    ENDDO
    !$ACC END PARALLEL LOOP

    ! !$ACC KERNELS
    ! !$ACC LOOP INDEPENDENT
    ! DO K = 1, NPTS(3)
    !  DO J = 1, NPTS(2)
    !    DO I = 1, NPTS(1)
    !      FIELD(I, J, K) = FIELD(I, J, K) * SCALING
    !    ENDDO
    !  ENDDO
    ! ENDDO
    ! !$ACC END KERNELS

    ! COPY "FIELD" BACK TO HOST AND FREE MEMORY ON DEVICE
    !$ACC END DATA

#ifdef TIMING
    CALL WATCH_STOP(TICTOC, MPI_COMM_SELF)

    STATS(6) = STATS(6) + TICTOC
#endif

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    ! COMPUTE VARIANCE AND MEAN OF RANDOM FIELD. IF DESIRED, REMOVE MEAN AND SET ACTUAL STANDARD DEVIATION TO DESIRED VALUE.
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    ALLOCATE(VAR(0:WORLD_SIZE - 1), MU(0:WORLD_SIZE - 1), NPOINTS(0:WORLD_SIZE - 1))

    ! COMPUTE VARIANCE AND MEAN OF RANDOM FIELD FOR EACH SINGLE PROCESS
    VAR(WORLD_RANK) = VARIANCE(FIELD)
    MU(WORLD_RANK)  = MEAN(FIELD)

    ! NUMBER OF GRID POINTS PER PROCESS
    NPOINTS(WORLD_RANK) = NPTS(1) * NPTS(2) * NPTS(3)

    ! SHARE RESULTS
    CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, VAR, 1, REAL_TYPE, MPI_COMM_WORLD, IERR)
    CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, MU, 1, REAL_TYPE, MPI_COMM_WORLD, IERR)
    CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, NPOINTS, 1, MPI_INTEGER, MPI_COMM_WORLD, IERR)

    ! COMPUTE TOTAL VARIANCE ("STATS(3)") AND MEAN ("STATS(4)")
    CALL PARALLEL_VARIANCE(VAR, MU, NPOINTS, STATS(3), STATS(4))

    ! RETURN STANDARD DEVIATION
    STATS(3) = SQRT(STATS(3))

    DEALLOCATE(VAR, MU, NPOINTS)

    ! APPLY TAPER/MUTE
    DO I = 1, SIZE(POI, 2)
      CALL TAPERING(DS, FS, FE, FIELD, POI(:, I), MUTE, TAPER)
    ENDDO

    CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

  END SUBROUTINE SCARF3D_STRUCTURED_SPEC

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  !=================================================================================================================================
  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --


END MODULE SCARFLIB_SPEC
