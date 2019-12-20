MODULE SCARFLIB_FFT2

  ! ALL VARIABLES AND SUBMODULE PROCEDURES ARE GLOBAL WITHIN THE SUBMODULE, BUT LIMITED TO IT.

!    USE, INTRINSIC     :: OMP_LIB

  USE, NON_INTRINSIC :: SCARFLIB

  IMPLICIT NONE

  INCLUDE 'fftw3.f03'

  PRIVATE

  PUBLIC :: SCARF3D_STRUCTURED

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  ! PROCEDURE POINTERS
  PROCEDURE(VK_PSDF), POINTER :: FUN

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  ! POINTER FOR FOURIER TRANSFORM
  COMPLEX(C_CPP),              DIMENSION(:),            POINTER :: CDUM => NULL()

  ! POINTER FOR FOURIER TRANSFORM
  REAL(C_FPP),                 DIMENSION(:),            POINTER :: RDUM => NULL()

  ! POINTER TO FFTW PLAN
  TYPE(C_PTR)                                                   :: PC

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE SCARF3D_STRUCTURED(LS, LE, DH, ACF, CL, SIGMA, HURST, SEED, POI, MUTE, TAPER, FIELD, TIME)

      ! GLOBAL INDICES AND COMMUNICATOR SUBGROUPPING COULD BE HANDLED ELEGANTLY IF VIRTUAL TOPOLOGIES ARE USED IN CALLING PROGRAM.
      ! HOWEVER WE ASSUME THAT THESE ARE NOT USED AND THEREFORE WE ADOPT A SIMPLER APPROACH BASED ON COLLECTIVE CALLS.

      INTEGER(IPP),                  DIMENSION(3),     INTENT(IN)    :: LS, LE         !< FIRST/LAST INDEX ALONG X, Y, Z
      REAL(FPP),                                       INTENT(IN)    :: DH             !< GRID-STEP
      CHARACTER(LEN=*),                                INTENT(IN)    :: ACF            !< AUTOCORRELATION FUNCTION: "VK" OR "GAUSS"
      REAL(FPP),                     DIMENSION(3),     INTENT(IN)    :: CL             !< CORRELATION LENGTH
      REAL(FPP),                                       INTENT(IN)    :: SIGMA          !< STANDARD DEVIATION
      REAL(FPP),                                       INTENT(IN)    :: HURST          !< HURST EXPONENT
      INTEGER(IPP),                                    INTENT(IN)    :: SEED           !< SEED NUMBER
      INTEGER(IPP),                  DIMENSION(:,:),   INTENT(IN)    :: POI            !< LOCATION OF POINT(S)-OF-INTEREST
      INTEGER(IPP),                                    INTENT(IN)    :: MUTE           !< NUMBER OF POINTS WHERE MUTING IS APPLIED
      INTEGER(IPP),                                    INTENT(IN)    :: TAPER          !< NUMBER OF POINTS WHERE TAPERING IS APPLIED
      REAL(FPP),                     DIMENSION(:,:,:), INTENT(INOUT) :: FIELD          !< VELOCITY FIELD TO BE PERTURBED
      REAL(FPP),                     DIMENSION(3),     INTENT(OUT)   :: TIME           !< TIMING FOR PERFORMANCE ANALYSIS
      COMPLEX(FPP),     ALLOCATABLE, DIMENSION(:,:,:)                :: SPEC           !< SPECTRUM/RANDOM FIELD
      INTEGER(IPP)                                                   :: IERR
      INTEGER(IPP)                                                   :: I, J, K        !< COUNTERS
      REAL(FPP)                                                      :: SCALING        !< SCALING FACTOR
      REAL(FPP)                                                      :: TICTOC
      REAL(FPP),                     DIMENSION(2)                    :: ET             !< DUMMY FOR ELAPSED TIME

      !-----------------------------------------------------------------------------------------------------------------------------

      ! GET RANK NUMBER
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, WORLD_RANK, IERR)

      ! GET AVAILABLE MPI PROCESSES
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, WORLD_SIZE, IERR)

      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

      ALLOCATE(GS(3, 0:WORLD_SIZE-1), GE(3, 0:WORLD_SIZE-1))

      ! STORE GLOBAL INDICES
      GS(:, WORLD_RANK) = LS
      GE(:, WORLD_RANK) = LE

      ! MAKE ALL PROCESSES AWARE OF GLOBAL INDICES ALONG EACH AXIS
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, GS, 3, MPI_INTEGER, MPI_COMM_WORLD, IERR)
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, GE, 3, MPI_INTEGER, MPI_COMM_WORLD, IERR)

      ! FIND TOTAL NUMBER OF POINTS ALONG EACH AXIS
      N = MAXVAL(GE, DIM = 2)

      ! ALLOCATE MEMORY FOR SPECTRUM (NOTE USE OF GLOBAL INDICES)
      !ALLOCATE(SPEC(LS(1):LE(1), LS(2):LE(2), LS(3):LE(3)))
      ALLOCATE(SPEC(LE(1)-LS(1)+1, LE(2)-LS(2)+1, LE(3)-LS(3)+1))

      ! COMPUTE SPECTRUM AND APPLY HERMITIAN SYMMETRY
      CALL COMPUTE_SPECTRUM(LS, LE, DH, ACF, CL, SIGMA, HURST, SEED, SPEC, ET)

      TIME(1:2) = ET

      ! START TIMER
      CALL WATCH_START(TICTOC)

      ! TRANSFORM ALONG EACH DIRECTION
      CALL TRANSFORM_ALONG_Z(SPEC)
      CALL TRANSFORM_ALONG_Y(SPEC)
      CALL TRANSFORM_ALONG_X(SPEC)

      CALL WATCH_STOP(TICTOC)

      TIME(3) = TICTOC

      ! SCALING PARAMETER
      SCALING = 1._FPP / SQRT(REAL(N(1), FPP) * REAL(N(2), FPP) * REAL(N(3), FPP) * DH**3)

      DO K = 1, SIZE(SPEC, 3)
        DO J = 1, SIZE(SPEC, 2)
          DO I = 1, SIZE(SPEC, 1)
            SPEC(I, J, K)  = SPEC(I, J, K)  * SCALING                                  !< APPLY SCALING FACTOR
            FIELD(I, J, K) = FIELD(I, J, K) * (1._FPP + REAL(SPEC(I, J, K), FPP))      !< ADD PERTURBATIONS
          ENDDO
        ENDDO
      ENDDO


      !CALL IO_WRITE_ONE(REAL(SPEC, FPP), 'random_struct.bin')
      !CALL IO_WRITE_SLICE(1, 50, REAL(SPEC, FPP), 'random_struct_slice.bin')

      DEALLOCATE(SPEC)

      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

    END SUBROUTINE SCARF3D_STRUCTURED

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE TRANSFORM_ALONG_Z(SPECTRUM)

      ! COMPUTE MANY IN-PLACE COMPLEX-TO-COMPLEX IFFT ALONG Z-DIRECTION ON INPUT SPECTRUM. PROCESSES ARE GATHERED IN Z-ORIENTED PENCILS
      ! BASED ON THE INPUT GEOMETRY. INSIDE EACH PENCIL, THE ALGORITHM WORKS ON Y-SLICES (TO REQUIRE AS LESS MEMORY AS POSSIBLE) TO
      ! EXCHANGE DATA AND COMPUTE THE IFFT.

      COMPLEX(FPP),              DIMENSION(:,:,:), INTENT(INOUT) :: SPECTRUM                   !< INPUT FIELD
      COMPLEX(FPP), ALLOCATABLE, DIMENSION(:)                    :: SENDBUF, RECVBUF           !< SENDER/RECEIVER BUFFER
      INTEGER(IPP)                                               :: I, J, K, P, C              !< COUNTERS
      INTEGER(IPP)                                               :: RANK, NTASKS, IERR         !< MPI STUFF
      INTEGER(IPP)                                               :: PENCIL
      INTEGER(IPP)                                               :: NX, NY, NZ                 !< TOTAL NUMBER OF POINTS FOR PENCIL ALONG EACH AXIS
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: COLOR                      !< ARRAY USED TO CREATE PENCIL COMMUNICATOR
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: I0, I1                     !< LOCAL I-INDEX FOR PROCESSES BELONGING TO PENCIL
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: K0, K1                     !< LOCAL K-INDEX FOR PROCESSES BELONGING TO PENCIL
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: LX, LZ                     !< LOCAL NUMBER OF POINTS ALONG X AND Z AXIS
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: SENDCOUNTS, RECVCOUNTS     !< MPI STUFF
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: SDISPLS, RDISPLS           !< MPI STUFF
      LOGICAL,                   DIMENSION(2)                    :: BOOL
      TYPE(C_PTR)                                                :: P1

      !-------------------------------------------------------------------------------------------------------------------------------

      ALLOCATE(COLOR(0:WORLD_SIZE-1))

      ! GROUP PROCESSES SPANNING SAME X/Y-AXIS IN STENCILS ORIENTED ALONG Z-AXIS
      DO I = 0, WORLD_SIZE - 1
        COLOR(I) = 0
        BOOL(1) = (GS(1, I) .EQ. GS(1, WORLD_RANK)) .AND. (GE(1, I) .EQ. GE(1, WORLD_RANK))
        BOOL(2) = (GS(2, I) .EQ. GS(2, WORLD_RANK)) .AND. (GE(2, I) .EQ. GE(2, WORLD_RANK))
        IF (ALL(BOOL .EQV. .TRUE.)) COLOR(I) = I + 1
      ENDDO

      COLOR(WORLD_RANK) = MAXVAL(COLOR, DIM = 1)

      ! CREATE COMMUNICATOR SUBGROUP
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, COLOR(WORLD_RANK), WORLD_RANK, PENCIL, IERR)

      DEALLOCATE(COLOR)

      ! PROCESS ID AND COMMUNICATOR SIZE
      CALL MPI_COMM_RANK(PENCIL, RANK, IERR)
      CALL MPI_COMM_SIZE(PENCIL, NTASKS, IERR)

      ! THESE ARRAYS WILL CONTAIN THE NUMBER OF POINTS LOCAL TO EACH PROCESS IN CURRENT PENCIL
      ALLOCATE(LX(0:NTASKS-1), LZ(0:NTASKS-1))

      ! NUMBER OF POINTS ALONG X AND Y. THESE ARE INDENTICAL FOR ALL PROCESSES IN CURRENT PENCIL.
      NX = SIZE(SPECTRUM, 1)
      NY = SIZE(SPECTRUM, 2)

      ! POINTS ALONG Z-AXIS VARY ACCORDING TO RANK
      LZ(RANK) = SIZE(SPECTRUM, 3)

      ! MAKE ALL PROCESSES IN THE COMMUNICATOR AWARE OF HOW MANY POINTS ALONG Z-AXIS EACH PROCESS HAS
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, LZ, 1, MPI_INTEGER, PENCIL, IERR)

      ! TOTAL NUMBER OF POINTS ALONG Z-AXIS
      NZ = SUM(LZ)

      ALLOCATE(K0(0:NTASKS-1), K1(0:NTASKS-1))

      ! FIRST/LAST K-INDEX FOR EACH PROCESS
      DO P = 0, NTASKS - 1
        IF (P .EQ. 0) THEN
          K0(P) = 1
          K1(P) = LZ(P)
        ELSE
          K0(P) = K1(P - 1) + 1
          K1(P) = K1(P - 1) + LZ(P)
        ENDIF
      ENDDO

      ALLOCATE(I0(0:NTASKS-1), I1(0:NTASKS-1))

      ! DISTRIBUTE POINTS ALONG X-AXIS BETWEEN PROCESSES
      CALL MPI_SPLIT_TASK(NX, NTASKS, I0, I1)

      ! NUMBER OF POINTS ALONG X-AXIS FOR EACH PROCESS
      DO P = 0, NTASKS - 1
        LX(P) = I1(P) - I0(P) + 1
      ENDDO

      ALLOCATE(SENDCOUNTS(0:NTASKS-1), SDISPLS(0:NTASKS-1), RECVCOUNTS(0:NTASKS-1), RDISPLS(0:NTASKS-1))

      ! SET COUNTS/DISPLACEMENT ARRAYS
      DO P = 0, NTASKS - 1

        SENDCOUNTS(P) = LZ(RANK) * LX(P)              !< POINTS SEND FROM PROCESS "RANK" TO PROCESS "P"
        RECVCOUNTS(P) = LZ(P) * LX(RANK)              !< POINTS RECEIVED FROM PROCESS "RANK" BY PROCESS "P"

        IF (P .EQ. 0) THEN
          SDISPLS(P) = 0
          RDISPLS(P) = 0
        ELSE
          SDISPLS(P) = SDISPLS(P - 1) + SENDCOUNTS(P - 1)
          RDISPLS(P) = RDISPLS(P - 1) + RECVCOUNTS(P - 1)
        ENDIF

      ENDDO

      ALLOCATE(SENDBUF(NX * LZ(RANK)))
      ALLOCATE(RECVBUF(NZ * LX(RANK)))

      ! BUFFER "CDUM" WILL BE USED FOR TRANSFORM AND IT IS ALLOCATED USING FFTW MEMORY ALLOCATION UTILITY
      PC = FFTW_ALLOC_COMPLEX(INT(LX(RANK) * NZ, C_SIZE_T))

      CALL C_F_POINTER(PC, CDUM, [LX(RANK) * NZ])

      ! PREPARE FFTW PLAN ALONG Z-AXIS ("LX(RANK)" TRANSFORMS, EACH HAVING "NZ" NUMBER OF POINTS)
      P1 = FFTW_PLAN_MANY_DFT(1, [NZ], LX(RANK), CDUM, [NZ], LX(RANK), 1, CDUM, [NZ], LX(RANK), 1, FFTW_BACKWARD, FFTW_ESTIMATE)

      ! WORK Y-SLICES ONE BY ONE
      DO J = 1, NY

        C = 0

        ! REARRANGE ELEMENTS OF "SPECTRUM"
        DO P = 0, NTASKS - 1
          DO K = 1, LZ(RANK)
            DO I = I0(P), I1(P)
              C          = C + 1
              SENDBUF(C) = SPECTRUM(I, J, K)
            ENDDO
          ENDDO
        ENDDO

        ! EXCHANGE DATA
        CALL MPI_ALLTOALLV(SENDBUF, SENDCOUNTS, SDISPLS, COMPLEX_TYPE, CDUM, RECVCOUNTS, RDISPLS, COMPLEX_TYPE, PENCIL, IERR)

        ! IFFT
        CALL FFTW_EXECUTE_DFT(P1, CDUM, CDUM)

        C = 0

        ! REARRANGE ELEMENTS OF "CDUM" INTO "RECVBUF"
        DO P = 0, NTASKS - 1
          DO I = 1, LX(RANK)
            DO K = K0(P), K1(P)
              C          = C + 1
              RECVBUF(C) = CDUM(I + (K - 1)*LX(RANK))
            ENDDO
          ENDDO
        ENDDO

        ! SEND TRANSFORMED DATA BACK
        CALL MPI_ALLTOALLV(RECVBUF, RECVCOUNTS, RDISPLS, COMPLEX_TYPE, SENDBUF, SENDCOUNTS, SDISPLS, COMPLEX_TYPE, PENCIL, IERR)

        ! REARRANGE DATA INTO ORIGINAL LAYOUT
        DO K = 1, LZ(RANK)
          DO I = 1, NX
            SPECTRUM(I, J, K) = SENDBUF(K + (I - 1)*LZ(RANK))
          ENDDO
        ENDDO

      ENDDO

      ! DESTROY PLAN
      CALL FFTW_DESTROY_PLAN(P1)

      ! RELEASE MEMORY
      NULLIFY(CDUM)
      CALL FFTW_FREE(PC)

      DEALLOCATE(SENDBUF, RECVBUF, SDISPLS, RDISPLS, SENDCOUNTS, RECVCOUNTS)
      DEALLOCATE(I0, I1, K0, K1)

      ! RELEASE COMMUNICATOR
      CALL MPI_COMM_FREE(PENCIL, IERR)

    END SUBROUTINE TRANSFORM_ALONG_Z

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE TRANSFORM_ALONG_Y(SPECTRUM)

      ! COMPUTE MANY IN-PLACE COMPLEX-TO-COMPLEX IFFT ALONG Y-DIRECTION ON INPUT SPECTRUM. PROCESSES ARE GATHERED IN Y-ORIENTED PENCILS
      ! BASED ON THE INPUT GEOMETRY. INSIDE EACH PENCIL, THE ALGORITHM WORKS ON Z-SLICES (TO REQUIRE AS LESS MEMORY AS POSSIBLE) TO
      ! EXCHANGE DATA AND COMPUTE THE IFFT.

      COMPLEX(FPP),              DIMENSION(:,:,:), INTENT(INOUT) :: SPECTRUM                   !< INPUT FIELD
      COMPLEX(FPP), ALLOCATABLE, DIMENSION(:)                 :: SENDBUF, RECVBUF           !< SENDER/RECEIVER BUFFER
      INTEGER(IPP)                                            :: I, J, K, P, C              !< COUNTERS
      INTEGER(IPP)                                            :: RANK, NTASKS, IERR         !< MPI STUFF
      INTEGER(IPP)                                            :: PENCIL
      INTEGER(IPP)                                            :: NX, NY, NZ                 !< TOTAL NUMBER OF POINTS FOR PENCIL ALONG EACH AXIS
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                 :: COLOR                      !< ARRAY USED TO CREATE PENCIL COMMUNICATOR
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                 :: I0, I1                     !< LOCAL I-INDEX FOR PROCESSES BELONGING TO PENCIL
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                 :: J0, J1                     !< LOCAL J-INDEX FOR PROCESSES BELONGING TO PENCIL
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                 :: LX, LY                     !< LOCAL NUMBER OF POINTS ALONG X AND Y AXIS
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                 :: SENDCOUNTS, RECVCOUNTS     !< MPI STUFF
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                 :: SDISPLS, RDISPLS           !< MPI STUFF
      LOGICAL,                   DIMENSION(2)                 :: BOOL
      TYPE(C_PTR)                                             :: P1

      !-------------------------------------------------------------------------------------------------------------------------------

      ALLOCATE(COLOR(0:WORLD_SIZE-1))

      ! GROUP PROCESSES SPANNING SAME X/Z-AXIS IN STENCILS ORIENTED ALONG Z-AXIS
      DO I = 0, WORLD_SIZE - 1
        COLOR(I) = 0
        BOOL(1) = (GS(1, I) .EQ. GS(1, WORLD_RANK)) .AND. (GE(1, I) .EQ. GE(1, WORLD_RANK))
        BOOL(2) = (GS(3, I) .EQ. GS(3, WORLD_RANK)) .AND. (GE(3, I) .EQ. GE(3, WORLD_RANK))
        IF (ALL(BOOL .EQV. .TRUE.)) COLOR(I) = I + 1
      ENDDO

      COLOR(WORLD_RANK) = MAXVAL(COLOR, DIM = 1)

      ! CREATE COMMUNICATOR SUBGROUP
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, COLOR(WORLD_RANK), WORLD_RANK, PENCIL, IERR)

      DEALLOCATE(COLOR)

      ! PROCESS ID AND COMMUNICATOR SIZE
      CALL MPI_COMM_RANK(PENCIL, RANK, IERR)
      CALL MPI_COMM_SIZE(PENCIL, NTASKS, IERR)

      ! THESE ARRAYS WILL CONTAIN THE NUMBER OF POINTS LOCAL TO EACH PROCESS IN CURRENT PENCIL
      ALLOCATE(LX(0:NTASKS-1), LY(0:NTASKS-1))

      ! NUMBER OF POINTS ALONG X AND Z. THESE ARE INDENTICAL FOR ALL PROCESSES IN CURRENT PENCIL.
      NX = SIZE(SPECTRUM, 1)
      NZ = SIZE(SPECTRUM, 3)

      ! POINTS ALONG Y-AXIS VARY ACCORDING TO RANK
      LY(RANK) = SIZE(SPECTRUM, 2)

      ! MAKE ALL PROCESSES IN THE COMMUNICATOR AWARE OF HOW MANY POINTS ALONG Y-AXIS EACH PROCESS HAS
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, LY, 1, MPI_INTEGER, PENCIL, IERR)

      ! TOTAL NUMBER OF POINTS ALONG Y-AXIS
      NY = SUM(LY)

      ALLOCATE(J0(0:NTASKS-1), J1(0:NTASKS-1))

      ! FIRST/LAST K-INDEX FOR EACH PROCESS
      DO P = 0, NTASKS - 1
        IF (P .EQ. 0) THEN
          J0(P) = 1
          J1(P) = LY(P)
        ELSE
          J0(P) = J1(P - 1) + 1
          J1(P) = J1(P - 1) + LY(P)
        ENDIF
      ENDDO

      ALLOCATE(I0(0:NTASKS-1), I1(0:NTASKS-1))

      ! DISTRIBUTE POINTS ALONG X-AXIS BETWEEN PROCESSES
      CALL MPI_SPLIT_TASK(NX, NTASKS, I0, I1)

      ! NUMBER OF POINTS ALONG X-AXIS FOR EACH PROCESS
      DO P = 0, NTASKS - 1
        LX(P) = I1(P) - I0(P) + 1
      ENDDO

      ALLOCATE(SENDCOUNTS(0:NTASKS-1), SDISPLS(0:NTASKS-1), RECVCOUNTS(0:NTASKS-1), RDISPLS(0:NTASKS-1))

      ! SET COUNTS/DISPLACEMENT ARRAYS
      DO P = 0, NTASKS - 1

        SENDCOUNTS(P) = LY(RANK) * LX(P)              !< POINTS SEND FROM PROCESS "RANK" TO PROCESS "P"
        RECVCOUNTS(P) = LY(P) * LX(RANK)              !< POINTS RECEIVED FROM PROCESS "RANK" BY PROCESS "P"

        IF (P .EQ. 0) THEN
          SDISPLS(P) = 0
          RDISPLS(P) = 0
        ELSE
          SDISPLS(P) = SDISPLS(P - 1) + SENDCOUNTS(P - 1)
          RDISPLS(P) = RDISPLS(P - 1) + RECVCOUNTS(P - 1)
        ENDIF

      ENDDO

      ALLOCATE(SENDBUF(NX * LY(RANK)))
      ALLOCATE(RECVBUF(NY * LX(RANK)))

      ! BUFFER "CDUM" WILL BE USED FOR TRANSFORM AND IT IS ALLOCATED USING FFTW MEMORY ALLOCATION UTILITY
      PC = FFTW_ALLOC_COMPLEX(INT(LX(RANK) * NY, C_SIZE_T))

      CALL C_F_POINTER(PC, CDUM, [LX(RANK) * NY])

      ! PREPARE FFTW PLAN ALONG Y-AXIS ("LX(RANK)" TRANSFORMS, EACH HAVING "NY" NUMBER OF POINTS)
      P1 = FFTW_PLAN_MANY_DFT(1, [NY], LX(RANK), CDUM, [NY], LX(RANK), 1, CDUM, [NY], LX(RANK), 1, FFTW_BACKWARD, FFTW_ESTIMATE)

      ! WORK Z-SLICES ONE BY ONE
      DO K = 1, NZ

        C = 0

        ! REARRANGE ELEMENTS OF "SPECTRUM"
        DO P = 0, NTASKS - 1
          DO J = 1, LY(RANK)
            DO I = I0(P), I1(P)
              C          = C + 1
              SENDBUF(C) = SPECTRUM(I, J, K)
            ENDDO
          ENDDO
        ENDDO

        ! EXCHANGE DATA
        CALL MPI_ALLTOALLV(SENDBUF, SENDCOUNTS, SDISPLS, COMPLEX_TYPE, CDUM, RECVCOUNTS, RDISPLS, COMPLEX_TYPE, PENCIL, IERR)

        ! IFFT
        CALL FFTW_EXECUTE_DFT(P1, CDUM, CDUM)

        C = 0

        ! REARRANGE ELEMENTS OF "CDUM" INTO "RECVBUF"
        DO P = 0, NTASKS - 1
          DO I = 1, LX(RANK)
            DO J = J0(P), J1(P)
              C          = C + 1
              RECVBUF(C) = CDUM(I + (J - 1)*LX(RANK))
            ENDDO
          ENDDO
        ENDDO

        ! SEND TRANSFORMED DATA BACK
        CALL MPI_ALLTOALLV(RECVBUF, RECVCOUNTS, RDISPLS, COMPLEX_TYPE, SENDBUF, SENDCOUNTS, SDISPLS, COMPLEX_TYPE, PENCIL, IERR)

        ! REARRANGE DATA INTO ORIGINAL LAYOUT
        DO J = 1, LY(RANK)
          DO I = 1, NX
            SPECTRUM(I, J, K) = SENDBUF(J + (I - 1)*LY(RANK))
          ENDDO
        ENDDO

      ENDDO

      ! DESTROY PLAN
      CALL FFTW_DESTROY_PLAN(P1)

      ! RELEASE MEMORY
      NULLIFY(CDUM)
      CALL FFTW_FREE(PC)

      DEALLOCATE(SENDBUF, RECVBUF, SDISPLS, RDISPLS, SENDCOUNTS, RECVCOUNTS)
      DEALLOCATE(I0, I1, J0, J1)

      ! RELEASE COMMUNICATOR
      CALL MPI_COMM_FREE(PENCIL, IERR)

    END SUBROUTINE TRANSFORM_ALONG_Y

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE TRANSFORM_ALONG_X(SPECTRUM)

      ! COMPUTE MANY IN-PLACE COMPLEX-TO-REAL IFFT ALONG X-DIRECTION ON INPUT SPECTRUM. PROCESSES ARE GATHERED IN X-ORIENTED PENCILS
      ! BASED ON THE INPUT GEOMETRY. INSIDE EACH PENCIL, THE ALGORITHM WORKS ON Z-SLICES (TO REQUIRE AS LESS MEMORY AS POSSIBLE) TO
      ! EXCHANGE DATA AND COMPUTE THE IFFT.

      COMPLEX(FPP),              DIMENSION(:,:,:), INTENT(INOUT) :: SPECTRUM                       !< INPUT FIELD
      COMPLEX(FPP), ALLOCATABLE, DIMENSION(:)                    :: SENDBUF                        !< SENDER BUFFER
      INTEGER(IPP)                                               :: I, J, K, P, C                  !< COUNTERS
      INTEGER(IPP)                                               :: RANK, NTASKS, IERR             !< MPI STUFF
      INTEGER(IPP)                                               :: NYQUIST, PENCIL
      INTEGER(IPP)                                               :: IMAX, P_IMAX
      INTEGER(IPP)                                               :: NX, NY, NZ                     !< TOTAL NUMBER OF POINTS FOR PENCIL ALONG EACH AXIS
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: OFFSET
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: COLOR                          !< ARRAY USED TO CREATE PENCIL COMMUNICATOR
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: I0, I1                         !< LOCAL I-INDEX FOR PROCESSES BELONGING TO PENCIL
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: J0, J1                         !< LOCAL J-INDEX FOR PROCESSES BELONGING TO PENCIL
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: LX, LY                         !< LOCAL NUMBER OF POINTS ALONG X AND Y AXIS
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: SENDCOUNTS, RECVCOUNTS         !< MPI STUFF
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: SDISPLS, RDISPLS               !< MPI STUFF
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: R_SENDCOUNTS, R_RECVCOUNTS     !< MPI STUFF
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                    :: R_SDISPLS, R_RDISPLS           !< MPI STUFF
      LOGICAL,                   DIMENSION(2)                    :: BOOL
      REAL(FPP),    ALLOCATABLE, DIMENSION(:)                    :: R_SENDBUF, RECVBUF
      TYPE(C_PTR)                                                :: P1

      !-------------------------------------------------------------------------------------------------------------------------------

      ! INDEX OF NYQUIST FREQUENCY ALONG X-AXIS
      NYQUIST = N(1) / 2 + 1

      ! PARAMETER USED TO DETERMINE ELEMENTS TO BE SENT/RECEIVED
      !OFFSET = GS(1, WORLD_RANK) - 1

      ALLOCATE(COLOR(0:WORLD_SIZE-1))

      ! GROUP PROCESSES SPANNING SAME Y/Z-AXIS IN STENCILS ORIENTED ALONG X-AXIS
      DO I = 0, WORLD_SIZE - 1
        COLOR(I) = 0
        BOOL(1) = (GS(2, I) .EQ. GS(2, WORLD_RANK)) .AND. (GE(2, I) .EQ. GE(2, WORLD_RANK))
        BOOL(2) = (GS(3, I) .EQ. GS(3, WORLD_RANK)) .AND. (GE(3, I) .EQ. GE(3, WORLD_RANK))
        IF (ALL(BOOL .EQV. .TRUE.)) COLOR(I) = I + 1
      ENDDO

      COLOR(WORLD_RANK) = MAXVAL(COLOR, DIM = 1)

      ! CREATE COMMUNICATOR SUBGROUP
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, COLOR(WORLD_RANK), WORLD_RANK, PENCIL, IERR)

      DEALLOCATE(COLOR)

      ! PROCESS ID AND COMMUNICATOR SIZE
      CALL MPI_COMM_RANK(PENCIL, RANK, IERR)
      CALL MPI_COMM_SIZE(PENCIL, NTASKS, IERR)

      ALLOCATE(OFFSET(0:NTASKS-1))

      OFFSET(RANK) = GS(1, WORLD_RANK) - 1

      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, OFFSET, 1, MPI_INTEGER, PENCIL, IERR)

      ! THESE ARRAYS WILL CONTAIN THE NUMBER OF POINTS LOCAL TO EACH PROCESS IN CURRENT PENCIL
      ALLOCATE(LX(0:NTASKS-1), LY(0:NTASKS-1))

      ! NUMBER OF POINTS ALONG Y AND Z. THESE ARE INDENTICAL FOR ALL PROCESSES IN CURRENT PENCIL.
      NY = SIZE(SPECTRUM, 2)
      NZ = SIZE(SPECTRUM, 3)

      ! POINTS ALONG X-AXIS VARY ACCORDING TO RANK
      LX(RANK) = SIZE(SPECTRUM, 1)

      ! MAKE ALL PROCESSES IN THE COMMUNICATOR AWARE OF HOW MANY POINTS ALONG X-AXIS EACH PROCESS HAS
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, LX, 1, MPI_INTEGER, PENCIL, IERR)

      ! TOTAL NUMBER OF POINTS ALONG X-AXIS
      NX = SUM(LX)

      ! MAX I-INDEX USED TO SEND DATA
      IMAX = MIN(LX(RANK) + OFFSET(RANK), NYQUIST) - OFFSET(RANK)

      IF (IMAX .LT. 1) IMAX = 0

      ALLOCATE(I0(0:NTASKS-1), I1(0:NTASKS-1))

      ! FIRST/LAST I-INDEX FOR EACH PROCESS
      DO P = 0, NTASKS - 1
        IF (P .EQ. 0) THEN
          I0(P) = 1
          I1(P) = LX(P)
        ELSE
          I0(P) = I1(P - 1) + 1
          I1(P) = I1(P - 1) + LX(P)
        ENDIF
      ENDDO

      ALLOCATE(J0(0:NTASKS-1), J1(0:NTASKS-1))

      ! DISTRIBUTE POINTS ALONG Y-AXIS BETWEEN PROCESSES
      CALL MPI_SPLIT_TASK(NY, NTASKS, J0, J1)

      ! NUMBER OF POINTS ALONG Y-AXIS FOR EACH PROCESS
      DO P = 0, NTASKS - 1
        LY(P) = J1(P) - J0(P) + 1
      ENDDO

      ALLOCATE(SENDCOUNTS(0:NTASKS-1), SDISPLS(0:NTASKS-1), RECVCOUNTS(0:NTASKS-1), RDISPLS(0:NTASKS-1))
      ALLOCATE(R_SENDCOUNTS(0:NTASKS-1), R_SDISPLS(0:NTASKS-1), R_RECVCOUNTS(0:NTASKS-1), R_RDISPLS(0:NTASKS-1))

      ! SET COUNTS/DISPLACEMENT ARRAYS
      DO P = 0, NTASKS - 1

        ! DEFINE "IMAX" FOR PROCESS "P"
        P_IMAX = MIN(LX(P) + OFFSET(P), NYQUIST) - OFFSET(P)

        IF (P_IMAX .LT. 1) P_IMAX = 0

        ! PROCESS COMPLETELY TO THE RIGHT OF "NYQUIST" WON'T SEND ANY ELEMENT ("SENDCOUNTS=0")
        SENDCOUNTS(P) = LY(P) * IMAX                                            !< POINTS SEND FROM PROCESS "RANK" TO PROCESS "P"

        ! A PROCESS WON'T RECEIVE ANYTHING FROM PROCESS "P" IF THE LATTER IS TO THE RIGHT OF "NYQUIST" ("RECVCOUNTS=0")
        RECVCOUNTS(P) = LY(RANK) * P_IMAX                                       !< POINTS RECEIVED FROM PROCESS "RANK" BY PROCESS "P"

        R_SENDCOUNTS(P) = LY(RANK) * LX(P)                                      !< POINTS SEND FROM PROCESS "RANK" TO PROCESS "P"
        R_RECVCOUNTS(P) = LY(P) * LX(RANK)                                      !< POINTS RECEIVED FROM PROCESS "RANK" BY PROCESS "P"

        IF (P .EQ. 0) THEN
          SDISPLS(P)   = 0
          RDISPLS(P)   = 0
          R_SDISPLS(P) = 0
          R_RDISPLS(P) = 0
        ELSE
          SDISPLS(P)   = SDISPLS(P - 1) + SENDCOUNTS(P - 1)
          RDISPLS(P)   = RDISPLS(P - 1) + RECVCOUNTS(P - 1)
          R_SDISPLS(P) = R_SDISPLS(P - 1) + R_SENDCOUNTS(P - 1)
          R_RDISPLS(P) = R_RDISPLS(P - 1) + R_RECVCOUNTS(P - 1)
        ENDIF

      ENDDO

      ALLOCATE(SENDBUF(NY * IMAX), R_SENDBUF(NX * LY(RANK)))
      ALLOCATE(RECVBUF(NY * LX(RANK)))

      ! ALLOCATE MEMORY USING FFTW ALLOCATION UTILITY
      PC = FFTW_ALLOC_COMPLEX(INT(LY(RANK) * NYQUIST, C_SIZE_T))

      ! FOR IN-PLACE TRANSFORMS, OUTPUT REAL ARRAY IS SLIGHTLY LONGER THAN ACTUAL PHYSICAL DIMENSION
      CALL C_F_POINTER(PC, CDUM, [LY(RANK) * NYQUIST])
      CALL C_F_POINTER(PC, RDUM, [LY(RANK) * NYQUIST * 2])

      ! PREPARE FFTW PLAN ALONG X-AXIS ("LY(RANK)" TRANSFORMS, EACH HAVING "N(1)" NUMBER OF POINTS)
      P1 = FFTW_PLAN_MANY_DFT_C2R(1, [N(1)], LY(RANK), CDUM, [N(1)], LY(RANK), 1, RDUM, [N(1)], LY(RANK), 1, FFTW_ESTIMATE)

      ! WORK Z-SLICES ONE BY ONE
      DO K = 1, NZ

        C = 0

        ! REARRANGE ELEMENTS OF "SPECTRUM"
        DO P = 0, NTASKS - 1
          DO I = 1, IMAX
            DO J = J0(P), J1(P)
              C          = C + 1
              SENDBUF(C) = SPECTRUM(I, J, K)
            ENDDO
          ENDDO
        ENDDO

        ! EXCHANGE DATA
        CALL MPI_ALLTOALLV(SENDBUF, SENDCOUNTS, SDISPLS, COMPLEX_TYPE, CDUM, RECVCOUNTS, RDISPLS, COMPLEX_TYPE, PENCIL, IERR)

        ! IFFT
        CALL FFTW_EXECUTE_DFT_C2R(P1, CDUM, RDUM)

        C = 0

        ! REARRANGE ELEMENTS OF "RDUM" INTO "RECVBUF"
        DO P = 0, NTASKS - 1
          DO J = 1, LY(RANK)
            DO I = I0(P), I1(P)
              C            = C + 1
              R_SENDBUF(C) = RDUM(J + (I - 1)*LY(RANK))
            ENDDO
          ENDDO
        ENDDO

        ! SEND TRANSFORMED DATA BACK
        CALL MPI_ALLTOALLV(R_SENDBUF, R_SENDCOUNTS, R_SDISPLS, REAL_TYPE, RECVBUF, R_RECVCOUNTS, R_RDISPLS, REAL_TYPE, PENCIL, IERR)

        ! REARRANGE DATA INTO ORIGINAL LAYOUT
        DO J = 1, NY
          DO I = 1, LX(RANK)
            SPECTRUM(I, J, K) = CMPLX(RECVBUF(I + (J - 1)*LX(RANK)), 0._FPP, FPP)
          ENDDO
        ENDDO

      ENDDO

      ! DESTROY PLAN
      CALL FFTW_DESTROY_PLAN(P1)

      ! RELEASE MEMORY
      NULLIFY(CDUM, RDUM)
      CALL FFTW_FREE(PC)

      DEALLOCATE(SENDBUF, RECVBUF, SDISPLS, RDISPLS, SENDCOUNTS, RECVCOUNTS)
      DEALLOCATE(R_SENDBUF, R_SDISPLS, R_RDISPLS, R_SENDCOUNTS, R_RECVCOUNTS)
      DEALLOCATE(I0, I1, J0, J1, OFFSET)

      ! RELEASE COMMUNICATOR
      CALL MPI_COMM_FREE(PENCIL, IERR)

    END SUBROUTINE TRANSFORM_ALONG_X

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE COMPUTE_SPECTRUM(LS, LE, DH, ACF, CL, SIGMA, HURST, SEED, SPEC, TIME)

      ! COMPUTE THE SPECTRUM OF A RANDOM FIELD CHARACTERISED BY A SPECIFIC AUTOCORRELATION FUNCTION, CORRELATION LENGTH, STANDARD DEVIATION,
      ! AND HURST EXPONENT (NOT USED FOR GAUSSIAN FIELDS). BASED ON THE FOURIER INTEGRAL METHOD OF PARDO-IGUZQUIZA AND CHICA-OLMO.

      INTEGER(IPP),     DIMENSION(3),                                   INTENT(IN)  :: LS, LE                !< GLOBAL INDICES
      REAL(FPP),                                                        INTENT(IN)  :: DH                    !< GRID-STEP
      CHARACTER(LEN=*),                                                 INTENT(IN)  :: ACF                   !< AUTOCORRELATION FUNCTION
      REAL(FPP),        DIMENSION(3),                                   INTENT(IN)  :: CL                    !< CORRELATION LENGTH (CAN BE ANISOTROPIC)
      REAL(FPP),                                                        INTENT(IN)  :: SIGMA                 !< STANDARD DEVIATION (SIGMA%/100)
      REAL(FPP),                                                        INTENT(IN)  :: HURST                 !< HURST EXPONENT
      INTEGER(IPP),                                                     INTENT(IN)  :: SEED                  !< INITIAL SEED NUMBER
      COMPLEX(FPP),     DIMENSION(LS(1):LE(1),LS(2):LE(2),LS(3):LE(3)), INTENT(OUT) :: SPEC                  !< SPECTRUM
      REAL(FPP),        DIMENSION(2),                                   INTENT(OUT) :: TIME                  !< ELAPSED TIME
      INTEGER(IPP)                                                                  :: I, J, K               !< INDICES
      LOGICAL,          DIMENSION(3)                                                :: BOOL
      REAL(FPP)                                                                     :: BUTTER, NUM, AMP      !< USED TO COMPUTE SPECTRUM
      REAL(FPP)                                                                     :: KNYQ, KC, KR          !< USED TO COMPUTE SPECTRUM
      REAL(FPP)                                                                     :: TICTOC                !< USED FOR TIMING
      REAL(FPP),        DIMENSION(3)                                                :: DK                    !< RESOLUTION IN WAVENUMBER DOMAIN
      REAL(FPP),        DIMENSION(N(1))                                             :: HARVEST               !< RANDOM VALUES
      REAL(FPP),        DIMENSION(N(1))                                             :: KX                    !< WAVENUMBER VECTOR (X)
      REAL(FPP),        DIMENSION(N(2))                                             :: KY                    !< WAVENUMBER VECTOR (Y)
      REAL(FPP),        DIMENSION(N(3))                                             :: KZ                    !< WAVENUMBER VECTOR (Z)

      !-----------------------------------------------------------------------------------------------------------------------------

      ! RESOLUTION IN WAVENUMBER DOMAIN ALONG EACH DIRECTION
      DO I = 1, 3
        DK(I) = 2._FPP * PI / (REAL(N(I), FPP) * DH)
      ENDDO

      ! NYQUIST WAVENUMBER
      KNYQ = PI / DH

      ! CORNER WAVENUMBER FOR SPECTRUM TAPERING TO AVOID ALIASING IS HALF NYQUIST WAVENUMBER
      KC = KNYQ / 2._FPP

      ! VECTORS GO FROM 0 TO NYQUIST AND THEN BACK AGAIN UNTIL DK
      KX = [[(I * DK(1), I = 0, N(1)/2)], [(I * DK(1), I = N(1)/2-1, 1, -1)]]
      KY = [[(J * DK(2), J = 0, N(2)/2)], [(J * DK(2), J = N(2)/2-1, 1, -1)]]
      KZ = [[(K * DK(3), K = 0, N(3)/2)], [(K * DK(3), K = N(3)/2-1, 1, -1)]]

      ! COMPUTE PART OF POWER SPECTRAL DENSITY OUTSIDE LOOP
      IF (ACF .EQ. 'VK') THEN
        NUM = 8._FPP * SQRT(PI**3) * GAMMA(HURST + 1.5_FPP) * SIGMA**2 * PRODUCT(CL) / GAMMA(HURST)
        FUN => VK_PSDF
      ELSEIF (ACF .EQ. 'GAUSS') THEN
        NUM = SIGMA**2 * PRODUCT(CL) * SQRT(PI**3)
        FUN => GS_PSDF
      ENDIF

      ! INITIALISE RANDOM GENERATOR
      CALL SET_STREAM(SEED)

      ! START TIMER
      CALL WATCH_START(TICTOC)

      ! COMPUTE SPECTRUM
      DO K = 1, N(3)

        BOOL(3) = (K .GE. LS(3)) .AND. (K .LE. LE(3))

        DO J = 1, N(2)

          BOOL(2) = (J .GE. LS(2)) .AND. (J .LE. LE(2))

          CALL RANDOM_NUMBER(HARVEST)

          DO I = 1, N(1)

            BOOL(1) = (I .GE. LS(1)) .AND. (I .LE. LE(1))

            IF (ALL(BOOL .EQV. .TRUE.)) THEN

              ! RADIAL WAVENUMBER
              !KR = SQRT(KX(I)**2 + KY(J)**2 + KZ(K)**2)

              ! FOURTH-ORDER LOW-PASS BUTTERWORTH FILTER
              !BUTTER = 1._FPP / SQRT(1._FPP + (KR / KC)**(2 * 4))

              ! NOW "KR" IS THE PRODUCT "K * CL"
              KR = (KX(I) * CL(1))**2 + (KY(J) * CL(2))**2 + (KZ(K) * CL(3))**2

              ! COMPLETE POWER SPECTRAL DENSITY AND GO TO AMPLITUDE SPECTRUM
              AMP = SQRT(NUM / FUN(KR, HURST))

              ! APPLY FILTER
              !AMP = AMP * BUTTER

              HARVEST(I) = HARVEST(I) * 2._FPP * PI

              ! COMBINE AMPLITUDE AND RANDOM PHASE
              SPEC(I, J, K) = CMPLX(COS(HARVEST(I)) * AMP, SIN(HARVEST(I)) * AMP)

            ENDIF

          ENDDO
        ENDDO
      ENDDO

      CALL WATCH_STOP(TICTOC)

      TIME(1) = TICTOC

      CALL WATCH_START(TICTOC)

      CALL ENFORCE_SYMMETRY(LS, LE, SPEC)

      CALL WATCH_STOP(TICTOC)

      TIME(2) = TICTOC

    END SUBROUTINE COMPUTE_SPECTRUM

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE MPI_SPLIT_TASK(N, NTASKS, I0, I1)

      ! DISTRIBUTE "N" POINTS AMONGST "NTASKS" PROCESSES AND, FOR EACH PROCESS, RETURN THE LOWEST AND HIGHEST INDEX.

      INTEGER(IPP),                        INTENT(IN)  :: N                          !< NUMBER OF POINTS TO BE SPLIT
      INTEGER(IPP),                        INTENT(IN)  :: NTASKS                     !< NUMBER OF MPI PROCESSES
      INTEGER(IPP), DIMENSION(0:NTASKS-1), INTENT(OUT) :: I0, I1                     !< 1ST/LAST INDEX
      INTEGER(IPP)                                     :: P                          !< COUNTER

      !------------------------------------------------------------------------------------------------------------------------------

      DO P = 0, NTASKS - 1
        I0(P) = 1 + INT( REAL(N, FPP) / REAL(NTASKS, FPP) * REAL(P, FPP) )
        I1(P) = INT( REAL(N, FPP) / REAL(NTASKS, FPP) * REAL(P + 1, FPP) )
      ENDDO

    END SUBROUTINE MPI_SPLIT_TASK

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(FPP) FUNCTION VK_PSDF(M, HURST)

      ! FUNCTION USED TO GENERATE VON KARMAN (AND EXPONENTIAL) RANDOM FIELDS

      REAL(FPP), INTENT(IN) :: M                  !< PRODUCT (WAVENUMBER*CL)**2
      REAL(FPP), INTENT(IN) :: HURST              !< HURST EXPONENT

      !-----------------------------------------------------------------------------------------------------------------------------

      VK_PSDF = (1._FPP + M)**(HURST + 1.5_FPP)

    END FUNCTION VK_PSDF

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(FPP) FUNCTION GS_PSDF(M, DUMMY)

      ! FUNCTION USED TO GENERATE GAUSSIAN RANDOM FIELDS

      REAL(FPP), INTENT(IN) :: M                  !< PRODUCT (WAVENUMBER*CL)**2
      REAL(FPP), INTENT(IN) :: DUMMY              !< THIS VARIABLE IS NOT USED

      !-----------------------------------------------------------------------------------------------------------------------------

      GS_PSDF = EXP(0.25_FPP * M)

    END FUNCTION GS_PSDF

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE ENFORCE_SYMMETRY(FS, FE, SPEC)

      ! MAKE SURE THE COMPLEX SPECTRUM HAS HERMITIAN SYMMETRY. SYMMETRY CONDITIONS ARE ENFORCED AT "I=1" AND "I=NYQUIST" ON THE YZ-PLANE;
      ! SYMMETRY ALONG X-DIRECTION IS IMPLICIT AS WE MAKE USE OF A COMPLEX-TO-REAL IFFT. THE ALGORITHM RUNS OVER EACH PROCESS CUT BY
      ! THE YZ-PLANE(S) STEP BY STEP, IN ORDER TO REQUIRE AS LESS MEMORY AS POSSIBLE (THIS RESULTS IN A PERFORMANCE PERNALTY SINCE MORE
      ! COMMUNICATIONS ARE NEEDED). SEE FIGURE 4 AND 5 OF PARDO-IGUZQUIZA AND CHICA-OLMO FOR A BETTER UNDERSTANDING OF THE SYMMETRY
      ! CONDITIONS.

      INTEGER(IPP),              DIMENSION(3),                                   INTENT(IN)    :: FS, FE                       !< LOCAL INDICES
      COMPLEX(FPP),              DIMENSION(FS(1):FE(1),FS(2):FE(2),FS(3):FE(3)), INTENT(INOUT) :: SPEC                         !< SPECTRUM
      COMPLEX(FPP), ALLOCATABLE, DIMENSION(:)                                                  :: BUFFER                       !< BUFFER FOR MESSAGING
      INTEGER(IPP)                                                                             :: C, I, J, K, P                !< COUNTERS
      INTEGER(IPP)                                                                             :: I0, J0, J1, K0, K1           !< VARIOUS INDICES
      INTEGER(IPP)                                                                             :: NTASKS, RANK, COLOR, SLICE   !< MPI STUFF
      INTEGER(IPP)                                                                             :: COMM, IERR, KEY              !< MPI STUFF
      INTEGER(IPP),              DIMENSION(3)                                                  :: NYQUIST                      !< INDEX NYQUIST WAVENUMBER
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                                                  :: SLICE2GLOBAL, JMAP, KMAP     !< ARRAYS FOR MESSAGING
      LOGICAL                                                                                  :: BOOL

      !-----------------------------------------------------------------------------------------------------------------------------

      ! INDEX OF NYQUIST WAVENUMBER
      DO I = 1, 3
        NYQUIST(I) = N(I) / 2 + 1
      ENDDO

      ! CONSIDER ONLY I-INDEX = 1 AND I-INDEX = NYQUIST(1)
      DO I = 1, NYQUIST(1), NYQUIST(1) - 1

        ! "TRUE" IF PROCESS CONTAINS POINT "ALPHA"
        BOOL = (FS(1) .EQ. 1) .AND. (FS(2) .EQ. 1) .AND. (FS(3) .EQ. 1)

        IF (BOOL) SPEC(1, 1, 1) = 0._FPP                                                         !< ALPHA MUST BE ZERO (ZERO-MEAN FIELD)

        ! "TRUE" IF PROCESS CONTAINS POINT "ALPHA"
        BOOL = (FS(1) .LE. I) .AND. (FE(1) .GE. I) .AND. (FS(2) .EQ. 1) .AND. (FS(3) .EQ. 1)

        IF (BOOL) SPEC(I, 1, 1) = REAL(SPEC(I, 1, 1), FPP)                                       !< ALPHA MUST BE REAL

        ! "TRUE" IF PROCESS CONTAINS POINT "PSI"
        BOOL = (FS(1) .LE. I) .AND. (FE(1) .GE. I) .AND. (FS(2) .LE. NYQUIST(2)) .AND. (FE(2) .GE. NYQUIST(2)) .AND.      &
               (FS(3) .LE. NYQUIST(3)) .AND. (FE(3) .GE. NYQUIST(3))

        IF (BOOL) SPEC(I, NYQUIST(2), NYQUIST(3)) = REAL(SPEC(I, NYQUIST(2), NYQUIST(3)), FPP)   !< PSI MUST BE REAL

        ! "TRUE" IF PROCESS CONTAINS POINT "XSI"
        BOOL = (FS(1) .LE. I) .AND. (FE(1) .GE. I) .AND. (FS(2) .EQ. 1) .AND. (FS(3) .LE. NYQUIST(3)) .AND. (FE(3) .GE. NYQUIST(3))

        IF (BOOL) SPEC(I, 1, NYQUIST(3)) = REAL(SPEC(I, 1, NYQUIST(3)), FPP)                     !< XSI MUST BE REAL

        ! "TRUE" IF PROCESS CONTAINS POINT "ETA"
        BOOL = (FS(1) .LE. I) .AND. (FE(1) .GE. I) .AND. (FS(2) .LE. NYQUIST(2)) .AND. (FE(2) .GE. NYQUIST(2)) .AND. (FS(3) .EQ. 1)

        IF (BOOL) SPEC(I, NYQUIST(2), 1) = REAL(SPEC(I, NYQUIST(2), 1), FPP)                     !< ETA MUST BE REAL

        ! "TRUE" ONLY FOR THOSE PROCESSES CONTAINING "I=1" OR "I=NYQUIST"
        BOOL = (FS(1) .LE. I) .AND. (FE(1) .GE. I)

        COLOR = 0

        IF (BOOL) COLOR = I

        ! FIRST SPLIT OF ORIGINAL COMMUNICATOR
        CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, COLOR, WORLD_RANK, SLICE, IERR)

        ! GLOBAL TO LOCAL MAPPING FOR I-INDEX
        I0 = I - FS(1) + 1

        ! ONLY PROCESSES BELONGING TO THE (YZ-)SLICE (I.E. WITH "COLOR=I") ENTER THIS IF STATEMENT
        IF (COLOR .EQ. I) THEN

          ! NUMBER OF PROCESSES IN THE NEW COMMUNICATOR
          CALL MPI_COMM_SIZE(SLICE, NTASKS, IERR)

          CALL MPI_COMM_RANK(SLICE, RANK, IERR)

          ALLOCATE(SLICE2GLOBAL(0:NTASKS-1))

          ! SLICE-RANK TO GLOBAL-RANK MAPPING
          SLICE2GLOBAL(RANK) = WORLD_RANK

          ! UPDATE GLOBAL-RANK MAPPING FOR ALL PROCESSES
          CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_INTEGER, SLICE2GLOBAL, 1, MPI_INTEGER, SLICE, IERR)

          ! LOOP OVER EACH PROCESS IN THE "SLICE" COMMUNICATOR: THOSE PROCESSES WHOSE DATA ARE NEEDED FOR COMPLEX CONJUGATION IN PROCESS
          ! "P" WILL BE PART OF A NEW COMMUNICATOR
          DO P = 0, NTASKS - 1

            ! ENTER STATEMENT IF PROCESS "P" HAS POINTS IN THE EPSILON*/THETA*/PHI* REGION
            IF ( (GE(2, SLICE2GLOBAL(P)) .GE. NYQUIST(2)) .AND. (GE(3, SLICE2GLOBAL(P)) .GE. NYQUIST(3)) ) THEN

              ! CORRESPONDING MIN/MAX INDICES IN THE EPSILON/THETA/PHI REGION
              J0 = N(2) - GE(2, SLICE2GLOBAL(P)) + 2
              J1 = N(2) - MAX(GS(2, SLICE2GLOBAL(P)), NYQUIST(2)) + 2

              K0 = N(3) - GE(3, SLICE2GLOBAL(P)) + 2
              K1 = N(3) - MAX(GS(3, SLICE2GLOBAL(P)), NYQUIST(3)) + 2

              ! "TRUE" IF CALLING PROCESS HAS DATA NEEDED BY PROCESS "P"
              BOOL = (( (FS(2) .LE. J0) .AND. (J0 .LE. FE(2)) ) .OR. ( (FS(2) .LE. J1) .AND. (J1 .LE. FE(2)) ) .OR.      &
                      ( (FS(2) .GE. J0) .AND. (FE(2) .LE. J1) ))                                               .AND.     &
                     (( (FS(3) .LE. K0) .AND. (K0 .LE. FE(3)) ) .OR. ( (FS(3) .LE. K1) .AND. (K1 .LE. FE(3)) ) .OR.      &
                      ( (FS(3) .GE. K0) .AND. (FE(3) .LE. K1) ))

              COLOR = 0

              ! CALLING PROCESS IS PART OF NEW COMMUNICATOR
              IF (BOOL) THEN
                COLOR = 1
                KEY   = RANK
              ENDIF

              ! ADD PROCESS "P" TO COMMUNICATOR AS IT WILL COLLECT DATA
              IF (P .EQ. RANK) THEN
                COLOR = 1
                KEY   = -1                                                  !< MAKE SURE IT GETS LOWEST (0) RANK VALUE
              ENDIF

              ! SPLIT "SLICE" COMMUNICATOR
              CALL MPI_COMM_SPLIT(SLICE, COLOR, KEY, COMM, IERR)

              ! ONLY PROCESSES IN THE NEW COMMUNICATOR EXECUTE STATEMENTS BELOW
              IF (COLOR .EQ. 1) THEN

                ! GATHER DATA FROM PROCESSES
                CALL GATHER_DATA(COMM, I0, SPEC, JMAP, KMAP, BUFFER)

                DO C = 1, SIZE(BUFFER)

                  ! "TRUE" IF CURRENT POINT IS INSIDE THE EPSILON/THETA/PHI REGION
                  BOOL = (JMAP(C) .GE. J0) .AND. (JMAP(C) .LE. J1) .AND. (KMAP(C) .GE. K0) .AND. (KMAP(C) .LE. K1)

                  ! INDICES WHERE COMPLEX-CONJUGATED VALUES MUST BE PLACED
                  J = N(2) - JMAP(C) + 2
                  K = N(3) - KMAP(C) + 2

                  ! UPDATE SPECTRUM (ONLY ROOT PROCESS IN COMMUNICATOR "COMM")
                  IF (BOOL .AND. (KEY .EQ. -1)) SPEC(I, J, K) = CONJG(BUFFER(C))

                ENDDO

              ENDIF

              CALL MPI_COMM_FREE(COMM, IERR)

            ENDIF
            ! END EPSILON*/THETA*/PHI* REGION

            ! ENTER STATEMENT IF PROCESS "P" HAS POINTS IN THE GAMMA* REGION
            IF ( (GS(2, SLICE2GLOBAL(P)) .LT. NYQUIST(2)) .AND. (GE(3, SLICE2GLOBAL(P)) .GT. NYQUIST(3)) ) THEN

              ! CORRESPONDING MIN/MAX INDICES IN THE GAMMA REGION
              J0 = N(2) - MIN(GE(2, SLICE2GLOBAL(P)), NYQUIST(2) - 1) + 2
              J1 = N(2) - GS(2, SLICE2GLOBAL(P)) + 2

              K0 = N(3) - GE(3, SLICE2GLOBAL(P)) + 2
              K1 = N(3) - MAX(GS(3, SLICE2GLOBAL(P)), NYQUIST(3) + 1) + 2

              ! "TRUE" IF CALLING PROCESS HAS DATA NEEDED BY PROCESS "P"
              BOOL = (( (FS(2) .LE. J0) .AND. (J0 .LE. FE(2)) ) .OR. ( (FS(2) .LE. J1) .AND. (J1 .LE. FE(2)) ) .OR.      &
                      ( (FS(2) .GE. J0) .AND. (FE(2) .LE. J1) ))                                               .AND.     &
                     (( (FS(3) .LE. K0) .AND. (K0 .LE. FE(3)) ) .OR. ( (FS(3) .LE. K1) .AND. (K1 .LE. FE(3)) ) .OR.      &
                      ( (FS(3) .GE. K0) .AND. (FE(3) .LE. K1) ))

              COLOR = 0

              ! CALLING PROCESS IS PART OF NEW COMMUNICATOR
              IF (BOOL) THEN
                COLOR = 1
                KEY   = RANK
              ENDIF

              ! ADD PROCESS "P" TO COMMUNICATOR AS IT WILL COLLECT DATA
              IF (P .EQ. RANK) THEN
                COLOR = 1
                KEY   = -1                                                  !< MAKE SURE IT GETS LOWEST (0) RANK VALUE
              ENDIF

              ! SPLIT "SLICE" COMMUNICATOR
              CALL MPI_COMM_SPLIT(SLICE, COLOR, KEY, COMM, IERR)

              ! ONLY PROCESSES IN THE NEW COMMUNICATOR EXECUTE STATEMENTS BELOW
              IF (COLOR .EQ. 1) THEN

                ! GATHER DATA FROM PROCESSES
                CALL GATHER_DATA(COMM, I0, SPEC, JMAP, KMAP, BUFFER)

                DO C = 1, SIZE(BUFFER)

                  ! "TRUE" IF CURRENT POINT IS INSIDE THE GAMMA REGION
                  BOOL = (JMAP(C) .GE. J0) .AND. (JMAP(C) .LE. J1) .AND. (KMAP(C) .GE. K0) .AND. (KMAP(C) .LE. K1)

                  ! INDICES WHERE COMPLEX-CONJUGATED VALUES MUST BE PLACED
                  J = N(2) - JMAP(C) + 2
                  K = N(3) - KMAP(C) + 2

                  ! UPDATE SPECTRUM (ONLY ROOT PROCESS IN COMMUNICATOR "COMM")
                  IF (BOOL .AND. (KEY .EQ. -1)) SPEC(I, J, K) = CONJG(BUFFER(C))

                ENDDO

              ENDIF

              CALL MPI_COMM_FREE(COMM, IERR)

            ENDIF
            ! END GAMMA* REGION

            ! ENTER STATEMENT IF PROCESS "P" HAS POINTS IN THE BETA* REGION
            IF ( (GE(2, SLICE2GLOBAL(P)) .GT. NYQUIST(2)) .AND. (GS(3, SLICE2GLOBAL(P)) .EQ. 1) ) THEN

              ! CORRESPONDING MIN/MAX INDICES IN THE BETA REGION
              J0 = N(2) - GE(2, SLICE2GLOBAL(P)) + 2
              J1 = N(2) - MAX(GS(2, SLICE2GLOBAL(P)), NYQUIST(2) + 1) + 2

              K0 = 1
              K1 = 1

              ! "TRUE" IF CALLING PROCESS HAS DATA NEEDED BY PROCESS "P"
              BOOL = (( (FS(2) .LE. J0) .AND. (J0 .LE. FE(2)) ) .OR. ( (FS(2) .LE. J1) .AND. (J1 .LE. FE(2)) ) .OR.      &
                      ( (FS(2) .GE. J0) .AND. (FE(2) .LE. J1) ))                                               .AND.     &
                     (( (FS(3) .LE. K0) .AND. (K0 .LE. FE(3)) ) .OR. ( (FS(3) .LE. K1) .AND. (K1 .LE. FE(3)) ) .OR.      &
                      ( (FS(3) .GE. K0) .AND. (FE(3) .LE. K1) ))

              COLOR = 0

              ! CALLING PROCESS IS PART OF NEW COMMUNICATOR
              IF (BOOL) THEN
                COLOR = 1
                KEY   = RANK
              ENDIF

              ! ADD PROCESS "P" TO COMMUNICATOR AS IT WILL COLLECT DATA
              IF (P .EQ. RANK) THEN
                COLOR = 1
                KEY   = -1                                                  !< MAKE SURE IT GETS LOWEST (0) RANK VALUE
              ENDIF

              ! SPLIT "SLICE" COMMUNICATOR
              CALL MPI_COMM_SPLIT(SLICE, COLOR, KEY, COMM, IERR)

              ! ONLY PROCESSES IN THE NEW COMMUNICATOR EXECUTE STATEMENTS BELOW
              IF (COLOR .EQ. 1) THEN

                ! GATHER DATA FROM PROCESSES
                CALL GATHER_DATA(COMM, I0, SPEC, JMAP, KMAP, BUFFER)

                DO C = 1, SIZE(BUFFER)

                  ! "TRUE" IF CURRENT POINT IS INSIDE THE BETA REGION
                  BOOL = (JMAP(C) .GE. J0) .AND. (JMAP(C) .LE. J1) .AND. (KMAP(C) .GE. K0) .AND. (KMAP(C) .LE. K1)

                  ! INDICES WHERE COMPLEX-CONJUGATED VALUES MUST BE PLACED
                  J = N(2) - JMAP(C) + 2
                  K = 1

                  ! UPDATE SPECTRUM (ONLY ROOT PROCESS IN COMMUNICATOR "COMM")
                  IF (BOOL .AND. (KEY .EQ. -1)) SPEC(I, J, K) = CONJG(BUFFER(C))

                ENDDO

              ENDIF

              CALL MPI_COMM_FREE(COMM, IERR)

            ENDIF
            ! END BETA* REGION

            ! ENTER STATEMENT IF PROCESS "P" HAS POINTS IN THE DELTA* REGION
            IF ( (GS(2, SLICE2GLOBAL(P)) .EQ. 1) .AND. (GE(3, SLICE2GLOBAL(P)) .GT. NYQUIST(3)) ) THEN

              ! CORRESPONDING MIN/MAX INDICES IN THE DELTA REGION
              J0 = 1
              J1 = 1

              K0 = N(3) - GE(3, SLICE2GLOBAL(P)) + 2
              K1 = N(3) - MAX(GS(3, SLICE2GLOBAL(P)), NYQUIST(3) + 1) + 2

              ! "TRUE" IF CALLING PROCESS HAS DATA NEEDED BY PROCESS "P"
              BOOL = (( (FS(2) .LE. J0) .AND. (J0 .LE. FE(2)) ) .OR. ( (FS(2) .LE. J1) .AND. (J1 .LE. FE(2)) ) .OR.      &
                      ( (FS(2) .GE. J0) .AND. (FE(2) .LE. J1) ))                                               .AND.     &
                     (( (FS(3) .LE. K0) .AND. (K0 .LE. FE(3)) ) .OR. ( (FS(3) .LE. K1) .AND. (K1 .LE. FE(3)) ) .OR.      &
                      ( (FS(3) .GE. K0) .AND. (FE(3) .LE. K1) ))

              COLOR = 0

              ! CALLING PROCESS IS PART OF NEW COMMUNICATOR
              IF (BOOL) THEN
                COLOR = 1
                KEY   = RANK
              ENDIF

              ! ADD PROCESS "P" TO COMMUNICATOR AS IT WILL COLLECT DATA
              IF (P .EQ. RANK) THEN
                COLOR = 1
                KEY   = -1                                                  !< MAKE SURE IT GETS LOWEST (0) RANK VALUE
              ENDIF

              ! SPLIT "SLICE" COMMUNICATOR
              CALL MPI_COMM_SPLIT(SLICE, COLOR, KEY, COMM, IERR)

              ! ONLY PROCESSES IN THE NEW COMMUNICATOR EXECUTE STATEMENTS BELOW
              IF (COLOR .EQ. 1) THEN

                ! GATHER DATA FROM PROCESSES
                CALL GATHER_DATA(COMM, I0, SPEC, JMAP, KMAP, BUFFER)

                DO C = 1, SIZE(BUFFER)

                  ! "TRUE" IF CURRENT POINT IS INSIDE THE BETA REGION
                  BOOL = (JMAP(C) .GE. J0) .AND. (JMAP(C) .LE. J1) .AND. (KMAP(C) .GE. K0) .AND. (KMAP(C) .LE. K1)

                  ! INDICES WHERE COMPLEX-CONJUGATED VALUES MUST BE PLACED
                  J = 1
                  K = N(3) - KMAP(C) + 2

                  ! UPDATE SPECTRUM (ONLY ROOT PROCESS IN COMMUNICATOR "COMM")
                  IF (BOOL .AND. (KEY .EQ. -1)) SPEC(I, J, K) = CONJG(BUFFER(C))

                ENDDO

              ENDIF

              CALL MPI_COMM_FREE(COMM, IERR)

            ENDIF
            ! END BETA* REGION


          ENDDO
          ! END LOOP OVER PROCESSES IN "SLICE" COMMUNICATOR

          DEALLOCATE(SLICE2GLOBAL)

        ENDIF
        ! END IF FOR COLOR = 1

        CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

        CALL MPI_COMM_FREE(SLICE, IERR)


      ENDDO
      ! END LOOP OVER INDEX "I"

    END SUBROUTINE ENFORCE_SYMMETRY

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE GATHER_DATA(COMM, I0, SPEC, JMAP, KMAP, BUFFER)

      ! GATHER DATA FROM ALL PROCESSES IN A COMMUNICATOR IN ROOT PROCESS. GLOBAL INDICES OF GATHERED DATA ARE RETURNED AS WELL.

      INTEGER(IPP),                                INTENT(IN)  :: COMM                             !< MPI COMMUNICATOR
      INTEGER(IPP),                                INTENT(IN)  :: I0                               !< LOCAL INDEX
      COMPLEX(FPP),              DIMENSION(:,:,:), INTENT(IN)  :: SPEC                             !< SPECTRUM
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:),     INTENT(OUT) :: JMAP, KMAP                       !< GLOBAL INDICES OF COLLECTED DATA
      COMPLEX(FPP), ALLOCATABLE, DIMENSION(:),     INTENT(OUT) :: BUFFER                           !< COLLECTED DATA
      COMPLEX(FPP), ALLOCATABLE, DIMENSION(:,:)                :: DUM                              !< LOCAL BUFFER
      INTEGER(IPP)                                             :: C, J, K, P                       !< COUNTERS
      INTEGER(IPP)                                             :: NPTS
      INTEGER(IPP)                                             :: NTASKS, RANK, IERR               !< MPI STUFF
      INTEGER(IPP), ALLOCATABLE, DIMENSION(:)                  :: COMM2GLOBAL, RECVCOUNTS, DISPLS  !< MPI STUFF

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL MPI_COMM_SIZE(COMM, NTASKS, IERR)

      CALL MPI_COMM_RANK(COMM, RANK, IERR)

      ALLOCATE(COMM2GLOBAL(0:NTASKS-1))

      ! COMM-RANK TO GLOBAL-RANK MAPPING
      COMM2GLOBAL(RANK) = WORLD_RANK

      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_INTEGER, COMM2GLOBAL, 1, MPI_INTEGER, COMM, IERR)

      ALLOCATE(RECVCOUNTS(0:NTASKS-1), DISPLS(0:NTASKS-1))

      RECVCOUNTS(RANK) = (GE(2, WORLD_RANK) - GS(2, WORLD_RANK) + 1) * (GE(3, WORLD_RANK) - GS(3, WORLD_RANK) + 1)

      ! NUMBER OF POINTS TO BE RECEIVED BY EACH PROCESS IN "COMM" COMMUNICATOR
      CALL MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_INTEGER, RECVCOUNTS, 1, MPI_INTEGER, COMM, IERR)

      ! TOTAL NUMBER OF POINTS TO BE RECEIVED
      NPTS = SUM(RECVCOUNTS)

      ALLOCATE(BUFFER(NPTS), JMAP(NPTS), KMAP(NPTS))

      ! DETERMINE THEIR RELATIVE DISPLACEMENT
      DISPLS(0) = 0

      DO P = 1, NTASKS - 1
        DISPLS(P) = DISPLS(P - 1) + RECVCOUNTS(P - 1)
      ENDDO

      ALLOCATE(DUM(SIZE(SPEC, 2), SIZE(SPEC, 3)))

      ! PREPARE DATA TO BE TRANSFERRED
      DO K = 1, SIZE(SPEC, 3)
        DO J = 1, SIZE(SPEC, 2)
          DUM(J, K) = SPEC(I0, J, K)
        ENDDO
      ENDDO

      ! COLLECT DATA INTO "BUFFER"
      CALL MPI_GATHERV(DUM, SIZE(DUM), COMPLEX_TYPE, BUFFER, RECVCOUNTS, DISPLS, COMPLEX_TYPE, 0, COMM, IERR)

      DEALLOCATE(DUM, DISPLS, RECVCOUNTS)

      C = 0

      ! CREATE MAPS FOR "J" AND "K" INDICES FOR ALL POINTS IN "BUFFER"
      DO P = 0, NTASKS - 1
        DO K = GS(3, COMM2GLOBAL(P)), GE(3, COMM2GLOBAL(P))
          DO J = GS(2, COMM2GLOBAL(P)), GE(2, COMM2GLOBAL(P))
            C       = C + 1
            JMAP(C) = J
            KMAP(C) = K
          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE(COMM2GLOBAL)

    END SUBROUTINE GATHER_DATA

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE SCARFLIB_FFT2
