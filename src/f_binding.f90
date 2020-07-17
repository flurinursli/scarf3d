MODULE m_scarflib_c_binding

  ! Purpose:
  !   To allow C wrappers to interface with the FORTRAN library (m_scarflib).
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   04/05/20                  original version
  !

  USE, INTRINSIC     :: iso_c_binding
  USE, NON_INTRINSIC :: m_scarflib
  USE, NON_INTRINSIC :: m_scarflib_common, only: c_real, f_int

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! make all interfaces accessible from the outside
  PUBLIC

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTEGER(f_int),               PRIVATE :: n                     !< number of grid nodes (process view, unstructured grid)
  INTEGER(f_int), DIMENSION(3), PRIVATE :: dims                  !< number of grid nodes (process view, structured grid)
  LOGICAL,                      PRIVATE :: structured            !< flag for structured/unstructured grid

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE struct_initialize(fs, fe, ds, acf, cl, sigma, method, hurst, dh, poi, np, mute, taper, rescale, pad, nc, fc, alpha, &
                                 beta, gamma) BIND(c, name="struct_initialize")

      ! Purpose:
      !   To allow C wrapper to interface with "scarf_initialize" written in FORTRAN (see m_scarflib for input arguments).
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !   13/07/20                  added rotation angles
      !

      INTEGER(c_int), DIMENSION(3),                     INTENT(IN) :: fs, fe
      REAL(c_real),                                     INTENT(IN) :: ds
      INTEGER(c_int),                                   INTENT(IN) :: acf
      REAL(c_real),   DIMENSION(3),                     INTENT(IN) :: cl
      REAL(c_real),                                     INTENT(IN) :: sigma
      INTEGER(c_int),                         OPTIONAL, INTENT(IN) :: method
      REAL(c_real),                           OPTIONAL, INTENT(IN) :: hurst
      REAL(c_real),                           OPTIONAL, INTENT(IN) :: dh
      TYPE(c_ptr),                            OPTIONAL, INTENT(IN) :: poi
      INTEGER(c_int),                         OPTIONAL, INTENT(IN) :: np
      REAL(c_real),                           OPTIONAL, INTENT(IN) :: mute, taper
      INTEGER(c_int),                         OPTIONAL, INTENT(IN) :: rescale, pad
      REAL(c_real),   DIMENSION(3),           OPTIONAL, INTENT(IN) :: nc, fc
      REAL(c_real),                           OPTIONAL, INTENT(IN) :: alpha, beta, gamma
      REAL(c_real),   DIMENSION(:,:), POINTER                      :: ptr => NULL()
      REAL(c_real),   DIMENSION(0,0), TARGET                       :: void

      !-----------------------------------------------------------------------------------------------------------------------------

      structured = .true.

      dims(:) = fe(:) - fs(:) + 1

      IF (PRESENT(poi)) THEN
        CALL c_f_pointer(poi, ptr, [3, np])
      ELSE
        ptr => void
      ENDIF

      CALL scarf_initialize(fs, fe, ds, acf, cl, sigma, method, hurst, dh, ptr, mute, taper, rescale, pad, nc, fc, alpha, beta, &
                            gamma)

      NULLIFY(ptr)

    END SUBROUTINE struct_initialize

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE unstruct_initialize(npts, x, y, z, dh, acf, cl, sigma, method, hurst, ds, poi, np, mute, taper, rescale, pad, nc,  &
                                   fc, alpha, beta, gamma) BIND(c, name="unstruct_initialize")

      ! Purpose:
      !   To allow C wrapper to interface with "scarf_initialize" written in FORTRAN (see m_scarflib for input arguments).
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !   13/07/20                  added rotation angles
      !

      INTEGER(c_int),                                    INTENT(IN) :: npts
      REAL(c_real),   DIMENSION(npts),                   INTENT(IN) :: x, y, z
      REAL(c_real),                                      INTENT(IN) :: dh
      INTEGER(c_int),                                    INTENT(IN) :: acf
      REAL(c_real),   DIMENSION(3),                      INTENT(IN) :: cl
      REAL(c_real),                                      INTENT(IN) :: sigma
      INTEGER(c_int),                          OPTIONAL, INTENT(IN) :: method
      REAL(c_real),                            OPTIONAL, INTENT(IN) :: hurst
      REAL(c_real),                            OPTIONAL, INTENT(IN) :: ds
      TYPE(c_ptr),                             OPTIONAL, INTENT(IN) :: poi
      INTEGER(c_int),                          OPTIONAL, INTENT(IN) :: np
      REAL(c_real),                            OPTIONAL, INTENT(IN) :: mute, taper
      INTEGER(c_int),                          OPTIONAL, INTENT(IN) :: rescale, pad
      REAL(c_real),   DIMENSION(3),            OPTIONAL, INTENT(IN) :: nc, fc
      REAL(c_real),                            OPTIONAL, INTENT(IN) :: alpha, beta, gamma
      REAL(c_real),   DIMENSION(:,:),  POINTER                      :: ptr => NULL()
      REAL(c_real),   DIMENSION(0,0),  TARGET                       :: void

      !-----------------------------------------------------------------------------------------------------------------------------

      structured = .false.

      n = npts

      IF (PRESENT(poi)) THEN
        CALL c_f_pointer(poi, ptr, [3, np])
      ELSE
        ptr => void
      ENDIF

      CALL scarf_initialize(x, y, z, dh, acf, cl, sigma, method, hurst, ds, ptr, mute, taper, rescale, pad, nc, fc, alpha, beta, &
                            gamma)

      NULLIFY(ptr)

    END SUBROUTINE unstruct_initialize

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE execute(seed, field, stats) BIND(c, name="execute")

      ! Purpose:
      !   To allow C wrapper to interface with "scarf_execute" written in FORTRAN (see m_scarflib for input/output arguments).
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(c_int),                           INTENT(IN)  :: seed
      TYPE(c_ptr),                              INTENT(OUT) :: field
      REAL(c_real),   DIMENSION(8),             INTENT(OUT) :: stats
      REAL(c_real),   DIMENSION(:),     POINTER             :: f_unstruct
      REAL(c_real),   DIMENSION(:,:,:), POINTER             :: f_struct

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (structured .eqv. .false.) THEN
        CALL c_f_pointer(field, f_unstruct, [n])
        CALL scarf_execute(seed, f_unstruct, stats)
        NULLIFY(f_unstruct)
      ELSE
        CALL c_f_pointer(field, f_struct, [dims])
        CALL scarf_execute(seed, f_struct, stats)
        NULLIFY(f_struct)
      ENDIF

    END SUBROUTINE execute

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE finalize() BIND(c, name="finalize")

      ! Purpose:
      !   To allow C wrapper to interface with "scarf_finalize" written in FORTRAN.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL scarf_finalize()

    END SUBROUTINE finalize

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE io_one(npts, field, fout, nwriters) BIND(c, name="io_one")

      ! Purpose:
      !   To allow C wrapper to interface with "scarf_io" written in FORTRAN (see m_scarflib for input arguments).
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(c_int),                DIMENSION(3),                       INTENT(IN) :: npts
      TYPE(c_ptr),                                                       INTENT(IN) :: field
      CHARACTER(c_char),             DIMENSION(*),                       INTENT(IN) :: fout
      INTEGER(c_int),                                          OPTIONAL, INTENT(IN) :: nwriters
      CHARACTER(:),      ALLOCATABLE                                                :: filename
      INTEGER(f_int)                                                                :: i
      REAL(c_real),                  DIMENSION(:,:,:), POINTER                      :: ptr

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL c_f_pointer(field, ptr, [dims])

      i = 1

      filename = ""

      DO WHILE(fout(i) .ne. c_null_char)
        filename = filename // fout(i)
        i        = i + 1
      ENDDO

      CALL scarf_io(npts, ptr, filename, nwriters)

      NULLIFY(ptr)

    END SUBROUTINE io_one

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE io_slice(n, axis, plane, field, fout) BIND(c, name="io_slice")

      ! Purpose:
      !   To allow C wrapper to interface with "scarf_io" written in FORTRAN (see m_scarflib for input arguments).
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(c_int),                DIMENSION(3),             INTENT(IN) :: n
      CHARACTER(c_char),                                       INTENT(IN) :: axis
      INTEGER(c_int),                                          INTENT(IN) :: plane
      TYPE(c_ptr),                                             INTENT(IN) :: field
      CHARACTER(c_char),             DIMENSION(*),             INTENT(IN) :: fout
      CHARACTER(len=1)                                                    :: fax
      CHARACTER(:),      ALLOCATABLE                                      :: filename
      INTEGER(f_int)                                                      :: i
      REAL(c_real),                  DIMENSION(:,:,:), POINTER            :: ptr

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL c_f_pointer(field, ptr, [dims])

      i = 1

      filename = ""

      DO WHILE(fout(i) .ne. c_null_char)
        filename = filename // fout(i)
        i        = i + 1
      ENDDO

      fax = axis

      CALL scarf_io(n, fax, plane, ptr, filename)

      NULLIFY(ptr)

    END SUBROUTINE io_slice

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE m_scarflib_c_binding
