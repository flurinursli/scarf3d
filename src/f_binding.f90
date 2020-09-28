MODULE m_scarflib_c_binding

  ! This file is part of SCARF3D, version: 2.4
  !
  ! SCARF3D is free software: you can redistribute it and/or modify
  ! it under the terms of the GNU General Public License as published by
  ! the Free Software Foundation, either version 3 of the License, or
  ! (at your option) any later version.
  !
  ! SCARF3D is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ! GNU General Public License for more details.
  !
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

  INTEGER(f_int),                            PRIVATE :: n                     !< number of grid nodes (process view, unstructured grid)
  INTEGER(f_int), ALLOCATABLE, DIMENSION(:), PRIVATE :: dims                  !< number of grid nodes (process view, structured grid)
  LOGICAL,                                   PRIVATE :: structured            !< flag for structured/unstructured grid
  LOGICAL,                                   PRIVATE :: three_dim

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE struct_initialize(nd, fs, fe, ds, acf, cl, sigma, method, hurst, dh, poi, np, mute, taper, rescale, pad, nc, fc,  &
                                 alpha, beta, gamma) BIND(c, name="struct_initialize")

      ! Purpose:
      !   To allow C wrapper to interface with "scarf_initialize" written in FORTRAN (see m_scarflib for input arguments).
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !   13/07/20                  added rotation angles
      !

      INTEGER(c_int),                                   INTENT(IN) :: nd                   !< define problem dimension (2 or 3D)
      INTEGER(c_int), DIMENSION(nd),                    INTENT(IN) :: fs, fe
      REAL(c_real),                                     INTENT(IN) :: ds
      INTEGER(c_int),                                   INTENT(IN) :: acf
      REAL(c_real),   DIMENSION(nd),                    INTENT(IN) :: cl
      REAL(c_real),                                     INTENT(IN) :: sigma
      INTEGER(c_int),                         OPTIONAL, INTENT(IN) :: method
      REAL(c_real),                           OPTIONAL, INTENT(IN) :: hurst
      REAL(c_real),                           OPTIONAL, INTENT(IN) :: dh
      TYPE(c_ptr),                            OPTIONAL, INTENT(IN) :: poi
      INTEGER(c_int),                         OPTIONAL, INTENT(IN) :: np
      REAL(c_real),                           OPTIONAL, INTENT(IN) :: mute, taper
      INTEGER(c_int),                         OPTIONAL, INTENT(IN) :: rescale, pad
      REAL(c_real),   DIMENSION(nd),          OPTIONAL, INTENT(IN) :: nc, fc
      REAL(c_real),                           OPTIONAL, INTENT(IN) :: alpha, beta, gamma
      REAL(c_real),   DIMENSION(:,:), POINTER                      :: ptr => NULL()
      REAL(c_real),   DIMENSION(0,0), TARGET                       :: void

      !-----------------------------------------------------------------------------------------------------------------------------

      structured = .true.

      IF (nd .eq. 3) THEN
        three_dim = .true.
      ELSE
        three_dim = .false.
      ENDIF

      ALLOCATE(dims(nd))

      dims(:) = fe(:) - fs(:) + 1

      IF (PRESENT(poi)) THEN
        CALL c_f_pointer(poi, ptr, [nd, np])
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

    SUBROUTINE unstruct_initialize(nd, npts, dh, acf, cl, sigma, x, y, z, method, hurst, ds, poi, np, mute, taper, rescale, pad,  &
                                   nc, fc, alpha, beta, gamma) BIND(c, name="unstruct_initialize")

      ! Purpose:
      !   To allow C wrapper to interface with "scarf_initialize" written in FORTRAN (see m_scarflib for input arguments).
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !   13/07/20                  added rotation angles
      !

      INTEGER(c_int),                                    INTENT(IN) :: nd                   !< define problem dimension (2 or 3D)
      INTEGER(c_int),                                    INTENT(IN) :: npts
      REAL(c_real),                                      INTENT(IN) :: dh
      INTEGER(c_int),                                    INTENT(IN) :: acf
      REAL(c_real),   DIMENSION(nd),                     INTENT(IN) :: cl
      REAL(c_real),                                      INTENT(IN) :: sigma
      REAL(c_real),   DIMENSION(npts),                   INTENT(IN) :: x, y
      REAL(c_real),   DIMENSION(npts),         OPTIONAL, INTENT(IN) :: z
      INTEGER(c_int),                          OPTIONAL, INTENT(IN) :: method
      REAL(c_real),                            OPTIONAL, INTENT(IN) :: hurst
      REAL(c_real),                            OPTIONAL, INTENT(IN) :: ds
      TYPE(c_ptr),                             OPTIONAL, INTENT(IN) :: poi
      INTEGER(c_int),                          OPTIONAL, INTENT(IN) :: np
      REAL(c_real),                            OPTIONAL, INTENT(IN) :: mute, taper
      INTEGER(c_int),                          OPTIONAL, INTENT(IN) :: rescale, pad
      REAL(c_real),   DIMENSION(nd),           OPTIONAL, INTENT(IN) :: nc, fc
      REAL(c_real),                            OPTIONAL, INTENT(IN) :: alpha, beta, gamma
      REAL(c_real),   DIMENSION(:,:),  POINTER                      :: ptr => NULL()
      REAL(c_real),   DIMENSION(0,0),  TARGET                       :: void

      !-----------------------------------------------------------------------------------------------------------------------------

      structured = .false.

      IF (nd .eq. 3) THEN
        three_dim = .true.
      ELSE
        three_dim = .false.
      ENDIF

      n = npts

      IF (PRESENT(poi)) THEN
        CALL c_f_pointer(poi, ptr, [nd, np])
      ELSE
        ptr => void
      ENDIF

      CALL scarf_initialize(dh, acf, cl, sigma, x, y, z, method, hurst, ds, ptr, mute, taper, rescale, pad, nc, fc, alpha, beta, &
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
      REAL(c_real),   DIMENSION(:,:),   POINTER             :: f_struct_2d
      REAL(c_real),   DIMENSION(:,:,:), POINTER             :: f_struct_3d

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (structured .eqv. .false.) THEN
        CALL c_f_pointer(field, f_unstruct, [n])
        CALL scarf_execute(seed, f_unstruct, stats)
        NULLIFY(f_unstruct)
      ELSE
        IF (three_dim) THEN
          CALL c_f_pointer(field, f_struct_3d, [dims])
          CALL scarf_execute(seed, f_struct_3d, stats)
          NULLIFY(f_struct_3d)
        ELSE
          CALL c_f_pointer(field, f_struct_2d, [dims])
          CALL scarf_execute(seed, f_struct_2d, stats)
          NULLIFY(f_struct_2d)
        ENDIF
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

      IF (ALLOCATED(dims)) DEALLOCATE(dims)

    END SUBROUTINE finalize

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE io_one(nd, npts, field, fout, nwriters) BIND(c, name="io_one")

      ! Purpose:
      !   To allow C wrapper to interface with "scarf_io" written in FORTRAN (see m_scarflib for input arguments).
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      INTEGER(c_int),                                                    INTENT(IN) :: nd
      INTEGER(c_int),                DIMENSION(nd),                      INTENT(IN) :: npts
      TYPE(c_ptr),                                                       INTENT(IN) :: field
      CHARACTER(c_char),             DIMENSION(*),                       INTENT(IN) :: fout
      INTEGER(c_int),                                          OPTIONAL, INTENT(IN) :: nwriters
      CHARACTER(:),      ALLOCATABLE                                                :: filename
      INTEGER(f_int)                                                                :: i
      REAL(c_real),                  DIMENSION(:,:),   POINTER                      :: ptr_2d
      REAL(c_real),                  DIMENSION(:,:,:), POINTER                      :: ptr_3d

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (three_dim) THEN
        CALL c_f_pointer(field, ptr_3d, [dims])
      ELSE
        CALL c_f_pointer(field, ptr_2d, [dims])
      ENDIF

      i = 1

      filename = ""

      DO WHILE(fout(i) .ne. c_null_char)
        filename = filename // fout(i)
        i        = i + 1
      ENDDO

      IF (three_dim) THEN
        CALL scarf_io(npts, ptr_3d, filename, nwriters)
        NULLIFY(ptr_3d)
      ELSE
        CALL scarf_io(npts, ptr_2d, filename, nwriters)
        NULLIFY(ptr_2d)
      ENDIF

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
