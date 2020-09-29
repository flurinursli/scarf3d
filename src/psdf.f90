MODULE m_psdf

  ! Copyright (c) 2020, Eidgenoessische Technische Hochschule Zurich, ETHZ.
  !
  ! Written by:
  ! Walter Imperatori (walter.imperatori@sed.ethz.ch)
  !
  ! All rights reserved.
  !
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
  !   To define power spectral density functions for the SCARF3D library.
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   04/05/20                  original version
  !

  USE, NON_INTRINSIC :: m_scarflib_common, only: f_int, f_real, f_dble, pi

  IMPLICIT NONE

  PUBLIC

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! this variable represents the hurst exponent term "nu + D/2"
  REAL(f_real) :: nu

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(f_real) FUNCTION fn_vk(x)

      ! Purpose:
      ! To compute the Von Karman (or exponential, if hurst = 0.5) PSDF
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_real), INTENT(IN) :: x                !< "k" or "sqrt(sum_{i=1}{D} k_i^2 * a_i^2)"

      !-----------------------------------------------------------------------------------------------------------------------------

      fn_vk = 1._f_real / (1._f_real + x**2)**nu

    END FUNCTION fn_vk

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(f_real) FUNCTION fn_gs(x)

      ! Purpose:
      ! To compute the Gaussian PSDF
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_real), INTENT(IN) :: x               !< "k" or "sqrt(sum_{i=1}{D} k_i^2 * a_i^2)"

      !-----------------------------------------------------------------------------------------------------------------------------

      fn_gs = EXP(-0.25_f_real * x**2)

    END FUNCTION fn_gs

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(f_real) FUNCTION fn_ud(x)

      ! Purpose:
      ! To compute a user-defined PSDF
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_real), INTENT(IN) :: x               !< "k" or "sqrt(sum_{i=1}{D} k_i^2 * a_i^2)"

      !-----------------------------------------------------------------------------------------------------------------------------

      !fn_ud = ...

    END FUNCTION fn_ud

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE m_psdf
