MODULE m_psdf

  USE, NON_INTRINSIC :: m_scarflib_common, only: f_int, f_real, f_dble, pi

  IMPLICIT NONE

  PUBLIC

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTEGER(f_int) :: d

  REAL(f_real)   :: nu

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(f_real) FUNCTION fim_vk(x)

      ! Purpose:
      ! To compute the Von Karman (or exponential, if hurst = 0.5) PSDF for the FIM
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_real), INTENT(IN) :: x

      !-----------------------------------------------------------------------------------------------------------------------------

      fim_vk = 1._f_real / (1._f_real + x)**nu

    END FUNCTION fim_vk

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(f_real) FUNCTION fim_gs(x)

      ! Purpose:
      ! To compute the Gaussian PSDF for the FIM
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_real), INTENT(IN) :: x

      !-----------------------------------------------------------------------------------------------------------------------------

      fim_gs = EXP(-0.25_f_real * x)

    END FUNCTION fim_gs

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(f_real) FUNCTION fim_ud(x)

      ! Purpose:
      ! To compute a user-defined PSDF for the FIM
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_real), INTENT(IN) :: x

      !-----------------------------------------------------------------------------------------------------------------------------

      !fim_ud = ...

    END FUNCTION fim_ud

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(f_dble) FUNCTION srm_vk(x)

      ! Purpose:
      ! To compute Von Karman (or exponential, if hurst = 0.5) PSDF for the SRM. Input and result must be in double precision as
      ! requested by numerical integration routines.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_dble), INTENT(IN) :: x
      REAL(f_real)             :: z

      !-----------------------------------------------------------------------------------------------------------------------------

      z = x**2

      srm_vk = x * fim_vk(z)

      IF (d .eq. 3) srm_vk = srm_vk * x

    END FUNCTION srm_vk

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(f_dble) FUNCTION srm_gs(x)

      ! Purpose:
      ! To compute 2D and 3D Gaussian PSDF for the SRM. Input and result must be in double precision as requested by numerical
      ! integration routines.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_dble), INTENT(IN) :: x
      REAL(f_real)             :: z

      !-----------------------------------------------------------------------------------------------------------------------------

      z = x**2

      srm_gs = x * fim_gs(z)

      ! sqrt(pi) is introduced in order to have srm always < 1
      IF (d .eq. 3) srm_gs = srm_gs * x / sqrt(pi)

    END FUNCTION srm_gs

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(f_dble) FUNCTION srm_ud(x)

      ! Purpose:
      ! To compute 2D and 3D user-defined PSDF for the SRM. Input and result must be in double precision as requested by numerical
      ! integration routines.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(f_dble), INTENT(IN) :: x

      !-----------------------------------------------------------------------------------------------------------------------------

      ! here we may use the PSDF form defined for the FIM
      !srm_ud = ...

      ! a scaling factor may be introduced such that psdf always < 1
      !IF (d .eq. 3) srm_ud = srm_ud * ...

    END FUNCTION srm_ud

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


END MODULE m_psdf
