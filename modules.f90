MODULE precisions
!
! Description:
!
!   Set precision to guarantee portability on different machines.
!
! Authors: W. Imperatori (walter.imperatori@sed.ethz.ch)
!
! Version: Created on August 2009, v1.0
!          Modified on January 2011, v1.1 - exploit Fortran2008 intrinsics
!

use, intrinsic :: iso_fortran_env

implicit none

! Real single precision
integer,parameter:: rsp = REAL32
! Real double precision
integer,parameter:: rdp = REAL64
! Integer single representation
integer,parameter:: isp = INT32
! Integer double representation 
integer,parameter:: idp = INT64

END MODULE precisions

!=========================================================================================

MODULE constants
!
! Description:
!
!   Define several constants used in the code
!
! Authors: W. Imperatori (walter.imperatori@sed.ethz.ch)
!
! Version: Created on August 2009, v1.0
!          Modified on January 2011, v1.1 - overloading for pi
!

use precisions

interface set_pi
   module procedure set_pi_r, set_pi_d
end interface set_pi

! Imaginary unit
complex(rsp),parameter:: zeta = cmplx(0.0,1.0)

contains

subroutine set_pi_r(pi)
   real(rsp),intent(out) :: pi 
   pi = 4. * atan(1.)
end subroutine set_pi_r

subroutine set_pi_d(pi)
   real(rdp),intent(out) :: pi
   pi = dble(4.) * atan(dble(1.))
end subroutine set_pi_d   


END MODULE constants

!=========================================================================================

MODULE interfaces

INTERFACE
   FUNCTION arth(first,increment,n)
      use precisions
      integer(isp),intent(in)  :: first,increment,n
      integer(isp),dimension(n):: arth
   END FUNCTION arth
END INTERFACE

INTERFACE 
   FUNCTION arth_d(first,increment,n)
      use precisions
      real(rdp), intent(in) :: first,increment 
      integer(isp), intent(in) :: n 
      real(rdp), dimension(n) :: arth_d
   END FUNCTION arth_d
END INTERFACE

INTERFACE four1d
   SUBROUTINE four1d(data,isign)
      use precisions
      implicit none
      complex(rsp),dimension(:),intent(inout):: data
      integer(isp),intent(in)                :: isign
   END SUBROUTINE four1d
END INTERFACE four1d

INTERFACE fourrow
   SUBROUTINE fourrow(data,isign)
      use precisions
      implicit none
      complex(rsp),dimension(:,:),intent(inout):: data
      integer(isp),intent(in)                  :: isign
   END SUBROUTINE fourrow
END INTERFACE fourrow

INTERFACE polint
   SUBROUTINE polint(xa,ya,x,y)
      use precisions
      implicit none
      real,dimension(:),intent(in):: xa,ya
      real,intent(in)             :: x
      real,intent(out)            :: y
   END SUBROUTINE polint
END INTERFACE polint

INTERFACE poly_interp
   SUBROUTINE poly_interp(x1a,x2a,x3a,ya,x1,x2,x3,y)
      use precisions
      implicit none	
      real,dimension(:),intent(in)    :: x1a,x2a,x3a  
      real,dimension(:,:,:),intent(in):: ya
      real,intent(in)                 :: x1,x2,x3
      real,intent(out)                :: y
   END SUBROUTINE poly_interp
END INTERFACE poly_interp 

INTERFACE swap
   SUBROUTINE swap(a,b)
      use precisions
      complex(rsp),dimension(:),intent(inout):: a,b
   END SUBROUTINE swap
END INTERFACE swap

END MODULE interfaces