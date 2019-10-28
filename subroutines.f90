SUBROUTINE four1d(data,isign)
!-----------------------------------------------------------------------------------------
!
! Description:
!
!   External subroutine from Numerical Recipes, prepares input for 
!   FAST FOURIER TRANSFORM (FFT) or inverse (IFFT) computation
!
! References:
!
!   Numerical Recipes in Fortran90, pag. 1239-1240
!
! Notes:
!
!   i_sign=1 -> FFT, i_sign=-1 -> IFFT (see Numerical Recipes convention)
!
! Warning:
!
!   This subroutine works ONLY with 2**n number of points
!
! Authors: W. Imperatori (walter.imperatori@sed.ethz.ch)
!
! Version: Created on August 2009, v1.0
!

use constants; use precisions; use interfaces, only: fourrow, arth

implicit none

! input-output vector where complex transform will be stored
complex(rsp),dimension(:),intent(inout) :: data
! flag for FFT/IFFT
integer(isp),intent(in)                 :: isign
! locals
complex(rsp),dimension(:,:),allocatable :: dat,temp
complex(rdp),dimension(:),allocatable   :: w,wp
real(rdp),dimension(:),allocatable      :: theta
integer(isp)                            :: n,m1,m2,j
real(rdp)                               :: pi

!-----------------------------------------------------------------------------------------

n=size(data)                            

if ( iand(n,n-1) /=0 ) then
   print*, 'error in four1d'
endif   

call set_pi(pi)
pi = dble(2.) * pi

m1=2**ceiling(0.5*log(real(n,rsp))/0.693147)  
m2=n/m1

allocate(dat(m1,m2),theta(m1),w(m1),wp(m1),temp(m2,m1))

dat=reshape(data,shape(dat))

call fourrow(dat,isign)    !first transform   

theta=arth(0,isign,m1)*pi/n
wp=cmplx(-2.0*sin(0.5*theta)**2,sin(theta),rdp)
w=cmplx(1.0,0.0,rdp)

do j=2,m2
   w=w*wp+w
   dat(:,j)=dat(:,j)*w
enddo

temp=transpose(dat)

call fourrow(temp,isign)   !second transform
                                           
data=reshape(temp,shape(data))
                                           
deallocate(dat,w,wp,theta,temp)

END SUBROUTINE four1d

!=========================================================================================

SUBROUTINE fourrow(data,isign)
!-----------------------------------------------------------------------------------------
!
! Description:
!
!   External subroutine from Numerical Recipes, calculates FAST FOURIER
!   TRANSFORM (FFT) and its inverse (IFFT) for single-precision complex
!   data
!
! References:
!
!   Numerical Recipes in Fortran90, pag. 1235
!
! Notes:
!
!   i_sign=1 -> FFT, i_sign=-1 -> IFFT (see Numerical Recipes convention)
!
! Warning:
!
!   This subroutine works ONLY with 2**n number of points
!
! Authors: W. Imperatori (walter.imperatori@sed.ethz.ch)
!
! Version: Created on August 2009, v1.0
!
	
use constants; use precisions; use interfaces, only: swap

implicit none

! input/output array
complex(rsp),dimension(:,:),intent(inout) :: data
! flag for FFT/IFFT
integer(isp),intent(in)                   :: isign
! locals
integer(isp)                              :: n,i,istep,j,m,mmax,n2
real(rdp)                                 :: theta
complex(rsp), dimension(size(data,1))     :: temp
complex(rdp)                              :: w,wp
complex(rsp)                              :: ws
real(rdp)                                 :: pi

!-----------------------------------------------------------------------------------------

n=size(data,2)                      
                 
if ( iand(n,n-1) /= 0 ) then
   print*, 'error in fourrow'
endif

call set_pi(pi)

n2=n/2

j=n2

do i=1,n-2

   if (j > i) call swap(data(:,j+1),data(:,i+1))

   m=n2

      do
	     if (m < 2 .or. j < m) exit
		 j=j-m
		 m=m/2
	  enddo
	  
	  j=j+m
enddo

mmax=1

do
   if (n <= mmax) exit
   istep=2*mmax
   theta=pi/(isign*mmax)
   wp=cmplx(-2.0*sin(0.5*theta)**2,sin(theta),rdp)
   w=cmplx(1.0,0.0,rdp)

   do m=1,mmax
      ws=w

	     do i=m,n,istep
		    j=i+mmax
			temp=ws*data(:,j)
			data(:,j)=data(:,i)-temp
			data(:,i)=data(:,i)+temp
		 enddo

		 w=w*wp+w

   enddo

   mmax=istep

enddo
	
END SUBROUTINE fourrow

!=========================================================================================

FUNCTION arth(first,increment,n)
!-----------------------------------------------------------------------------------------
!
! Description:
!
!   Return an arithmetic progression, working with integer numbers
!
! References:
!
!   Numerical Recipes in Fortran90, pag. 1371
!
! Authors: W. Imperatori (walter.imperatori@sed.ethz.ch)
!
! Version: Created on August 2009, v1.0
!

use precisions

implicit none

integer(isp),intent(in)  :: first,increment,n
integer(isp),dimension(n):: arth
! local variables
integer(isp)             :: k,k2,temp
integer(isp),parameter   :: NPAR_ARTH=16,NPAR2_ARTH=8

!-----------------------------------------------------------------------------------------
                                  
if (n > 0) arth(1)=first

if (n <= NPAR_ARTH) then

   do k=2,n
      arth(k)=arth(k-1)+increment
   enddo

else

   do k=2,NPAR2_ARTH
      arth(k)=arth(k-1)+increment
   enddo
		
   temp=increment*NPAR2_ARTH
   k=NPAR2_ARTH
	
   do
      if (k >= n) exit
   	  k2=k+k
	  arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
	  temp=temp+temp
	  k=k2
   enddo
	
endif
	
END FUNCTION arth

!=========================================================================================

FUNCTION arth_d(first,increment,n) 
!-----------------------------------------------------------------------------------------
!
! Description:
!
!   Return an arithmetic progression, working with double precision numbers
!
! References:
!
!   Numerical Recipes in Fortran90, pag. 1371
!
! Authors: W. Imperatori (walter.imperatori@sed.ethz.ch)
!
! Version: Created on August 2009, v1.0
!

use precisions

implicit none

real(rdp), intent(in) :: first,increment 
integer(isp), intent(in) :: n 
real(rdp), dimension(n) :: arth_d 
integer(isp) :: k,k2 
real(rdp) :: temp 
integer(isp),parameter   :: NPAR_ARTH=16,NPAR2_ARTH=8

!-----------------------------------------------------------------------------------------

if (n > 0) arth_d(1)=first 
if (n <= NPAR_ARTH) then
   do k=2,n 
      arth_d(k)=arth_d(k-1)+increment
   enddo 
else
   do k=2,NPAR2_ARTH 
      arth_d(k)=arth_d(k-1)+increment
   enddo 
   temp=increment*NPAR2_ARTH 
   k=NPAR2_ARTH 
   do
      if (k >= n) exit 
      k2=k+k 
      arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k)) 
      temp=temp+temp 
      k=k2
   enddo 
endif

END FUNCTION arth_d

!=========================================================================================

SUBROUTINE swap(a,b)
!-----------------------------------------------------------------------------------------
!
! Description:
!
!   Swap the content of a and b
!
! References:
!
!   Numerical Recipes in Fortran90, pag. 1366-1367
!
! Authors: W. Imperatori (walter.imperatori@sed.ethz.ch)
!
! Version: Created on August 2009, v1.0
!

use precisions

implicit none

! arrays to be swapped
complex(rsp),dimension(:),intent(inout):: a,b
! local variable
complex(rsp),dimension(size(a))        :: dum

!-----------------------------------------------------------------------------------------
                                  
dum=a
a=b
b=dum
	
END SUBROUTINE swap

!=========================================================================================

FUNCTION gammaln(xx) 
!-----------------------------------------------------------------------------------------
!
! Description:
!
!   Compute gamma
!
! References:
!
!   Numerical Recipes in Fortran90
!
! Authors: W. Imperatori (walter.imperatori@sed.ethz.ch)
!
! Version: Created on August 2009, v1.0

use precisions; use interfaces, only: arth_d

implicit none

real(rsp), intent(in) :: xx 
real(rsp) :: gammaln
real(rdp) :: tmp,x
real(rdp) :: stp = 2.5066282746310005 
real(rdp), dimension(6) :: coef = (/dble(76.18009172947146),&
                                            dble(-86.50532032941677),dble(24.01409824083091),& 
                                            dble(-1.231739572450155),dble(0.1208650973866179e-2),& 
                                            dble(-0.5395239384953e-5)/)

!-----------------------------------------------------------------------------------------

x=xx 
tmp=x+dble(5.5) 
tmp=(x+dble(0.5))*log(tmp)-tmp 
gammaln=tmp+log(stp*(dble(1.000000000190015)+&
               sum(coef(:)/arth_d(x+dble(1.0),dble(1.0),size(coef))))/x) 
               
gammaln = exp(gammaln)               
                 
END FUNCTION gammaln

!=========================================================================================

SUBROUTINE poly_interp(x1a,x2a,x3a,ya,x1,x2,x3,y)
!------------------------------------------------------------------------
!
! Description:
!
!   Returns an interpolated value using polynomial interpolation (3D). If
!   input function is [2x2x2], bilinear interpolation in 3D is performed
!
! References:
!
!   Numerical Recipes in Fortran90
!
! Authors: W. Imperatori (walter.imperatori@sed.ethz.ch)
!
! Version: Created on August 2009, v1.0
!

use interfaces, only: polint

implicit none

! arrays for x,y coordinates	
real,dimension(:),intent(in)    :: x1a,x2a,x3a
! values of function at (x,y,z)
real,dimension(:,:,:),intent(in):: ya
! coordinates of interpolated point
real,intent(in)                 :: x1,x2,x3
! interpolated point
real,intent(out)                :: y
! counter, vector length 
integer                         :: k,j,m,w
! temporary arrays 
real,dimension(size(x1a))       :: yntmp
real,dimension(size(x2a))       :: ymtmp
real,dimension(size(x3a))       :: ywtmp

!-----------------------------------------------------------------------

! check for data 
if ( size(x1a) /= size(ya,1) ) print*,'error in POLY_INTERP'
if ( size(x2a) /= size(ya,2) ) print*,'error in POLY_INTERP'
if ( size(x3a) /= size(ya,3) ) print*,'error in POLY_INTERP'

m=size(x2a); w=size(x3a)

do k=1,w
   do j=1,m
      yntmp = ya(:,j,k)
      call polint(x1a,yntmp,x1,ymtmp(j))
   enddo
   call polint(x2a,ymtmp,x2,ywtmp(k))
enddo

call polint(x3a,ywtmp,x3,y)   


END SUBROUTINE poly_interp

!=========================================================================================

SUBROUTINE polint(xa,ya,x,y)
!-----------------------------------------------------------------------------------------
!
! Description:
!
!   Returns value through polynomial interpolation (1D), that means
!   the returned value is y=P(x), where P is polynomial of degree
!   N-1 (N is the number of points in xa)
!
! References:
!
!   Numerical Recipes in Fortran90
!
! Authors: W. Imperatori (walter.imperatori@sed.ethz.ch)
!
! Version: Created on August 2009, v1.0
!

implicit none

! arrays for x coordinate and y values
real,dimension(:),intent(in):: xa,ya
! coordinate of interpolated value
real,intent(in)             :: x
! interpolated value and error estimation
real,intent(out)            :: y
! counter and dummies (locals)
real                        :: dy
integer                     :: m,n,ns
real, dimension(size(xa))   :: c,d,den,ho
integer,dimension(1)        :: imin

!-----------------------------------------------------------------------------------------

! check for data consistency 
if ( size(xa) /= size(ya) ) print*,'error in POLINT'

n=size(xa)   

! initialize tableau
c=ya
d=ya

ho=xa-x

! find index closest to table entry
imin=minloc(abs(x-xa))
ns=imin(1)   

! initial approximation
y=ya(ns)
ns=ns-1

do m=1,n-1
   
   den(1:n-m)=ho(1:n-m)-ho(1+m:n)      
   
   ! abort if two xa are identical
   if (any(den(1:n-m) == 0.0))  print*,'error in POLINT'

   den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
   d(1:n-m)=ho(1+m:n)*den(1:n-m)
   c(1:n-m)=ho(1:n-m)*den(1:n-m)

   if (2*ns < n-m) then
      dy=c(ns+1)
   else
      dy=d(ns)
      ns=ns-1
   endif

   y=y+dy

enddo

END SUBROUTINE polint

!=========================================================================================