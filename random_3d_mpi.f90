PROGRAM random_3d_RAM_MPI
!=========================================================================================
!
! Description:
!
!   This program generates a 3D medium with randomly distributed velocity heterogeneities.
!   Calculations occur in the spectral domain following the approach of Pardo-Iguzquiza et
!   Chica-Olmo (1994).
!   The random field can be characterised by a Gaussian, exponential or von Karman
!   autocorrelation function. The random field can be tapered and/or muted at the sides,
!   except the free-surface (z=0). Eventually, a reference point can be defined: the random
!   field can be tapered/muted around this point or it can be tapered/muted everywhere except
!   around the point.
!   The random field can be either written to disk or can be superimposed to a background
!   velocity model (1D or 3D). In the latter case, the background velocity model is resampled
!   on the grid of the random field.
!
!
!   MKS system of units is assumed.
!
!   A right-hand reference system is assumed:
!
!          *------------------*
!         /|                 /|
!      X / |                / |
!       /  |     Y         /  |   O indicates axes origin
!      O------------------*   |
!      |   *--------------|---*
!      |  /               |  /
!    Z | /                | /
!      |/                 |/
!      *------------------*
!
!   If a background velocity model is provided, it must cover the entire grid where the
!   random field is computed. If a 1D velocity model is provided, it must cover the whole
!   vertical axis. The code returns an error message if this condition is not met.
!
!   In a first step, the X axis is split amongst the MPI tasks. Spectrum and IFFTs are computed
!   on Y-Z planes. In a second step, data are redistributed on X-Y planes and final IFFTs
!   are performed along X direction. Normally, each single X-Y plane is assigned to a different
!   MPI task; however, if not enough RAM is available, the Y axis is divided into N chunks
!   and a MPI task owns the subregion X-Y/N.
!
!   The desired continuous standard deviation, specified in input, differs from the discrete
!   one (i.e. the standard deviation of the output random field) because of spectral truncation
!   (see Frenje et Julin, 2000, for further details). Furthermore, we apply a taper between
!   Knyquist/2 and Knyquist to avoid aliasing.
!
!   Input parameters:
!
!      acf       [char]        --> autocorrelation function, either "VK" (Von Karman) or "GS"
!                                  (Gaussian).
!      acf_cl    [3*real]      --> correlation length.
!      sigma     [real]        --> desired continuous standard deviation.
!      v         [real]        --> Hurst exponent. If acf=VK and v=0.5, an exponential auto-
!                                  correlation function is produced.
!      iseed     [int]         --> seed number to initialise the random generator.
!      output_N  [3*real]      --> number of points along X,Y,Z to be saved to disk.
!      output_dx [real]        --> grid-spacing of the grid where the random field is computed.
!                                  This is also the grid-spacing used for output.
!      output_Vp [char]        --> name of output file containing Vp values.
!      output_Vs [char]        --> name of output file containing Vs values.
!      RAM       [real]        --> memory available to each MPI task.
!      backgr_Vp [char]        --> background velocity model (Vp values).
!      backgr_Vs [char]        --> background velocity model (Vs values).
!      backgr_mode [char]      --> kind of background velocity: 't' for text, 'b' for binary,
!                                  'n' for none.
!      backgr_N  [3*int]       --> number of points along X,Y,Z in the background velocity model.
!                                  Used only for binary files (3D velocity models).
!      backgr_dx [3*real] 	   --> grid-spacing of the grid where the background velocity model
!                                  is defined. Used only for binary files (3D velocity models).
!      taper_side [2*int]      --> grid points to be tapered (1) and muted (2).
!      ref_point_xyz [3*real]  --> position of reference point in absolute coordinates.
!      ref_point_taper [2*int] --> grid points to be tapered (1) and muted (2).
!      ref_point_mode [int]    --> define how the tapering/muting is applied: 0 - around the
!                                  the reference point; 1 - everywhere except around the
!                                  reference point.
!
!   Output:
!
!      A single file with Vs values (either random field or random field superimposed to a
!      background velocity model) or two single files with Vp and Vs values (random field
!      superimposed to a background velocity model)
!
! Author: W. Imperatori (walter.imperatori@sed.ethz.ch)
!
! Version: Created on August 2011, v1.0
!          Modified on December 2011, v1.1 - 3D background velocity model added
!          Modified on February 2012, v1.2 - Gaussian ACF added
!          Modified on December 2017, v1.3 - Structure of input file modified, several check-
!                                            points introduced, more comments added
!          Modified on August 2018, v1.4   - Removed bug leading to wrong random field amplitude
!                                            when "output_N" is not power of 2
!
! Notes:
!       A) With OpenMPI 3.0, an IO error prevents writing results correctly to disk if the
!       Y axis is split into a large number of chunks. Adding the option "--mca io romio314"
!       to the mpirun command avoids such error (see also https://github.com/open-mpi/ompi/
!       issues/4336)
!
!       B) Please report any bug to walter.imperatori@sed.ethz.ch
!
!=========================================================================================

use mpi
use omp_lib

use precisions
use constants, only: set_pi
use interfaces, only: four1d

implicit none

!include "mpif.h"

! constants
real(rsp) :: pi

! spectral stuff
real(rsp),dimension(3)                    :: dk, K_nyq, K_max
real(rsp)                                 :: sigma, v, exponent, hanning, num, gamma, Kr
real(rsp)                                 :: gammaln, Kr_nyq, Kr_max, scaling
real(rdp)                                 :: PWR  ! DO WE NEED DOUBLE PRECISION FOR GAUSSIAN ACF??
real(rsp),allocatable,dimension(:)        :: Kx, Ky, Kz
integer(isp),dimension(3)                 :: N
complex(rsp),allocatable,dimension(:,:)   :: SPEC, SPEC_Z
complex(rsp),allocatable,dimension(:,:,:) :: SPEC_STORE

! random generator stuff
integer(isp)                          :: iseed, seed_size
integer(isp),allocatable,dimension(:) :: tmp_seed
real(rsp),allocatable,dimension(:,:)  :: R

! correlation function stuff
character(len=2)       :: acf
real(rsp),dimension(3) :: acf_cl

! background model stuff
character(len=99)         :: backgr_Vp, backgr_Vs
character(len=1)          :: backgr_mode
integer(isp)              :: n_layers
integer(isp),dimension(3) :: backgr_N
real(rsp),dimension(3)    :: backgr_dx
real(rsp),allocatable,dimension(:) :: depth1D, vs1D, vp1D

! output random model stuff
character(len=99)         :: output_Vp, output_Vs
!integer(isp)              :: output_mode
integer(isp),dimension(3) :: output_N
real(rsp)                 :: output_dx

! reference point stuff
integer(isp)              :: ref_point_mode
integer(isp),dimension(2) :: ref_point_taper
real(rsp),dimension(3)    :: ref_point_xyz

! tapering stuff
integer(isp),dimension(2) :: taper_side

! MPI stuff
integer(isp)                            :: rank, ntasks, err, errcode, y_chunks, blocks
integer(isp)                            :: c_block, filep, status(MPI_STATUS_SIZE), nwriters
integer(isp)                            :: MCW, MCW_GR, WRT_GR, WRT_COMM, mypoints
integer(isp),allocatable,dimension(:)   :: counts, displacement, z_level, writers, npts_ychunk
integer(isp),allocatable,dimension(:)   :: proc, recvcounts, displ
integer(isp),allocatable,dimension(:,:) :: y_level

! statistics stuff
real(rsp)                          :: loc_min, loc_max, glob_min, glob_max, compos_Vp_min
real(rsp)                          :: compos_Vp_max, compos_Vs_min, compos_Vs_max
real(rdp)                          :: var, mu, npts
real(rdp),allocatable,dimension(:) :: varar, muar, nar

! counters
integer(isp) :: i, j, k, l, u, i_min, i_max, i_size, j_first, j_last

! timing stuff
real(rdp) :: exec_time, tic, tac, tic2, tac2, io_tic, io_tac

! IO stuff
logical                              :: lwrite = .false.
logical,dimension(5)                 :: lprint
real(rsp),allocatable,dimension(:,:) :: IOar

! system stuff
real(rsp) :: RAM, task_memory, By2Gy

!-----------------------------------------------------------------------------------------

! Initialise MPI
call MPI_INIT(err)

! Duplicate communicator (i.e. MPI_COMM_WORLD into MCW)
call MPI_COMM_DUP(MPI_COMM_WORLD,MCW,err)

! Get rank number
call MPI_COMM_RANK(MCW,rank,err)

! Get number of tasks
call MPI_COMM_SIZE(MCW,ntasks,err)

! Set barrier
call MPI_BARRIER(MCW,err)

exec_time = MPI_Wtime()

! Default values
backgr_Vp = 'none'

! pi
call set_pi(pi)

! set factor to convert from Byte to GByte
By2Gy = 1. / 1024. / 1024. / 1024.

! Master rank reads input file
if (rank .eq. 0) then

   open(1,file='input.inp',status='old')
   read(1,*); read(1,*) acf                ! auto-correlation function             2 char
   read(1,*); read(1,*) acf_cl             ! correlation length, 					3 real
   read(1,*); read(1,*) sigma              ! standard deviation, 					1 real
   read(1,*); read(1,*) v                  ! Hurst exponent,     					1 real
   read(1,*); read(1,*) iseed              ! seed number,        					1 integer
   read(1,*); read(1,*) output_N           ! points to save,     		 			3 integer
   read(1,*); read(1,*) output_dx          ! grid spacing,       					1 real
   read(1,*); read(1,*) output_Vp          ! filename for Vs,    		  			99 char
   read(1,*); read(1,*) output_Vs          ! filename for Vp,    		   			99 char
   !read(1,*) output_mode        ! describe output mode,			        1 integer
   read(1,*); read(1,*) RAM                ! available memory per task, 			1 real
   read(1,*); read(1,*) backgr_Vp          ! b.v.m. for Vp,                        99  char
   read(1,*); read(1,*) backgr_Vs          ! b.v.m. for Vs,                        99  char
   read(1,*); read(1,*) backgr_mode        ! specify b.v.m. type,                  1  char
   read(1,*); read(1,*) backgr_N           ! points in background velocity model,  3 integer
   read(1,*); read(1,*) backgr_dx          ! grid spacing of b.v.m.,               3 real
   read(1,*); read(1,*) taper_side         ! specify tapering of sides,            2 integer
   read(1,*); read(1,*) ref_point_xyz      ! position of r.p. and tapering,        3 real
   read(1,*); read(1,*) ref_point_taper
   read(1,*); read(1,*) ref_point_mode

   close(1)

   print*,'Summary of input parameters'
   print*,'***************************************************'
   print*,'Auto-correlation function: ',acf
   print*,'Correlation length (m): ',acf_cl
   print*,'Std. dev.: ',sigma
   print*,'Hurst exponent (v): ',v
   print*,'Seed number: ',iseed
   print*,'Points in output model: ',output_N
   print*,'Grid-step (m): ',output_dx
   print*,'Output file for Vp: ',trim(output_Vp)
   print*,'Output file for Vs: ',trim(output_Vs)
   !print*,'Output type: ',output_mode
   print*,'Available RAM: ',RAM
   print*,'Backgr model for Vp: ',trim(backgr_Vp)
   print*,'Backgr model for Vs: ',trim(backgr_Vs)
   print*,'Backgr model type: ',backgr_mode
   print*,'Point in backgr. model: ',backgr_N
   print*,'Grid-step in backgr model: ',backgr_dx
   print*,'Taper and Muting: ',taper_side
   print*,'Reference point location: ',ref_point_xyz
   if (any(ref_point_xyz .ne. -1.)) then
      print*,'Tapering and muting for ref-point: ', ref_point_taper
      print*,'Type of tapering/muting: ',ref_point_mode
   endif
   print*,'***************************************************'
   print*,''

endif

call MPI_BARRIER(MCW,err)

! Broadcast input-file parameters
call MPI_BCAST(acf,2,MPI_CHARACTER,0,MCW,err)
call MPI_BCAST(acf_cl,3,MPI_REAL,0,MCW,err)
call MPI_BCAST(output_dx,1,MPI_REAL,0,MCW,err)
call MPI_BCAST(sigma,1,MPI_REAL,0,MCW,err)
call MPI_BCAST(v,1,MPI_REAL,0,MCW,err)
call MPI_BCAST(iseed,1,MPI_INTEGER,0,MCW,err)
call MPI_BCAST(output_N,3,MPI_INTEGER,0,MCW,err)
call MPI_BCAST(output_Vp,99,MPI_CHARACTER,0,MCW,err)
call MPI_BCAST(output_Vs,99,MPI_CHARACTER,0,MCW,err)
call MPI_BCAST(RAM,1,MPI_REAL,0,MCW,err)
!call MPI_BCAST(output_mode,1,MPI_INTEGER,0,MCW,err)
call MPI_BCAST(backgr_Vp,99,MPI_CHARACTER,0,MCW,err)
call MPI_BCAST(backgr_Vs,99,MPI_CHARACTER,0,MCW,err)
call MPI_BCAST(backgr_mode,1,MPI_INTEGER,0,MCW,err)
call MPI_BCAST(backgr_N,3,MPI_INTEGER,0,MCW,err)
call MPI_BCAST(backgr_dx,3,MPI_REAL,0,MCW,err)
call MPI_BCAST(taper_side,2,MPI_INTEGER,0,MCW,err)
call MPI_BCAST(ref_point_xyz,3,MPI_REAL,0,MCW,err)
call MPI_BCAST(ref_point_taper,2,MPI_INTEGER,0,MCW,err)
call MPI_BCAST(ref_point_mode,1,MPI_INTEGER,0,MCW,err)

! here below we do some checking

! minimum wavenumber must be smaller than corner wavenumber, i.e. N*dx > 2*pi*acf_cl
if (any(output_dx*output_N .le. 2.*pi*acf_cl)) then
   if (rank .eq. 0) then
      print*,'Minimum wavenumber must be smaller than corner wavenumber'
      print*,'Either decrease acf_cl or increase N (or dx)'
      print*,'Aborting'
   endif
   call MPI_BARRIER(MCW,err)
   !call MPI_ABORT(MCW,errcode,err)
endif

! check grid spacing. We follow the suggestion of Frenje & Juhlin (2000), eq. 8. However,
! since we apply a taper between Knyquist/2 and Knyquist, we apply a somewhat more
! conservative assumption (i.e. dx ≤ a/2 instead of dx ≤ a)
if (any(output_dx .gt. acf_cl/2.)) then
   if (rank .eq. 0) then
      print*,'Grid-spacing must be smaller than acf_cl/2'
      print*,'Either increase acf_cl or decrease dx'
      print*,'Aborting'
   endif
   call MPI_BARRIER(MCW,err)
   !call MPI_ABORT(MCW,errcode,err)
endif

! if a background model is given, check whether it is large enough (3D) or load it into
! memory (1D)
if (backgr_mode .ne. 'n') call check_domain

! set number of points to work in spectral domain (FFT)
N = output_N

! check power of 2 requirement for FFT. Increase N if necessary
do i = 1,3
   exponent = log(real(N(i))) / log(2.0)
   if (exponent /= nint(exponent))  N(i) = 2**(ceiling(exponent))
enddo

! compute delta-k in each direction
dk = (1. / real(N)) * (2.*pi / output_dx)

! compute Nyquist wavenumber
K_nyq = pi/output_dx

! set maximum wavenumber to half of Nyquist wavenumber. Spectrum values above this will
! be tapered
K_max = K_nyq / 2.

! compute radial Nyquist and maximum wavenumbers
Kr_nyq = sqrt(K_nyq(1)**2 + K_nyq(2)**2 + K_nyq(3)**2)
Kr_max = Kr_nyq / 2.

allocate(Kx(N(1)),Ky(N(2)),Kz(N(3)))

! define wavenumber vectors
Kx = (/0.0,(dk(1)*i,i=1,N(1)/2)/); Kx = (/Kx(1:N(1)/2+1),(dk(1)*i,i=N(1)/2-1,1,-1)/)
Ky = (/0.0,(dk(2)*i,i=1,N(2)/2)/); Ky = (/Ky(1:N(2)/2+1),(dk(2)*i,i=N(2)/2-1,1,-1)/)
Kz = (/0.0,(dk(3)*i,i=1,N(3)/2)/); Kz = (/Kz(1:N(3)/2+1),(dk(3)*i,i=N(3)/2-1,1,-1)/)

! print some more info
if(rank .eq. 0) then
   print*,'Minimum wavenumbers (dk): ',dk
   print*,'Corner wavenumbers: ',1/acf_cl
   print*,'Nyquist wavenumbers: ',K_nyq
   print*,'Upper wavenumbers before tapering: ',K_max
endif

! Split loop over x-axis amongst MPI tasks. Some tasks may get larger x-intervals than
! others
i_min = 1 + nint( real(N(1)/2+1) / real(ntasks) * real(rank) )
i_max = nint( real(N(1)/2+1) / real(ntasks) * real(rank+1) )
i_size = i_max - i_min + 1

! Rough estimate of max memory consumption (GB) per task
! Add complex arrays for spectrum (SPEC_STORE and SPEC)
task_memory = 2. * kind(SPEC) * (output_N(2)*output_N(3)*i_size + N(2)*N(3))

! Add real array for random phase (R)
task_memory = task_memory + kind(R) * N(2)*N(3)

! memory needed in GB
task_memory = task_memory * By2Gy

! stop execution if we are above RAM
if (task_memory .gt. RAM) then
   if (rank .eq. 0) print*,'Program requires much memory than available... stopping'
   if (rank .eq. 0) print*,'Required memory (GB): ',task_memory
   call MPI_BARRIER(MCW,err)
   call MPI_ABORT(MCW,errcode,err)
endif

if(rank .eq. 0) print*,''
if(rank .eq. 0) write(*,'(X,A,F5.1)') 'Estimated memory per MPI task (GB): ',task_memory
if(rank .eq. 0) print*,''
if(rank .eq. 0) print*,'Loop for computing spectrum on 2D array (Y-Z)'

! compute part of power spectrum independent on wavenumber
! VON KARMAN (or EXPONENTIAL if v=0.5)
if (acf .eq. 'VK') then
   gamma = gammaln(v)
   num = (sigma**2) * (acf_cl(1)*acf_cl(2)*acf_cl(3)) * ((2.*sqrt(pi))**3) * gammaln(v+3./2.)

! GAUSSIAN
elseif (acf .eq. 'GS') then
   num = (sigma**2) * (acf_cl(1)*acf_cl(2)*acf_cl(3)) * sqrt(pi**3)
endif

! compute scaling factor
scaling = 1. / sqrt(N(1) * N(2) * N(3) * output_dx**3)

! prepare random numbers generator
call random_seed(size=seed_size)
allocate(tmp_seed(seed_size))

! Timing for computing spectrum
tic = MPI_Wtime()

! allocate array for random numbers
allocate(R(N(2),N(3)))

! allocate array for spectrum
allocate(SPEC(N(2),N(3)))

! Loop along x-axis to fill 2D (Y-Z) array
do i = i_min,i_max

   ! Update seed with x-level index and initialise random generator accordingly. This way
   ! results do not depend on the number of tasks
   tmp_seed = iseed + i
   call random_seed(put=tmp_seed)

   ! Rough estimate of total time needed to compute spectrum. Only values for master task
   ! are considered
   if(i .eq. (i_min + 1)) then
      tac = MPI_Wtime(); tac = tac - tic
      if(rank .eq. 0) print*,'Loop will be terminated in about: ',nint(tac * i_size),' sec'
      if(rank .eq. 0) print*,''
   endif

   !if(rank .eq. 0) print*,'Working on YZ array, x-index: ',i - i_min + 1,' of: ',i_max - i_min + 1

   ! generate random numbers
   call random_number(R)

   tic2 = MPI_Wtime()

   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,j,Kr,PWR,hanning)
   do k = 1,N(3)
      do j = 1,N(2)

         Kr = sqrt( (Kx(i)**2 * acf_cl(1)**2) + (Ky(j)**2 * acf_cl(2)**2) + (Kz(k)**2 * acf_cl(3)**2) )

         ! VON KARMAN (or EXPONENTIAL if v=0.5)
         if (acf .eq. 'VK') then
            PWR = num / ( gamma * (1 + Kr**2)**(v + 3./2.) )
         ! GAUSSIAN
         elseif (acf .eq. 'GS') then
            PWR = exp(0.25*Kr**2); PWR = num / PWR
         endif

         ! now we apply a separable window taper (Hann) to reduce risk of aliasing at high
         ! wave-numbers. Tapering occurs between Knyq/2 (i.e. Kmax) and Knyq

         ! x-direction
         if (Kx(i) .gt. K_max(1)) then
            hanning = 0.5 * (1. - cos(2. * pi * Kx(i) / K_nyq(1)))
         else
            hanning = 1.
         endif

         ! y-direction
         if (Ky(j) .gt. K_max(2)) then
            hanning = 0.5 * (1. - cos(2. * pi * Ky(j) / K_nyq(2))) * hanning
         else
            hanning = 1. * hanning
         endif

         ! z-direction
         if (Kz(k) .gt. K_max(3)) then
            hanning = 0.5 * (1. - cos(2. * pi * Kz(k) / K_nyq(3))) * hanning
         else
            hanning = 1. * hanning
         endif

         ! apply taper
         PWR = PWR * hanning

         ! putting together amplitude spectrum and random phase
         SPEC(j,k) = cmplx(sqrt(real(PWR))*cos(R(j,k)*2*pi),sqrt(real(PWR))*sin(R(j,k)*2*pi))

      enddo
   enddo
   !$OMP END PARALLEL DO

   tac2 = MPI_Wtime(); tac2 = tac2 - tic2
   if ( (rank .eq. 0) .and. (i .eq. i_min) ) print*,'{Spectrum computed in ',real(tac2),' sec}'


   ! define symmetry conditions and nugget effect
   if ( (i .eq. 1) .or. (i .eq. N(1)/2+1) ) then
      SPEC(1,1) = 0.0
      SPEC(N(2)/2+1,1) = real(SPEC(N(2)/2+1,1))
      SPEC(1,N(3)/2+1) = real(SPEC(1,N(3)/2+1))
      SPEC(N(2)/2+1,N(3)/2+1) = real(SPEC(N(2)/2+1,N(3)/2+1))

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k)
      do k=2,N(3)/2
         SPEC(1,N(3)-k+2) = conjg(SPEC(1,k))
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j)
      do j=2,N(2)/2
         SPEC(N(2)-j+2,1) = conjg(SPEC(j,1))
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k)
      do k=2,N(3)/2
         SPEC(N(2)/2+1,N(3)-k+2) = conjg(SPEC(N(2)/2+1,k))
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j)
      do j=2,N(2)/2
         SPEC(N(2)-j+2,N(3)/2+1) = conjg(SPEC(j,N(3)/2+1))
      enddo
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,k)
      do k=2,N(3)/2
         do j=2,N(2)/2
            SPEC(N(2)-j+2,N(3)-k+2) = conjg(SPEC(j,k))
            SPEC(N(2)-j+2,N(3)/2-k+2) = conjg(SPEC(j,N(3)/2+k))
         enddo
      enddo
      !$OMP END PARALLEL DO
   endif

   tic2 = MPI_Wtime()

   ! Inverse FFT along Y (no cache misses -> fast)
   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k)
   do k = 1,N(3)
      call four1d(SPEC(:,k),-1)
   enddo
   !$OMP END PARALLEL DO

   ! Inverse FFT along Z (cache misses -> slow)
   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j)
   do j = 1,N(2)
      call four1d(SPEC(j,:),-1)
   enddo
   !$OMP END PARALLEL DO

   tac2 = MPI_Wtime(); tac2 = tac2 - tic2
   if ( (rank .eq. 0) .and. (i .eq. i_min) ) print*,'{IFFTs computed in ',real(tac2),' sec}'


   if (.not.allocated(SPEC_STORE)) allocate(SPEC_STORE(output_N(2),output_N(3),i_size))

   SPEC_STORE(:,:,i - i_min + 1) = SPEC(1:output_N(2),1:output_N(3))

enddo

call MPI_BARRIER(MCW,err)

! deallocate some arrays
deallocate(R,SPEC)

! elapsed time for computing spectrum
tac = MPI_Wtime(); tac = tac - tic

call MPI_REDUCE(tac,tic,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MCW,err)

if(rank .eq. 0) print*,'Loop to fill 2D (Y-Z) array executed in: ',real(tic/real(ntasks)),' sec'
if(rank .eq. 0) print*,''

! deallocate some arrays
deallocate(Kx,Ky,Kz,tmp_seed)

! Here below we compute how many pieces ('y_chunks') the Y-axis should be divided into. The
! total number of blocks (output_N(3)*y_chunks) are assigned to each MPI task in a round-robin
! fashion

! allocate arrays
allocate(counts(0:ntasks-1),displacement(0:ntasks-1))

! store in 'counts' the number of Y-Z arrays (i.e. X-slices) for each MPI process
call MPI_ALLGATHER(i_size,1,MPI_INTEGER,counts,1,MPI_INTEGER,MCW,err)

! compute stride between each X-block. Basically, we are measuring the distance of each
! block from X=0
displacement(0) = 0
do i = 1,ntasks - 1
   displacement(i) = counts(i - 1) + displacement(i - 1)
enddo

! remove R...
task_memory = task_memory - kind(R) * N(2)*N(3) * By2Gy
! ...and SPEC from memory requirement
task_memory = task_memory - 2. * kind(SPEC) * N(2)*N(3) * By2Gy

! compute 'y_chunks', i.e. the amount of pieces the Y axis is split into in order to not
! exceed the available memory
y_chunks = ceiling(2. * kind(SPEC) * output_N(2)*(N(1) + i_size) * By2Gy / (RAM - task_memory))
!y_chunks = 8

if (rank .eq. 0) then
   write(*,'(X,A,X,I2,X,A)'),'Y-axis will be divided into',y_chunks,'chunk(s)'
   print*,''
endif

! total number of blocks in the Y-Z plane
blocks = y_chunks * output_N(3)

! compute points along Y for each block
allocate(npts_ychunk(blocks))

do i = 1,blocks
   npts_ychunk(i) = output_N(2) / y_chunks
   if (mod(i,y_chunks) .eq. 0) npts_ychunk(i) = output_N(2) / y_chunks + mod(output_N(2),y_chunks)
enddo

! assign each block to each task in a round-robin fashion
allocate(proc(blocks))
do j = 0,ntasks - 1
   do i = j,blocks - 1,ntasks
      proc(i + 1) = j
   enddo
enddo

! allocate vectors where we'll store Z-level and Y-level for each MPI task
allocate(z_level(0:ntasks - 1),y_level(0:ntasks - 1,2))

! initialise variables for statistics
mu = 0.; var = 0.; npts = 0.
loc_min = huge(1.); loc_max = -huge(1.)
glob_min = huge(1.); glob_max = -huge(1.)
compos_Vs_min = huge(1.); compos_Vs_max = -huge(1.)
compos_Vp_min = huge(1.); compos_Vp_max = -huge(1.)

! initialise number of points per task
mypoints = 0

lprint(:) = .true.

! start timing
tic = MPI_Wtime()

! Now we loop over all blocks to perform missing IFFT along X-direction.
! Note: some MPI tasks may not write data to disk, depending on 'proc'
do l = 1,output_N(3)

   ! Partial timing: get a rough estimate to code completion
   if ( (rank .eq. 0) .and. lprint(1) .and. (lprint(3) .eqv. .false.) ) then
      tac = MPI_Wtime(); tac = tac - tic
      print*,'Program will complete in about: ',nint(tac*blocks/c_block),' sec'
      print*,''
      lprint(1) = .false.
   endif

   tic2 = MPI_Wtime()

   j_first = 0; j_last = 0

   ! Loop over the pieces the Y-axis has been divided into
   do u = 1,y_chunks

      ! compute current block number
      c_block = (l - 1) * y_chunks + u

      ! numer of points along Y for current task
      if (rank .eq. proc(c_block)) mypoints = npts_ychunk(c_block)

      ! first index along Y
      !j_first = (u - 1) * npts_ychunk(c_block) + 1
      j_first = j_last + 1

      ! last index along Y
      !j_last = u * npts_ychunk(c_block)
      j_last = j_last + npts_ychunk(c_block)

      ! allocate array where results from previous calculations will be stored. Shape of
      ! 'SPEC_Z' is irrelevant to all tasks except to task owing current block
      if ( .not.allocated(SPEC_Z) .and. (rank .ne. proc(c_block)) ) then
         allocate(SPEC_Z(npts_ychunk(c_block),N(1)))
         SPEC_Z = cmplx(0.,0.)
      endif

      ! correct allocation for current task
      if ( .not.allocated(SPEC_Z) .and. (rank .eq. proc(c_block)) ) then
         allocate(SPEC_Z(mypoints,N(1)))
         SPEC_Z = cmplx(0.,0.)
      endif

      ! allocate temporary array to store task-specific partial results
      allocate(SPEC(npts_ychunk(c_block),i_size))

      ! read previous partial results
      do i = 1,i_size
         SPEC(:,i) = SPEC_STORE(j_first:j_last,l,i)
      enddo

      ! elements received
      allocate(recvcounts(0:ntasks - 1))
      recvcounts = counts * npts_ychunk(c_block)

      ! displacement relative to SPEC_Z
      allocate(displ(0:ntasks - 1))
      displ = displacement * npts_ychunk(c_block)

      ! Gather partial results to a specific MPI task
      call MPI_GATHERV(SPEC,size(SPEC),MPI_COMPLEX,SPEC_Z,recvcounts,displ,MPI_COMPLEX,  &
                       proc(c_block),MCW,err)

      deallocate(SPEC)
      deallocate(recvcounts, displ)

      ! At this point 'SPEC_Z' contains the whole [npts_ychunk(c_block),N(1)] data.

      ! Here we store Z-level and Y-segment (Y first/last index) for each MPI task. They
      ! will be used during interpolation
      z_level(proc(c_block))   = l
      y_level(proc(c_block),1) = j_first
      y_level(proc(c_block),2) = j_last


      ! Once all MPI tasks have 'SPEC_Z' loaded with data or there are no blocks left,
      ! continue with calculations (otherwise keep looping over blocks to load data). Note
      ! that all tasks enter the instruction set below
      if ( (proc(c_block) .eq. (ntasks - 1)) .or. (c_block .eq. blocks) ) then

         ! Timing
         tac2 = MPI_Wtime(); tac2 = tac2 - tic2
         if ( (rank .eq. 0) .and. lprint(2) ) then
             print*,'{Data for IO gathered in ',real(tac2),' sec}'
             print*,''
             lprint(2) = .false.
         endif

         forall (i = 1:mypoints, j = 2:N(1)/2)
            SPEC_Z(i,N(1) - j + 2) = conjg(SPEC_Z(i,j))
         end forall

         ! IFFT over X-direction
         !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
         do i = 1,mypoints
            call four1d(SPEC_Z(i,:),-1)
         enddo
         !$OMP END PARALLEL DO

         ! apply scaling factor. The random field is now in the real part of 'SPEC_Z'
         do j = 1,output_N(1)
            SPEC_Z(:,j) = cmplx( real(SPEC_Z(:,j)) * scaling, 0. )
         enddo

         ! Do we need to apply a taper and/or mute sides of the random field?
         if(any(taper_side .ne. 0)) then
            !if(rank .eq. 0) print*,'Tapering and muting edges'
            !tic2 = MPI_Wtime()
            call tapering
            !tac2 = MPI_Wtime(); tac2 = tac2 - tic2
            !if(rank .eq. 0) print*,'{Taper applied in ',real(tac2),' sec}'
         endif

         ! Do we need to apply a taper and/or mute a region around a PoI?
         if(all(ref_point_xyz(1:3) .ge. 0)) then
            !if(rank .eq. 0) print*,'Tapering and muting point of interest'
            !tic2 = MPI_Wtime()
            call taper_region
            !tac2 = MPI_Wtime(); tac2 = tac2 - tic2
            !if(rank .eq. 0) print*,'{Taper PoI applied in ',real(tac2),' sec}'
         endif

         ! Complete previous writing process (only tasks in WRT_COMM enter this IF)
         if(lwrite) then

            ! finish writing results and close file
            call MPI_FILE_WRITE_ORDERED_END(filep,IOar,status,err)
            call MPI_FILE_CLOSE(filep,err)

            io_tac = MPI_Wtime(); io_tac = io_tac - io_tic
            if ( (rank .eq. 0) .and. lprint(3) ) then
               print*,'{Data written in ',real(io_tac),' sec}'
               print*,''
               lprint(3) = .false.
            endif

            ! deallocate array used to write result to disk
            deallocate(IOar)

            ! Release group and communicator created for writing
            call MPI_GROUP_FREE(WRT_GR,err)
            call MPI_COMM_FREE(WRT_COMM,err)

            ! reset flag
            lwrite = .false.
         endif

         ! Here below we create a sub-group and communicator for MPI I/O, i.e. we setup a
         ! group made by only those tasks with data to be written to disk

         ! find number of tasks that must be write to disk. Normally, this is equal to 'ntasks',
         ! unless we are toward the end of the blocks
         nwriters = proc(c_block) + 1

         ! allocate array with rank of IO tasks
         allocate(writers(nwriters))
         do i = 0,proc(c_block)
            writers(i + 1) = i
         enddo

         ! create group
         call MPI_COMM_GROUP(MCW,MCW_GR,err)
         call MPI_GROUP_INCL(MCW_GR,nwriters,writers,WRT_GR,err)
         call MPI_COMM_CREATE(MCW,WRT_GR,WRT_COMM,err)

         ! tasks that do not participate to IO can deallocate array
         if (rank .gt. proc(c_block)) deallocate(SPEC_Z)

         ! Only tasks belonging to the group can execute the instructions below
         if(rank .le. proc(c_block)) then

            ! do some statistics on current block of data
            do j = 1,output_N(1)
               do i = 1,mypoints
                  call online_variance(real(SPEC_Z(i,j)), var, mu, npts)
               enddo
            enddo

            ! min / max
            loc_min = min( minval(real(SPEC_Z(:,1:output_N(1)))), loc_min )
            loc_max = max( maxval(real(SPEC_Z(:,1:output_N(1)))), loc_max )

            ! Read, resample and add values from a background velocity model
            if(backgr_mode .ne. 'n') then
               !if(rank .eq. 0) print*,'Resampling background model'
               tic2 = MPI_Wtime()
               call resample_deterministic_model
               tac2 = MPI_Wtime(); tac2 = tac2 - tic2
               if ( (rank .eq. 0) .and. lprint(4) ) then
                  print*,'{Background model read&resampled in ',real(tac2),' sec}'
                  print*,''
                  lprint(4) = .false.
               endif
            endif

            !if (rank .eq. 0) print*,'Writing files to disk'

            ! allocate array for writing to disk
            allocate(IOar(output_N(1),mypoints))

            io_tic = MPI_Wtime()

            ! start non-blocking IO
            call write_to_disk

            ! complete writing to disk immediately if we are at last block
            if (c_block .eq. blocks) then

               call MPI_FILE_WRITE_ORDERED_END(filep,IOar,status,err)
               call MPI_FILE_CLOSE(filep,err)

               ! deallocate array used to write result to disk
               deallocate(IOar)

               ! release group and communicator created for writing
               call MPI_GROUP_FREE(WRT_GR,err)
               call MPI_COMM_FREE(WRT_COMM,err)

            endif

            deallocate(SPEC_Z)

            !tac2 = MPI_Wtime(); tac2 = tac2 - tic2
            !if (rank .eq. 0) print*,'{Files written in: ',real(tac2),'sec}'

         endif  ! end block for writing to disk

         deallocate(writers)

         ! reset number of points for current task
         mypoints = 0

      endif  ! end block for assigning blocks to MPI tasks

      ! make sure 'SPEC_Z' is deallocated
      if ( allocated(SPEC_Z) .and. (mypoints .eq. 0) ) deallocate(SPEC_Z)

   enddo ! end loop over Y-chunks

enddo ! end loop over Z-levels

! deallocate pending arrays
deallocate(SPEC_STORE)

! deallocate arrays for 1D velocity model
if (backgr_mode .eq. 't') deallocate(depth1D,vs1D,vp1D)

! retrieve info above min/max values of random field
call MPI_REDUCE(loc_min,glob_min,1,MPI_REAL,MPI_MIN,0,MCW,err)
call MPI_REDUCE(loc_max,glob_max,1,MPI_REAL,MPI_MAX,0,MCW,err)

if (rank .eq. 0) then
   print*,'Minimum value of random field: ',glob_min
   print*,'Maximum value of random field: ',glob_max
endif

! retrieve info above min/max values of composed velocity model
if (backgr_mode .ne. 'n') then

   call MPI_REDUCE(compos_Vs_min,glob_min,1,MPI_REAL,MPI_MIN,0,MCW,err)
   call MPI_REDUCE(compos_Vs_max,glob_max,1,MPI_REAL,MPI_MAX,0,MCW,err)

   if (rank .eq. 0) then
      print*,'Minimum value for Vs: ',glob_min
      print*,'Maximum value for Vs: ',glob_max
   endif

   if (backgr_Vp .ne. 'none') then

      call MPI_REDUCE(compos_Vp_min,glob_min,1,MPI_REAL,MPI_MIN,0,MCW,err)
      call MPI_REDUCE(compos_Vp_max,glob_max,1,MPI_REAL,MPI_MAX,0,MCW,err)

      if (rank .eq. 0) then
         print*,'Minimum value for Vp: ',glob_min
         print*,'Maximum value for Vp: ',glob_max
      endif

   endif
endif

! compute variance for each task
var = var / (npts - 1.)

! allocate stat arrays
allocate(varar(ntasks),muar(ntasks),nar(ntasks))

! gather stat parameters
call MPI_ALLGATHER(var, 1, MPI_DOUBLE_PRECISION, varar, 1, MPI_DOUBLE_PRECISION, MCW, err)
call MPI_ALLGATHER(mu, 1, MPI_DOUBLE_PRECISION, muar, 1, MPI_DOUBLE_PRECISION, MCW, err)
call MPI_ALLGATHER(npts, 1, MPI_DOUBLE_PRECISION, nar, 1, MPI_DOUBLE_PRECISION, MCW, err)

! assign starting values
var = varar(1); mu = muar(1); npts = nar(1)

! final statistics
do i = 2,ntasks
   call parallel_variance(varar(i), muar(i), nar(i), var, mu, npts)
enddo

if (rank .eq. 0) then
   print*,'Average of random field: ',real(mu)
   print*,'Std. dev. of random field: ',real(sqrt(var))
   print*,''
endif

! deallocate stat arrays
deallocate(varar, muar, nar)

tac = MPI_Wtime(); exec_time = tac - exec_time

if (rank .eq. 0) print*,'Total elapsed time: ',real(exec_time),' sec'

if(rank .eq. 0) then
   print*,''
   print*,'Progam completed'
endif

call MPI_BARRIER(MCW,err)

call MPI_FINALIZE(err)

CONTAINS

!=========================================================================================

SUBROUTINE check_domain
!
! Description:
!
!    This subroutine accomplish the following tasks:
!
!    1) It verifies that the background velocity model is larger(3D case) / deeper(1D case)
!       than the desired output composite model
!
!    2) Load and broadcast velocity data (1D case)
!
!    The subroutine is not called if no background velocity model is specified (i.e. if
!    'backgr_mode' = 'n')
!
! Author: W. Imperatori (walter.imperatori@sed.ethz.ch)
!
! Version: Created on August 2011, v1.0
!          Modified on December 2017, v1.3 - Some variables were renamed to reflect changes
!                                            in the main program
!

implicit none

integer(isp) :: i
real(rsp)    :: xmax, ymax, zmax, x, y, z

!-------------------------------------------------------

! Check whether 1D deterministic velocity model is deep enough to encompass the output grid
if (trim(backgr_mode) .eq. 't') then

   ! only master task reads input file
   if(rank .eq. 0) then
      open(1,file=trim(backgr_Vs),status='old',form='formatted')
      read(1,*) n_layers
   endif

   ! broadcast info
   call MPI_BCAST(n_layers,1,MPI_INTEGER,0,MCW,err)

   if (n_layers .lt. 2) then
      if (rank .eq. 0) then
         print*,'ERROR, at least two layers are required'
         print*,'Number of layers ',n_layers
      endif
      call MPI_BARRIER(MCW,err)
      call MPI_ABORT(MCW,errcode,err)
   endif

   ! allocate some arrays
   allocate(depth1D(n_layers),vs1D(n_layers),vp1D(n_layers))

   if(rank .eq. 0) then
      do i = 1,n_layers
         read(1,*) depth1D(i),vp1D(i),vs1D(i)
      enddo
      close(1)
   endif

   ! broadcast medium parameters
   call MPI_BCAST(depth1D,n_layers,MPI_REAL,0,MCW,err)
   call MPI_BCAST(vp1D,n_layers,MPI_REAL,0,MCW,err)
   call MPI_BCAST(vs1D,n_layers,MPI_REAL,0,MCW,err)

   zmax = (output_N(3) - 1) * output_dx

   ! check whether background model is deep enough
   if ( depth1D(n_layers) .lt. zmax ) then
      if (rank .eq. 0) then
         print*,'Maximum depth of input 1D model must be larger than ', zmax
         print*,'Execution aborted'
      endif
      call MPI_BARRIER(MCW,err)
      call MPI_ABORT(MCW,errcode,err)
   endif

endif

! check whether background model encompasses output grid
if (backgr_mode .eq. 'b') then

   xmax = (output_N(1) - 1) * output_dx
   ymax = (output_N(2) - 1) * output_dx
   zmax = (output_N(3) - 1) * output_dx

   x = (backgr_N(1) - 1) * backgr_dx(1)
   y = (backgr_N(2) - 1) * backgr_dx(2)
   z = (backgr_N(3) - 1) * backgr_dx(3)

   if( (xmax .ge. x) .or. (ymax .ge. y) .or. (zmax .ge. z) ) then
      if(rank .eq. 0) then
         print*,'3D deterministic velocity model must be encompass output grid'
         print*,'Execution aborted'
      endif
      call MPI_BARRIER(MCW,err)
      call MPI_ABORT(MCW,errcode,err)
   endif

endif

END SUBROUTINE check_domain

!=========================================================================================

SUBROUTINE tapering
!
! Description:
!
!    This subroutine taper and mute the random field close to model sides, except at Z=0
!
! Author: W. Imperatori (walter.imperatori@sed.ethz.ch)
!
! Version: Created on August 2011, v1.0
!          Modified on December 2017, v1.3 - Some variables were renamed to reflect changes
!                                            in the main program
!

implicit none

integer(isp) :: taper, mute, i, k, l
real(rsp)    :: ratio, hanning

!-----------------------------------------------------------------------------------------

mute = taper_side(2)

! adjust input taper value
taper = taper_side(1) + 1

! mute in the X-direction, both sides
do k = 1,mute
   SPEC_Z(:,k) = 0.
enddo
l = output_N(1) - mute + 1
do k = l,output_N(1)
   SPEC_Z(:,k) = 0.
enddo

! taper in the X-direction, both sides
do k = (mute + 1),(taper + mute)
   ratio = real(k - mute) / real(taper)
   hanning = 1. - cos(0.5 * pi * ratio)
   SPEC_Z(:,k) = real(SPEC_Z(:,k)) * hanning
enddo
l = output_N(1) - mute
do k = (l - taper + 1),l
   ratio = real(k - l) / real(taper)
   hanning = 1. - cos(0.5*pi*ratio)
   SPEC_Z(:,k) = real(SPEC_Z(:,k)) * hanning
enddo

! mute along the Y-direction, both sides. Muting occurs only if chunk is at either ends of
! Y-axis
do i = y_level(rank,1),y_level(rank,2)
   if (i .le. mute) SPEC_Z(i - y_level(rank,1) + 1,:) = 0.
enddo
l = output_N(2) - mute + 1
do i = y_level(rank,1),y_level(rank,2)
   if(i .ge. l) SPEC_Z(i - y_level(rank,1) + 1,:) = 0.
enddo

! taper along the Y-direction, both sides. Tapering occurs only if chunk is at either ends
! of Y-axis
do i = y_level(rank,1),y_level(rank,2)
   if( (i .ge. (mute + 1)) .and. (i .le. (taper + mute)) ) then
      ratio = real(i - mute) / real(taper)
      hanning = 1. - cos(0.5*pi*ratio)
      SPEC_Z(i - y_level(rank,1) + 1,:) = real(SPEC_Z(i - y_level(rank,1) + 1,:)) * hanning
   endif
enddo
l = output_N(2) - mute
do i = y_level(rank,1),y_level(rank,2)
   if( (i .ge. (l - taper + 1)) .and. (i .le. l) ) then
      ratio = real(i - l) / real(taper)
      hanning = 1. - cos(0.5*pi*ratio)
      SPEC_Z(i - y_level(rank,1) + 1,:) = real(SPEC_Z(i - y_level(rank,1) + 1,:)) * hanning
   endif
enddo

! mute and taper along Z-direction
l = output_N(3) - mute
if( (z_level(rank) .ge. (l - taper + 1)) .and. (z_level(rank) .le. l) ) then
   ratio = real(z_level(rank) - l) / real(taper)
   hanning = 1. - cos(0.5*pi*ratio)
   SPEC_Z(:,:) = real(SPEC_Z(:,:)) * hanning
elseif( (z_level(rank) .ge. (l + 1)) .and. (z_level(rank) .le. output_N(3)) ) then
   SPEC_Z(:,:) = 0.
endif

END SUBROUTINE tapering

!=========================================================================================

SUBROUTINE taper_region
!
! Description:
!
!    This subroutine taper and mute the random field around a specific point as indicated by
!    ref_point_xyz, ref_point_taper, ref_point_mode.
!
! Author: W. Imperatori (walter.imperatori@sed.ethz.ch)
!
! Version: Created on August 2011, v1.0
!          Modified on December 2017, v1.3 - Some variables were renamed to reflect changes
!                                            in the main program
!

implicit none

integer(isp)              :: taper, mute, i, k, reverse
integer(isp),dimension(3) :: point
real(rsp)                 :: ratio, hanning, dist, dist_bndr

!-----------------------------------------------------------------------------------------

! determine PoI position in grid-points
point(1) = nint(ref_point_xyz(1) / output_dx) + 1
point(2) = nint(ref_point_xyz(2) / output_dx) + 1
point(3) = nint(ref_point_xyz(3) / output_dx) + 1

! radius (in grid points) of the homogeneous/heterogeneous region around the PoI
mute = ref_point_taper(2)

! size of taper around the PoI
taper = ref_point_taper(1)

! flag to choose between taper and muting a region around point of interest (0) or
! the whole medium except a region around point of interest (1)
reverse = ref_point_mode

if (reverse .eq. 0) then

   do k = y_level(rank,1),y_level(rank,2)

      do i = 1,output_N(1)

         dist = sqrt( real ((i - point(1))**2 + (k - point(2))**2 + (z_level(rank)- point(3))**2) )

         if( (nint(dist) .gt. mute) .and. (nint(dist) .le. (mute + taper)) ) then

            ratio = (dist - mute) / real(taper)
            hanning = 1. - cos(0.5*pi*ratio)

            SPEC_Z(k - y_level(rank,1) + 1,i) = real(SPEC_Z(k - y_level(rank,1) + 1,i)) * hanning

         elseif (nint(dist) .le. mute) then

            SPEC_Z(k - y_level(rank,1) + 1,i) = 0.

         endif

      enddo
   enddo

elseif(reverse .eq. 1) then

   do k = y_level(rank,1),y_level(rank,2)

      do i = 1,output_N(1)

         dist = sqrt( real ((i - point(1))**2 + (k - point(2))**2 + (z_level(rank)- point(3))**2) )

         if( (nint(dist) .gt. mute) .and. (nint(dist) .le. (mute + taper)) ) then

            ratio = (dist - mute) / real(taper)
            hanning = 1. - (1. - cos(0.5*pi*ratio))

            SPEC_Z(k - y_level(rank,1) + 1,i) = real(SPEC_Z(k - y_level(rank,1) + 1,i)) * hanning

         elseif (nint(dist) .gt. (mute + taper)) then

            SPEC_Z(k - y_level(rank,1) + 1,i) = 0.

         endif

      enddo
   enddo

endif

END SUBROUTINE taper_region

!=========================================================================================

SUBROUTINE resample_deterministic_model
!
! Description:
!
!    This subroutine resample the deterministic velocity model (either 1D or 3D) on the same
!    grid of the random field.
!    Resampling occurs also for Vp values (if provided via backgr_Vp)
!
! Author: W. Imperatori (walter.imperatori@sed.ethz.ch)
!
! Version: Created on August 2011, v1.0
!          Modified on December 2017, v1.3 - Some variables were renamed to reflect changes
!                                            in the main program.
!

use interfaces, only: poly_interp,polint

implicit none

! background velocity model stuff
integer(isp)                       :: j_min, j_max, j_npts, k_min, k_max, i_min, i_max
integer(isp)                       :: i_npts, ij_npts, ijk_npts

! counters
integer(isp) :: i, j

! MPI stuff
integer(isp)                  :: fileq, rl, filetype
integer(kind=mpi_offset_kind) :: offset, disp
integer(isp),dimension(3)     :: starts, subsizes


! velocity model resampling stuff
real(rsp)                              :: Vs, Vp, xo, yo, zo, vpi, vsi
real(rsp),dimension(2)                 :: x, y, z, vsl, vpl

real(rsp),dimension(n_layers)          :: tmp
real(rsp),allocatable,dimension(:,:,:) :: VsB, VpB
real(rsp),pointer,dimension(:,:,:)     :: V
target                                 :: VsB, VpB

! timing stuff
real(rdp) :: tic, tac

!-----------------------------------------------------------------------------------------

! depth of actual Z-level
zo = (z_level(rank) - 1) * output_dx

! 3D deterministic velocity model provided as binary file
if(backgr_mode .eq. 'b') then

   ! find first/last index of the background Y-axis encompassing the current
   ! Y-axis chunk
   j_min = floor( (y_level(rank,1)-1) * output_dx / backgr_dx(2) ) + 1
   j_max = ceiling( (y_level(rank,2)-1) * output_dx / backgr_dx(2) ) + 1

   ! ...corresponding number of points
   j_npts = j_max - j_min + 1

   ! find first/last index of the background Z-axis encompassing the current Z-level
   k_min = floor( zo / backgr_dx(3) ) + 1
   k_max = k_min + 1

   ! handle the case when we are at end of the grid (and both input and output grids have
   ! identical step)
   if (k_max .gt. backgr_N(3)) then
      k_min = k_min - 1
      k_max = k_max - 1
   endif

   ! ...corresponding vector
   z = (/ (k_min - 1) * backgr_dx(3), (k_max - 1) * backgr_dx(3) /)

   ! find first/last index of the background X-axis encompassing output X-axis
   i_min = 1
   i_max = ceiling( (output_N(1)-1) * output_dx / backgr_dx(1) ) + 1

   ! ...corresponding number of points
   i_npts = i_max - i_min + 1

   ! number of points along XY
   ij_npts = i_npts * j_npts

   ! total points
   ijk_npts = ij_npts * 2

   ! Here below each MPI task (belonging to WRT_COMM) reads part of the background velocity
   ! model.

   ! allocate array where Vs values of background model are stored
   allocate(VsB(i_npts,j_npts,2))

   tic = MPI_Wtime()

   ! now we open and read Vs values for the background velocity model.
   call MPI_FILE_OPEN(WRT_COMM,trim(backgr_Vs),MPI_MODE_RDONLY,MPI_INFO_NULL,fileq,err)

   ! define dimension and starting points of subarray
   subsizes = (/i_npts,j_npts,2/); starts = (/i_min - 1,j_min - 1,k_min - 1/)

   call MPI_TYPE_CREATE_SUBARRAY(3,backgr_N,subsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL,filetype,err)
   call MPI_TYPE_COMMIT(filetype,err)
   disp = 0
   call MPI_FILE_SET_VIEW(fileq,disp,MPI_REAL,filetype,'native',MPI_INFO_NULL,err)

   call MPI_FILE_READ(fileq,VsB,ijk_npts,MPI_REAL,status,err)

   call MPI_FILE_CLOSE(fileq,err)

   tac = MPI_Wtime(); tac = tac - tic

   if ( (rank .eq. 0) .and. lprint(5) ) then
      print*,'{Background model read in ',real(tac),' sec}'
      print*,''
      lprint(5) = .false.
   endif

   ! do the same for Vp
   if (trim(backgr_Vp) .ne. 'none') then

      allocate(VpB(i_npts,j_npts,2))

      call MPI_FILE_OPEN(WRT_COMM,trim(backgr_Vp),MPI_MODE_RDONLY,MPI_INFO_NULL,fileq,err)

      call MPI_TYPE_CREATE_SUBARRAY(3,backgr_N,subsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL,  &
                                    filetype,err)
      call MPI_TYPE_COMMIT(filetype,err)
      disp = 0
      call MPI_FILE_SET_VIEW(fileq,disp,MPI_REAL,filetype,'native',MPI_INFO_NULL,err)

      call MPI_FILE_READ(fileq,VpB,ijk_npts,MPI_REAL,status,err)

      call MPI_FILE_CLOSE(fileq,err)

   endif

   ! Now interpolate Vs/Vp values onto current grid

   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,i_min,i_max,xo,x,j_min,j_max,yo,y,V,Vs,Vp) &
   !$OMP REDUCTION(MAX:compos_Vs_max,compos_Vp_max) REDUCTION(MIN:compos_Vs_min,compos_Vp_min)
   do i = 1,output_N(1)

      ! position along X of current point
      xo = (i - 1) * output_dx

      ! corresponding indexes of background velocity grid
      i_min = floor( xo / backgr_dx(1) ) + 1

      ! handle the case when we are at end of the grid (and both input and output grids
      ! have identical step)
      if (i_min .eq. size(VsB,1)) i_min = i_min - 1

      i_max = i_min + 1

      x = (/(i_min - 1) * backgr_dx(1), (i_max - 1) * backgr_dx(1)/)

      !do j = 1,npts_ychunk(c_block)
      do j = 1,mypoints

         ! position along Y of current point
         yo = (j - 1) * output_dx

         ! corresponding indexes of background velocity grid
         j_min = floor( yo / backgr_dx(2) ) + 1

         ! handle the case when we are at end of the grid (and both input and output grids
         ! have identical step)
         if (j_min .eq. size(VsB,2)) j_min = j_min - 1

         j_max = j_min + 1

         y = (/(j_min - 1) * backgr_dx(2), (j_max - 1) * backgr_dx(2)/)

         ! alias to part of background model
         V => VsB(i_min:i_max,j_min:j_max,:)

         call poly_interp(x,y,z,V,xo,yo,zo,Vs)

         if (backgr_Vp .ne. 'none') then

            V => VpB(i_min:i_max,j_min:j_max,:)

            call poly_interp(x,y,z,V,xo,yo,zo,Vp)

         endif

         ! Add random field to background velocity. The random field must be considered as
         ! perturbation in % (we express input sigma in %)
         Vs = Vs * (1. + real(SPEC_Z(j,i))/100.)

         Vp = Vp * (1. + real(SPEC_Z(j,i))/100.)

         if (trim(backgr_Vp) .eq. 'none') Vp = 0.

         ! store results into existing array 'SPEC_Z'
         SPEC_Z(j,i) = cmplx(Vs,Vp)

         ! store min/max values
         compos_Vs_min = min(compos_Vs_min,Vs)
         compos_Vs_max = max(compos_Vs_max,Vs)

         compos_Vp_min = min(compos_Vp_min,Vp)
         compos_Vp_max = max(compos_Vp_max,Vp)

         nullify(V)

      enddo
   enddo
   !$OMP END PARALLEL DO

   ! deallocate arrays
   deallocate(VsB)
   if (trim(backgr_Vp) .ne. 'none') deallocate(VpB)

endif

! 1D deterministic velocity model provided as text file. In this case we assume the velocity
! model has been already loaded in memory
if(backgr_mode .eq. 't') then

   ! make a copy of 1D z-vector
   tmp = depth1D

   ! find layers containing current point
   k_min = minloc(abs(depth1D - zo),dim = 1)

   if (zo .lt. depth1D(k_min)) tmp(k_min:n_layers) = -huge(1.)
   if (zo .gt. depth1D(k_min)) tmp(1:k_min)        = -huge(1.)
   if (zo .eq. depth1D(k_min)) tmp(k_min)          = -huge(1.)

   k_max = minloc(abs(tmp - zo),dim = 1)

   ! make sure k_min < k_max
   if (k_min .gt. k_max) then
      i = k_max
      k_max = k_min
      k_min = i
   endif

   z = (/depth1D(k_min),depth1D(k_max)/)

   vsl = (/vs1D(k_min),vs1D(k_max)/)
   vpl = (/vp1D(k_min),vp1D(k_max)/)

   ! Linear interpolation is used to introduce velocity gradients in 1D models
   call polint(z,vsl,zo,vsi)
   call polint(z,vpl,zo,vpi)

   !do j = 1,npts_ychunk(c_block)
   do j = 1,mypoints
      do i = 1,output_N(1)

         ! Add random field to background velocity. The random field must be considered as
         ! perturbation in % (we express input sigma in %)
         Vs = vsi * (1. + real(SPEC_Z(j,i))/100.)

         Vp = vpi * (1. + real(SPEC_Z(j,i))/100.)

         ! store results into existing array 'SPEC_Z'
         SPEC_Z(j,i) = cmplx(Vs,Vp)

         ! store min/max values
         compos_Vs_min = min(compos_Vs_min,Vs)
         compos_Vs_max = max(compos_Vs_max,Vs)

         !print*,compos_Vs_min,vsi,Vs,i,j,l,u

         compos_Vp_min = min(compos_Vp_min,Vp)
         compos_Vp_max = max(compos_Vp_max,Vp)

      enddo
   enddo

endif

END SUBROUTINE resample_deterministic_model

!=========================================================================================

SUBROUTINE write_to_disk
!
! Description:
!
!   Write result contained in array SPEC_Z to disk. This subroutine invokes non-blocking
!   IO command. Output is forced to complete in the main program.
!
! Author: W. Imperatori (walter.imperatori@sed.ethz.ch)
!
! Version: Created on August 2011, v1.0
!          Modified on December 2017, v1.3 - Some variables were renamed to reflect changes
!                                            in the main program
!

use precisions

implicit none

! local variables
integer(kind=mpi_offset_kind) :: offset
integer(isp)                  :: N

!-----------------------------------------------------------------------------------------

! number of points to be written to disk
!N = output_N(1) * npts_ychunk(c_block)
N = output_N(1) * mypoints

! first, write VP values to disk if available. In this case all tasks belonging to the
! communicator write synchronously to disk
if ( (trim(backgr_Vp) .ne. 'none') .and. (backgr_mode .ne. 'n') ) then

   call MPI_FILE_OPEN(WRT_COMM,trim(output_Vp),MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL, &
&                     filep,err)

   ! Always write at the end of the file (i.e. append)
   offset = 0
   call MPI_FILE_SEEK_SHARED(filep,offset,MPI_SEEK_END,err)

   ! copy values into output array
   IOar = transpose(aimag(SPEC_Z(:,1:output_N(1))))

   call MPI_FILE_WRITE_ORDERED(filep,IOar,N,MPI_REAL,status,err)

   call MPI_FILE_CLOSE(filep,err)

endif

! then, write Vs values
call MPI_FILE_OPEN(WRT_COMM,trim(output_Vs),MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL, &
&                  filep,err)

! Always write at the end of the file (i.e. append)
offset = 0
call MPI_FILE_SEEK_SHARED(filep,offset,MPI_SEEK_END,err)

! copy values into output array
IOar = transpose(real(SPEC_Z(:,1:output_N(1))))

call MPI_FILE_WRITE_ORDERED_BEGIN(filep,IOar,N,MPI_REAL,err)

! update status flag for writing
lwrite = .true.

END SUBROUTINE write_to_disk

!=========================================================================================

END PROGRAM random_3d_RAM_MPI

!=========================================================================================

SUBROUTINE online_variance(val, var, mu, n)
!
! Description:
!
!    Algorithm of Welford to compute online variance and average. Output var must be divided
!    by n-1 to obtain variance
!
! Author: W. Imperatori (walter.imperatori@sed.ethz.ch)
!
! Version: Created on August 2011, v1.0
!          Modified on December 2017, v1.3 - Some variables were renamed to reflect changes
!                                            in the main program. Switched to double precision
!                                            to avoid inaccurate calculations
!

use precisions

implicit none

real(rsp),intent(in)    :: val
real(rdp),intent(inout) :: mu, var, n

! local variables
real(rdp)    :: delta, delta2

!-----------------------------------------------------------------------------------------

! update number of points
n = n + 1

delta = val - mu

mu = mu + delta / n

delta2 = val - mu

var = var + delta * delta2

END SUBROUTINE online_variance

!=========================================================================================

SUBROUTINE parallel_variance(var1,mu1,n1,var2,mu2,n2)
!
! Description:
!
!    Algorithm of Chan et al. to compute variance and average of a partitioned dataset
!
! Author: W. Imperatori (walter.imperatori@sed.ethz.ch)
!
! Version: Created on August 2011, v1.0
!          Modified on December 2017, v1.3 - Some variables were renamed to reflect changes
!                                            in the main program. Switched to double precision
!                                            to avoid inaccurate calculations
!

use precisions

implicit none

real(rdp),intent(in)    :: var1, mu1, n1
real(rdp),intent(inout) :: var2, mu2, n2

! local variables
real(rdp) :: m_1, m_2, M2, delta

!-----------------------------------------------------------------------------------------

delta = mu2 - mu1

m_1 = var1 * (n1 - 1.)

m_2 = var2 * (n2 - 1.)

M2 = m_1 + m_2 + delta**2 * n1 * n2 / (n1 + n2)

M2 = M2 / (n1 + n2 - 1.)

! assign output mean
mu2 = mu1 + delta * n2 / (n1 + n2)

! assign output variance
var2 = M2

! assign output points
n2 = n1 + n2

END SUBROUTINE parallel_variance
