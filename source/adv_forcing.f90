!!====================================================================
!!
!!  Forcing scheme for statistically steady turbulent flows
!!
!!====================================================================

subroutine ADV_FORCING

!!====================================================================
!!
!! The forcing scheme, proposed by Eswaran & Pope, J. Comp. & Fluids
!! (1988) is based on a stochastic force added at low-wavenumbers.
!!
!! The stochastic force is given by a Ornstein-Uhlenbeck process,
!! parameterized by a timescale (TIME_FORCE) and a variancce
!! (SIGMA_FORCE). The range of modified wavenumber is controled by
!! KFORCING_MIN and KFORCING_MAX.
!!
!! 
!!   KFORCING_MIN: minimum forced wavenumber
!!   KFORCING_MAX: maximum forced wavenumber
!!   TIME_FORCE: timescale of forcing
!!   SIGMA_FORCE: variance of forcing
!!
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE
use FLUID_VARIABLE
use PARAM_PHYS
use FORCING
use CHECK_CPU
use RANDOM_VAR

use P3DFFT

implicit none

!!--------------------------------------------------------------------
!! ARRAYS STATEMENT
!!--------------------------------------------------------------------

!!--------------------------------------------------------------------
real(kind=8) :: XAMPF
real(kind=8) :: ALPHA, BETA

!- Random number generator
real(kind=8) :: GAUSS

real(kind=8) :: XAMPL, XE

!- number of forced waves
integer :: NF

!- Wavenumber
real(kind=8) :: KXFULL, KYFULL, KZFULL, KF

!- Imaginary and real part of forcing
real(kind=8) :: FORCE_REAL, FORCE_IMAG

!- Divergence
double complex :: DIVC


!- Statistics of the droped field
real(kind=8) :: STATX, STATY, STATZ

!- index
integer :: I, J, K, MJ, MK

!---------------------------------------------------------------------

!!====================================================================
!! 1. Compute stochastic force
!!====================================================================

if(MYID==0) then

 ALPHA =                     DTIME / TIME_FORCE
 BETA  = (2.0*SIGMA_FORCE**2*DTIME / TIME_FORCE)**0.5

 NF = 0
 do K = 1, NZ
  do J = 1, NY
   do I = 1, NX/2+1

!!- Wavenumbers
    KXFULL = (I-1)*KFIRST
    if(J<=NY/2+1) then
     KYFULL = (J-1)*KFIRST
    else
     KYFULL = (J-1-NY)*KFIRST
    end if
    if(K<=NZ/2+1) then
     KZFULL = (K-1)*KFIRST
    else
     KZFULL = (K-1-NZ)*KFIRST
    end if

!! KF = (KX(I)**2 + KY_FULL(J)**2 + KZ_FULL(K)**2)**0.5
    KF = (KXFULL**2 + KYFULL**2 + KZFULL**2)**0.5

    if ((KF>=KFORCE_MIN).and.(KF<=KFORCE_MAX)) then

     NF = NF + 1

     XAMPF = (KFORCE_MAX - KF) / 0.20 / KFORCE_MAX
     XAMPL = (exp(2.0*XAMPF) - 1.0) / (exp(2.0*XAMPF) + 1.0)

!!-------------------------------------------------------------------- 
!! 1.1. x-velocity
!!-------------------------------------------------------------------- 
!- real part
     XE = GAUSS(IDFORCE)
!  XE = GAUSS2(XRANDF)
 
     XE = XE * XAMPL
     FORCE_REAL = real(FORCING_UFOU(NF))*(1.0-ALPHA) + XE*BETA

!- complex part
     XE = GAUSS(IDFORCE)
!  XE = GAUSS2(XRANDF)

     XE = XE * XAMPL
     FORCE_IMAG = aimag(FORCING_UFOU(NF))*(1.0-ALPHA) + XE*BETA
 
!- Stochastic acceleration
     FORCING_UFOU(NF) = cmplx(FORCE_REAL,FORCE_IMAG) 


!!-------------------------------------------------------------------- 
!! 1.2. y-velocity
!!-------------------------------------------------------------------- 
!- real part
     XE = GAUSS(IDFORCE)
!  XE = GAUSS2(YRANDF)

     XE = XE * XAMPL
     FORCE_REAL = real(FORCING_VFOU(NF))*(1.0-ALPHA) + XE*BETA

!- complex part
     XE = GAUSS(IDFORCE)
!  XE = GAUSS2(YRANDF)

     XE = XE * XAMPL
     FORCE_IMAG = aimag(FORCING_VFOU(NF))*(1.0-ALPHA) + XE*BETA 

!- Stochastic acceleration
     FORCING_VFOU(NF) = cmplx(FORCE_REAL,FORCE_IMAG) 

!!-------------------------------------------------------------------- 
!! 1.3. z-velocity
!!-------------------------------------------------------------------- 

!- real part
     XE = GAUSS(IDFORCE)
!  XE = GAUSS2(ZRANDF)

     XE = XE * XAMPL
     FORCE_REAL = real(FORCING_WFOU(NF))*(1.0-ALPHA) + XE*BETA

!- complex part
     XE = GAUSS(IDFORCE)
!  XE = GAUSS2(ZRANDF)

     XE = XE * XAMPL
     FORCE_IMAG = aimag(FORCING_WFOU(NF))*(1.0-ALPHA) + XE*BETA

 
!- Stochastic acceleration
     FORCING_WFOU(NF) = cmplx(FORCE_REAL,FORCE_IMAG) 

    end if

   end do
  end do
 end do 

end if !!- End if MYID==0


!!- Broadcast to all task
 call MPI_BCAST(FORCING_UFOU,NFORCE_FULL,MPI_COMPLEX16,0,MPI_COMM_WORLD,IERR)
 call MPI_BCAST(FORCING_VFOU,NFORCE_FULL,MPI_COMPLEX16,0,MPI_COMM_WORLD,IERR)
 call MPI_BCAST(FORCING_WFOU,NFORCE_FULL,MPI_COMPLEX16,0,MPI_COMM_WORLD,IERR)



!STATX = ZERO
!STATY = ZERO
!STATZ = ZERO
!do I=1,NFORCE_FULL
! STATX = STATX + real(FORCING_UFOU(I)*conjg(FORCING_UFOU(I)))
! STATY = STATY + real(FORCING_VFOU(I)*conjg(FORCING_VFOU(I)))
! STATZ = STATZ + real(FORCING_VFOU(I)*conjg(FORCING_WFOU(I)))
!end do
!
!call RSUMCPU(STATX,STATX)
!call RSUMCPU(STATY,STATY)
!call RSUMCPU(STATZ,STATZ)
!
!if(MYID==0) write(*,*) '       <fx> = ',STATX / NFORCE_FULL
!if(MYID==0) write(*,*) '       <fy> = ',STATY / NFORCE_FULL
!if(MYID==0) write(*,*) '       <fz> = ',STATZ / NFORCE_FULL



end subroutine ADV_FORCING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function Gauss(ID)
!
! use RANDOM_VAR
!
!  implicit none
!
!  integer, intent(inout) :: ID
!
!  real(kind=8) :: Gauss
!  real(kind=8) :: FAC,RSQ,V1,V2
!  real(kind=8) :: XRAND
!  real(kind=8), external :: RAN1
!
!  intrinsic LOG, SQRT
!
!  if (ISET /= 0) then
!    Gauss = GSET
!    ISET = 0
!    return
!  end if
!  do
!    ! Generateur aleatoire tire de Numerical receipes (routine interne)
!    call random_number(XRAND)
!    V1 = 2.*XRAND - 1.
!    call random_number(XRAND)
!    V2 = 2.*XRAND - 1.
!
!!    V1 = 2.*RAN1(ID) - 1.
!!    V2 = 2.*RAN1(ID) - 1.
!    !
!    RSQ = V1**2 + V2**2
!    if (RSQ<1. .and. RSQ/=0.) exit
!  end do
!  FAC = SQRT(-2.*LOG(RSQ)/RSQ)
!  GSET = V1 * FAC
!  Gauss = V2 * FAC
!  ISET = 1
!end function Gauss
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function Gauss2(XRAND)
!
! use RANDOM_VAR
!
!  implicit none
!
!!!  integer, intent(inout) :: ID
!
!  real(kind=8) :: Gauss2
!  real(kind=8) :: FAC,RSQ,V1,V2
!  real(kind=8),intent(inout) :: XRAND
!
!!!  real(kind=8), external :: RAN1
!
!  intrinsic LOG, SQRT
!
!  if (ISET /= 0) then
!    Gauss2 = GSET
!    ISET = 0
!    return
!  end if
!  do
!    ! Generateur aleatoire tire de Numerical receipes (routine interne)
!    call random_number(XRAND)
!    V1 = 2.*XRAND - 1.
!
!    call random_number(XRAND)
!    V2 = 2.*XRAND - 1.
!
!!    V1 = 2.*RAN1(ID) - 1.
!!    V2 = 2.*RAN1(ID) - 1.
!    !
!    RSQ = V1**2 + V2**2
!    if (RSQ<1. .and. RSQ/=0.) exit
!  end do
!  FAC = SQRT(-2.*LOG(RSQ)/RSQ)
!  GSET = V1 * FAC
!  Gauss2 = V2 * FAC
!  ISET = 1
!end function Gauss2
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function RAN1(IDUM)

  use RANDOM_VAR

  implicit none

  integer,      parameter :: IA = 16807
  integer,      parameter :: IM = 2147483647
  real(kind=8), parameter :: AM = 1./IM
  integer,      parameter :: IQ = 127773
  integer,      parameter :: IR = 2836
  integer,      parameter :: NTAB = 32
  integer,      parameter :: NDIV = 1+(IM-1)/NTAB
  real(kind=8), parameter :: EPS = 1.2e-7
  real(kind=8), parameter :: RNMX = 1.-EPS

  integer, intent(inout) :: IDUM

  real(kind=8) :: RAN1

  integer :: J,K

!  intrinsic MAX, MIN

!  common /GENE2/ IV,IY
!  integer, dimension(32) :: IV
!  integer :: IY = 0

!!  data IV/NTAB*0/


  if (IDUM<=0 .or. IY==0) then
    IDUM = MAX(-IDUM,1)
    do J = NTAB+8,1,-1
      K = IDUM / IQ
      IDUM = IA*(IDUM-K*IQ) - IR*K
      if (IDUM < 0) then
        IDUM = IDUM + IM
      end if
      if (J <= NTAB) then
        IV(J) = IDUM
      end if
    end do
    IY = IV(1)
  end if
  K = IDUM / IQ
  IDUM = IA*(IDUM-K*IQ) - IR*K
  if (IDUM < 0) then
    IDUM = IDUM + IM
  end if
  J = 1 + IY/NDIV
  IY = IV(J)
  IV(J) = IDUM
  ran1 = MIN(AM*IY,RNMX)

end function RAN1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Gauss(ID)

 use RANDOM_VAR

 implicit none

 integer, intent(inout) :: ID

 real(kind=8) :: Gauss
 real(kind=8) :: FAC,RSQ,V1,V2

 real(kind=8), external :: RAN1

! intrinsic LOG, SQRT

! common /GENE1/ ISET,GSET

! integer :: ISET = 0
! real(kind=8) :: GSET

  if (ISET /= 0) then
    Gauss = GSET
    ISET = 0
    return
  end if
  do
    ! Generateur aleatoire tire de Numerical receipes (routine interne)
    V1 = 2.*RAN1(ID) - 1.
    V2 = 2.*RAN1(ID) - 1.
    !
    RSQ = V1**2 + V2**2
    if (RSQ<1. .and. RSQ/=0.) exit
  end do
  FAC = SQRT(-2.*LOG(RSQ)/RSQ)
  GSET = V1 * FAC
  Gauss = V2 * FAC
  ISET = 1
end function Gauss

!
