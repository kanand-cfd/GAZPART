!!====================================================================
!!
!!     3-dimensional turbulent spectrum and statistics
!!
!!====================================================================

subroutine SPEC3D(NOUT)

!!====================================================================
!!
!!
!!====================================================================

use dns_dim            !- Dimension
use param_phys         !- Physical & numerical parameters
use fluid_variable     !- Fluid velocity
use geometric_variable !- Mesh & Wavenumbers
use forcing, only:TIME_FORCE, SIGMA_FORCE
use statistics, only:EPS_FLU_SPEC, MEAN_DT, MAX_DT, MIN_DT


implicit none


!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Global arrays
!---------------------------------------------------------------------

!- 
integer, intent(in) :: NOUT

!---------------------------------------------------------------------
!- Local arrays
!---------------------------------------------------------------------
!- 3d spectrum
real(kind=8), dimension(:), allocatable :: SPECTRUM, SPECLOC

!- File name
character (len=15) :: FILESPEC

!- 
real(kind=8) :: KXGLOB, KAPPA2, KAPPA

!- Kinetic energy
real(kind=8) :: KF_SPEC

!- Dissipation
real(kind=8) :: EPS_SPEC

!- Eddy Lifetime
real(kind=8) :: TE_SPEC

!- Integral length scale
real(kind=8) :: LF_SPEC

!- Taylor length scale
real(kind=8) :: LAMBDA_SPEC

!- Reynolds number
real(kind=8) :: RE_SPEC

!- Kolmogorov length scale
real(kind=8) :: ETA_KOL

!- Kolmogorov Time scale
real(kind=8) :: TIME_KOL

!- Time step CFL
real(kind=8) :: DTIME_CFL

!-Forcing
real(kind=8) :: EPS_STAR
real(kind=8) :: RE_STAR
real(kind=8) :: TL_STAR
real(kind=8) :: T0_STAR



!- Index
integer :: I, J, K, N, KMAX, IK
!------------------------------------------------------------------


!!====================================================================
!! 1. initiation
!!====================================================================
KMAX = NX/2 + 1

allocate(SPECTRUM(KMAX))
allocate(SPECLOC(KMAX))


SPECTRUM(:) = 0.
SPECLOC(:) = 0.


!!====================================================================
!! 2. 3d spectrum
!!====================================================================
do K = FSTART(3),FEND(3)
 do J = FSTART(2),FEND(2)
  do I = FSTART(1),FEND(1)

    KAPPA2 = (KX(I)/KFIRST)**2 + (KY(J)/KFIRST)**2 + (KZ(K)/KFIRST)**2

    IK = sqrt(KAPPA2) + 1

    if(IK>KMAX) IK = NX-IK

    SPECLOC(IK) = SPECLOC(IK)                          &
                + real(UFOU(I,J,K)*conjg(UFOU(I,J,K))) &
                + real(VFOU(I,J,K)*conjg(VFOU(I,J,K))) &
                + real(WFOU(I,J,K)*conjg(WFOU(I,J,K))) 

  enddo
 enddo
enddo


call MPI_ALLREDUCE(SPECLOC,SPECTRUM, &
              KMAX,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

!!--------------------------------------------------------------------
!! 1.2. Normalisation and printing 
!!--------------------------------------------------------------------
do N = 1, KMAX
 SPECTRUM(N) = SPECTRUM(N)/KFIRST
end do


!!====================================================================
!! 3. Statistics
!!====================================================================

!!--------------------------------------------------------------------
!! 3.0. Spectrum statistics
!!--------------------------------------------------------------------
if(MYID==0) then

 KF_SPEC = 0.
EPS_SPEC = 0.
 LF_SPEC = 0.

do N = 2, KMAX
   KXGLOB = (N-1)*KFIRST
  KF_SPEC =  KF_SPEC +                   SPECTRUM(N)       *KFIRST
 EPS_SPEC = EPS_SPEC + 2.*VISC*KXGLOB**2*SPECTRUM(N)       *KFIRST
  LF_SPEC =  LF_SPEC +                   SPECTRUM(N)/KXGLOB*KFIRST
end do

!!- Save dissipation for spec3d_scalar.f90
EPS_FLU_SPEC = EPS_SPEC

ETA_KOL = (VISC**3/EPS_SPEC)**0.25
TIME_KOL = (VISC/EPS_SPEC)**0.5

LF_SPEC = LF_SPEC*PPI/(4.*KF_SPEC/3.)
RE_SPEC = LF_SPEC*sqrt(2.*KF_SPEC/3.)/VISC
TE_SPEC = LF_SPEC/sqrt(2.*KF_SPEC/3.)

!!--------------------------------------------------------------------
!! 2.1. Print in file
!!--------------------------------------------------------------------
if(NOUT<0) then
 !- Case for interpolation checking
 FILESPEC = 'spec_int.dat'
elseif(NOUT<10) then
 write(FILESPEC,4001)NOUT
elseif(NOUT<100) then
 write(FILESPEC,4002)NOUT
else
 write(FILESPEC,4003)NOUT
end if
open(unit=200, file=trim(FILESPEC))


 write(200,*)'# k, k*eta_K, E(k), &
&E(k)/(eps*nu^5)^0.25, k*E(k)/kf,2*nu*k^3*E(k)/eps'

do N=2, KMAX-1
 KXGLOB = (N-1)*KFIRST
 write (200,10000)KXGLOB, KXGLOB*ETA_KOL,                          &
                           SPECTRUM(N),                          &
                           SPECTRUM(N)/(EPS_SPEC*VISC**5)**0.25, &
                           KXGLOB*SPECTRUM(N)/KF_SPEC,            &
                           2.*VISC*KXGLOB**3.*SPECTRUM(N)/EPS_SPEC
end do
close (200)


!!- Time step for CFL restriction
!!  Restriction from Pope's book
DTIME_CFL = 1./20.*DX/KF_SPEC**0.5


deallocate(SPECTRUM)
deallocate(SPECLOC)


EPS_STAR = SIGMA_FORCE**2.0*TIME_FORCE
RE_STAR = EPS_STAR**(1./3.)*KFIRST**(-4./3.)/VISC
TL_STAR = TIME_FORCE*EPS_STAR**0.5*KFIRST**(2./3.)
T0_STAR = TE_SPEC*(EPS_STAR*KFIRST**2)**(2./3.)



!!--------------------------------------------------------------------
!! 2.2. Print in "stat.info" the last spectrum statistics
!!--------------------------------------------------------------------
if((NOUT == FOUT1).or.(NOUT == 99)) then
 I=1
 write(UNIT_INFO(I),*)
 write(UNIT_INFO(I),*)'====================================================================='
 write(UNIT_INFO(I),*)'FLUID STATISTICS FROM SPECTRUM'
 write(UNIT_INFO(I),*)'====================================================================='
 write(UNIT_INFO(I),*)
 write(UNIT_INFO(I),10603)'Kinetic energy:               kf = ',KF_SPEC,' m2/s2'
 write(UNIT_INFO(I),10603)'Longitudinal length scale:    Lf = ',LF_SPEC,' m'
 write(UNIT_INFO(I),10604)'                           Lx/Lf = ',LXMAX/LF_SPEC
 write(UNIT_INFO(I),10603)'Eddy turnover:          Te=Lf/up = ',TE_SPEC,' s'
 write(UNIT_INFO(I),10604)'Reynolds number:             ReL = ',RE_SPEC
 write(UNIT_INFO(I),*)
 write(UNIT_INFO(I),10603)'Dissipation:                 eps = ',EPS_SPEC,' m2/s3'
 write(UNIT_INFO(I),10603)'Kolmogorov length scale:     eta = ',ETA_KOL,' m'
 write(UNIT_INFO(I),10603)'Kolmogorov   time scale:   tau_k = ',TIME_KOL,' s'
 write(UNIT_INFO(I),10604)'               eta.(kmax*Rtrunc) = ',ETA_KOL*KLAST*RTRUNC
 write(UNIT_INFO(I),*)
 write(UNIT_INFO(I),10603)' Time step for optimal CFL stability = ',DTIME_CFL,' s'
 write(UNIT_INFO(I),*)
 if(STEADY) then
 write(UNIT_INFO(I),*)' Check forcing parameters'
 write(UNIT_INFO(I),*)' ------------------------'
 write(UNIT_INFO(I),10604)'     Te/Tf =',TE_SPEC/TIME_FORCE
 write(UNIT_INFO(I),10604)'     <dt>/Tf =',MEAN_DT/TIME_FORCE
 write(UNIT_INFO(I),10604)'  min(dt)/Tf =',MIN_DT/TIME_FORCE
 write(UNIT_INFO(I),10604)'  max(dt)/Tf =',MAX_DT/TIME_FORCE
 write(UNIT_INFO(I),10604)'  Tf/tau_k =',TIME_FORCE/TIME_KOL
 write(UNIT_INFO(I),10603)' epsilon_* = ',EPS_STAR,' m2/s'
 write(UNIT_INFO(I),10602)'      Re_* = ',RE_STAR
 write(UNIT_INFO(I),10602)'      Tf_* = ',TL_STAR
 write(UNIT_INFO(I),10602)'      T0_* = ',T0_STAR
 end if
end if

end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
4001 format ('spec_0',i1,'.stat')
4002 format ('spec_',i2,'.stat')
4003 format ('spec_',i3,'.stat')

10000 format (10(e17.7))

10601 format (2x,A,E13.6,A,E13.6)
10602 format (2x,A,E13.6)
10603 format (2x,A,E13.6,A)
10604 format (2x,A,f8.4)

end subroutine SPEC3D
