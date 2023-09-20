!!====================================================================
!!
!! Direct Numerical Simulation of turbulent fluid flow
!!
!!====================================================================

subroutine ADV_FLUID(NCYCLE)

!!====================================================================
!!
!!
!!====================================================================

use DNS_DIM             !- Dimension
use PARAM_PHYS          !- Physical & numerical parameters
use FORCING             !- Forcing
use FLUID_VARIABLE      !- Fluid velocity
use SCALAR_VARIABLE    
use GEOMETRIC_VARIABLE 
use STATISTICS         
use WORK_ARRAYS
use CHECK_CPU       

use MPI_STRUCTURES

use P3DFFT

implicit none

!!====================================================================
!! ARRAYS STATEMENT
!!--------------------------------------------------------------------
!- cycle number
integer, intent(in) :: NCYCLE

double complex, dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)) :: RHS_UFOU
double complex, dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)) :: RHS_VFOU
double complex, dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)) :: RHS_WFOU

double complex, dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)) :: UFOU_NP1
double complex, dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)) :: VFOU_NP1
double complex, dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)) :: WFOU_NP1

!- Runge-Kutta Step
integer :: RKSTEP
!- squared wavenumber
real(kind=8) :: KAPPA2


!- Time control variable
real(kind=8) :: TIME_START, TIME_END
real(kind=8) :: TIME_START2

!!real(kind=8) :: SHIFT
double complex :: SHIFT

integer :: I, J, K
!!
!!====================================================================

!!- CPU check
if(MYID == 0) then
 TIME_START = MPI_WTIME()
 TIME_START2 = TIME_START
end if




!!====================================================================
!! 1. Transform the variable from physical space to Fourier space
!!====================================================================

!!--------------------------------------------------------------------
!! 1.1. Fluid velocity from physical space to Fourier space
!!--------------------------------------------------------------------
TMPPHY(:,:,:) = UFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))
call P3DFFT_FTRAN_R2C(TMPPHY,UFOU,FFTFLAG) !- x-component of fluid velocity

TMPPHY(:,:,:) = VFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))
call P3DFFT_FTRAN_R2C(TMPPHY,VFOU,FFTFLAG) !- y-component of fluid velocity

TMPPHY(:,:,:) = WFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))
call P3DFFT_FTRAN_R2C(TMPPHY,WFOU,FFTFLAG) !- z-component of fluid velocity

!- FFT Normalization
UFOU(:,:,:) = UFOU(:,:,:)*FACTOR
VFOU(:,:,:) = VFOU(:,:,:)*FACTOR
WFOU(:,:,:) = WFOU(:,:,:)*FACTOR


if(MYID == 0) then
 TIME_END = MPI_WTIME()
 CPU_FLUID(8) = CPU_FLUID(8) +  TIME_END - TIME_START
end if


!!====================================================================
!! 1. Runge-Kutta integration scheme
!!====================================================================
!!- Advance the forcing scheme
if(STEADY) call ADV_FORCING


!!--------------------------------------------------------------------
!! 1.1. First step
!!--------------------------------------------------------------------
RKSTEP = 1

RHS_UFOU(:,:,:) = cmplx(ZERO,ZERO)
RHS_VFOU(:,:,:) = cmplx(ZERO,ZERO)
RHS_WFOU(:,:,:) = cmplx(ZERO,ZERO)

!!- add forcing
if(STEADY) call ADD_FORCING3(RHS_UFOU,RHS_VFOU,RHS_WFOU)

!!- add viscous terms
! call VISCOUS_TERMS(    UFOU,    VFOU,    WFOU, &
!                    RHS_UFOU,RHS_VFOU,RHS_WFOU)

!!- add non-linear terms
 call NONLINEAR_TERMS(    UFLU,    VFLU,    WFLU, &
                          UFOU,    VFOU,    WFOU, &
                      RHS_UFOU,RHS_VFOU,RHS_WFOU, RKSTEP )




do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)

   KAPPA2 = KX(I)**2 + KY(J)**2 + KZ(K)**2

!- Store Un+1 = Un + dt/2*k1
   UFOU_NP1(I,J,K) = (UFOU(I,J,K) + DTIME/2*RHS_UFOU(I,J,K))*exp(-VISC*KAPPA2*DTIME)
   VFOU_NP1(I,J,K) = (VFOU(I,J,K) + DTIME/2*RHS_VFOU(I,J,K))*exp(-VISC*KAPPA2*DTIME)
   WFOU_NP1(I,J,K) = (WFOU(I,J,K) + DTIME/2*RHS_WFOU(I,J,K))*exp(-VISC*KAPPA2*DTIME)



   !- Aproximate fluid velocity as un+1
   UFOU(I,J,K) = (UFOU(I,J,K) + DTIME * RHS_UFOU(I,J,K))*exp(-VISC*KAPPA2*DTIME)
   VFOU(I,J,K) = (VFOU(I,J,K) + DTIME * RHS_VFOU(I,J,K))*exp(-VISC*KAPPA2*DTIME)
   WFOU(I,J,K) = (WFOU(I,J,K) + DTIME * RHS_WFOU(I,J,K))*exp(-VISC*KAPPA2*DTIME)
  end do
 end do
end do


!!--------------------------------------------------------------------
!! 2.2. Second step
!!--------------------------------------------------------------------
RKSTEP = 2



!- Back in physical space
TMPFOU = UFOU
call P3DFFT_BTRAN_C2R(TMPFOU,UFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)       
TMPFOU = VFOU 
call P3DFFT_BTRAN_C2R(TMPFOU,VFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)
TMPFOU = WFOU
call P3DFFT_BTRAN_C2R(TMPFOU,WFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)       



RHS_UFOU(:,:,:) = cmplx(ZERO,ZERO)
RHS_VFOU(:,:,:) = cmplx(ZERO,ZERO)
RHS_WFOU(:,:,:) = cmplx(ZERO,ZERO)

!!- add forcing
if(STEADY) call ADD_FORCING3(RHS_UFOU,RHS_VFOU,RHS_WFOU)

!!- add viscous terms
! call VISCOUS_TERMS(    UFOU,    VFOU,    WFOU, &
!                    RHS_UFOU,RHS_VFOU,RHS_WFOU)

!!- add non-linear terms
 call NONLINEAR_TERMS(    UFLU,    VFLU,    WFLU, &
                          UFOU,    VFOU,    WFOU, &
                      RHS_UFOU,RHS_VFOU,RHS_WFOU, RKSTEP )



!- update the fluid velocity by the remaining part dt/2*k2
UFOU = UFOU_NP1 + DTIME/2*RHS_UFOU
VFOU = VFOU_NP1 + DTIME/2*RHS_VFOU
WFOU = WFOU_NP1 + DTIME/2*RHS_WFOU




!!====================================================================
!! 4. Finalizing fluid velocity
!!====================================================================

!- Projection on soloneidal basis
 call PROJ_DIVFREE



!!====================================================================
!! 5. Aliasing control
!!====================================================================

!!- Apply the filter for aliasing control on non-linear terms
do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
  UFOU(I,J,K) = UFOU(I,J,K) * FILTER(I,J,K)
  VFOU(I,J,K) = VFOU(I,J,K) * FILTER(I,J,K)
  WFOU(I,J,K) = WFOU(I,J,K) * FILTER(I,J,K)
  end do
 end do
end do



!!- CPU check
if(MYID == 0) then
 TIME_START = MPI_WTIME()
end if


!- Back in physical space
TMPFOU = UFOU
call P3DFFT_BTRAN_C2R(TMPFOU,UFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)       

TMPFOU = VFOU 
call P3DFFT_BTRAN_C2R(TMPFOU,VFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)

TMPFOU = WFOU
call P3DFFT_BTRAN_C2R(TMPFOU,WFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)       







!!====================================================================
!! X. Passive scalar equation
!!====================================================================
!if(SOLVE_SCALAR) then
!
! call ADV_SCALAR
!
! TMPFOU = THETAFOU
! call P3DFFT_BTRAN_C2R(TMPFOU,THETA(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)       
!
!end if







!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()
 CPU_FLUID(8) = CPU_FLUID(8) + TIME_END - TIME_START

 CPU_FLUID(1) = CPU_FLUID(1) + TIME_END - TIME_START2
end if




end subroutine ADV_FLUID
