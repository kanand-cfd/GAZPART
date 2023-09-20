!!====================================================================
!!
!!  Forcing scheme for statistically steady turbulent flows
!!
!!====================================================================

subroutine ADD_FORCING3(RHS_UFOU,RHS_VFOU,RHS_WFOU)

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
use PARAM_PHYS
use FORCING
use CHECK_CPU

use P3DFFT

implicit none

!!--------------------------------------------------------------------
!! ARRAYS STATEMENT
!!--------------------------------------------------------------------
!- Right-Hand-Side
double complex, &
   dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)), intent(inout) :: RHS_UFOU
double complex, &
   dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)), intent(inout) :: RHS_VFOU
double complex, &
   dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)), intent(inout) :: RHS_WFOU


!!--------------------------------------------------------------------

!- number of forced waves
integer :: NF

!- Wavenumber
real(kind=8) :: KF


!- index
integer :: I, J, K, MJ, MK

!---------------------------------------------------------------------

do NF = 1, NFORCE_FULL

 if(     IFORCING(NF)>=FSTART(1).and.(IFORCING(NF)<=FEND(1)) &
    .and.JFORCING(NF)>=FSTART(2).and.(JFORCING(NF)<=FEND(2)) &
    .and.KFORCING(NF)>=FSTART(3).and.(KFORCING(NF)<=FEND(3)) ) then

    RHS_UFOU(IFORCING(NF),JFORCING(NF),KFORCING(NF)) = &
         RHS_UFOU(IFORCING(NF),JFORCING(NF),KFORCING(NF)) + FORCING_UFOU(NF)

    RHS_VFOU(IFORCING(NF),JFORCING(NF),KFORCING(NF)) = &
         RHS_VFOU(IFORCING(NF),JFORCING(NF),KFORCING(NF)) + FORCING_VFOU(NF)

    RHS_WFOU(IFORCING(NF),JFORCING(NF),KFORCING(NF)) = &
         RHS_WFOU(IFORCING(NF),JFORCING(NF),KFORCING(NF)) + FORCING_WFOU(NF)

  end if
end do



end subroutine ADD_FORCING3
!
