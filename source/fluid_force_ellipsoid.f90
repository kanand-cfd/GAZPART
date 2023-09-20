subroutine FLUID_FORCE_ELLIPSOID(IG, IDP, &
                                UPART, VPART, WPART, &
                                UFLUID, VFLUID, WFLUID, &
                                quat, &
                                FDRAG, FLIFT)

use DNS_DIM               
use PARAM_PHYS
use mod_quaternion

implicit none

integer, intent(in) :: IG, IDP

real(kind=8), intent(in) :: UPART, VPART, WPART
real(kind=8), intent(in) :: UFLUID, VFLUID, WFLUID

type(quaternion), intent(in) :: quat

real(kind=8), dimension(ndim, 1), intent(inout) :: FDRAG, FLIFT

!real(kind=8), dimension(ndim, 1) :: FDRAG, FLIFT

real(kind=8) :: CDRAG, CLIFT
real(kind=8) :: CDRAG_PHI0, CDRAG_PHI90

!!! Fitting Parameters for Drag Coefficient !!!
real(kind=8) :: a0, a1, a2, a3, a4, a5, a6, a7, a8

!!! Fitting Parameters for Lift Coefficient !!!
real(kind=8) :: b1, b2, b3, b4, b5, b6, b7, b8, b9, b10

real(kind=8), dimension(ndim, 1) :: UREL            ! relative velocity
real(kind=8), dimension(ndim, 1) :: UFLU_PFRAME     ! fluid velocity in particle frame
real(kind=8), dimension(ndim, 1) :: principal_axis  ! 

real(kind=8) :: phi      ! Angle of Incidence
real(kind=8) :: dequiv   ! Equivalent Diameter
real(kind=8) :: Rep      ! Particle Reynolds number
real(kind=8) :: UREL_mag ! Magnitude of relative velocity
real(kind=8) :: UFLU_mag ! Magnitude of fluid velocity in particle frame
real(kind=8) :: PMASS

!initialise FDRAG, FLIFT
FDRAG(:,:) = ZERO; FLIFT(:,:) = ZERO
!!==================================================================================!!
 
!! Fitting parameters taken from Zastawny et al. (2012) for Ellipsoid 2 (a/b=1.25) !!
!! Drag 
a0 = 1.95
a1 = 18.12; a2 = 1.023; a3 = 4.26; a4 = 0.384
a5 = 21.52; a6 = 0.990; a7 = 2.86; a8 = 0.260

!! Lift
b1 = 0.083; b2 = -0.21; b3 = 1.582; b4 = 0.851; b5  = 1.842
b6 =-0.802; b7 =-0.006; b8 = 0.874; b9 = 0.009; b10 = 0.570
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Mass of the particle
PMASS = RHOP(IG, IDP)*(4.0/3.0)*PPI*(EMAJ_PART(IG,IDP)**3.0)/((APR_PART(IG)**(2.0)))

! Relative velocity of particle relative to fluid 
UREL(1,1) = UFLUID - UPART
UREL(2,1) = VFLUID - VPART
UREL(3,1) = WFLUID - WPART

! Magnitude of Relative velocity
UREL_mag = sqrt(UREL(1,1)**2 + UREL(2,1)**2 + UREL(3,1)**2)

! Equivalent Diameter !! Diameter of a sphere (pi/6*dequiv^3) with same volume as ellipsoid (4/3*pi*a^3/beta^2)
dequiv = (2.0*EMAJ_PART(IG,IDP))/(APR_PART(IG)**(2.0/3.0))

! Particle Reynolds Number
Rep = (UREL_mag*dequiv)/VISC

! Transform the fluid velocity from the fixed frame to the particle frame
! Initialize fluid velocity vector with fluid velocity from fixed frame
UFLU_PFRAME(1,1) = UFLUID - UPART
UFLU_PFRAME(2,1) = VFLUID - VPART
UFLU_PFRAME(3,1) = WFLUID - WPART

call transform_basis(UFLU_PFRAME, conj_q(quat), shape(UFLU_PFRAME))


phi = abs(ATAN((sqrt(UFLU_PFRAME(2,1)**2 + UFLU_PFRAME(3,1)**2))/(UFLU_PFRAME(1,1))))

UFLU_mag = sqrt(UFLU_PFRAME(1,1)**2 + UFLU_PFRAME(2,1)**2 + UFLU_PFRAME(3,1)**2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Calculation of Drag Coefficient and Lift Coefficient using correlations     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! C_D
CDRAG_PHI0  = (a1 / Rep**a2) + (a3 / Rep**a4)

CDRAG_PHI90 = (a5 / Rep**a6) + (a7 / Rep**a8) 

CDRAG = CDRAG_PHI0 + (CDRAG_PHI90 - CDRAG_PHI0) * (sin(phi))**a0

! C_L
CLIFT = ( (b1 / Rep**b2) + (b3 / Rep**b4)) * (sin(phi))**(b5 + b6*(Rep**b7)) * (cos(phi))**(b8 + b9*(Rep**b10))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if(APR_PART(IG) == 1.0) then

!    CLIFT = 0.0

!    if(Rep < 1000.0) then

!        CDRAG = (24.0/Rep)*(1.0+0.15*Rep**0.687)

!    else

!        CDRAG = 0.44

!    end if

!end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                          Calculation of Drag and Lift Force                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FDRAG(1,1) = (1.0/2.0) * (RHOF/PMASS) * (1.0/4.0) * (PPI * dequiv**2) * CDRAG * UFLU_mag * UFLU_PFRAME(1,1)

FDRAG(2,1) = (1.0/2.0) * (RHOF/PMASS) * (1.0/4.0) * (PPI * dequiv**2) * CDRAG * UFLU_mag * UFLU_PFRAME(2,1)

FDRAG(3,1) = (1.0/2.0) * (RHOF/PMASS) * (1.0/4.0) * (PPI * dequiv**2) * CDRAG * UFLU_mag * UFLU_PFRAME(3,1)


FLIFT(1,1) = (1.0/2.0) * (RHOF/PMASS) * (1.0/4.0) * (PPI * dequiv**2) * CLIFT &
            * UFLU_mag**2 * sin(phi) * sign(1.0, -UFLU_PFRAME(1,1))

FLIFT(2,1) = (1.0/2.0) * (RHOF/PMASS) * (1.0/4.0) * (PPI * dequiv**2) * CLIFT &
            * UFLU_mag**2 * cos(phi) * (UFLU_PFRAME(2,1)/sqrt(UFLU_PFRAME(2,1)**2 + UFLU_PFRAME(3,1)**2))

FLIFT(3,1) = (1.0/2.0) * (RHOF/PMASS) * (1.0/4.0) * (PPI * dequiv**2) * CLIFT &
            * UFLU_mag**2 * cos(phi) * (UFLU_PFRAME(3,1)/sqrt(UFLU_PFRAME(2,1)**2 + UFLU_PFRAME(3,1)**2))



!! Total Fluid Force on Particle !!
!FFLUID(1,1)  = FDRAG(1, 1) + FLIFT(1,1)
!FFLUID(2,1)  = FDRAG(2, 1) + FLIFT(2,1)
!FFLUID(3,1)  = FDRAG(3, 1) + FLIFT(3,1) + GRAVITY(IG)
!write(*,*) 'FFLUID ', FFLUID


! Transform force back to fixed frame
call transform_basis(FDRAG, quat, shape(FDRAG))

call transform_basis(FLIFT, quat, shape(FLIFT))


if(CDRAG /= CDRAG) then

    write(*,*) UPART, VPART, WPART, UFLUID, VFLUID, WFLUID
    write(*,*) UFLUID - UPART, VFLUID - VPART, WFLUID - WPART
    write(*,*) phi, Rep, UFLU_PFRAME
    write(*,*) quat

    stop 'In Fluid force'

end if 

return
end subroutine FLUID_FORCE_ELLIPSOID