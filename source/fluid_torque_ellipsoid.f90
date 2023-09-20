subroutine FLUID_TORQUE_ELLIPSOID(I, IG, IDP, &
                                  OX, OY, OZ, &
                                  QUAT      , &
                                  TORQUE_PITCHING, TORQUE_ROTATION)

use DNS_DIM
use PARAM_PHYS
use PARTICLE_PARALLEL
use FLUID_VARIABLE
use mod_quaternion
use GEOMETRIC_VARIABLE

implicit none

integer, intent(in) :: I, IG, IDP

real(kind=8), intent(in) :: OX, OY, OZ

type(quaternion), intent(in) :: QUAT

real(kind=8), dimension(ndim, 1), intent(inout) :: TORQUE_PITCHING, TORQUE_ROTATION

real(kind=8) :: CTPITCH
real(kind=8) :: CTROTATION, CTROTATION_2

!! Fitting Parameters for Pitching torque coefficient !!
real(kind=8) :: c1, c2, c3, c4, c5, c6, c7, c8, c9, c10

!! Fitting Parameters for Rotation torque coefficient !!
real(kind=8) :: r1, r2, r3, r4
real(kind=8) :: r1_, r2_, r3_, r4_

real(kind=8), dimension(ndim, 1) :: UREL            ! relative velocity
real(kind=8), dimension(ndim, 1) :: UFLU_PFRAME     ! fluid velocity in particle frame
real(kind=8), dimension(ndim, 1) :: OREL            ! Relative rotational velocity
real(kind=8), dimension(ndim, 1) :: curl_UFLU       ! curl of fluid velocity

real(kind=8) :: phi      ! Angle of Incidence
real(kind=8) :: dequiv   ! Equivalent Diameter
real(kind=8) :: Rep      ! Particle Reynolds number
real(kind=8) :: UREL_mag ! Magnitude of relative velocity
real(kind=8) :: UFLU_mag ! Magnitude of fluid velocity in particle frame
real(kind=8) :: Re_R     ! Rotational Particle Reynolds Number
real(kind=8) :: OREL_MAG ! Magnitude of Relative Rotational Velocity
real(kind=8) :: OREL_YZ
real(kind=8) :: PMASS

!!==================================================================================!!
real(kind=8), dimension(ndim, 1) :: FDRAG, FLIFT, XCP

real(kind=8) :: CDRAG, CLIFT
real(kind=8) :: CDRAG_PHI0, CDRAG_PHI90

!!! Fitting Parameters for Drag Coefficient !!!
real(kind=8) :: a0, a1, a2, a3, a4, a5, a6, a7, a8

!!! Fitting Parameters for Lift Coefficient !!!
real(kind=8) :: b1, b2, b3, b4, b5, b6, b7, b8, b9, b10
!!==================================================================================!!

!initialise FTORQUE 
!FTORQUE(:,:) = ZERO
!!==================================================================================!!
!! Fitting parameters taken from Zastawny et al. (2012) for Ellipsoid 2 (a/b=1.25) !!
! Drag 
a0 = 1.95
a1 = 18.12; a2 = 1.023; a3 = 4.26; a4 = 0.384
a5 = 21.52; a6 = 0.990; a7 = 2.86; a8 = 0.260

!! Lift
b1 = 0.083; b2 = -0.21; b3 = 1.582; b4 = 0.851; b5  = 1.842
b6 =-0.802; b7 =-0.006; b8 = 0.874; b9 = 0.009; b10 = 0.570

!! Pitching Torque
c1 = 0.935; c2 = 0.146; c3 = -0.469; c4 = 0.145; c5  = 0.116
c6 = 0.748; c7 = 0.041; c8 =  0.221; c9 = 0.657; c10 = 0.044

!! Rotation Torque (Mode 1)
r1 = 0.573; r2 =-0.154; r3 = 116.61; r4 = 1.0

!! Rotation Torque (Mode 2)
r1_ = 1.244; r2_ = 0.239; r3_ = 378.12; r4_ = 0.789
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Mass of the particle
PMASS = RHOP(IG, IDP)*(4.0/3.0)*PPI*(EMAJ_PART(IG,IDP)**3.0)/((APR_PART(IG)**(2.0)))

! Relative velocity of particle relative to fluid 
UREL(1,1) = PART(I,IG)%UFAP - PART(I,IG)%UP
UREL(2,1) = PART(I,IG)%VFAP - PART(I,IG)%VP
UREL(3,1) = PART(I,IG)%WFAP - PART(I,IG)%WP

! Magnitude of Relative velocity
UREL_mag = sqrt(UREL(1,1)**2 + UREL(2,1)**2 + UREL(3,1)**2)

! Equivalent Diameter !! Diameter of a sphere (pi/6*dequiv^3) with same volume as ellipsoid (4/3*pi*a^3/beta^2)
dequiv = 2.0*EMAJ_PART(IG,IDP)/(APR_PART(IG)**(2.0/3.0))

! Particle Reynolds Number
Rep = (UREL_mag*dequiv)/VISC

! Transform the fluid velocity from the fixed frame to the particle frame
! Initialize fluid velocity vector with fluid velocity from fixed frame
UFLU_PFRAME(1,1) = PART(I,IG)%UFAP - PART(I,IG)%UP
UFLU_PFRAME(2,1) = PART(I,IG)%VFAP - PART(I,IG)%VP
UFLU_PFRAME(3,1) = PART(I,IG)%WFAP - PART(I,IG)%WP

call transform_basis(UFLU_PFRAME, conj_q(QUAT), shape(UFLU_PFRAME))


phi = abs(ATAN((sqrt((UFLU_PFRAME(2,1))**2 + UFLU_PFRAME(3,1)**2))/(UFLU_PFRAME(1,1))))


UFLU_mag = sqrt(UFLU_PFRAME(1,1)**2 + UFLU_PFRAME(2,1)**2 + UFLU_PFRAME(3,1)**2)

!! Relative Rotational Velocity !!
curl_UFLU(1,1) = PART(I,IG)%VORTXAP
curl_UFLU(2,1) = PART(I,IG)%VORTYAP
curl_UFLU(3,1) = PART(I,IG)%VORTZAP 

call transform_basis(curl_UFLU, conj_q(QUAT), shape(curl_UFLU))

OREL(1,1) = 0.5*curl_UFLU(1,1) - OX
OREL(2,1) = 0.5*curl_UFLU(2,1) - OY
OREL(3,1) = 0.5*curl_UFLU(3,1) - OZ

OREL_MAG = sqrt(OREL(1,1)**2 + OREL(2,1)**2 + OREL(3,1)**2)

OREL_YZ = sqrt(OREL(3,1)**2 + OREL(2,1)**2)

!! Rotational Reynolds Number !!
Re_R = (OREL_MAG * dequiv**2)/(VISC)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   Calculation of Pitching & Rotational Torque Coefficient using correlations    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CTPITCH !!
CTPITCH = (c1/Rep**c2 + c3/Rep**c4) * ((sin(phi))**(c5 + c6*(Rep**c7))) * ((cos(phi))**(c8 + c9*(Rep**c10)))

!! CTROTATION !!
CTROTATION = r1*(Re_R**r2) + r3/(Re_R**r4)

CTROTATION_2 = r1_ * (Re_R ** r2_) + r3_/(Re_R ** r4_) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

TORQUE_PITCHING(1,1) = 0.0

TORQUE_PITCHING(2,1) = (1.0/4.0) * (RHOF/PMASS) * (1.0/4.0) * PPI * dequiv**3 * CTPITCH &
                     * (UFLU_mag*UFLU_mag) * (abs(UFLU_PFRAME(3,1))/sqrt(UFLU_PFRAME(2,1)**2 + UFLU_PFRAME(3,1)**2)) &
                     * sign(1.0, UFLU_PFRAME(1,1)*UFLU_PFRAME(3,1))

TORQUE_PITCHING(3,1) = (1.0/4.0) * (RHOF/PMASS) * (1.0/4.0) * PPI * dequiv**3 * CTPITCH &
                     * (UFLU_mag*UFLU_mag) * (abs(UFLU_PFRAME(2,1))/sqrt(UFLU_PFRAME(2,1)**2 + UFLU_PFRAME(3,1)**2)) &
                     * sign(1.0, UFLU_PFRAME(1,1)*UFLU_PFRAME(2,1))


!TORQUE_PITCHING(:,1) = 0.0

TORQUE_ROTATION(1,1) = CTROTATION * (0.5*RHOF/PMASS) * (dequiv/2.0)**5 * abs(OREL(1,1)) * OREL(1,1)

TORQUE_ROTATION(2,1) = CTROTATION_2 * (0.5*RHOF/PMASS) * (dequiv/2.0)**5 * OREL_YZ * OREL(2,1)

TORQUE_ROTATION(3,1) = CTROTATION_2 * (0.5*RHOF/PMASS) * (dequiv/2.0)**5 * OREL_YZ * OREL(3,1)


!TORQUE_ROTATION(:,1) = 0.0

!! Total Fluid Torque
!FTORQUE = (TORQUE_PITCHING + TORQUE_ROTATION)
!if(APR_PART(IG)== 1.0) then

!    TORQUE_PITCHING(:,1) = 0.0

!    TORQUE_ROTATION(:,1) = 0.0

!end if


if(OX /= OX) then
    write(*,*) I
    write(*,*) 'FLUID_TORQUE_ELLIPSOID', TORQUE_ROTATION, TORQUE_PITCHING
    write(*,*) CTROTATION, CTPITCH
    write(*,*) 'OREL', OREL
    write(*,*) ' Re_R ', Re_R
    write(*,*) 'phi', phi
    write(*,*) 'curl_UFLU', curl_UFLU
    write(*,*) 'Omega', OX, OY, OZ
    stop
end if 

return 
end subroutine FLUID_TORQUE_ELLIPSOID