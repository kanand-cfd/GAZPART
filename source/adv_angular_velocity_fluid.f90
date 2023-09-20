subroutine EULER_INTEGRATION_FLUID(I, IG, IDP)

use PARAM_PHYS
use DNS_DIM
!use ellipsoid_particle
use PARTICLE_PARALLEL

implicit none

!=========== INPUT VARIABLES =================!

integer, intent(in) :: I, IG, IDP

!type(PARTTYPE), intent(inout) :: particle

!=============================================!
real(kind=8) :: dt

!!============== RK4 Method =================!!
real(kind=8) :: k1_x, k2_x, k3_x, k4_x
real(kind=8) :: k1_y, k2_y, k3_y, k4_y
real(kind=8) :: k1_z, k2_z, k3_z, k4_z

real(kind=8) :: w0_x, w0_y, w0_z
real(kind=8) :: w1_x, w1_y, w1_z
!=============================================!

real(kind=8), dimension(ndim, 1) :: FTORQUE, FTORQUE_ROT, FTORQUE_PITCH

!=============================================!

dt = DTIME

! Initial 
w0_x = PART(I,IG)%OMEGAX
w0_y = PART(I,IG)%OMEGAY
w0_z = PART(I,IG)%OMEGAZ

call FLUID_TORQUE_ELLIPSOID(I, IG, IDP, &
                            w0_x, w0_y, w0_z, &
                            FTORQUE_PITCH, FTORQUE_ROT) !FTORQUE(:,1) = 0.0

FTORQUE = FTORQUE_PITCH + FTORQUE_ROT

! First step
k1_x = dt*(FTORQUE(1,1) + w0_y*w0_z*(IPYY(IG,IDP) - IPZZ(IG,IDP))) / IPXX(IG,IDP) ! omega_x
k1_y = dt*(FTORQUE(2,1) + w0_z*w0_x*(IPZZ(IG,IDP) - IPXX(IG,IDP))) / IPYY(IG,IDP) ! omega_y
k1_z = dt*(FTORQUE(3,1) + w0_x*w0_y*(IPXX(IG,IDP) - IPYY(IG,IDP))) / IPZZ(IG,IDP) ! omega_z

call FLUID_TORQUE_ELLIPSOID(I, IG, IDP, &
                            (w0_x + 0.5*k1_x), (w0_y + 0.5*k1_y), (w0_z + 0.5*k1_z), &
                            FTORQUE_PITCH, FTORQUE_ROT) !FTORQUE(:,1) = 0.0

FTORQUE = FTORQUE_PITCH + FTORQUE_ROT

!Second step
k2_x = dt*(FTORQUE(1,1) + (w0_y + 0.5*k1_y)*(w0_z + 0.5*k1_z)*(IPYY(IG, IDP) - IPZZ(IG, IDP))) / IPXX(IG, IDP) ! omega_x
k2_y = dt*(FTORQUE(2,1) + (w0_z + 0.5*k1_z)*(w0_x + 0.5*k1_x)*(IPZZ(IG, IDP) - IPXX(IG, IDP))) / IPYY(IG, IDP) ! omega_y
k2_z = dt*(FTORQUE(3,1) + (w0_x + 0.5*k1_x)*(w0_y + 0.5*k1_y)*(IPXX(IG, IDP) - IPYY(IG, IDP))) / IPZZ(IG, IDP) ! omega_z

call FLUID_TORQUE_ELLIPSOID(I, IG, IDP, &
                            (w0_x + 0.5*k2_x), (w0_y + 0.5*k2_y), (w0_z + 0.5*k2_z), &
                            FTORQUE_PITCH, FTORQUE_ROT) !FTORQUE(:,1) = 0.0

FTORQUE = FTORQUE_PITCH + FTORQUE_ROT

! Third Step
k3_x = dt*(FTORQUE(1,1) + (w0_y + 0.5*k2_y)*(w0_z + 0.5*k2_z)*(IPYY(IG, IDP) - IPZZ(IG, IDP))) / IPXX(IG, IDP) ! omega_x
k3_y = dt*(FTORQUE(2,1) + (w0_z + 0.5*k2_z)*(w0_x + 0.5*k2_x)*(IPZZ(IG, IDP) - IPXX(IG, IDP))) / IPYY(IG, IDP) ! omega_y
k3_z = dt*(FTORQUE(3,1) + (w0_x + 0.5*k2_x)*(w0_y + 0.5*k2_y)*(IPXX(IG, IDP) - IPYY(IG, IDP))) / IPZZ(IG, IDP) ! omega_z

call FLUID_TORQUE_ELLIPSOID(I, IG, IDP, &
                            (w0_x + k3_x), (w0_y + k3_y), (w0_z + k3_z), &
                            FTORQUE_PITCH, FTORQUE_ROT) !FTORQUE(:,1) = 0.0

FTORQUE = FTORQUE_PITCH + FTORQUE_ROT

! Fourth Step
k4_x = dt*(FTORQUE(1,1) + (w0_y + k3_y)*(w0_z + k3_z)*(IPYY(IG, IDP) - IPZZ(IG, IDP))) / IPXX(IG, IDP) ! omega_x
k4_y = dt*(FTORQUE(2,1) + (w0_z + k3_z)*(w0_x + k3_x)*(IPZZ(IG, IDP) - IPXX(IG, IDP))) / IPYY(IG, IDP) ! omega_y
k4_z = dt*(FTORQUE(3,1) + (w0_x + k3_x)*(w0_y + k3_y)*(IPXX(IG, IDP) - IPYY(IG, IDP))) / IPZZ(IG, IDP) ! omega_z

! Final Step
w1_x = w0_x + (k1_x + 2.0*k2_x + 2.0*k3_x + k4_x)/6.0
w1_y = w0_y + (k1_y + 2.0*k2_y + 2.0*k3_y + k4_y)/6.0
w1_z = w0_z + (k1_z + 2.0*k2_z + 2.0*k3_z + k4_z)/6.0


PART(I,IG)%OMEGAX = w1_x
PART(I,IG)%OMEGAY = w1_y
PART(I,IG)%OMEGAZ = w1_z


end subroutine EULER_INTEGRATION_FLUID