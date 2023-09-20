subroutine pair_depth(IG, IDP1, IDP2, &
					  x, b1, IQUAT, &
					  b2, JQUAT, &
					  cvec1, pt1, pt2)

use PARTICLE_PARALLEL
use DNS_DIM               
use PARAM_PHYS
use GEOMETRIC_VARIABLE
use COLLISION_VARIABLE
use mod_quaternion

implicit none

integer, intent(in) :: IG, IDP1, IDP2

! Optimisation Parameters
real(kind=8), intent(in) :: x(2)

! Particle 1
real(kind=8), dimension(ndim, 1), intent(in) :: b1 ! Ellipsoid centre
type(quaternion), intent(in) :: IQUAT

! Particle 2
real(kind=8), dimension(ndim, 1), intent(in) :: b2 ! Ellipsoid centre
type(quaternion), intent(in) :: JQUAT

! Output
real(kind=8), dimension(ndim, 1), intent(inout) :: cvec1, pt1, pt2

!=================== Local  Variables =================!
real(kind=8), dimension(ndim, 1):: cvec2

! Ellipsoid matrix after rotation
real(kind=8), dimension(ndim, ndim) :: E1_matrix, E2_matrix

! Scalar normal parameter
real(kind=8):: lambda1, lambda2

! Temporary vectors
real(kind=8), dimension(ndim, 1) :: tmp1, tmp2
!======================================================!

cvec1(1,1) = cos(x(1))*cos(x(2))
cvec1(2,1) = sin(x(1))*cos(x(2))
cvec1(3,1) = sin(x(2))

cvec2 = - cvec1

! Ellipsoid matrix 1 after rotation
call calculate_Ellipsoid_matrix(IG, IDP1, E1_matrix, IQUAT)

! Ellipsoid matrix 2 after rotation
call calculate_Ellipsoid_matrix(IG, IDP2, E2_matrix, JQUAT)

tmp1 = matmul(E1_matrix, cvec1)
tmp2 = matmul(E2_matrix, cvec2)

lambda1 = 0.25*(cvec1(1,1)*tmp1(1,1) + cvec1(2,1)*tmp1(2,1) + cvec1(3,1)*tmp1(3,1))
lambda1 = sqrt(lambda1)

lambda2 = 0.25*(cvec2(1,1)*tmp2(1,1) + cvec2(2,1)*tmp2(2,1) + cvec2(3,1)*tmp2(3,1))
lambda2 = sqrt(lambda2)

! Contact Point 1
pt1 = (1.0/(2.0*lambda1))*matmul(E1_matrix, cvec1) + b1

! Contact Point 2
pt2 = (1.0/(2.0*lambda2))*matmul(E2_matrix, cvec2) + b2

! Penetration Depth
!depth = pt1 - pt2

return
end subroutine pair_depth