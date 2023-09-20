subroutine STAT_PARTICLE_COLLISION(IG, IDP1,            &
								   VRELX, VRELY, VRELZ, &
								   quat1, quat2,        &
								   normal)

use PARTICLE_PARALLEL
use DNS_DIM               
use PARAM_PHYS
use GEOMETRIC_VARIABLE
use COLLISION_VARIABLE
use STATISTICS
use mod_quaternion

implicit none

!========================================================!
integer, intent(in) :: IG, IDP1

real(kind=8), intent(in) :: VRELX, VRELY, VRELZ

type(quaternion), intent(in) :: quat1, quat2

real(kind=8), dimension(ndim, 1), intent(in) :: normal
!========================================================!

real(kind=8) :: NRM_VREL, NVREL, THETA

real(kind=8), dimension(NSTAT,NIG,POLYDISPMAX, NPDF) :: MEAN_PART_COL_PDF

real(kind=8) :: MIN_THETA, MAX_THETA
real(kind=8) :: MIN_NR, MAX_NR, MIN_NV, MAX_NV 
real(kind=8) :: DPDF_NV, DPDF_NR, DPDF_THETA

real(kind=8), dimension(ndim, 1) :: principal_axis_p1, principal_axis_p2, principal_axis_p12

real(kind=8) :: mod_principal_axis1, mod_principal_axis2, mod_principal_axis12

real(kind=8) :: dcos_x1, dcos_y1, dcos_z1
real(kind=8) :: dcos_x2, dcos_y2, dcos_z2
real(kind=8) :: dcos_x12, dcos_y12, dcos_z12

real(kind=8) :: MIN_DCOS, MAX_DCOS, DPDF_DCOS

integer :: IPDF
!========================================================!

MEAN_PART_COL_PDF(:,:,:,:) = ZERO

MIN_NV = 0.0
MAX_NV = 5.0
DPDF_NV = (MAX_NV - MIN_NV)/real(NPDF)

MIN_NR = 0.0
MAX_NR = 5.0
DPDF_NR = (MAX_NR - MIN_NR)/real(NPDF)

MIN_THETA = 0.0
MAX_THETA = 1.5*PPI
DPDF_THETA = (MAX_THETA - MIN_THETA)/real(NPDF)

MIN_DCOS = -1.0
MAX_DCOS =  1.0
DPDF_DCOS = (MAX_DCOS - MIN_DCOS)/real(NPDF)
!========================================================!

! Dot Product of Relative Velocity and Vector linking encroachment points)
NVREL = (VRELX*normal(1,1) + VRELY*normal(2,1) + VRELZ*normal(3,1))

! Modulus of relative velocity
NRM_VREL = sqrt(VRELX*VRELX + VRELY*VRELY + VRELZ*VRELZ)

! Collision angle
THETA = acos(NVREL/NRM_VREL)


! Radial Relative Velocity
IPDF = int((abs(NVREL) - MIN_NV)/DPDF_NV) + 1
if(IPDF < 1) IPDF = 1
if(IPDF > NPDF) IPDF = NPDF
MEAN_PART_COL_PDF(1, IG, IDP1, IPDF) =  MEAN_PART_COL_PDF(1, IG, IDP1, IPDF) + 1.0 


! Modulus of Relative Velocity 
IPDF = int((NRM_VREL - MIN_NR)/DPDF_NR) + 1
if(IPDF < 1) IPDF = 1
if(IPDF > NPDF) IPDF = NPDF
MEAN_PART_COL_PDF(2, IG, IDP1, IPDF) =  MEAN_PART_COL_PDF(2, IG, IDP1, IPDF) + 1.0 


! Radial Relative Velocity
IPDF = int((THETA - MIN_THETA)/DPDF_THETA) + 1
if(IPDF < 1) IPDF = 1
if(IPDF > NPDF) IPDF = NPDF
MEAN_PART_COL_PDF(3, IG, IDP1, IPDF) =  MEAN_PART_COL_PDF(3, IG, IDP1, IPDF) + 1.0 

!========================================================!

! Orientation of paritcle 1
principal_axis_p1(1,1) = 1.0
principal_axis_p1(2,1) = 0.0
principal_axis_p1(3,1) = 0.0

call transform_basis(principal_axis_p1, quat1, shape(principal_axis_p1))


mod_principal_axis1 = sqrt(principal_axis_p1(1,1)**2 + principal_axis_p1(2,1)**2 + principal_axis_p1(3,1)**2)

dcos_x1 = (principal_axis_p1(1,1)/mod_principal_axis1)
dcos_y1 = (principal_axis_p1(2,1)/mod_principal_axis1)
dcos_z1 = (principal_axis_p1(3,1)/mod_principal_axis1)

IPDF = int((dcos_x1 - MIN_DCOS)/DPDF_DCOS) + 1
if(IPDF < 1) IPDF = 1
if(IPDF > NPDF) IPDF = NPDF
MEAN_PART_COL_PDF(4, IG, IDP1, IPDF) =  MEAN_PART_COL_PDF(4, IG, IDP1, IPDF) + 1.0 

IPDF = int((dcos_y1 - MIN_DCOS)/DPDF_DCOS) + 1
if(IPDF < 1) IPDF = 1
if(IPDF > NPDF) IPDF = NPDF
MEAN_PART_COL_PDF(5, IG, IDP1, IPDF) =  MEAN_PART_COL_PDF(5, IG, IDP1, IPDF) + 1.0


IPDF = int((dcos_z1 - MIN_DCOS)/DPDF_DCOS) + 1
if(IPDF < 1) IPDF = 1
if(IPDF > NPDF) IPDF = NPDF
MEAN_PART_COL_PDF(6, IG, IDP1, IPDF) =  MEAN_PART_COL_PDF(6, IG, IDP1, IPDF) + 1.0


! Orientation of paritcle 2
principal_axis_p2(1,1) = 1.0
principal_axis_p2(2,1) = 0.0
principal_axis_p2(3,1) = 0.0

call transform_basis(principal_axis_p2, quat2, shape(principal_axis_p2))

mod_principal_axis2 = sqrt(principal_axis_p2(1,1)**2 + principal_axis_p2(2,1)**2 + principal_axis_p2(3,1)**2)

dcos_x2 = (principal_axis_p2(1,1)/mod_principal_axis2)
dcos_y2 = (principal_axis_p2(2,1)/mod_principal_axis2)
dcos_z2 = (principal_axis_p2(3,1)/mod_principal_axis2)


IPDF = int((dcos_x2 - MIN_DCOS)/DPDF_DCOS) + 1
if(IPDF < 1) IPDF = 1
if(IPDF > NPDF) IPDF = NPDF
MEAN_PART_COL_PDF(7, IG, IDP1, IPDF) =  MEAN_PART_COL_PDF(7, IG, IDP1, IPDF) + 1.0 

IPDF = int((dcos_y2 - MIN_DCOS)/DPDF_DCOS) + 1
if(IPDF < 1) IPDF = 1
if(IPDF > NPDF) IPDF = NPDF
MEAN_PART_COL_PDF(8, IG, IDP1, IPDF) =  MEAN_PART_COL_PDF(8, IG, IDP1, IPDF) + 1.0


IPDF = int((dcos_z2 - MIN_DCOS)/DPDF_DCOS) + 1
if(IPDF < 1) IPDF = 1
if(IPDF > NPDF) IPDF = NPDF
MEAN_PART_COL_PDF(9, IG, IDP1, IPDF) =  MEAN_PART_COL_PDF(9, IG, IDP1, IPDF) + 1.0


! Relative Orientation of paritcle 2
principal_axis_p12 = principal_axis_p1 - principal_axis_p2

mod_principal_axis12 = sqrt(principal_axis_p12(1,1)**2 + principal_axis_p12(2,1)**2 + principal_axis_p12(3,1)**2)

dcos_x12 = (principal_axis_p12(1,1)/mod_principal_axis12)
dcos_y12 = (principal_axis_p12(2,1)/mod_principal_axis12)
dcos_z12 = (principal_axis_p12(3,1)/mod_principal_axis12)

IPDF = int((dcos_x12 - MIN_DCOS)/DPDF_DCOS) + 1
if(IPDF < 1) IPDF = 1
if(IPDF > NPDF) IPDF = NPDF
MEAN_PART_COL_PDF(10, IG, IDP1, IPDF) =  MEAN_PART_COL_PDF(10, IG, IDP1, IPDF) + 1.0 

IPDF = int((dcos_y12 - MIN_DCOS)/DPDF_DCOS) + 1
if(IPDF < 1) IPDF = 1
if(IPDF > NPDF) IPDF = NPDF
MEAN_PART_COL_PDF(11, IG, IDP1, IPDF) =  MEAN_PART_COL_PDF(11, IG, IDP1, IPDF) + 1.0


IPDF = int((dcos_z12 - MIN_DCOS)/DPDF_DCOS) + 1
if(IPDF < 1) IPDF = 1
if(IPDF > NPDF) IPDF = NPDF
MEAN_PART_COL_PDF(12, IG, IDP1, IPDF) =  MEAN_PART_COL_PDF(12, IG, IDP1, IPDF) + 1.0



NEVEN_COLL_PDF = NEVEN_COLL_PDF + 1
MEAN_TIME_PART_COL_PDF = MEAN_TIME_PART_COL_PDF + MEAN_PART_COL_PDF

    
end subroutine STAT_PARTICLE_COLLISION