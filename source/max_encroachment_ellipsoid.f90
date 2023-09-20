subroutine max_encroachment_ellipsoid(XP, YP, ZP, quat, &
                                      normal, point,    &
                                      IG, IDP)

use param_phys
use DNS_DIM
use mod_quaternion

implicit none

!=====================================================!
real(kind=8), intent(in) :: XP, YP, ZP

type(quaternion), intent(in) :: quat

real(kind=8), dimension(ndim, 1), intent(in) :: normal

real(kind=8), dimension(ndim, 1), intent(inout) :: point

integer, intent(in) :: IG, IDP
!=====================================================!

! Ellipsoid Center Position Vector
real(kind=8), dimension(ndim, 1) :: ELLC 

! Ellipsoid Matrix and it inverse
real(kind=8), dimension(ndim, ndim) :: ELL, inv_ELL

! Transpose of wall normal vector
real(kind=8), dimension(1, ndim) :: tr_normal

! Product of inversed ELL and wall normal
real(kind=8), dimension(ndim, 1) :: inv_ELL_N

real(kind=8), dimension(1,1) :: trN_invELL_N

real(kind=8) :: sqrt_trN_invELL_N


! Margin function
real(kind=8), dimension(ndim,1) :: vec

real(kind=8), dimension(ndim, 1) :: ELL_

real(kind=8), dimension(1, 1) :: margin_vec

real(kind=8) :: margin
!==============================================================!
ELL_A = EMAJ_PART(IG,IDP)

lambda = APR_PART(IG)

! Semi Minor axes
ELL_B = ELL_A/lambda
ELL_C = ELL_B

! Initialize the Ellipsoid Center
ELLC(1,1) = XP
ELLC(2,1) = YP
ELLC(3,1) = ZP

! Initialize the ellipsoid matrix
ELL(:,:) = 0.0
ELL(1,1) = 1.0/(ELL_A*ELL_A)
ELL(2,2) = 1.0/(ELL_B*ELL_B)
ELL(3,3) = 1.0/(ELL_C*ELL_C)

! Transform the Ellipsoid Matrix acc. to Orientation
call transform_basis(ELL, quat, shape(ELL))

! Invert ELL
call invert_ndim3_matrix(ELL, inv_ELL)

! Transpose of wall normal vector
tr_normal = transpose(normal)

! [ELL^(-1)]N
inv_ELL_N = matmul(inv_ELL, normal)

! trN [ELL^(-1)] N
trN_invELL_N = matmul(tr_normal, inv_ELL_N)

sqrt_trN_invELL_N = sqrt(trN_invELL_N(1,1))

point = ELLC - (inv_ELL_N/sqrt_trN_invELL_N)



!!=== Calculate margin function ===!!
vec = point - ELLC

ELL_ = matmul(ELL, vec)


margin_vec = matmul(transpose(vec), ELL_)


margin = margin_vec(1,1)

!write(*,*) " "
!write(*,*) 'Ellipsoid Encroachment Point :', point
!write(*,*) 'Ellipsoid Center Position :', ELLC
!write(*,*) 'Margin :', margin

!stop    
end subroutine max_encroachment_ellipsoid