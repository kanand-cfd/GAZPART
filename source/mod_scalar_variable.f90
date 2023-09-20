!!=====================================================================
!!
!! Scalar variable
!!
!!=====================================================================
module SCALAR_VARIABLE

implicit none

!- Scalar field in physical space
real(kind=8), dimension(:,:,:), allocatable :: THETA

!- Complex scalar field
double complex, dimension(:,:,:), allocatable :: THETAFOU


!- Right-Hand-Side term of Passive scalar equation
double complex, dimension(:,:,:,:), allocatable :: RHS_SCL

end module SCALAR_VARIABLE
