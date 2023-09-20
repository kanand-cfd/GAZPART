!!=====================================================================
!!
!! Fluid variable
!!
!!=====================================================================
module fluid_variable

implicit none

!- Fluid velocity at tn
real(kind=8), dimension(:,:,:), allocatable :: UFLU
real(kind=8), dimension(:,:,:), allocatable :: VFLU
real(kind=8), dimension(:,:,:), allocatable :: WFLU

!- Complex velocity
double complex, dimension(:,:,:), allocatable :: UFOU
double complex, dimension(:,:,:), allocatable :: VFOU
double complex, dimension(:,:,:), allocatable :: WFOU


!- Fluid Vorticity Field in Physical space -- Used for Ellipsoid
real(kind=8), dimension(:,:,:), allocatable :: VORTICITY_X
real(kind=8), dimension(:,:,:), allocatable :: VORTICITY_Y
real(kind=8), dimension(:,:,:), allocatable :: VORTICITY_Z


end module fluid_variable
