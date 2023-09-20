!!====================================================================
!!
!! Forcing
!!
!!====================================================================
module FORCING

implicit none

!- Number of forced wavenumber
integer :: NFORCE_FULL

!- Number of forced wavenumber for each CPU
integer :: NFORCE_CPU

!- Integer
integer :: ISAVEFORCE


!- Wavenumber of forcing
real(kind=8) :: KFORCE_MIN
real(kind=8) :: KFORCE_MAX

!- Intensity of forcing
real(kind=8) :: SIGMA_FORCE

!- Time of forcing
real(kind=8) :: TIME_FORCE


!- Forcing force
double complex, dimension(:), allocatable :: FORCING_UFOU
double complex, dimension(:), allocatable :: FORCING_VFOU
double complex, dimension(:), allocatable :: FORCING_WFOU





integer, dimension(:), allocatable :: IFORCING
integer, dimension(:), allocatable :: JFORCING
integer, dimension(:), allocatable :: KFORCING
!
!!- Wave in the full domain
!real(kind=8), dimension(:), allocatable :: KY_FULL
!real(kind=8), dimension(:), allocatable :: KZ_FULL
!
!!- Index for Hermitian symetry
!integer, dimension(:), allocatable :: NHERM



end module FORCING
