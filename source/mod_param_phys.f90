!!====================================================================
!!
!! Physical parameters
!!
!!====================================================================
module PARAM_PHYS

implicit none


!!--------------------------------------------------------------------
!!- Fluid physical parameter
!!--------------------------------------------------------------------
real(kind=8) :: VISC      !- Molecular viscosity
real(kind=8) :: RHOF      !- Fluid density


integer :: FROZEN_FLOW    !- Flag for Frozen Flow with fluid velocity
integer :: INIT_FLUID_VELOCITY
!- Fluid Initiation
! = 0 : Uf = 0
! = 1 : Random value
! = 2 : Uniform fluid velocity (defined in fortran file)
! = 3 : Read fluid velocity field read from thi generating code
! = 4 : Read fluid velocity field read from stored files

real(kind=8), dimension(3) :: UREF !- Initial mean velocity field


!!--------------------------------------------------------------------
!!- Scalar parameter
!!--------------------------------------------------------------------
real(kind=8) :: DIFF_SCL  !- Diffusivity of Passive SCalar
real(kind=8) :: CP_SCL    !- "Specific heat" of scalar
real(kind=8) :: GRAD_SCL  !- Gradient of scalar

integer :: INIT_SCALAR

!!--------------------------------------------------------------------
!! Numerical parameter
!!--------------------------------------------------------------------
real(kind=8) :: DTIME       !- Time step
real(kind=8) :: DTIME_USER  !- User time step
integer      :: NCYCLEMAX   !- Maximum of cycle number

!- Forced turbulence
logical :: STEADY

!!--------------------------------------------------------------------
!! Output parameters
!!--------------------------------------------------------------------
integer :: FOUT0       !- Screen printing
integer :: FOUT1       !- Fluid Solution printing
integer :: FOUT2       !- Stat printing
integer :: FOUT3       !- Particles Solution printing
integer :: FOUT4       !- To define
integer :: FOUT5       !- To define

integer :: ENSIGHT_OUT !- Ensight's format



!!--------------------------------------------------------------------
!! Particle parameters
!!--------------------------------------------------------------------
!!- 
!!----ELLIPSOIDAL PARTICLE----!!
! Dimensions of the ellipsoidal particles
real(kind=8):: ELL_A
real(kind=8):: ELL_B
real(kind=8):: ELL_C

! Aspect Ratio
real(kind=8):: lambda 

!!----------------------------!!

logical :: WALL_BOUNDARY      !- WALL B.C. FLAG 
logical :: PERIODICITY        !- Flag for Periodicity
integer :: INTERP_SCHEME      !- Interpolation scheme
integer :: INIT_PART_POSITION !- Particle position initiation
integer :: INIT_PART_VELOCITY !- Particle velocity initiation
integer :: INIT_PART_SCALAR   !- Particle scalar initiation

integer,      dimension(:), allocatable :: PARTDEF !- Particle kind
real(kind=8), dimension(:), allocatable :: RHOP_USER    !- Particle density
!real(kind=8), dimension(:), allocatable :: DPART_USER   !- Particle diameter

!! Ellipsoidal Particle Dimension Parameters
real(kind=8), dimension(:), allocatable :: EMAJ_PART_USER !-  Major Axis of Ellipsoidal Particle
real(kind=8), dimension(:), allocatable :: APR_PART !- Aspect Ratio

real(kind=8), dimension(:), allocatable :: GRAVITY !- Gravity exeprienced by particle
real(kind=8), dimension(:), allocatable :: CP_PART !- Specific heat of particles
integer,      dimension(:), allocatable :: COLLIDEF  !- Flag for collision
integer,      dimension(:), allocatable :: POLYDISP!- Flag for polydisperse flows
real(kind=8), dimension(:), allocatable :: RHOQ_RHOP  !- rho_q / rho_p
real(kind=8), dimension(:), allocatable :: DQ_DP      !- dq / dp
real(kind=8), dimension(:), allocatable :: NQ_NP      !- Np / Nq 


integer, dimension(:,:), allocatable :: NPCLASS   !- Particle diameter forpolydisperse
!real(kind=8), dimension(:,:), allocatable :: DPART   !- Particle diameter for polydisperse

real(kind=8), dimension(:,:), allocatable :: EMAJ_PART !- Particle major axis for polydisperse

! Moment of Inertia
real(kind=8), dimension(:,:), allocatable :: IPXX
real(kind=8), dimension(:,:), allocatable :: IPYY
real(kind=8), dimension(:,:), allocatable :: IPZZ

real(kind=8), dimension(:,:), allocatable :: RHOP   !- Particle density for polysolid
real(kind=8), dimension(:,:), allocatable :: MPART ! Particle Mass for polysolid

integer                                 :: POLYDISPMAX





end module param_phys
