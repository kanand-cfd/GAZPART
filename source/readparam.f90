!!====================================================================
!!
!! Routine reading the simulation parameters from the ASCII file
!!                             'param.in'
!!
!!====================================================================

subroutine READPARAM

!!====================================================================
!! 
!!====================================================================

use DNS_DIM
use PARAM_PHYS 
use FORCING
use GEOMETRIC_VARIABLE
use COLLISION_VARIABLE

implicit none

!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- index
integer :: I
!---------------------------------------------------------------------


!- Open file
open(unit=200,file='param.in',status='old')

read(200,*) ! comment line
read(200,*) ! comment line
read(200,*) ! comment line
read(200,*) ! comment line


!!====================================================================
!! 1. NUMERICS
!!====================================================================
read(200,*) ! comment line
read(200,*) ! comment line
read(200,*) ! comment line
read(200,*) SOLVE_FLUID
read(200,*) ! comment line

read(200,*) NX         !- Mesh size
read(200,*) NY         !- Mesh size
read(200,*) NZ         !- Mesh size
read(200,*)            ! comment line
read(200,*) DTIME_USER !- Time step
read(200,*) NCYCLEMAX  !- Maximum of cycle


!!====================================================================
!! 2. Fluid numerics
!!====================================================================
read(200,*)                     ! comment line
read(200,*) FROZEN_FLOW			!- Flag for Frozen Flow with fluid velocity
read(200,*) INIT_FLUID_VELOCITY !- Fluid initiation
read(200,*) STEADY              !- Flag for forcing
read(200,*) SIGMA_FORCE         !- Intensity
read(200,*) TIME_FORCE          !- Time of forcing
read(200,*) KFORCE_MIN          !- Wavenumber
read(200,*) KFORCE_MAX          !- Wavenumber

!!====================================================================
!! 3. Scalar numerics
!!====================================================================
read(200,*) ! comment line
read(200,*) SOLVE_SCALAR
read(200,*) INIT_SCALAR


!!====================================================================
!! 4. Particle numerics
!!====================================================================
read(200,*)            		! comment line
read(200,*) SOLVE_PART 		!- Particle tracking
!read(200,*) DRY_CHANNEL
read(200,*) NIG             !- Number of particle's classes
read(200,*) NPART_FULL 		!- Number of particle
read(200,*) WALL_BOUNDARY 	!- WALL B.C. FLAG
!- Interpolation scheme
! INTERP_SCHEME = 1 : Lagrangian polynomial 1st order
!               = 2 : Lagrangian polynomial 2nd order
!               = 3 : Lagrangian polynomial 3rd order
!               = 4 : Lagrangian polynomial 4th order
read(200,*) INTERP_SCHEME      !- Interpolation scheme

read(200,*) INIT_PART_POSITION !- Particle position initiation
read(200,*) INIT_PART_VELOCITY !- Particle velocity initiation
read(200,*) INIT_PART_SCALAR   !- Particle scalar initiation




allocate(RHOP_USER(NIG))
!allocate(DPART_USER(NIG))

!- Ellipsoid Dimension Parameter
allocate(EMAJ_PART_USER(NIG)) !- ELLIPSOID MAJOR AXIS
allocate(APR_PART(NIG))       !- ELLIPSOID ASPECT RATIO

allocate(PARTDEF(NIG))
allocate(GRAVITY(NIG))
allocate(CP_PART(NIG))
allocate(POLYDISP(NIG))
allocate(NQ_NP(NIG))
allocate(RHOQ_RHOP(NIG))
allocate(DQ_DP(NIG))

allocate(COLLIDEF(NIG))
allocate(ECP(NIG))
allocate(MUP(NIG))


!!====================================================================
!! 5. Outputing
!!====================================================================
read(200,*) ! comment line
read(200,*) ! comment line
read(200,*) ! comment line


read(200,*) READSTAT     !- read statistics from previous run
read(200,*) STAT_TIME    !- Flag to perform averaging over the time
read(200,*)              ! comment line
                           !- Level of statistics for the fluid
read(200,*) LEVEL_STFLU !- =1: One point statistics (mean, Reynolds stress, 3rd & 4th-order)
                        !- =2: Divergence, dissipation
                        !- =3: Skewness, Flatness, PDF of fluid velocity & gradients
read(200,*)              ! comment line
read(200,*) LEVEL1_STSCL !- One point statistics (mean, Reynolds stress, 3rd & 4th-order)
read(200,*) LEVEL2_STSCL !- Dissipation and gradients
read(200,*)              ! comment line
                         !- Level of statistics for the particle
read(200,*) LEVEL1_STPAR !- One point statistics (mean, Reynolds stress, 3rd & 4th-order)
read(200,*) LEVEL2_STPAR !- Lagrangian functiond
read(200,*) LEVEL3_STPAR !- Spatial distribution
read(200,*) LEVELX_STPAR !- Statistics along X
read(200,*) LEVEL1_STPDF !- Calculation of PDFS

!- Outputing parameters
!---------------------
read(200,*) ! comment line
read(200,*) ! comment line
read(200,*) ! comment line


read(200,*) FOUT0       !- Screen printing
read(200,*) FOUT1       !- Fluid Solution printing
read(200,*) FOUT2       !- Stat printing
read(200,*) FOUT3       !- Part Solution printing
read(200,*) FOUT4       !- Printing X-statistics
read(200,*) FOUT5       !- To be defined
read(200,*) ENSIGHT_OUT !- Ensight's format velocity field
read(200,*) NSLICE      !- Number of Slices for X-statistics
read(200,*) REFINE_MESH !- Geometric ratio for slice spacing


!!====================================================================
!! 6. Fluid properties
!!====================================================================
read(200,*) ! comment line
read(200,*) ! comment line
read(200,*) ! comment line


read(200,*) LXMAX !- Computational domain size
read(200,*) LYMAX !- Computational domain size
read(200,*) LZMAX !- Computational domain size
read(200,*) VISC  !- Molecular viscosity
read(200,*) RHOF  !-  Fluid density


!!====================================================================
!! 7. Scalar field properties
!!====================================================================
read(200,*)          ! comment line
read(200,*)          ! comment line
read(200,*)          ! comment line
read(200,*) DIFF_SCL !- Scalar diffusivity
read(200,*) CP_SCL   !- Specific heat of scalar
read(200,*) GRAD_SCL !- Mean gradient of scalar


!!====================================================================
!! 8. Particle properties
!!====================================================================
read(200,*) ! comment line
read(200,*) ! comment line
read(200,*) ! comment line

do I=1, NIG

    read(200,*) PARTDEF(I), COLLIDEF(I) !- Particle definition
    !read(200,*) DPART_USER(I)   !- Diameter

    !- Ellipsoid Dimension Parameter
    read(200,*) EMAJ_PART_USER(I) !- ELLIPSOID MAJOR AXIS
    read(200,*) APR_PART(I)       !- ELLIPSOID ASPECT RATIO

    read(200,*) RHOP_USER(I)    !- Density
    read(200,*) ECP(I), MUP(I)    !- Particle Coefficient of Restitution
    read(200,*) GRAVITY(I) !- Gravity
    read(200,*) CP_PART(I) !- "Specific heat"
    read(200,*) POLYDISP(I)!- Number of Polydisperse class for stat
    read(200,*) NQ_NP(I), DQ_DP(I), RHOQ_RHOP(I) 
    read(200,*)            ! comment line


    if(POLYDISP(I)==1) then
        NQ_NP(I) = 0.0
        DQ_DP(I) = 1.0
        RHOQ_RHOP(I) = 1.0
    end if

end do

!!====================================================================
!! 8. Wall properties
!!====================================================================
read(200,*) ECW     !- Restitution Coefficient of Wall
read(200,*) MUW     !- Friction Coefficient of Wall
read(200,*) BETAW   !- Tangential Restitution Coefficient of Wall

read(200,*) 
read(200,*) FILTERING, KCUT

close(200)



!- Forbid collision for motionless and fluid elements
do I=1, NIG

 if(PARTDEF(I)<2) COLLIDEF(I) = 0

end do

!- Activate collision solver if at least one class have collision
SOLVE_COLLISION = .false.
if(sum(COLLIDEF(:))>0) SOLVE_COLLISION = .true.

PERIODICITY = .true.
if (WALL_BOUNDARY) PERIODICITY = .false.


!- Solution printing
if(FOUT1>NCYCLEMAX) FOUT1 = NCYCLEMAX

!- Statistic printing
if(FOUT2>NCYCLEMAX) FOUT2 = NCYCLEMAX

!- Statistic printing
if(FOUT3>NCYCLEMAX) FOUT3 = NCYCLEMAX

!- maximum number of particle class
if(NIG>20) then
 write(UNIT_INFO(1),*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 write(UNIT_INFO(1),*)'!!            ERROR'
 write(UNIT_INFO(1),*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 write(UNIT_INFO(1),*)'!!    NIG=',NIG,' > 20'
 write(UNIT_INFO(1),*)'!!'
 write(UNIT_INFO(1),*)'!! Problem for files managing see openclose.f90'
 write(UNIT_INFO(1),*)'!!'
 write(UNIT_INFO(1),*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 stop
end if


if(MYID==0) write(*,*) 'Parameters reading --> OK'


end subroutine READPARAM
