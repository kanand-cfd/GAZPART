!!------------------------PARTICLE_PARALLEL----------------------------!!
!!  This module defines the particle type and all the variables which  !!
!!  are common to all particles (position, velocity, etc.) through  a  !!
!!  particle structure PARTTYPE; In this version of GASPART, we  have  !!
!!  modified this module for working with ELLIPSOIDAL particles, thus  !!
!!  a few more variables like orientation (via quaternion) & ang. vel. !!  
!!---------------------------------------------------------------------!!

module PARTICLE_PARALLEL

use DNS_DIM
use mod_quaternion
use PARAM_PHYS

implicit none

type PARTTYPE

!!- Particle position
 real(kind=8) :: XP
 real(kind=8) :: YP
 real(kind=8) :: ZP

!!- Particle translational velocity
 real(kind=8) :: UP
 real(kind=8) :: VP
 real(kind=8) :: WP

!!- Fluid Velocity at Particle Position
 real(kind=8) :: UFAP
 real(kind=8) :: VFAP
 real(kind=8) :: WFAP

 !!- Fluid Vorticity at Particle Position
 real(kind=8) :: VORTXAP
 real(kind=8) :: VORTYAP
 real(kind=8) :: VORTZAP

 real(kind=8) :: INVTAUP

!!- Particle Orientation via Quaternions
 type(quaternion):: ELLQUAT 

!!- Particle angular velocity 
 real(kind=8):: OMEGAX
 real(kind=8):: OMEGAY
 real(kind=8):: OMEGAZ

!!- Global Angular Velocity
 real(kind=8) :: GOMEGAX
 real(kind=8) :: GOMEGAY
 real(kind=8) :: GOMEGAZ

!!- Particle class with respect to the diameter
 integer :: IDP

!!- Flag for wall rebound in the previous time step
 logical :: COLL_BIDISP
 logical :: COLL_FLAG

!!- Stored particle velocity for Lagrangian correlation
 real(kind=8),dimension(NBLGRMAX) :: UPT0
 real(kind=8),dimension(NBLGRMAX) :: VPT0
 real(kind=8),dimension(NBLGRMAX) :: WPT0

!!- Stored fluid velocity for Lagrangian correlation
 real(kind=8),dimension(NBLGRMAX) :: UFAPT0
 real(kind=8),dimension(NBLGRMAX) :: VFAPT0
 real(kind=8),dimension(NBLGRMAX) :: WFAPT0

!!- Color of the particles
 real(kind=8) :: COLOR

!!====================================================
!!- Specific variable for investigating the subgrid
!!  turbulent scales
 real(kind=8) :: DUFAP
 real(kind=8) :: DVFAP
 real(kind=8) :: DWFAP
 real(kind=8),dimension(NBLGRMAX) :: DUFAPT0
 real(kind=8),dimension(NBLGRMAX) :: DVFAPT0
 real(kind=8),dimension(NBLGRMAX) :: DWFAPT0
!!====================================================

 integer :: MYIDP
! integer :: ID 
end type PARTTYPE

integer :: NDOUBLE = 20 + NBLGRMAX*9 

integer :: NINTEGER = 2
!integer :: NINTEGER = 1

!! MPI TYPE FOR PARTICLE EXCHANGE
integer :: MPI_PARTICLETYPE
integer :: TAG


!- Number of particle exchanged
integer, dimension(:), allocatable :: NBR_EXCHANGE

!- Particle structure
type(PARTTYPE), dimension(:,:), allocatable :: PART

!- Function to calculate if a given point in space is inside or outside the ellipsoid
public :: fmargin  !- Global definition

private :: margin_function !- Private definition

!- Linking the global and private definition
interface fmargin 
    module procedure margin_function 
end interface


contains

  subroutine MPI_PART_TYPE

    integer, dimension(0:1) ::  oldtypes, blockcounts, offsets(0:1)
!    integer(kind=MPI_OFFSET_KIND) :: extent
    integer :: extent

!   Setup description of the NDOUBLE MPI_DOUBLE_PRECISION fields in PARTTYPE  
    offsets(0) = 0 
    oldtypes(0) = MPI_DOUBLE_PRECISION
    blockcounts(0) = NDOUBLE

!  Setup description of the NINTEGER MPI_INTEGER fields in PARTYPE
!  Need to first figure offset by getting size of MPI_DOUBLE_PRECISION 
    call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, extent, ierr) 
    offsets(1) = NDOUBLE * extent 
    oldtypes(1) = MPI_INTEGER 
    blockcounts(1) = NINTEGER 

!  Now define structured type and commit it  
    call MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes, MPI_PARTICLETYPE, ierr) 
    call MPI_TYPE_COMMIT(MPI_PARTICLETYPE, ierr) 


  end subroutine MPI_PART_TYPE


!!- Function definition for margin function procedure
  function margin_function(IG, IDP, point, center, quat) result(m)

    integer, intent(in) :: IG, IDP

    real(kind=8) :: m

    real(kind=8), dimension(1, 1) :: margin_vec

    real(kind=8) , dimension(ndim, 1), intent(in) :: point, center

    type(quaternion), intent(in) :: quat

    real(kind=8), dimension(ndim, ndim) :: ELL, R, Trn_R

    real(kind=8), dimension(ndim, 1) :: ELL_, vec

    !! 
    ELL_A = EMAJ_PART(IG, IDP)
    ELL_B = ELL_A / APR_PART(IG)
    ELL_C = ELL_B

    ! Initialize the ellipsoid matrix
    ELL(:,:) = 0.0
    ELL(1,1) = 1.0/(ELL_A*ELL_A)
    ELL(2,2) = 1.0/(ELL_B*ELL_B)
    ELL(3,3) = 1.0/(ELL_C*ELL_C)

    !write(*,*) "Quaternion : ", quat%a, quat%b, quat%c, quat%d
    !write(*,*) "Norm :", norm_q(q_p)

    ! Initialise Rotation Matrix
    R(1,1) = 1.0 - 2.0*(quat%c*quat%c + quat%d*quat%d)
    R(1,2) = 2.0*(quat%b*quat%c - quat%a*quat%d)
    R(1,3) = 2.0*(quat%b*quat%d + quat%a*quat%c)

    R(2,1) = 2.0*(quat%b*quat%c + quat%a*quat%d)
    R(2,2) = 1.0 - 2.0*(quat%b*quat%b + quat%d*quat%d)
    R(2,3) = 2.0*(quat%c*quat%d - quat%a*quat%b)

    R(3,1) = 2.0*(quat%b*quat%d - quat%a*quat%c)
    R(3,2) = 2.0*(quat%c*quat%d + quat%a*quat%b)
    R(3,3) = 1.0 - 2.0*(quat%b*quat%b + quat%c*quat%c)

    ! Transpose of R 
    Trn_R = TRANSPOSE(R) ! R is orthogonal so Trn_R = Inverse(R)

    ELL = matmul(R, matmul(ELL, Trn_R))

    vec = point - center

    ELL_ = matmul(ELL, vec)

    margin_vec = matmul(transpose(vec), ELL_)

    m = margin_vec(1,1)

end function margin_function


end module PARTICLE_PARALLEL


