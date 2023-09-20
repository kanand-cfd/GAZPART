subroutine WALL_BOUNDARY_CONDITION_NEW(NCYCLE, I, IG, IDP,  &
                                   XP1, YP1, ZP1,           &
                                   quat,                    &
                                   wall_normal, impact_arm, &
                                   delta_x,                 &
                                   isIntersect)

use DNS_DIM
use param_phys
use mod_quaternion
use PARTICLE_PARALLEL
use GEOMETRIC_VARIABLE

implicit none

!=================================================!
! 		            Input Variables               !
!=================================================!
integer, intent(in) :: I, IG, IDP, NCYCLE
real(kind=8), intent(in) :: XP1, YP1, ZP1
type(quaternion), intent(in) :: quat
real(kind=8), dimension(ndim,1), intent(inout) :: wall_normal, impact_arm
real(kind=8), intent(inout) :: delta_x 
logical, intent(inout) :: isIntersect 
!=================================================!

! Wall Parameters
!real(kind=8), dimension(ndim,1) :: wall_normal ! normal vector
real(kind=8), dimension(ndim,1) :: wall_pos ! Position vector

! Closest point on plane from Ellipsoid
real(kind=8), dimension(ndim,1) :: plane_pt

! Maximum Encroachment point on Ellipsoid
real(kind=8), dimension(ndim,1) :: ellip_pt

! Margin function for Ellipsoid
real(kind=8) :: margin

! Impact Arm
!real(kind=8), dimension(ndim,1) :: impact_arm
!=================================================!
!=================================================!

! Call subroutine for closest point between ellipsoid and plane !
! by passing the plane location and particle location. Check if !
! the closest point on plane is inside or outside the ellipsoid !
! based on the margin function (i.e. if margin < 1 -> inside).  !
! We the point in the ellipsoid which has maximum encroachment  !
! Update particle position and  Apply collision forces.  Repeat !
! the procedure for all the 6 walls.							!

impact_arm(:,:) = 0.0
isIntersect = .false.


if(XP1 < LXMAX/2.0) then
!!=====================WALL#1======================!
! Wall at (0,0,0) with normal (1,0,0)

    wall_pos(1,1) = 0.0
    wall_pos(2,1) = YP1 !0.0
    wall_pos(3,1) = ZP1 !0.0

    wall_normal(1,1) = 1.0
    wall_normal(2,1) = 0.0
    wall_normal(3,1) = 0.0

else

!!=====================WALL#2======================!
! Wall at (LX,0,0) with normal (-1,0,0)

    wall_pos(1,1) = LXMAX
    wall_pos(2,1) = YP1 !0.0
    wall_pos(3,1) = ZP1 !0.0

    wall_normal(1,1) = -1.0
    wall_normal(2,1) =  0.0
    wall_normal(3,1) =  0.0

end if



! call closest point
call closest_point_plane(XP1, YP1, ZP1, &
                         quat,wall_pos,wall_normal,     &
                         plane_pt,margin, IG, IDP) 


if(margin < 1.0) then

    isIntersect = .true.

    call max_encroachment_ellipsoid(XP1, YP1, ZP1, &
                                    quat,wall_normal,ellip_pt,     &
                                    IG, IDP)


    ! Calculate impact arm
    impact_arm(1,1) = (XP1 - ellip_pt(1,1))
    impact_arm(2,1) = (YP1 - ellip_pt(2,1))
    impact_arm(3,1) = (ZP1 - ellip_pt(3,1))

    delta_x = wall_pos(1,1) - ellip_pt(1,1)


end if     
!!=================================================!
return

end subroutine WALL_BOUNDARY_CONDITION_NEW