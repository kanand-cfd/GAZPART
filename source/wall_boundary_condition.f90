subroutine WALL_BOUNDARY_CONDITION(NCYCLE, I, IG, IDP)

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
!type(ELL_PART), intent(inout) ::  PART(I,IG)
!=================================================!

! Wall Parameters
real(kind=8), dimension(ndim,1) :: wall_normal ! normal vector
real(kind=8), dimension(ndim,1) :: wall_pos ! Position vector

! Closest point on plane from Ellipsoid
real(kind=8), dimension(ndim,1) :: plane_pt

! Maximum Encroachment point on Ellipsoid
real(kind=8), dimension(ndim,1) :: ellip_pt

! Margin function for Ellipsoid
real(kind=8) :: margin

! Impact Arm
real(kind=8), dimension(ndim,1) :: impact_arm

logical :: flag
!=================================================!

! Variables for updating the position once collision is detected
real(kind=8) :: delta_t

real(kind=8) :: delta_x, delta_y, delta_z
real(kind=8) :: x_wall, y_wall, z_wall
!=================================================!

! Call subroutine for closest point between ellipsoid and plane !
! by passing the plane location and particle location. Check if !
! the closest point on plane is inside or outside the ellipsoid !
! based on the margin function (i.e. if margin < 1 -> inside).  !
! We the point in the ellipsoid which has maximum encroachment  !
! Update particle position and  Apply collision forces.  Repeat !
! the procedure for all the 6 walls.							!


!!=====================WALL#1======================!
! Wall at (0,0,0) with normal (1,0,0)

wall_pos(1,1) = 0.0
wall_pos(2,1) = PART(I,IG)%YP !0.0
wall_pos(3,1) = PART(I,IG)%ZP !0.0

wall_normal(1,1) = 1.0
wall_normal(2,1) = 0.0
wall_normal(3,1) = 0.0

! call closest point
call closest_point_plane(PART(I,IG)%XP, PART(I,IG)%YP, PART(I,IG)%ZP, &
                         PART(I,IG)%ELLQUAT,wall_pos,wall_normal,     &
                         plane_pt,margin, IG, IDP) 


if(margin < 1.0) then

    call max_encroachment_ellipsoid(PART(I,IG)%XP, PART(I,IG)%YP, PART(I,IG)%ZP, &
                                    PART(I,IG)%ELLQUAT,wall_normal,ellip_pt,     &
                                    IG, IDP)


    ! Calculate impact arm
    impact_arm(1,1) = (PART(I,IG)%XP - ellip_pt(1,1))
    impact_arm(2,1) = (PART(I,IG)%YP - ellip_pt(2,1))
    impact_arm(3,1) = (PART(I,IG)%ZP - ellip_pt(3,1))

    !write(*,*) 'wall_pos(1,1)', wall_pos(1,1), 'ellip_pt(1,1)', ellip_pt(1,1)
    !write(*,*) 'delta_x', delta_x
    !write(*,*) 'x_wall', x_wall

    !PART(I,IG)%XP = x_wall


    call closest_point_plane(PART(I,IG)%XP, PART(I,IG)%YP, PART(I,IG)%ZP, &
                             PART(I,IG)%ELLQUAT,wall_pos,wall_normal,     &
                             plane_pt,margin, IG, IDP)

    !write(*,*) "margin ", margin

    !if(abs(margin - 1) > 1.0e-6) then

    !    write(*,*) "Ellipsoid still intersecting wall", PART(I,IG)%XP, PART(I,IG)%YP, PART(I,IG)%ZP
    !    write(*,*) "margin ", margin, abs(margin - 1)

    !    stop

    !end if

    !write(*,*) "Position before collision",    PART(I,IG)%XP,  PART(I,IG)%YP,  PART(I,IG)%ZP
    !write(*,*) "Orientation :",  PART(I,IG)%ELLQUAT 
    !write(*,*) "Wall at (0,0,0) with normal (1,0,0)"

    ! calculate impulse
    call wall_ellipsoid_rebound(PART(I,IG)%UP, PART(I,IG)%VP, PART(I,IG)%WP,             &
                                PART(I,IG)%OMEGAX, PART(I,IG)%OMEGAY, PART(I,IG)%OMEGAZ, &
                                PART(I,IG)%ELLQUAT,impact_arm,wall_normal,               &
                                IG, IDP, NCYCLE, flag)

    !write(*,*) "  "
!    if(flag) then

        ! Update the particle position
!        delta_x = (ellip_pt(1,1) - wall_pos(1,1))!(wall_pos(1,1) - ellip_pt(1,1))

!        x_wall =  PART(I,IG)%XP - delta_x

!    end if


end if     
!!=================================================!

!stop
!!=====================WALL#2======================!
! Wall at (LX,0,0) with normal (-1,0,0)

wall_pos(1,1) = LXMAX
wall_pos(2,1) = PART(I,IG)%YP !0.0
wall_pos(3,1) = PART(I,IG)%ZP !0.0

wall_normal(1,1) = -1.0
wall_normal(2,1) =  0.0
wall_normal(3,1) =  0.0

call closest_point_plane(PART(I,IG)%XP, PART(I,IG)%YP, PART(I,IG)%ZP, &
                         PART(I,IG)%ELLQUAT,wall_pos,wall_normal,     &
                         plane_pt,margin, IG, IDP) 


if(margin < 1.0) then

    call max_encroachment_ellipsoid(PART(I,IG)%XP, PART(I,IG)%YP, PART(I,IG)%ZP, &
                                    PART(I,IG)%ELLQUAT,wall_normal,ellip_pt,     &
                                    IG, IDP)

    ! Calculate impact arm
    impact_arm(1,1) = (PART(I,IG)%XP - ellip_pt(1,1))
    impact_arm(2,1) = (PART(I,IG)%YP - ellip_pt(2,1))
    impact_arm(3,1) = (PART(I,IG)%ZP - ellip_pt(3,1))


    

    !write(*,*) 'wall_pos(1,1)', wall_pos(1,1), 'ellip_pt(1,1)', ellip_pt(1,1)
    !write(*,*) 'delta_x', delta_x
    !write(*,*) 'x_wall', x_wall

    !PART(I,IG)%XP = x_wall


    call closest_point_plane(PART(I,IG)%XP, PART(I,IG)%YP, PART(I,IG)%ZP, &
                             PART(I,IG)%ELLQUAT,wall_pos,wall_normal,     &
                             plane_pt,margin, IG, IDP) 

    !write(*,*) "margin ", margin

    !if(abs(margin - 1) > 1.0e-6) then

    !    write(*,*) "Ellipsoid still intersecting wall", PART(I,IG)%XP, PART(I,IG)%YP, PART(I,IG)%ZP
    !    write(*,*) "margin ", margin

    !    stop

    !end if

    !write(*,*) "Position before collision",  I,  PART(I,IG)%XP,  PART(I,IG)%YP,  PART(I,IG)%ZP
    !write(*,*) "Orientation :",  PART(I,IG)%ELLQUAT 
    !write(*,*) "Wall at (LX,0,0) with normal (-1,0,0)"

    ! calculate impulse
    call wall_ellipsoid_rebound(PART(I,IG)%UP, PART(I,IG)%VP, PART(I,IG)%WP,             &
                                PART(I,IG)%OMEGAX, PART(I,IG)%OMEGAY, PART(I,IG)%OMEGAZ, &
                                PART(I,IG)%ELLQUAT,impact_arm,wall_normal,               &
                                IG, IDP, NCYCLE, flag)

!    if(flag) then

        ! Update the particle position
!        delta_x = (ellip_pt(1,1) - wall_pos(1,1))

!        x_wall =  PART(I,IG)%XP - delta_x

!    end if

    !write(*,*) " "
    
end if
!!=================================================!!

end subroutine WALL_BOUNDARY_CONDITION