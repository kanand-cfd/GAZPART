!!====================================================================
!!
!!          Particle TRACKING
!!
!!====================================================================

subroutine ADV_PARTICLE_ELLP(NCYCLE)

!!====================================================================
!!
!! The particle position and velocity equations are time-advanced
!! by a Runge-Kutta scheme.
!!
!! 
!!    dup,i(t)                
!!   --------- = - a.up,i(t) + b 
!!       dt                   
!! where
!!           1              ufa@p,i(t)
!!     a = -----       b = ------------  + gi
!!         tau_p               tau_p      
!!--------------------------------------------------------------------
!! Algorithm
!!----------
!!
!!  1. First RK2 step
!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!
!!  1.1. Fluid velocity interpolation at particle position xp(t)
!!  -------------------------------------------------------------
!!      uf@p[t(n)] = interp{uf[x[t(n)],t(n)]}
!!
!!  1.2. First slope
!!  ----------------
!!               upi(t)
!!      k1i = - ------ + bi(t)
!!              tau_p
!!
!!  1.3. Approximation of particle position and velocity
!!  ----------------------------------------------------
!!
!!      ~upi = upi(t) + k1i*dt
!!
!!      ~xpi = xpi(t) + ~upi*dt
!!
!!
!!  2. Second RK2 step
!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!
!!  2.1. Fluid velocity interpolation at particle position ~xp
!!  ----------------------------------------------------------
!!      ~uf@p = interp{uf[~xp]}
!!
!!  2.2. Second slope
!!  ----------------
!!                ~upi
!!      k2i = - ------ + ~bi
!!              ~tau_p
!!
!!  3. Finalizing RK2 predictions
!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!
!!
!!  4. Apply boundary conditions
!! +++++++++++++++++++++++++++++
!!
!!====================================================================
!! Warning !!
!!!!!!!!!!!!!
!! In the present version the fluid velocity at the particle
!! position is kept constant. It can be imporved later!
!!
!!====================================================================

use STATISTICS 
use PARTICLE_PARALLEL
use MPI_STRUCTURES
use FLUID_VARIABLE
use DNS_DIM               
use PARAM_PHYS
use CHECK_CPU
use mod_quaternion

implicit none

!!====================================================================
!! ARRAYS
!!====================================================================
!!- cycle number
integer, intent(in) :: NCYCLE


!!
!real(kind=8) :: INVTAUP, REP, FCORR
real(kind=8) :: VRNRM !, VRX, VRY, VRZ

real(kind=8), dimension(ndim,1) :: K1_XP, K1_UP, K2_XP, K2_UP
real(kind=8), dimension(ndim,1) :: K1_OMG, K2_OMG

real(kind=8), dimension(ndim,1) :: FFLUID, FDRAG, FLIFT
real(kind=8), dimension(ndim, 1) :: FTORQUE, FTORQUE_ROT, FTORQUE_PITCH

real(kind=8), dimension(ndim,1) :: impact_arm, wall_normal, wall_pos, plane_pt

real(kind=8), dimension(ndim, 1) :: U_P, OMG_lab, VPC

real(kind=8) :: XP1, YP1, ZP1
real(kind=8) :: UP1, VP1, WP1
real(kind=8) :: OMGX1, OMGY1, OMGZ1 

real(kind=8) :: x_wall, delta_x, delt

real(kind=8) :: XP2, YP2, ZP2, margin

type(quaternion) :: QUAT1, QUAT2

!!- dummy variable for subroutine argument
!integer      :: IDUMMY

!- Time control variable
real(kind=8) :: TIME_START, TIME_END

!!- Index
integer :: I, J, IDP ! ,K, NP
logical :: flag, isIntersect
!!====================================================================

!!- CPU check
if(MYID == 0) then
 TIME_START=MPI_WTIME()
end if




!!====================================================================
!! 1. First Step
!!====================================================================
do J = 1, NIG


    if(SOLVE_FLUID>0 .or. FROZEN_FLOW>1) then 

        !!--------------------------------------------------------------------
        !! 1.1 Fluid velocity Interpolation
        !!--------------------------------------------------------------------
        call INTERPH(INTERP_SCHEME,              &
                     XMESH, YMESH, ZMESH,        &
                     UFLU,                       &
                     NPART_LOC(J),               &
                     PART(1:NPART_LOC(J),J)%XP,  &
                     PART(1:NPART_LOC(J),J)%YP,  &
                     PART(1:NPART_LOC(J),J)%ZP,  &
                     PART(1:NPART_LOC(J),J)%UFAP )


        call INTERPH(INTERP_SCHEME,              &
                     XMESH, YMESH, ZMESH,        &
                     VFLU,                       &
                     NPART_LOC(J),               &
                     PART(1:NPART_LOC(J),J)%XP,  &
                     PART(1:NPART_LOC(J),J)%YP,  &
                     PART(1:NPART_LOC(J),J)%ZP,  &
                     PART(1:NPART_LOC(J),J)%VFAP )


        call INTERPH(INTERP_SCHEME,              &
                     XMESH, YMESH, ZMESH,        &
                     WFLU,                       &
                     NPART_LOC(J),               &
                     PART(1:NPART_LOC(J),J)%XP,  &
                     PART(1:NPART_LOC(J),J)%YP,  &
                     PART(1:NPART_LOC(J),J)%ZP,  &
                     PART(1:NPART_LOC(J),J)%WFAP )

        !!--------------------------------------------------------------------
        !! 1.2 Fluid Vorticity Interpolation
        !!--------------------------------------------------------------------
        call INTERPH(INTERP_SCHEME,                &
                     XMESH, YMESH, ZMESH,          &
                     VORTICITY_X,                  &
                     NPART_LOC(J),                 &
                     PART(1:NPART_LOC(J),J)%XP,    &
                     PART(1:NPART_LOC(J),J)%YP,    &
                     PART(1:NPART_LOC(J),J)%ZP,    &
                     PART(1:NPART_LOC(J),J)%VORTXAP)


        call INTERPH(INTERP_SCHEME,                &
                     XMESH, YMESH, ZMESH,          &
                     VORTICITY_Y,                  &
                     NPART_LOC(J),                 &
                     PART(1:NPART_LOC(J),J)%XP,    &
                     PART(1:NPART_LOC(J),J)%YP,    &
                     PART(1:NPART_LOC(J),J)%ZP,    &
                     PART(1:NPART_LOC(J),J)%VORTYAP)


        call INTERPH(INTERP_SCHEME,                &
                     XMESH, YMESH, ZMESH,          &
                     VORTICITY_Z,                  &
                     NPART_LOC(J),                 &
                     PART(1:NPART_LOC(J),J)%XP,    &
                     PART(1:NPART_LOC(J),J)%YP,    &
                     PART(1:NPART_LOC(J),J)%ZP,    &
                     PART(1:NPART_LOC(J),J)%VORTZAP)


    else if(FROZEN_FLOW == 1) then

        PART(1:NPART_LOC(J),J)%UFAP = 0.0
        PART(1:NPART_LOC(J),J)%VFAP = 0.0
        PART(1:NPART_LOC(J),J)%WFAP = 0.0

        PART(1:NPART_LOC(J),J)%VORTXAP = 0.0
        PART(1:NPART_LOC(J),J)%VORTYAP = 0.0
        PART(1:NPART_LOC(J),J)%VORTZAP = 0.0

    else

        PART(1:NPART_LOC(J),J)%UFAP = ZERO
        PART(1:NPART_LOC(J),J)%VFAP = ZERO
        PART(1:NPART_LOC(J),J)%WFAP = ZERO


    end if !- end if SOLVE_FLUID>1

end do !- end loop  J = 1, NIG


!!====================================================================
!! 2. Advance particle velocities
!!====================================================================

!!--------------------------------------------------------------------
!! 2.1 Save variables
!!--------------------------------------------------------------------
do J = 1, NIG


    !!- Motionless particles
    if(PARTDEF(J)==0) then

        PART(1:NPART_LOC(J),J)%UP = PART(1:NPART_LOC(J),J)%UFAP
        PART(1:NPART_LOC(J),J)%VP = PART(1:NPART_LOC(J),J)%VFAP
        PART(1:NPART_LOC(J),J)%WP = PART(1:NPART_LOC(J),J)%WFAP


    elseif(PARTDEF(J)==1) then

        PART(1:NPART_LOC(J),J)%UP = PART(1:NPART_LOC(J),J)%UFAP
        PART(1:NPART_LOC(J),J)%VP = PART(1:NPART_LOC(J),J)%VFAP
        PART(1:NPART_LOC(J),J)%WP = PART(1:NPART_LOC(J),J)%WFAP

        do I = 1, NPART_LOC(J)
            
            PART(I,J)%XP = PART(I,J)%XP + DTIME*PART(I,J)%UP 
            PART(I,J)%YP = PART(I,J)%YP + DTIME*PART(I,J)%VP
            PART(I,J)%ZP = PART(I,J)%ZP + DTIME*PART(I,J)%WP

        end do


    elseif(PARTDEF(J)==2) then

        if(SOLVE_FLUID>0 .or. (FROZEN_FLOW>0)) then

            do I = 1, NPART_LOC(J)

                !Quaternion integration
                !call QUATERNION_INTEGRATION(PART(I,J)%ELLQUAT, PART(I,J)%OMEGAX, PART(I,J)%OMEGAY, PART(I,J)%OMEGAZ, DTIME)

                IDP = PART(I,J)%IDP
                flag = .false.

                !! RK2 =========================== 1st Step =========================== RK2 !!
                !!==========================================================================!!
                K1_XP(1,1) = PART(I,J)%UP
                K1_XP(2,1) = PART(I,J)%VP
                K1_XP(3,1) = PART(I,J)%WP


                call FLUID_FORCE_ELLIPSOID(J, IDP, &
                                           PART(I,J)%UP, PART(I,J)%VP, PART(I,J)%WP, &
                                           PART(I,J)%UFAP, PART(I,J)%VFAP, PART(I,J)%WFAP, &
                                           PART(I,J)%ELLQUAT, &
                                           FDRAG, FLIFT)

                FFLUID = FDRAG + FLIFT

                call FLUID_TORQUE_ELLIPSOID(I, J, IDP,        &
                                           PART(I,J)%OMEGAX,  &
                                           PART(I,J)%OMEGAY,  &
                                           PART(I,J)%OMEGAZ,  &
                                           PART(I,J)%ELLQUAT, & 
                                           FTORQUE_PITCH,     &
                                           FTORQUE_ROT) 

                FTORQUE = FTORQUE_PITCH + FTORQUE_ROT

                K1_UP(1,1) = FFLUID(1,1)
                K1_UP(2,1) = FFLUID(2,1)
                K1_UP(3,1) = FFLUID(3,1) + GRAVITY(J)

                K1_OMG(1,1) = (FTORQUE(1,1) + PART(I,J)%OMEGAY*PART(I,J)%OMEGAZ*(IPYY(J,IDP) - IPZZ(J,IDP))) / IPXX(J,IDP) ! omega_x
                K1_OMG(2,1) = (FTORQUE(2,1) + PART(I,J)%OMEGAZ*PART(I,J)%OMEGAX*(IPZZ(J,IDP) - IPXX(J,IDP))) / IPYY(J,IDP) ! omega_y
                K1_OMG(3,1) = (FTORQUE(3,1) + PART(I,J)%OMEGAX*PART(I,J)%OMEGAY*(IPXX(J,IDP) - IPYY(J,IDP))) / IPZZ(J,IDP) ! omega_z


                XP1 = PART(I,J)%XP + DTIME*K1_XP(1,1)
                YP1 = PART(I,J)%YP + DTIME*K1_XP(2,1)
                ZP1 = PART(I,J)%ZP + DTIME*K1_XP(3,1)

                UP1 = PART(I,J)%UP + DTIME*K1_UP(1,1)
                VP1 = PART(I,J)%VP + DTIME*K1_UP(2,1)
                WP1 = PART(I,J)%WP + DTIME*K1_UP(3,1)

                OMGX1 = PART(I,J)%OMEGAX + DTIME*K1_OMG(1,1)
                OMGY1 = PART(I,J)%OMEGAY + DTIME*K1_OMG(2,1)
                OMGZ1 = PART(I,J)%OMEGAZ + DTIME*K1_OMG(3,1)

                !!!!!!! Integrate Quaternion !!!!!!
                QUAT1 = PART(I,J)%ELLQUAT

                call QUATERNION_INTEGRATION(QUAT1, OMGX1, OMGY1, OMGZ1, DTIME)

                !! =========================== End of 1st Step ======================== RK2 !!
                !!==========================================================================!!

                ! Check for wall boundary condition
                call WALL_BOUNDARY_CONDITION_NEW(NCYCLE, I, J, IDP,      &
                                                XP1, YP1, ZP1,           &
                                                QUAT1,                   &
                                                wall_normal, impact_arm, &
                                                delta_x,                 &
                                                isIntersect)

                U_P(1,1) = UP1; U_P(2,1) = VP1; U_P(3,1) = WP1
                OMG_lab(1,1) = OMGX1; OMG_lab(2,1) = OMGY1; OMG_lab(3,1) = OMGZ1

                ! Transform to global frame
                call transform_basis(OMG_lab, QUAT1, shape(OMG_lab))

                ! Velocity of contact point 
                VPC(1,1) = U_P(1,1) + (OMG_lab(2,1)*impact_arm(3,1) - OMG_lab(3,1)*impact_arm(2,1))
                VPC(2,1) = U_P(2,1) + (OMG_lab(3,1)*impact_arm(1,1) - OMG_lab(1,1)*impact_arm(3,1))
                VPC(3,1) = U_P(3,1) + (OMG_lab(1,1)*impact_arm(2,1) - OMG_lab(2,1)*impact_arm(1,1))

                VRNRM = VPC(1,1)*wall_normal(1,1) + VPC(2,1)*wall_normal(2,1) + VPC(3,1)*wall_normal(3,1)

                x_wall =  XP1 + delta_x
                    
                if(isIntersect) XP1 = x_wall

                !!==========================================================================!!
                !!              If Intersection & Particle is moving towards wall           !!
                !!==========================================================================!!
                if(isIntersect .and. VRNRM<0.0) then

                    !delt = abs(delta_x)/abs(VPC(1,1)) !DTIME * (abs(PART(I,J)%XP - x_wall)/abs(PART(I,J)%XP - XP1)) !abs(delta_x)/abs(PART(I,J)%UP)

                    !write(*,*) ' Before Update ', XP1, YP1, ZP1, UP1, VP1, WP1
                    !Update the particle position
                    !YP1 = PART(I,J)%YP + delt*K1_XP(2,1)
                    !ZP1 = PART(I,J)%ZP + delt*K1_XP(3,1)

                    !UP1 = PART(I,J)%UP + delt*K1_UP(1,1)
                    !VP1 = PART(I,J)%VP + delt*K1_UP(2,1)
                    !WP1 = PART(I,J)%WP + delt*K1_UP(3,1)

                    !OMGX1 = PART(I,J)%OMEGAX + delt*K1_OMG(1,1)
                    !OMGY1 = PART(I,J)%OMEGAY + delt*K1_OMG(2,1) 
                    !OMGZ1 = PART(I,J)%OMEGAZ + delt*K1_OMG(3,1)

                    !!=================================================================================!!
                    !! RK2 ========================= Particle Wall Bouncing ====================== RK2 !!
                    !!=================================================================================!!

                    !write(*,*) 'Before Impulse ', UP1, VP1, WP1
                    call WALL_ELLIPSOID_REBOUND(UP1, VP1, WP1,          &
                                               OMGX1, OMGY1, OMGZ1,     &
                                               QUAT1,                   &
                                               impact_arm, wall_normal, &
                                               J, IDP, NCYCLE,          &
                                               flag)
                    !write(*,*) 'After Impulse ', UP1, VP1, WP1

                    !!=================================================================================!!
                    !! RK2 ==================== 2nd Step for Bounced Particle ==================== RK2 !!
                    !!=================================================================================!!
                    K2_XP(1,1) = UP1
                    K2_XP(2,1) = VP1
                    K2_XP(3,1) = WP1

                    call FLUID_FORCE_ELLIPSOID(J, IDP, &
                                               UP1, VP1, WP1, &
                                               PART(I,J)%UFAP, PART(I,J)%VFAP, PART(I,J)%WFAP, &
                                               QUAT1, &
                                               FDRAG, FLIFT)

                    FFLUID = FDRAG + FLIFT

                    call FLUID_TORQUE_ELLIPSOID(I, J, IDP,           &
                                                OMGX1, OMGY1, OMGZ1, &
                                                QUAT1,               &
                                                FTORQUE_PITCH, FTORQUE_ROT) 

                    FTORQUE = FTORQUE_PITCH + FTORQUE_ROT


                    K2_UP(1,1) = FFLUID(1,1)
                    K2_UP(2,1) = FFLUID(2,1)
                    K2_UP(3,1) = FFLUID(3,1) + GRAVITY(J)

                    K2_OMG(1,1) = (FTORQUE(1,1) + OMGY1*OMGZ1*(IPYY(J,IDP) - IPZZ(J,IDP))) / IPXX(J,IDP) ! omega_x
                    K2_OMG(2,1) = (FTORQUE(2,1) + OMGZ1*OMGX1*(IPZZ(J,IDP) - IPXX(J,IDP))) / IPYY(J,IDP) ! omega_y
                    K2_OMG(3,1) = (FTORQUE(3,1) + OMGX1*OMGY1*(IPXX(J,IDP) - IPYY(J,IDP))) / IPZZ(J,IDP) ! omega_z

                    PART(I,J)%XP = XP1 + (DTIME)*K2_XP(1,1)
                    PART(I,J)%YP = YP1 + (DTIME)*K2_XP(2,1)
                    PART(I,J)%ZP = ZP1 + (DTIME)*K2_XP(3,1)

                    PART(I,J)%UP = UP1 + (DTIME)*K2_UP(1,1)
                    PART(I,J)%VP = VP1 + (DTIME)*K2_UP(2,1)
                    PART(I,J)%WP = WP1 + (DTIME)*K2_UP(3,1)

                    PART(I,J)%OMEGAX = OMGX1 + (DTIME)*K2_OMG(1,1)
                    PART(I,J)%OMEGAY = OMGY1 + (DTIME)*K2_OMG(2,1)
                    PART(I,J)%OMEGAZ = OMGZ1 + (DTIME)*K2_OMG(3,1)

                    !!!!!!! Integrate Quaternion !!!!!!
                    PART(I,J)%ELLQUAT = QUAT1

                    call QUATERNION_INTEGRATION(PART(I,J)%ELLQUAT, PART(I,J)%OMEGAX, PART(I,J)%OMEGAY, PART(I,J)%OMEGAZ, DTIME)


                    !if(wall_normal(1,1) == 1.0) then

                    !    wall_pos(1,1) = 0.0
                    !    wall_pos(2,1) = PART(I,J)%YP
                    !    wall_pos(3,1) = PART(I,J)%ZP

                    !else

                    !    wall_pos(1,1) = LXMAX
                    !    wall_pos(2,1) = PART(I,J)%YP
                    !    wall_pos(3,1) = PART(I,J)%ZP

                    !end if

                    !call closest_point_plane(PART(I,J)%XP, PART(I,J)%YP, PART(I,J)%ZP, &
                    !                        PART(I,J)%ELLQUAT,                         &
                    !                        wall_pos,wall_normal,                      &
                    !                        plane_pt,margin, J, IDP) 

                    !if((1.0 - margin) > 1.0e-3) then 
                    !    write(*,*) ' '
                    !    write(*,*) 'Still intersecting after rebound', NCYCLE, I, margin
                    !    write(*,*) PART(I,J)%XP, PART(I,J)%YP, PART(I,J)%ZP, PART(I,J)%UP, PART(I,J)%VP, PART(I,J)%WP 
                    !    write(*,*) DTIME, delt
                    !    if(wall_normal(1,1)==1.0) write(*,*)  'Wall at X=0'
                    !    if(wall_normal(1,1)==-1.0) write(*,*) 'Wall at X=LXMAX'
                        !stop

                    !end if 

                    !!=================================================================================!!
                    !! RK2 ==================== 2nd Step for Bounced Particle ==================== RK2 !!
                    !!=================================================================================!!

                else
                    !!=================================================================================!!
                    !! RK2 =================== 2nd Step for Unbounced Particle =================== RK2 !!
                    !!=================================================================================!!

                    K2_XP(1,1) = UP1
                    K2_XP(2,1) = VP1
                    K2_XP(3,1) = WP1

                    call FLUID_FORCE_ELLIPSOID(J, IDP, &
                                               UP1, VP1, WP1, &
                                               PART(I,J)%UFAP, PART(I,J)%VFAP, PART(I,J)%WFAP, &
                                               QUAT1, &
                                               FDRAG, FLIFT)

                    FFLUID = FDRAG + FLIFT

                    call FLUID_TORQUE_ELLIPSOID(I, J, IDP,           &
                                                OMGX1, OMGY1, OMGZ1, &
                                                QUAT1,               &
                                                FTORQUE_PITCH, FTORQUE_ROT) 

                    FTORQUE = FTORQUE_PITCH + FTORQUE_ROT


                    K2_UP(1,1) = FFLUID(1,1)
                    K2_UP(2,1) = FFLUID(2,1)
                    K2_UP(3,1) = FFLUID(3,1) + GRAVITY(J)

                    K2_OMG(1,1) = (FTORQUE(1,1) + OMGY1*OMGZ1*(IPYY(J,IDP) - IPZZ(J,IDP))) / IPXX(J,IDP) ! omega_x
                    K2_OMG(2,1) = (FTORQUE(2,1) + OMGZ1*OMGX1*(IPZZ(J,IDP) - IPXX(J,IDP))) / IPYY(J,IDP) ! omega_y
                    K2_OMG(3,1) = (FTORQUE(3,1) + OMGX1*OMGY1*(IPXX(J,IDP) - IPYY(J,IDP))) / IPZZ(J,IDP) ! omega_z

                    !!==========================================================================!!
                    !!          Check if Particle is outside the wall after second step         !!
                    !!==========================================================================!!

                    XP2 = PART(I,J)%XP + 0.5*DTIME*(K1_XP(1,1) + K2_XP(1,1))
                    YP2 = PART(I,J)%YP + 0.5*DTIME*(K1_XP(2,1) + K2_XP(2,1))
                    ZP2 = PART(I,J)%ZP + 0.5*DTIME*(K1_XP(3,1) + K2_XP(3,1))


                    call WALL_BOUNDARY_CONDITION_NEW(NCYCLE, I, J, IDP,      &
                                                    XP2, YP2, ZP2,           &
                                                    QUAT1,                   &
                                                    wall_normal, impact_arm, &
                                                    delta_x,                 &
                                                    isIntersect)


                    ! If the particle is still intersecting the wall
                    
                    if(isIntersect) then 

                        PART(I,J)%XP = PART(I,J)%XP + DTIME*K1_XP(1,1)
                        PART(I,J)%YP = PART(I,J)%YP + DTIME*K1_XP(2,1)
                        PART(I,J)%ZP = PART(I,J)%ZP + DTIME*K1_XP(3,1)

                        PART(I,J)%UP = PART(I,J)%UP + DTIME*K1_UP(1,1)
                        PART(I,J)%VP = PART(I,J)%VP + DTIME*K1_UP(2,1)
                        PART(I,J)%WP = PART(I,J)%WP + DTIME*K1_UP(3,1)

                        PART(I,J)%OMEGAX = PART(I,J)%OMEGAX + DTIME*K1_OMG(1,1)
                        PART(I,J)%OMEGAY = PART(I,J)%OMEGAY + DTIME*K1_OMG(2,1)
                        PART(I,J)%OMEGAZ = PART(I,J)%OMEGAZ + DTIME*K1_OMG(3,1)

                        !!!!!!! Integrate Quaternion !!!!!!
                        PART(I,J)%ELLQUAT = QUAT1

                        call QUATERNION_INTEGRATION(PART(I,J)%ELLQUAT, PART(I,J)%OMEGAX, PART(I,J)%OMEGAY, PART(I,J)%OMEGAZ, DTIME)

                        !if(wall_normal(1,1) == 1.0) then

                        !    wall_pos(1,1) = 0.0
                        !    wall_pos(2,1) = PART(I,J)%YP
                        !    wall_pos(3,1) = PART(I,J)%ZP

                        !else

                        !    wall_pos(1,1) = LXMAX
                        !    wall_pos(2,1) = PART(I,J)%YP
                        !    wall_pos(3,1) = PART(I,J)%ZP

                        !end if

                        !call closest_point_plane(PART(I,J)%XP, PART(I,J)%YP, PART(I,J)%ZP, &
                        !                        PART(I,J)%ELLQUAT,                         &
                        !                        wall_pos,wall_normal,                      &
                        !                        plane_pt,margin, J, IDP) 

                        !if((1.0 - margin) > 1.0e-6) then 
                        !    write(*,*) 'Still intersecting after not having a rebound ', NCYCLE, I, margin
                        !    write(*,*) 'Position ', PART(I,J)%XP, PART(I,J)%YP, PART(I,J)%ZP
                        !    write(*,*) 'Velocity ', PART(I,J)%UP, PART(I,J)%VP, PART(I,J)%WP
                        !    write(*,*) VRNRM
                        !    if(wall_normal(1,1)==1.0) write(*,*)  'Wall at X=0'
                        !    if(wall_normal(1,1)==-1.0) write(*,*) 'Wall at X=LXMAX'
                        !    !stop
                        !end if 

                    ! If the particle is inside the domain
                    else

                        PART(I,J)%XP = PART(I,J)%XP + 0.5*DTIME*(K1_XP(1,1) + K2_XP(1,1))
                        PART(I,J)%YP = PART(I,J)%YP + 0.5*DTIME*(K1_XP(2,1) + K2_XP(2,1))
                        PART(I,J)%ZP = PART(I,J)%ZP + 0.5*DTIME*(K1_XP(3,1) + K2_XP(3,1))

                        PART(I,J)%UP = PART(I,J)%UP + 0.5*DTIME*(K1_UP(1,1) + K2_UP(1,1))
                        PART(I,J)%VP = PART(I,J)%VP + 0.5*DTIME*(K1_UP(2,1) + K2_UP(2,1))
                        PART(I,J)%WP = PART(I,J)%WP + 0.5*DTIME*(K1_UP(3,1) + K2_UP(3,1))

                        PART(I,J)%OMEGAX = PART(I,J)%OMEGAX + 0.5*DTIME*(K1_OMG(1,1) + K2_OMG(1,1))
                        PART(I,J)%OMEGAY = PART(I,J)%OMEGAY + 0.5*DTIME*(K1_OMG(2,1) + K2_OMG(2,1))
                        PART(I,J)%OMEGAZ = PART(I,J)%OMEGAZ + 0.5*DTIME*(K1_OMG(3,1) + K2_OMG(3,1))


                        !!!!!!! Integrate Quaternion !!!!!!
                        PART(I,J)%ELLQUAT = QUAT1

                        call QUATERNION_INTEGRATION(PART(I,J)%ELLQUAT, PART(I,J)%OMEGAX, PART(I,J)%OMEGAY, PART(I,J)%OMEGAZ, DTIME)


                    end if

                end if

                !!=================================================================================!!
                !! RK2 =========================== End of 2nd Step =========================== RK2 !!
                !!=================================================================================!!
            

            end do


        else 

            !if(MYID==0) write(*,*) "Dry Granular Flow"

            do I = 1, NPART_LOC(J)

                IDP = PART(I,J)%IDP

                K1_XP(1,1) = DTIME*PART(I,J)%UP
                K1_XP(2,1) = DTIME*PART(I,J)%VP
                K1_XP(3,1) = DTIME*PART(I,J)%WP

                K2_XP(1,1) = DTIME*(PART(I,J)%UP + K1_XP(1,1))
                K2_XP(2,1) = DTIME*(PART(I,J)%VP + K1_XP(2,1))
                K2_XP(3,1) = DTIME*(PART(I,J)%WP + K1_XP(3,1))

                PART(I,J)%XP = PART(I,J)%XP + K2_XP(1,1)
                PART(I,J)%YP = PART(I,J)%YP + K1_XP(2,1)
                PART(I,J)%ZP = PART(I,J)%ZP + K2_XP(3,1)

                !Quaternion integration
                call QUATERNION_INTEGRATION(PART(I,J)%ELLQUAT, PART(I,J)%OMEGAX, PART(I,J)%OMEGAY, PART(I,J)%OMEGAZ, DTIME)

                ! Euler integration for angular velocity
                call EULER_INTEGRATION(J, IDP, PART(I,J)%OMEGAX, PART(I,J)%OMEGAY, PART(I,J)%OMEGAZ)

            end do

        end if 

    end if !!- end if (PARTDEF(J)==2)

end do !!- end loop: do J = 1, NIG


!!====================================================================
!! 4. Particle boundary conditions
!!====================================================================
if(FROZEN_FLOW == 3) call BOUNDARY_PARTICLE_Y_Z
if(FROZEN_FLOW < 3 .and. PERIODICITY) call BOUNDARY_PARTICLE

!! call periodic boundary 


!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()
 CPU_PART(1) = CPU_PART(1) + TIME_END - TIME_START
end if




!!--------------------------------------------------------------------
10000 format (30(e17.7))
10703 format (2x,' Part. in Class ',i3,' :  Max=',i7,' Min=',i7)

end subroutine ADV_PARTICLE_ELLP

