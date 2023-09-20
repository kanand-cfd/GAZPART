!!====================================================================
!!
!!          Particle TRACKING
!!
!!====================================================================

subroutine ADV_PARTICLE(NCYCLE)

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
real(kind=8) :: INVTAUP, REP, FCORR
real(kind=8) :: VRNRM, VRX, VRY, VRZ

real(kind=8), dimension(ndim,1) :: K1_XP, K1_UP, K2_XP, K2_UP
real(kind=8), dimension(ndim,1) :: FFLUID, FDRAG, FLIFT

real(kind=8) :: XP1, YP1, ZP1
real(kind=8) :: UP1, VP1, WP1

real(kind=8) :: XP2, YP2, ZP2
real(kind=8) :: UP2, VP2, WP2

!!- dummy variable for subroutine argument
integer      :: IDUMMY

!- Time control variable
real(kind=8) :: TIME_START, TIME_END

!!- Index
integer :: I, J ,K, NP, IDP
!!====================================================================

!!- CPU check
if(MYID == 0) then
 TIME_START=MPI_WTIME()
end if




!!====================================================================
!! 1. First RK2 Step
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

        !!--------------------------------------------------------------------
        !! Fluid Vorticity Interpolation
        !!--------------------------------------------------------------------
!        call INTERPH(INTERP_SCHEME,                &
!                     XMESH, YMESH, ZMESH,          &
!                     VORTICITY_X,                  &
!                     NPART_LOC(J),                 &
!                     PART(1:NPART_LOC(J),J)%XP,    &
!                     PART(1:NPART_LOC(J),J)%YP,    &
!                     PART(1:NPART_LOC(J),J)%ZP,    &
!                     PART(1:NPART_LOC(J),J)%VORTXAP)
!
!
!       call INTERPH(INTERP_SCHEME,                &
!                     XMESH, YMESH, ZMESH,          &
!                     VORTICITY_Y,                  &
!                     NPART_LOC(J),                 &
!                     PART(1:NPART_LOC(J),J)%XP,    &
!                     PART(1:NPART_LOC(J),J)%YP,    &
!                     PART(1:NPART_LOC(J),J)%ZP,    &
!                     PART(1:NPART_LOC(J),J)%VORTYAP)
!
!
!        call INTERPH(INTERP_SCHEME,                &
!                     XMESH, YMESH, ZMESH,          &
!                     VORTICITY_Z,                  &
!                     NPART_LOC(J),                 &
!                     PART(1:NPART_LOC(J),J)%XP,    &
!                     PART(1:NPART_LOC(J),J)%YP,    &
!                     PART(1:NPART_LOC(J),J)%ZP,    &
!                     PART(1:NPART_LOC(J),J)%VORTZAP)

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

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!         To be done for ellipsoids            !!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do I = 1, NPART_LOC(J)

                IDP = PART(I,J)%IDP

                !! RK2 - 1st Step
                K1_XP(1,1) = PART(I,J)%UP
                K1_XP(2,1) = PART(I,J)%VP
                K1_XP(3,1) = PART(I,J)%WP


                call FLUID_FORCE_ELLIPSOID(J, IDP, &
                                           PART(I,J)%UP, PART(I,J)%VP, PART(I,J)%WP, &
                                           PART(I,J)%UFAP, PART(I,J)%VFAP, PART(I,J)%WFAP, &
                                           PART(I,J)%ELLQUAT, &
                                           FDRAG, FLIFT)


                
                FFLUID = FDRAG + FLIFT

                K1_UP(1,1) = FFLUID(1,1)
                K1_UP(2,1) = FFLUID(2,1)
                K1_UP(3,1) = FFLUID(3,1) + GRAVITY(J)

                PART(I,J)%XP = PART(I,J)%XP + 0.5*DTIME*K1_XP(1,1)
                PART(I,J)%YP = PART(I,J)%YP + 0.5*DTIME*K1_XP(2,1)
                PART(I,J)%ZP = PART(I,J)%ZP + 0.5*DTIME*K1_XP(3,1)


                PART(I,J)%UP = PART(I,J)%UP + 0.5*DTIME*K1_UP(1,1)
                PART(I,J)%VP = PART(I,J)%VP + 0.5*DTIME*K1_UP(2,1)
                PART(I,J)%WP = PART(I,J)%WP + 0.5*DTIME*K1_UP(3,1)


                ! Check for wall boundary condition
                if(WALL_BOUNDARY) call WALL_BOUNDARY_CONDITION(NCYCLE, I, J, IDP)

                !! RK2 - 2nd Step

                K2_XP(1,1) = PART(I,J)%UP
                K2_XP(2,1) = PART(I,J)%UP
                K2_XP(3,1) = PART(I,J)%UP

                call FLUID_FORCE_ELLIPSOID(J, IDP, &
                                           PART(I,J)%UP, PART(I,J)%VP, PART(I,J)%WP, &
                                           PART(I,J)%UFAP, PART(I,J)%VFAP, PART(I,J)%WFAP, &
                                           PART(I,J)%ELLQUAT, &
                                           FDRAG, FLIFT)

                FFLUID = FDRAG + FLIFT


                K2_UP(1,1) = FFLUID(1,1)
                K2_UP(2,1) = FFLUID(2,1)
                K2_UP(3,1) = FFLUID(3,1) + GRAVITY(J)

                PART(I,J)%XP = PART(I,J)%XP + 0.5*DTIME*K2_XP(1,1)
                PART(I,J)%YP = PART(I,J)%YP + 0.5*DTIME*K2_XP(2,1)
                PART(I,J)%ZP = PART(I,J)%ZP + 0.5*DTIME*K2_XP(3,1)


                PART(I,J)%UP = PART(I,J)%UP + 0.5*DTIME*K2_UP(1,1)
                PART(I,J)%VP = PART(I,J)%VP + 0.5*DTIME*K2_UP(2,1)
                PART(I,J)%WP = PART(I,J)%WP + 0.5*DTIME*K2_UP(3,1)


                !Quaternion integration
                call QUATERNION_INTEGRATION(PART(I,J)%ELLQUAT, PART(I,J)%OMEGAX, PART(I,J)%OMEGAY, PART(I,J)%OMEGAZ)

                ! Euler integration for angular velocity
                call EULER_INTEGRATION_FLUID(I, J, IDP)

                !== The scheme below is unstable ==!
                !call QUATERNION_ANGULAR_VELOCITY_INTGERATION(I, J)

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
                call QUATERNION_INTEGRATION(PART(I,J)%ELLQUAT, PART(I,J)%OMEGAX, PART(I,J)%OMEGAY, PART(I,J)%OMEGAZ)

                ! Euler integration for angular velocity
                call EULER_INTEGRATION(J, IDP, PART(I,J)%OMEGAX, PART(I,J)%OMEGAY, PART(I,J)%OMEGAZ)

            end do

        end if 

    end if !!- end if (PARTDEF(J)==2)

end do !!- end loop: do J = 1, NIG


!!====================================================================
!! 4. Particle boundary conditions
!!====================================================================
if(FROZEN_FLOW < 3 .and. PERIODICITY) call BOUNDARY_PARTICLE
if(FROZEN_FLOW == 3) call BOUNDARY_PARTICLE_Y_Z

!! call periodic boundary 


!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()
 CPU_PART(1) = CPU_PART(1) + TIME_END - TIME_START
end if




!!--------------------------------------------------------------------
10000 format (30(e17.7))
10703 format (2x,' Part. in Class ',i3,' :  Max=',i7,' Min=',i7)

end subroutine ADV_PARTICLE

