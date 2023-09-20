!!===================================================================
!!
!! This subroutine calculates all the statistics along the x-span of 
!! the channel.
!!
!!===================================================================
subroutine STAT_PART_CHANNEL(NCYCLE, TIME)
!!-------------------------------------------------------------------
!! Statistics (temporaly and spatially averaged) are performed using 
!! macro array called MEAN_TIME_PART_CHAN. The averaging is done as
!! --> MEAN_TIME_PART_CHAN = MEAN_TIME_PART_CHAN + MEAN_PART_CHAN
!!
!!===================================================================

use DNS_DIM            !- Dimension
use PARAM_PHYS         !- Physical & numerical parameters
use PARTICLE_PARALLEL
use GEOMETRIC_VARIABLE !- Mesh parameters
use STATISTICS         !- Statistics
use CHECK_CPU 
use FLUID_VARIABLE        

implicit none

!---------------------------------------------------------------------
!- Global Variables
!---------------------------------------------------------------------
!- Curent time
real(kind=8), intent(in) :: NCYCLE, TIME
!--------------------------------------------------------------------
!	ARRAYS STATEMENT
!--------------------------------------------------------------------
!- Local arrays
!---------------------------------------------------------------------
!- Arrays containing statistics
real(kind=8), dimension(NSTAT,NIG,POLYDISPMAX, NSLICE) :: MEAN_PART_CHAN_LOC
real(kind=8), dimension(NSTAT,NIG,POLYDISPMAX, NSLICE) :: MEAN_PART_CHAN

! Fluid velocty at particle position
real(kind=8) :: UF_P, VF_P, WF_P

! Fluid force
real(kind=8), dimension(ndim,1) :: FDRAG, FLIFT
real(kind=8), dimension(ndim, 1) :: FTORQUE_PITCH, FTORQUE_ROT

real(kind=8) :: CDRAG, CLIFT
real(kind=8) :: CDRAG_PHI0, CDRAG_PHI90

real(kind=8) :: CTPITCH
real(kind=8) :: CTROTATION, CTROTATION_2

!!! Fitting Parameters for Drag Coefficient !!!
real(kind=8) :: a0, a1, a2, a3, a4, a5, a6, a7, a8

!!! Fitting Parameters for Lift Coefficient !!!
real(kind=8) :: b1, b2, b3, b4, b5, b6, b7, b8, b9, b10

!! Fitting Parameters for Pitching torque coefficient !!
real(kind=8) :: c1, c2, c3, c4, c5, c6, c7, c8, c9, c10

!! Fitting Parameters for Rotation torque coefficient !!
real(kind=8) :: r1, r2, r3, r4

real(kind=8) :: r1_, r2_, r3_, r4_

real(kind=8), dimension(ndim, 1) :: UREL            ! relative velocity
real(kind=8), dimension(ndim, 1) :: UFLU_PFRAME     ! fluid velocity in particle frame
real(kind=8), dimension(ndim, 1) :: OREL            ! Relative rotational velocity
real(kind=8), dimension(ndim, 1) :: curl_UFLU       ! curl of fluid velocity

real(kind=8) :: phi      ! Angle of Incidence
real(kind=8) :: dequiv   ! Equivalent Diameter
real(kind=8) :: Rep      ! Particle Reynolds number
real(kind=8) :: UREL_mag ! Magnitude of relative velocity
real(kind=8) :: UFLU_mag ! Magnitude of fluid velocity in particle frame
real(kind=8) :: Re_R     ! Rotational Particle Reynolds Number
real(kind=8) :: OREL_MAG ! Magnitude of Relative Rotational Velocity
real(kind=8) :: OREL_YZ
real(kind=8) :: PMASS


!- Index of particle class for polydisperse stat
integer :: IDP

!- Variables for initialising the slice index
real(kind=8) :: DELX, NSM1, NSP1, somme

!- Time control variable
real(kind=8) :: TIME_START, TIME_END

integer :: I,J,K, IFLAG1, IPART, JPART, KPART, NS

!- file number
integer :: NUMFILE

!- Collision flag
logical :: CFLAG_A, CFLAG_B 

!---------------------------------------------------------------------

!!==================================================================================!!
!! Fitting parameters taken from Zastawny et al. (2012) for Ellipsoid 2 (a/b=1.25) !!
!! Drag 
a0 = 1.95
a1 = 18.12; a2 = 1.023; a3 = 4.26; a4 = 0.384
a5 = 21.52; a6 = 0.990; a7 = 2.86; a8 = 0.260
!! Lift
b1 = 0.083; b2 = -0.21; b3 = 1.582; b4 = 0.851; b5  = 1.842
b6 =-0.802; b7 =-0.006; b8 = 0.874; b9 = 0.009; b10 = 0.570

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Fitting parameters taken from Zastawny et al. (2012) for Ellipsoid 2 (a/b=1.25) !!
!! Pitching Torque
c1 = 0.935; c2 = 0.146; c3 = -0.469; c4 = 0.145; c5  = 0.116
c6 = 0.748; c7 = 0.041; c8 =  0.221; c9 = 0.657; c10 = 0.044

!! Rotation Torque (Mode 1)
r1 = 0.573; r2 =-0.154; r3 = 116.61; r4 = 1.0

!! Rotation Torque (Mode 2)
r1_ = 1.244; r2_ = 0.239; r3_ = 378.12; r4_ = 0.789
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!- CPU check
if(MYID == 0) then
 TIME_START=MPI_WTIME()
end if

somme = ZERO
MEAN_PART_CHAN_LOC(:,:,:,:) = ZERO
MEAN_PART_CHAN(:,:,:,:) = ZERO

!--------------------------------------------------------------------
! Fluid Velocity Interpolation
!--------------------------------------------------------------------

! do J = 1, NIG


!     if(SOLVE_FLUID>0 .or. FROZEN_FLOW>0) then 


!         call INTERPH(INTERP_SCHEME,              &
!                      XMESH, YMESH, ZMESH,        &
!                      UFLU,                       &
!                      NPART_LOC(J),               &
!                      PART(1:NPART_LOC(J),J)%XP,  &
!                      PART(1:NPART_LOC(J),J)%YP,  &
!                      PART(1:NPART_LOC(J),J)%ZP,  &
!                      PART(1:NPART_LOC(J),J)%UFAP )

!         call INTERPH(INTERP_SCHEME,              &
!                      XMESH, YMESH, ZMESH,        &
!                      VFLU,                       &
!                      NPART_LOC(J),               &
!                      PART(1:NPART_LOC(J),J)%XP,  &
!                      PART(1:NPART_LOC(J),J)%YP,  &
!                      PART(1:NPART_LOC(J),J)%ZP,  &
!                      PART(1:NPART_LOC(J),J)%VFAP )

!         call INTERPH(INTERP_SCHEME,              &
!                      XMESH, YMESH, ZMESH,        &
!                      WFLU,                       &
!                      NPART_LOC(J),               &
!                      PART(1:NPART_LOC(J),J)%XP,  &
!                      PART(1:NPART_LOC(J),J)%YP,  &
!                      PART(1:NPART_LOC(J),J)%ZP,  &
!                      PART(1:NPART_LOC(J),J)%WFAP )

!     else

!         PART(1:NPART_LOC(J),J)%UFAP = ZERO
!         PART(1:NPART_LOC(J),J)%VFAP = ZERO
!         PART(1:NPART_LOC(J),J)%WFAP = ZERO

!     end if !- end if SOLVE_FLUID>1


! end do !- end loop  J = 1, NIG



do J = 1, NIG

    do I = 1, NPART_LOC(J)

        IDP = PART(I,J)%IDP


        call WALL_INTERP_LAG3(XMESH, YMESH, ZMESH, &
                              UFLU,                &
                              PART(I,J)%XP, PART(I,J)%YP, PART(I,J)%ZP, &
                              UF_P)

        call WALL_INTERP_LAG3(XMESH, YMESH, ZMESH, &
                              VFLU,                &
                              PART(I,J)%XP, PART(I,J)%YP, PART(I,J)%ZP, &
                              VF_P)

        call WALL_INTERP_LAG3(XMESH, YMESH, ZMESH, &
                              WFLU,                &
                              PART(I,J)%XP, PART(I,J)%YP, PART(I,J)%ZP, &
                              WF_P)

!--------------------------------------------------------------------
! If particle lies in the vicinity (+-delx) of slice
!--------------------------------------------------------------------
        call LOCATE_PART(J, PART(I,J), 2.0*EMAJ_PART(J,IDP)/APR_PART(J), NS)

!!--------------------------------------------------------------------
!!- Particle number for each class
!!--------------------------------------------------------------------
!!- Np
        MEAN_PART_CHAN_LOC(44,J,IDP,NS) = MEAN_PART_CHAN_LOC(44,J,IDP,NS) + 1.0      
        
!!--------------------------------------------------------------------
!!- Particle velocities
!!--------------------------------------------------------------------
!!- <up>
        MEAN_PART_CHAN_LOC(1,J,IDP,NS) = MEAN_PART_CHAN_LOC(1,J,IDP,NS) + PART(I,J)%UP

!!- <vp>
        MEAN_PART_CHAN_LOC(2,J,IDP,NS) = MEAN_PART_CHAN_LOC(2,J,IDP,NS) + PART(I,J)%VP

!!- <wp>
        MEAN_PART_CHAN_LOC(3,J,IDP,NS) = MEAN_PART_CHAN_LOC(3,J,IDP,NS) + PART(I,J)%WP

!!--------------------------------------------------------------------
!!- Fluid velocities at particle position
!!--------------------------------------------------------------------
!!- <uf@p>
        MEAN_PART_CHAN_LOC(10,J,IDP,NS) = MEAN_PART_CHAN_LOC(10,J,IDP,NS) + UF_P

!!- <vf@p>
        MEAN_PART_CHAN_LOC(11,J,IDP,NS) = MEAN_PART_CHAN_LOC(11,J,IDP,NS) + VF_P

!!- <wf@p>
        MEAN_PART_CHAN_LOC(12,J,IDP,NS) = MEAN_PART_CHAN_LOC(12,J,IDP,NS) + WF_P

!!--------------------------------------------------------------------
!! Particle kinetic stress
!!--------------------------------------------------------------------
!!- <up*up>
        MEAN_PART_CHAN_LOC(4,J,IDP,NS) = MEAN_PART_CHAN_LOC(4,J,IDP,NS) + PART(I,J)%UP*PART(I,J)%UP

!!- <vp*vp>
        MEAN_PART_CHAN_LOC(5,J,IDP,NS) = MEAN_PART_CHAN_LOC(5,J,IDP,NS) + PART(I,J)%VP*PART(I,J)%VP

!!- <wp*wp>
        MEAN_PART_CHAN_LOC(6,J,IDP,NS) = MEAN_PART_CHAN_LOC(6,J,IDP,NS) + PART(I,J)%WP*PART(I,J)%WP

!!- <up*vp>
        MEAN_PART_CHAN_LOC(7,J,IDP,NS) = MEAN_PART_CHAN_LOC(7,J,IDP,NS) + PART(I,J)%UP*PART(I,J)%VP

!!- <up*wp>
        MEAN_PART_CHAN_LOC(8,J,IDP,NS) = MEAN_PART_CHAN_LOC(8,J,IDP,NS) + PART(I,J)%UP*PART(I,J)%WP

!!- <vp*wp>
        MEAN_PART_CHAN_LOC(9,J,IDP,NS) = MEAN_PART_CHAN_LOC(9,J,IDP,NS) + PART(I,J)%VP*PART(I,J)%WP

!!--------------------------------------------------------------------
!! Fluid at Particle position kinetic stress
!!--------------------------------------------------------------------
!!- <uf@p*uf@p>
        MEAN_PART_CHAN_LOC(13,J,IDP,NS) = MEAN_PART_CHAN_LOC(13,J,IDP,NS) + UF_P*UF_P

!!- <vf@p*vf@p
        MEAN_PART_CHAN_LOC(14,J,IDP,NS) = MEAN_PART_CHAN_LOC(14,J,IDP,NS) + VF_P*VF_P

!!- <wf@p*wf@p>
        MEAN_PART_CHAN_LOC(15,J,IDP,NS) = MEAN_PART_CHAN_LOC(15,J,IDP,NS) + WF_P*WF_P

!!-------------------------------------------------------------------
!! Particle Angular Velocity
!!-------------------------------------------------------------------
!!- <omega_px>
        MEAN_PART_CHAN_LOC(16,J,IDP,NS) = MEAN_PART_CHAN_LOC(16,J,IDP,NS) + PART(I,J)%OMEGAX

!!- <omega_py>
        MEAN_PART_CHAN_LOC(17,J,IDP,NS) = MEAN_PART_CHAN_LOC(17,J,IDP,NS) + PART(I,J)%OMEGAY

!!- <omega_pz>
        MEAN_PART_CHAN_LOC(18,J,IDP,NS) = MEAN_PART_CHAN_LOC(18,J,IDP,NS) + PART(I,J)%OMEGAZ

!!--------------------------------------------------------------------
!! Particle Rotational kinetic stress
!!--------------------------------------------------------------------
!!- <omega_px.omega_px>
        MEAN_PART_CHAN_LOC(19,J,IDP,NS) = MEAN_PART_CHAN_LOC(19,J,IDP,NS) + PART(I,J)%OMEGAX*PART(I,J)%OMEGAX

!!- <omega_py.omega_py>
        MEAN_PART_CHAN_LOC(20,J,IDP,NS) = MEAN_PART_CHAN_LOC(20,J,IDP,NS) + PART(I,J)%OMEGAY*PART(I,J)%OMEGAY

!!- <omega_pz.omega_pz>
        MEAN_PART_CHAN_LOC(21,J,IDP,NS) = MEAN_PART_CHAN_LOC(21,J,IDP,NS) + PART(I,J)%OMEGAZ*PART(I,J)%OMEGAZ

!!- <omega_px.omega_py>
        MEAN_PART_CHAN_LOC(22,J,IDP,NS) = MEAN_PART_CHAN_LOC(22,J,IDP,NS) + PART(I,J)%OMEGAX*PART(I,J)%OMEGAY

!!- <omega_px.omega_pz>
        MEAN_PART_CHAN_LOC(23,J,IDP,NS) = MEAN_PART_CHAN_LOC(23,J,IDP,NS) + PART(I,J)%OMEGAX*PART(I,J)%OMEGAZ

!!- <omega_py.omega_pz>
        MEAN_PART_CHAN_LOC(24,J,IDP,NS) = MEAN_PART_CHAN_LOC(24,J,IDP,NS) + PART(I,J)%OMEGAY*PART(I,J)%OMEGAZ

!!-------------------------------------------------------------------
!! Force on the particle
!!-------------------------------------------------------------------
        call FLUID_FORCE_ELLIPSOID(J, IDP, &
                                   PART(I,J)%UP, PART(I,J)%VP, PART(I,J)%WP, &
                                   UF_P, VF_P, WF_P, &
                                   PART(I,J)%ELLQUAT, &
                                   FDRAG, FLIFT)

!!- <F Drag X>
        MEAN_PART_CHAN_LOC(25,J,IDP,NS) = MEAN_PART_CHAN_LOC(25,J,IDP,NS) + FDRAG(1,1)
!!- <F Drag Y>
        MEAN_PART_CHAN_LOC(26,J,IDP,NS) = MEAN_PART_CHAN_LOC(26,J,IDP,NS) + FDRAG(2,1)
!!- <F Drag Z>
        MEAN_PART_CHAN_LOC(27,J,IDP,NS) = MEAN_PART_CHAN_LOC(27,J,IDP,NS) + FDRAG(3,1)

!!- <F Lift X>
        MEAN_PART_CHAN_LOC(28,J,IDP,NS) = MEAN_PART_CHAN_LOC(28,J,IDP,NS) + FLIFT(1,1)
!!- <F Lift Y>
        MEAN_PART_CHAN_LOC(29,J,IDP,NS) = MEAN_PART_CHAN_LOC(29,J,IDP,NS) + FLIFT(2,1)
!!- <F Lift Z>
        MEAN_PART_CHAN_LOC(30,J,IDP,NS) = MEAN_PART_CHAN_LOC(30,J,IDP,NS) + FLIFT(3,1)

!!-------------------------------------------------------------------
!! Torque on the particle
!!-------------------------------------------------------------------
        call FLUID_TORQUE_ELLIPSOID(I, J, IDP, &
                            PART(I,J)%OMEGAX, PART(I,J)%OMEGAY, PART(I,J)%OMEGAZ, &
                            PART(I,J)%ELLQUAT, &
                            FTORQUE_PITCH, FTORQUE_ROT)

!!- <Pitching Torque X>
        MEAN_PART_CHAN_LOC(31,J,IDP,NS) = MEAN_PART_CHAN_LOC(31,J,IDP,NS) + FTORQUE_PITCH(1,1)

!!- <Pitching Torque Y>
        MEAN_PART_CHAN_LOC(32,J,IDP,NS) = MEAN_PART_CHAN_LOC(32,J,IDP,NS) + FTORQUE_PITCH(2,1)

!!- <Pitching Torque Z>
        MEAN_PART_CHAN_LOC(33,J,IDP,NS) = MEAN_PART_CHAN_LOC(33,J,IDP,NS) + FTORQUE_PITCH(3,1)

!!- <Rotation Torque X>
        MEAN_PART_CHAN_LOC(34,J,IDP,NS) = MEAN_PART_CHAN_LOC(34,J,IDP,NS) + FTORQUE_ROT(1,1)

!!- <Rotation Torque Y>
        MEAN_PART_CHAN_LOC(35,J,IDP,NS) = MEAN_PART_CHAN_LOC(35,J,IDP,NS) + FTORQUE_ROT(2,1)

!!- <Rotation Torque Z>
        MEAN_PART_CHAN_LOC(36,J,IDP,NS) = MEAN_PART_CHAN_LOC(36,J,IDP,NS) + FTORQUE_ROT(3,1) 

!!-------------------------------------------------------------------
!! Drag Lift & Torque Coefficients for the particle
!!-------------------------------------------------------------------
        ! Relative velocity of particle relative to fluid 
        UREL(1,1) = UF_P - PART(I,J)%UP
        UREL(2,1) = VF_P - PART(I,J)%VP
        UREL(3,1) = WF_P - PART(I,J)%WP

        ! Magnitude of Relative velocity
        UREL_mag = sqrt(UREL(1,1)**2 + UREL(2,1)**2 + UREL(3,1)**2)

        ! Equivalent Diameter !! Diameter of a sphere (pi/6*dequiv^3) with same volume as ellipsoid (4/3*pi*a^3/beta^2)
        dequiv = 2.0*EMAJ_PART(J,IDP)/(APR_PART(J)**(2.0/3.0))

        ! Particle Reynolds Number
        Rep = (UREL_mag*dequiv)/VISC

        ! Transform the fluid velocity from the fixed frame to the particle frame
        ! Initialize fluid velocity vector with fluid velocity from fixed frame
        UFLU_PFRAME(1,1) = UF_P - PART(I,J)%UP
        UFLU_PFRAME(2,1) = VF_P - PART(I,J)%VP
        UFLU_PFRAME(3,1) = WF_P - PART(I,J)%WP

        call transform_basis(UFLU_PFRAME, conj_q(PART(I,J)%ELLQUAT), shape(UFLU_PFRAME))

        phi = abs(ATAN((sqrt((UFLU_PFRAME(2,1))**2 + UFLU_PFRAME(3,1)**2))/(UFLU_PFRAME(1,1))))

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !       Calculation of Drag Coefficient and Lift Coefficient using correlations     !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! C_D
        CDRAG_PHI0  = a1/Rep**a2 + a3/Rep**a4

        CDRAG_PHI90 = a5/Rep**a6 + a7/Rep**a8 

        CDRAG = CDRAG_PHI0 + (CDRAG_PHI90 - CDRAG_PHI0)*(sin(phi))**a0

        ! C_L
        CLIFT = (b1/Rep**b2 + b3/Rep**b4) * (sin(phi))**(b5 + b6*(Rep**b7)) * (cos(phi))**(b8 + b9*(Rep**b10))
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!- <Rep>
        MEAN_PART_CHAN_LOC(37,J,IDP,NS) = MEAN_PART_CHAN_LOC(37,J,IDP,NS) + Rep

!!- <phi>
        MEAN_PART_CHAN_LOC(38,J,IDP,NS) = MEAN_PART_CHAN_LOC(38,J,IDP,NS) + phi

!!- <CDRAG>
        MEAN_PART_CHAN_LOC(39,J,IDP,NS) = MEAN_PART_CHAN_LOC(39,J,IDP,NS) + CDRAG

!!- <CLIFT>
        MEAN_PART_CHAN_LOC(40,J,IDP,NS) = MEAN_PART_CHAN_LOC(40,J,IDP,NS) + CLIFT

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Relative Rotational Velocity !!
        curl_UFLU(1,1) = PART(I,J)%VORTXAP
        curl_UFLU(2,1) = PART(I,J)%VORTYAP
        curl_UFLU(3,1) = PART(I,J)%VORTZAP 


        call transform_basis(curl_UFLU, conj_q(PART(I,J)%ELLQUAT), shape(curl_UFLU))

        OREL(1,1) = 0.5*curl_UFLU(1,1) - PART(I,J)%OMEGAX
        OREL(2,1) = 0.5*curl_UFLU(2,1) - PART(I,J)%OMEGAY
        OREL(3,1) = 0.5*curl_UFLU(3,1) - PART(I,J)%OMEGAZ

        OREL_MAG = sqrt(OREL(1,1)**2 + OREL(2,1)**2 + OREL(3,1)**2)

        OREL_YZ = sqrt(OREL(3,1)**2 + OREL(2,1)**2)

        !! Rotational Reynolds Number !!
        Re_R = (OREL_MAG * dequiv**2)/(VISC)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!   Calculation of Pitching & Rotational Torque Coefficient using correlations    !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! CTPITCH !!
        CTPITCH = (c1/Rep**c2 + c3/Rep**c4) * (sin(phi))**(c5 + c6*(Rep**c7)) * (cos(phi))**(c8 + c9*(Rep**c10))

        !! CTROTATION !!
        CTROTATION = r1*(Re_R**r2) + r3/(Re_R**r4)

        CTROTATION_2 = r1_ * (Re_R ** r2_) + r3_/(Re_R ** r4_)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!- <CTPITCH>
        MEAN_PART_CHAN_LOC(41,J,IDP,NS) = MEAN_PART_CHAN_LOC(41,J,IDP,NS) + CTPITCH

!!- <CTROTATION>
        MEAN_PART_CHAN_LOC(42,J,IDP,NS) = MEAN_PART_CHAN_LOC(42,J,IDP,NS) + CTROTATION

!!- <CTROTATION_2>
        MEAN_PART_CHAN_LOC(43,J,IDP,NS) = MEAN_PART_CHAN_LOC(43,J,IDP,NS) + CTROTATION_2
                      

!!-------------------------------------------------------------------
!! Triple correlations for particles
!!-------------------------------------------------------------------


!!-------------------------------------------------------------------
!! Moment of the drag force
!!-------------------------------------------------------------------

    end do ! do I = 1, NPART_LOC(J)

    somme = somme + sum(MEAN_PART_CHAN_LOC(44,J,:,:))

end do ! do J = 1, NIG




if (somme .ne. NPART_LOC(1)) then
    
    write(*,*) somme, NPART_LOC(1)
    write(*,*) 'All Particles have not been accounted for'
    stop

end if




!!====================================================================
!! Summation overall domain and normalization
!!====================================================================
!if (NPROC>1) then
do J = 1, NIG

        do IDP = 1, POLYDISP(J)

	       do NS = 1, NSLICE

                call MPI_ALLREDUCE(MEAN_PART_CHAN_LOC(:,J,IDP,NS),MEAN_PART_CHAN(:,J,IDP,NS),NSTAT,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

               end do

        end do

end do


!!======================================================================
!! Time averaging if so
!!======================================================================
if(STAT_TIME) then
 MEAN_TIME_PART_CHAN = MEAN_TIME_PART_CHAN + MEAN_PART_CHAN
end if 



!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()
 CPU_PART(6) = CPU_PART(6) + TIME_END - TIME_START
end if


10000 format (30(e17.7))


end subroutine STAT_PART_CHANNEL
