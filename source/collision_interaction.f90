!!====================================================================
!!
!!          Particle-particle interaction
!!
!!====================================================================

subroutine COLLISION_INTERACTION(TIME,NCYCLE,IG,NCLOSE,HOME,NEAR)

!!====================================================================
!!
!!
!!====================================================================

use PARTICLE_PARALLEL
use DNS_DIM               
use PARAM_PHYS
use GEOMETRIC_VARIABLE
use COLLISION_VARIABLE
use STATISTICS
use mod_quaternion

implicit none

!!====================================================================
!! ARRAYS
!!====================================================================
!- Curent time
real(kind=8),               intent(in) :: TIME

!!- Cycle number
integer,                    intent(in) :: NCYCLE

!!- particle specy
integer,                    intent(in) :: IG

!!- Number of particle "close"
integer,                    intent(inout) :: NCLOSE

!!- Index of particle "home"
integer, dimension(NCLOSE), intent(inout) :: HOME

!!- Index of particle "near"
integer, dimension(NCLOSE), intent(inout) :: NEAR

!!====================================================================

!- Collision diameter
real(kind=8) :: DCOLL

!- Impact parameter
real(kind=8) :: VKX, VKY, VKZ
real(kind=8) :: NRM_VK

!- Relative velocity
real(kind=8) :: WX, WY, WZ
real(kind=8) :: NRM_W

!- Scalar product between k and w
real(kind=8) :: WK
real(kind=8) :: THETA

real(kind=8):: M1, M2, LI, BI, DELTAT


real(kind=8):: RDUMMY

real(kind=8) :: XNEAR, YNEAR, ZNEAR
real(kind=8) :: VNEAR

logical :: flag

! Vector linking particle centers
real(kind=8), dimension(ndim, 1) :: VPC, P2P1

! Output of the contact detection algorithm
! Point on Ellipsoid 1
real(kind=8), dimension(ndim, 1) :: EPT1

! Point on Ellipsoid 2
real(kind=8), dimension(ndim, 1) :: EPT2

! Common normal
real(kind=8), dimension(ndim, 1) :: cnormal

! Ellipsoid 1
real(kind=8) :: UP1, VP1, WP1 ! Velocity of contact points
real(kind=8), dimension(ndim, 1) :: RX1 ! Dist. of ct Pt from center
real(kind=8), dimension(ndim, 1) :: OMG1 ! Comoving frame angular velocity

! Ellipsoid 2
real(kind=8) :: UP2, VP2, WP2 ! Velocity of contact points
real(kind=8), dimension(ndim, 1) :: RX2! Dist. of ct Pt from center
real(kind=8), dimension(ndim, 1) :: OMG2 ! Comoving frame angular velocity

! Relative velocity of particles
real(kind=8) :: VRELX, VRELY, VRELZ, NRM_VREL

real(kind=8) :: NVREL

real(kind=8), dimension(ndim, 1) :: U_P, R_impact1


! Overlap of particles
real(kind=8) :: overlap_depth


!!--------------------------------------------------------------------
real(kind=8), dimension(NCLOSE) :: DIST
integer,      dimension(NCLOSE) :: IHOME, INEAR
!!--------------------------------------------------------------------
integer :: NCYCLESTAT


!!- Index
integer :: N, NP1, NP2, NCLOSE_SORT
integer :: I, IDP1, IDP2, NS, IPDF
integer :: NUMFILE
integer :: NCYCLE_TRACK, NP1_TRACK, NP2_TRACK

!!====================================================================
! Local Statistics
real(kind=8), dimension(NSTAT, POLYDISPMAX, POLYDISPMAX) :: MEAN_PART_COL

!!====================================================================
!! In case of parallel computation the particles are sorted according 
!! to the distance between neighboring particles (from smallest to
!! the largest).
!!====================================================================
!if(NPROC>1) then

! NCLOSE_SORT = 0

! do N=1, NCLOSE

!  NP1 = HOME(N)
!  NP2 = NEAR(N)

!  IDP1 = PART(NP1,IG)%IDP
!  IDP2 = PART(NP2,IG)%IDP 

  !- Collision diameter
!  DCOLL = 0.5*(DPART(IG,IDP1) + DPART(IG,IDP2))

  !- Impact parameter
!  VKX = PART(NP1,IG)%XP - PART(NP2,IG)%XP
!  VKY = PART(NP1,IG)%YP - PART(NP2,IG)%YP
!  VKZ = PART(NP1,IG)%ZP - PART(NP2,IG)%ZP

  !- Periodicity only in x-direction
!  if((abs(VKX) > LXMAX/2.0) .and. PERIODICITY) VKX = VKX - SIGN(LXMAX,VKX)

!  NRM_VK = sqrt(VKX*VKX + VKY*VKY + VKZ*VKZ)

!  if((PART(NP1,IG)%MYIDP==MYID.or.PART(NP2,IG)%MYIDP==MYID).and.NRM_VK<=DCOLL) then

!   NCLOSE_SORT = NCLOSE_SORT+1
!   IHOME(NCLOSE_SORT) = NP1
!   INEAR(NCLOSE_SORT) = NP2
!   DIST(NCLOSE_SORT) = sqrt(VKX**2+VKY**2+VKZ**2)
!  end if
! end do

! call QUICKSORT(NCLOSE_SORT,IHOME(1:NCLOSE_SORT),INEAR(1:NCLOSE_SORT),DIST(1:NCLOSE_SORT))

! HOME(1:NCLOSE_SORT)=IHOME(1:NCLOSE_SORT)
! NEAR(1:NCLOSE_SORT)=INEAR(1:NCLOSE_SORT)

! NCLOSE = NCLOSE_SORT
!end if

MEAN_PART_COL(:,:,:) = ZERO
!MEAN_PART_COL_PDF(:,:,:,:) = ZERO
!!====================================================================
!!====================================================================
!! Collision Statistics Initialisation (Not needed Now)
!!====================================================================
!!====================================================================
!MIN_PDF = -5.0
!MAX_PDF = 5.0
!MIN_THETA = -2.0*PPI
!MAX_THETA = 2.0*PPI


!DPDF = (MAX_PDF - MIN_PDF)/real(NPDF)
!DPDF_THETA = (MAX_THETA - MIN_THETA)/real(NPDF)

!!-----------------------!!-----------------------
!!-----------------------!!-----------------------

PART(:,IG)%COLOR = 1
PART(:,IG)%COLOR = 1


do N = 1, NCLOSE

    NP1 = HOME(N)

    NP2 = NEAR(N)

    IDP1 = PART(NP1,IG)%IDP

    IDP2 = PART(NP2,IG)%IDP 

    flag = .false.

    !!- Allow only one collision per particle
    if(PART(NP1,IG)%COLOR==1 .and. PART(NP2,IG)%COLOR==1) then

        !!====================================================================
        !! 1. Compute the impact vector
        !!====================================================================
        !! The impact vector is the vector linking the centres of the two
        !! neighboring particles.
        !!--------------------------------------------------------------------
        
        XNEAR = PART(NP1,IG)%XP
        YNEAR = PART(NP1,IG)%YP
        ZNEAR = PART(NP1,IG)%ZP

        !- Impact parameters
        VKX = PART(NP2,IG)%XP - XNEAR
        VKY = PART(NP2,IG)%YP - YNEAR
        VKZ = PART(NP2,IG)%ZP - ZNEAR

        !- Periodicity in x-direction
        if((abs(VKX) > LXMAX/2) .and. PERIODICITY) VKX = VKX - SIGN(LXMAX,VKX)

        !- Periodicity in y- and z-direction only if scalar computation
        if(NPROC==1) then

            !- Periodicity only in y-direction
            if(abs(VKY) > LYMAX/2.) VKY = VKY - SIGN(LYMAX,VKY)

            !- Periodicity only in z-direction
            if(abs(VKZ) > LZMAX/2.) VKZ = VKZ - SIGN(LZMAX,VKZ)

        end if


        call ellipsoid_contact_detection(IG, IDP1, IDP2,         &
                                         PART(NP1, IG)%XP,       &
                                         PART(NP1, IG)%YP,       &
                                         PART(NP1, IG)%ZP,       &
                                         PART(NP1, IG)%ELLQUAT,  &
                                         PART(NP1, IG)%XP + VKX, & 
                                         PART(NP1, IG)%YP + VKY, &
                                         PART(NP1, IG)%ZP + VKZ, &
                                         PART(NP2, IG)%ELLQUAT,  &
                                         flag, EPT1, EPT2, cnormal)


        !!====================================================================
        !! 2. Overlapped particles
        !!====================================================================
        if(flag) then

            ! For 2 overlapped ellipsoids we need to figure out whether they are approaching
            ! or going away, this can be accurately done by considering the point on the  
            ! ellipsoid surface having the maximum encroachment into the other ellipsoid

            ! Calculate the vector linking the max. encroachment points                
            VPC(1,1) = EPT1(1,1) - EPT2(1,1) ! (EPT1(1,1) + EPT2(1,1))/2.0 !
            VPC(2,1) = EPT1(2,1) - EPT2(2,1) ! (EPT1(2,1) + EPT2(2,1))/2.0 ! 
            VPC(3,1) = EPT1(3,1) - EPT2(3,1) ! (EPT1(3,1) + EPT2(3,1))/2.0 ! 


            ! Calculate the velocity of the max. encroachment points
            ! Particle 1
            RX1(1,1) = (EPT1(1,1) - PART(NP1, IG)%XP) ! (VPC(1,1) - PART(NP1, IG)%XP)  ! 
            RX1(2,1) = (EPT1(2,1) - PART(NP1, IG)%YP) ! (VPC(2,1) - PART(NP1, IG)%YP)  ! 
            RX1(3,1) = (EPT1(3,1) - PART(NP1, IG)%ZP) ! (VPC(3,1) - PART(NP1, IG)%ZP)  ! 

            OMG1(1,1) = PART(NP1, IG)%OMEGAX
            OMG1(2,1) = PART(NP1, IG)%OMEGAY
            OMG1(3,1) = PART(NP1, IG)%OMEGAZ

            call transform_basis(OMG1, PART(NP1, IG)%ELLQUAT, shape(OMG1))

            UP1 = PART(NP1, IG)%UP + (OMG1(2,1)*RX1(3,1) - OMG1(3,1)*RX1(2,1))
            VP1 = PART(NP1, IG)%VP + (OMG1(3,1)*RX1(1,1) - OMG1(1,1)*RX1(3,1))
            WP1 = PART(NP1, IG)%WP + (OMG1(1,1)*RX1(2,1) - OMG1(2,1)*RX1(1,1))

            ! Particle 2
            RX2(1,1) = (EPT2(1,1) - (PART(NP1, IG)%XP + VKX)) ! (VPC(1,1) - (PART(NP1, IG)%XP + VKX)) ! 
            RX2(2,1) = (EPT2(2,1) - (PART(NP1, IG)%YP + VKY)) ! (VPC(2,1) - (PART(NP1, IG)%YP + VKY)) ! 
            RX2(3,1) = (EPT2(3,1) - (PART(NP1, IG)%ZP + VKZ)) ! (VPC(3,1) - (PART(NP1, IG)%ZP + VKZ)) ! 


            OMG2(1,1) = PART(NP2, IG)%OMEGAX
            OMG2(2,1) = PART(NP2, IG)%OMEGAY
            OMG2(3,1) = PART(NP2, IG)%OMEGAZ

            call transform_basis(OMG2, PART(NP2, IG)%ELLQUAT, shape(OMG2))


            UP2 = PART(NP2, IG)%UP + (OMG2(2,1)*RX2(3,1) - OMG2(3,1)*RX2(2,1))
            VP2 = PART(NP2, IG)%VP + (OMG2(3,1)*RX2(1,1) - OMG2(1,1)*RX2(3,1))
            WP2 = PART(NP2, IG)%WP + (OMG2(1,1)*RX2(2,1) - OMG2(2,1)*RX2(1,1))

            ! Calculate the relative velocity
            VRELX = UP2 - UP1
            VRELY = VP2 - VP1
            VRELZ = WP2 - WP1

            ! Modulus of relative velocity
            NRM_VREL = sqrt(VRELX*VRELX + VRELY*VRELY + VRELZ*VRELZ)


            !P2P1(1,1) = EPT1(1,1) - EPT2(1,1) !
            !P2P1(2,1) = EPT1(2,1) - EPT2(2,1) !
            !P2P1(3,1) = EPT1(3,1) - EPT2(3,1) !

            !if(sqrt(P2P1(1,1)**2 + P2P1(2,1)**2 + P2P1(3,1)**2) == 0.0) then

            !    P2P1 = cnormal

            !else 

            !    P2P1 = P2P1/sqrt(P2P1(1,1)**2 + P2P1(2,1)**2 + P2P1(3,1)**2)

            !end if        

            ! Dot Product of Relative Velocity and Vector linking encroachment points)
            NVREL = VRELX*cnormal(1,1) + VRELY*cnormal(2,1) + VRELZ*cnormal(3,1)
            !NVREL = VRELX*P2P1(1,1) + VRELY*P2P1(2,1) + VRELZ*P2P1(3,1)

            ! Collision angle
            THETA = acos(NVREL/NRM_VREL)

            if(mod(NCYCLE, FOUT2)==0 .and. LEVEL1_STPAR) then

                ! Number of overlapped pairs 
                MEAN_PART_COL(1, IDP1, IDP2) = MEAN_PART_COL(1, IDP1, IDP2) + 1.0

                ! Overlap depth
                overlap_depth = sqrt(VPC(1,1)**2 + VPC(2,1)**2 + VPC(3,1)**2)
                MEAN_PART_COL(4, IDP1, IDP2) = MEAN_PART_COL(4, IDP1, IDP2) + overlap_depth

                MEAN_PART_COL(6, IDP1, IDP2) = MEAN_PART_COL(6, IDP1, IDP2) + NVREL ! Radial relative velocity

                MEAN_PART_COL(7, IDP1, IDP2) = MEAN_PART_COL(7, IDP1, IDP2) + NRM_VREL ! Modulus of relative velocity

                MEAN_PART_COL(8, IDP1, IDP2) = MEAN_PART_COL(8, IDP1, IDP2) + THETA ! Collision Angle

            end if

            
            if(NVREL < ZERO) then ! Particles are approaching !!

                !! Collsion Response
                
                !write(*,*) '   '
                !write(*,*) 'Before Collision '
                
                !write(*,*)'1 Translational ', PART(NP1, IG)%UP, PART(NP1, IG)%VP, PART(NP1, IG)%WP
                !write(*,*)'1 Rotational ', PART(NP1, IG)%OMEGAX, PART(NP1, IG)%OMEGAY, PART(NP1, IG)%OMEGAZ

                !write(*,*)'2 Translational ', PART(NP2, IG)%UP, PART(NP2, IG)%VP, PART(NP2, IG)%WP
                !write(*,*)'2 Rotational ', PART(NP2, IG)%OMEGAX, PART(NP2, IG)%OMEGAY, PART(NP2, IG)%OMEGAZ

                if(LEVEL1_STPDF) then

                    call STAT_PARTICLE_COLLISION(IG, IDP1,             & 
                                                VRELX, VRELY,  VRELZ,  &
                                                PART(NP1, IG)%ELLQUAT, &
                                                PART(NP2, IG)%ELLQUAT, &
                                                cnormal)


                    if(FOUT3>0 .and. NEVEN_COLL_PDF<200002) then

                        call print_colliding_pair(IG, NP1, NP2, NEVEN_COLL_PDF, &
                                                  RX1, RX2, cnormal) 

                    end if

                end if
                

                ! Calculate the impulse on the particles and Update the velocities
                call  COLLISION_RESPONSE(IG, IDP1, &
                                        PART(NP1, IG),   &
                                        PART(NP2, IG),   &
                                        RX1, RX2, cnormal)

                !write(*,*) '   '
                !write(*,*) 'After Collision '
                
                !write(*,*)'1 Translational ', PART(NP1, IG)%UP, PART(NP1, IG)%VP, PART(NP1, IG)%WP
                !write(*,*)'1 Rotational ', PART(NP1, IG)%OMEGAX, PART(NP1, IG)%OMEGAY, PART(NP1, IG)%OMEGAZ

                !write(*,*)'2 Translational ', PART(NP2, IG)%UP, PART(NP2, IG)%VP, PART(NP2, IG)%WP
                !write(*,*)'2 Rotational ', PART(NP2, IG)%OMEGAX, PART(NP2, IG)%OMEGAY, PART(NP2, IG)%OMEGAZ
                !write(*,*) '   '
                

                PART(NP1,IG)%COLOR = -1.0
                PART(NP2,IG)%COLOR = -1.0
                PART(NP1,IG)%COLL_FLAG = .true.
                PART(NP2,IG)%COLL_FLAG = .true.
                
                if (IDP1 .ne. IDP2) then
                    PART(NP1, IG)%COLL_BIDISP = .true.
                    PART(NP2, IG)%COLL_BIDISP = .true.
                end if 


                ! Number of collision
                if(mod(NCYCLE, FOUT2)==0 .and. LEVEL1_STPAR)  MEAN_PART_COL(2, IDP1, IDP2) = MEAN_PART_COL(2, IDP1, IDP2) + 1.0 

                !stop

            end if !- end if(NVREL<ZERO)

        end if !- end if(flag)

    end if !- end if color==1

end do !- end loop do N = 1, NCLOSE

!stop

if(mod(NCYCLE, FOUT2)==0 .and. LEVEL1_STPAR) then

    do IDP1 = 1, POLYDISP(IG)
        
        do IDP2 = 1, POLYDISP(IG)

            MEAN_PART_COL(3, IDP1, IDP2) = 2.0*MEAN_PART_COL(2, IDP1, IDP2)/DTIME
            MEAN_PART_COL(5, IDP1, IDP2) = MEAN_PART_COL(4, IDP1, IDP2)/MEAN_PART_COL(1, IDP1, IDP2)

            NUMFILE = 560 
            NUMFILE = NUMFILE + sum(POLYDISP(1:IG)**2)-POLYDISP(IG)**2
            NUMFILE = NUMFILE + (IDP1-1)*POLYDISP(IG)+IDP2

            write(NUMFILE, 10000)       TIME, &
                MEAN_PART_COL(1, IDP1, IDP2), & ! Pairs of Overlapped particles
                MEAN_PART_COL(2, IDP1, IDP2), & ! Number of collision
                MEAN_PART_COL(3, IDP1, IDP2), & ! Frequency of collision
                MEAN_PART_COL(4, IDP1, IDP2)/EMAJ_PART(IG,IDP1), & ! Overlap
                MEAN_PART_COL(5, IDP1, IDP2), & ! Overlap normalized by number of pair of particles 
                MEAN_PART_COL(6, IDP1, IDP2), & ! Radial relative velocity
                MEAN_PART_COL(7, IDP1, IDP2), & ! Modulus of relative velocity
                MEAN_PART_COL(8, IDP1, IDP2)    ! Collision Angle

        end do

    end do

end if


if(mod(NCYCLE, FOUT2)==0 .and. LEVEL1_STPAR .and. STAT_TIME) MEAN_TIME_PART_COL = MEAN_TIME_PART_COL + MEAN_PART_COL
                    

!!--------------------------------------------------------------------
10000 format (30(e17.7))
10703 format (2x,' Part. in Class ',i3,' :  Max=',i7,' Min=',i7)

end subroutine COLLISION_INTERACTION