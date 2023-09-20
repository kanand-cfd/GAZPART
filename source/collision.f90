!!====================================================================
!!
!!          Particle-particle interaction
!!
!!====================================================================

subroutine COLLISION(NCYCLE,TIME)

!!====================================================================
!!
!!
!!====================================================================

use STATISTICS 
use PARTICLE_PARALLEL
use MPI_STRUCTURES
use DNS_DIM               
use PARAM_PHYS
use CHECK_CPU
use COLLISION_VARIABLE

implicit none

!!====================================================================
!! ARRAYS
!!====================================================================
!!- cycle number
integer, intent(in) :: NCYCLE

!- Curent time
real(kind=8), intent(in) :: TIME
!!====================================================================


integer, parameter :: NPCLOSEMAX = 100000000
integer, dimension(NPCLOSEMAX) :: HOME
integer, dimension(NPCLOSEMAX) :: NEAR

integer :: NCLOSE

integer :: NP1, NP2
!!====================================================================

!!- Index
integer :: I, J ,K
integer :: NS, NF 

!!====================================================================


if(STAT_TIME .and. LEVEL0_STPAR .and. LEVELX_STPAR .and. mod(NCYCLE, FOUT2)==0) call BEFORE_COLLISION

do J = 1, NIG

  do I = 1, NPART_LOC(J)

    PART(I,J)%MYIDP = MYID

  end do

end do


!- Loop over particle class
do J = 1, NIG

  !- Only some particle classes are allowed to have collisions

  if(COLLIDEF(J)==1) then

    !!====================================================================
    !! 1. Manage multi-task collision algorithm
    !!====================================================================

    !!--------------------------------------------------------------------
    !! 1.2 Gather all particles from neighbours
    !!--------------------------------------------------------------------

    if(NPROC>1) then

      !!- Counting 
      call PARTICLES_COUNTING_COLL(J)

      !!- Merge
      call EXCHANGE_P_COLL(J,NPNEIGH_COLL)

      NS = NPART_LOC(J)+1
      NF = NPART_LOC(J)+NPNEIGH_COLL

      !!--------------------------------------------------------------------
      !! 1.3 For all ghost particles shift for periodicity in y- and z-direction
      !!--------------------------------------------------------------------
      do I = NS, NF

        if(YMESH(ISTART(2)) <= ZERO) then

          if (PART(I,J)%YP > LYMAX - DELTA_COLL .and. PART(I,J)%YP <= LYMAX) then

            PART(I,J)%YP = PART(I,J)%YP - LYMAX

          end if

        end if



        if(YMESH(IEND(2))+DY >= LYMAX) then

          if (PART(I,J)%YP >= ZERO .and. PART(I,J)%YP < DELTA_COLL) then

            PART(I,J)%YP = PART(I,J)%YP + LYMAX

          end if

        end if

        

        if(ZMESH(ISTART(3)) <= ZERO) then

          if (PART(I,J)%ZP > LZMAX - DELTA_COLL .and. PART(I,J)%ZP <= LZMAX) then

            PART(I,J)%ZP = PART(I,J)%ZP - LZMAX

          end if

        end if

        

        if(ZMESH(IEND(3))+DZ >= LZMAX) then

          if (PART(I,J)%ZP >= ZERO .and. PART(I,J)%ZP < DELTA_COLL) then

            PART(I,J)%ZP = PART(I,J)%ZP + LZMAX

          end if

        end if

      

      end do !- ending I = NS, NF

    
    else 

      NPNEIGH_COLL = 0
      NS = NPART_LOC(J)
      NF = NPART_LOC(J)

    end if !!- end if NPROC>1

    !!====================================================================
    !! 2. Detection of neighboring particles
    !!====================================================================
    call COLLISION_DETECTION_CONNECTIVITY(NF, &
                                          PART(1:NF,J)%XP, &
                                          PART(1:NF,J)%YP, &
                                          PART(1:NF,J)%ZP, &
                                          NCLOSE, & 
                                          NPCLOSEMAX, &
                                          HOME, &
                                          NEAR  )

    !!====================================================================
    !! 4. Particle-particle interaction
    !!====================================================================
    call COLLISION_INTERACTION(TIME,NCYCLE,J,NCLOSE,HOME,NEAR)


  end if !- end if(COLLIDEF(J)==1)

end do !- end do J = 1, NIG


if(STAT_TIME .and. LEVEL0_STPAR .and. LEVELX_STPAR .and. mod(NCYCLE, FOUT2)==0) call AFTER_COLLISION


!if(STAT_TIME .and. LEVELX_STPAR) then

!  NEVEN_COLL_PDF = NEVEN_COLL_PDF + 1

!  NEVEN_PAIR_PDF = NEVEN_PAIR_PDF + 1

!end if 


!!====================================================================
!! 3. Periodic boundary conditions
!!====================================================================
!!- Boundary conditions if particles move
! call BOUNDARY_PARTICLE


!!--------------------------------------------------------------------
10000 format (30(e17.7))
10703 format (2x,' Part. in Class ',i3,' :  Max=',i7,' Min=',i7)

end subroutine COLLISION

