!!====================================================================
!!
!!   3-dimensional periodical boundary for particles
!!
!!====================================================================

subroutine EXCHANGE_PART

!!====================================================================
!!
!!
!!====================================================================

use MPI_structures
use particle_parallel

implicit none

!- Index
integer :: I, J

!-
integer :: IDUMMY

!- 
integer :: NPCHECK

!- Loop on each class of particles
do J=1,NIG

 ! Count number of particles for each proc (staying + leaving)
 call PARTICLES_COUNTING(J)


 ! NPART_LOC before exchanging particles sum staying + leaving
 NPART_LOC(J) = COUNT_STAY + SUM(COUNTER(:,MYID)) 

!!write(*,*)MYID,J,NPART_LOC(J),COUNT_STAY,SUM(COUNTER(:,MYID))

 ! Total particles number checking
 call MPI_ALLREDUCE(NPART_LOC(J),NPCHECK,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
  if (MYID==0) then
    if (NPCHECK/=NPART_FULL) then
      write(*,*) 'Particles total number ',NPCHECK,' is different from the initial number of particles ',NPART_FULL
      stop
    end if
  end if



 ! Particles processor exchange 
 call EXCHANGE_P(J)




 call IMAXCPU(maxval(COUNTER(:,MYID)),IDUMMY)
 NBR_EXCHANGE(J) = max(NBR_EXCHANGE(J),IDUMMY)


 call IMAXCPU(NPART_LOC(J),IDUMMY)
 NPMAX_CPU(J) = max(NPMAX_CPU(J),IDUMMY)

 call IMINCPU(NPART_LOC(J),IDUMMY)
 NPMIN_CPU(J) = min(NPMIN_CPU(J),IDUMMY)

end do

!------------------------------------------------------------------
end subroutine EXCHANGE_PART
