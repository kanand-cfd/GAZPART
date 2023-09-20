!!====================================================================
!!
!!   3-dimensional periodical boundary for particles
!!
!!====================================================================

subroutine BOUNDARY_PARTICLE

!!====================================================================
!!
!!
!!====================================================================

use DNS_DIM
use PARTICLE_PARALLEL
use GEOMETRIC_VARIABLE
use PARAM_PHYS 
use CHECK_CPU

implicit none



!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Time control variable
real(kind=8) :: TIME_START, TIME_END

!- Index
integer :: I, J
!---------------------------------------------------------------------

!!- CPU check
if(MYID == 0) then
 TIME_START=MPI_WTIME()
end if


!! exchange of particles leaving their processor
if (NPROC > 1) call EXCHANGE_PART



do J = 1, NIG

 do I = 1, NPART_LOC(J)

!!- x-component
  if(PART(I,J)%XP > LXMAX) PART(I,J)%XP = PART(I,J)%XP - LXMAX
  if(PART(I,J)%XP < 0.   ) PART(I,J)%XP = PART(I,J)%XP + LXMAX

!!- y-component
  if(PART(I,J)%YP > LYMAX) PART(I,J)%YP = PART(I,J)%YP - LYMAX
  if(PART(I,J)%YP < 0.   ) PART(I,J)%YP = PART(I,J)%YP + LYMAX

!!- z-component
  if(PART(I,J)%ZP > LZMAX) PART(I,J)%ZP = PART(I,J)%ZP - LZMAX
  if(PART(I,J)%ZP < 0.   ) PART(I,J)%ZP = PART(I,J)%ZP + LZMAX


 end do

end do


!! exchange of particles leaving their processor
if (NPROC > 1) call EXCHANGE_PART




!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()
 CPU_PART(5) = CPU_PART(5) + TIME_END - TIME_START
end if


!------------------------------------------------------------------
end subroutine BOUNDARY_PARTICLE
