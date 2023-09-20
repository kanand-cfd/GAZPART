!!====================================================================
!!
!!   Check remaining time
!!
!!====================================================================

subroutine TREMAIN(ANSWER)

use DNS_DIM
use CPUTIME_CONTROL

implicit none

logical, intent(out) :: ANSWER

real(kind=8) :: CPU_CURRENT
!!--------------------------------------------------------------------

ANSWER = .false.

CPU_CURRENT = MPI_WTIME()


!!if (rang==0) write(*,*) 'Time', CPU_CURRENT - CPU_INIT

if (WALLTIME-(CPU_CURRENT - CPU_INIT) < 240.) then
 if (MYID == 0)write(*,*) 'Warning Remainig CPU Time = ',WALLTIME-(CPU_CURRENT - CPU_INIT)
  ANSWER = .true.
end if 


end subroutine TREMAIN
