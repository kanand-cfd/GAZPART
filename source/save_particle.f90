!!====================================================================
!!
!!
!!====================================================================

subroutine SAVE_PARTICLE(NCYCLE)

!!====================================================================
!!
!!
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE
use PARAM_PHYS
use PARTICLE_PARALLEL
use CHECK_CPU
use MPI_STRUCTURES


implicit none


!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Cycle number
integer, intent(in) :: NCYCLE

!- File name 
character(len=40) :: FILENAME

integer :: IDUMMY
integer :: RECSIZE, SIZE_INT, SIZE_REAL, NVARIABLE

!- Variable for MPI I/O
integer :: FILETYPE
integer :: LOC_SIZE
integer :: DESCRIPTEUR


!- temporary array
real(kind=8), dimension(NPART_FULL,14) :: TMPVAR

!- number of particles of class J for all procs
integer, dimension(NPROC,NIG) :: NPART_LOC_ALL

!Various variables for ALL_GATHERV
integer, dimension(NPROC) :: COUNTS
integer, dimension(NPROC) :: DISPLS
integer                   :: BLOCK


!- Time control variable
real(kind=8) :: TIME_START, TIME_END

!- Index
integer :: I, J, NP, P
!---------------------------------------------------------------------


!!- CPU check
if(MYID == 0) then
 TIME_START=MPI_WTIME()
end if



!!====================================================================
!! 1. Multiple binary files
!!====================================================================
if(ISAVEPART == 1) then 

!- Define file name
if(NCYCLE>=0) then
 write(FILENAME,10403)'PART',trim(FILE_EXT),'_t',NCYCLE,'.bin'
else
 write(FILENAME,10303)'PART',trim(FILE_EXT),'.end'
end if 



!- Open file containing the last particle position and velocity
open(unit = 150, file = trim(FILENAME),form='unformatted')

do J = 1, NIG

 write(150)NPART_LOC(J)
 write(150)(PART(I,J)%IDP,I=1,NPART_LOC(J))
 write(150)(PART(I,J)%XP,I=1,NPART_LOC(J))
 write(150)(PART(I,J)%YP,I=1,NPART_LOC(J))
 write(150)(PART(I,J)%ZP,I=1,NPART_LOC(J))

 write(150)(PART(I,J)%UP,I=1,NPART_LOC(J))
 write(150)(PART(I,J)%VP,I=1,NPART_LOC(J))
 write(150)(PART(I,J)%WP,I=1,NPART_LOC(J))
 
 write(150)(PART(I,J)%UFAP,I=1,NPART_LOC(J))
 write(150)(PART(I,J)%VFAP,I=1,NPART_LOC(J))
 write(150)(PART(I,J)%WFAP,I=1,NPART_LOC(J))

end do !!- Loop:J = 1, NIG


!- Close file
close(150)


if(NCYCLE>0) then 
 if(MYID==0)write(*,*) 'Particle position and velocity file dropped --> ',trim(FILENAME)
end if

!!====================================================================
!! 2. One binary file per class
!!====================================================================
elseif(ISAVEPART == 2) then

if(NCYCLE<0) then 
 if(MYID==0)write(*,*) 'Final particle position and velocity --> Binary File'
else
 if(MYID==0)write(*,*) 'Particle position and velocity file dropped --> Binary file'
end if


do J = 1, NIG


 call MPI_ALLGATHER(NPART_LOC(J), 1, MPI_INT, NPART_LOC_ALL(:,J), 1, MPI_INT, MPI_COMM_WORLD, IERR)

 if (sum(NPART_LOC_ALL(:,J)) /= NPART_FULL) stop 'MPI_ALLGATHER is wrong!!!'

 COUNTS = [ (NPART_LOC_ALL(P,J), P = 1, NPROC) ]
 DISPLS = [ (sum(NPART_LOC_ALL(:P-1,J)), P=1, NPROC) ]




 call MPI_TYPE_VECTOR(NPART_LOC(J), 1, 1, MPI_DOUBLE_PRECISION, BLOCK, IERR)
 call MPI_TYPE_COMMIT(BLOCK,IERR)     


!- Drop x-position
 call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%XP       , 1, BLOCK, TMPVAR(:,1),  COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
 call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%YP       , 1, BLOCK, TMPVAR(:,2),  COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
 call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%ZP       , 1, BLOCK, TMPVAR(:,3),  COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
 call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%UP       , 1, BLOCK, TMPVAR(:,4),  COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
 call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%VP       , 1, BLOCK, TMPVAR(:,5),  COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
 call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%WP       , 1, BLOCK, TMPVAR(:,6),  COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
 call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%ELLQUAT%a, 1, BLOCK, TMPVAR(:,7),  COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
 call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%ELLQUAT%b, 1, BLOCK, TMPVAR(:,8),  COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
 call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%ELLQUAT%c, 1, BLOCK, TMPVAR(:,9),  COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
 call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%ELLQUAT%d, 1, BLOCK, TMPVAR(:,10), COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
 call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%OMEGAX   , 1, BLOCK, TMPVAR(:,11), COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
 call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%OMEGAY   , 1, BLOCK, TMPVAR(:,12), COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
 call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%OMEGAZ   , 1, BLOCK, TMPVAR(:,13), COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
 call MPI_ALLGATHERV(real(PART(1:NPART_LOC(J),J)%IDP), 1, BLOCK, TMPVAR(:,14), COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)

!- Define file name
 if(NCYCLE>=0) then
  write(FILENAME,10102)'part_p',J,'_t',NCYCLE,'.bin'
 else
  write(FILENAME,10100)'part_p',J,'.end'
 end if 

 if(MYID==0)write(*,*) '   + file:',trim(FILENAME)


 open(unit=150,file=trim(FILENAME),status='replace',form='unformatted')

 write(150)NPART_FULL
 write(150)NPCLASS(J,:)
 write(150)RHOP(J,:)
 write(150)EMAJ_PART(J,:)
 write(150)APR_PART(J)
 write(150)(TMPVAR(I,1) ,I=1,NPART_FULL)
 write(150)(TMPVAR(I,2) ,I=1,NPART_FULL)
 write(150)(TMPVAR(I,3) ,I=1,NPART_FULL)
 write(150)(TMPVAR(I,4) ,I=1,NPART_FULL)
 write(150)(TMPVAR(I,5) ,I=1,NPART_FULL)
 write(150)(TMPVAR(I,6) ,I=1,NPART_FULL)
 write(150)(TMPVAR(I,7) ,I=1,NPART_FULL)
 write(150)(TMPVAR(I,8) ,I=1,NPART_FULL)
 write(150)(TMPVAR(I,9) ,I=1,NPART_FULL)
 write(150)(TMPVAR(I,10),I=1,NPART_FULL)
 write(150)(TMPVAR(I,11),I=1,NPART_FULL)
 write(150)(TMPVAR(I,12),I=1,NPART_FULL)
 write(150)(TMPVAR(I,13),I=1,NPART_FULL)
 write(150)(TMPVAR(I,14),I=1,NPART_FULL)
 close(150)

end do !- end do IG = 1, NIG


if(DEBUG) then
 do J=1, NIG
  call ISUMCPU(NPART_LOC(J),IDUMMY)
  if (MYID==0) write(*,*) 'Save Part -> Class:',J,' Full number of particles ',IDUMMY
 end do
end if




!!====================================================================
!! 3. MPI I/O
!!====================================================================
elseif(ISAVEPART == 3) then

FILENAME = 'traj.end'

call SAVE_PART_MPIIO(PART,FILENAME)


if(MYID==0)write(*,*) 'Final particle position and velocity --> Saved'
if(MYID==0)write(*,*) '    + MPI I/O'
if(MYID==0)write(*,10300) '     + File name = ',trim(FILENAME)


end if !!- If: ISAVEPART



!!--------------------------------------------------------------------
!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()
 CPU_PART(7) = CPU_PART(7) + TIME_END - TIME_START
end if

!!--------------------------------------------------------------------
10000 format (15(e17.7))

10100 format (A,I2.2,A)
10101 format (A,I2.2,A,A)
10102 format (A,I2.2,A,I8.8,A)

10200 format (A,F8.2,A)
10300 format (A,A)

10303 format (A,A,A)
10403 format (A,A,A,I8.8,A)
20000 format (I2,2x,I2,2x,I2,10(e17.7))

end subroutine SAVE_PARTICLE
