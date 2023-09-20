!!====================================================================
!!
!!
!!====================================================================

subroutine INITIATION_SCALAR

!!====================================================================
!!
!!
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE
use SCALAR_VARIABLE
use PARAM_PHYS
use MPI_STRUCTURES

implicit none


!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------
!- File name 
character(len=40) :: FILENAME
!- 
integer, dimension(3) :: ISIZE_READ, ISTART_READ, IEND_READ

real(kind=8) :: UNIF_VALUE

real(kind=8) :: XRAND

!- File descriptor
integer :: DESCRIPTEUR


!- Statistics of the droped field
real(kind=8), dimension(4) :: STATSCL

!- Size of the record
integer :: RECSIZE

!- Index
integer :: I, J, K

!---------------------------------------------------------------------

UNIF_VALUE = 1.0



!!====================================================================
!! 0. FLUID VELOCITY EQUAL TO ZERO
!!====================================================================
if(INIT_SCALAR <= 1) then

 THETA(:,:,:) = ZERO

 if(MYID==0) write(*,*)'Scalar initiation --> Theta(t=0)=0'



!!====================================================================
!! 2. RANDOM FLUID VELOCITY
!!====================================================================
elseif(INIT_SCALAR == 2) then

 do K = ISTART(3), IEND(3)
  do J = ISTART(2), IEND(2)
   do I = ISTART(1), IEND(1)

    call random_number(XRAND)
    THETA(I,J,K) = (0.5 - XRAND)*UNIF_VALUE

   end do
  end do
 end do

 if(MYID==0) write(*,*)'Scalar initiation --> Random field'



!!====================================================================
!! 3. Restart file
!!====================================================================
elseif(INIT_SCALAR == 3) then


!!--------------------------------------------------------------------
!! 3.1 Multiple Binary Files
!!--------------------------------------------------------------------
if(ISAVEFLUID == 1) then 

!!--------------------------------------------------------------------
!! 3.1.1. Fluid velocity
!!--------------------------------------------------------------------
 !- File name
 write(FILENAME,10101)'scl',trim(FILE_EXT),'.ini'

 !- Open file containing the initial scalar field
 open(unit = 150, file=trim(FILENAME), form='unformatted')

 !- Read size of stored file
 read(150)ISIZE_READ(1),ISTART_READ(1),IEND_READ(1)
 read(150)ISIZE_READ(2),ISTART_READ(2),IEND_READ(2)
 read(150)ISIZE_READ(3),ISTART_READ(3),IEND_READ(3)


 !- Check size of the stored field 
 if(    (ISIZE_READ(1)/=ISIZE(1)) &
    .or.(ISIZE_READ(2)/=ISIZE(2)) &
    .or.(ISIZE_READ(3)/=ISIZE(3)) ) then
  write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write(*,*)'!!                     ERROR                    !!'
  write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write(*,*)'!!'
  write(*,*)'!! file     : init_scalar.f90'
  write(*,*)'!!'
  write(*,*)'!! PROBLEM DESCRIPTION: '
  write(*,*)'!! initial scalar field has wrong dimensions' 
  write(*,*)'!!'
  write(*,*)'!! NX_READ=',ISIZE_READ(1),'  local NX=',ISIZE(1)
  write(*,*)'!! NY_READ=',ISIZE_READ(2),'  local NY=',ISIZE(2)
  write(*,*)'!! NZ_READ=',ISIZE_READ(3),'  local NZ=',ISIZE(3)
  write(*,*)'!!'
  write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  stop
 end if  

 read(150)(((THETA(I,J,K),I=ISTART(1),IEND(1)),J=ISTART(2),IEND(2)),K=ISTART(3),IEND(3))

 !- Close file
 close(150)

 if(MYID==0) write(*,*)'Scalar initiation --> Read from file'
 if(MYID==0) write(*,*)'       + Multiple Binary Files'





!!--------------------------------------------------------------------
!! 3.2 MPI I/O 
!!--------------------------------------------------------------------
elseif(ISAVEFLUID == 2) then

!- Scalar field
 FILENAME='scl.ini'
 call READ_MPIIO(THETA,FILENAME)

if(MYID==0) write(*,*)'Scalar initiation --> Read from file'
if(MYID==0) write(*,*)'       + MPI I/O'

end if

end if



if(MYID==0) write(*,*) 'Scalar initiation --> OK'


STATSCL(:) = ZERO

!!- 1st order moment
call RSUMCPU(SUM(THETA(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))   ),STATSCL(1))

!!- 2nd order moment
call RSUMCPU(SUM(THETA(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))**2),STATSCL(2))

!!- 3rd order moment
call RSUMCPU(SUM(THETA(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))**3),STATSCL(3))

!!- 4th order moment
call RSUMCPU(SUM(THETA(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))**4),STATSCL(4))


if(MYID==0) write(*,*) '      <T  > = ',STATSCL(1) / NGLOB
if(MYID==0) write(*,*) '      <T^2> = ',STATSCL(2) / NGLOB
if(MYID==0) write(*,*) '      <T^3> = ',STATSCL(3) / NGLOB
if(MYID==0) write(*,*) '      <T^4> = ',STATSCL(4) / NGLOB
if(MYID==0) write(*,*)


!!====================================================================
10200 format (A,A,A,I8.8,A)
10201 format (A,i1)
10202 format (A,i2)
10203 format (A,i3)
10204 format (A,i4)
10101 format (A,A,A)

10500 format (A,I3,A,I3,A,I3,A,I3,A,E13.6)

end subroutine INITIATION_SCALAR
