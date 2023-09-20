!!====================================================================
!!
!!
!!====================================================================

subroutine INITIATION_FLUID

!!====================================================================
!!
!!
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE
use FLUID_VARIABLE
use PARAM_PHYS
use MPI_STRUCTURES

implicit none


!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------
real(kind=8) :: X0, Y0, Z0
real(kind=8) :: RADIUS, U0, V0, W0, OMEGA, PSI

!- random number
real(kind=8) :: XRAND

!- File name 
character(len=40) :: FILENAME

!- 
integer, dimension(3) :: ISIZE_READ, ISTART_READ, IEND_READ

!!- Read from NTMIX
real(kind=8) :: UNTMIX, LNTMIX, TNTMIX
real(kind=8) :: RDUMMY

!- Statistics of the droped field
real(kind=8), dimension(4) :: STATX, STATY, STATZ

!- File descriptor
integer :: DESCRIPTEUR

integer :: RECSIZE

!- Index
integer :: I, J, K, IJK

!---------------------------------------------------------------------

UREF(:) = ZERO

!!====================================================================
!! 0. FLUID VELOCITY EQUAL TO ZERO
!!====================================================================
if(INIT_FLUID_VELOCITY == 0) then

 UFLU(:,:,:) = ZERO 
 VFLU(:,:,:) = ZERO
 WFLU(:,:,:) = ZERO

 if(MYID==0) write(*,*)'Fluid velocity initiation: Uf(t=0)=0 --> OK'


!!====================================================================
!! 1. SINGLE EDDY
!!====================================================================
elseif(INIT_FLUID_VELOCITY == 1) then

 RADIUS = LXMAX/10.
 OMEGA = 0.1 !- 1/s

 U0 = 1.
 V0 = ZERO
 W0 = ZERO

 X0 = LXMAX/2.
 Y0 = LYMAX/2.
 Z0 = ZERO

 do K = ISTART(3), IEND(3)
  do J = ISTART(2), IEND(2)
   do I = ISTART(1), IEND(1)

    PSI = OMEGA*exp(-((XMESH(I)-X0)**2+(YMESH(J)-Y0)**2)/(2.*RADIUS**2))

    UFLU(I,J,K) = U0 - (YMESH(J)-Y0)/RADIUS**2*PSI
    VFLU(I,J,K) = V0 + (XMESH(I)-X0)/RADIUS**2*PSI
    WFLU(I,J,K) = ZERO

!    RADIUS = sqrt((XMESH(I)-LXMAX/2)**2 + (YMESH(J)-LYMAX/2)**2)
!    call random_number(XRAND)
!    if(RADIUS < LXMAX/5) then
!      UFLU(I,J,K) = 10. + XRAND
!    else
!      UFLU(I,J,K) = 0.
!    end if     


   end do
  end do
 end do

 UREF(1) = U0
 UREF(2) = V0
 UREF(3) = W0


 if(MYID==0) write(*,*)'Fluid velocity initiation: 2d eddy --> OK'



!!====================================================================
!! 2. RANDOM FLUID VELOCITY
!!====================================================================
elseif(INIT_FLUID_VELOCITY == 2) then

 do K = ISTART(3), IEND(3)
  do J = ISTART(2), IEND(2)
   do I = ISTART(1), IEND(1)

    call random_number(XRAND)
    UFLU(I,J,K) = 0.5 - XRAND

    call random_number(XRAND)
    VFLU(I,J,K) = 0.5 - XRAND

    call random_number(XRAND)
    WFLU(I,J,K) = 0.5 - XRAND

   end do
  end do
 end do


 if(MYID==0) write(*,*)'Fluid velocity initiation: Random --> OK'


!!====================================================================
!! 3. Restart file
!!====================================================================
elseif(INIT_FLUID_VELOCITY >= 3) then

!!--------------------------------------------------------------------
!! 3.1. Multiple Binary Files
!!--------------------------------------------------------------------
if(ISAVEFLUID == 1) then 

!!--------------------------------------------------------------------
!!- Uf
!!--------------------------------------------------------------------
!- Define file name
write(FILENAME,10101)'uf',trim(FILE_EXT),'.ini'

!- Open file containing the initial fluid velocity field
open(unit = 120, file=trim(FILENAME), form='unformatted')

!- Read size of stored file
read(120)ISIZE_READ(1),ISTART_READ(1),IEND_READ(1)
read(120)ISIZE_READ(2),ISTART_READ(2),IEND_READ(2)
read(120)ISIZE_READ(3),ISTART_READ(3),IEND_READ(3)

!- Check size of the stored field 
if(	(ISIZE_READ(1)/=ISIZE(1)) &
   .or.(ISIZE_READ(2)/=ISIZE(2)) &
   .or.(ISIZE_READ(3)/=ISIZE(3)) ) then
 write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 write(*,*)'!! 		    ERROR		     !!'
 write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 write(*,*)'!!'
 write(*,*)'!! file	 : init_fluid.f90'
 write(*,*)'!!'
 write(*,*)'!! PROBLEM DESCRIPTION: '
 write(*,*)'!! initial fluid velocity field has'
 write(*,*)'!! wrong dimensions' 
 write(*,*)'!!'
 write(*,*)'!! NX_READ=',ISIZE_READ(1),'  local NX=',ISIZE(1)
 write(*,*)'!! NY_READ=',ISIZE_READ(2),'  local NY=',ISIZE(2)
 write(*,*)'!! NZ_READ=',ISIZE_READ(3),'  local NZ=',ISIZE(3)
 write(*,*)'!!'
 write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 stop
end if  

read(120)(((UFLU(I,J,K),I=ISTART(1),IEND(1)),J=ISTART(2),IEND(2)),K=ISTART(3),IEND(3))

!- Close file
close(120)



!!--------------------------------------------------------------------
!!- Vf
!!--------------------------------------------------------------------
!- Define file name
write(FILENAME,10101)'vf',trim(FILE_EXT),'.ini'

!- Open file containing the initial fluid velocity field
open(unit = 120, file=trim(FILENAME), form='unformatted')

!- Read size of stored file
read(120)ISIZE_READ(1),ISTART_READ(1),IEND_READ(1)
read(120)ISIZE_READ(2),ISTART_READ(2),IEND_READ(2)
read(120)ISIZE_READ(3),ISTART_READ(3),IEND_READ(3)

!- Check size of the stored field 
if(	(ISIZE_READ(1)/=ISIZE(1)) &
   .or.(ISIZE_READ(2)/=ISIZE(2)) &
   .or.(ISIZE_READ(3)/=ISIZE(3)) ) then
 write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 write(*,*)'!! 		    ERROR		     !!'
 write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 write(*,*)'!!'
 write(*,*)'!! file	 : init_fluid.f90'
 write(*,*)'!!'
 write(*,*)'!! PROBLEM DESCRIPTION: '
 write(*,*)'!! initial fluid velocity field has'
 write(*,*)'!! wrong dimensions' 
 write(*,*)'!!'
 write(*,*)'!! NX_READ=',ISIZE_READ(1),'  local NX=',ISIZE(1)
 write(*,*)'!! NY_READ=',ISIZE_READ(2),'  local NY=',ISIZE(2)
 write(*,*)'!! NZ_READ=',ISIZE_READ(3),'  local NZ=',ISIZE(3)
 write(*,*)'!!'
 write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 stop
end if  

read(120)(((VFLU(I,J,K),I=ISTART(1),IEND(1)),J=ISTART(2),IEND(2)),K=ISTART(3),IEND(3))

!- Close file
close(120)


!!--------------------------------------------------------------------
!!- Wf
!!--------------------------------------------------------------------
!- Define file name
write(FILENAME,10101)'wf',trim(FILE_EXT),'.ini'

!- Open file containing the initial fluid velocity field
open(unit = 120, file=trim(FILENAME), form='unformatted')

!- Read size of stored file
read(120)ISIZE_READ(1),ISTART_READ(1),IEND_READ(1)
read(120)ISIZE_READ(2),ISTART_READ(2),IEND_READ(2)
read(120)ISIZE_READ(3),ISTART_READ(3),IEND_READ(3)

!- Check size of the stored field 
if(	(ISIZE_READ(1)/=ISIZE(1)) &
   .or.(ISIZE_READ(2)/=ISIZE(2)) &
   .or.(ISIZE_READ(3)/=ISIZE(3)) ) then
 write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 write(*,*)'!! 		    ERROR		     !!'
 write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 write(*,*)'!!'
 write(*,*)'!! file	 : init_fluid.f90'
 write(*,*)'!!'
 write(*,*)'!! PROBLEM DESCRIPTION: '
 write(*,*)'!! initial fluid velocity field has'
 write(*,*)'!! wrong dimensions' 
 write(*,*)'!!'
 write(*,*)'!! NX_READ=',ISIZE_READ(1),'  local NX=',ISIZE(1)
 write(*,*)'!! NY_READ=',ISIZE_READ(2),'  local NY=',ISIZE(2)
 write(*,*)'!! NZ_READ=',ISIZE_READ(3),'  local NZ=',ISIZE(3)
 write(*,*)'!!'
 write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 stop
end if  

read(120)(((WFLU(I,J,K),I=ISTART(1),IEND(1)),J=ISTART(2),IEND(2)),K=ISTART(3),IEND(3))

!- Close file
close(120)

!!--------------------------------------------------------------------
!! 3.2 MPI I/O 
!!--------------------------------------------------------------------
elseif(ISAVEFLUID == 2) then

!!--------------------------------------------------------------------
!!- Uf
!!--------------------------------------------------------------------
 FILENAME='uf.ini'
 call READ_MPIIO(UFLU,FILENAME)
 
!!--------------------------------------------------------------------
!!- Wf
!!--------------------------------------------------------------------
 FILENAME='vf.ini'
 call READ_MPIIO(VFLU,FILENAME)
 
!!--------------------------------------------------------------------
!!- Wf
!!--------------------------------------------------------------------
 FILENAME='wf.ini'
 call READ_MPIIO(WFLU,FILENAME)
 
 if(MYID==0) write(*,*)'Fluid velocity initiation: Read from file --> OK'
 
!!--------------------------------------------------------------------
!! 3.3. Statistics of the read file
!!--------------------------------------------------------------------

STATX(:) = ZERO
STATY(:) = ZERO
STATZ(:) = ZERO

!!- 1st order moment
call RSUMCPU(SUM(UFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))),STATX(1))
call RSUMCPU(SUM(VFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))),STATY(1))
call RSUMCPU(SUM(WFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))),STATZ(1))

!!- 2nd order moment
call RSUMCPU(SUM(UFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))**2),STATX(2))
call RSUMCPU(SUM(VFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))**2),STATY(2))
call RSUMCPU(SUM(WFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))**2),STATZ(2))

!!- 3rd order moment
call RSUMCPU(SUM(UFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))**3),STATX(3))
call RSUMCPU(SUM(VFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))**3),STATY(3))
call RSUMCPU(SUM(WFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))**3),STATZ(3))

!!- 4th order moment
call RSUMCPU(SUM(UFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))**4),STATX(4))
call RSUMCPU(SUM(VFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))**4),STATY(4))
call RSUMCPU(SUM(WFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))**4),STATZ(4))


if(MYID==0) write(*,*) '       <Ux> = ',STATX(1) / NGLOB
if(MYID==0) write(*,*) '       <Uy> = ',STATY(1) / NGLOB
if(MYID==0) write(*,*) '       <Uz> = ',STATZ(1) / NGLOB
if(MYID==0) write(*,*) '     <Ux^2> = ',STATX(2) / NGLOB
if(MYID==0) write(*,*) '     <Uy^2> = ',STATY(2) / NGLOB
if(MYID==0) write(*,*) '     <Uz^2> = ',STATZ(2) / NGLOB
if(MYID==0) write(*,*) '     <Ux^3> = ',STATX(3) / NGLOB
if(MYID==0) write(*,*) '     <Uy^3> = ',STATY(3) / NGLOB
if(MYID==0) write(*,*) '     <Uz^3> = ',STATZ(3) / NGLOB
if(MYID==0) write(*,*) '     <Ux^4> = ',STATX(4) / NGLOB
if(MYID==0) write(*,*) '     <Uy^4> = ',STATY(4) / NGLOB
if(MYID==0) write(*,*) '     <Uz^4> = ',STATZ(4) / NGLOB
if(MYID==0) write(*,*) 




end if

end if



if(MYID==0) write(*,*) 'Fluid initiation --> OK'


!!====================================================================
10201 format (A,I1)
10202 format (A,I2)
10203 format (A,I3)
10204 format (A,I4)
10101 format (A,A,A)

end subroutine INITIATION_FLUID
