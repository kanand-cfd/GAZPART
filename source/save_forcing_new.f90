!!====================================================================
!!
!!                   Save forcing coefficient
!!
!!====================================================================

subroutine SAVE_FORCING_NEW

!!====================================================================
!!
!! ISAVEFLUID = 1 : Multiple binary files
!!            = 2 : MPI I/O
!!
!!====================================================================

use DNS_DIM
use PARAM_PHYS
use FORCING
use RANDOM_VAR
use MPI_STRUCTURES


implicit none


!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- File name 
character(len=40) :: FILENAME

!- Record size
integer :: RECSIZE

!- Wavenumber modulus
real(kind=8) :: KF

!- Statistics of the droped field
real(kind=8) :: STATX , STATY , STATZ

!- Index
integer :: NF
integer :: I, J ,K
!---------------------------------------------------------------------


if(MYID ==0) then

 !- Define file name
 FILENAME='forcing_new.end'

 !- Open file containing the forcing
 open(unit = 120, file = trim(FILENAME),form='unformatted')

 !- Number of forced waves
 write(120) NFORCE_FULL

 !- Random seed
 write(120) IDFORCE
 write(120) ISET, GSET
 write(120) XRANDF, YRANDF, ZRANDF

 !- Forcing coefficients
 write(120)(FORCING_UFOU(I),I=1,NFORCE_FULL)
 write(120)(FORCING_VFOU(I),I=1,NFORCE_FULL)
 write(120)(FORCING_WFOU(I),I=1,NFORCE_FULL)

 !- Close file
 close(120)


!!====================================================================
!! 3. Satistics of the forcing
!!====================================================================

 STATX = ZERO
 STATY = ZERO
 STATZ = ZERO

 do I=1,NFORCE_FULL
  STATX = STATX + real(FORCING_UFOU(I)*conjg(FORCING_UFOU(I)))
  STATY = STATY + real(FORCING_VFOU(I)*conjg(FORCING_VFOU(I)))
  STATZ = STATZ + real(FORCING_WFOU(I)*conjg(FORCING_WFOU(I)))
 end do


 if(MYID==0) write(*,*) 'Final forcing coefficients --> Saved'

 if(MYID==0) write(*,*) '       <fx> = ',STATX / real(NFORCE_FULL)
 if(MYID==0) write(*,*) '       <fy> = ',STATY / real(NFORCE_FULL)
 if(MYID==0) write(*,*) '       <fz> = ',STATZ / real(NFORCE_FULL)



end if !- end if MYID == 0


!!====================================================================
10100 format (A,A,A)
10200 format (A,A,A,I8.8,A)
10300 format (A,A)
10400 format (A,F8.2,A)


end subroutine SAVE_FORCING_NEW
