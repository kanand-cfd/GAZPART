!!====================================================================
!!
!!                   Save forcing coefficient
!!
!!====================================================================

subroutine SAVE_FORCING

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

!!====================================================================
!! 1. Save with multiple binary
!!====================================================================
 if(ISAVEFORCE==1) then

 !- Define file name
 FILENAME='forcing.end'

 write(FILENAME,10100)'forcing',trim(FILE_EXT),'.end'

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
!! 2. Save direct access file
!!====================================================================
 elseif(ISAVEFORCE==2) then

  !- RECSIZE = 8*2*NFORCE_FULL
  RECSIZE = (4 + 4 + 8 + 4 + 32*4) + 8*2*NFORCE_FULL

  FILENAME= 'forcing.end'

  !- Open files
  open(unit=120,file=trim(FILENAME),access='direct',recl=RECSIZE,form='unformatted'   )

  !- Print data
  write(unit=120,rec=3*MYID+1)IDFORCE,ISET,GSET,IY,IV,FORCING_UFOU
  write(unit=120,rec=3*MYID+2)IDFORCE,ISET,GSET,IY,IV,FORCING_VFOU
  write(unit=120,rec=3*MYID+3)IDFORCE,ISET,GSET,IY,IV,FORCING_WFOU

  close(120)

!!====================================================================
!! 3. Save direct access file without random seed
!!====================================================================
 elseif (ISAVEFORCE==3) then

  RECSIZE = 8*2*3

  FILENAME= 'forcing.end'

  !- Open files
  open(unit=120,file=trim(FILENAME),access='direct',recl=RECSIZE,form='unformatted'   )

  !- Print data
  NF = 0
  do K = FSTART(3), FEND(3)
   do J = FSTART(2), FEND(2)
    do I = FSTART(1), FEND(1)

     KF = (KX(I)**2 + KY(J)**2 + KZ(K)**2)**0.5

     if ((KF>=KFORCE_MIN).and.(KF<=KFORCE_MAX)) then
      NF = NF +1
    
      write(unit=120,rec=NF+1)FORCING_UFOU(NF), FORCING_VFOU(NF), FORCING_WFOU(NF)

     end if

    end do
   end do
  end do

  close(120)

 end if



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

 if(MYID==0) write(*,*) '       <fx> = ',STATX / NFORCE_FULL
 if(MYID==0) write(*,*) '       <fy> = ',STATY / NFORCE_FULL
 if(MYID==0) write(*,*) '       <fz> = ',STATZ / NFORCE_FULL



end if !- end if MYID == 0


!!====================================================================
10100 format (A,A,A)
10200 format (A,A,A,I8.8,A)
10300 format (A,A)
10400 format (A,F8.2,A)


end subroutine SAVE_FORCING
