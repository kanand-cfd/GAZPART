!!====================================================================
!!
!!  Initiation of Forcing scheme for statistically steady 
!!  turbulent flows
!!
!!====================================================================

subroutine INITIATION_FORCING_NEW

!!====================================================================
!! This routine works with ADD_FORCING
!!--------------------------------------------------------------------
!! The forcing scheme, proposed by Eswaran & Pope, J. Comp. & Fluids
!! (1988) is based on a stochastic force added at low-wavenumbers.
!!
!! The stochastic force is given by a Ornstein-Uhlenbeck process,
!! parameterized by a timescale (TIME_FORCE) and a variancce
!! (SIGMA_FORCE). The range of modified wavenumber is controled by
!! KFORCE_MIN and KFORCE_MAX.
!!
!! 
!!   KFORCE_MIN: minimum forced wavenumber
!!   KFORCE_MAX: maximum forced wavenumber
!!   TIME_FORCE: timescale of forcing
!!   SIGMA_FORCE: variance of forcing
!!
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE
use PARAM_PHYS
use FORCING
use RANDOM_VAR

use P3DFFT

implicit none



!!--------------------------------------------------------------------
!! ARRAYS STATEMENT
!!--------------------------------------------------------------------

!- Wavenumber
real(kind=8) :: KXFULL, KYFULL, KZFULL, KF

!- 
integer :: NF, NFM

!- Record size
integer :: RECSIZE

!- File name 
character(len=40) :: FILENAME


!- Statistics of the droped field
real(kind=8) :: STATX, STATY, STATZ

!- index
integer :: I, J, K
!---------------------------------------------------------------------


NFORCE_CPU = 0

do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)

   KF = (KX(I)**2 + KY(J)**2 + KZ(K)**2)**0.5

   if ((KF>=KFORCE_MIN).and.(KF<=KFORCE_MAX)) then
    NFORCE_CPU = NFORCE_CPU + 1
   end if

  end do
 end do
end do



!!- Sum over the domain
call ISUMCPU(NFORCE_CPU,NFORCE_FULL)


if(MYID==0) write(*,10800) 'Forcing initiation --> Nwave= ',NFORCE_FULL


!!- allocation arrays containing stochastic acceleration
allocate(FORCING_UFOU(NFORCE_FULL))
allocate(FORCING_VFOU(NFORCE_FULL))
allocate(FORCING_WFOU(NFORCE_FULL))

allocate(IFORCING(NFORCE_FULL))
allocate(JFORCING(NFORCE_FULL))
allocate(KFORCING(NFORCE_FULL))


NF = 0
do K = 1, NZ
 do J = 1, NY
  do I = 1, NX/2+1

!!- Wavenumbers
   KXFULL = (I-1)*KFIRST

   if(J<=NY/2+1) then
    KYFULL = (J-1)*KFIRST
   else
    KYFULL = (J-1-NY)*KFIRST
   end if
   if(K<=NZ/2+1) then
    KZFULL = (K-1)*KFIRST
   else
    KZFULL = (K-1-NZ)*KFIRST
   end if

!! KF = (KX(I)**2 + KY_FULL(J)**2 + KZ_FULL(K)**2)**0.5
   KF = (KXFULL**2 + KYFULL**2 + KZFULL**2)**0.5

   if ((KF>=KFORCE_MIN).and.(KF<=KFORCE_MAX)) then

    NF = NF + 1

    if(     (I>=FSTART(1)).and.(I<=FEND(1)) &
       .and.(J>=FSTART(2)).and.(J<=FEND(2)) &
       .and.(K>=FSTART(3)).and.(K<=FEND(3)) ) then
    
      IFORCING(NF) = I
      JFORCING(NF) = J
      KFORCING(NF) = K
    else
      IFORCING(NF) = 0
      JFORCING(NF) = 0
      KFORCING(NF) = 0
    end if

   end if
  end do
 end do
end do 


!- Barrier
 call MPI_BARRIER(MPI_COMM_WORLD,IERR)

!do NF=1, NFORCE_FULL
!write(870+MYID,*)NF,IFORCING(NF),JFORCING(NF),KFORCING(NF)
!end do

if(MYID==0) then

!!====================================================================
!! 
!!====================================================================
 if(INIT_FLUID_VELOCITY <= 2 .or.(INIT_FLUID_VELOCITY == 4)) then

  FORCING_UFOU(:) = cmplx(ZERO,ZERO)
  FORCING_VFOU(:) = cmplx(ZERO,ZERO)
  FORCING_WFOU(:) = cmplx(ZERO,ZERO)

  if(MYID==0) write(*,*)'Forcing coefficients initiation --> F=0'


!!====================================================================
!! 2. Read from file
!!====================================================================
 elseif(INIT_FLUID_VELOCITY == 3) then

  !- Define file name
  FILENAME='forcing_new.ini'

  !- Open file containing the forcing coefficients
  open(unit =120, file=trim(FILENAME), form='unformatted')

  !- Number of forced waves
  read(120) NFORCE_FULL

  !- Random seed
  read(120) IDFORCE
  read(120) ISET, GSET
  read(120) XRANDF, YRANDF, ZRANDF

  !- Forcing coefficients
  read(120)(FORCING_UFOU(I),I=1,NFORCE_FULL)
  read(120)(FORCING_VFOU(I),I=1,NFORCE_FULL)
  read(120)(FORCING_WFOU(I),I=1,NFORCE_FULL)

  !- Close file
  close(120)

  if(MYID==0) write(*,*) 'Forcing coefficients initiation --> Read from file'


 end if !- end if INIT_FLUID_VELOCITY == XXX

 STATX = ZERO
 STATY = ZERO
 STATZ = ZERO

 do I=1, NFORCE_FULL
  STATX = STATX + real(FORCING_UFOU(I)*conjg(FORCING_UFOU(I)))
  STATY = STATY + real(FORCING_VFOU(I)*conjg(FORCING_VFOU(I)))
  STATZ = STATZ + real(FORCING_WFOU(I)*conjg(FORCING_WFOU(I)))
 end do

 if(MYID==0) write(*,*) '       <fx> = ',STATX / real(NFORCE_FULL)
 if(MYID==0) write(*,*) '       <fy> = ',STATY / real(NFORCE_FULL)
 if(MYID==0) write(*,*) '       <fz> = ',STATZ / real(NFORCE_FULL)


end if !- end if MYID ==0


!- Broadcast forcing coefficient to all tasks
 call MPI_BCAST(FORCING_UFOU,NFORCE_FULL,MPI_COMPLEX16,0,MPI_COMM_WORLD,IERR)
 call MPI_BCAST(FORCING_VFOU,NFORCE_FULL,MPI_COMPLEX16,0,MPI_COMM_WORLD,IERR)
 call MPI_BCAST(FORCING_WFOU,NFORCE_FULL,MPI_COMPLEX16,0,MPI_COMM_WORLD,IERR)




!!--------------------------------------------------------------------
10800 format (1x,A,I3,A)
10101 format (A,A,A)


end subroutine INITIATION_FORCING_NEW

