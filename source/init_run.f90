!!====================================================================
!!
!! This routine defines the constants for the simulation
!!
!!====================================================================

subroutine INIT_RUN

!!====================================================================
!!
!!
!!====================================================================

use DNS_DIM
use PARAM_PHYS 
use FORCING
use GEOMETRIC_VARIABLE
use ENSIGHT_VAR
use STATISTICS
use RANDOM_VAR


implicit none


!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------

integer :: I, N
!---------------------------------------------------------------------

!!====================================================================
!! 1. Basics
!!====================================================================


!! In case of (i.e. Stokes flow) scalar and DPS are switch off
if(SOLVE_FLUID == 2) then
 SOLVE_SCALAR = .false.
 SOLVE_PART = .false.

!!- no forcing because no turbulence !
 STEADY = .false.
end if


if(SOLVE_FLUID == 0) then

 SOLVE_SCALAR = .false.

!!- No-fluid -> no-stat
 LEVEL_STFLU = 0

!!- no forcing because no turbulence !
 STEADY = .false.
end if




!!- Ghost cells
NGHTCELL = 0
if(SOLVE_PART.and. (SOLVE_FLUID>0 .or. FROZEN_FLOW>0)) then
 NGHTCELL = 2
end if


!- File containing all informations
if(MYID==0) then
 UNIT_INFO(1) = 151
 UNIT_INFO(2) = 152
 UNIT_INFO(3) = 153
 UNIT_INFO(4) = 154
 UNIT_INFO(5) = 155

 !- Opening file
 open(unit=UNIT_INFO(1), file='run.info', action='write', status='replace')

 !- Opening file
 open(unit=UNIT_INFO(2), file='stat.info', action='write', status='replace')

 !- Opening file
 open(unit=UNIT_INFO(3), file='cpu.info', action='write', status='replace')


!- CPU loading  
if(SOLVE_PART) then
 open(unit=UNIT_INFO(4), file='partcpu_load.stat', action='write', status='replace')
 write(UNIT_INFO(4),*)'# t  [Np max/CPU] [Np min/CPU] [var(Np/CPU)]  '
end if


 !- Opening file
 open(unit=UNIT_INFO(5), file='numerics.info', action='write', status='replace')


end if



!!- Method for saving restart file
!!  This variable is filled in subroutine INIT_RUN
!! 
!! ISAVEFLUID = 1: Multiple binary files
!!            = 2: Using of MPI I/O
ISAVEFLUID = 2

!! 
!! ISAVEPART = 1: Multiple binary files
!!           = 2: Single binary file 
ISAVEPART = 2

!! ISAVEFORCE= 1: Multiple binary files
!!           = 2: Direct access file
!!           = 3: Direct access file without random seed generator
ISAVEFORCE = 2


!!====================================================================
!! 2. Physics
!!====================================================================
PPI = 4.*atan(1.D0) !- Pi
TWOPI = 2.*PPI      !- 2.*Pi

ICMPL = CMPLX(0.D0,1.D0) !!- complex number

KFIRST = TWOPI/LXMAX   !!- First wavenumber
KLAST = PPI/(LXMAX/NX) !!- Last wavenumber

!!------------------------------------------------------------------
!! 2.1 Time step initiation
!!------------------------------------------------------------------
!- Default time step is User's time step
DTIME = DTIME_USER


!!- Time step statistics
NEVEN_DT = 0
MAX_DT = ZERO
MIN_DT = INFINITY
MEAN_DT = ZERO






!!- Forced wavenumber
if(STEADY) then
 KFORCE_MIN = KFORCE_MIN*KFIRST
 KFORCE_MAX = KFORCE_MAX*KFIRST
end if



!!- Seed for random number generator
!! Must be the same for all process for the forcing
IDFORCE = -2
ISET = 0
IY = 0
IV(:) = 0


!!- Size of particle's array
if(SOLVE_PART)  then

NPCPU_UNIF = int(NPART_FULL/NPROC)



if(NPCPU_UNIF <= 100) then

! NPMAX_LOC = 2*NPART_FULL
 NPEXCH_MAX = NPART_FULL

elseif(NPCPU_UNIF <= 1000) then

! NPMAX_LOC = int(3*NPCPU_UNIF)
 NPEXCH_MAX = int(2*NPCPU_UNIF)

elseif(NPCPU_UNIF <= 10000) then

! NPMAX_LOC = int(2.5*NPCPU_UNIF)
 NPEXCH_MAX = int(1.5*NPCPU_UNIF)

elseif(NPCPU_UNIF <= 100000) then

! NPMAX_LOC = int(1.5*NPCPU_UNIF)
 NPEXCH_MAX = int(1.0*NPCPU_UNIF)

else

! NPMAX_LOC = int(1.5*NPCPU_UNIF)
 NPEXCH_MAX = int(0.75*NPCPU_UNIF)

end if

NPMAX_LOC = NPCPU_UNIF + NPEXCH_MAX




!!- Array size of echanged particles limited to 5%
!!NPEXCH_MAX = int(0.05*NPMAX_LOC)
!!NPEXCH_MAX = NPMAX_LOC - NPCPU_UNIF



if(MYID==0) write(*,*) ' NPCPU_UNIF =',NPCPU_UNIF
if(MYID==0) write(*,*) '  NPMAX_LOC =',NPMAX_LOC
if(MYID==0) write(*,*) ' NPEXCH_MAX =',NPEXCH_MAX

!!- Scalar case !!!AP - Attention, ce n'est plus vrai pour les collisions!!!!
!if(NPROC == 1) then
! NPMAX_LOC = NPART_FULL
! NPEXCH_MAX = 0
!end if

end if



!!FILTERING = .false.
!!KCUT = 2*KLAST
if(FILTERING) then
!! KCUT = 8*KFIRST
 if(MYID==0) write(*,*) ' Filtering activated'
 if(MYID==0) write(*,*) '     + Kcut/K0 =',KCUT/KFIRST
end if



!!====================================================================
!! 3. Outputing
!!====================================================================
ENSIGHT_BIN = .false.
ENSIGHT_DIM = 3
SLCNUM = NZ/2


NFILEOUT = 0

write(FILE_EXT,10205)'.p',MYID




if(ENSIGHT_OUT > 0) then
CASEFLU = 'fluid'
NBNODE = 2

allocate(NODELIST(NBNODE))
allocate(NODETYPE(NBNODE))

end if



!!====================================================================
!! 4. Statistics
!!====================================================================
if(.not.SOLVE_PART)  then
 LEVEL0_STPAR = .false.
 LEVEL1_STPAR = .false.
 LEVEL2_STPAR = .false.
end if

if(.not.SOLVE_SCALAR)  then
 LEVEL0_STSCL = .false.
 LEVEL1_STSCL = .false.
end if



if(LEVEL1_STPAR.or.LEVEL2_STPAR.or.LEVEL3_STPAR) then
 LEVEL0_STPAR=.true.
else
 LEVEL0_STPAR=.false.
end if

if(LEVEL1_STSCL.or.LEVEL2_STSCL) then
 LEVEL0_STSCL=.true.
else
 LEVEL0_STSCL=.false.
end if


if(SOLVE_FLUID==0) then
 LEVEL_STFLU = 0
 LEVEL0_STSCL = .false.
 
! LEVEL1_STPAR = .false.
! LEVEL2_STPAR = .false.
end if


if(LEVEL_STFLU>0.or.LEVEL0_STPAR.or.LEVEL0_STSCL) then
 NSTAT = 100
else
 NSTAT = 1
end if

!!====================================================================
!! Paticles (in this case Ellipsoidal)
!!====================================================================
! Semi Minor axes
!ELL_B = ELL_A/lambda
!ELL_C = ELL_B

! Moment of Inertia of the Ellipsoid
!IPXX = (ELL_B**2 + ELL_C**2)/5.0
!IPYY = (ELL_A**2 + ELL_C**2)/5.0
!IPZZ = (ELL_A**2 + ELL_B**2)/5.0


!- Maximum of polydisperse species
POLYDISPMAX = maxval(POLYDISP(:))

!!- Allocate particle density
allocate(RHOP(NIG,POLYDISPMAX))
!allocate(DPART(NIG,POLYDISPMAX))

allocate(EMAJ_PART(NIG, POLYDISPMAX))
allocate(IPXX(NIG, POLYDISPMAX))
allocate(IPYY(NIG, POLYDISPMAX))
allocate(IPZZ(NIG, POLYDISPMAX))

allocate(NPCLASS(NIG,POLYDISPMAX))
allocate(MPART(NIG, POLYDISPMAX))

!- Default is user's value
do I=1, NIG

    RHOP(I,:) = RHOP_USER(I)
    
    !DPART(I,:) = DPART_USER(I)
    EMAJ_PART(I,:) = EMAJ_PART_USER(I)

    NPCLASS(I,:) = NPART_FULL

end do 

!!====================================================================


!!--------------------------------------------------------------------
!! Lagrangian correlation
!!--------------------------------------------------------------------
!~ FOUT3 = 1
DIMLGR = 1

!!- Cycle where the Lagrangian function are computed
!! Be careful the Lagrangian stat will starts
!! at NT0(i)*FOUT2 (frequency of statistics ..)
!NT0(1) = 1
!NT0(2) = 10
!NT0(3) = 25
!NT0(4) = 50
!NT0(5) = 75
!NT0(6) = 100


NT0(1) = 1
NT0(2) = 200
NT0(3) = 400
NT0(4) = 600
NT0(5) = 800
NT0(6) = 1000
NT0(7) = 1200
NT0(8) = 1400
NT0(9) = 1600
NT0(10) = 1800
NT0(11) = 2000
NT0(12) = 2200
NT0(13) = 2400
NT0(14) = 2600
NT0(15) = 2800
NT0(16) = 3000


if(LEVEL2_STPAR) then

 DIMLGR = 1200

!!- Overestimation of length of lagrangian function
 DIMLGR = int(NCYCLEMAX/FOUT2)
 

 if(DIMLGR>NCYCLEMAX) DIMLGR = NCYCLEMAX


 if(MYID==0) write(*,*)'NBLGRMAX=',NBLGRMAX,' DIMLGR=',DIMLGR
 
!- Time of Lagrangian correlation
 allocate(TIME_LGR(DIMLGR,NBLGRMAX))
 TIME_LGR(:,:) = ZERO

 allocate(RPX_LOC(DIMLGR,NIG,NBLGRMAX,POLYDISPMAX))
 allocate(RPY_LOC(DIMLGR,NIG,NBLGRMAX,POLYDISPMAX))
 allocate(RPZ_LOC(DIMLGR,NIG,NBLGRMAX,POLYDISPMAX))
 !allocate(RFAPX_LOC(DIMLGR,NIG,NBLGRMAX,POLYDISPMAX))
 !allocate(RFAPY_LOC(DIMLGR,NIG,NBLGRMAX,POLYDISPMAX))
 !allocate(RFAPZ_LOC(DIMLGR,NIG,NBLGRMAX,POLYDISPMAX))

 RPX_LOC(:,:,:,:) = ZERO
 RPY_LOC(:,:,:,:) = ZERO
 RPZ_LOC(:,:,:,:) = ZERO

 !RFAPX_LOC(:,:,:,:) = ZERO
 !RFAPY_LOC(:,:,:,:) = ZERO
 !RFAPZ_LOC(:,:,:,:) = ZERO


end if





if(MYID==0 .and. DEBUG) then
 write(*,*)
 write(*,*)'-------------------'
 write(*,*)'--> DEBUG MODE <---'
 write(*,*)'-------------------'
 write(*,*)
end if



!!--------------------------------------------------------------------
10201 format(A,I1)
10202 format(A,I2)
10203 format(A,I3)
10204 format(A,I4)

10205 format(A,I4.4)

10600 format (A,I3,A,I4,A,E13.6)


end subroutine INIT_RUN
