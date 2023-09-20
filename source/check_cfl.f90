!!====================================================================
!! 
!! Check CFL condition for time-integration stability
!!
!!====================================================================

subroutine CHECK_CFL(NCYCLE)

!!====================================================================
!!
!!
!!====================================================================

use DNS_DIM
use FLUID_VARIABLE
use GEOMETRIC_VARIABLE, only: DX, DY, DZ
use PARAM_PHYS,         only: DTIME, DTIME_USER, VISC, FOUT0
use STATISTICS,         only: MEAN_DT, MAX_DT, MIN_DT, NEVEN_DT
use PARTICLE_PARALLEL

implicit none



!------------------------------------------------------------------
! ARRAYS STATEMENT
!------------------------------------------------------------------
integer, intent(in) :: NCYCLE

!- Maximum velocity
real(kind=8) :: UMAX, VMAX, WMAX

!- Time step according to CFL condition
real(kind=8) :: DT_FOURIER
real(kind=8) :: CFLMAX, FOUMAX
real(kind=8) :: DT_COURANT,DT_UP, DT_VP, DT_WP

real(kind=8), dimension(NIG) :: MINTAUP 

real(kind=8) :: DTIME_MIN


integer :: IFLAG1

integer :: I, J

!-----------------------------------------------------------------
DTIME_MIN = DTIME


CFLMAX = 0.05
FOUMAX = 1.0 

!!====================================================================
!! 1. Time step limitationfor the fluid
!!====================================================================
UMAX = ZERO
VMAX = ZERO
WMAX = ZERO

if(NPROC>1) then

 call RMAXCPU(maxval(abs(UFLU)),UMAX)
 call RMAXCPU(maxval(abs(VFLU)),VMAX)
 call RMAXCPU(maxval(abs(WFLU)),WMAX)

else

 UMAX=maxval(abs(UFLU))
 VMAX=maxval(abs(VFLU))
 WMAX=maxval(abs(WFLU))

end if

if(UMAX>ZERO) then
 DT_UP  = CFLMAX*DX/UMAX
else
 DT_UP = INFINITY
end if
if(VMAX>ZERO) then
 DT_VP  = CFLMAX*DY/VMAX
else
 DT_VP = INFINITY
end if
if(WMAX>ZERO) then
 DT_WP  = CFLMAX*DZ/WMAX
else
 DT_WP = INFINITY
end if



DT_COURANT = INFINITY
DT_COURANT = min(DT_COURANT,DT_UP)
DT_COURANT = min(DT_COURANT,DT_VP)
DT_COURANT = min(DT_COURANT,DT_WP)

DT_FOURIER = INFINITY
DT_FOURIER = min(DT_FOURIER,FOUMAX*DX**2/VISC)
DT_FOURIER = min(DT_FOURIER,FOUMAX*DY**2/VISC)
DT_FOURIER = min(DT_FOURIER,FOUMAX*DZ**2/VISC)


!!====================================================================
!! 2. Time step limitationfor the particles
!!====================================================================
!INVTAUPMAX = ZERO
!do J=1, NIG
!
! if(PARTDEF(J) == 2) then
!
!  if(NPROC>1) then
!   call RMAXCPU(PART(1:NPART_LOC(J),J)%INVTAUP,RDUMMY)
!  else
!   RDUMMY=maxval(PART(1:NPART_LOC(J),J)%INVTAUP)
!  end if
!
!  INVTAUPMAX = max(RDUMMY,INVTAUPMAX)
!
!
! end if !- end if on PARTDEF
!
!end do






DTIME = min(DT_COURANT,DT_FOURIER,DTIME_USER)

!if((mod(NCYCLE,FOUT0) == 0).and.(MYID==0)) then
! write(*,10603)'Dt_Cou=',DT_COURANT,' Dt_Fou=',DT_FOURIER,' Dt=',DTIME
!end if

NEVEN_DT = NEVEN_DT + 1.0
MAX_DT = max(DTIME,MAX_DT)
MIN_DT = min(DTIME,MIN_DT)
MEAN_DT = MEAN_DT + DTIME




!!======================================================================
!! Print in run.info
!!======================================================================
if((mod(NCYCLE,FOUT0) == 0).and.(MYID==0)) then
 IFLAG1 = 5

 write(UNIT_INFO(IFLAG1),*)'-----------------------'
 write(UNIT_INFO(IFLAG1),*)'Ncycle=',NCYCLE
 write(UNIT_INFO(IFLAG1),10602)'max(|uf|) =',UMAX,' Courant =',UMAX*DTIME/DX
 write(UNIT_INFO(IFLAG1),10602)'max(|vf|) =',VMAX,' Courant =',VMAX*DTIME/DY
 write(UNIT_INFO(IFLAG1),10602)'max(|wf|) =',WMAX,' Courant =',WMAX*DTIME/DZ
 write(UNIT_INFO(IFLAG1),10601)'Dt = ',DTIME
 write(UNIT_INFO(IFLAG1),10601)'Dt/Dt(Courant) =',DTIME/DT_COURANT
 write(UNIT_INFO(IFLAG1),10601)'Dt/Dt(Fourier) =',DTIME/DT_FOURIER
 write(UNIT_INFO(IFLAG1),10601)'Dt/Dt(User   ) =',DTIME/DTIME_USER
 write(UNIT_INFO(IFLAG1),*)
end if

!-----------------------------------------------------------------------
10601 format (2x,A,1x,E13.6)
10602 format (2x,A,1x,E13.6,A,2x,E13.6)
10603 format (2x,A,1x,E13.6,A,2x,E13.6,A,2x,E13.6)

end subroutine CHECK_CFL
