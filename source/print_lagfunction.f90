!!====================================================================
!!
!!
!!====================================================================
!!
subroutine PRINT_LAGFUNCTION
!!
!!====================================================================
!! The normalized autocorrelation function is computed as
!!
!!                 <u(t0).u(t0+t)>
!! R(t) = ---------------------------------
!!        sqrt(<u(t0)^2>).sqrt(<u(t0+t)^2>)
!!
!!====================================================================

use dns_dim
use param_phys
use particle_parallel
use statistics

implicit none

!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!- Local arrays
!---------------------------------------------------------------------
!- Lagrangian function
real(kind=8), dimension(DIMLGR) ::   RPX,   RPY,   RPZ
real(kind=8), dimension(DIMLGR) :: RFAPX, RFAPY, RFAPZ

!- True time step (accouting for frequency)
real(kind=8) :: DTIME0

!- File name
character(len=40) :: FILENAME

integer :: NUNIT

!- Index
integer :: I, J, LGR, K, IDP
!---------------------------------------------------------------------

if(MYID==0)write(*,*)' Lag. Func.: NBLGRMAX =',NBLGRMAX


DTIME0 = DTIME*FOUT2

!!- loop on lagrangian function
do K =1, NBLGRMAX


if(MYID==0)write(*,10700)' Lag. Func. #',K,' starts at N=',NT0(K)*FOUT2,' Time =',TIME_LGR(1,K)


 !!- Loop on particle kind
 do J = 1, NIG
  do IDP = 1, POLYDISP(J)


!!======================================================================
!! 1. Add all Lagrangian function from each CPU
!!======================================================================
 call MPI_ALLREDUCE(RPX_LOC(:,J,K,IDP),RPX,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
 call MPI_ALLREDUCE(RPY_LOC(:,J,K,IDP),RPY,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
 call MPI_ALLREDUCE(RPZ_LOC(:,J,K,IDP),RPZ,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

 !call MPI_ALLREDUCE(RFAPX_LOC(:,J,K,IDP),RFAPX,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
 !call MPI_ALLREDUCE(RFAPY_LOC(:,J,K,IDP),RFAPY,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
 !call MPI_ALLREDUCE(RFAPZ_LOC(:,J,K,IDP),RFAPZ,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

!!- Normalization
! RPX = RPX / real(NPROC)
! RPY = RPY / real(NPROC)
! RPZ = RPZ / real(NPROC)

! RFAPX = RFAPX / real(NPROC)
! RFAPY = RFAPY / real(NPROC)
! RFAPZ = RFAPZ / real(NPROC)
 

!!======================================================================
!! 2. Print in file
!!======================================================================
 if(MYID==0) then

!!----------------------------------------------------------------------
!!- 3.1. Define filename
!!----------------------------------------------------------------------
   write(FILENAME,10600)'part_l2_Rp_p',J,'_c',IDP,'_f',K,'.stat'
   open(unit=600, file=trim(FILENAME), status='replace')
   write(600,20104)
!!----------------------------------------------------------------------
!!- 3.2. Print in file
!!----------------------------------------------------------------------
   do I = 1,DIMLGR-1
      write(600,10000) TIME_LGR(I,K),  &
                        RPX(I), &
                        RPY(I), &
                        RPZ(I), &
                        (RPX(I)+RPY(I)+RPZ(I))/3.0D0 !, &
                      !RFAPX(I), &
                      !RFAPY(I), &
                      !RFAPZ(I), &
                      !(RFAPX(I)+RFAPY(I)+RFAPZ(I))/3.0D0
   end do 

 close(600)

 
 end if !- If: (MYID==0)

  end do !!- Loop: IDP=1,

end do  !!- Loop: J = 1, NIG

end do !!- Loop: K= 1, NBLGR




!!--------------------------------------------------------------------
10600 format (A,I2.2,A,I1,A,I2.2,A)
10700 format (A,I2.2,A,I6.6,A,F12.3)
20104 format('# t, Rpx, Rpy, Rpz, Rp, Rf@px, Rf@py, Rf@pz, Rf@p')


10000 format (30(e17.7))
10001 format (50(e17.7))

end subroutine PRINT_LAGFUNCTION 
