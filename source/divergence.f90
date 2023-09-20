!!====================================================================
!!
!!
!!====================================================================

subroutine DIVERGENCE(DIVMAX,DIVMEAN,DIVMEAN2)

!!====================================================================
!!
!!
!!====================================================================

use dns_dim
use fluid_variable     !- Fluid velocity
use geometric_variable
use param_phys
use work_arrays

use p3dfft

implicit none


!------------------------------------------------------------------
! ARRAYS STATEMENT
!------------------------------------------------------------------
!- Global arrays
!---------------

!- Mean of divergence
real(kind=8), intent(out) :: DIVMEAN

!- Mean of absolute divergence
real(kind=8), intent(out) :: DIVMEAN2

!- Maximum of divergence
real(kind=8), intent(out) :: DIVMAX

!---------------------------------------------------------------------
!- Local arrays
!---------------------------------------------------------------------

!- Index
integer :: I, J, K, IJK
!---------------------------------------------------------------------


DIVMEAN2 = 0.
DIVMEAN = 0.
DIVMAX = 0.

!!====================================================================
!! 1. Compute divergence in Fourier sapce
!!====================================================================
do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)

   TMPFOU(I,J,K) = ICMPL*KX(I)*UFOU(I,J,K) &
                 + ICMPL*KY(J)*VFOU(I,J,K) &
                 + ICMPL*KZ(K)*WFOU(I,J,K)
  end do
 end do
end do


!!====================================================================
!! 2. Transform the divergence in Physical space
!!====================================================================
!- Synchronize all the process
call MPI_BARRIER(MPI_COMM_WORLD,IERR)

call P3DFFT_BTRAN_C2R(TMPFOU,TMPPHY,FFTFLAG)


!!====================================================================
!! 3. Statistic on divergence
!!====================================================================
if (NPROC>1) then
! Collective communications to compute max and mean values 

!!----------------------------------------------------------------------
!! 3.1. Max 
!!----------------------------------------------------------------------
!!- ( |d uf_i/d x_i| )_max
 call RMAXCPU(maxval(ABS(TMPPHY)),DIVMAX)

!  call MPI_ALLREDUCE(MAXVAL(ABS(TMPPHY)),DIVMAX, &
!                     1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERR)


!!----------------------------------------------------------------------
!! 3.2. Mean divergence
!!----------------------------------------------------------------------
!!- < d uf_i/d x_i >
 call RSUMCPU(SUM(TMPPHY),DIVMEAN)

!  call MPI_ALLREDUCE(SUM(TMPPHY),DIVMEAN, &
!                     1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)


!!----------------------------------------------------------------------
!! 3.3. Mean of squared divegrence
!!----------------------------------------------------------------------
!!- < (d uf_i/d x_i)^2 >
 call RSUMCPU(SUM(TMPPHY**2),DIVMEAN2)

!  call MPI_ALLREDUCE(SUM(TMPPHY**2),DIVMEAN2, &
!                     1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

else

!!----------------------------------------------------------------------
!! 3.1. Max
!!----------------------------------------------------------------------
!!- < d uf_i/d x_i >
 DIVMAX = MAXVAL(ABS(TMPPHY))

!!----------------------------------------------------------------------
!! 3.2. Mean divergence
!!----------------------------------------------------------------------
!!- < d uf_i/d x_i >
 DIVMEAN = SUM(TMPPHY)

!!----------------------------------------------------------------------
!! 3.3. Mean of squared divegrence
!!----------------------------------------------------------------------
!!- < (d uf_i/d x_i)^2 >
 DIVMEAN2 = SUM(TMPPHY**2)

end if

!!====================================================================
!! 2. Statistics normalization
!!====================================================================
DIVMEAN = DIVMEAN / NGLOB
DIVMEAN2 = DIVMEAN2 / NGLOB


end subroutine DIVERGENCE
