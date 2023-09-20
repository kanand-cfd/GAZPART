!!====================================================================
!!
!! Sum over all cpu for real
!!
!!====================================================================
subroutine RSUMCPU(VARIN,VAROUT)

use dns_dim            !- Dimension

implicit none

!!--------------------------------------------------------------------
real(kind=8), intent(in) :: VARIN
real(kind=8), intent(out) :: VAROUT
!!--------------------------------------------------------------------

VAROUT = VARIN

if (NPROC>1) call MPI_ALLREDUCE(VARIN,VAROUT,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

end subroutine RSUMCPU




!!====================================================================
!!
!! Sum over all cpu for integer
!!
!!====================================================================
subroutine ISUMCPU(VARIN,VAROUT)

use dns_dim            !- Dimension

implicit none

!!--------------------------------------------------------------------
integer, intent(in) :: VARIN
integer, intent(out) :: VAROUT
!!--------------------------------------------------------------------


VAROUT = VARIN

if (NPROC>1) call MPI_ALLREDUCE(VARIN,VAROUT,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)


end subroutine ISUMCPU




!!====================================================================
!!
!! Maximum over all cpu for real
!!
!!====================================================================
subroutine RMAXCPU(VARIN,VAROUT)

use dns_dim            !- Dimension

implicit none

!!--------------------------------------------------------------------
real(kind=8), intent(in) :: VARIN
real(kind=8), intent(out) :: VAROUT
!!--------------------------------------------------------------------

VAROUT = VARIN

if (NPROC>1) call MPI_ALLREDUCE(VARIN,VAROUT,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERR)


end subroutine RMAXCPU



!!====================================================================
!!
!! Maximum over all cpu for integer
!!
!!====================================================================
subroutine IMAXCPU(VARIN,VAROUT)

use dns_dim            !- Dimension

implicit none

!!--------------------------------------------------------------------
integer, intent(in) :: VARIN
integer, intent(out) :: VAROUT
!!--------------------------------------------------------------------

VAROUT = VARIN

if (NPROC>1) call MPI_ALLREDUCE(VARIN,VAROUT,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERR)

end subroutine IMAXCPU


!!====================================================================
!!
!! Minumum over all cpu for real
!!
!!====================================================================
subroutine RMINCPU(VARIN,VAROUT)

use dns_dim            !- Dimension

implicit none

!!--------------------------------------------------------------------
real(kind=8), intent(in) :: VARIN
real(kind=8), intent(out) :: VAROUT
!!--------------------------------------------------------------------

VAROUT = VARIN

if (NPROC>1) call MPI_ALLREDUCE(VARIN,VAROUT,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERR)


end subroutine RMINCPU



!!====================================================================
!!
!! Maximum over all cpu for integer
!!
!!====================================================================
subroutine IMINCPU(VARIN,VAROUT)

use dns_dim            !- Dimension

implicit none

!!--------------------------------------------------------------------
integer, intent(in) :: VARIN
integer, intent(out) :: VAROUT
!!--------------------------------------------------------------------

VAROUT = VARIN

if (NPROC>1) call MPI_ALLREDUCE(VARIN,VAROUT,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,IERR)

end subroutine IMINCPU
