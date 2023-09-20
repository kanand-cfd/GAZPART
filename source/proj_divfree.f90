!!====================================================================
!!
!!
!!====================================================================

subroutine PROJ_DIVFREE

!!====================================================================
!!
!!
!!====================================================================

use DNS_DIM	       !- Dimension
use PARAM_PHYS         !- Physical & numerical parameters
use FLUID_VARIABLE     !- Fluid velocity
use GEOMETRIC_VARIABLe !- 
use STATISTICS         !- Statistics
use CHECK_CPU


use P3DFFT

implicit none


!------------------------------------------------------------------
! ARRAYS STATEMENT
!------------------------------------------------------------------
!- Divergence in Fourier space
double complex :: DIVC

!- Squared wavenumber
real(kind=8) :: KAPPA2
 
!- Time control variable
real(kind=8) :: TIME_START, TIME_END

!- index
integer :: I, J, K, IJK
!-----------------------------------------------------------------


!!- CPU check
if(MYID == 0) then
 TIME_START = MPI_WTIME()
end if




do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)

   KAPPA2 = KX(I)**2 + KY(J)**2 +  KZ(K)**2

   if(KAPPA2 <=ZERO) KAPPA2 = INFINITY

   !- Compute the divergence
   DIVC = KX(I)*UFOU(I,J,K) + KY(J)*VFOU(I,J,K) + KZ(K)*WFOU(I,J,K)

   UFOU(I,J,K) = UFOU(I,J,K) - KX(I)/KAPPA2*DIVC
   VFOU(I,J,K) = VFOU(I,J,K) - KY(J)/KAPPA2*DIVC
   WFOU(I,J,K) = WFOU(I,J,K) - KZ(K)/KAPPA2*DIVC

  end do
 end do
end do





!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()
 CPU_FLUID(5) = CPU_FLUID(5)  + TIME_END - TIME_START
end if


end subroutine PROJ_DIVFREE
