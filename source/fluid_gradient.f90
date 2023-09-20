!!====================================================================
!!
!!
!!====================================================================

subroutine FLUID_GRADIENT

!!====================================================================
!!
!!
!!====================================================================

use DNS_DIM
use FLUID_VARIABLE     !- Fluid velocity
use GEOMETRIC_VARIABLE
use PARAM_PHYS
use STATISTICS
use WORK_ARRAYS

use P3DFFT

implicit none


!------------------------------------------------------------------
! ARRAYS STATEMENT
!------------------------------------------------------------------

!- Index
integer :: I, J, K, IJK
!---------------------------------------------------------------------

DUIDXJ(:,:,:) = 0.

!!====================================================================
!! 1. (du/dx)
!!====================================================================
do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
   TMPFOU(I,J,K) = ICMPL*KX(I)*UFOU(I,J,K)
  end do
 end do
end do

!- Synchronize all the process
call MPI_BARRIER(MPI_COMM_WORLD,IERR)

call P3DFFT_BTRAN_C2R(TMPFOU,TMPPHY,FFTFLAG)       

do K = ISTART(3), IEND(3)
 do J = ISTART(2), IEND(2)
  do I = ISTART(1), IEND(1)
    
   !- < (dui/dxi)>
   DUIDXJ(1,1,1)= DUIDXJ(1,1,1) + TMPPHY(I,J,K)

   !- < (dui/dxi)^2>
   DUIDXJ(1,1,2)= DUIDXJ(1,1,2) + TMPPHY(I,J,K)**2

   !- < (dui/dxi)^3>
   DUIDXJ(1,1,3)= DUIDXJ(1,1,3) + TMPPHY(I,J,K)**3

   !- < (dui/dxi)^4>
   DUIDXJ(1,1,4)= DUIDXJ(1,1,4) + TMPPHY(I,J,K)**4

  end do
 end do
end do



!!====================================================================
!! 2. (dv/dy)
!!====================================================================
do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
   TMPFOU(I,J,K) = ICMPL*KY(J)*VFOU(I,J,K)
  end do
 end do
end do

!- Synchronize all the process
call MPI_BARRIER(MPI_COMM_WORLD,IERR)

call P3DFFT_BTRAN_C2R(TMPFOU,TMPPHY,FFTFLAG)

do K = ISTART(3), IEND(3)
 do J = ISTART(2), IEND(2)
  do I = ISTART(1), IEND(1)
    
   !- < (dui/dxi)>
   DUIDXJ(2,2,1)= DUIDXJ(2,2,1) + TMPPHY(I,J,K)

   !- < (dui/dxi)^2>
   DUIDXJ(2,2,2)= DUIDXJ(2,2,2) + TMPPHY(I,J,K)**2

   !- < (dui/dxi)^3>
   DUIDXJ(2,2,3)= DUIDXJ(2,2,3) + TMPPHY(I,J,K)**3

   !- < (dui/dxi)^4>
   DUIDXJ(2,2,4)= DUIDXJ(2,2,4) + TMPPHY(I,J,K)**4

  end do
 end do
end do


!!====================================================================
!! 3. (dw/dz)
!!====================================================================
do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
   TMPFOU(I,J,K) = ICMPL*KZ(K)*WFOU(I,J,K)
  end do
 end do
end do

!- Synchronize all the process
call MPI_BARRIER(MPI_COMM_WORLD,IERR)

call P3DFFT_BTRAN_C2R(TMPFOU,TMPPHY,FFTFLAG)

do K = ISTART(3), IEND(3)
 do J = ISTART(2), IEND(2)
  do I = ISTART(1), IEND(1)
    
   !- < (dui/dxi)>
   DUIDXJ(3,3,1)= DUIDXJ(3,3,1) + TMPPHY(I,J,K)

   !- < (dui/dxi)^2>
   DUIDXJ(3,3,2)= DUIDXJ(3,3,2) + TMPPHY(I,J,K)**2

   !- < (dui/dxi)^3>
   DUIDXJ(3,3,3)= DUIDXJ(3,3,3) + TMPPHY(I,J,K)**3

   !- < (dui/dxi)^4>
   DUIDXJ(3,3,4)= DUIDXJ(3,3,4) + TMPPHY(I,J,K)**4

  end do
 end do
end do



!!====================================================================
!! 4. (du/dy)
!!====================================================================
do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
   TMPFOU(I,J,K) = ICMPL*KY(J)*UFOU(I,J,K)
  end do
 end do
end do

!- Synchronize all the process
call MPI_BARRIER(MPI_COMM_WORLD,IERR)

call P3DFFT_BTRAN_C2R(TMPFOU,TMPPHY,FFTFLAG)

do K = ISTART(3), IEND(3)
 do J = ISTART(2), IEND(2)
  do I = ISTART(1), IEND(1)
    
   !- < (dui/dxi)>
   DUIDXJ(1,2,1)= DUIDXJ(1,2,1) + TMPPHY(I,J,K)

   !- < (dui/dxi)^2>
   DUIDXJ(1,2,2)= DUIDXJ(1,2,2) + TMPPHY(I,J,K)**2

   !- < (dui/dxi)^3>
   DUIDXJ(1,2,3)= DUIDXJ(1,2,3) + TMPPHY(I,J,K)**3

   !- < (dui/dxi)^4>
   DUIDXJ(1,2,4)= DUIDXJ(1,2,4) + TMPPHY(I,J,K)**4

  end do
 end do
end do



!!====================================================================
!! 5. (du/dz)
!!====================================================================
do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
   TMPFOU(I,J,K) = ICMPL*KZ(K)*UFOU(I,J,K)
  end do
 end do
end do

!- Synchronize all the process
call MPI_BARRIER(MPI_COMM_WORLD,IERR)

call P3DFFT_BTRAN_C2R(TMPFOU,TMPPHY,FFTFLAG)

do K = ISTART(3), IEND(3)
 do J = ISTART(2), IEND(2)
  do I = ISTART(1), IEND(1)
    
   !- < (dui/dxi)>
   DUIDXJ(1,3,1)= DUIDXJ(1,3,1) + TMPPHY(I,J,K)

   !- < (dui/dxi)^2>
   DUIDXJ(1,3,2)= DUIDXJ(1,3,2) + TMPPHY(I,J,K)**2

   !- < (dui/dxi)^3>
   DUIDXJ(1,3,3)= DUIDXJ(1,3,3) + TMPPHY(I,J,K)**3

   !- < (dui/dxi)^4>
   DUIDXJ(1,3,4)= DUIDXJ(1,3,4) + TMPPHY(I,J,K)**4

  end do
 end do
end do



!!====================================================================
!! 6. (dv/dx)
!!====================================================================
do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
   TMPFOU(I,J,K) = ICMPL*KX(I)*VFOU(I,J,K)
  end do
 end do
end do

!- Synchronize all the process
call MPI_BARRIER(MPI_COMM_WORLD,IERR)

call P3DFFT_BTRAN_C2R(TMPFOU,TMPPHY,FFTFLAG)

do K = ISTART(3), IEND(3)
 do J = ISTART(2), IEND(2)
  do I = ISTART(1), IEND(1)
    
   !- < (dui/dxi)>
   DUIDXJ(2,1,1)= DUIDXJ(2,1,1) + TMPPHY(I,J,K)

   !- < (dui/dxi)^2>
   DUIDXJ(2,1,2)= DUIDXJ(2,1,2) + TMPPHY(I,J,K)**2

   !- < (dui/dxi)^3>
   DUIDXJ(2,1,3)= DUIDXJ(2,1,3) + TMPPHY(I,J,K)**3

   !- < (dui/dxi)^4>
   DUIDXJ(2,1,4)= DUIDXJ(2,1,4) + TMPPHY(I,J,K)**4

  end do
 end do
end do



!!====================================================================
!! 7. (dv/dz)
!!====================================================================
do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
   TMPFOU(I,J,K) = ICMPL*KZ(K)*VFOU(I,J,K)
  end do
 end do
end do

!- Synchronize all the process
call MPI_BARRIER(MPI_COMM_WORLD,IERR)

call P3DFFT_BTRAN_C2R(TMPFOU,TMPPHY,FFTFLAG)

do K = ISTART(3), IEND(3)
 do J = ISTART(2), IEND(2)
  do I = ISTART(1), IEND(1)
    
   !- < (dui/dxi)>
   DUIDXJ(2,3,1)= DUIDXJ(2,3,1) + TMPPHY(I,J,K)

   !- < (dui/dxi)^2>
   DUIDXJ(2,3,2)= DUIDXJ(2,3,2) + TMPPHY(I,J,K)**2

   !- < (dui/dxi)^3>
   DUIDXJ(2,3,3)= DUIDXJ(2,3,3) + TMPPHY(I,J,K)**3

   !- < (dui/dxi)^4>
   DUIDXJ(2,3,4)= DUIDXJ(2,3,4) + TMPPHY(I,J,K)**4

  end do
 end do
end do



!!====================================================================
!! 8. (dw/dx)
!!====================================================================
do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
   TMPFOU(I,J,K) = ICMPL*KX(I)*WFOU(I,J,K)
  end do
 end do
end do

!- Synchronize all the process
call MPI_BARRIER(MPI_COMM_WORLD,IERR)

call P3DFFT_BTRAN_C2R(TMPFOU,TMPPHY,FFTFLAG)

do K = ISTART(3), IEND(3)
 do J = ISTART(2), IEND(2)
  do I = ISTART(1), IEND(1)
    
   !- < (dui/dxi)>
   DUIDXJ(3,1,1)= DUIDXJ(3,1,1) + TMPPHY(I,J,K)

   !- < (dui/dxi)^2>
   DUIDXJ(3,1,2)= DUIDXJ(3,1,2) + TMPPHY(I,J,K)**2

   !- < (dui/dxi)^3>
   DUIDXJ(3,1,3)= DUIDXJ(3,1,3) + TMPPHY(I,J,K)**3

   !- < (dui/dxi)^4>
   DUIDXJ(3,1,4)= DUIDXJ(3,1,4) + TMPPHY(I,J,K)**4

  end do
 end do
end do



!!====================================================================
!! 9. (dw/dy)
!!====================================================================
do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
   TMPFOU(I,J,K) = ICMPL*KY(J)*WFOU(I,J,K)
  end do
 end do
end do

!- Synchronize all the process
call MPI_BARRIER(MPI_COMM_WORLD,IERR)

call P3DFFT_BTRAN_C2R(TMPFOU,TMPPHY,FFTFLAG)

do K = ISTART(3), IEND(3)
 do J = ISTART(2), IEND(2)
  do I = ISTART(1), IEND(1)
    
   !- < (dui/dxi)>
   DUIDXJ(3,2,1)= DUIDXJ(3,2,1) + TMPPHY(I,J,K)

   !- < (dui/dxi)^2>
   DUIDXJ(3,2,2)= DUIDXJ(3,2,2) + TMPPHY(I,J,K)**2

   !- < (dui/dxi)^3>
   DUIDXJ(3,2,3)= DUIDXJ(3,2,3) + TMPPHY(I,J,K)**3

   !- < (dui/dxi)^4>
   DUIDXJ(3,2,4)= DUIDXJ(3,2,4) + TMPPHY(I,J,K)**4

  end do
 end do
end do



!!====================================================================
!! 10. Normalization
!!====================================================================

 DUIDXJ = DUIDXJ / NGLOB


end subroutine FLUID_GRADIENT
