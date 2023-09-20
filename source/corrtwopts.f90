!!====================================================================
!!
!!     Spatial correlation function computation
!!
!!====================================================================

subroutine CORRTWOPTS(NI,NJ,NK,U,UM,CORR)

!!====================================================================
!!
!!
!!====================================================================

implicit none

!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Global arrays
!---------------
!- Size
integer,                              intent(in) :: NI, NJ, NK

!- Fluid velocity
real(kind=8), dimension(NI,NJ,NK), intent(in) :: U

!- Mean fluid velocity
real(kind=8),                          intent(in) :: UM


!- Spatial correlation
real(kind=8), dimension(NI/2+1),   intent(out) :: CORR

!---------------------------------------------------------------------
!- Local arrays
!--------------
!- Correlation
real(kind=8), dimension(NI/2+1) ::  ULRS
real(kind=8), dimension(NI/2+1) :: UULRS

real(kind=8), dimension(NI,NI/2+1) ::  ULR
real(kind=8), dimension(NI,NI/2+1) :: UULR


!- Index
integer :: I, J, K, IJ, IK, JK, LR, M, N
!---------------------------------------------------------------------


!!====================================================================
!! 1. z-direction
!!====================================================================
do I = 1,NI

 do LR = 1, NI/2+1

  ULR(I,LR) = 0.
  UULR(I,LR) = 0.

  M = I + LR - 1
  if (M > NI) M = I + LR - 1 - NI

  N = I - LR + 1
  if (N < 1 ) N = I - LR + 1 + NI

  do JK = 1,NJ*NK
   K = INT((JK-1)/NJ) + 1
   J = JK - NJ*(K-1)

    ULR(I,LR) =  ULR(I,LR) +          (U(M,J,K)+U(N,J,K))/2./NJ/NK
   UULR(I,LR) = UULR(I,LR) + U(I,J,K)*(U(M,J,K)+U(N,J,K))/2./NJ/NK

  end do
 end do
end do

ULRS(:) = 0.
UULRS(:) = 0.

do LR = 1,NI/2+1
 do I = 1,NI
   ULRS(LR) =  ULRS(LR) +  ULR(I,LR)/NI
  UULRS(LR) = UULRS(LR) + UULR(I,LR)/NI
 end do
end do

do LR = 1, NI/2+1
 CORR(LR) = UULRS(LR) - UM*ULRS(LR)
end do


end subroutine CORRTWOPTS
