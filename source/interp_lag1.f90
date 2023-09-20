!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     
!                   Lagrangian Polynomial INTERPOLATION               
!                                                                     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     
! The velocity field "UAIN" located on the mesh "XAIN,YAIN,ZAIN" is   
! interpolated on the mesh "XINT,YINT,ZINT" and the result is "UINT". 
!                                                                     
!=====================================================================
subroutine INTERP_LAG1( XAIN,YAIN,ZAIN, &
                        UAIN,           &
                        NINT,           &
                        XINT,YINT,ZINT, &
                        UINT            )
!=====================================================================
!
!
!                  8--------7
!     k+1         /|       /|
!                5--------6 |
!                | |      | |
!                | 4------|-3
!                |/       |/
!      k         1--------2
!
!---------------------------------------------------------------------
!                             P. FEDE     --  I.M.F.T. --  31/03/2011
!---------------------------------------------------------------------

use dns_dim

implicit none

!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------

!=====================================================================
! Input Arrays
!=============
!- Mesh of data for interpolation
real(kind=8), dimension(ISTART(1):IEND(1)), intent(in) :: XAIN
real(kind=8), dimension(ISTART(2):IEND(2)), intent(in) :: YAIN
real(kind=8), dimension(ISTART(3):IEND(3)), intent(in) :: ZAIN

  
!- Data field for interpolation
real(kind=8),                                         &
   dimension(ISTART(1)         :IEND(1)               &
            ,ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL      &
            ,ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL ) ,  &
   intent(in) :: UAIN
  
  
!- Interpolated array size
integer, intent(in) :: NINT
  
!- Interpolated velocity field
real(kind=8), dimension(NINT), intent(out) :: UINT
   
!- Positions for interpolation 
real(kind=8), dimension(NINT), intent(in) :: XINT
real(kind=8), dimension(NINT), intent(in) :: YINT
real(kind=8), dimension(NINT), intent(in) :: ZINT
  
                               
!!- Normalized distance to the reference node
real(kind=8) :: ALFA, BETA, GAMA


!!- Lagrange's Polynom
real(kind=8) :: FX0, FX1
real(kind=8) :: FY0, FY1
real(kind=8) :: FZ0, FZ1

!!- Mesh step
real(kind=8) :: DXAINT
real(kind=8) :: DYAINT
real(kind=8) :: DZAINT

!!- Index
integer :: NP

integer :: I0,  J0,  K0  
integer :: IP1, JP1, KP1

intrinsic int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DXAINT = XAIN(ISTART(1)+1) - XAIN(ISTART(1))
DYAINT = YAIN(ISTART(2)+1) - YAIN(ISTART(2))
DZAINT = ZAIN(ISTART(3)+1) - ZAIN(ISTART(3))



do NP = 1,NINT


!!====================================================================
!! 1. Interpolation nodes location
!!====================================================================
 I0 = int(XINT(NP)/ DXAINT) + 1
 J0 = int(YINT(NP)/ DYAINT) + 1
 K0 = int(ZINT(NP)/ DZAINT) + 1


! x-direction --> periodic boundary condition
!----------------------------------------------
 IP1 = I0 + 1
 if(I0 == IEND(1)) IP1 = 1


! y-direction
!----------------
 JP1 = J0 + 1

! z-direction
!----------------
 KP1 = K0 + 1


!!====================================================================
!! 2. Distance to the first node location
!!====================================================================
 ALFA = (XINT(NP)-XAIN(I0)) / DXAINT
 BETA = (YINT(NP)-YAIN(J0)) / DYAINT
 GAMA = (ZINT(NP)-ZAIN(K0)) / DZAINT


!!====================================================================
!! 3. Lagrange polynomial functions
!!====================================================================

 FX0 = 1. - ALFA
 FX1 =      ALFA

 FY0 = 1. - BETA
 FY1 =      BETA

 FZ0 = 1. - GAMA
 FZ1 =      GAMA


!!====================================================================
!! 4.Interpolation
!!====================================================================
 UINT(NP) = 0.

 UINT(NP) = FZ0*(  FY0*(  FX0 * UAIN(I0, J0, K0 )     &
                        + FX1 * UAIN(IP1,J0, K0 ) )   &
                 + FY1*(  FX0 * UAIN(I0, JP1,K0 )     &
                        + FX1 * UAIN(IP1,JP1,K0 ) ) ) &
          + FZ1*(  FY0*(  FX0 * UAIN(I0, J0, KP1)     &
                        + FX1 * UAIN(IP1,J0, KP1) )   &
                 + FY1*(  FX0 * UAIN(I0, JP1,KP1)     &
                        + FX1 * UAIN(IP1,JP1,KP1) ) ) 

end do !- on NP = 1, NPMAX




10000 format (8(A,f6.2,4x))


end subroutine INTERP_LAG1
