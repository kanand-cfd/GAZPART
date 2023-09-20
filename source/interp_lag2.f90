!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      
!                   Lagrangian Polynomial INTERPOLATION                
!                                                                      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
! The velocity field "UAIN" located on the mesh "XAIN,YAIN,ZAIN" is     
! interpolated on the mesh "XINT,YINT,ZINT" and the result is "UINT".   
!                                                                       
!=====================================================================
subroutine INTERP_LAG2(XAIN,YAIN,ZAIN, &
                       UAIN,           &
                       NINT,           &
                       XINT,YINT,ZINT, &
                       UINT            )
!=====================================================================
!
!
!                25-------26-------27
!               /|       /        /|           
!              22-------23-------24|
!             /  |     /        /  |
!     k+2    19-------20-------21  |
!            |   |             |   |
!            |   16-------17---|---18
!            |  /|       /     |  /|           
!            | 13-------14-----|-15|
!            |/  |     /       |/  |
!     k+1    10-------11-------12  |
!            |   |             |   |
!            |   7--------8----|---9
!            |  /        /     |  /            
!            | 4--------5------|-6
!            |/        /       |/
!     k      1--------2--------3
!                                                                       !
!-----------------------------------------------------------------------!
!                             P. FEDE     --  I.M.F.T. --  31/03/2011   !
!-----------------------------------------------------------------------!

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
  
                               
 
integer :: IM1, JM1, KM1
integer :: I0,  J0,  K0  
integer :: IP1, JP1, KP1
                                

real(kind=8) :: ALFA, BETA, GAMA

real(kind=8) :: FXM1, FX0, FXP1 
real(kind=8) :: FYM1, FY0, FYP1 
real(kind=8) :: FZM1, FZ0, FZP1

real(kind=8) :: DXAINT
real(kind=8) :: DYAINT
real(kind=8) :: DZAINT

!!- Index
integer :: NP
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


!---------------------------------------------------------------------
! x-direction --> periodic boundary condition
!---------------------------------------------------------------------
 IP1 = I0 + 1
 if(I0 == IEND(1)) IP1 = 1

 IM1 = I0 - 1
 if(I0 == ISTART(1)) IM1 = IEND(1)

!---------------------------------------------------------------------
! y-direction
!---------------------------------------------------------------------
 JP1 = J0 + 1
 JM1 = J0 - 1


!---------------------------------------------------------------------
! z-direction
!---------------------------------------------------------------------
 KP1 = K0 + 1
 KM1 = K0 - 1



!!====================================================================
!! 2. Distance to the first node location
!!====================================================================
 ALFA = (XINT(NP)-XAIN(I0)) / DXAINT
 BETA = (YINT(NP)-YAIN(J0)) / DYAINT
 GAMA = (ZINT(NP)-ZAIN(K0)) / DZAINT


!!====================================================================
!! 3. Lagrange polynomial functions
!!====================================================================
 FXM1 =  0.5*ALFA*(ALFA - 1.)
 FX0  = -         (ALFA - 1.)*(ALFA + 1.)
 FXP1 =  0.5*ALFA            *(ALFA + 1.)

 FYM1 =  0.5*BETA*(BETA - 1.)
 FY0  = -         (BETA - 1.)*(BETA + 1.)
 FYP1 =  0.5*BETA            *(BETA + 1.)

 FZM1 =  0.5*GAMA*(GAMA - 1.)
 FZ0  = -         (GAMA - 1.)*(GAMA + 1.)
 FZP1 =  0.5*GAMA            *(GAMA + 1.)


!!====================================================================
!! 4.Interpolation
!!====================================================================

 UINT(NP) = 0.


 UINT(NP) =                                        &
   FZM1*(  FYM1*(   FXM1 * UAIN(IM1,JM1,KM1)       &
                  + FX0  * UAIN(I0 ,JM1,KM1)       &
                  + FXP1 * UAIN(IP1,JM1,KM1)  )    &
         + FY0 *(   FXM1 * UAIN(IM1,J0 ,KM1)       &
		  + FX0  * UAIN(I0 ,J0 ,KM1)       &
              	  + FXP1 * UAIN(IP1,J0 ,KM1)  )    &
         + FYP1*(   FXM1 * UAIN(IM1,JP1,KM1)       &
		  + FX0  * UAIN(I0 ,JP1,KM1)       &
		  + FXP1 * UAIN(IP1,JP1,KM1)  )  ) &
!
 + FZ0 *(  FYM1*(   FXM1 * UAIN(IM1,JM1,K0)        &
                  + FX0  * UAIN(I0 ,JM1,K0)        &
                  + FXP1 * UAIN(IP1,JM1,K0)   )    &
         + FY0 *(   FXM1 * UAIN(IM1,J0 ,K0)        &
		  + FX0  * UAIN(I0 ,J0 ,K0)        &
              	  + FXP1 * UAIN(IP1,J0 ,K0)   )    &
         + FYP1*(   FXM1 * UAIN(IM1,JP1,K0)        &
		  + FX0  * UAIN(I0 ,JP1,K0)        &
		  + FXP1 * UAIN(IP1,JP1,K0)   )  ) &
!
 + FZP1*(  FYM1*(   FXM1 * UAIN(IM1,JM1,KP1)       &
                  + FX0  * UAIN(I0 ,JM1,KP1)       &
                  + FXP1 * UAIN(IP1,JM1,KP1)    )  &
         + FY0 *(   FXM1 * UAIN(IM1,J0 ,KP1)       &
		  + FX0  * UAIN(I0 ,J0 ,KP1)       &
              	  + FXP1 * UAIN(IP1,J0 ,KP1)    )  &
         + FYP1*(   FXM1 * UAIN(IM1,JP1,KP1)       &
		  + FX0  * UAIN(I0 ,JP1,KP1)       &
		  + FXP1 * UAIN(IP1,JP1,KP1)    )  ) 

end do !- on NP = 1, NINT



10000 format (8(A,f6.2,4x))


end subroutine INTERP_LAG2
