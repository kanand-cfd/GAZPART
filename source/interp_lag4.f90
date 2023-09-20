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
subroutine INTERP_LAG4(XAIN,YAIN,ZAIN, &
                       UAIN,           &
                       NINT,           &
                       XINT,YINT,ZINT, &
                       UINT            )
!=====================================================================
!
!  No sketch, too boring to do ...
!                                                                       
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
  
                               
 
integer :: IM2, JM2, KM2
integer :: IM1, JM1, KM1
integer :: I0,  J0,  K0  
integer :: IP1, JP1, KP1
integer :: IP2, JP2, KP2
      

real(kind=8) :: ALFA, BETA, GAMA


real(kind=8) :: FXM2, FXM1, FX0, FXP1, FXP2 
real(kind=8) :: FYM2, FYM1, FY0, FYP1, FYP2 
real(kind=8) :: FZM2, FZM1, FZ0, FZP1, FZP2

real(kind=8) :: DXAINT
real(kind=8) :: DYAINT
real(kind=8) :: DZAINT


!!--------------------------------------------------------------------
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

 IP2 = IP1 + 1
 if(IP1 == IEND(1)) IP2 = 1

 IM1 = I0 - 1
 if(I0 == ISTART(1)) IM1 = IEND(1)

 IM2 = IM1 - 1
 if(IM1 == ISTART(1)) IM2 = IEND(1)


!---------------------------------------------------------------------
! y-direction
!---------------------------------------------------------------------
 JM2 = J0 - 2
 JM1 = J0 - 1
 JP1 = J0 + 1
 JP2 = J0 + 2


!---------------------------------------------------------------------
! z-direction
!---------------------------------------------------------------------
 KM2 = K0 - 2
 KM1 = K0 - 1
 KP1 = K0 + 1
 KP2 = K0 + 2



!!====================================================================
!! 2. Distance to the first node location
!!====================================================================
 ALFA = (XINT(NP)-XAIN(I0)) / DXAINT
 BETA = (YINT(NP)-YAIN(J0)) / DYAINT
 GAMA = (ZINT(NP)-ZAIN(K0)) / DZAINT


!!====================================================================
!! 3. Lagrange polynomial functions
!!====================================================================
 FXM2 =  1./24.*(ALFA - 2.)*(ALFA - 1.)*ALFA*(ALFA + 1.)
 FXM1 = -1./6. *(ALFA - 2.)*(ALFA - 1.)*ALFA            *(ALFA + 2.)
 FX0  =  1./4. *(ALFA - 2.)*(ALFA - 1.)     *(ALFA + 1.)*(ALFA + 2.)
 FXP1 = -1./6. *(ALFA - 2.)            *ALFA*(ALFA + 1.)*(ALFA + 2.)
 FXP2 =  1./24.            *(ALFA - 1.)*ALFA*(ALFA + 1.)*(ALFA + 2.)


 FYM2 =  1./24.*(BETA - 2.)*(BETA - 1.)*BETA*(BETA + 1.)
 FYM1 = -1./6. *(BETA - 2.)*(BETA - 1.)*BETA            *(BETA + 2.)
 FY0  =  1./4. *(BETA - 2.)*(BETA - 1.)     *(BETA + 1.)*(BETA + 2.)
 FYP1 = -1./6. *(BETA - 2.)            *BETA*(BETA + 1.)*(BETA + 2.)
 FYP2 =  1./24.            *(BETA - 1.)*BETA*(BETA + 1.)*(BETA + 2.)


 FZM2 =  1./24.*(GAMA - 2.)*(GAMA - 1.)*GAMA*(GAMA + 1.)
 FZM1 = -1./6. *(GAMA - 2.)*(GAMA - 1.)*GAMA            *(GAMA + 2.)
 FZ0  =  1./4. *(GAMA - 2.)*(GAMA - 1.)     *(GAMA + 1.)*(GAMA + 2.)
 FZP1 = -1./6. *(GAMA - 2.)            *GAMA*(GAMA + 1.)*(GAMA + 2.)
 FZP2 =  1./24.            *(GAMA - 1.)*GAMA*(GAMA + 1.)*(GAMA + 2.)


!!====================================================================
!! 4.Interpolation
!!====================================================================

 UINT(NP) = 0.

 UINT(NP) =                                        &
   FZM2*(  FYM2*(   FXM2 * UAIN(IM2,JM2,KM2)       &
                  + FXM1 * UAIN(IM1,JM2,KM2)       &
                  + FX0  * UAIN(I0 ,JM2,KM2)       &
                  + FXP1 * UAIN(IP1,JM2,KM2)       &
                  + FXP2 * UAIN(IP2,JM2,KM2)  )    &
         + FYM1*(   FXM2 * UAIN(IM2,JM1,KM2)       &
                  + FXM1 * UAIN(IM1,JM1,KM2)       &
                  + FX0  * UAIN(I0 ,JM1,KM2)       &
                  + FXP1 * UAIN(IP1,JM1,KM2)       &
                  + FXP2 * UAIN(IP2,JM1,KM2)  )    &
         + FY0 *(   FXM2 * UAIN(IM2,J0 ,KM2)       &
                  + FXM1 * UAIN(IM1,J0 ,KM2)       &
                  + FX0  * UAIN(I0 ,J0 ,KM2)       &
                  + FXP1 * UAIN(IP1,J0 ,KM2)       &
                  + FXP2 * UAIN(IP2,J0 ,KM2)  )    &
         + FYP1*(   FXM2 * UAIN(IM2,JP1,KM2)       &
                  + FXM1 * UAIN(IM1,JP1,KM2)       &
                  + FX0  * UAIN(I0 ,JP1,KM2)       &
                  + FXP1 * UAIN(IP1,JP1,KM2)       &
                  + FXP2 * UAIN(IP2,JP1,KM2)  )    &
         + FYP2*(   FXM2 * UAIN(IM2,JP1,KM2)       &
                  + FXM1 * UAIN(IM1,JP1,KM2)       &
                  + FX0  * UAIN(I0 ,JP1,KM2)       &
                  + FXP1 * UAIN(IP1,JP1,KM2)       &
                  + FXP2 * UAIN(IP2,JP1,KM2)  ) )  &
!
 + FZM1*(  FYM2*(   FXM2 * UAIN(IM2,JM2,KM1)       &
                  + FXM1 * UAIN(IM1,JM2,KM1)       &
                  + FX0  * UAIN(I0 ,JM2,KM1)       &
                  + FXP1 * UAIN(IP1,JM2,KM1)       &
                  + FXP2 * UAIN(IP2,JM2,KM1)  )    &
         + FYM1*(   FXM2 * UAIN(IM2,JM1,KM1)       &
                  + FXM1 * UAIN(IM1,JM1,KM1)       &
                  + FX0  * UAIN(I0 ,JM1,KM1)       &
                  + FXP1 * UAIN(IP1,JM1,KM1)       &
                  + FXP2 * UAIN(IP2,JM1,KM1)  )    &
         + FY0 *(   FXM2 * UAIN(IM2,J0 ,KM1)       &
                  + FXM1 * UAIN(IM1,J0 ,KM1)       &
                  + FX0  * UAIN(I0 ,J0 ,KM1)       &
                  + FXP1 * UAIN(IP1,J0 ,KM1)       &
                  + FXP2 * UAIN(IP2,J0 ,KM1)  )    &
         + FYP1*(   FXM2 * UAIN(IM2,JP1,KM1)       &
                  + FXM1 * UAIN(IM1,JP1,KM1)       &
                  + FX0  * UAIN(I0 ,JP1,KM1)       &
                  + FXP1 * UAIN(IP1,JP1,KM1)       &
                  + FXP2 * UAIN(IP2,JP1,KM1)  )    &
         + FYP2*(   FXM2 * UAIN(IM2,JP1,KM1)       &
                  + FXM1 * UAIN(IM1,JP1,KM1)       &
                  + FX0  * UAIN(I0 ,JP1,KM1)       &
                  + FXP1 * UAIN(IP1,JP1,KM1)       &
                  + FXP2 * UAIN(IP2,JP1,KM1)  ) )  &
!
 + FZ0 *(  FYM2*(   FXM2 * UAIN(IM2,JM2,K0 )       &
                  + FXM1 * UAIN(IM1,JM2,K0 )       &
                  + FX0  * UAIN(I0 ,JM2,K0 )       &
                  + FXP1 * UAIN(IP1,JM2,K0 )       &
                  + FXP2 * UAIN(IP2,JM2,K0 )  )    &
         + FYM1*(   FXM2 * UAIN(IM2,JM1,K0 )       &
                  + FXM1 * UAIN(IM1,JM1,K0 )       &
                  + FX0  * UAIN(I0 ,JM1,K0 )       &
                  + FXP1 * UAIN(IP1,JM1,K0 )       &
                  + FXP2 * UAIN(IP2,JM1,K0 )  )    &
         + FY0 *(   FXM2 * UAIN(IM2,J0 ,K0 )       &
                  + FXM1 * UAIN(IM1,J0 ,K0 )       &
                  + FX0  * UAIN(I0 ,J0 ,K0 )       &
                  + FXP1 * UAIN(IP1,J0 ,K0 )       &
                  + FXP2 * UAIN(IP2,J0 ,K0 )  )    &
         + FYP1*(   FXM2 * UAIN(IM2,JP1,K0 )       &
                  + FXM1 * UAIN(IM1,JP1,K0 )       &
                  + FX0  * UAIN(I0 ,JP1,K0 )       &
                  + FXP1 * UAIN(IP1,JP1,K0 )       &
                  + FXP2 * UAIN(IP2,JP1,K0 )  )    &
         + FYP2*(   FXM2 * UAIN(IM2,JP1,K0 )       &
                  + FXM1 * UAIN(IM1,JP1,K0 )       &
                  + FX0  * UAIN(I0 ,JP1,K0 )       &
                  + FXP1 * UAIN(IP1,JP1,K0 )       &
                  + FXP2 * UAIN(IP2,JP1,K0 )  ) )  &
!
 + FZP1*(  FYM2*(   FXM2 * UAIN(IM2,JM2,KP1)       &
                  + FXM1 * UAIN(IM1,JM2,KP1)       &
                  + FX0  * UAIN(I0 ,JM2,KP1)       &
                  + FXP1 * UAIN(IP1,JM2,KP1)       &
                  + FXP2 * UAIN(IP2,JM2,KP1)  )    &
         + FYM1*(   FXM2 * UAIN(IM2,JM1,KP1)       &
                  + FXM1 * UAIN(IM1,JM1,KP1)       &
                  + FX0  * UAIN(I0 ,JM1,KP1)       &
                  + FXP1 * UAIN(IP1,JM1,KP1)       &
                  + FXP2 * UAIN(IP2,JM1,KP1)  )    &
         + FY0 *(   FXM2 * UAIN(IM2,J0 ,KP1)       &
                  + FXM1 * UAIN(IM1,J0 ,KP1)       &
                  + FX0  * UAIN(I0 ,J0 ,KP1)       &
                  + FXP1 * UAIN(IP1,J0 ,KP1)       &
                  + FXP2 * UAIN(IP2,J0 ,KP1)  )    &
         + FYP1*(   FXM2 * UAIN(IM2,JP1,KP1)       &
                  + FXM1 * UAIN(IM1,JP1,KP1)       &
                  + FX0  * UAIN(I0 ,JP1,KP1)       &
                  + FXP1 * UAIN(IP1,JP1,KP1)       &
                  + FXP2 * UAIN(IP2,JP1,KP1)  )    &
         + FYP2*(   FXM2 * UAIN(IM2,JP1,KP1)       &
                  + FXM1 * UAIN(IM1,JP1,KP1)       &
                  + FX0  * UAIN(I0 ,JP1,KP1)       &
                  + FXP1 * UAIN(IP1,JP1,KP1)       &
                  + FXP2 * UAIN(IP2,JP1,KP1)  ) )  &
!
 + FZP2*(  FYM2*(   FXM2 * UAIN(IM2,JM2,KP2)       &
                  + FXM1 * UAIN(IM1,JM2,KP2)       &
                  + FX0  * UAIN(I0 ,JM2,KP2)       &
                  + FXP1 * UAIN(IP1,JM2,KP2)       &
                  + FXP2 * UAIN(IP2,JM2,KP2)  )    &
         + FYM1*(   FXM2 * UAIN(IM2,JM1,KP2)       &
                  + FXM1 * UAIN(IM1,JM1,KP2)       &
                  + FX0  * UAIN(I0 ,JM1,KP2)       &
                  + FXP1 * UAIN(IP1,JM1,KP2)       &
                  + FXP2 * UAIN(IP2,JM1,KP2)  )    &
         + FY0 *(   FXM2 * UAIN(IM2,J0 ,KP2)       &
                  + FXM1 * UAIN(IM1,J0 ,KP2)       &
                  + FX0  * UAIN(I0 ,J0 ,KP2)       &
                  + FXP1 * UAIN(IP1,J0 ,KP2)       &
                  + FXP2 * UAIN(IP2,J0 ,KP2)  )    &
         + FYP1*(   FXM2 * UAIN(IM2,JP1,KP2)       &
                  + FXM1 * UAIN(IM1,JP1,KP2)       &
                  + FX0  * UAIN(I0 ,JP1,KP2)       &
                  + FXP1 * UAIN(IP1,JP1,KP2)       &
                  + FXP2 * UAIN(IP2,JP1,KP2)  )    &
         + FYP2*(   FXM2 * UAIN(IM2,JP1,KP2)       &
                  + FXM1 * UAIN(IM1,JP1,KP2)       &
                  + FX0  * UAIN(I0 ,JP1,KP2)       &
                  + FXP1 * UAIN(IP1,JP1,KP2)       &
                  + FXP2 * UAIN(IP2,JP1,KP2)  )  )  

end do !- on NP = 1, NINT



10000 format (8(A,f6.2,4x))

end subroutine INTERP_LAG4
