!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      
!                   Lagrangian Polynomial INTERPOLATION                
!                                                                      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       
! The velocity field "UAIN" located on the mesh "XAIN(NX),YAIN(NY),ZAIN(NZ)" is     
! interpolated on the mesh "XINT(NINT),YINT(NINT),ZINT(NINT)" and the result is "UINT".   
!                                                                       
!=====================================================================
subroutine FROZEN_INTERP_LAG3(XAIN,YAIN,ZAIN, &
                              NX,NY,NZ,       &
                              UAIN,           &
                              NINT,           &
                              XINT,YINT,ZINT, &
                              UINT            )
!=====================================================================
!
!                  61-------62-------63-------64
!                 /|       /        /        /|
!                57-------58-------59-------60|
!               /  |     /        /        /  |  
!              53-------54-------55-------56  |
!             /    |   /        /        /    |
!     k+3    49-------50-------51-------52    |
!            |     |                    |     |
!            |     45-------46-------47-|-----48
!            |    /|       /        /   |    /|
!            |   41-------42-------43---|---44|
!            |  /  |     /        /     |  /  |  
!            | 37-------38-------39-----|-40  |
!            |/    |   /        /       |/    |
!     k+2    33-------34-------35-------36    |
!            |     |                    |     |
!            |     29-------30-------31-|-----32
!            |    /|       /        /   |    /|
!            |   25-------26-------27---|---28|
!            |  /  |     /        /     |  /  |  
!            | 21-------22-------23-----|-24  |
!            |/    |   /        /       |/    |
!     k+1    17-------18-------19-------20    |
!            |     |                    |     |
!            |     13-------14-------15-|-----16
!            |    /        /        /   |    /
!            |   9--------10-------11---|---12
!            |  /        /        /     |  /    
!            | 5--------6--------7------|-8
!            |/        /        /       |/
!     k      1--------2--------3--------4
!                                                                       !
!-----------------------------------------------------------------------!
!                             P. FEDE     --  I.M.F.T. --  31/03/2011   !
!-----------------------------------------------------------------------!

! use dns_dim
! use param_phys

implicit none

!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------

!=====================================================================
! Input Arrays
!=============
!- Mesh of data for interpolation
integer, intent(in) :: NX, NY, NZ
real(kind=8), dimension(NX), intent(in) :: XAIN
real(kind=8), dimension(NY), intent(in) :: YAIN
real(kind=8), dimension(NZ), intent(in) :: ZAIN
  
!- Data field for interpolation
real(kind=8), dimension(NX,NY,NZ), intent(in) :: UAIN
  
  
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
integer :: IP2, JP2, KP2
      

real(kind=8) :: ALFA, BETA, GAMA

real(kind=8) :: FXM1, FX0, FXP1, FXP2 
real(kind=8) :: FYM1, FY0, FYP1, FYP2 
real(kind=8) :: FZM1, FZ0, FZP1, FZP2

real(kind=8) :: DXAINT
real(kind=8) :: DYAINT
real(kind=8) :: DZAINT


!!--------------------------------------------------------------------
!!- Index
integer :: NP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DXAINT = XAIN(2) - XAIN(1)
DYAINT = YAIN(2) - YAIN(1)
DZAINT = ZAIN(2) - ZAIN(1)


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
 if(I0 == NX) IP1 = NX - 1
 
 IP2 = IP1 + 1
 if(IP1 == NX) IP2 = NX - 1
 

 IM1 = I0 - 1
 if(I0 == 1) IM1 = 2
 
!---------------------------------------------------------------------
! y-direction
!---------------------------------------------------------------------
 JM1 = J0 - 1
 if (J0 == 1) JM1 = 2 
 JP1 = J0 + 1
 if (J0 == NY) JP1 = NY - 1
 JP2 = JP1 + 1
 if (JP1 == NZ) JP2 = NY - 1

!---------------------------------------------------------------------
! z-direction
!---------------------------------------------------------------------
 KM1 = K0 - 1
 if (K0 == 1) KM1 = 2
 KP1 = K0 + 1
 if (K0 == NZ) KP1 = NZ - 1
 KP2 = KP1 + 1
 if (KP1 == NZ) KP2 = NZ - 1


!!====================================================================
!! 2. Distance to the first node location
!!====================================================================
 ALFA = (XINT(NP)-XAIN(I0)) / DXAINT
 BETA = (YINT(NP)-YAIN(J0)) / DYAINT
 GAMA = (ZINT(NP)-ZAIN(K0)) / DZAINT


!!====================================================================
!! 3. Lagrange polynomial functions
!!====================================================================
 FXM1 = -1./6.*(ALFA - 2.)*(ALFA - 1.)*ALFA
 FX0  =  1./2.*(ALFA - 2.)*(ALFA - 1.)     *(ALFA + 1.)
 FXP1 = -1./2.*(ALFA - 2.)            *ALFA*(ALFA + 1.)
 FXP2 =  1./6.            *(ALFA - 1.)*ALFA*(ALFA + 1.)


 FYM1 = -1./6.*(BETA - 2.)*(BETA - 1.)*BETA
 FY0  =  1./2.*(BETA - 2.)*(BETA - 1.)     *(BETA + 1.)
 FYP1 = -1./2.*(BETA - 2.)            *BETA*(BETA + 1.)
 FYP2 =  1./6.            *(BETA - 1.)*BETA*(BETA + 1.)


 FZM1 = -1./6.*(GAMA - 2.)*(GAMA - 1.)*GAMA
 FZ0  =  1./2.*(GAMA - 2.)*(GAMA - 1.)     *(GAMA + 1.)
 FZP1 = -1./2.*(GAMA - 2.)            *GAMA*(GAMA + 1.)
 FZP2 =  1./6.            *(GAMA - 1.)*GAMA*(GAMA + 1.)


!!====================================================================
!! 4.Interpolation
!!====================================================================
 UINT(NP) = 0.


 UINT(NP) =                                        &
   FZM1*(  FYM1*(   FXM1 * UAIN(IM1,JM1,KM1)       &
                  + FX0  * UAIN(I0 ,JM1,KM1)       &
                  + FXP1 * UAIN(IP1,JM1,KM1)       &
                  + FXP2 * UAIN(IP2,JM1,KM1)  )    &
         + FY0 *(   FXM1 * UAIN(IM1,J0 ,KM1)       &
          + FX0  * UAIN(I0 ,J0 ,KM1)       &
                  + FXP1 * UAIN(IP1,J0 ,KM1)       &
                  + FXP2 * UAIN(IP2,J0 ,KM1)  )    &
         + FYP1*(   FXM1 * UAIN(IM1,JP1,KM1)       &
          + FX0  * UAIN(I0 ,JP1,KM1)       &
          + FXP1 * UAIN(IP1,JP1,KM1)       &
          + FXP2 * UAIN(IP2,JP1,KM1)  )    &
         + FYP2*(   FXM1 * UAIN(IM1,JP2,KM1)       &
          + FX0  * UAIN(I0 ,JP2,KM1)       &
          + FXP1 * UAIN(IP1,JP2,KM1)       &
          + FXP2 * UAIN(IP2,JP2,KM1)  )  ) &

!
 + FZ0 *(  FYM1*(   FXM1 * UAIN(IM1,JM1,K0 )       &
                  + FX0  * UAIN(I0 ,JM1,K0 )       &
                  + FXP1 * UAIN(IP1,JM1,K0 )       &
                  + FXP2 * UAIN(IP2,JM1,K0 )  )    &
     + FY0 *(   FXM1 * UAIN(IM1,J0 ,K0 )       &
              + FX0  * UAIN(I0 ,J0 ,K0 )       &
                  + FXP1 * UAIN(IP1,J0 ,K0 )       &
              + FXP2 * UAIN(IP2,J0 ,K0 )  )    &
     + FYP1*(   FXM1 * UAIN(IM1,JP1,K0 )       &
          + FX0  * UAIN(I0 ,JP1,K0 )       &
          + FXP1 * UAIN(IP1,JP1,K0 )       &
          + FXP2 * UAIN(IP2,JP1,K0 )  )    &
     + FYP2*(   FXM1 * UAIN(IM1,JP2,K0 )       &
          + FX0  * UAIN(I0 ,JP2,K0 )       &
          + FXP1 * UAIN(IP1,JP2,K0 )       &
          + FXP2 * UAIN(IP2,JP2,K0 )  )  ) &
!
 + FZP1*(  FYM1*(   FXM1 * UAIN(IM1,JM1,KP1)       &
                  + FX0  * UAIN(I0 ,JM1,KP1)       &
                  + FXP1 * UAIN(IP1,JM1,KP1)       &
                  + FXP2 * UAIN(IP2,JM1,KP1)  )    &
     + FY0 *(   FXM1 * UAIN(IM1,J0 ,KP1)       &
              + FX0  * UAIN(I0 ,J0 ,KP1)       &
                  + FXP1 * UAIN(IP1,J0 ,KP1)       &
              + FXP2 * UAIN(IP2,J0 ,KP1)  )    &
     + FYP1*(   FXM1 * UAIN(IM1,JP1,KP1)       &
          + FX0  * UAIN(I0 ,JP1,KP1)       &
          + FXP1 * UAIN(IP1,JP1,KP1)       &
          + FXP2 * UAIN(IP2,JP1,KP1)  )    &
     + FYP2*(   FXM1 * UAIN(IM1,JP2,KP1)       &
          + FX0  * UAIN(I0 ,JP2,KP1)       &
          + FXP1 * UAIN(IP1,JP2,KP1)       &
          + FXP2 * UAIN(IP2,JP2,KP1)  )  ) &
!
 + FZP2*(  FYM1*(   FXM1 * UAIN(IM1,JM1,KP2)       &
                  + FX0  * UAIN(I0 ,JM1,KP2)       &
                  + FXP1 * UAIN(IP1,JM1,KP2)       &
                  + FXP2 * UAIN(IP2,JM1,KP2)  )    &
     + FY0 *(   FXM1 * UAIN(IM1,J0 ,KP2)       &
              + FX0  * UAIN(I0 ,J0 ,KP2)       &
                  + FXP1 * UAIN(IP1,J0 ,KP2)       &
              + FXP2 * UAIN(IP2,J0 ,KP2)  )    &
     + FYP1*(   FXM1 * UAIN(IM1,JP1,KP2)       &
          + FX0  * UAIN(I0 ,JP1,KP2)       &
          + FXP1 * UAIN(IP1,JP1,KP2)       &
          + FXP2 * UAIN(IP2,JP1,KP2)  )    &
     + FYP2*(   FXM1 * UAIN(IM1,JP2,KP2)       &
          + FX0  * UAIN(I0 ,JP2,KP2)       &
          + FXP1 * UAIN(IP1,JP2,KP2)       &
          + FXP2 * UAIN(IP2,JP2,KP2)  )  ) 

end do !- on NP = 1, NINT




10000 format (8(A,f6.2,4x))


end subroutine FROZEN_INTERP_LAG3
