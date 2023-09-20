module COLLISION_VARIABLE

 implicit none 

 real (kind=8) :: DELTA_COLL
 integer :: NPNEIGH_COLL

!- Inter-particle restitution coefficient
 real(kind=8), dimension(:), allocatable :: ECP
 real(kind=8), dimension(:), allocatable :: MUP 

!- Particle Wall Coefficient of Restitution
real(kind=8) :: ECW 
!- Particle Wall Coefficient of Friction Normal 
real(kind=8) :: MUW
!- Particle Wall Coefficient of Friction Tangential
real(kind=8) :: BETAW

!- Number of detection boxes
 integer :: NXBOX
 integer :: NYBOX
 integer :: NZBOX

 integer :: NBOX

!- Size of the domain including the halo
 real(kind=8) :: XMINBOX, XMAXBOX, DXBOX
 real(kind=8) :: YMINBOX, YMAXBOX, DYBOX
 real(kind=8) :: ZMINBOX, ZMAXBOX, DZBOX

!- Connectivity for the detection boxes
 type cell_box
  integer ::  NCON
  integer,dimension(27) :: CON
 end type cell_box

!- Structure for collision detection
 type(cell_box), dimension(:), allocatable ::  COLBOX


!- Statistics of inter-particle collision
 integer, parameter :: NPDFMAX = 40
 real(kind=8), dimension(:,:,:,:,:), allocatable :: PDFCOL
 real(kind=8), dimension(5) :: DPDFCOL
 real(kind=8)  :: NMOM
 real(kind=8), dimension(:,:,:), allocatable :: MEAN_TIME_PART_COL

 contains
  recursive subroutine QUICKSORT(NBROW,HOME,NEAR,D)

   !!====================================================================
   !!
   !!          QUICKSORT
   !!
   !! We would like to sort the Home and Near structures using the distance between the particles 
   !! to be sure that if a processor X does a collision between a particle A (which have the ID X)
   !! and a particle B (which have the ID Y) (the particles are in the halo), then the processor Y
   !! will also make the collision between the particles A and B. If the Home and Near structures
   !! aren't sorted, the processor Y can do a collision between the particle B and an other 
   !! particle C first. Considering the fact that we allowed only one collision per particle 
   !! per time step, the Y processor will not do the A-B collision because the B particle had
   !! already collide with an other particle. So the X processor have done the A-B collision but 
   !! the Y processor didn't ! In consequence, the B particle velocity is not changed, 
   !! only the A particle velocity is changed.
   !!====================================================================

   implicit none 

   !!====================================================================

   !!- Number of particle "close" (nb de lignes du tableau)(dans les appels récursifs,
                                                        ! il ne s'agira que d'une 
                                                        ! partie des particule proches)
   integer,                    intent(in) :: NBROW

   !data structure D (column 1 : index particle Home
   !                     column 2 : index particle Near
   !                     column 3 : distance between the particle Home and the particle Near)
   integer,      dimension(NBROW), intent(inout) :: HOME, NEAR
   real(kind=8), dimension(NBROW), intent(inout) :: D

   ! Pivot 
   real(kind=8) :: PIVOT

   ! temporary variable for exchanges
   integer :: ITMP1
   integer :: ITMP2
   real(kind=8) :: RTMP

   ! compteur début pour les échanges
   integer :: J
   ! compteur fin pour les échanges
   integer :: I
   
   integer :: K ! INDEX
   !!====================================================================

   ! vérifications :
   !if (NBROW==1) then
   !  write(*,*)'NBROW = ',NBROW
   !end if

   !vérifications
   !write(*,*)'NBROW = ',NBROW,'MYID = ',MYID
   !write(*,*)'distances : ',D(:,3)

    if (NBROW >= 2) then
     
     !!! Cas NBROW > 2 !!!
     !--------------------  
     if (NBROW > 2) then
        !!====================================================================
        !! PRELIMINARY
        !!====================================================================

        !On prend pour pivot la première distance (le tableau n'est à priori pas partiellement trié)
        PIVOT = D(1)

        !On place la ligne du pivot à la fin
        ITMP1 = HOME(NBROW)
        ITMP2 = NEAR(NBROW)
        RTMP = D(NBROW)

        D(NBROW) = D(1)
        HOME(NBROW) = HOME(1)
        NEAR(NBROW) = NEAR(1)
        
        HOME(1)=ITMP1
        NEAR(1)=ITMP2
        D(1)=RTMP

        J = 1  ! compteur début pour les échanges
        I = NBROW-1  ! compteur fin pour les échanges

        !!====================================================================
        !! SORT
        !!====================================================================

        do K=1, NBROW-1

          if (D(J)>=PIVOT) then
            ITMP1=HOME(I)
            ITMP2=NEAR(I)
            RTMP = D(I)

            HOME(I) = HOME(J)
            NEAR(I) = NEAR(J)
            D(I) = D(J)

            HOME(J) = ITMP1
            NEAR(J) = ITMP2
            D(J) = RTMP

            I = I-1
          else
            J = J+1
          end if
        end do

        !!====================================================================
        !! REPLACEMENT DU PIVOT
        !!====================================================================
        ITMP1 = HOME(NBROW)
        ITMP2 = NEAR(NBROW)
        RTMP = D(NBROW)

        HOME(NBROW) = HOME(I+1)
        NEAR(NBROW) = NEAR(I+1)
        D(NBROW) = D(I+1)  

        HOME(I+1) = ITMP1
        NEAR(I+1) = ITMP2
        D(I+1) = RTMP

        
        !vérifications 
        !write(*,*)'distances : '
        !write(*,*)D(:,3)
        !!====================================================================
        !! RECURSIVITY
        !!====================================================================
        if ((I+1)>1) then
        
          call QUICKSORT(I+1,HOME(:I+1),NEAR(:I+2),D(:I+1)) ! partie avant le pivot (pivot inclu)  
        
        end if

        if (NBROW-(I+1)>1) then
        
          call QUICKSORT(NBROW-(I+1),HOME(I+2:),NEAR(I+2:),D(I+2:)) ! partie après le pivot (pivot exclu) 
        
        end if 
     
     end if ! end if NBROW > 2

     !!! Cas NBROW = 2 !!!
     !-------------------- 
     if (NBROW==2) then
      
      if (D(1) /= D(2)) then

        !!====================================================================
        !! PRELIMINARY
        !!====================================================================

        !On prend pour pivot la première distance (le tableau n'est à priori pas partiellement trié)
        PIVOT = D(1)

        !On place la ligne du pivot à la fin
        ITMP1 = HOME(NBROW)
        ITMP2 = NEAR(NBROW)
        RTMP = D(NBROW)

        D(NBROW) = D(1)
        HOME(NBROW) = HOME(1)
        NEAR(NBROW) = NEAR(1)
        
        HOME(1)=ITMP1
        NEAR(1)=ITMP2
        D(1)=RTMP


        J = 1  ! compteur début pour les échanges
        I = NBROW-1  ! compteur fin pour les échanges

        !!====================================================================
        !! SORT
        !!====================================================================

        do K=1, NBROW-1
          
          if (D(J)>=PIVOT) then
            
            ITMP1=HOME(I)
            ITMP2=NEAR(I)
            RTMP = D(I)

            HOME(I) = HOME(J)
            NEAR(I) = NEAR(J)
            D(I) = D(J)

            HOME(J) = ITMP1
            NEAR(J) = ITMP2
            D(J) = RTMP


            I = I-1
          else
            
            J = J+1
          
          end if
        
        end do

        !!====================================================================
        !! REPLACEMENT DU PIVOT
        !!====================================================================
        ITMP1 = HOME(NBROW)
        ITMP2 = NEAR(NBROW)
        RTMP = D(NBROW)

        HOME(NBROW) = HOME(I+1)
        NEAR(NBROW) = NEAR(I+1)
        D(NBROW) = D(I+1)  

        HOME(I+1) = ITMP1
        NEAR(I+1) = ITMP2
        D(I+1) = RTMP

    
        !vérifications 
        !write(*,*)'distances : '
        !write(*,*)D(:,3)
        !!====================================================================
        !! RECURSIVITY
        !!====================================================================
        if ((I+1)>1) then
          
          call QUICKSORT(I+1,HOME(:I+1),NEAR(:I+1),D(:I+1)) ! partie avant le pivot (pivot inclu)  
        
        end if

        if (NBROW-(I+1)>1) then
          
          call QUICKSORT(NBROW-(I+1),HOME(I+2:),NEAR(I+2:),D(I+2:)) ! partie après le pivot (pivot exclu) 
        
        end if 

      end if   ! end if (D(1,3) /= D(2,3))
     
     end if   ! end if NBROW==2 
        
    end if  ! end if NBROW >=2

  end subroutine QUICKSORT


end module COLLISION_VARIABLE


