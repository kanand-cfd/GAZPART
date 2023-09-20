!!====================================================================
!!
!!          Particle-particle detection grid connectivity
!!
!!====================================================================

subroutine COLLISION_BOX_CONNECTIVITY

!!====================================================================
!!
!!====================================================================
use DNS_DIM, only: MYID, NPROC
use PARAM_PHYS, only: PERIODICITY, WALL_BOUNDARY
use COLLISION_VARIABLE

implicit none


!!====================================================================
! ARRAYS STATEMENT
!---------------------------------------------------------------------

!!--------------------------------------------------------------------
!!- Size of each box
character(len=20) :: FILENAME

!- Index
integer :: N
integer :: IJK, IJKV
integer :: I0, J0, K0
integer :: IP1, JP1, KP1 
integer :: IM1, JM1, KM1

integer, dimension(27) :: COUNTV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!- 
allocate(COLBOX(NBOX))


do I0 = 1, NXBOX

  do J0 = 1, NYBOX

    do K0 = 1, NZBOX

      IJK = I0 + (J0-1)*NXBOX + (K0-1)*NXBOX*NYBOX

      IM1 = I0 - 1
      IP1 = I0 + 1

      JM1 = J0 - 1
      JP1 = J0 + 1

      KM1 = K0 - 1
      KP1 = K0 + 1

      if(IM1<1  .and. PERIODICITY) IM1 = NXBOX
      if(IP1>NXBOX .and. PERIODICITY) IP1 = 1

      
      if(NPROC==1) then

        if(JM1<1    ) JM1 = NYBOX
        
        if(JP1>NYBOX) JP1 = 1

        if(KM1<1    ) KM1 = NZBOX

        if(KP1>NZBOX) KP1 = 1

      end if

      !!-   
      N = 0

      !!===============================================
      !! Cells located KM1
      !!===============================================
      if(KM1>=1) then
        !  if(JM1>=1) then
        !!- Box 1 -> IM1, JM1, KM1
        !   N = N + 1
        !   IJKV = IM1 + (JM1-1)*NXBOX + (KM1-1)*NXBOX*NYBOX
        !   COLBOX(IJK)%CON(N) = IJKV
        
        !!- Box 2 -> I0, JM1, KM1
        !   N = N + 1
        !   IJKV = I0  + (JM1-1)*NXBOX + (KM1-1)*NXBOX*NYBOX
        !   COLBOX(IJK)%CON(N) = IJKV

        !!- Box 3 -> IP1, JM1, KM1
        !   N = N + 1
        !   IJKV = IP1 + (JM1-1)*NXBOX + (KM1-1)*NXBOX*NYBOX
        !   COLBOX(IJK)%CON(N) = IJKV
        !  end if

        !!- Box 4 -> IM1, J0, KM1
        !  N = N + 1
        !  IJKV = IM1 + (J0-1)*NXBOX + (KM1-1)*NXBOX*NYBOX
        !  COLBOX(IJK)%CON(N) = IJKV

        !!- Box 5 -> I0, J0, KM1
        !  N = N + 1
        !  IJKV = I0  + (J0-1)*NXBOX + (KM1-1)*NXBOX*NYBOX
        !  COLBOX(IJK)%CON(N) = IJKV

        !!- Box 6 -> IP1, J0, KM1
        !  N = N + 1
        !  IJKV = IP1 + (J0-1)*NXBOX + (KM1-1)*NXBOX*NYBOX
        !  COLBOX(IJK)%CON(N) = IJKV

        !!- Box 7 -> IM1, JP1, KM1
        if(JP1<=NYBOX) then

          !   N = N + 1
          !   IJKV = IM1 + (JP1-1)*NXBOX + (KM1-1)*NXBOX*NYBOX
          !   COLBOX(IJK)%CON(N) = IJKV

          !!- Box 8 -> I0, JP1, KM1
          !   N = N + 1
          !   IJKV = I0  + (JP1-1)*NXBOX + (KM1-1)*NXBOX*NYBOX
          !   COLBOX(IJK)%CON(N) = IJKV

          !!- Box 9 -> IP1, JP1, KM1
          if(IP1<=NXBOX) then
            
            N = N + 1
            IJKV = IP1 + (JP1-1)*NXBOX + (KM1-1)*NXBOX*NYBOX
            COLBOX(IJK)%CON(N) = IJKV

          end if

        end if

      end if !- end if(KM1>=1)

      !!===============================================
      !! Cells located K0 --> 8 cells
      !!===============================================
      if(JM1>=1) then

      !!- Box 10 -> IM1, JM1, K0
      !   N = N + 1
      !   IJKV = IM1 + (JM1-1)*NXBOX + (K0-1)*NXBOX*NYBOX
      !   COLBOX(IJK)%CON(N) = IJKV

      !!- Box 11 -> I0, JM1, K0
      !   N = N + 1
      !   IJKV = I0 + (JM1-1)*NXBOX + (K0-1)*NXBOX*NYBOX
      !   COLBOX(IJK)%CON(N) = IJKV

      !!- Box 12 -> IP1, JM1, K0
        
        if(IP1<=NXBOX) then

          N = N + 1
          IJKV = IP1 + (JM1-1)*NXBOX + (K0-1)*NXBOX*NYBOX
          COLBOX(IJK)%CON(N) = IJKV
        
        end if
      
      end if

      !!- Box 13 -> IM1, J0, K0
      !  N = N + 1
      !  IJKV = IM1 + (J0-1)*NXBOX + (K0-1)*NXBOX*NYBOX
      !  COLBOX(IJK)%CON(N) = IJKV

      !!- Box 14 -> I0, J0, K0
      !  N = N + 1
      !  IJKV = I0 + (J0-1)*NXBOX + (K0-1)*NXBOX*NYBOX
      !  COLBOX(IJK)%CON(N) = IJKV

      !!- Box 15 -> IP1, J0, K0
      if(IP1<=NXBOX) then

        N = N + 1
        
        IJKV = IP1 + (J0-1)*NXBOX + (K0-1)*NXBOX*NYBOX

        COLBOX(IJK)%CON(N) = IJKV

      end if

      if(JP1<=NYBOX) then

        !!- Box 16 -> IM1, JP1, K0
        !   N = N + 1
        !   IJKV = IM1 + (JP1-1)*NXBOX + (K0-1)*NXBOX*NYBOX
        !   COLBOX(IJK)%CON(N) = IJKV

        !!- Box 17 -> I0, JP1, K0
        N = N + 1
        IJKV = I0 + (JP1-1)*NXBOX + (K0-1)*NXBOX*NYBOX
        COLBOX(IJK)%CON(N) = IJKV

        !!- Box 18 -> IP1, JP1, K0

        if(IP1<=NXBOX) then

          N = N + 1
          
          IJKV = IP1 + (JP1-1)*NXBOX + (K0-1)*NXBOX*NYBOX

          COLBOX(IJK)%CON(N) = IJKV

        end if

      end if

      !!===============================================
      !! Cells located KP1 --> 9 cells
      !!===============================================

      if(KP1<=NZBOX) then

        if(JM1>=1) then

          !!- Box 19 -> IM1, JM1, KP1
          !   N = N + 1
          !   IJKV = IM1 + (JM1-1)*NXBOX + (KP1-1)*NXBOX*NYBOX
          !   COLBOX(IJK)%CON(N) = IJKV

          !!- Box 20 -> I0, JM1, KP1
          N = N + 1

          IJKV = I0  + (JM1-1)*NXBOX + (KP1-1)*NXBOX*NYBOX

          COLBOX(IJK)%CON(N) = IJKV

          !!- Box 21 -> IP1, JM1, KP1

          if(IP1<=NXBOX) then

            N = N + 1

            IJKV = IP1 + (JM1-1)*NXBOX + (KP1-1)*NXBOX*NYBOX

            COLBOX(IJK)%CON(N) = IJKV

          end if

        end if

        !!- Box 22 -> IM1, J0, KP1

        if(IM1>=1) then

          N = N + 1

          IJKV = IM1 + (J0 -1)*NXBOX + (KP1-1)*NXBOX*NYBOX

          COLBOX(IJK)%CON(N) = IJKV

        end if

        !!- Box 23 -> I0, J0, KP1
        N = N + 1

        IJKV = I0  + (J0 -1)*NXBOX + (KP1-1)*NXBOX*NYBOX

        COLBOX(IJK)%CON(N) = IJKV

        !!- Box 24 -> IP1, J0, KP1
        if(IP1<=NXBOX) then

          N = N + 1

          IJKV = IP1 + (J0 -1)*NXBOX + (KP1-1)*NXBOX*NYBOX

          COLBOX(IJK)%CON(N) = IJKV

        end if

        !!- Box 25 -> IM1, JP1, KP1

        if(JP1<=NYBOX) then

          if(IM1>=1) then

            N = N + 1

            IJKV = IM1 + (JP1-1)*NXBOX + (KP1-1)*NXBOX*NYBOX

            COLBOX(IJK)%CON(N) = IJKV

          end if

          !!- Box 26 -> I0, JP1, KP1
          N = N + 1

          IJKV = I0  + (JP1-1)*NXBOX + (KP1-1)*NXBOX*NYBOX

          COLBOX(IJK)%CON(N) = IJKV

          !!- Box 27 -> IP1, JP1, KP1

          if(IP1<=NXBOX) then

            N = N + 1

            IJKV = IP1 + (JP1-1)*NXBOX + (KP1-1)*NXBOX*NYBOX

            COLBOX(IJK)%CON(N) = IJKV

          end if

        end if

      end if !- end if(KP1<=NZBOX)

      !!==================================================

      COLBOX(IJK)%NCON = N

    end do

  end do

end do

COUNTV(:) = 0

do IJK = 1, NXBOX*NYBOX*NZBOX

  do N=1, 27

    if(N==COLBOX(IJK)%NCON)  COUNTV(N) = COUNTV(N) +1

  end do

end do


if(MYID==0) then

  write(*,*)' Connectivity for collision detection grid created'

  write(*,*)'  + Full number of boxes = ',NBOX
  
  do N = 1, 27

    if(COUNTV(N)/=0)write(*,'(A,I3,A,I8)')' + Number of cell with ',N,' Neighbors =',COUNTV(N)

  end do

  write(*,*)' + Number of cell with Neighbors =',sum(COUNTV(:))
  write(*,*)

end if


!!- 
!write(FILENAME,10101)'gridcol_',MYID,'.dat'
!open(unit=300,file=trim(FILENAME))
!write(300,2000)
!write(300,2010)NXBOX,NYBOX,NZBOX,1.0
!do K0 = 1, NZBOX
! do J0 = 1, NYBOX
!  do I0 = 1, NXBOX
!  IJK = I0 + (J0-1)*NXBOX + (K0-1)*NXBOX*NYBOX
!  write(300,10000)XMINBOX+DXBOX/2*I0, &
!                  YMINBOX+DYBOX/2*J0, &
!                  ZMINBOX+DZBOX/2*K0, &
!                  real(COLBOX(IJK)%NCON)
!  end do
! end do
!end do
!close(300)

!!--------------------------------------------------------------------
10000 format (10(e17.7))
10100 format (A,I2.2,A)
10101 format (A,I2.2,A,A)
10102 format (A,I2.2,A,I8)
10300 format (A,A)
10303 format (A,A,A)

2000 format ('VARIABLES = "x", "y", "z", "nbv"')
2010 format ('ZONE F=POINT I=',i4,' J=',i4,' K=',i4,' SOLUTIONTIME=',E13.6)

end subroutine COLLISION_BOX_CONNECTIVITY
