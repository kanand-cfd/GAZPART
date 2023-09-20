!!====================================================================
!!
!!          Particle-particle detection
!!
!!====================================================================

subroutine COLLISION_DETECTION_CONNECTIVITY(NPART, &
                               XPART, &
                               YPART, &
                               ZPART, &
                              NCLOSE, &
                           NCLOSEMAX, &
                                HOME, &
                                NEAR  )

!!====================================================================
!!
!!====================================================================
use DNS_DIM, only: MYID, NPROC

use COLLISION_VARIABLE

implicit none

!!====================================================================
! ARRAYS STATEMENT
!---------------------------------------------------------------------
integer,                        intent(in) :: NPART
real(kind=8), dimension(NPART), intent(in) :: XPART
real(kind=8), dimension(NPART), intent(in) :: YPART
real(kind=8), dimension(NPART), intent(in) :: ZPART


!!- Number of particle close
integer,                     intent(out) :: NCLOSE
integer,                     intent(in ) :: NCLOSEMAX

!!- index of home particle
integer, dimension(NCLOSEMAX), intent(out) :: HOME

!!- index of near particle
integer, dimension(NCLOSEMAX), intent(out) :: NEAR

!!--------------------------------------------------------------------

integer, dimension(NBOX) :: NPBOX
integer, dimension(:,:), allocatable :: PARTBOX

integer :: NPCON

!- Index
integer :: NP, I, J, K, IJK, I0, J0, K0, NP0, NP1, NP2
integer :: NV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

HOME(:) = 0
NEAR(:) = 0

NCLOSE = 0



!!--------------------------------------------------------------------
!!- Number of particle in each box
!!--------------------------------------------------------------------
NPBOX(:) = 0

do NP = 1, NPART

  !!- index of the box containing the particle
  !write(*,*) DXBOX, DYBOX, DZBOX
  !write(*,*) NP, XPART(NP), XMINBOX
 
  I0 = int((XPART(NP) - XMINBOX)/DXBOX) + 1
  J0 = int((YPART(NP) - YMINBOX)/DYBOX) + 1 
  K0 = int((ZPART(NP) - ZMINBOX)/DZBOX) + 1

  IJK = I0 + (J0-1)*NXBOX + (K0-1)*NXBOX*NYBOX


  if(I0<1.or.I0>NXBOX.or.J0<1.or.J0>NYBOX.or.K0<1.or.K0>NZBOX) then

    write(*,*)'I0=',I0,' J0=', J0,' K0=',K0,' IJK=',IJK
    write(*,*)'np=',NP
    write(*,*)'xp=',XPART(NP),YPART(NP),ZPART(NP)
    write(*,*)'dx=',DXBOX, DYBOX, DZBOX
    write(*,*)'xmi=',XMINBOX, YMINBOX,ZMINBOX

  end if

  NPBOX(IJK) = NPBOX(IJK) + 1

end do

NCLOSE = maxval(NPBOX(:))

allocate(PARTBOX(NBOX,NCLOSE))

!write(*,*)'Id=',MYID,' PARTBOX allocated with:',NCLOSE

NCLOSE = 0
NPBOX(:) = 0
PARTBOX(:,:) = 0

do NP = 1, NPART

  !!- index of the box containing the particle
  I0 = int((XPART(NP) - XMINBOX)/DXBOX) + 1
  J0 = int((YPART(NP) - YMINBOX)/DYBOX) + 1 
  K0 = int((ZPART(NP) - ZMINBOX)/DZBOX) + 1

  IJK = I0 + (J0-1)*NXBOX + (K0-1)*NXBOX*NYBOX

  NPBOX(IJK) = NPBOX(IJK) + 1

  PARTBOX(IJK,NPBOX(IJK)) = NP

end do

!do I=1,NXBOX
! do J=1,NYBOX
!  do K=1,NZBOX
!   write(*,*)I,J,K,NPBOX(I,J,K),(PARTBOX(I,J,K,NP),NP=1,NPBOX(I,J,K))
!  end do
! end do
!end do
!write(*,'(2(A,I5))')'MYId=',MYID,' Detection     Np =',NPART
!write(*,'(2(A,I5))')'MYId=',MYID,' Detection sum(Np)=',sum(NPBOX(:))
!write(*,'(2(A,I5))')'MYId=',MYID,' Detection max(Np)=',maxval(NPBOX(:))
!write(*,'(2(A,I5))')'MYId=',MYID,' Detection min(Np)=',minval(NPBOX(:))

!do IJK=1, NBOX
!write(*,'(A,I4,A,I4,A,10(I4))')'IJK=',IJK,' Np=',NPBOX(IJK),' PARTBOX = ',PARTBOX(IJK,:)
!end do


!write(*,*)
!write(*,*)
!write(*,*)


!!--------------------------------------------------------------------
!!- List of particles in each box
!!--------------------------------------------------------------------
NP0 = 0

do IJK = 1, NBOX

  NPCON = 0

  NPCON = NPCON + NPBOX(IJK)

  do NV = 1, COLBOX(IJK)%NCON

    NPCON = NPCON + NPBOX(COLBOX(IJK)%CON(NV))

  end do

  !write(*,'(3(A,I4))')'MYID=',MYID,' IJK=',IJK,' NPCON=',NPCON

  if(NPBOX(IJK) /= 0) then

    !write(*,'(A,I4,A,I4)')'IJK=',IJK,' Np=',NPBOX(IJK)

    !!- all particles in box IJK is a potential partner
    do NP1 = 1, NPBOX(IJK)

      !- partners are  the others particle in the cell
      do NP2 = NP1+1, NPBOX(IJK)

        NP0 = NP0 + 1
        
        HOME(NP0) = PARTBOX(IJK,NP1)

        NEAR(NP0) = PARTBOX(IJK,NP2)

        !write(MYID+700,'(A,I4,A,I4,A,I4)')'NP0=',NP0,' HOME=',HOME(NP0),' NEAR=',NEAR(NP0)

      end do

      !- and the one in the neighboring cells
      do NV = 1, COLBOX(IJK)%NCON

        !  write(*,'(A,I4,A,I4)')'V: IJKV=',COLBOX(IJK)%CON(NV),' Np=',NPBOX(COLBOX(IJK)%CON(NV))
        if(NPBOX(COLBOX(IJK)%CON(NV)) /= 0)  then

          do NP2 = 1, NPBOX(COLBOX(IJK)%CON(NV))

            NP0 = NP0 + 1
            
            !     write(*,*) PARTBOX(IJK,NP1), NP0, HOME(NP0)
            HOME(NP0) = PARTBOX(IJK,NP1)
            NEAR(NP0) = PARTBOX(COLBOX(IJK)%CON(NV),NP2)

            !write(MYID+700,'(A,I4,A,I4,A,I4)')'NP0=',NP0,' HOME=',HOME(NP0),' NEAR=',NEAR(NP0)

          end do

        end if

      end do !- loop  NV = 1, COLBOX(IJK)%NCON

    end do !- loop NP1 = 1, NPBOX(IJK)

  end if !- if  NPBOX(I0,J0,K0) /= 0

end do !- loop do IJK = 1, NBOX


NCLOSE = NP0


!write(*,*)'MYID=',MYID,' NCLOSE=',NCLOSE
!do NP=1,NCLOSE
!write(*,'(A,I4,A,I4,A,I4)')'Np=',NP,' home=',HOME(NP),' near=',NEAR(NP)
!end do
!pause

deallocate(PARTBOX)




!!--------------------------------------------------------------------
10100 format (A,I2.2,A)
10101 format (A,I2.2,A,A)
10102 format (A,I2.2,A,I8)
10300 format (A,A)
10303 format (A,A,A)


end subroutine COLLISION_DETECTION_CONNECTIVITY
