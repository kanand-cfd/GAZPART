!!====================================================================
!!
!!          Initiation of particle-particle collision
!!
!!====================================================================

subroutine COLLISION_INITIATION

!!====================================================================
!!
!! Here we defined the detection boxes and the connectivity of each
!!
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE
use PARAM_PHYS
use COLLISION_VARIABLE
use PARTICLE_PARALLEL


implicit none


!!====================================================================
! ARRAYS STATEMENT
!---------------------------------------------------------------------

!!- Parameter of grid initiation
!!   INIT_GRID = 1 : DELTA_COLL and NXBOX both imposed
!!             = 2 : DELTA_COLL = DX and DX imposed
integer :: INIT_GRID


real(kind=8) :: ALPHAP, ALPHAP_TOT, DPMAX

!!- Index
!integer :: IM1, JM1, KM1
!integer :: I0,  J0,  K0
!integer :: IP1, JP1, KP1
!integer :: IJK, IJKL

integer :: I, J, IDP


!!====================================================================


!- Maximum diameter
DPMAX = 2.0*maxval(EMAJ_PART(:,:))

!!====================================================================
!! 1. Define the geometry of the collision domain
!!====================================================================

INIT_GRID = 2

if(INIT_GRID ==1) then

 !- Distance of overlaped domain
 DELTA_COLL = 2.0*DPMAX

 !- Boxes number
 NXBOX = 42

else

 !- Desired size of each box
 DXBOX = 2.0*DPMAX

 !- Boxes number
 NXBOX = int(LXMAX/DXBOX)

 !- Compute the true size
 DXBOX = LXMAX/real(NXBOX)

 !- Overlap is one boxe
 DELTA_COLL = DXBOX

end if 


if(NPROC==1) DELTA_COLL = ZERO

!- Size of the domain
XMINBOX = XMESH(ISTART(1))   ! -DELTA_COLL
XMAXBOX = XMESH(  IEND(1))+DX! +DELTA_COLL

YMINBOX = YMESH(ISTART(2))   -DELTA_COLL
YMAXBOX = YMESH(  IEND(2))+DY+DELTA_COLL

ZMINBOX = ZMESH(ISTART(3))   -DELTA_COLL
ZMAXBOX = ZMESH(  IEND(3))+DZ+DELTA_COLL


NYBOX = int((YMAXBOX-YMINBOX)/DXBOX)
NZBOX = int((ZMAXBOX-ZMINBOX)/DXBOX)

DYBOX = (YMAXBOX - YMINBOX)/real(NYBOX)
DZBOX = (ZMAXBOX - ZMINBOX)/real(NZBOX)

NBOX = NXBOX*NYBOX*NZBOX



if(MYID==0) then

    write(*,*)'Collision initiation:'
    write(*,'(2(A,E13.6))')' -> delta_col=',DELTA_COLL,' delta_col/max(dp)=',DELTA_COLL/DPMAX
    write(*,'(3(A,E13.6))')' -> min(xcol)=',XMINBOX, &
                           ' max(xcol)=',XMAXBOX, &
                           ' Lx_col/Lx=',(XMAXBOX-XMINBOX)/LXMAX
    write(*,'(3(A,E13.6))')' -> min(ycol)=',YMINBOX, &
                           ' max(ycol)=',YMAXBOX, &
                           ' Ly_col/Ly=',(YMAXBOX-YMINBOX)/LYMAX
    write(*,'(3(A,E13.6))')' -> min(zcol)=',ZMINBOX, &
                           ' max(zcol)=',ZMAXBOX, &
                           ' Lz_col/Lz=',(ZMAXBOX-ZMINBOX)/LZMAX
    write(*,'(A,I3,A,E13.6,A,E13.6)')' -> Box: Nx=',NXBOX, &
                                     ' dx=',(XMAXBOX-XMINBOX)/NXBOX, &
                                     ' dx/max(dp)=',(XMAXBOX-XMINBOX)/NXBOX/DPMAX
    write(*,'(A,I3,A,E13.6,A,E13.6)')' -> Box: Ny=',NYBOX, &
                                     ' dy=',(YMAXBOX-YMINBOX)/NYBOX, &
                                     ' dy/max(dp)=',(YMAXBOX-YMINBOX)/NYBOX/DPMAX
    write(*,'(A,I3,A,E13.6,A,E13.6)')' -> Box: Nz=',NZBOX,&
                                     ' dz=',(ZMAXBOX-ZMINBOX)/NZBOX,&
                                     ' dz/max(dp)=',(ZMAXBOX-ZMINBOX)/NZBOX/DPMAX

end if



do J = 1, NIG

    ALPHAP = ZERO
    
    do I = 1, NPART_LOC(J)

        IDP = PART(I,J)%IDP
        ALPHAP = ALPHAP + (4.0*PPI/3.0)*(EMAJ_PART(J,IDP)**3)/(LXMAX*LYMAX*LZMAX * APR_PART(J)**2)

    end do

    call RSUMCPU(ALPHAP,ALPHAP_TOT)
    if(MYID==0) write(*,'(1(A,I2.2,A,F5.2,A,E13.6))')' -> C',J,' ec=',ECP(J),' ap=',ALPHAP_TOT

end do

!if(MYID==0) write(*,*)



!!====================================================================
!! 2. Grid detection initiation
!!====================================================================
call COLLISION_BOX_CONNECTIVITY



!!====================================================================
!! 3. Collision statistics initiation
!!====================================================================
if(LEVEL1_STPAR) then 
    
    ! DPDFCOL(1) = .5*PPI / real(NPDFMAX)
    ! DPDFCOL(2) =  3.    / real(NPDFMAX)
    ! DPDFCOL(3) =  3.    / real(NPDFMAX)

    if(STAT_TIME) then

        allocate(MEAN_TIME_PART_COL(NSTAT,POLYDISPMAX,POLYDISPMAX))
        MEAN_TIME_PART_COL(:,:,:) = ZERO

    end if
    
    ! allocate(PDFCOL(NIG,POLYDISPMAX,POLYDISPMAX,NPDFMAX,5))
    ! PDFCOL(:,:,:,:,:) = ZERO

end if




!!--------------------------------------------------------------------
10100 format (A,I2.2,A)
10101 format (A,I2.2,A,A)
10102 format (A,I2.2,A,I8)
10300 format (A,A)
10303 format (A,A,A)


end subroutine COLLISION_INITIATION
