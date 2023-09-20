!!====================================================================
!!
!!
!!====================================================================

subroutine MESHING

!!====================================================================
!!
!!
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE
use PARAM_PHYS

implicit none


!---------------------------------------------------------------------
!- Local arrays
!---------------------------------------------------------------------
integer :: I, J, K

real(kind=8) :: K2, K2SCL, DAMPING
!---------------------------------------------------------------------

!!- Flog for the method of aliasing control
!!  ALIASCONTROL = 1 : Sharp aliasing 
!!               = 2 : Smooth aliasing control 
ALIASCONTROL = 2

!- Truncation wavenumber
!RTRUNC =  2.0/3.0   
!RTRUNC =  2.0*sqrt(2.0)/3.0   
RTRUNC =  8.0/9.0   

!!====================================================================
!! 1. COMPUTATIONNAL GRID (physical space)
!!====================================================================

!- Space step
DX = LXMAX / real(NX-1)
DY = LYMAX / real(NY-1)
DZ = LZMAX / real(NZ-1)

!!--------------------------------------------------------------------
!! 1.1. x-direction
!!--------------------------------------------------------------------
!! I=     1       2             NX-1     NX
!!        !-------!---- .... ----!-------!
!! X=     0      DX            L-2DX    L-DX
!!--------------------------------------------------------------------
do I = ISTART(1), IEND(1)
 XMESH(I) = (I-1)*DX
end do

!!--------------------------------------------------------------------
!! 1.2. y-direction
!!--------------------------------------------------------------------
!! J=     1       2             NY-1     NY
!!        !-------!---- .... ----!-------!
!! Y=     0      DY            L-2DY    L-DY
!!--------------------------------------------------------------------
do J = ISTART(2), IEND(2)
 YMESH(J) = DY*(J-1)
end do

!!--------------------------------------------------------------------
!! 1.3. z-direction
!!--------------------------------------------------------------------
!!
!! K=     1       2             NZ-1     NZ
!!        !-------!---- .... ----!-------!
!! Z=     0      DZ            L-2DZ    L-DZ
!!
!!--------------------------------------------------------------------
!! In physical space the slab are oriented along z
!!--------------------------------------------------------------------
do K = ISTART(3), IEND(3)
 ZMESH(K) = DZ*(K-1)
end do



!!====================================================================
!! 2. DEFINE WAVENUMBERS
!!====================================================================
if(SOLVE_FLUID>0) then

!!--------------------------------------------------------------------
!! 2.1. x-direction
!!--------------------------------------------------------------------
!!
!!  I = 1       2       NX/2+1
!!      !-------!---...---!
!! KX = 0       1       NX/2
!!
!!--------------------------------------------------------------------
do I = FSTART(1), FEND(1)
 KX(I) = (I-1)*KFIRST
end do

!!--------------------------------------------------------------------
!! 2.2. y-direction
!!--------------------------------------------------------------------
!!
!!  J = 1       2       NY/2+1  NY/2+2         NY-1     NY
!!      !-------!---...---!-------!-----...-----!-------!
!! Ky = 0       1       NY/2   -NY/2+1         -2      -1
!
!!--------------------------------------------------------------------
do J = FSTART(2), FEND(2)
 if(J<=NY/2+1) then
  KY(J) = (J-1)*KFIRST
 else
  KY(J) = (J-1-NY)*KFIRST
 end if
end do


!!--------------------------------------------------------------------
!! 2.3. z-direction
!!--------------------------------------------------------------------
!!
!!  K = 1       2       NZ/2+1  NZ/2+2         NZ-1     NZ
!!      !-------!---...---!-------!-----...-----!-------!
!! KZ = 0       1       NZ/2   -NZ/2+1         -2      -1
!!
!!--------------------------------------------------------------------
do K = FSTART(3), FEND(3)
 if(K<=NZ/2+1) then
  KZ(K) = (K-1)*KFIRST
 else
  KZ(K) = (K-1-NZ)*KFIRST
 end if
enddo


end if !!- - end if SOLVE_FLUID>0)

!!--------------------------------------------------------------------
!! 3. Mask for aliasing control
!!--------------------------------------------------------------------
if (SOLVE_FLUID ==1) then

!!- Initiation
FILTER(:,:,:) = ZERO

    
do K = FSTART(3),FEND(3)
 do J = FSTART(2),FEND(2)
  do I = FSTART(1),FEND(1)


   K2 = KX(I)**2 + KY(J)**2 + KZ(K)**2

   K2SCL= (KX(I)/KLAST)**2 + (KY(J)/KLAST)**2 + (KZ(K)/KLAST)**2

!!- Default no aliasing control
   FILTER(I,J,K) = 1.0

!!- Sharp aliasing control
   if(ALIASCONTROL == 1) then

    if( K2SCL > RTRUNC**2) FILTER(I,J,K) = ZERO

!!- Smooth aliasing control
   elseif(ALIASCONTROL == 2) then

    if( (K2SCL > RTRUNC**2) .and. (K2SCL <= 1.0) ) then
     FILTER(I,J,K) = cos( (PPI/2.0)*(sqrt(K2SCL)-RTRUNC)/(1.0-RTRUNC) )**4
    elseif(K2SCL > 1.0) then
     FILTER(I,J,K) = ZERO
    end if

   end if !- end if ALIASCONTROL 


  enddo
 enddo
enddo

end if

if(MYID==0) write(*,*)'Mesh construction --> OK'
if(ALIASCONTROL == 1) then
if(MYID==0) write(*,*)' Aliasing control --> Sharp'
if(MYID==0) write(*,*)'           Rtrunk --> ',RTRUNC
elseif(ALIASCONTROL == 2) then
if(MYID==0) write(*,*)' Aliasing control --> Smooth'
if(MYID==0) write(*,*)'           Rtrunk --> ',RTRUNC
else
if(MYID==0) write(*,*)' Aliasing control --> No dealiasing'
end if



end subroutine MESHING
