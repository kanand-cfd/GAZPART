subroutine LOCATE_PART(J,PARTICLE,DPARTICLE,NS)

use DNS_DIM            !- Dimension
use PARAM_PHYS         !- Physical & numerical parameters
use PARTICLE_PARALLEL
use GEOMETRIC_VARIABLE !- Mesh parameters
use STATISTICS         !- Statistics
use CHECK_CPU 
use FLUID_VARIABLE        


implicit none

!------------------------------------------
!------------------------------------------
integer, intent(in) :: J
type(PARTTYPE), intent(in) :: PARTICLE
real(kind=8), intent(in) :: DPARTICLE
integer     , intent(out) :: NS
!------------------------------------------
!------------------------------------------

real(kind=8) :: DELX, XPART

real(kind=8) :: ZZ1, ZZ2

integer:: N, IDP

!------------------------------------------

DELX = (LXMAX-DPARTICLE)/real(NSLICE-1)
XPART = PARTICLE%XP

IDP = PARTICLE%IDP

!!- Uniform mesh
if(REFINE_MESH .le. 1.0) then

    NS = int((XPART-0.5*DPARTICLE)/DELX) + 1

else
!!- Non-uniform mesh
    ZZ1 = 0.5*DPARTICLE
    ZZ2 = LXMAX - 0.5*DPARTICLE

    if (XPART<=((ZZ1+ZZ2)/2)) then

        NS=1+int(log(1+(REFINE_MESH-1)*(XPART-ZZ1)/(XSTAT(J,IDP,2)-ZZ1))/log(REFINE_MESH))

    else

        NS=floor(NSLICE+1-log(1+(REFINE_MESH-1)*(ZZ2-XPART)/(XSTAT(J,IDP,2)-ZZ1))/log(REFINE_MESH))

    end if

    if (NS == 0) NS=1
    if (NS == NSLICE+1) NS = NSLICE

end if



end subroutine LOCATE_PART
