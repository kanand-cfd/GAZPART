subroutine INIT_CHANNEL_STAT

use DNS_DIM            !- Dimension
use PARAM_PHYS         !- Physical & numerical parameters
use PARTICLE_PARALLEL
use GEOMETRIC_VARIABLE !- Mesh parameters
use STATISTICS         !- Statistics
use CHECK_CPU 
use FLUID_VARIABLE        


implicit none
!-----------------------------------------------------------
integer :: I, J, NS, IDP
real(kind=8) :: ZZ1, ZZ2
real(kind=8) :: DELX, PDU, PDV, PDW
real(kind=8) :: DPART
!-----------------------------------------------------------


!REFINE_MESH = 1.15


do J = 1, NIG

   do IDP = 1, POLYDISP(J)

        DPART = 2.0*EMAJ_PART(J,IDP)/APR_PART(J)
!!=================================
!!- Non-uniform mesh
!!=================================
		if(REFINE_MESH>1.0) then

        ZZ1 = 0.5*DPART
        ZZ2 = LXMAX - 0.5*DPART
        XSTAT(J,IDP,1)=ZZ1
        XSTAT(J,IDP,2)=XSTAT(J,IDP,1)+(ZZ2-ZZ1)/2*(REFINE_MESH-1)/(REFINE_MESH**(NSLICE/2)-1)

            do NS = 1,(NSLICE/2+1)
            
                XSTAT(J,IDP,NS)=(XSTAT(J,IDP,2)-ZZ1)/(REFINE_MESH-1)*(REFINE_MESH**(NS-1)-1)+ZZ1
                XSTAT(J,IDP,NSLICE+2-NS)=ZZ2-XSTAT(J,IDP,NS)+ZZ1
                !write(*,*) NS, XSTAT(J,IDP,NS)
                !write(*,*) NSLICE+2-NS, XSTAT(J,IDP,NSLICE+2-NS)

		    end do

!!================================
!! Uniform mesh
!!================================
		else

            DELX = (LXMAX-DPART)/real(NSLICE-1)
 
            XSTAT(J,IDP,1) = 0.5*DPART

            do NS = 2 , NSLICE+1
                XSTAT(J,IDP,NS)=XSTAT(J,IDP,NS-1) + DELX
 			end do

		end if

        do NS = 1,NSLICE

            DEL_X(J,IDP,NS) = XSTAT(J,IDP,NS+1) - XSTAT(J,IDP,NS)
            !write(*,*) XSTAT(J,IDP,NS) + 0.5*DEL_X(J,IDP,NS)

        end do

	end do

end do


! Initiation of PDF wall
! DPDFWALL(1) = 3.5 / (0.5*real(NPDF_MAX))
! DPDFWALL(2) = 3.5 / (0.5*real(NPDF_MAX))
! DPDFWALL(3) = 3.5 / (0.5*real(NPDF_MAX))
! allocate(PDFWALL(NIG,POLYDISPMAX,NPDF_MAX,5))
! allocate(MEAN_PDFWALL(NIG,POLYDISPMAX,NPDF_MAX,5))
! allocate(DPDFWALL(NSTAT,NPDF_MAX))

! PDFWALL(:,:,:,:) = ZERO
! MEAN_PDFWALL(:,:,:,:) = ZERO

! PDU = 7.0/real(NPDF_MAX)
! PDV = 3.0/real(NPDF_MAX)
! PDW = 20.0/real(NPDF_MAX)

! DPDFWALL(1,1) = -3.5
! DPDFWALL(2,1) = 0.0
! DPDFWALL(3,1) = 0.0

! do I = 2, NPDF_MAX

!     DPDFWALL(1,I) = DPDFWALL(1,I-1) + PDU
!     DPDFWALL(2,I) = DPDFWALL(2,I-1) + PDV
!     DPDFWALL(3,I) = DPDFWALL(3,I-1) + PDW    
    
! end do


end subroutine INIT_CHANNEL_STAT