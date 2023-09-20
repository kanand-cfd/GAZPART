!!=====================================================================
!!
!!   Arrays allocation
!!
!!=====================================================================
subroutine ALLOCATE_ARRAYS
!!=====================================================================
!!
!!
!!=====================================================================

use DNS_DIM
use FLUID_VARIABLE
use GEOMETRIC_VARIABLE
use PARAM_PHYS 
use FORCING
use STATISTICS
use WORK_ARRAYS
use SCALAR_VARIABLE

use PARTICLE_PARALLEL

implicit none


!---------------------------------------------------------------------
integer :: K
!---------------------------------------------------------------------

!!=====================================================================
!! 1. Allocate arrays for the Fluid DNS
!!=====================================================================

!!---------------------------------------------------------------------
!! 1.1. Mesh and Wavenumber
!!---------------------------------------------------------------------

!- Allocate arrays for mesh
allocate(XMESH(ISTART(1):IEND(1)))
allocate(YMESH(ISTART(2):IEND(2)))
allocate(ZMESH(ISTART(3):IEND(3)))


if(SOLVE_FLUID >0) then
 !- Wavenumbers
 allocate(KX(FSTART(1):FEND(1)))
 allocate(KY(FSTART(2):FEND(2)))
 allocate(KZ(FSTART(3):FEND(3)))
end if

!!---------------------------------------------------------------------
!! 1.2. Arrays for the fluid
!!---------------------------------------------------------------------
if(SOLVE_FLUID >0) then

!- Fluid velocity in real space
allocate(UFLU(ISTART(1)         :IEND(1)             &
             ,ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL    &
             ,ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL )  )

allocate(VFLU(ISTART(1)         :IEND(1)             &
             ,ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL    &
             ,ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL )  )

allocate(WFLU(ISTART(1)         :IEND(1)             &
             ,ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL    &
             ,ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL)   )

!- Fluid velocity in Fourier space
allocate( UFOU(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)),  &
          VFOU(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)),  &
          WFOU(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3))   )


allocate(TMPPHY(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)))
allocate(TMPFOU(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)))



!!- Case of Navier-Stokes equations
if(SOLVE_FLUID==1) then

!- Integrating factor
 allocate(FILTER(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)))


end if !- end if: SOLVE_FLUID==1



end if !- end if: SOLVE_FLUID >0

if(FROZEN_FLOW > 0) then

  allocate(UFLU(ISTART(1)         :IEND(1)             &
             ,ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL    &
             ,ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL )  )

  allocate(VFLU(ISTART(1)         :IEND(1)             &
             ,ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL    &
             ,ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL )  )

  allocate(WFLU(ISTART(1)         :IEND(1)             &
             ,ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL    &
             ,ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL)   )


  allocate(VORTICITY_X(ISTART(1):IEND(1)                    &
                    ,ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL    &
                    ,ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL)   )

  allocate(VORTICITY_Y(ISTART(1):IEND(1)    &
                    ,ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL    &
                    ,ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL)   )

  allocate(VORTICITY_Z(ISTART(1):IEND(1)    &
                    ,ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL    &
                    ,ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL)   )

end if

!!---------------------------------------------------------------------
!! 1.4. Arrays for statistics
!!---------------------------------------------------------------------
if(LEVEL_STFLU>0) then
 allocate(MEAN_FLUID(NSTAT))
 if(STAT_TIME) allocate(MEAN_TIME_FLUID(NSTAT))
end if

if(LEVEL0_STSCL) then
 allocate(MEAN_SCL(NSTAT))
 if(STAT_TIME) allocate(MEAN_TIME_SCL(NSTAT))
end if

!!--------------------------------------------------------------------
!! Spatial correlation
!!--------------------------------------------------------------------
if(STAT_TIME.and.LEVEL_STFLU==4) then
 DIMSCOR = ISIZE(1)/2+1
 
 allocate(MEAN_RUXLOC(DIMSCOR))
 allocate(MEAN_RVXLOC(DIMSCOR))
 allocate(MEAN_RWXLOC(DIMSCOR))

 MEAN_RUXLOC(:) = ZERO
 MEAN_RVXLOC(:) = ZERO
 MEAN_RWXLOC(:) = ZERO
 
end if

!!=====================================================================
!! 2. Allocate arrays for the Particles
!!=====================================================================
!! The size of particle's array 10% larger than uniform distribution of
!! particle for each CPU.
!! To avoid a crash a specific subroutine check the number particles 
!! per CPU before MPI exchange.
!!---------------------------------------------------------------------

!!---------------------------------------------------------------------
!! 2.1. Arrays for point particles
!!---------------------------------------------------------------------
if(SOLVE_PART) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(NPART_LOC(NIG))
allocate(PART(NPMAX_LOC,NIG))
allocate(NBR_EXCHANGE(NIG))
NBR_EXCHANGE(:) = 0

allocate(NPMAX_CPU(NIG))
NPMAX_CPU(:) = 0

allocate(NPMIN_CPU(NIG))
NPMIN_CPU(:) = 1000000
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end if

!!---------------------------------------------------------------------
!! X.X. Arrays for statistics
!!---------------------------------------------------------------------
if(LEVEL0_STPAR) then

 if(STAT_TIME) then

  allocate(MEAN_TIME_PART(NSTAT,NIG,POLYDISPMAX))
  allocate(MEAN_TIME_PART_PDF(NSTAT, NIG, POLYDISPMAX, NPDF))
  allocate(MEAN_TIME_PART_COL_PDF(NSTAT,NIG,POLYDISPMAX, NPDF))
  allocate(MEAN_TIME_PARTFLUID(NSTAT,NIG,POLYDISPMAX))

 end if


 if (STAT_TIME .and. LEVELX_STPAR) then
  allocate(MEAN_TIME_PART_CHAN(NSTAT,NIG,POLYDISPMAX,NSLICE))
  allocate(MEAN_TIME_COLL_CHAN(NSTAT,NIG,POLYDISPMAX,POLYDISPMAX,NSLICE))
  
 ! allocate(MEAN_TIME_PDF_CHAN(NSTAT,NIG,POLYDISPMAX,NSLICE,NPDF))
 ! allocate(MEAN_TIME_PDF_RELV_CHAN(NSTAT,NIG,POLYDISPMAX,POLYDISPMAX,NSLICE,NPDF))

  allocate(B_COLL(NSTAT,NIG,NPMAX_LOC,POLYDISPMAX,POLYDISPMAX))
  allocate(A_COLL(NSTAT,NIG,NPMAX_LOC,POLYDISPMAX,POLYDISPMAX))

  allocate(XSTAT(NIG,POLYDISPMAX,NSLICE+1))
  allocate(DEL_X(NIG,POLYDISPMAX,NSLICE))
  
 end if
 
 if(LEVEL3_STPAR) then
  allocate(RMAX(NIG))
  allocate(RMIN(NIG))
  allocate(DR(NIG))

  !maximum radius of RDF
  NTEST = NPROC                         !number of test particles
  RMIN  = EMAJ_PART_USER(:)/APR_PART(:) !DPART_USER(:) !minimum distance of RDF is particle diameter
  RMAX  = LXMAX/2.                      !maximum distance of RDF can go up to L/2.
  NR    = 100.                          !number of increments in [DPART,RMAX] interval
  DR    = (RMAX-RMIN)/real(NR)          !interval of RDF (spherical shell thickness)

! RDFs for all particles
  allocate(RDF_LOC(NR,NIG,POLYDISPMAX,POLYDISPMAX))
  allocate(RDFWR_LOC(NR,NIG,POLYDISPMAX,POLYDISPMAX))
  allocate(RDFTHR_LOC(NR,NIG,POLYDISPMAX,POLYDISPMAX))

  RDF_LOC(:,:,:,:)     = ZERO
  RDFWR_LOC(:,:,:,:)   = ZERO
  RDFTHR_LOC(:,:,:,:)  = ZERO

  ! RDFs for approaching particles
  allocate(RDFIN_LOC(NR,NIG,POLYDISPMAX,POLYDISPMAX))
  allocate(RDFWRIN_LOC(NR,NIG,POLYDISPMAX,POLYDISPMAX))
  allocate(RDFTHRIN_LOC(NR,NIG,POLYDISPMAX,POLYDISPMAX))
  RDFIN_LOC(:,:,:,:)    = ZERO
  RDFWRIN_LOC(:,:,:,:)  = ZERO
  RDFTHRIN_LOC(:,:,:,:) = ZERO

  ! RDFs for departing particles
  allocate(RDFOUT_LOC(NR,NIG,POLYDISPMAX,POLYDISPMAX))
  allocate(RDFWROUT_LOC(NR,NIG,POLYDISPMAX,POLYDISPMAX))
  allocate(RDFTHROUT_LOC(NR,NIG,POLYDISPMAX,POLYDISPMAX))
  RDFOUT_LOC(:,:,:,:)    = ZERO
  RDFWROUT_LOC(:,:,:,:)  = ZERO
  RDFTHROUT_LOC(:,:,:,:) = ZERO
 end if
  
end if




!!=====================================================================
!! 3. Allocate arrays for scalar
!!=====================================================================
if(SOLVE_SCALAR.and.(SOLVE_FLUID==1)) then

!- Scalar field in real space
allocate(THETA(ISTART(1)         :IEND(1)             &
              ,ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL    &
              ,ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL )  )

!- Scalar in Fourier space
allocate(THETAFOU(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3))   )

!- Right-Hand-Side of Scalar transport equation
allocate(RHS_SCL(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3),3))



end if



!!=====================================================================
if(MYID==0) write(*,*)'Arrays allocation --> OK'

end subroutine ALLOCATE_ARRAYS
