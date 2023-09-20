!!=====================================================================
!!
!! Statistics
!!
!!=====================================================================
module STATISTICS

implicit none


!- Stored number of event
integer :: NEVEN

integer :: NCYCLELAST

!- Global array containing statistics on the fluid
real(kind=8), dimension(:), allocatable :: MEAN_FLUID

!- Time-averaged statistics on the fluid
real(kind=8), dimension(:), allocatable :: MEAN_TIME_FLUID

!- Global array containing statistics on the scalar
real(kind=8), dimension(:), allocatable :: MEAN_SCL

!- Time-averaged statistics on the scalar
real(kind=8), dimension(:), allocatable :: MEAN_TIME_SCL

!- Dissipation
real(kind=8) :: EPS_FLU

!- Dissipation computed from spectrum
real(kind=8) :: EPS_FLU_SPEC

!- Mean gradients
real(kind=8), dimension(3,3,4) :: DUIDXJ


!- Time step statistics
real(kind=8) :: MEAN_DT, MAX_DT, MIN_DT, NEVEN_DT

integer, parameter :: NPDF = 100 

!- Time-averaged statistics on the particles
real(kind=8), dimension(:,:,:), allocatable :: MEAN_TIME_PART
real(kind=8), dimension(:,:,:,:), allocatable :: MEAN_TIME_PART_PDF
real(kind=8), dimension(:,:,:,:), allocatable :: MEAN_TIME_PART_COL_PDF

real(kind=8), dimension(:,:,:), allocatable :: MEAN_TIME_PARTFLUID

!integer :: NBLGRMAX

!- 
!real(kind=8), dimension(:), allocatable :: UPT0_MEAN
!real(kind=8), dimension(:), allocatable :: VPT0_MEAN
!real(kind=8), dimension(:), allocatable :: WPT0_MEAN
!
!real(kind=8), dimension(:), allocatable :: UUPT0_MEAN
!real(kind=8), dimension(:), allocatable :: VVPT0_MEAN
!real(kind=8), dimension(:), allocatable :: WWPT0_MEAN
!
!real(kind=8), dimension(:), allocatable :: UFAPT0_MEAN
!real(kind=8), dimension(:), allocatable :: VFAPT0_MEAN
!real(kind=8), dimension(:), allocatable :: WFAPT0_MEAN
!
!real(kind=8), dimension(:), allocatable :: UUFAPT0_MEAN
!real(kind=8), dimension(:), allocatable :: VVFAPT0_MEAN
!real(kind=8), dimension(:), allocatable :: WWFAPT0_MEAN
!

!- Time of Lagrangian correlation (for variable timestep)
real(kind=8), dimension(:,:), allocatable :: TIME_LGR

real(kind=8), dimension(:,:,:,:), allocatable :: RPX_LOC
real(kind=8), dimension(:,:,:,:), allocatable :: RPY_LOC
real(kind=8), dimension(:,:,:,:), allocatable :: RPZ_LOC

real(kind=8), dimension(:,:,:,:), allocatable :: RFAPX_LOC
real(kind=8), dimension(:,:,:,:), allocatable :: RFAPY_LOC
real(kind=8), dimension(:,:,:,:), allocatable :: RFAPZ_LOC


!!- Fluid velocity spatial correlation
integer :: DIMSCOR
real(kind=8), dimension(:), allocatable :: MEAN_RUXLOC
real(kind=8), dimension(:), allocatable :: MEAN_RVXLOC
real(kind=8), dimension(:), allocatable :: MEAN_RWXLOC


!- size of array for radial distribution function
integer :: NTEST
integer :: NR
real(kind=8), dimension(:), allocatable :: DR
real(kind=8), dimension(:), allocatable :: RMAX
real(kind=8), dimension(:), allocatable :: RMIN



!- Radial distribution functions of neighbouring particles
real(kind=8), dimension(:,:,:,:), allocatable :: RDF_LOC
real(kind=8), dimension(:,:,:,:), allocatable :: RDFWR_LOC
real(kind=8), dimension(:,:,:,:), allocatable :: RDFTHR_LOC

!- Radial distribution functions of approaching particles
real(kind=8), dimension(:,:,:,:), allocatable :: RDFIN_LOC
real(kind=8), dimension(:,:,:,:), allocatable :: RDFWRIN_LOC
real(kind=8), dimension(:,:,:,:), allocatable :: RDFTHRIN_LOC

!- Radial distribution functions of departing particles
real(kind=8), dimension(:,:,:,:), allocatable :: RDFOUT_LOC
real(kind=8), dimension(:,:,:,:), allocatable :: RDFWROUT_LOC
real(kind=8), dimension(:,:,:,:), allocatable :: RDFTHROUT_LOC


!Statistics for Channel
integer :: XFLAG, NEVEN_CHAN, NEVEN_COLL_CHAN, NEVEN_PDF, NEVEN_PAIR_PDF, NEVEN_COLL_PDF

real(kind=8), dimension(:,:,:), allocatable :: XSTAT, DEL_X
real(kind=8), dimension(:,:,:,:,:), allocatable :: B_COLL, A_COLL

real(kind=8), dimension(:,:,:,:,:), allocatable :: MEAN_TIME_COLL_CHAN
real(kind=8), dimension(:,:,:,:), allocatable :: MEAN_TIME_PART_CHAN
real(kind=8), dimension(:,:,:,:,:), allocatable :: MEAN_TIME_PDF_CHAN
real(kind=8), dimension(:,:,:,:,:,:), allocatable :: MEAN_TIME_PDF_RELV_CHAN

real(kind=8), dimension(2) :: MEAN_PC_NORM
real(kind=8), dimension(2) :: MEAN_PC_X
real(kind=8), dimension(2) :: MEAN_PC_Y
real(kind=8), dimension(2) :: MEAN_PC_Z
integer :: NEVEN_L, NEVEN_R, NEVEN_WALL

! real(kind=8), dimension(:,:,:,:), allocatable :: PDFWALL
! real(kind=8), dimension(:,:,:,:), allocatable :: MEAN_PDFWALL
! real(kind=8), dimension(:,:), allocatable :: DPDFWALL

end module STATISTICS
