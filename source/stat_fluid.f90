!!====================================================================
!!
!! This fortran file compute all statistics on the fluid velocity
!! field.
!!
!!====================================================================

subroutine STAT_FLUID(NCYCLE,TIME)

!!====================================================================
!! We compute 4 levels of statistics, with graduate computational cost.
!!
!! LEVEL_STFLU = 1
!!     + Mean velocity
!!     + Reynolds stresses
!!     + Kolmogorov's scales
!!     + Taylor's scale (Tg)
!!
!! LEVEL_STFLU = 2
!!     + Divergence
!!     + Dissipation
!!
!! LEVEL_STFLU = 3
!!     + Flatness coefficient
!!     + Skewness coefficient
!!     + Fluid velocity PDF
!!     + Fluid velocity gradient PDF
!!     + Taylor's scales
!!
!! LEVEL_STFLU = 4
!!     + Autocorrelation function f & g
!!     + Integral scales Lf & Lg
!!
!!--------------------------------------------------------------------
!! Time-averaged statistics are performed using a macro array
!! called MEAN_FLUID_TIME. The averaging is done as
!!   --> MEAN_FLUID_TIME = MEAN_FLUID_TIME + MEAN_FLUID
!!
!!--------------------------------------------------------------------
!! MEAN_FLUID( 1): <uf>
!!             2 : <vf>
!!             3 : <wf>
!!             4 : <uf.uf>
!!             5 : <vf.vf>
!!             6 : <wf.wf>
!!             7 : <uf.vf>
!!             8 : <uf.wf>
!!             9 : <vf.wf>
!!            10 :  kf = (<uf.uf>+<vf.vf>+<wf.wf>)/2
!!            11 : epsilon
!!            12 : eta_K
!!            13 : tau_K
!!            14 : max(|dui/dxi|)
!!            15 :  < dui/dxi >
!!            16 : <(dui/dxi)^2>
!!            17 : v_K
!!            18 : Sk (skewness)
!!            19 : Tk (flatness)
!!            20 : <(dui/dxi)>
!!            21 : <(dui/dxi)^2>
!!            22 : <(dui/dxi)^3>
!!            23 : <(dui/dxi)^4>
!!            24 : <(dui/dxj)^2>
!!            25 : Taylor f
!!            26 : Taylor g
!!====================================================================

use dns_dim            !- Dimension
use param_phys         !- Physical & numerical parameters
use fluid_variable     !- Fluid velocity
use geometric_variable !- 
use statistics         !- Statistics
use check_cpu          !- CPU time checks

implicit none


!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Global arrays
!---------------------------------------------------------------------
!- cycle number
integer, intent(in) :: NCYCLE

!- Curent time
real(kind=8), intent(in) :: TIME


real(kind=8), dimension(DIMSCOR) :: RUXLOC, RVXLOC, RWXLOC
!---------------------------------------------------------------------
!- Local arrays
!---------------------------------------------------------------------
!- 
real(kind=8) :: RDUMMY

!- Taylor microscale
real(kind=8) :: LAMBDA_F, LAMBDA_G

!- Time control variable
real(kind=8) :: TIME_START, TIME_END


integer :: IFLAG1, I
!---------------------------------------------------------------------

!!- Check CPU time
if(MYID == 0) then
 TIME_START=MPI_WTIME()
end if


!!- Initiation
MEAN_FLUID(:) = 0.


!!====================================================================
!! 1. 1st level of  statistic
!!====================================================================
if(LEVEL_STFLU >= 1) then

!!--------------------------------------------------------------------
!! 1.1. Space-Averaging
!!--------------------------------------------------------------------
!!- <uf>
 call RSUMCPU(SUM(UFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))),MEAN_FLUID(1))

!!- <vf>
 call RSUMCPU(SUM(VFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))),MEAN_FLUID(2))

!!- <wf>
 call RSUMCPU(SUM(WFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))),MEAN_FLUID(3))

!!- <uf*uf>
 call RSUMCPU(SUM(UFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))**2),MEAN_FLUID(4))

!!- <vf*vf>
 call RSUMCPU(SUM(VFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))**2),MEAN_FLUID(5))

!!- <wf*wf>
 call RSUMCPU(SUM(WFLU(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))**2),MEAN_FLUID(6))

!!- <uf*vf>
 call RSUMCPU(SUM(UFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))&
                 *VFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))),MEAN_FLUID(7))

!!- <uf*wf>
 call RSUMCPU(SUM(UFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3)) &
                 *WFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))),MEAN_FLUID(8))

!!- <wf*vf>
 call RSUMCPU(SUM(VFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3)) &
                 *WFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))),MEAN_FLUID(9))

!!- qf = 0.5*(uf*uf + vf*vf + wf*wf)
 MEAN_FLUID(10) = 0.5*(MEAN_FLUID(4) + MEAN_FLUID(5) + MEAN_FLUID(6))

!!- Dissipation
 if(NCYCLE>1) call RSUMCPU(EPS_FLU,MEAN_FLUID(11))



!- Normalization by the full number of points
 MEAN_FLUID = MEAN_FLUID / NGLOB




 if(MEAN_FLUID(11) > ZERO) then
!!- Kolmogorov length scale
 MEAN_FLUID(12) = (VISC**3/MEAN_FLUID(11))**0.25

!!- Kolmogorov time scale
 MEAN_FLUID(13) = (VISC/MEAN_FLUID(11))**0.5

!!- Kolmogorov velocity scale
 MEAN_FLUID(17) = (VISC*MEAN_FLUID(11))**0.25

!!- Taylor g
 MEAN_FLUID(25) = (15.*VISC*2./3.*MEAN_FLUID(10)/MEAN_FLUID(11))**0.5

 else
 MEAN_FLUID(12) = 0.
 MEAN_FLUID(13) = 0.
 MEAN_FLUID(17) = 0.
 MEAN_FLUID(25) = 0.
 end if
!!--------------------------------------------------------------------
!! 1.2. Print in file
!!--------------------------------------------------------------------
 if(MYID==0) write(301,10000)TIME,       &
                   MEAN_FLUID(1 ), & !- <uf>
                   MEAN_FLUID(2 ), & !- <vf>
                   MEAN_FLUID(3 ), & !- <wf>
                   MEAN_FLUID(10), & !-  kf
                   MEAN_FLUID(11), & !-  epsilon_f
                   MEAN_FLUID(4 ), & !- <uf.uf>
                   MEAN_FLUID(5 ), & !- <vf.vf>
                   MEAN_FLUID(6 ), & !- <wf.wf>
                   MEAN_FLUID(7 ), & !- <uf.vf>
                   MEAN_FLUID(8 ), & !- <uf.wf>
                   MEAN_FLUID(9 ), & !- <vf.wf>
                   MEAN_FLUID(12), & !- eta_K
                   MEAN_FLUID(13), & !- tau_K
                   MEAN_FLUID(17), & !- v_K
                   MEAN_FLUID(25)    !- Taylor g

end if


!!====================================================================
!! 2. Divergence and so on
!!====================================================================
if(LEVEL_STFLU >= 2) then

 call DIVERGENCE(MEAN_FLUID(14),MEAN_FLUID(15),MEAN_FLUID(16))

!!--------------------------------------------------------------------
!! 2.2. Print in file
!!--------------------------------------------------------------------
 if(MYID==0) write(302,10000)TIME, &
                   MEAN_FLUID(14), & !- max(|dui/dxi|)
                   MEAN_FLUID(15), & !- < dui/dxi >
                   MEAN_FLUID(16)    !-<(dui/dxi)^2>
end if



!!======================================================================
!! . Skewness, Flatness, PDF of fluid velocity & gradients
!!======================================================================
!! These statistics are time consuming because of the large number
!! of FFT required for the gradient computation
!!======================================================================
if(LEVEL_STFLU >= 3) then

!!----------------------------------------------------------------------
!! 3.1.  Fluid gradient
!!---------------------------------------------------------------------
!!-  using FFT
 call FLUID_GRADIENT

!!-  using finite difference
!!  call FLUID_GRADIENT_FD

!!----------------------------------------------------------------------
!! 3.2. CPU-averaging
!!---------------------------------------------------------------------
!!- <(dui/dxi)^2>
 RDUMMY = (DUIDXJ(1,1,2) + DUIDXJ(2,2,2) + DUIDXJ(3,3,2))/3.
 call RSUMCPU(RDUMMY,MEAN_FLUID(21))

!!- <(dui/dxi)^3>
 RDUMMY = (DUIDXJ(1,1,3) + DUIDXJ(2,2,3) + DUIDXJ(3,3,3))/3.
 call RSUMCPU(RDUMMY,MEAN_FLUID(22))

!!- <(dui/dxi)^4>
 RDUMMY = (DUIDXJ(1,1,4) + DUIDXJ(2,2,4) + DUIDXJ(3,3,4))/3.
 call RSUMCPU(RDUMMY,MEAN_FLUID(23))

!!- <(dui/dxj)^2>
 RDUMMY = (   DUIDXJ(1,2,2) + DUIDXJ(1,3,2) &
            + DUIDXJ(2,1,2) + DUIDXJ(2,3,2) &
            + DUIDXJ(3,1,2) + DUIDXJ(3,2,2) )/6.
 call RSUMCPU(RDUMMY,MEAN_FLUID(24))


!!- Skewness
 MEAN_FLUID(18) = MEAN_FLUID(22)/ MEAN_FLUID(21)**1.5

!!- Flatness
 MEAN_FLUID(19) = MEAN_FLUID(23)/ MEAN_FLUID(21)**2

!!- Taylor f
 MEAN_FLUID(24) = (4./3.*MEAN_FLUID(10)/MEAN_FLUID(21))**0.5



!!----------------------------------------------------------------------
!! 3.3. Print in file
!!---------------------------------------------------------------------
 if(MYID==0) write(303,10000)TIME, &
                   MEAN_FLUID(18), & !- Skewness
                   MEAN_FLUID(19), & !- Flatness
                   MEAN_FLUID(24), & !- Taylor f
                   MEAN_FLUID(25), & !- Taylor g
                   MEAN_FLUID(24)/MEAN_FLUID(25)/sqrt(2.)
end if


!!======================================================================
!! 4. Two points statistics (Spatial correlation)
!!======================================================================
if(LEVEL_STFLU >= 4.and.STAT_TIME) then
!
!
 call CORRTWOPTS(ISIZE(1),  &
                   ISIZE(2),  &
		   ISIZE(3),  &
		   UFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),&
		   MEAN_FLUID(1),&
		   RUXLOC)
		   
 call CORRTWOPTS(ISIZE(1),  &
                   ISIZE(2),  &
		   ISIZE(3),  &
		   VFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),  &
		   MEAN_FLUID(2),&
		   RVXLOC)
 
 call CORRTWOPTS(ISIZE(1),  &
                   ISIZE(2),  &
		   ISIZE(3),  &
		   WFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),  &
		   MEAN_FLUID(3),  &
		   RWXLOC)
 !
 MEAN_RUXLOC(:) = MEAN_RUXLOC(:)  + RUXLOC(:)
 MEAN_RVXLOC(:) = MEAN_RVXLOC(:)  + RVXLOC(:)
 MEAN_RWXLOC(:) = MEAN_RWXLOC(:)  + RWXLOC(:)


end if





!!======================================================================
!! 5. Time averaging if so
!!======================================================================
if(STAT_TIME) then
 MEAN_TIME_FLUID = MEAN_TIME_FLUID + MEAN_FLUID
end if 

!!======================================================================
!! 6. Print statistics in "info" files
!!======================================================================

!- Print in file "stat.info"
if((mod(NCYCLE,FOUT0) == 0).and.(MYID==0)) then
 IFLAG1 = 2
 write(UNIT_INFO(IFLAG1),*)
 write(UNIT_INFO(IFLAG1),*)'====================================='
 write(UNIT_INFO(IFLAG1),*)'Statistic for cycle = ',NCYCLE
 write(UNIT_INFO(IFLAG1),*)'====================================='
 write(UNIT_INFO(IFLAG1),10601)'             <uf>= ',MEAN_FLUID(1)
 write(UNIT_INFO(IFLAG1),10601)'             <vf>= ',MEAN_FLUID(2)
 write(UNIT_INFO(IFLAG1),10601)'             <wf>= ',MEAN_FLUID(3)   
 write(UNIT_INFO(IFLAG1),10601)'--'
 write(UNIT_INFO(IFLAG1),10601)'  <uf.uf>/(2kf/3)= ',MEAN_FLUID(4)/MEAN_FLUID(10)
 write(UNIT_INFO(IFLAG1),10601)'  <vf.vf>/(2kf/3)= ',MEAN_FLUID(5)/MEAN_FLUID(10)
 write(UNIT_INFO(IFLAG1),10601)'  <wf.wf>/(2kf/3)= ',MEAN_FLUID(6)/MEAN_FLUID(10)   
 write(UNIT_INFO(IFLAG1),*)
 write(UNIT_INFO(IFLAG1),10601)'  <uf.vf>/(2kf/3)= ',MEAN_FLUID(7)/MEAN_FLUID(10)
 write(UNIT_INFO(IFLAG1),10601)'  <uf.wf>/(2kf/3)= ',MEAN_FLUID(8)/MEAN_FLUID(10)
 write(UNIT_INFO(IFLAG1),10601)'  <vf.wf>/(2kf/3)= ',MEAN_FLUID(9)/MEAN_FLUID(10)   
 write(UNIT_INFO(IFLAG1),*)
 write(UNIT_INFO(IFLAG1),10601)'               kf=',MEAN_FLUID(10)
 write(UNIT_INFO(IFLAG1),10601)'             epsf=',MEAN_FLUID(11)
 if(LEVEL_STFLU>=2) then
 write(UNIT_INFO(IFLAG1),10601)'--'
 write(UNIT_INFO(IFLAG1),10601)'    max(|div(U)|)= ',MEAN_FLUID(14)
 write(UNIT_INFO(IFLAG1),10601)'        < div(U)>= ',MEAN_FLUID(15)
 write(UNIT_INFO(IFLAG1),10601)'       <div(U)^2>= ',MEAN_FLUID(16)
 end if
 if(LEVEL_STFLU>=3) then
 write(UNIT_INFO(IFLAG1),10601)'             epsf=',15.*VISC*MEAN_FLUID(21)
 write(UNIT_INFO(IFLAG1),10601)'--'
 write(UNIT_INFO(IFLAG1),10601)'                  Sk=',MEAN_FLUID(18)
 write(UNIT_INFO(IFLAG1),10601)'                  Tk=',MEAN_FLUID(19)
 write(UNIT_INFO(IFLAG1),10601)'    <(duf_i/dx_i)^2>=',MEAN_FLUID(21)
 write(UNIT_INFO(IFLAG1),10601)'Taylor            lf=',MEAN_FLUID(24)
 write(UNIT_INFO(IFLAG1),10601)'Taylor            lg=',MEAN_FLUID(25)
 write(UNIT_INFO(IFLAG1),10601)'Taylor lf/(lg*2^0.5)=',MEAN_FLUID(24)/MEAN_FLUID(25)/sqrt(2.)
 end if
 write(UNIT_INFO(IFLAG1),10601)'---------------------------------'
end if





!!- CPU check
if(MYID == 0) then
 TIME_END=MPI_WTIME()
 CPU_FLUID(6) = CPU_FLUID(6) + TIME_END - TIME_START
end if


!!----------------------------------------------------------------------
10000 format (20(e17.7))
10601 format (2x,A,E13.6)
10602 format (2x,A,E13.6,1x,A,E13.6)
10603 format (2x,A,E13.6,2x,A,E13.6,2x,A,E13.6)

end subroutine STAT_FLUID
