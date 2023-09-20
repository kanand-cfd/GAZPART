!!====================================================================
!!
!! This fortran file compute all statistics on the paticles
!!
!!====================================================================

subroutine STAT_PARTICLE(NCYCLE,TIME)

!!====================================================================
!! We compute 4 levels of statistics, with graduate computational cost.
!!
!! LEVEL1_STPAR: + Mean velocity
!!               + Kinetic stress
!!
!! LEVEL2_STPAR: + Lagrangian autocorrelation functions
!!
!!
!! LEVEL3_STPAR: + Two-particle correlation
!!
!!
!!--------------------------------------------------------------------
!! Time-averaged statistics are performed using a macro array
!! called MEAN_TIME_PART. The averaging is done as
!!   --> MEAN_TIME_PART = MEAN_TIME_PART + MEAN_PART
!!
!!--------------------------------------------------------------------
!!  MEAN_PART( 1): <up>
!!             2 : <vp>
!!             3 : <wp>
!!             4 : <up.up>
!!             5 : <vp.vp>
!!             6 : <wp.wp>
!!             7 : <up.vp>
!!             8 : <up.wp>
!!             9 : <vp.wp>
!!            10 :  qp = (<up.up>+<vp.vp>+<wp.wp>)/2
!!            11 : <uf@p>
!!            12 : <vf@p>
!!            13 : <wf@p>
!!            14 : <uf@p.uf@p>
!!            15 : <vf@p.vf@p>
!!            16 : <wf@p.wf@p>
!!            17 : <uf@p.vf@p>
!!            18 : <uf@p.wf@p>
!!            19 : <vf@p.wf@p>
!!            20 :  qf@p = (<uf@p.uf@p>+<vf@p.vf@p>+<wf@p.wf@p>)/2
!!            21 : <up.uf@p>
!!            22 : <vp.vf@p>
!!            23 : <wp.wf@p>
!!            24 : <up.vf@p>
!!            25 : <up.wf@p>
!!            26 : <vp.wf@p>
!!            27 :  qfp = <up.uf@p>+<vp.vf@p>+<wp.wf@p>
!!            28 : <1/tfpF>
!!            29 : <Vr>
!!            30 : <Cd>
!!            31 : <Rep>
!!            32 : Np_class
!!            33 : <dp>
!!            34 : <dp^2>
!!            35 : <dp^3>
!!
!!            59 : 
!!            60 : 
!!            61 : <np>
!!            62 : <np'^2>
!!            63 : <np^2>
!!            64 : 
!!            65 : 
!!            66 : 
!!            67 : 
!!
!!====================================================================

use DNS_DIM	       !- Dimension
use PARAM_PHYS         !- Physical & numerical parameters
use PARTICLE_PARALLEL
use GEOMETRIC_VARIABLE !- 
use STATISTICS         !- Statistics
use CHECK_CPU	       !- CPU time checks

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

!---------------------------------------------------------------------
!- Local arrays
!---------------------------------------------------------------------
!- Arrays containing statistics
real(kind=8), dimension(NSTAT,NIG,POLYDISPMAX) :: MEAN_PART_LOC
real(kind=8), dimension(NSTAT,NIG,POLYDISPMAX) :: MEAN_PART
real(kind=8), dimension(NSTAT, NIG, POLYDISPMAX, NPDF) :: MEAN_PART_PDF



!- Number of particle per class
real(kind=8), dimension(POLYDISPMAX) :: NPM

!- Mean velocities
real(kind=8), dimension(POLYDISPMAX) :: UPM, VPM, WPM
real(kind=8), dimension(POLYDISPMAX) :: OXM, OYM, OZM, GOXM, GOYM, GOZM

!- Fluctuating particle velocities
real(kind=8) :: UPFLC, VPFLC, WPFLC

!- Fluctuating particle angular velocities
real(kind=8) :: OXFLC, OYFLC, OZFLC, GOXFLC, GOYFLC, GOZFLC

real(kind=8), dimension(ndim, 1) :: global_angular_velocity, principal_axis !, omega_fluct, Amom

real(kind=8), dimension(ndim, ndim) :: I_p


real(kind=8) :: MIN_PDF, MAX_PDF, DPDF, MIN_PX, MAX_PX, DPX
real(kind=8) :: MIN_OMEGA_PDF, MAX_OMEGA_PDF, DPDF_OMEGA


!- file number
integer :: NUMFILE


!- Particle concentration
real(kind=8) :: DXCP, DYCP,DZCP, CMEAN


!- 
integer :: NTLGR 

!- 
real(kind=8) :: RDUMMY

integer :: IFLAG

!- Index of particle class for polydisperse stat
integer :: IDP


!- Time control variable
real(kind=8) :: TIME_START, TIME_END

integer :: I,J,K, IFLAG1, IPART, JPART, KPART, IPDF, NLGR, NX0
!---------------------------------------------------------------------

!!- CPU check
if(MYID == 0) then

  TIME_START=MPI_WTIME()

end if


MEAN_PART(:,:,:) = ZERO
MEAN_PART_LOC(:,:,:) = ZERO
MEAN_PART_PDF(:,:,:,:) = ZERO

!!====================================================================
!! 1. 1st level of  statistic
!!====================================================================
if(LEVEL1_STPAR) then

  do J = 1, NIG

    !!--------------------------------------------------------------------
    !! 1.1. Mean variable
    !!--------------------------------------------------------------------
    do I = 1, NPART_LOC(J)

      !!- index of PSD
      IDP = PART(I,J)%IDP

      !!--------------------------------------------------------------------
      !!- Particle number for each class
      !!--------------------------------------------------------------------
      !!- Np
      MEAN_PART_LOC(23,J,IDP) = MEAN_PART_LOC(23,J,IDP) + 1.0


      !!--------------------------------------------------------------------
      !!- Particle velocities
      !!--------------------------------------------------------------------
      !!- <up>
      MEAN_PART_LOC(1,J,IDP) = MEAN_PART_LOC(1,J,IDP) + PART(I,J)%UP

      !!- <vp>
      MEAN_PART_LOC(2,J,IDP) = MEAN_PART_LOC(2,J,IDP) + PART(I,J)%VP

      !!- <wp>
      MEAN_PART_LOC(3,J,IDP) = MEAN_PART_LOC(3,J,IDP) + PART(I,J)%WP


      !!--------------------------------------------------------------------
      !!- Particle Rotational velocities (particle frame)
      !!--------------------------------------------------------------------
      !!- <omega_x,p>
      MEAN_PART_LOC(11,J,IDP) = MEAN_PART_LOC(11,J,IDP) + PART(I,J)%OMEGAX

      !!- <omega_y,p>
      MEAN_PART_LOC(12,J,IDP) = MEAN_PART_LOC(12,J,IDP) + PART(I,J)%OMEGAY      

      !!- <omega_z,p>
      MEAN_PART_LOC(13,J,IDP) = MEAN_PART_LOC(13,J,IDP) + PART(I,J)%OMEGAZ


      !!--------------------------------------------------------------------
      !!- Particle Rotational velocities (comoving frame)
      !!--------------------------------------------------------------------
      global_angular_velocity(1, 1) = PART(I,J)%OMEGAX
      global_angular_velocity(2, 1) = PART(I,J)%OMEGAY
      global_angular_velocity(3, 1) = PART(I,J)%OMEGAZ

      call transform_basis(global_angular_velocity, PART(I,J)%ELLQUAT, shape(global_angular_velocity))

      MEAN_PART_LOC(14, J, IDP) = MEAN_PART_LOC(14, J, IDP) + global_angular_velocity(1,1)
      MEAN_PART_LOC(15, J, IDP) = MEAN_PART_LOC(15, J, IDP) + global_angular_velocity(2,1)
      MEAN_PART_LOC(16, J, IDP) = MEAN_PART_LOC(16, J, IDP) + global_angular_velocity(3,1)


    end do


    do IDP = 1, POLYDISP(J)


      call RSUMCPU(MEAN_PART_LOC(23,J,IDP),NPM(IDP))

      !!- Compute the mean velocities over the whole domain
      call RSUMCPU(MEAN_PART_LOC( 1,J,IDP),UPM(IDP))
      call RSUMCPU(MEAN_PART_LOC( 2,J,IDP),VPM(IDP))
      call RSUMCPU(MEAN_PART_LOC( 3,J,IDP),WPM(IDP))
  

      UPM(IDP) = UPM(IDP) / NPM(IDP)
      VPM(IDP) = VPM(IDP) / NPM(IDP)
      WPM(IDP) = WPM(IDP) / NPM(IDP)


      call RSUMCPU(MEAN_PART_LOC(11, J, IDP), OXM(IDP))
      call RSUMCPU(MEAN_PART_LOC(12, J, IDP), OYM(IDP))
      call RSUMCPU(MEAN_PART_LOC(13, J, IDP), OZM(IDP))

      OXM(IDP) = OXM(IDP) / NPM(IDP)
      OYM(IDP) = OYM(IDP) / NPM(IDP)
      OZM(IDP) = OZM(IDP) / NPM(IDP)

      call RSUMCPU(MEAN_PART_LOC(14, J, IDP), GOXM(IDP))
      call RSUMCPU(MEAN_PART_LOC(15, J, IDP), GOYM(IDP))
      call RSUMCPU(MEAN_PART_LOC(16, J, IDP), GOZM(IDP))

      GOXM(IDP) = GOXM(IDP) / NPM(IDP)
      GOYM(IDP) = GOYM(IDP) / NPM(IDP)
      GOZM(IDP) = GOZM(IDP) / NPM(IDP)


    end do

    !!--------------------------------------------------------------------
    !! 1.2. Particle kinetic stress
    !!--------------------------------------------------------------------

    do I = 1, NPART_LOC(J)

      !!- index of PSD
      IDP = PART(I,J)%IDP

      !!- u'p = up-<up>
      UPFLC = PART(I,J)%UP - UPM(IDP)

      !!- v'p = vp-<vp>
      VPFLC = PART(I,J)%VP - VPM(IDP)

      !!- w'p = wp-<wp>
      WPFLC = PART(I,J)%WP - WPM(IDP)


      !!- <up*up>
      MEAN_PART_LOC(4,J,IDP) = MEAN_PART_LOC(4,J,IDP) + UPFLC*UPFLC

      !!- <vp*vp>
      MEAN_PART_LOC(5,J,IDP) = MEAN_PART_LOC(5,J,IDP) + VPFLC*VPFLC

      !!- <wp*wp>
      MEAN_PART_LOC(6,J,IDP) = MEAN_PART_LOC(6,J,IDP) + WPFLC*WPFLC

      !!- <up*vp>
      MEAN_PART_LOC(7,J,IDP) = MEAN_PART_LOC(7,J,IDP) + UPFLC*VPFLC

      !!- <up*wp>
      MEAN_PART_LOC(8,J,IDP) = MEAN_PART_LOC(8,J,IDP) + UPFLC*WPFLC

      !!- <vp*wp>
      MEAN_PART_LOC(9,J,IDP) = MEAN_PART_LOC(9,J,IDP) + VPFLC*WPFLC

    end do

    !!- qp - Translational
    do IDP = 1, POLYDISP(J)

      MEAN_PART_LOC(10,J,IDP) = 0.5*(MEAN_PART_LOC(4,J,IDP) &
                                   + MEAN_PART_LOC(5,J,IDP) &
                                   + MEAN_PART_LOC(6,J,IDP) )

    end do

    !!--------------------------------------------------------------------
    !! 1.2. Particle Rotational kinetic stress
    !!--------------------------------------------------------------------

    do I = 1, NPART_LOC(J)

      !!- index of PSD
      IDP = PART(I,J)%IDP

      !!---------------------- Particle Frame --------------------!!

      !!- omega_x' = omega_x - <omega_x>
      OXFLC = PART(I,J)%OMEGAX - OXM(IDP)

      !!- omega_x' = omega_x - <omega_x>
      OYFLC = PART(I,J)%OMEGAY - OYM(IDP)

      !!- omega_x' = omega_x - <omega_x>
      OZFLC = PART(I,J)%OMEGAZ - OZM(IDP)


      !!-<omega_x*omega_x>
      MEAN_PART_LOC(17, J, IDP) = MEAN_PART_LOC(17, J, IDP) + OXFLC*OXFLC
      !!-<omega_y*omega_y>
      MEAN_PART_LOC(18, J, IDP) = MEAN_PART_LOC(18, J, IDP) + OYFLC*OYFLC
      !!-<omega_z*omega_z>
      MEAN_PART_LOC(19, J, IDP) = MEAN_PART_LOC(19, J, IDP) + OZFLC*OZFLC


      !!-<omega_x*omega_y>
      MEAN_PART_LOC(20, J, IDP) = MEAN_PART_LOC(20, J, IDP) + OXFLC*OYFLC
      !!-<omega_x*omega_z>
      MEAN_PART_LOC(21, J, IDP) = MEAN_PART_LOC(21, J, IDP) + OXFLC*OZFLC
      !!-<omega_y*omega_z>
      MEAN_PART_LOC(22, J, IDP) = MEAN_PART_LOC(22, J, IDP) + OYFLC*OZFLC      


      !!-<I_p,xx*omega_x*omega_x>
      MEAN_PART_LOC(24, J, IDP) = MEAN_PART_LOC(24, J, IDP) + IPXX(J, IDP)*OXFLC*OXFLC
      !!-<I_p,yy*omega_y*omega_y>
      MEAN_PART_LOC(25, J, IDP) = MEAN_PART_LOC(25, J, IDP) + IPYY(J, IDP)*OYFLC*OYFLC
      !!-<I_p,zz*omega_z*omega_z>
      MEAN_PART_LOC(26, J, IDP) = MEAN_PART_LOC(26, J, IDP) + IPZZ(J, IDP)*OZFLC*OZFLC


      !!--------------------- Comoving Frame -------------------!!
      ! Initialization of Inertia Matrix
      I_p(:,:) = 0.0
      I_p(1,1) = IPXX(J, IDP)
      I_p(2,2) = IPYY(J, IDP)
      I_p(3,3) = IPZZ(J, IDP)

      global_angular_velocity(1, 1) = PART(I,J)%OMEGAX
      global_angular_velocity(2, 1) = PART(I,J)%OMEGAY
      global_angular_velocity(3, 1) = PART(I,J)%OMEGAZ

      call transform_basis(global_angular_velocity, PART(I,J)%ELLQUAT, shape(global_angular_velocity))

      call transform_basis(I_p, PART(I,J)%ELLQUAT, shape(I_p))

      !! !! !! !! !! !! !! !! !! !! !! !! !! !! !! !!

      GOXFLC = global_angular_velocity(1,1) - GOXM(IDP)
      GOYFLC = global_angular_velocity(2,1) - GOYM(IDP)
      GOZFLC = global_angular_velocity(3,1) - GOZM(IDP)


      MEAN_PART_LOC(27, J, IDP) = MEAN_PART_LOC(27, J, IDP) + GOXFLC*GOXFLC

      MEAN_PART_LOC(28, J, IDP) = MEAN_PART_LOC(28, J, IDP) + GOYFLC*GOYFLC    

      MEAN_PART_LOC(29, J, IDP) = MEAN_PART_LOC(29, J, IDP) + GOZFLC*GOZFLC


      MEAN_PART_LOC(30, J, IDP) = MEAN_PART_LOC(30, J, IDP) + GOXFLC*GOYFLC

      MEAN_PART_LOC(31, J, IDP) = MEAN_PART_LOC(31, J, IDP) + GOXFLC*GOZFLC

      MEAN_PART_LOC(32, J, IDP) = MEAN_PART_LOC(32, J, IDP) + GOYFLC*GOZFLC


      MEAN_PART_LOC(33, J, IDP) = MEAN_PART_LOC(33, J, IDP) + I_p(1,1)*GOXFLC*GOXFLC

      MEAN_PART_LOC(34, J, IDP) = MEAN_PART_LOC(34, J, IDP) + I_p(2,2)*GOYFLC*GOYFLC

      MEAN_PART_LOC(35, J, IDP) = MEAN_PART_LOC(35, J, IDP) + I_p(3,3)*GOZFLC*GOZFLC


      MEAN_PART_LOC(36, J, IDP) = MEAN_PART_LOC(36, J, IDP) + I_p(1,2)*GOXFLC*GOYFLC

      MEAN_PART_LOC(37, J, IDP) = MEAN_PART_LOC(37, J, IDP) + I_p(1,3)*GOXFLC*GOZFLC

      MEAN_PART_LOC(38, J, IDP) = MEAN_PART_LOC(38, J, IDP) + I_p(2,3)*GOYFLC*GOZFLC 


    end do


    !!-qp - Rotational
    do IDP = 1, POLYDISP(J)

      MEAN_PART_LOC(39, J, IDP) = 0.5*(IPXX(J, IDP)*MEAN_PART_LOC(17, J, IDP) &
                                     + IPYY(J, IDP)*MEAN_PART_LOC(18, J, IDP) &
                                     + IPZZ(J, IDP)*MEAN_PART_LOC(19, J, IDP) )


      MEAN_PART_LOC(40, J, IDP) = 0.5*(MEAN_PART_LOC(24, J, IDP) &
                                     + MEAN_PART_LOC(25, J, IDP) &
                                     + MEAN_PART_LOC(26, J, IDP) )

      !!- Total
      MEAN_PART_LOC(41, J, IDP) = MEAN_PART_LOC(10, J, IDP) + MEAN_PART_LOC(40, J, IDP)

    end do  

 

  end do!!- End loop over NIG

  !!====================================================================
  !! Summation overall domain and normalization
  !!====================================================================
  !if (NPROC>1) then
  do J = 1, NIG

    do IDP = 1, POLYDISP(J)

      call MPI_ALLREDUCE(MEAN_PART_LOC(:,J,IDP),MEAN_PART(:,J,IDP),NSTAT,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

      !!- Normalization by the full number of particle
      MEAN_PART(:,J,IDP) = MEAN_PART(:,J,IDP) / NPM(IDP)

      !!- Rewrite the number of particles
      MEAN_PART(23,J,IDP) = NPM(IDP)

    end do

  end do

!!--------------------------------------------------------------------
!! 1.4. Print in file
!!--------------------------------------------------------------------
if(MYID==0)  then

  do J = 1, NIG

    do IDP = 1, POLYDISP(J)

      NUMFILE = 500+(J-1)*POLYDISP(J)+IDP

      write(NUMFILE,10000)   &
        TIME, MEAN_PART( 1,J,IDP),  & !- <up>
              MEAN_PART( 2,J,IDP),  & !- <vp>
              MEAN_PART( 3,J,IDP),  & !- <wp>
  
              MEAN_PART( 4,J,IDP),  & !- <up.up>
              MEAN_PART( 5,J,IDP),  & !- <vp.vp>
              MEAN_PART( 6,J,IDP),  & !- <wp.wp>

              MEAN_PART( 7,J,IDP),  & !- <up.vp>
              MEAN_PART( 8,J,IDP),  & !- <up.wp>
              MEAN_PART( 9,J,IDP),  & !- <vp.wp>

              MEAN_PART(10,J,IDP)     !- qp - Translational


      NUMFILE = 520+(J-1)*POLYDISP(J)+IDP


      write(NUMFILE,10000) &
        TIME, MEAN_PART(11,J,IDP),  & !- <omega_x>
              MEAN_PART(12,J,IDP),  & !- <omega_y>
              MEAN_PART(13,J,IDP),  & !- <omega_z>

              MEAN_PART(14,J,IDP),  & !- global <omega_x>
              MEAN_PART(15,J,IDP),  & !- global <omega_y>
              MEAN_PART(16,J,IDP),  & !- global <omega_z>

              !!Particle Frame !!
              MEAN_PART(17,J,IDP),  & !! <oxp.oxp>
              MEAN_PART(18,J,IDP),  & !! <oyp.oyp>
              MEAN_PART(19,J,IDP),  & !! <ozp.ozp>

              MEAN_PART(20,J,IDP),  & !! <oxp.oyp>
              MEAN_PART(21,J,IDP),  & !! <oxp.ozp>
              MEAN_PART(22,J,IDP),  & !! <oyp.ozp>

              MEAN_PART(24,J,IDP),  & !! Ixx<oxp.oxp>
              MEAN_PART(25,J,IDP),  & !! Iyy<oyp.oyp>
              MEAN_PART(26,J,IDP),  & !! Izz<ozp.ozp>

              !!Comoving Frame !!
              MEAN_PART(27,J,IDP),  & !! Global-<oxp.oxp>
              MEAN_PART(28,J,IDP),  & !! Global-<oyp.oyp>
              MEAN_PART(29,J,IDP),  & !! Global-<ozp.ozp>

              MEAN_PART(30,J,IDP),  & !! Global-<oxp.oyp>
              MEAN_PART(31,J,IDP),  & !! Global-<oxp.ozp>
              MEAN_PART(32,J,IDP),  & !! Global-<oyp.ozp>

              MEAN_PART(33,J,IDP),  & !! Global-<Ixx.oxp.oxp>
              MEAN_PART(34,J,IDP),  & !! Global-<Iyy.oyp.oyp>
              MEAN_PART(35,J,IDP),  & !! Global-<Izz.ozp.ozp>

              MEAN_PART(36,J,IDP),  & !! Global-<Ixy.oxp.oyp>
              MEAN_PART(37,J,IDP),  & !! Global-<Ixz.oxp.ozp>
              MEAN_PART(38,J,IDP),  & !! Global-<Iyz.oyp.ozp>

              !!Energy!!
              MEAN_PART(39,J,IDP),  & !- qp - Rotational - 1
              MEAN_PART(40,J,IDP),  & !- qp - Rotational - 2
              MEAN_PART(41,J,IDP),  & !- qp - Total
                
              MEAN_PART(23,J,IDP)     !- <Np>

    end do !!- End loop over POLYDISP

  end do !!- End loop over NIG

end if !- End uf MYID==0

end if

!!======================================================================
!! Calculation of PDF
!!======================================================================
if((mod(NCYCLE, 500)==0) .and. LEVEL1_STPDF) then

  MAX_PDF =  2.0
  MIN_PDF = -2.0
  DPDF = (MAX_PDF - MIN_PDF)/real(NPDF)

  MIN_PX = - 1.0
  MAX_PX =   1.0
  DPX = (MAX_PX - MIN_PX)/real(NPDF)


  MAX_OMEGA_PDF = 2000
  MIN_OMEGA_PDF = -2000
  DPDF_OMEGA = (MAX_OMEGA_PDF - MIN_OMEGA_PDF)/real(NPDF)



  do J = 1, NIG

    do I = 1, NPART_LOC(J)


      IDP = PART(I,J)%IDP

      ! Velocity PDFs
      IPDF = int((PART(I,J)%UP - MIN_PDF)/DPDF) + 1
      if(IPDF < 1) IPDF = 1
      if(IPDF > NPDF) IPDF = NPDF
      MEAN_PART_PDF(1, J, IDP, IPDF) =  MEAN_PART_PDF(1, J, IDP, IPDF) + 1.0

      IPDF = int((PART(I,J)%VP - MIN_PDF)/DPDF) + 1
      if(IPDF < 1) IPDF = 1
      if(IPDF > NPDF) IPDF = NPDF
      MEAN_PART_PDF(2, J, IDP, IPDF) =  MEAN_PART_PDF(2, J, IDP, IPDF) + 1.0    

      IPDF = int((PART(I,J)%WP - MIN_PDF)/DPDF) + 1
      if(IPDF < 1) IPDF = 1
      if(IPDF > NPDF) IPDF = NPDF
      MEAN_PART_PDF(3, J, IDP, IPDF) =  MEAN_PART_PDF(3, J, IDP, IPDF) + 1.0


      !! Angular Velocity PDFs
      IPDF = int((PART(I,J)%OMEGAX - MIN_OMEGA_PDF)/DPDF_OMEGA) + 1
      if(IPDF < 1) IPDF = 1
      if(IPDF > NPDF) IPDF = NPDF
      MEAN_PART_PDF(4, J, IDP, IPDF) =  MEAN_PART_PDF(4, J, IDP, IPDF) + 1.0

      IPDF = int((PART(I,J)%OMEGAY - MIN_OMEGA_PDF)/DPDF_OMEGA) + 1
      if(IPDF < 1) IPDF = 1
      if(IPDF > NPDF) IPDF = NPDF
      MEAN_PART_PDF(5, J, IDP, IPDF) =  MEAN_PART_PDF(5, J, IDP, IPDF) + 1.0

      IPDF = int((PART(I,J)%OMEGAZ - MIN_OMEGA_PDF)/DPDF_OMEGA) + 1
      if(IPDF < 1) IPDF = 1
      if(IPDF > NPDF) IPDF = NPDF
      MEAN_PART_PDF(6, J, IDP, IPDF) =  MEAN_PART_PDF(6, J, IDP, IPDF) + 1.0

      !! Global angular velocity pdf
      global_angular_velocity(1, 1) = PART(I,J)%OMEGAX
      global_angular_velocity(2, 1) = PART(I,J)%OMEGAY
      global_angular_velocity(3, 1) = PART(I,J)%OMEGAZ

      call transform_basis(global_angular_velocity, PART(I,J)%ELLQUAT, shape(global_angular_velocity))

      IPDF = int((global_angular_velocity(1, 1) - MIN_OMEGA_PDF)/DPDF_OMEGA) + 1
      if(IPDF < 1) IPDF = 1
      if(IPDF > NPDF) IPDF = NPDF
      MEAN_PART_PDF(7, J, IDP, IPDF) =  MEAN_PART_PDF(7, J, IDP, IPDF) + 1.0

      IPDF = int((global_angular_velocity(2, 1) - MIN_OMEGA_PDF)/DPDF_OMEGA) + 1
      if(IPDF < 1) IPDF = 1
      if(IPDF > NPDF) IPDF = NPDF
      MEAN_PART_PDF(8, J, IDP, IPDF) =  MEAN_PART_PDF(8, J, IDP, IPDF) + 1.0

      IPDF = int((global_angular_velocity(3, 1) - MIN_OMEGA_PDF)/DPDF_OMEGA) + 1
      if(IPDF < 1) IPDF = 1
      if(IPDF > NPDF) IPDF = NPDF
      MEAN_PART_PDF(9, J, IDP, IPDF) =  MEAN_PART_PDF(9, J, IDP, IPDF) + 1.0

      !! Orientation
      principal_axis(1,1) = 1.0
      principal_axis(2,1) = 0.0
      principal_axis(3,1) = 0.0

      call transform_basis(principal_axis, PART(I,J)%ELLQUAT, shape(principal_axis))

      IPDF = int((principal_axis(1,1) - MIN_PX)/DPX) + 1
      if(IPDF < 1) IPDF = 1
      if(IPDF > NPDF) IPDF = NPDF
      MEAN_PART_PDF(10, J, IDP, IPDF) =  MEAN_PART_PDF(10, J, IDP, IPDF) + 1.0

      IPDF = int((principal_axis(2,1) - MIN_PX)/DPX) + 1
      if(IPDF < 1) IPDF = 1
      if(IPDF > NPDF) IPDF = NPDF
      MEAN_PART_PDF(11, J, IDP, IPDF) =  MEAN_PART_PDF(11, J, IDP, IPDF) + 1.0

      IPDF = int((principal_axis(3,1) - MIN_PX)/DPX) + 1
      if(IPDF < 1) IPDF = 1
      if(IPDF > NPDF) IPDF = NPDF
      MEAN_PART_PDF(12, J, IDP, IPDF) =  MEAN_PART_PDF(12, J, IDP, IPDF) + 1.0


    end do

  end do


  MEAN_TIME_PART_PDF = MEAN_TIME_PART_PDF + MEAN_PART_PDF
  NEVEN_PDF = NEVEN_PDF + 1

end if


if((NCYCLE == NCYCLEMAX) .and. LEVEL1_STPDF) then


  MEAN_TIME_PART_PDF = MEAN_TIME_PART_PDF/real(NEVEN_PDF)

  do J = 1, NIG

      do IDP = 1, POLYDISP(J)

        NUMFILE = 600+(J-1)*POLYDISP(J)+IDP

        do IPDF = 1, NPDF

          write(NUMFILE,10000) real(IPDF), &
          MEAN_TIME_PART_PDF(1, J, IDP, IPDF), & !- up
          MEAN_TIME_PART_PDF(2, J, IDP, IPDF), & !- vp
          MEAN_TIME_PART_PDF(3, J, IDP, IPDF), & !- wp
          MEAN_TIME_PART_PDF(4, J, IDP, IPDF), & !- oxp
          MEAN_TIME_PART_PDF(5, J, IDP, IPDF), & !- oyp
          MEAN_TIME_PART_PDF(6, J, IDP, IPDF), & !- ozp
          MEAN_TIME_PART_PDF(7, J, IDP, IPDF), & !- global.oxp         
          MEAN_TIME_PART_PDF(8, J, IDP, IPDF), & !- global.oyp
          MEAN_TIME_PART_PDF(9, J, IDP, IPDF), & !- global.ozp
          MEAN_TIME_PART_PDF(10,J, IDP, IPDF), & !- px
          MEAN_TIME_PART_PDF(11,J, IDP, IPDF), & !- py
          MEAN_TIME_PART_PDF(12,J, IDP, IPDF)    !- pz

        end do

      end do

  end do


  MEAN_TIME_PART_COL_PDF =  MEAN_TIME_PART_COL_PDF/real(NEVEN_COLL_PDF)

  do J = 1, NIG

    do IDP = 1, POLYDISP(J)

      NUMFILE = 650+(J-1)*POLYDISP(J)+IDP

      do IPDF = 1, NPDF

        write(NUMFILE,10000) real(IPDF), &
        MEAN_TIME_PART_COL_PDF(1, J, IDP, IPDF),  & !
        MEAN_TIME_PART_COL_PDF(2, J, IDP, IPDF),  & !
        MEAN_TIME_PART_COL_PDF(3, J, IDP, IPDF),  & !
        MEAN_TIME_PART_COL_PDF(4, J, IDP, IPDF),  & !
        MEAN_TIME_PART_COL_PDF(5, J, IDP, IPDF),  & !
        MEAN_TIME_PART_COL_PDF(6, J, IDP, IPDF),  & !
        MEAN_TIME_PART_COL_PDF(7, J, IDP, IPDF),  & !
        MEAN_TIME_PART_COL_PDF(8, J, IDP, IPDF),  & !
        MEAN_TIME_PART_COL_PDF(9, J, IDP, IPDF),  & !
        MEAN_TIME_PART_COL_PDF(10, J, IDP, IPDF), & !
        MEAN_TIME_PART_COL_PDF(11, J, IDP, IPDF), & !
        MEAN_TIME_PART_COL_PDF(12, J, IDP, IPDF)    !


      end do

    end do

  end do


end if



!!======================================================================
!! 4. Lagrangian correlation function
!!======================================================================
!!
!! The normalized autocorrelation function is computed as
!!
!!                 <u(t0).u(t0+t)>
!! R(t) = ---------------------------------
!!        sqrt(<u(t0)^2>).sqrt(<u(t0+t)^2>)
!!
!!====================================================================
!!- The following procedure works only if 
!!  the end of computation is normal

if(LEVEL2_STPAR) then


 do NLGR = 1, NBLGRMAX
 
!- Index
   NTLGR = int( (NCYCLE - NT0(NLGR)*FOUT2)/FOUT2) + 1
  
  
  if(NCYCLE >= NT0(NLGR)*FOUT2) then

!!--------------------------------------------------------------------
!! 4.1 Storage Lagrangian fluctuating velocity at t=t0
!!--------------------------------------------------------------------
   if(NCYCLE == NT0(NLGR)*FOUT2) then

   if(MYID==0) write(*,*)'Lagrangian function computing starts at NCYCLE= ',NT0(NLGR)*FOUT2

    do J = 1, NIG
     do I = 1, NPART_LOC(J)

      !- index of PSD
      IDP = PART(I,J)%IDP

      !- Particle velocity 
      PART(I,J)%UPT0(NLGR) = PART(I,J)%UP - MEAN_PART(1,J,IDP)
      PART(I,J)%VPT0(NLGR) = PART(I,J)%VP - MEAN_PART(2,J,IDP)
      PART(I,J)%WPT0(NLGR) = PART(I,J)%WP - MEAN_PART(3,J,IDP)

      !- Fluid velocity at particle position
      !PART(I,J)%UFAPT0(NLGR) = PART(I,J)%UFAP - MEAN_PART(11,J,IDP)
      !PART(I,J)%VFAPT0(NLGR) = PART(I,J)%VFAP - MEAN_PART(12,J,IDP)
      !PART(I,J)%WFAPT0(NLGR) = PART(I,J)%WFAP - MEAN_PART(13,J,IDP)
 
     end do
    end do !->  J = 1, NIG

   end if !-> (NCYCLE == NT0*FOUT2)

!!--------------------------------------------------------------------
!! 4.2 correlation at time t
!!--------------------------------------------------------------------

   TIME_LGR(NTLGR,NLGR) = TIME

   do J = 1, NIG
    do I = 1, NPART_LOC(J)

!!- index of PSD
     IDP = PART(I,J)%IDP

!!-------------------------------------------------------------------
!! 4.1 Particle-particle
!!-------------------------------------------------------------------
!- <up(t0).up(t0+t)>
     RPX_LOC(NTLGR,J,NLGR,IDP) = RPX_LOC(NTLGR,J,NLGR,IDP) &
              + (PART(I,J)%UP - MEAN_PART(1,J,IDP))*PART(I,J)%UPT0(NLGR) /NPM(IDP)

!- <vp(t0).vp(t0+t)>
     RPY_LOC(NTLGR,J,NLGR,IDP) = RPY_LOC(NTLGR,J,NLGR,IDP) &
              + (PART(I,J)%VP - MEAN_PART(2,J,IDP))*PART(I,J)%VPT0(NLGR) /NPM(IDP)

!- <wp(t0).wp(t0+t)>
     RPZ_LOC(NTLGR,J,NLGR,IDP) = RPZ_LOC(NTLGR,J,NLGR,IDP) &
              + (PART(I,J)%WP - MEAN_PART(3,J,IDP))*PART(I,J)%WPT0(NLGR) /NPM(IDP)


!!-------------------------------------------------------------------
!! 4.2 Fluid-particle
!!-------------------------------------------------------------------
!- <ufap(t0).ufap(t0+t)>
!     RFAPX_LOC(NTLGR,J,NLGR,IDP) = RFAPX_LOC(NTLGR,J,NLGR,IDP) &
!              + (PART(I,J)%UFAP-MEAN_PART(11,J,IDP))*PART(I,J)%UFAPT0(NLGR) /NPM(IDP)

!- <vfap(t0).vfap(t0+t)>
!     RFAPY_LOC(NTLGR,J,NLGR,IDP) = RFAPY_LOC(NTLGR,J,NLGR,IDP) &
!              + (PART(I,J)%VFAP-MEAN_PART(12,J,IDP))*PART(I,J)%VFAPT0(NLGR) /NPM(IDP)

!- <wfap(t0).wfap(t0+t)>
!     RFAPZ_LOC(NTLGR,J,NLGR,IDP) = RFAPZ_LOC(NTLGR,J,NLGR,IDP) &
!              + (PART(I,J)%WFAP-MEAN_PART(13,J,IDP))*PART(I,J)%WFAPT0(NLGR) /NPM(IDP)


    end do !-> I = 1, NPART_LOC(J)
   end do !->  J = 1, NIG

  end if !!- end if(NCYCLE >= NT0(NLGR)*FOUT2)
 end do !!- end do NLGR=1, NBLGRMAX

end if !!-> (LEVEL2_STPAR)


!!======================================================================
!! 5. Two-particles statistics
!!======================================================================
!if(LEVEL3_STPAR)  call RDF_PARTICLE




!!======================================================================
!! 7. Time averaging if so
!!======================================================================
if(STAT_TIME) then
 MEAN_TIME_PART = MEAN_TIME_PART + MEAN_PART
end if 



!!======================================================================
!! 7. Print statistics in "info" files
!!======================================================================
!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()
 CPU_PART(6) = CPU_PART(6) + TIME_END - TIME_START
end if

!!----------------------------------------------------------------------
10000 format (40(e17.7))
10601 format (2x,A,E13.6)
10602 format (2x,A,E13.6,1x,A,E13.6)
10603 format (2x,A,E13.6,2x,A,E13.6,2x,A,E13.6)

end subroutine STAT_PARTICLE
