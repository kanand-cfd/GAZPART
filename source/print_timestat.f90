!!====================================================================
!!
!!
!!====================================================================

subroutine PRINT_TIMESTAT(NCYCLE)

!!====================================================================
!!
!!
!!====================================================================
use DNS_DIM            !- Dimension
use STATISTICS         !- Statistics
use PARAM_PHYS
use GEOMETRIC_VARIABLE
use COLLISION_VARIABLE
use PARTICLE_PARALLEL


implicit none

!!--------------------------------------------------------------------
!! ARRAYS STATEMENT
!!--------------------------------------------------------------------
integer, intent(in) :: NCYCLE

!!--------------------------------------------------------------------
!!- Local arrays
!!--------------------------------------------------------------------

!!- Fluid spatial correlation function
real(kind=8), dimension(DIMSCOR) :: RUX, RVX, RWX

!!- Radial distribution function
real(kind=8), dimension(NR) :: RDF,    RDFIN,    RDFOUT
real(kind=8), dimension(NR) :: RDFWR,  RDFWRIN,  RDFWROUT
real(kind=8), dimension(NR) :: RDFTHR, RDFTHRIN, RDFTHROUT
real(kind=8) :: VSHELL

!- Solid volume fraction and collision timescale from theory
real(kind=8) :: TKIN, FKIN

real(kind=8), dimension(POLYDISPMAX) :: NP

real(kind=8) :: TCOL, FCOL, DCOL

real(kind=8) :: COEF1, COEF2

!- File name
character(len=40) :: FILENAME

! File number
integer :: NUMFILE

!- index
integer :: IDP1, IDP2

!- Index
integer :: I, J, LGR, M ,N, K
!---------------------------------------------------------------------


!!====================================================================
!! 1. Fluid statistics normalization
!!====================================================================
if(LEVEL_STFLU>0) then

!!--------------------------------------------------------------------
!! 1.1. Time-averaged statistics
!!--------------------------------------------------------------------
 MEAN_TIME_FLUID = MEAN_TIME_FLUID / real(NEVEN)

!!--------------------------------------------------------------------
!! 1.2. Spatial Correlation functions
!!--------------------------------------------------------------------
if(LEVEL_STFLU==4) then
 call MPI_ALLREDUCE(MEAN_RUXLOC,RUX, &
                    DIMSCOR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
 call MPI_ALLREDUCE(MEAN_RVXLOC,RVX, &
                    DIMSCOR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
 call MPI_ALLREDUCE(MEAN_RWXLOC,RWX, &
                    DIMSCOR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
!
!
 if(MYID==0) then
  do I = 1, DIMSCOR 
   write(304,10000)XMESH(I),    RUX(I)                         &
                            ,0.5*(RVX(I)+RWX(I))               &
                            ,     RUX(I)        / real(NEVEN)  &
                            ,0.5*(RVX(I)+RWX(I))/ real(NEVEN)  &
                            ,     RUX(I)        /RUX(1)        &
                            ,(RVX(I)+RWX(I))/(RVX(1)+RWX(1))
  end do
 end if

end if !- endi f (LEVEL_STFLU==4)

end if !- end if(LEVEL_STFLU>0)

!!====================================================================
!! 2. Scalar statistics normalization
!!====================================================================
if(LEVEL0_STSCL) then

!!--------------------------------------------------------------------
!! 2.1. Time-averaged statistics
!!--------------------------------------------------------------------
 MEAN_TIME_SCL = MEAN_TIME_SCL / real(NEVEN)

end if


!!====================================================================
!! 3. Particle statistics normalization
!!====================================================================
if(LEVEL0_STPAR) then

!!--------------------------------------------------------------------
!! 3.1. One point statistics
!!--------------------------------------------------------------------
 MEAN_TIME_PART = MEAN_TIME_PART / real(NEVEN)

 !if(FROZEN_FLOW>0) MEAN_TIME_PARTFLUID = MEAN_TIME_PARTFLUID / real(NEVEN)

!- Statistics of collision
 if(SOLVE_COLLISION)  MEAN_TIME_PART_COL = MEAN_TIME_PART_COL / real(NEVEN) 


!!--------------------------------------------------------------------
!! 3.2. Radial Distribution Function
!!--------------------------------------------------------------------
 if(LEVEL3_STPAR) then
 
  do J = 1, NIG

   if ((PARTDEF(J)==2).or.((PARTDEF(J)==1).and.(SOLVE_FLUID>0))) then

!- Number of particle per class
!    NP_LOC(:) = 0
!    do I = 1, NPART_LOC(J)
!     NP_LOC(PART(I,J)%IDP) = NP_LOC(PART(I,J)%IDP) + 1.0
!    end do
!    call MPI_ALLREDUCE(NP_LOC,NP,POLYDISPMAX,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)

    do IDP1 = 1, POLYDISP(J)
     do IDP2 = 1, POLYDISP(J)

     call MPI_ALLREDUCE(RDF_LOC(:,J,IDP1,IDP2),    RDF,    NR, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
     call MPI_ALLREDUCE(RDFWR_LOC(:,J,IDP1,IDP2),  RDFWR,  NR, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
     call MPI_ALLREDUCE(RDFTHR_LOC(:,J,IDP1,IDP2), RDFTHR, NR, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)

     call MPI_ALLREDUCE(RDFIN_LOC(:,J,IDP1,IDP2),    RDFIN,    NR, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
     call MPI_ALLREDUCE(RDFWRIN_LOC(:,J,IDP1,IDP2),  RDFWRIN,  NR, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
     call MPI_ALLREDUCE(RDFTHRIN_LOC(:,J,IDP1,IDP2), RDFTHRIN, NR, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)

     call MPI_ALLREDUCE(RDFOUT_LOC(:,J,IDP1,IDP2),    RDFOUT,    NR, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
     call MPI_ALLREDUCE(RDFWROUT_LOC(:,J,IDP1,IDP2),  RDFWROUT,  NR, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)
     call MPI_ALLREDUCE(RDFTHROUT_LOC(:,J,IDP1,IDP2), RDFTHROUT, NR, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, IERR)

    !normalizing by number of particles within the respective shell
     RDFWR  = RDFWR /RDF
     RDFTHR = RDFTHR/RDF
    !normalizing by number of approaching particles within the respective shell
     RDFWRIN  = RDFWRIN /RDFIN
     RDFTHRIN = RDFTHRIN/RDFIN
    !normalizing by number of departing particles within the respective shell
     RDFWROUT  = RDFWROUT /RDFOUT
     RDFTHROUT = RDFTHROUT/RDFOUT

     COEF1 = NTEST
     COEF2 = NPART_FULL
     if(POLYDISP(J) ==2) then
       COEF1 = NTEST
       COEF2 = NPCLASS(J,IDP2)
     end if
   
    !normalizing by total particle number fraction 
    RDF    = RDF   /real(COEF2-1.0)*LXMAX*LYMAX*LZMAX
    RDFIN  = RDFIN /real(COEF2-1.0)*LXMAX*LYMAX*LZMAX
    RDFOUT = RDFOUT/real(COEF2-1.0)*LXMAX*LYMAX*LZMAX

    !normalizing by number of test particles and statistical events
    RDF     = RDF    /real(COEF1)/real(NEVEN)
    RDFIN   = RDFIN  /real(COEF1)/real(NEVEN)
    RDFOUT  = RDFOUT /real(COEF1)/real(NEVEN)

    write(FILENAME,10102)'part_l3_rdf_p',J,'_c',IDP1,IDP2,'.stat'
    open(unit=600, file=trim(FILENAME), status='replace')
    write(600,40001)

    do I = 1, NR
     VSHELL = 4./3.*PPI*((RMIN(J)+I*DR(J))**3 - (RMIN(J)+(I-1)*DR(J))**3)
     RDF(I)    = RDF(I)   /VSHELL
     RDFIN(I)  = RDFIN(I) /VSHELL
     RDFOUT(I) = RDFOUT(I)/VSHELL
     write(600,10000) (RMIN(J)+(I-0.5)*DR(J)), RDF(I), RDFIN(I), RDFOUT(I), RDFWR(I), &
                            RDFWRIN(I), RDFWROUT(I), RDFTHR(I), RDFTHRIN(I), RDFTHROUT(I)
    end do
    close(600)

    end do !- end do IDP2 =
    end do !- end do IDP1 =


   end if !!- end if ((PARTDEF(J)==2).or. ...

  end do !!- end do J = 1, NIG

 end if !!- end if(LEVEL3_STPAR)

end if !!- end if(LEVEL0_STPAR)


!!---------------------------------------------------
!! Probability Distribution Function for Collision
!!---------------------------------------------------
!if (LEVEL1_STPAR .and. SOLVE_COLLISION) then
!  do J = 1, NIG
!    do IDP1 = 1, POLYDISP(J)
!      do IDP2 = 1, POLYDISP(J)
!        NUMFILE = 7000 + 100*J + 10*IDP1 + IDP2 
!        do I = 1, NPDFMAX
!          if (MYID == 0) write(NUMFILE, 10000)real(I), PDFCOL(J,IDP1,IDP2,I,1), PDFCOL(J,IDP1,IDP2,I,2), PDFCOL(J,IDP1,IDP2,I,3)
!        end do
!      end do  
!    end do
!  end do
!end if


!!====================================================================
!! 3. Print in "stat.scilab" for matlab or scilab post-processing
!!====================================================================
!!call PRINT_SCILAB



!!====================================================================
!! 4. Print in "run.info"
!!====================================================================
if(MYID==0) then 

I = 1
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'====================================================================='
write(UNIT_INFO(I),*)'TIME AVERAGED STATISTICS '
write(UNIT_INFO(I),*)'====================================================================='
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10600)'Number of computed cycle      =',real(NCYCLE)
write(UNIT_INFO(I),10600)'Number of event for statistics=',real(NEVEN)
write(UNIT_INFO(I),*)

if(LEVEL_STFLU>0) then
write(UNIT_INFO(I),*)'---------------------------------------------------------------------'
write(UNIT_INFO(I),*)' FLUID'
write(UNIT_INFO(I),*)'---------------------------------------------------------------------'
write(UNIT_INFO(I),*)
end if

if(LEVEL_STFLU>=1) then
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Mean fluid Velocity'
write(UNIT_INFO(I),*)'-------------------'
write(UNIT_INFO(I),10601)'<uf>f = ',MEAN_TIME_FLUID(1),' [m/s]'
write(UNIT_INFO(I),10601)'<vf>f = ',MEAN_TIME_FLUID(2),' [m/s]'
write(UNIT_INFO(I),10601)'<wf>f = ',MEAN_TIME_FLUID(3),' [m/s]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Reynolds stress tensor'
write(UNIT_INFO(I),*)'----------------------'
write(UNIT_INFO(I),10601)'<uf.uf>f = ',MEAN_TIME_FLUID(4),' [m2/s2]'
write(UNIT_INFO(I),10601)'<vf.vf>f = ',MEAN_TIME_FLUID(5),' [m2/s2]'
write(UNIT_INFO(I),10601)'<wf.wf>f = ',MEAN_TIME_FLUID(6),' [m2/s2]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'    <kf> = ',MEAN_TIME_FLUID(10),' [m2/s2]'
write(UNIT_INFO(I),10601)'  <epsf> = ',MEAN_TIME_FLUID(11),' [m2/s3]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*    )'   u_rms=(2/3*kf)**0.5'
write(UNIT_INFO(I),10601)'   u_rms = ',sqrt(2./3.*MEAN_TIME_FLUID(10)),' [m/s]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'<uf.vf>f = ',MEAN_TIME_FLUID(7),' [m2/s2]'
write(UNIT_INFO(I),10601)'<uf.wf>f = ',MEAN_TIME_FLUID(8),' [m2/s2]'
write(UNIT_INFO(I),10601)'<vf.wf>f = ',MEAN_TIME_FLUID(9),' [m2/s2]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'   eta_K = ',MEAN_TIME_FLUID(12),' [m]'
write(UNIT_INFO(I),10601)'   tau_K = ',MEAN_TIME_FLUID(13),' [s]'
write(UNIT_INFO(I),10601)'     v_K = ',MEAN_TIME_FLUID(17),' [m/s]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'Taylor      lg = ',MEAN_TIME_FLUID(25),' [m]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Reynolds number'
write(UNIT_INFO(I),*)'---------------'
write(UNIT_INFO(I),10602)' Re_Tay = ',sqrt(2.*MEAN_TIME_FLUID(10)/3.)*MEAN_TIME_FLUID(25)/VISC
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)
end if

if(LEVEL_STFLU>=2) then
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Divergence'
write(UNIT_INFO(I),*)'----------'
write(UNIT_INFO(I),10601)'   < div(U)  > = ',MEAN_TIME_FLUID(15),' [1/s]'
write(UNIT_INFO(I),10601)'   < div(U)^2> = ',MEAN_TIME_FLUID(16),' [1/s2]'
write(UNIT_INFO(I),10601)'max(|div(U)| ) = ',MEAN_TIME_FLUID(14),' [1/s]'
write(UNIT_INFO(I),*)
end if

if(LEVEL_STFLU>=3) then
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Velocity gradients'
write(UNIT_INFO(I),*)'------------------'
write(UNIT_INFO(I),10601)' <(dui/dxi)^2> = ',MEAN_TIME_FLUID(21),' [1/s2]'
write(UNIT_INFO(I),10601)' <(dui/dxj)^2> = ',MEAN_TIME_FLUID(31),' [1/s2]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10602)'Skewness:   Sk = ',MEAN_TIME_FLUID(18)
write(UNIT_INFO(I),10602)'Flatness:   Tk = ',MEAN_TIME_FLUID(19)
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'Taylor      lf = ',MEAN_TIME_FLUID(24),' [m]'
write(UNIT_INFO(I),10601)'Taylor      lg = ',MEAN_TIME_FLUID(25),' [m]'
write(UNIT_INFO(I),10602)'   lf/lg/2^0.5 = ',MEAN_TIME_FLUID(24)/MEAN_TIME_FLUID(25)/sqrt(2.)
write(UNIT_INFO(I),*)
end if

write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)


if(LEVEL0_STPAR) then
write(UNIT_INFO(I),*)'---------------------------------------------------------------------'
write(UNIT_INFO(I),*)' PARTICLE'
write(UNIT_INFO(I),*)'---------------------------------------------------------------------'
write(UNIT_INFO(I),*)
end if

if(LEVEL1_STPAR) then

  call TIME_AVERAGED_STATISTICS

  do J = 1, NIG
    write(UNIT_INFO(I),*)
    write(UNIT_INFO(I),*)'---------------------------------------------------------------------'

    if(PARTDEF(J) == 0) then
      write(UNIT_INFO(I),10603)'Particle:',J,'  ==> Motionless'

    elseif(PARTDEF(J) == 1) then
      write(UNIT_INFO(I),10603)'Particle:',J,'  ==> Fluid element'

    elseif(PARTDEF(J) == 2) then
      write(UNIT_INFO(I),10603)'Particle:',J,'  ==> Solid particle'

    end if

    !- Particle number density
    do IDP1 = 1, POLYDISP(J)
    
      NP(IDP1) = NPCLASS(J,IDP1)/LXMAX/LYMAX/LZMAX
    
    end do

    do IDP1 = 1, POLYDISP(J)

      write(UNIT_INFO(I),*)'---------------------------------------------------------------------'
      write(UNIT_INFO(I),*)' Class ',IDP1
      write(UNIT_INFO(I),*)'---------------------------------------------------------------------'
      write(UNIT_INFO(I),*)
      write(UNIT_INFO(I),10601)'  Np=',MEAN_TIME_PART(23,J,IDP1),' [-]'
      write(UNIT_INFO(I),10603)'  Np=',NPCLASS(J,IDP1),' [-]'
      !write(UNIT_INFO(I),10601)'  dp=',DPART(J,IDP1),' [m]'
      write(UNIT_INFO(I),10601)'Major Axis=',EMAJ_PART(J, IDP1),' [m]'
      write(UNIT_INFO(I),10601)'Aspect Ratio=',APR_PART(J),' [-]'
      write(UNIT_INFO(I),10601)'rhop=',RHOP(J,IDP1),' [kg/m3]'
      write(UNIT_INFO(I),*)
      write(UNIT_INFO(I),*)'Mean Velocity'
      write(UNIT_INFO(I),*)'-------------'
      write(UNIT_INFO(I),10601)'<up>p = ',MEAN_TIME_PART(1,J,IDP1),' [m/s]'
      write(UNIT_INFO(I),10601)'<vp>p = ',MEAN_TIME_PART(2,J,IDP1),' [m/s]'
      write(UNIT_INFO(I),10601)'<wp>p = ',MEAN_TIME_PART(3,J,IDP1),' [m/s]'
      write(UNIT_INFO(I),*)
      write(UNIT_INFO(I),10601)'<uf>p = ',MEAN_TIME_PART(11,J,IDP1),' [m/s]'
      write(UNIT_INFO(I),10601)'<vf>p = ',MEAN_TIME_PART(12,J,IDP1),' [m/s]'
      write(UNIT_INFO(I),10601)'<wf>p = ',MEAN_TIME_PART(13,J,IDP1),' [m/s]'
      write(UNIT_INFO(I),*)
      write(UNIT_INFO(I),*)'Fluctuating motion'
      write(UNIT_INFO(I),*)'------------------'
      write(UNIT_INFO(I),10601)'  qp = ',MEAN_TIME_PART(10,J,IDP1),' [m2/s2]'
      write(UNIT_INFO(I),10601)'qf@p = ',MEAN_TIME_PART(20,J,IDP1),' [m2/s2]'
      write(UNIT_INFO(I),10601)' qfp = ',MEAN_TIME_PART(27,J,IDP1),' [m2/s2]'
      write(UNIT_INFO(I),*)
      write(UNIT_INFO(I),*)
      write(UNIT_INFO(I),*)'Particle-Particle kinetic stress tensor'
      write(UNIT_INFO(I),*)'--------------------------------------'
      write(UNIT_INFO(I),10601)'<upup>p = ',MEAN_TIME_PART(4,J,IDP1),' [m2/s2]'
      write(UNIT_INFO(I),10601)'<vpvp>p = ',MEAN_TIME_PART(5,J,IDP1),' [m2/s2]'
      write(UNIT_INFO(I),10601)'<wpwp>p = ',MEAN_TIME_PART(6,J,IDP1),' [m2/s2]'
      write(UNIT_INFO(I),*)
      write(UNIT_INFO(I),10601)'<upvp>p = ',MEAN_TIME_PART(7,J,IDP1),' [m2/s2]'
      write(UNIT_INFO(I),10601)'<upwp>p = ',MEAN_TIME_PART(8,J,IDP1),' [m2/s2]'
      write(UNIT_INFO(I),10601)'<vpwp>p = ',MEAN_TIME_PART(9,J,IDP1),' [m2/s2]'
      write(UNIT_INFO(I),*)
      
      write(UNIT_INFO(I),10601)'<upup>p/(2/3qp) = ',&
                               MEAN_TIME_PART(4,J,IDP1)/(2./3.*MEAN_TIME_PART(10,J,IDP1)),' [-]'
  
      write(UNIT_INFO(I),10601)'<vpvp>p/(2/3qp) = ',&
                               MEAN_TIME_PART(5,J,IDP1)/(2./3.*MEAN_TIME_PART(10,J,IDP1)),' [-]'
  
      write(UNIT_INFO(I),10601)'<wpwp>p/(2/3qp) = ',&
                               MEAN_TIME_PART(6,J,IDP1)/(2./3.*MEAN_TIME_PART(10,J,IDP1)),' [-]'
  
      write(UNIT_INFO(I),*)
      
      write(UNIT_INFO(I),10601)'<upvp>p/(2/3qp) = ',&
                               MEAN_TIME_PART(7,J,IDP1)/(2./3.*MEAN_TIME_PART(10,J,IDP1)),' [-]'
  
      write(UNIT_INFO(I),10601)'<upwp>p/(2/3qp) = ',&
                               MEAN_TIME_PART(8,J,IDP1)/(2./3.*MEAN_TIME_PART(10,J,IDP1)),' [-]'
  
      write(UNIT_INFO(I),10601)'<vpwp>p/(2/3qp) = ',&
                               MEAN_TIME_PART(9,J,IDP1)/(2./3.*MEAN_TIME_PART(10,J,IDP1)),' [-]'
  
      write(UNIT_INFO(I),*)
      write(UNIT_INFO(I),*)
  
  
      if(PARTDEF(J)>1) then

      
        if(SOLVE_COLLISION) then
          write(UNIT_INFO(I),*)
          write(UNIT_INFO(I),*)'Collision statistics'
          write(UNIT_INFO(I),*)'--------------------------------------------------'

          do IDP2 = 1, POLYDISP(J)

            write(UNIT_INFO(I),10604)'        Nover[',IDP1,',',IDP2,'] = ',MEAN_TIME_PART_COL(1,IDP1,IDP2),' [-]'
            write(UNIT_INFO(I),10604)'        N_col[',IDP1,',',IDP2,'] = ',MEAN_TIME_PART_COL(2,IDP1,IDP2),' [-]'
            !write(UNIT_INFO(I),10604)'   <|w.k|>col[',IDP1,',',IDP2,'] = ',MEAN_TIME_PART_COL(4,IDP1,IDP2),' [m/s]'
            !write(UNIT_INFO(I),10604)'     <|w|>col[',IDP1,',',IDP2,'] = ',MEAN_TIME_PART_COL(5,IDP1,IDP2),' [m/s]'
            !write(UNIT_INFO(I),10604)'   <theta>col[',IDP1,',',IDP2,'] = ',MEAN_TIME_PART_COL(6,IDP1,IDP2),' [rad]'
            !write(UNIT_INFO(I),10604)' <theta^2>col[',IDP1,',',IDP2,'] = ',MEAN_TIME_PART_COL(7,IDP1,IDP2),' [rad]'
            !write(UNIT_INFO(I),10604)'       qp_col[',IDP1,',',IDP2,'] = ',MEAN_TIME_PART_COL(8,IDP1,IDP2),' [m2/s2]'
            !write(UNIT_INFO(I),10604)'    delta[qp][',IDP1,',',IDP2,'] = ',MEAN_TIME_PART_COL(9,IDP1,IDP2),' [m2/s2]'
            write(UNIT_INFO(I),10604)'        f_col[',IDP1,',',IDP2,'] = ',MEAN_TIME_PART_COL(3,IDP1,IDP2)/MEAN_TIME_PART(23,J,IDP1),' [1/s]'
            write(UNIT_INFO(I),10604)'      tau_col[',IDP1,',',IDP2,'] = ',MEAN_TIME_PART(23,J,IDP1)/MEAN_TIME_PART_COL(3,IDP1,IDP2),' [s]'
            write(UNIT_INFO(I),10604)'      Overlap[',IDP1,',',IDP2,'] = ',MEAN_TIME_PART_COL(4,IDP1,IDP2),' [-]'
            write(UNIT_INFO(I),10604)'Overlap/Nover[',IDP1,',',IDP2,'] = ',MEAN_TIME_PART_COL(5,IDP1,IDP2),' [-]'
            !write(UNIT_INFO(I),10604)'      tau_col[',IDP1,',',IDP2,'] = ',1./MEAN_TIME_PART_COL(10,IDP1,IDP2),' [1/s]'
            !write(UNIT_INFO(I),10604)'      tau_kin[',IDP1,',',IDP2,'] = ',TAUK,' [1/s]'
            !write(UNIT_INFO(I),10604)'           g0[',IDP1,',',IDP2,'] = ',TAUK*MEAN_TIME_PART_COL(10,IDP1,IDP2),' [-]'
            !write(UNIT_INFO(I),10604)'        dl/dp[',IDP1,',',IDP2,'] = ',1.5*sqrt(2.0*PPI/3.0*MEAN_TIME_PART(10,J,IDP1))*DTIME/DPART(J),' [1/s]'
            !write(UNIT_INFO(I),*)

          end do

        end if !- end if(SOLVE_COLLISION)

      end if !- end if(PARTDEF(J)>1)

!!- Spatial distribution statistics


end do !- end IDP1

end do !!- end loop J = 1, NIG

end if

!if(LEVELX_STPAR) then

!write(*,*) ' '
!write(*,*) 'Before Normalisation '
!write(*,*) 'Impulse at Left Wall (X=0) ', MEAN_PC_X(1), MEAN_PC_Y(1), MEAN_PC_Z(1), MEAN_PC_NORM(1)
!write(*,*) 'Impulse at Right Wall (X=LX) ', MEAN_PC_X(2), MEAN_PC_Y(2), MEAN_PC_Z(2), MEAN_PC_NORM(2)
!write(*,*) 'NEVEN_WALL', NEVEN_WALL, ' NEVEN_L ', NEVEN_L, 'NEVEN_R ', NEVEN_R
!write(*,*) ' '
!write(*,*) 'After Normalisation by NEVEN_L = ', NEVEN_L
!write(*,*) 'Impulse at Left Wall (X=0) ', MEAN_PC_X(1)/NEVEN_L, MEAN_PC_Y(1)/NEVEN_L, MEAN_PC_Y(1)/NEVEN_L, MEAN_PC_NORM(1)/NEVEN_L
!write(*,*) ' '
!write(*,*) 'After Normalisation by NEVEN_R = ', NEVEN_R
!write(*,*) 'Impulse at Right Wall (X=LX) ', MEAN_PC_X(2)/NEVEN_R, MEAN_PC_Y(2)/NEVEN_R, MEAN_PC_Y(2)/NEVEN_R, MEAN_PC_NORM(2)/NEVEN_R
!write(*,*) ' '

!write(*,*) 'After Normalisation by NEVEN_WALL'
!write(*,*) 'Impulse at Left Wall ', MEAN_PC_X(1)/NEVEN_WALL, MEAN_PC_Y(1)/NEVEN_WALL, MEAN_PC_Y(1)/NEVEN_WALL, MEAN_PC_NORM(1)/NEVEN_WALL
!write(*,*) 'Impulse at Right Wall', MEAN_PC_X(2)/NEVEN_WALL, MEAN_PC_Y(2)/NEVEN_WALL, MEAN_PC_Y(2)/NEVEN_WALL, MEAN_PC_NORM(2)/NEVEN_WALL
!write(*,*) ' '

!end if

!!====================================================================
!! Print statistic for postprocessing
!!====================================================================
!open(unit=900, file='data.mean', status='replace')
!write(900,*)'#   Particle class 1 2 3 ....'
!if(LEVEL1_STPAR) then
!write(900,10001)(MEAN_TIME_PART( 1,J),J=1,NIG) !- <up>p
!write(900,10001)(MEAN_TIME_PART( 2,J),J=1,NIG) !- <vp>p
!write(900,10001)(MEAN_TIME_PART( 3,J),J=1,NIG) !- <wp>p
!write(900,10001)(MEAN_TIME_PART(11,J),J=1,NIG) !- <uf>p
!write(900,10001)(MEAN_TIME_PART(12,J),J=1,NIG) !- <vf>p
!write(900,10001)(MEAN_TIME_PART(13,J),J=1,NIG) !- <wf>p
!write(900,10001)(MEAN_TIME_PART(10,J),J=1,NIG) !-  qp
!write(900,10001)(MEAN_TIME_PART(20,J),J=1,NIG) !-  qf@p
!write(900,10001)(MEAN_TIME_PART(27,J),J=1,NIG) !-  qfp
!write(900,10001)(MEAN_TIME_PART( 4,J),J=1,NIG) !- <upup>p
!write(900,10001)(MEAN_TIME_PART( 5,J),J=1,NIG) !- <vpvp>p
!write(900,10001)(MEAN_TIME_PART( 6,J),J=1,NIG) !- <wpwp>p
!write(900,10001)(MEAN_TIME_PART( 7,J),J=1,NIG) !- <upvp>p
!write(900,10001)(MEAN_TIME_PART( 8,J),J=1,NIG) !- <upwp>p
!write(900,10001)(MEAN_TIME_PART( 9,J),J=1,NIG) !- <vpwp>p
!write(900,10001)(MEAN_TIME_PART(14,J),J=1,NIG) !- <ufuf>p
!write(900,10001)(MEAN_TIME_PART(15,J),J=1,NIG) !- <vfvf>p
!write(900,10001)(MEAN_TIME_PART(16,J),J=1,NIG) !- <wfwf>p
!write(900,10001)(MEAN_TIME_PART(17,J),J=1,NIG) !- <ufvf>p
!write(900,10001)(MEAN_TIME_PART(18,J),J=1,NIG) !- <ufwf>p
!write(900,10001)(MEAN_TIME_PART(19,J),J=1,NIG) !- <vfwf>p
!write(900,10001)(MEAN_TIME_PART(21,J),J=1,NIG) !- <ufup>p
!write(900,10001)(MEAN_TIME_PART(22,J),J=1,NIG) !- <vfvp>p
!write(900,10001)(MEAN_TIME_PART(23,J),J=1,NIG) !- <wfwp>p
!write(900,10001)(MEAN_TIME_PART(24,J),J=1,NIG) !- <ufvp>p
!write(900,10001)(MEAN_TIME_PART(25,J),J=1,NIG) !- <ufwp>p
!write(900,10001)(MEAN_TIME_PART(26,J),J=1,NIG) !- <vfwp>p
!write(900,10001)(MEAN_TIME_PART(28,J),J=1,NIG) !- <1/tfp_F>p
!write(900,10001)(MEAN_TIME_PART(31,J),J=1,NIG) !-  <Rep>p
!write(900,10001)(MEAN_TIME_PART(29,J),J=1,NIG) !-  <Vfp>p
!write(900,10001)(MEAN_TIME_PART(30,J),J=1,NIG) !-   <Cd>p

!if(LEVEL3_STPAR) then
!write(900,10001)(MEAN_TIME_PART(61,J),J=1,NIG) !- <C>
!write(900,10001)(MEAN_TIME_PART(63,J),J=1,NIG) !-<C^2> 
!write(900,10001)(MEAN_TIME_PART(62,J),J=1,NIG) !-<Cp^2>
!write(900,10001)((sqrt(MEAN_TIME_PART(62,J))-sqrt(MEAN_TIME_PART(61,J)))/MEAN_TIME_PART(61,J),J=1,NIG) !- sigma
!end if
!end if !!- end if LEVEL1_STPAR
!close(900)




end if !!- MyID==0









!!====================================================================
10000 format (40(e17.7))
10001 format (15(2x,e17.7))

10100 format (A,I2.2,A)
10102 format (A,I2.2,A,I1,I1,A)


10600 format (2x,A,f8.2)
10601 format (2x,A,E13.6,A)
10602 format (2x,A,E13.6)
10603 format (2x,A,I2,A)
10604 format (2x,A,I2,A,I2,A,E13.6,A)
10605 format (2x,A,I2,A,E13.6,A)

40001 format('# R, RDF, RDF_in, RDF_out, RDFWR, RDFWR_in, RDFWR_out, RDFTHR, RDFTHR_in, RDFTHR_out')


end subroutine PRINT_TIMESTAT
