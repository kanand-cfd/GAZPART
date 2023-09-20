!!====================================================================
!! 
!!  
!!====================================================================

program gaspart

!!====================================================================
!! Numerical code solving the 3-dimensionnal Navier-Stokes equations
!! using a spectral decomposition. The code is coupled with a 
!! Lagrangian particle tracking module.
!!
!! The numerical features are:
!!  * Fluid: 
!!       - 3d Navier-Stokes equation
!!       - Incompressible flow
!!       - 3rd order Adam-Bashforth (time-advancing)
!!       - Integrating factor (viscous terms)
!!
!!  * Scalar: 
!!       - Transport equations
!!       - 3rd order Adam-Bashforth (time-advancing)
!!       - Integrating factor (diffusive terms)
!!
!!  * Particles: 
!!       - fixed points, fluid elements, inertial particles
!!       - 2nd Order Runge-Kutta (time-advancing)
!!       - Integrating factor (diffusive terms)
!!====================================================================

use DNS_DIM            !- Dimension
use PARAM_PHYS         !- Physical & numerical parameters
use FORCING            !- Forcing
use FLUID_VARIABLE     !- Fluid velocity
use SCALAR_VARIABLE
use GEOMETRIC_VARIABLE !- 
use STATISTICS         !- Statistics
use WORK_ARRAYS
use CHECK_CPU          !- Variable for cpu time checking
use CPUTIME_CONTROL

use PARTICLE_PARALLEL

use MPI_STRUCTURES

use P3DFFT

implicit none


!---------------------------------------------------------------------
! ARGUMENT STATEMENT
!---------------------------------------------------------------------
integer :: IARGC,NB_ARGS
character(len=10) :: ARG

!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Current time
real(kind=8) :: TIME

!- Integer flag for subroutine argument
integer      :: IFLAG1, IFLAG2
logical      :: LFLAG1
real(kind=8) :: RDUMMY
integer      :: IDUMMY, IDUMMY2

!- Time control variable
real(kind=8) :: TIME_START, TIME_END

!- Time measure variable
real(kind=8) :: MEASURE_START, MEASURE_END

real ETIME          ! Declare the type of etime()
real ELAPSED(2)     ! For receiving user and system time
real ELAPSED2(2)     ! For receiving user and system time
real TOTAL          ! For receiving total time
real TOTAL2          ! For receiving total time

!- Stop flag
logical :: CONT, CONT_CPU

!- Index
integer :: I, J, K, NCYCLE

real(kind=8) :: TIME0, TIME1

!-
integer :: NOUT2, NOUT3

!- Temporary index of switching for Time integration scheme
integer :: ISAVE
!---------------------------------------------------------------------




!!====================================================================
!! Get walltime
!!====================================================================

! lecture eventuelle du decoupage sur la ligne de commande
NB_ARGS = IARGC()


if ( NB_ARGS /= 0 ) then
 call getarg(1,ARG)
 read(ARG,*) WALLTIME
else 
 WALLTIME = INFINITY
end if


DEBUG = .false.


!!====================================================================
!! 1. INITIATION
!!====================================================================

!!--------------------------------------------------------------------
!! 1.1 MPI World initiation
!!--------------------------------------------------------------------
call MPI_INIT(IERR)
call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,IERR)
call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)

if (MYID==0) write(*,*) 'WallTime =', WALLTIME



CPU_INIT = MPI_WTIME()


!- CPU
CPU_FLUID(:)=0.
CPU_ELAPSED = 0.
CPU_CYCLE = 0.
CPU_INITIATION = 0.


!!- Check cpu time
TIME_START = MPI_WTIME()


!!--------------------------------------------------------------------
!! 1.2 Read parameter
!!--------------------------------------------------------------------
call READPARAM

!!--------------------------------------------------------------------
!! 1.3 Run initiation
!!--------------------------------------------------------------------
call INIT_RUN


if(SOLVE_PART) CPU_PART(:)=ZERO


!!- Define the particle data structure
call MPI_PART_TYPE




NCYCLE = 1 !- Initiation of time cycle
TIME = 0.  !- Initiation of time


!!--------------------------------------------------------------------
!! 1.4 Domain splitting
!!--------------------------------------------------------------------
!- Dimentionality of the cpu decomposition (=1: slab, =2: squared)
NDIMS = 2

if(NDIMS == 1) then
 DIMS(1) = 1
 DIMS(2) = NPROC
else if(NDIMS == 2) then
 if (MYID==0) print *, 'Creating proc. grid with mpi_dims_create'
 DIMS(1) = 0
 DIMS(2) = 0
 call MPI_DIMS_CREATE(NPROC,2,DIMS,IERR)
 if(DIMS(1) > DIMS(2)) then
  DIMS(1) = DIMS(2)
  DIMS(2) = NPROC / DIMS(1)
 endif
endif

IPROC = DIMS(1)
JPROC = DIMS(2)

if(MYID == 0)write(*,*)'Using processor grid ',iproc,' x ',jproc


!!- FFT flag 
FFTFLAG = 'fft'


!!- FFt Initiation
call P3DFFT_SETUP(DIMS,NX,NY,NZ,MPI_COMM_WORLD,NX,NY,NZ,.TRUE.)


!!- Split the geometry
call P3DFFT_GET_DIMS(ISTART,IEND,ISIZE,1)
call P3DFFT_GET_DIMS(FSTART,FEND,FSIZE,2)


!!- Dimension initiation
NTOT = FSIZE(1)*FSIZE(2)*FSIZE(3)
NGLOB = NX * NY * NZ
FACTOR = 1.0D0/real(NGLOB)



!!- Create ghost cell if needed
if(NGHTCELL>0.or.SOLVE_PART) then

! Create MPI structures for ghost cells MPI-exchange
 call CREATE_MPI_VECTOR
! Find neighbouring for each processor
 call NEIGHBOURING

end if




!!--------------------------------------------------------------------
!! 1.5. Allocate arrays
!!--------------------------------------------------------------------
call ALLOCATE_ARRAYS

!!- Initiation particle CPU passing if needed
if(SOLVE_PART) call INITIATION_EXCHANGE

!!--------------------------------------------------------------------
!! 1.6. Mesh, wavenumbers, integrating factor and aliasing control
!!--------------------------------------------------------------------
call MESHING



!!--------------------------------------------------------------------
!! 1.7. Initiation of Physics
!!--------------------------------------------------------------------
!!- Fluid Initiation
if(SOLVE_FLUID == 1) call INITIATION_FLUID

if(FROZEN_FLOW>0) call INITIATION_FROZEN_FLOW


!- Update Ghost cells
if(NGHTCELL>0) then
 call FLUIDCOMM(UFLU)
 call FLUIDCOMM(VFLU)
 call FLUIDCOMM(WFLU)
end if



!!- Initiation forcing
!!if(STEADY) call INITIATION_FORCING
if(STEADY) call INITIATION_FORCING_NEW



!!- Scalar initiation
!if(SOLVE_SCALAR) then

! call INITIATION_SCALAR

 !- Update Ghost cells
! if(NGHTCELL>0) call FLUIDCOMM(THETA)

!end if

!!- Particle initiation
if(SOLVE_PART) then

!- Positions
 call INITIATION_PARTICLE_POSITION

!- Velocities
 call INITIATION_PARTICLE_VELOCITY

!- Collision solver is necessary
 if(SOLVE_COLLISION) call COLLISION_INITIATION

 if(STAT_TIME .and. LEVELX_STPAR) call INIT_CHANNEL_STAT

end if


!!- Time-averaging initiation
if(STAT_TIME) then
 if(LEVEL_STFLU>0) MEAN_TIME_FLUID = ZERO
 if(LEVEL0_STPAR) MEAN_TIME_PART = ZERO
 if(LEVEL0_STPAR) MEAN_TIME_PART_PDF = ZERO

 if(LEVEL0_STPAR .and. FROZEN_FLOW > 0) MEAN_TIME_PARTFLUID  = ZERO

 if(LEVEL0_STPAR .and. SOLVE_COLLISION) MEAN_TIME_PART_COL_PDF = ZERO

 if(LEVEL0_STSCL) MEAN_TIME_SCL = ZERO
 
 if(LEVELX_STPAR) then
    MEAN_TIME_PART_CHAN = ZERO
    MEAN_PC_NORM(:) = ZERO
    MEAN_PC_X(:) = ZERO
    MEAN_PC_Y(:) = ZERO
    MEAN_PC_Z(:) = ZERO
    NEVEN_L = 0
    NEVEN_R = 0
    MEAN_TIME_COLL_CHAN = ZERO
    !MEAN_TIME_PDF_CHAN = ZERO
    !MEAN_TIME_PDF_RELV_CHAN = ZERO
 end if
 
 NEVEN = 0
 NEVEN_CHAN = 0
 NEVEN_COLL_CHAN = 0
 NEVEN_PDF = 0
 NEVEN_COLL_PDF = 0
 
end if


!!--------------------------------------------------------------------
!! 1.8. Opening files
!!--------------------------------------------------------------------
LFLAG1 = .true.
call OPENCLOSE(LFLAG1)


!!--------------------------------------------------------------------
!! 1.9. Print info about the numerical simulation
!!--------------------------------------------------------------------
if(MYID==0) call INFORUN


!!- End initiation
TIME_END = MPI_WTIME()
CPU_INITIATION = TIME_END - TIME_START

CPU_ELAPSED = CPU_ELAPSED + CPU_INITIATION

if (MYID == 0) then

  write(*,10700)100.*NCYCLE/NCYCLEMAX,CPU_ELAPSED

end if




!!- print the initial particle position
!if(SOLVE_PART) then
!  NOUT2 = 0
!  call SAVE_PARTICLE(NOUT2)
!end if

CONT     = .true.
CONT_CPU = .true.

NOUT2 = 1
NOUT3 = 1
XFLAG = 0


if(FOUT3>0.and.MYID==0) then
 call system('rm -rf postprocessing')
 call system('mkdir postprocessing')
end if


!== Time loop ==!== Time loop ==!== Time loop ==!== Time loop ==!
!== Time loop ==!== Time loop ==!== Time loop ==!== Time loop ==!
!== Time loop ==!== Time loop ==!== Time loop ==!== Time loop ==!
!== Time loop ==!== Time loop ==!== Time loop ==!== Time loop ==!
!== Time loop ==!== Time loop ==!== Time loop ==!== Time loop ==!
!== Time loop ==!== Time loop ==!== Time loop ==!== Time loop ==!

!- Synchronize all the process
call MPI_BARRIER(MPI_COMM_WORLD,IERR)


if(MYID==0)write(*,*)
if(MYID==0)write(*,*) 'Start the time-loop !!'
if(MYID==0)write(*,*)


do while(CONT)


!!- CPU check
if(MYID == 0) TIME_START = MPI_WTIME()


!!- Print particles for post-processing
if((FOUT3>0).and.(NCYCLEMAX>=FOUT3).and.(mod(NCYCLE,NCYCLEMAX/FOUT3)==0)) then

  call PRINT_PARTICLE(TIME,NCYCLE)
  !NFILEOUT = NFILEOUT + 1

end if


!!====================================================================
!! 2. STATISTICS
!!====================================================================
!! The statistic are compute at the begining of the time loop
!! because all variables are at the same cycle.
!!--------------------------------------------------------------------
if(mod(NCYCLE,FOUT2) == 0) then

  !!--------------------------------------------------------------------
  !! 2.1 Fluid Statistics 
  !!--------------------------------------------------------------------
  if(LEVEL_STFLU>0) call STAT_FLUID(NCYCLE,TIME)

  !!--------------------------------------------------------------------
  !! 2.2 Scalar Statistics 
  !!--------------------------------------------------------------------
  !! if(LEVEL0_STSCL) call STAT_SCALAR(NCYCLE,TIME)

  !!--------------------------------------------------------------------
  !! 2.3 Particle Statistics 
  !!--------------------------------------------------------------------
  !!- Here all Lagrangian variables are at n
  if(LEVEL0_STPAR) call STAT_PARTICLE(NCYCLE,TIME)

  if(LEVEL0_STPAR .and. FROZEN_FLOW > 0) call STAT_PARTFLUID(NCYCLE, TIME)

  !!- Event count for time-averaged statistics
  if((LEVEL_STFLU>0.or.LEVEL0_STSCL.or.LEVEL0_STPAR).and.STAT_TIME)  NEVEN = NEVEN + 1


end if

!!====================================================================
!! 4. Discrete Particle Simulation
!!====================================================================
if(SOLVE_PART) then

  !! Channel Statistics and Event Count
  !! Statistics called after Collisional Phase
  if((LEVELX_STPAR .and. STAT_TIME) .and. (mod(NCYCLE,FOUT2) == 0)) then

    call  STAT_PART_CHANNEL(NCYCLE, TIME)
    NEVEN_CHAN = NEVEN_CHAN + 1
    
    call STAT_COLL_PART_CHAN(NCYCLE, TIME)
    NEVEN_COLL_CHAN = NEVEN_COLL_CHAN + 1
    
    !call STAT_PART_PDF_CHANNEL(NCYCLE, TIME)
    !NEVEN_PDF = NEVEN_PDF + 1

  end if   

!!--------------------------------------------------------------------
!! 4.1. Particle tracking
!!--------------------------------------------------------------------
    !call ADV_PARTICLE(NCYCLE)
    call ADV_PARTICLE_ELLP(NCYCLE)


!! Channel Statistics and Event Count
!! Statistics called after Particle Advance and Rebound Phase
!! Before Collision 
  if((LEVELX_STPAR .and. STAT_TIME) .and. (mod(NCYCLE,FOUT2) == 0)) then

    call  STAT_PART_CHANNEL(NCYCLE, TIME)
      
    NEVEN_CHAN = NEVEN_CHAN + 1

  end if   

!!--------------------------------------------------------------------
!! 4.2. Particle-particle collisions
!!--------------------------------------------------------------------
  call COLLISION(NCYCLE,TIME)


end if


!!====================================================================
!! 3. Fluid flow prediction
!!====================================================================

!!--------------------------------------------------------------------
!! 3.1 Direct numerical simulation
!!--------------------------------------------------------------------
if(SOLVE_FLUID == 1) then


  call ADV_FLUID(NCYCLE)
 
!- Update Ghost cell
  if(NGHTCELL>0) then

    call FLUIDCOMM(UFLU)
    call FLUIDCOMM(VFLU)
    call FLUIDCOMM(WFLU)

  end if


!!--------------------------------------------------------------------
!! 3.2 Solve Stokes equation --> Blaise your job !
!!--------------------------------------------------------------------
elseif(SOLVE_FLUID == 2) then


!!--------------------------------------------------------------------
!! 3.3 Channel flow --> Pascal your job !
!!--------------------------------------------------------------------
elseif(SOLVE_FLUID == 3) then

end if







!!====================================================================
!! X. Manage time stepping
!!====================================================================
!! Two possibilities:
!! + First, the run reaches the user defined maximum step number
!! + The elapsed time reaches the limit imposed on the computer.
!!--------------------------------------------------------------------
NCYCLE = NCYCLE + 1
TIME = TIME + DTIME


!!- Check CFL limitation
if (SOLVE_FLUID > 0 ) call CHECK_CFL(NCYCLE)


if(NCYCLE == NCYCLEMAX) CONT = .false.


!!- Check Wall time 
call TREMAIN(CONT_CPU)
if(CONT_CPU) CONT = .false.
!!--------------------------------------------------------------------





!- CPU check
TIME_END = MPI_WTIME()
!!- Full CPU time elapsed
CPU_ELAPSED = CPU_ELAPSED + TIME_END - TIME_START

!!- Averaged CPU cycle time
CPU_CYCLE = CPU_CYCLE + TIME_END - TIME_START




!!- Print percentage of simulation accomplished
if((NCYCLEMAX>=FOUT0).and.(mod(NCYCLE,NCYCLEMAX/FOUT0)==0).and.(MYID==0)) then

  write(*,10700)100.*NCYCLE/NCYCLEMAX,CPU_ELAPSED

end if



end do




if(MYID==0)write(*,*)
if(MYID==0)write(*,*) ' Time loop ended !!!'
 


!!====================================================================
!! 5. STATISTICS FINALIZING
!!====================================================================
!- Synchronize all the process
!!call MPI_BARRIER(MPI_COMM_WORLD,IERR)


!!--------------------------------------------------------------------
!! 5.1 Last cycle
!!--------------------------------------------------------------------
if(LEVEL_STFLU>0) call STAT_FLUID(NCYCLEMAX,TIME)


!! if(LEVEL0_STSCL) call STAT_SCALAR(NCYCLEMAX,TIME)

if(LEVEL0_STPAR .and. (mod(NCYCLE, FOUT2) == 0) .and. STAT_TIME) call STAT_PARTICLE(NCYCLEMAX,TIME)


!!- Event count for time-averaged statistics
if((LEVEL_STFLU>0.or.LEVEL0_STSCL.or.LEVEL0_STPAR).and.STAT_TIME)  NEVEN = NEVEN + 1


! if(LEVEL_STFLU>0.or.LEVEL0_STPAR.or.LEVEL0_STSCL) then
!   if(MYID==0) call PRINT_LASTCYCLE(NCYCLE)
! end if

if(LEVELX_STPAR .and. STAT_TIME) then
  call STAT_PART_CHANNEL(NCYCLEMAX, TIME)
  NEVEN_CHAN = NEVEN_CHAN + 1
end if
!!--------------------------------------------------------------------
!! 5.2. Print time-averaged statistics
!!--------------------------------------------------------------------

!!--------------------------------------------------------------------
!! 5.2.1 Print time-averaged statistics on time step
!!--------------------------------------------------------------------
if(MYID == 0 .and. SOLVE_FLUID > 0) then

  MEAN_DT = MEAN_DT / NEVEN_DT

  IFLAG1=1
  write(UNIT_INFO(IFLAG1),*)
  write(UNIT_INFO(IFLAG1),*)'====================================================================='
  write(UNIT_INFO(IFLAG1),*)'TIME STEP STATISTIC'
  write(UNIT_INFO(IFLAG1),*)'====================================================================='
  write(UNIT_INFO(IFLAG1),*)
  write(UNIT_INFO(IFLAG1),10605)'<dt> = ',MEAN_DT,' s'
  write(UNIT_INFO(IFLAG1),10605)'min(dt) = ',MIN_DT,' s'
  write(UNIT_INFO(IFLAG1),10605)'max(dt) = ',MAX_DT,' s'

end if

!!--------------------------------------------------------------------
!! 5.2.2 Print time-averaged statistics
!!--------------------------------------------------------------------
if(STAT_TIME) call PRINT_TIMESTAT(NCYCLE)


if(LEVEL2_STPAR) call PRINT_LAGFUNCTION

if(LEVELX_STPAR .and. STAT_TIME) call PRINT_TIMESTAT_CHANNEL
!!--------------------------------------------------------------------
!! 5.3. Compute and print last spectrum
!!--------------------------------------------------------------------
IFLAG1 = 99
if (SOLVE_FLUID == 1) call SPEC3D(IFLAG1)


!IFLAG1 = 99
!if(SOLVE_SCALAR) call SPEC3D_SCALAR(IFLAG1)

!!--------------------------------------------------------------------
!! 5.5. File closing
!!--------------------------------------------------------------------
LFLAG1 = .false.
call OPENCLOSE(LFLAG1)


!!====================================================================
!! 6. SAVE SOLUTION FOR RESTART
!!====================================================================

!!--------------------------------------------------------------------
!! 6.1. Fluid velocity field 
!!--------------------------------------------------------------------
!!- Print last fluid solution for restart
IDUMMY = -99
if(SOLVE_FLUID>0 .or. (FROZEN_FLOW>0)) call SAVE_FLUID(IDUMMY)


if(STEADY) call SAVE_FORCING
if(STEADY) call SAVE_FORCING_NEW


!!--------------------------------------------------------------------
!! 6.2. Scalar field
!!--------------------------------------------------------------------
!!- Print last scalar solution for restart
!IDUMMY = -99
!!ISAVEFLUID = 1
!if(SOLVE_SCALAR) call SAVE_SCALAR(IDUMMY)
!!ISAVEFLUID = 4



!!--------------------------------------------------------------------
!! 6.3. Save particle position and velocity
!!--------------------------------------------------------------------
!!- Print final particle solution for restart
IDUMMY = -99
if(SOLVE_PART) call SAVE_PARTICLE(IDUMMY)



!!--------------------------------------------------------------------
!! 6.4. Particle velocities and orientations
!!--------------------------------------------------------------------
if(SOLVE_PART) then
do J=1, NIG
 if(MYID==0) write(*,10702)J,NBR_EXCHANGE(J),NPMAX_CPU(J)
end do

!!- check the number of particles
if(SOLVE_PART .and. DEBUG) then
 do J = 1, NIG
  call ISUMCPU(NPART_LOC(J),IDUMMY)
  if(MYID==0) write(*,10701)J,IDUMMY
 end do
end if

end if !- If:(SOLVE_PART) 




if(MYID==0) then
write(*,*)
write(*,*) '**************************************'
write(*,*) '       END COMPUTATION FOR NCYCLE'
write(*,*) '**************************************'
write(*,*) '           Cycle =', NCYCLE
write(*,*) '      Ending time=', TIME
write(*,*) '**************************************'
write(*,*)
end if



!!- runtime CPU control
if(MYID==0) call CPUTIME_INFO(NCYCLE)



!- Synchronize all the process
!call MPI_BARRIER(MPI_COMM_WORLD,IERR)



if(MYID==0) then
 close(UNIT_INFO(1))
 close(UNIT_INFO(2))
 close(UNIT_INFO(3))
 close(UNIT_INFO(4))
end if



!- Clean P3DFFT
 call P3DFFT_CLEAN

!- Free MPI environement
 call MPI_FINALIZE (ierr)


!!- runtime CPU control
if(MYID==0) write(*,*)' This is the END ... '

!!----------------------------------------------------------------------

10601 format (2x,A,  1x,E13.6 )
10602 format (2x,A,2(1x,E13.6))
10603 format (2x,A,3(3x,E13.6))
10604 format (6(A,E13.6,2x))
10605 format (2x,A,E13.6,A)

10000 format (30(e17.7))

!10700 format (2x,' Computation at ',f6.2,' %')
10700 format (2x,' Computation at ',f6.2,' %, Elapsed time:',f12.3,' s')
10701 format (2x,' Class ',i3,' has ',i7)
10702 format (2x,' C',i2.2,' has exchanged ',i5,' part. and occupied ',i6,' in one cpu')

end program gaspart

