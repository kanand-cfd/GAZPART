!!====================================================================
!!
!!          Particle positions initiation
!!
!!====================================================================

subroutine INITIATION_PARTICLE_POSITION

!!====================================================================
!! Here the particle position and velocity are initiated according
!! to the variable INIT_PART_POSITION specified by user in the
!! parameter file: 'param.in'.
!! Note that the fluid velocity at the particle position is computed
!! at the end of the subroutine because it is needed by time-advancing
!! numerical scheme.
!!====================================================================
!! Particle position initiation: 
!!------------------------------
!! INIT_PART_POSITION=
!!
!!  0: The particle positions are spefied in the Fortran file.
!!     The default position correspond of one particle in a cell of
!!     fluid mesh with a small shift from the cell center. This position
!!     distribution allows to check the accuracy of the interpolation
!!     scheme of fluid velocity at the particle position.
!!
!!  1: Random particle distribution in the box. Take care that a
!!     factor 0.99 im used in order to ensure that the particles are
!!     inside the box.
!!
!!  2: Particle positions are read from the binary file 'POSPART.ini'
!!
!!  3: Particle injection at one edge. This is like an injector.
!!
!!--------------------
!!
!! Warning: The particle number in each CPU is defined with the
!!          particle position
!!
!!--------------------
!!
!! Warning: Particle velocities and positions are read simultaneously
!!
!!====================================================================

use DNS_DIM

use GEOMETRIC_VARIABLE
use PARAM_PHYS
use PARTICLE_PARALLEL
use FLUID_VARIABLE
use MPI_STRUCTURES
use mod_quaternion

implicit none


!!====================================================================
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!---------------------------------------------------------------
!- Read data variables
!---------------------------------------------------------------
!- File name 
character(len=40) :: FILENAME

integer :: IDUMMY

real(kind=8) :: RDUMMY

integer :: RECSIZE, SIZE_INT, SIZE_REAL, NVARIABLE

!- 
integer :: NP_READ, ND_READ

!---------------------------------------------------------------
!- Variable for random and uniform initiation
!---------------------------------------------------------------
real(kind=8) :: X0, Y0, Z0
real(kind=8) :: U0, V0, W0
real(kind=8) :: XMAX, YMAX, ZMAX

!- 
real(kind=8) :: XRAND, YRAND, ZRAND, GAUSS
integer :: ID

real(kind=8) :: DD

real(kind=8) :: u, x, y, z, alpha2, alpha1

real(kind=8), dimension(3,1) :: principal_axis

!---------------------------------------------------------------
!- Injection parameter
!---------------------------------------------------------------
!- Position of injection center
real(kind=8) :: XINJ, YINJ 

!- Radius and width of injection
real(kind=8) :: RINJ, DRINJ

!- Angle and velocity of injection
real(kind=8) :: PHINJ, UINJ

!- 
real(kind=8) :: RADIUS, THETA

real(kind=8) :: NPART_LOC_UNIF

real(kind=8) :: DPMIN, DPMAX, DDP

integer, dimension(POLYDISPMAX) :: NPCLASS_LOC

! real(kind=8), dimension(NIG) :: NQ_NP   !- Nq/Np
! real(kind=8), dimension(NIG) :: DQ_DP   !- dq/dp
! real(kind=8), dimension(NIG) :: RHOQ_RHOP   !- rhoq/rhop

real(kind=8), dimension(:,:), allocatable :: TMPVAR

!- Shift
integer :: SHIFT

!- Index
integer :: I, J, K, L, M, N, INP_NQ, IDP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!- User -!- User -!- User -!- User -!- User
!!- User -!- User -!- User -!- User -!- User
!!- User -!- User -!- User -!- User -!- User
!NQ_NP(:) = 1.00


! NQ_NP(3) = 0.25
! NQ_NP(4) = 0.50
! NQ_NP(5) = 0.75
! NQ_NP(6) = 1.00
! NQ_NP(7) = 1.25
! NQ_NP(8) = 2.00

! DQ_DP(:) = 1.0
! RHOQ_RHOP(:) = 2.5

!!- User -!- User -!- User -!- User -!- User
!!- User -!- User -!- User -!- User -!- User
!!- User -!- User -!- User -!- User -!- User


!!- Number of particle in each CPU
do J = 1, NIG
    NPART_LOC(J) = int(NPART_FULL/NPROC)
end do


!!- Variable number of particle per proc
if(INIT_PART_POSITION==2) then

    NPART_LOC_UNIF = int(NPART_FULL/NPROC)

    if(MYID<NPROC-1) then
        do J = 1, NIG
            NPART_LOC(J) = NPART_LOC_UNIF
        end do
    else !- Last CPU get the remaining particles
        do J = 1, NIG

            NPART_LOC(J) = NPART_LOC_UNIF + (NPART_FULL - NPART_LOC_UNIF*NPROC)

        end do

    end if

    !!=======================================================================
    !!- Initiation of particle diameter (parameters read in param.in)
    !!=======================================================================
    !! FIX THIS PART (KARAN)
    do J = 1, NIG

        if(POLYDISP(J)>1) then

            INP_NQ = int(NPART_LOC(J)/(1+NQ_NP(J)))

            PART(1:INP_NQ,J)%IDP = 1

            PART(INP_NQ+1:NPART_LOC(J),J)%IDP = 2

            if(PARTDEF(J) ==0 .or. PARTDEF(J)==1) then

                RHOP(J,:) = RHOP_USER(J)
                !DPART(J,:) = DPART_USER(J)
                EMAJ_PART(J,:) = EMAJ_PART_USER(J)

            else

                if(POLYDISP(J) ==1) then

                    RHOP(J,1) = RHOP_USER(J)
                    !DPART(J,1) = DPART_USER(J)
                    EMAJ_PART(J,1) = EMAJ_PART_USER(J)

                else

                    RHOP(J,1) = RHOP_USER(J)
                    RHOP(J,2) = RHOP_USER(J)/RHOQ_RHOP(J)

                    !DPART(J,1) = DPART_USER(J)
                    !DPART(J,2) = DPART_USER(J)/DQ_DP(J)

                    EMAJ_PART(J,1) = EMAJ_PART_USER(J)
                    EMAJ_PART(J,2) = EMAJ_PART_USER(J)


                    if(MYID==0) write(*,'(A,I2.2,A,F8.3,A,F8.3)')' -> C',J, &
                    ' rho_p/rho_q=',RHOP(J,1)/RHOP(J,2),' dp/dq=',EMAJ_PART(J,1)/EMAJ_PART(J,2)

                end if

            end if

        else

            PART(1:NPART_LOC(J),J)%IDP = 1

        end if

    end do

end if



!!====================================================================
!! 1. Specific position given by user
!!====================================================================
if(INIT_PART_POSITION == 0 .or. INIT_PART_POSITION == 1) then

!- Shift 
! SHIFT = 1: 1 cell  between two particles
!       = 2: 2 cells between two particles 
    SHIFT = 2

    X0 = XMESH(ISTART(1)) + DX*0.5D0
    Y0 = YMESH(ISTART(2)) + DY*0.5D0
    Z0 = ZMESH(ISTART(3)) + DZ*0.5D0

    L = 0
    M = 0
    N = 0

    do J = 1, NIG

        do I = 1, NPART_LOC(J)

            PART(I,J)%XP = X0 + L*DX
            
            PART(I,J)%YP = Y0 + M*DY
            
            PART(I,J)%ZP = Z0 + N*DZ

            L = L + SHIFT

            if((PART(I,J)%XP + DX) > XMESH(IEND(1))) then

                L = 0
                M = M + SHIFT

                if ((PART(I,J)%YP + DY) > YMESH(IEND(2))) then

                    M = 0

                    N = N + SHIFT

                    if ((PART(I,J)%ZP + DZ) > ZMESH(IEND(3))) then

                        N = 0

                    end if

                end if

            end if

            !! Uniform initiation for Particle Orientation
            PART(I,J)%ELLQUAT%a = 1.0
            PART(I,J)%ELLQUAT%b = 0.0
            PART(I,J)%ELLQUAT%c = 0.0
            PART(I,J)%ELLQUAT%d = 0.0


        end do


    end do


    if(MYID==0) write(*,*)'Particle position initiation --> Uniform'

!!====================================================================
!! 2. Random position
!!====================================================================
elseif(INIT_PART_POSITION == 2) then

    do J = 1, NIG

        do I = 1, NPART_LOC(J)

            IDP = PART(I,J)%IDP

            call random_number(XRAND)
            PART(I,J)%XP = XRAND*(LXMAX - 2.0*EMAJ_PART(J,IDP)) + EMAJ_PART(J,IDP) !XMESH(ISTART(1)) + XRAND*(XMESH(IEND(1))+DX-XMESH(ISTART(1)))

            call random_number(XRAND)
            PART(I,J)%YP = XRAND*(LYMAX - 2.0*EMAJ_PART(J,IDP)) + EMAJ_PART(J,IDP)!YMESH(ISTART(2)) + XRAND*(YMESH(IEND(2))+DY-YMESH(ISTART(2)))

            call random_number(XRAND)
            PART(I,J)%ZP = XRAND*(LZMAX - 2.0*EMAJ_PART(J,IDP)) + EMAJ_PART(J,IDP)!ZMESH(ISTART(3)) + XRAND*(ZMESH(IEND(3))+DZ-ZMESH(ISTART(3)))

   
            !! Random initialisation of the quaternion
            theta = 0.0*(ppi/180.0) !ppi/2.0

            call random_number(XRAND)
            PART(I,J)%ELLQUAT%a = 2.0*XRAND - 1.0 ! cos(theta/2.0)
 
            call random_number(XRAND)
            PART(I,J)%ELLQUAT%b = 2.0*XRAND - 1.0 ! sin(theta) 0.0 !

            call random_number(XRAND)
            PART(I,J)%ELLQUAT%c = 2.0*XRAND - 1.0 !0.0 sin(theta/2.0) !

            call random_number(XRAND)
            PART(I,J)%ELLQUAT%d = 2.0*XRAND - 1.0 !0.0 0.0 !

            ! Normalize the Quaternion
            PART(I,J)%ELLQUAT = unit_quat(PART(I,J)%ELLQUAT)

            principal_axis(:,:) = 0.0; principal_axis(1,1) = 1.0

            !call transform_basis(principal_axis, PART(I,J)%ELLQUAT, shape(principal_axis))
            !write(*,*) 'Transformed p. axis = ', principal_axis

            PART(I,J)%COLOR = MYID
            PART(I,J)%COLOR = 1.0

        end do

    end do

    if(MYID==0) write(*,*)'Particle position initiation: Uniformly randomized --> OK'


!!====================================================================
!! 3. Stored particle position 
!!====================================================================
elseif(INIT_PART_POSITION == 3) then

    !!--------------------------------------------------------------------
    !! 3.1 Multiple binary files
    !!--------------------------------------------------------------------
    if(ISAVEPART == 1) then


        if(MYID==0) write(*,*)'Particle position initiation: Read from Multiple Binary Files'


        !- Define file name
        write(FILENAME,10303)'PART',trim(FILE_EXT),'.ini'

        !- Open file containing the last particle position and velocity
        open(unit=150, file=trim(FILENAME), status='old', form='unformatted')


        do J = 1, NIG

            read(150)NP_READ

            if(NP_READ>NPMAX_LOC) then

                write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                write(*,*)'!!                     ERROR                    !!'
                write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                write(*,*)'!!'
                write(*,*)'!! file     : initiation_particle_position.f90'
                write(*,*)'!!'
                write(*,*)'!!   NP_READ=',NP_READ
                write(*,*)'!! NPMAX_LOC=', NPMAX_LOC
                write(*,*)'!!'
                write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                stop

            end if  

            !!write(*,*)'ID:',MYID,' Class:',J,'  NPREAD=',NP_READ

            NPART_LOC(J) = NP_READ

            !- Read block data 
            read(150)(PART(I,J)%IDP,I=1,NP_READ)
            read(150)(RDUMMY,I=1,NP_READ)
            read(150)(PART(I,J)%XP,I=1,NP_READ)
            read(150)(PART(I,J)%YP,I=1,NP_READ)
            read(150)(PART(I,J)%ZP,I=1,NP_READ)

            read(150)(PART(I,J)%UP,I=1,NP_READ)
            read(150)(PART(I,J)%VP,I=1,NP_READ)
            read(150)(PART(I,J)%WP,I=1,NP_READ) 
            
            read(150)(PART(I,J)%UFAP,I=1,NP_READ)
            read(150)(PART(I,J)%VFAP,I=1,NP_READ)
            read(150)(PART(I,J)%WFAP,I=1,NP_READ)

            read(150)(PART(I,J)%ELLQUAT%a, I=1, NP_READ)
            read(150)(PART(I,J)%ELLQUAT%b, I=1, NP_READ)
            read(150)(PART(I,J)%ELLQUAT%c, I=1, NP_READ)
            read(150)(PART(I,J)%ELLQUAT%d, I=1, NP_READ)

            read(150)(PART(I,J)%OMEGAX, I=1, NP_READ)
            read(150)(PART(I,J)%OMEGAY, I=1, NP_READ)
            read(150)(PART(I,J)%OMEGAZ, I=1, NP_READ)

        end do

        !- Close file
        close(150)


    !!--------------------------------------------------------------------
    !! 3.2 Single binary file
    !!--------------------------------------------------------------------
    elseif(ISAVEPART == 2) then


        if(MYID==0) write(*,*)'Particle position initiation: Read from single binary file'

        allocate(TMPVAR(NPART_FULL,14))


        do J = 1, NIG

            !- Define file name
            write(FILENAME,10100)'part_p',J,'.ini'
            open(unit=150,file=trim(FILENAME),status='old',form='unformatted')

            read(150)NPART_FULL
            read(150)NPCLASS(J,:)
            read(150)RHOP(J,:)
            !read(150)DPART(J,:)
            read(150)EMAJ_PART(J,:)
            read(150)APR_PART(J)

            read(150)(TMPVAR(I,1),I=1,NPART_FULL) !- x-position
            read(150)(TMPVAR(I,2),I=1,NPART_FULL) !- y-position
            read(150)(TMPVAR(I,3),I=1,NPART_FULL) !- z-position
            
            read(150)(TMPVAR(I,4),I=1,NPART_FULL) !- x-velocity
            read(150)(TMPVAR(I,5),I=1,NPART_FULL) !- y-velocity
            read(150)(TMPVAR(I,6),I=1,NPART_FULL) !- z-velocity

            !- Orientation
            read(150)(TMPVAR(I,7),I=1,NPART_FULL) !- Quat a
            read(150)(TMPVAR(I,8),I=1,NPART_FULL) !- Quat b
            read(150)(TMPVAR(I,9),I=1,NPART_FULL) !- Quat c
            read(150)(TMPVAR(I,10),I=1,NPART_FULL)!- Quat d

            !- Angular Velocity
            read(150)(TMPVAR(I,11),I=1,NPART_FULL)!- omega x
            read(150)(TMPVAR(I,12),I=1,NPART_FULL)!- omega y
            read(150)(TMPVAR(I,13),I=1,NPART_FULL)!- omega z

            !  read(150)(RDUMMY,I=1,NPART_FULL) !- Particle diameter
            read(150)(TMPVAR(I,14),I=1,NPART_FULL) !- Index for polydisperse

            !- Particles are distributed on each processor
            NP_READ = 0

            do I = 1, NPART_FULL
                if(     TMPVAR(I,3)<=ZMESH(IEND(3))+DZ.and.TMPVAR(I,3)>=ZMESH(ISTART(3)) &
                    .and.TMPVAR(I,2)<=YMESH(IEND(2))+DY.and.TMPVAR(I,2)>=YMESH(ISTART(2)) ) then

                    NP_READ = NP_READ +1

                    !- Load the variable is the particles is in the domain managed by the task
                    PART(NP_READ,J)%XP        = TMPVAR(I,1)
                    PART(NP_READ,J)%YP        = TMPVAR(I,2)
                    PART(NP_READ,J)%ZP        = TMPVAR(I,3)
                    
                    PART(NP_READ,J)%UP        = TMPVAR(I,4)
                    PART(NP_READ,J)%VP        = TMPVAR(I,5)
                    PART(NP_READ,J)%WP        = TMPVAR(I,6)

                    PART(NP_READ,J)%ELLQUAT%a = TMPVAR(I,7)
                    PART(NP_READ,J)%ELLQUAT%b = TMPVAR(I,8)
                    PART(NP_READ,J)%ELLQUAT%c = TMPVAR(I,9)
                    PART(NP_READ,J)%ELLQUAT%d = TMPVAR(I,10)

                    PART(NP_READ,J)%OMEGAX    = TMPVAR(I,11)
                    PART(NP_READ,J)%OMEGAY    = TMPVAR(I,12)
                    PART(NP_READ,J)%OMEGAZ    = TMPVAR(I,13)

                    PART(NP_READ,J)%IDP   = int(TMPVAR(I,14))

                end if

            end do


            close(150)


            NPART_LOC(J) = NP_READ

            if(MYID==0) write(*,*)' -> file: ',trim(FILENAME)

        end do !- end do IG = 1, NIG


        if(DEBUG) then

            ! Total particles number checking

            do J=1, NIG

                call ISUMCPU(NPART_LOC(J),IDUMMY)

                if (MYID==0) write(*,*) 'Read Part -> Class:',J,' Full number of particles ',IDUMMY

            end do

        end if


        deallocate(TMPVAR)


    !!--------------------------------------------------------------------
    !! 3.2 MPI I/O
    !!--------------------------------------------------------------------    
    elseif(ISAVEPART == 3) then


        FILENAME = 'traj.ini'


        call READ_PART_MPIIO(PART,FILENAME)

        if(MYID==0) write(*,*)'Particle position initiation: Read from file --> OK'
        if(MYID==0)write(*,*) '    + MPI I/O'
        if(MYID==0)write(*,10300) '     + File name = ',trim(FILENAME)


    end if !!- If: ISAVEPART

    !!- Switch particle velocity initiation flag
    INIT_PART_VELOCITY = 3



end if


!!- Proc ID managing the particles
!PART(:,:)%PROC_ID = MYID




!!- Boundary conditions
if (PERIODICITY) call BOUNDARY_PARTICLE


!!=========================================================
!! 5. Check particle number
!!=========================================================
if(MYID==0) write(*,*) 'Particle position initiation --> OK'

do J=1, NIG

    call ISUMCPU(NPART_LOC(J),IDUMMY)

    if (MYID==0) write(*,10102) ' -> C',J,' Total Nbr of Part. = ',IDUMMY


    if(POLYDISP(J)>1) then

        NPCLASS_LOC(:) = 0

        do I = 1, NPART_LOC(J)
            NPCLASS_LOC(PART(I,J)%IDP) = NPCLASS_LOC(PART(I,J)%IDP) + 1
        end do

        do IDP=1,POLYDISP(J)

            call ISUMCPU(NPCLASS_LOC(IDP),NPCLASS(J,IDP))

            if(MYID==0) then
            
                write(*,'(A,I2.2,A,I2.2,A,I8)')' -> C',J,'-',IDP,' Np=',NPCLASS(J,IDP)

            end if

        end do !- en do IDP = 

        if(MYID==0)write(*,'(A,I2.2,A,I8)')' -> C',J,' Np+Nq=',NPCLASS(J,1)+NPCLASS(J,2)

        if(MYID==0)write(*,'(A,I2.2,A,F6.3)')' -> C',J,' Np/Nq=',real(NPCLASS(J,1))/real(NPCLASS(J,2))

    end if

    do IDP = 1, POLYDISP(J)

        ELL_A = EMAJ_PART(J,IDP)

        lambda = APR_PART(J)

        ! Semi Minor axes
        ELL_B = ELL_A/lambda
        ELL_C = ELL_B

        ! Moment of Inertia of the Ellipsoid
        IPXX(J,IDP) = (ELL_B**2 + ELL_C**2)/5.0
        IPYY(J,IDP) = (ELL_A**2 + ELL_C**2)/5.0
        IPZZ(J,IDP) = (ELL_A**2 + ELL_B**2)/5.0

        if(MYID==0)write(*,'(A,ES14.7)') "Semi-Major Axis = ", ELL_A
        if(MYID==0)write(*,'(A,ES14.7,ES14.7)') "Semi-Minor Axis = ", ELL_B, ELL_c
        if(MYID==0)write(*,'(A,ES14.7)') "Aspect Ratio = ", lambda

        if(MYID==0)write(*,*)'Moment of Inertia'
        
        if(MYID==0)write(*,'(A,I2.2,A,I2.2,A,ES14.7)')'-> C',J,' - ',IDP, ' I_xx = ', IPXX(J,IDP)
        if(MYID==0)write(*,'(A,I2.2,A,I2.2,A,ES14.7)')'-> C',J,' - ',IDP, ' I_yy = ', IPYY(J,IDP)
        if(MYID==0)write(*,'(A,I2.2,A,I2.2,A,ES14.7)')'-> C',J,' - ',IDP, ' I_zz = ', IPZZ(J,IDP)

    end do

end do



!!--------------------------------------------------------------------
10100 format (A,I2.2,A)
10101 format (A,I2.2,A,A)
10102 format (A,I2.2,A,I8)
10300 format (A,A)
10303 format (A,A,A)


end subroutine INITIATION_PARTICLE_POSITION
