!!====================================================================
!!
!!
!!====================================================================

subroutine INITIATION_FROZEN_FLOW

!!====================================================================
!!
!!
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE
use FLUID_VARIABLE
use PARAM_PHYS
use MPI_STRUCTURES

implicit none

!---------------------------------------------------------------------
! 	VARIABLE INITIALIZATION
!---------------------------------------------------------------------

!- Index
integer :: I,J,K, N !flag2, COUNT

!- Parameters for Velocity profile
real(kind=8) :: dy_frozen, dz_frozen
character(len=40) :: FILENAME


!logical :: flag

real(kind=8), dimension(21) :: x_in, y_in, z_in, Ug_in, Vorty_in 
real(kind=8), dimension(21,21,21) :: UTMP, VTMP, WTMP
real(kind=8), dimension(21,21,21) :: VORTX_TMP, VORTY_TMP, VORTZ_TMP 
real(kind=8), dimension(ISTART(1):IEND(1)) :: WINT
real(kind=8), dimension(ISTART(2):IEND(2)) :: VORTY_INT
!---------------------------------------------------------------------

VORTICITY_X(:,:,:) = 0.0
VORTICITY_Y(:,:,:) = 0.0
VORTICITY_Z(:,:,:) = 0.0

! COUNT= 1
! flag = .true.
! flag2 = 2

! Case of Zero Velocity
if (FROZEN_FLOW == 1) then


    UFLU(:,:,:) = 0.0 
    VFLU(:,:,:) = 0.0
    WFLU(:,:,:) = 0.0

    if(MYID==0) write(*,*)'Frozen Fluid Velocity Initiation: Imposed Uniform Flow Uf=0'


else if (FROZEN_FLOW == 2) then


    UFLU(:,:,:) = 0.0
    VFLU(:,:,:) = 0.0

    do I = ISTART(1), IEND(1)

        WFLU(I,:,:) = 2.0*(DX*(I-1))/LXMAX - 1.0

    end do


    if(MYID==0) write(*,*)'Frozen Fluid Velocity Initiation: Imposed Shear Flow'


else if (FROZEN_FLOW == 3) then


    N = 21

    FILENAME = 'Ug_mean.in'
    open(unit=310, file=trim(FILENAME), action='read')

    FILENAME = 'vorticity_y.txt'
    open(unit=290, file=trim(FILENAME), action='read')


    ! Read the data from file
    do I = 1, N 
        read(310, *) x_in(I), Ug_in(I)
        read(290, *) Vorty_in(I)
    end do

    !dx = LXMAX / real(N - 1)
    dy_frozen = LYMAX / real(N - 1)
    dz_frozen = LZMAX / real(N - 1)

    ! Initialize the frozen fluid velocity
    do I = 1, N

       !x_in(I) = dx*(I-1)    
        y_in(I) = dy_frozen*(I-1)
        z_in(I) = dz_frozen*(I-1)

    end do


    do K = 1, N
	   do J = 1, N
		  do I = 1, N

            WTMP(I,J,K) = Ug_in(I)
            VORTY_TMP(I,J,K) = Vorty_in(I)

		  end do 
	   end do 
    end do




    call FROZEN_INTERP_LAG3(x_in, y_in, z_in     &
                            ,N, N, N             &
                            ,WTMP                &
                            ,NX                  &
                            ,XMESH,YMESH,ZMESH   &
                            ,WINT)


    call FROZEN_INTERP_LAG3(x_in, y_in, z_in     &
                            ,N, N, N             &
                            ,VORTY_TMP           &
                            ,NY                  &
                            ,XMESH, YMESH, ZMESH &
                            ,VORTY_INT)



    do I = ISTART(1), IEND(1)

        UFLU(I,:,:) = ZERO
        VFLU(I,:,:) = ZERO
        WFLU(I,:,:) = WINT(I)

    end do


    do J = ISTART(2), IEND(2)

        VORTICITY_X(:,J,:) = ZERO
        VORTICITY_Y(:,J,:) = VORTY_INT(J)
        VORTICITY_Z(:,J,:) = ZERO

    end do

    if(MYID == 0) write(*,*) 'Frozen Fluid Velocity Initiation --> Imposed Channel flow Profile'
    if(MYID == 0) write(*,*) 'Frozen Vorticity Field Initiation --> Imposed OK'


end if  ! Type of Frozen Flow


!write(*,*) DX, DY, DZ


!! Calculation of vorticity
!do I = ISTART(1)+1, IEND(1)-1

!    do J = ISTART(2)+1, IEND(2)-1

!        do K = ISTART(3)+1, IEND(3)-1


!            VORTICITY_X(I,J,K) = (WFLU(I,J+1,K) - WFLU(I,J-1,K))/(2*DY) - (VFLU(I,J,K+1) - VFLU(I,J,K-1))/(2*DZ)

!            VORTICITY_Y(I,J,K) = (UFLU(I,J,K+1) - UFLU(I,J,K-1))/(2*DZ) - (WFLU(I+1,J,K) - WFLU(I-1,J,K))/(2*DX)

!            VORTICITY_Z(I,J,K) = (VFLU(I+1,J,K) - VFLU(I-1,J,K))/(2*DX) - (UFLU(I,J+1,K) - UFLU(I,J-1,K))/(2*DY)


!        end do

!    end do

!end do


 !FILENAME = 'fluid_vel_profile.dat'
 !open(unit=310, file=trim(FILENAME), status='replace')
close(310)
close(290)

 !do I = ISTART(1), IEND(1)
 !    if(MYID==0) write(310,10000) VORTICITY_X(I,64,64),VORTICITY_X(I,64,1),VORTICITY_X(I,64,32)
 !end do

10000 format (40(e17.7))

end subroutine INITIATION_FROZEN_FLOW