!!=====================================================================
!!
!!
!!  Print particle position and velocity
!!
!!
!!=====================================================================

subroutine PRINT_PARTICLE(TIME, NCYCLE)

!!- Commenting -PF - 29/02/2012
!!subroutine PRINT_PARTICLE(TIME,NFILEOUT)

!!=====================================================================
!!
!!
!!=====================================================================

use PARTICLE_PARALLEL
use PARAM_PHYS
use mod_quaternion
use DNS_DIM

implicit none

!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Time
integer, intent(in) :: NCYCLE
real(kind=8), intent(in) :: TIME

!- File name
character(len=50) :: FILENAME, FILENAME2
character(len=8) :: NUMFILE

real(kind=8), dimension(ndim, 1) :: orientation_vec

!- temporary array
real(kind=8), dimension(NPART_FULL,14) :: TMPVAR

!- number of particles of class J for all procs
integer, dimension(NPROC,NIG) :: NPART_LOC_ALL

!Various variables for ALL_GATHERV
integer, dimension(NPROC) :: COUNTS
integer, dimension(NPROC) :: DISPLS
integer                   :: BLOCK


!- Index
integer :: I, J, K, P, NP
!---------------------------------------------------------------------

!do J = 1, NIG
!
!
!-Print filename
!write(FILENAME1,10102)'postprocessing/part_',J,'_t_',NCYCLE,'.dat'
!write(FILENAME2,10102)'post_processing/part_quat_', J,'_t_',NCYCLE,'.dat'
!write(NUMFILE, '(i8)') NCYCLE
!
!
!- ASCII
!open(unit=121,file=trim(FILENAME1))
!open(unit=141,file=trim(FILENAME2))
!
!- Ecriture ASCII
!write(121,2000)
!write(121,2001)NPART_LOC(J),1,1,TIME

!do I = 1,NPART_LOC(J)
!
!
!    orientation_vec(1,1) = 0.0
!    orientation_vec(2,1) = 0.0
!    orientation_vec(3,1) = 1.0
!
!
!    call transform_basis(orientation_vec, PART(I,J)%ELLQUAT, shape(orientation_vec))
!
!    write(121,'(20(e17.7))')PART(I,J)%XP, PART(I,J)%YP, PART(I,J)%ZP, &
!                           PART(I,J)%UP, PART(I,J)%VP, PART(I,J)%WP, &
!                           orientation_vec(1,1), orientation_vec(2,1), orientation_vec(3,1)
!
!                        PART(I,J)%UFAP, PART(I,J)%VFAP, PART(I,J)%WFAP, &
!                        REAL(PART(I,J)%ID), REAL(PART(I,J)%PROC_ID)  
!
!
!end do
!
!
!write(141, 2004) TIME, NPART_LOC(J)
!
!do I = 1, NPART_LOC(J)
!
!    write(141, '(10(e17.7))') PART(I,J)%ELLQUAT%a, &
!                              PART(I,J)%ELLQUAT%b, &
!                              PART(I,J)%ELLQUAT%c, &
!                              PART(I,J)%ELLQUAT%d
!    
!end do
!
!
!- close file
!close(121)
!close(141)
!
!end do


!if(MYID==0) write(*,*) 'Particle position and velocity Binary file dropped --> ', NCYCLE, ' cycle'


do J = 1, NIG


    call MPI_ALLGATHER(NPART_LOC(J), 1, MPI_INT, NPART_LOC_ALL(:,J), 1, MPI_INT, MPI_COMM_WORLD, IERR)

    if (sum(NPART_LOC_ALL(:,J)) /= NPART_FULL) stop 'MPI_ALLGATHER is wrong!!!'

    COUNTS = [ (NPART_LOC_ALL(P,J), P = 1, NPROC) ]
    DISPLS = [ (sum(NPART_LOC_ALL(:P-1,J)), P=1, NPROC) ]

    call MPI_TYPE_VECTOR(NPART_LOC(J), 1, 1, MPI_DOUBLE_PRECISION, BLOCK, IERR)
    call MPI_TYPE_COMMIT(BLOCK,IERR)     

    !- Drop x-position
    call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%XP       , 1, BLOCK, TMPVAR(:,1),  COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
    call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%YP       , 1, BLOCK, TMPVAR(:,2),  COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
    call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%ZP       , 1, BLOCK, TMPVAR(:,3),  COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
    call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%UP       , 1, BLOCK, TMPVAR(:,4),  COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
    call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%VP       , 1, BLOCK, TMPVAR(:,5),  COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
    call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%WP       , 1, BLOCK, TMPVAR(:,6),  COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
    call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%ELLQUAT%a, 1, BLOCK, TMPVAR(:,7),  COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
    call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%ELLQUAT%b, 1, BLOCK, TMPVAR(:,8),  COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
    call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%ELLQUAT%c, 1, BLOCK, TMPVAR(:,9),  COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
    call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%ELLQUAT%d, 1, BLOCK, TMPVAR(:,10), COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
    call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%OMEGAX   , 1, BLOCK, TMPVAR(:,11), COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
    call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%OMEGAY   , 1, BLOCK, TMPVAR(:,12), COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
    call MPI_ALLGATHERV(PART(1:NPART_LOC(J),J)%OMEGAZ   , 1, BLOCK, TMPVAR(:,13), COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)
    call MPI_ALLGATHERV(real(PART(1:NPART_LOC(J),J)%IDP), 1, BLOCK, TMPVAR(:,14), COUNTS, DISPLS, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, IERR)

    !- Define file name
    if(NCYCLE>=0) write(FILENAME,10102)'postprocessing/part_p',J,'_t',NCYCLE,'.bin'

    !if(MYID==0)write(*,*) '   + file:',trim(FILENAME)

    open(unit=150,file=trim(FILENAME),status='replace',form='unformatted')


    write(150)(TMPVAR(I,1) ,I=1,NPART_FULL)
    write(150)(TMPVAR(I,2) ,I=1,NPART_FULL)
    write(150)(TMPVAR(I,3) ,I=1,NPART_FULL)
    write(150)(TMPVAR(I,4) ,I=1,NPART_FULL)
    write(150)(TMPVAR(I,5) ,I=1,NPART_FULL)
    write(150)(TMPVAR(I,6) ,I=1,NPART_FULL)
    write(150)(TMPVAR(I,7) ,I=1,NPART_FULL)
    write(150)(TMPVAR(I,8) ,I=1,NPART_FULL)
    write(150)(TMPVAR(I,9) ,I=1,NPART_FULL)
    write(150)(TMPVAR(I,10),I=1,NPART_FULL)
    write(150)(TMPVAR(I,11),I=1,NPART_FULL)
    write(150)(TMPVAR(I,12),I=1,NPART_FULL)
    write(150)(TMPVAR(I,13),I=1,NPART_FULL)
    write(150)(TMPVAR(I,14),I=1,NPART_FULL)
    close(150)

end do !- end do IG = 1, NIG


!===========================================================================
2000 format ('#VARIABLES = "x", "y", "z", "u", "v", "w"')
2001 format ('#ZONE F=POINT I=',i6,' J=',i4,' K=',i4,' SOLUTIONTIME=',E13.6)

10101 format (A,I2.2,A,I2.2,A,I2.2,A)
10205 format(A,I4.4)
2004 format ('# Time =', e17.7,'Particles=',i7)
10102 format (A,I2.2,A,I8.8,A)

end subroutine PRINT_PARTICLE
