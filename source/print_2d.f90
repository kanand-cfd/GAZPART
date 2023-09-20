!!=====================================================================
!!
!!
!!  Print particle 2d Fluid and particle position and velocity
!!
!!
!!=====================================================================

subroutine PRINT_2D(TIME)

!!subroutine PRINT_PARTICLE(TIME,NFILEOUT)

!!=====================================================================
!!
!!
!!=====================================================================

use PARTICLE_PARALLEL
use FLUID_VARIABLE
use GEOMETRIC_VARIABLE
use STATISTICS
use PARAM_PHYS
use DNS_DIM
use P3DFFT
use CHECK_CPU

implicit none

!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------

!- Vorticity
double complex, &
   dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)) :: VORT_FOU

real(kind=8),   &
   dimension(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)) :: VORT_REAL


!- Time
real(kind=8), intent(in) :: TIME

!- File name
character(len=50) :: FILENAME

!- Slice position along x-direction (index)
integer :: ISLC

!- With of the slice for particle position
real(kind=8) :: DSLC, DSLC1, DSLC2


!- Index
integer :: I, J, K, P, IDP, NP_SLC, NP_SLC1
!---------------------------------------------------------------------

ISLC = int(NX-1)
DSLC = DX*0.5

if(MYID==0) write(*,*)' Drop 2d file in postprocessing'

if(SOLVE_FLUID>0) then 

 !-Print filename
 write(FILENAME,10102)'postprocessing/2dflu',trim(FILE_EXT),'_t',NFILEOUT,'.bin'

 !- open binary file
 open(unit=300,file=trim(FILENAME),form='unformatted')
 !- 
 write(300)TIME
 write(300)ISIZE(2),ISTART(2),IEND(2)
 write(300)ISIZE(3),ISTART(3),IEND(3)
 write(300)((UFLU(ISLC,J,K),J=ISTART(2),IEND(2)),K=ISTART(3),IEND(3))
 write(300)((VFLU(ISLC,J,K),J=ISTART(2),IEND(2)),K=ISTART(3),IEND(3))
 write(300)((WFLU(ISLC,J,K),J=ISTART(2),IEND(2)),K=ISTART(3),IEND(3))


 !! [vort_x]f = i*ky*[wf]fou - i*kz*[vf]fou
 do K = FSTART(3), FEND(3)
  do J = FSTART(2), FEND(2)
   do I = FSTART(1), FEND(1)
    VORT_FOU(I,J,K) = ICMPL*WFOU(I,J,K)*KY(J) - ICMPL*VFOU(I,J,K)*KZ(K)
   end do
  end do
 end do

 !- Back in physical space
 call P3DFFT_BTRAN_C2R(VORT_FOU,VORT_REAL,FFTFLAG)

 write(300)((VORT_REAL(ISLC,J,K),J=ISTART(2),IEND(2)),K=ISTART(3),IEND(3))

 !! [vort_y]f = i*kz*[uf]fou - i*kx*[wf]fou
 do K = FSTART(3), FEND(3)
  do J = FSTART(2), FEND(2)
   do I = FSTART(1), FEND(1)
    VORT_FOU(I,J,K) = ICMPL*UFOU(I,J,K)*KZ(K) - ICMPL*WFOU(I,J,K)*KX(I)
   end do
  end do
 end do

 !- Back in physical space
 call P3DFFT_BTRAN_C2R(VORT_FOU,VORT_REAL,FFTFLAG)

 write(300)((VORT_REAL(ISLC,J,K),J=ISTART(2),IEND(2)),K=ISTART(3),IEND(3))

 !! [vort_z]f = i*kx*[vf]fou - i*kx*[uf]fou
 do K = FSTART(3), FEND(3)
  do J = FSTART(2), FEND(2)
   do I = FSTART(1), FEND(1)
    VORT_FOU(I,J,K) = ICMPL*VFOU(I,J,K)*KX(I) - ICMPL*UFOU(I,J,K)*KY(J)
   end do
  end do
 end do

 !- Back in physical space
 call P3DFFT_BTRAN_C2R(VORT_FOU,VORT_REAL,FFTFLAG)

 write(300)((VORT_REAL(ISLC,J,K),J=ISTART(2),IEND(2)),K=ISTART(3),IEND(3))



 !- close file
 close(300)

end if !- end if SOLVE_FLUID



if(SOLVE_PART) then

!=====================================================================
! 1. Save 2d fluid slice
!=====================================================================
!-Print filename
!write(FILENAME,10101)'2dflu_t',NFILEOUT,trim(FILE_EXT),'.dat'
!
!!- ASCII
!open(unit=300,file=trim(FILENAME))
!
!!- Ecriture ASCII
!write(300,2000)
!write(300,2001)1,ISIZE(2),ISIZE(3),TIME
!
!do K = ISTART(3), IEND(3)
! do J = ISTART(2), IEND(2)
!   write(300,10000)XMESH(ISLC), YMESH(J), ZMESH(K), &
!                            UFLU(ISLC,J,K), &
!                            VFLU(ISLC,J,K), &
!                            WFLU(ISLC,J,K)
! end do
!end do
!
!
!!- close file
!close(300)



!=====================================================================
! 2. Save 2d particle slice slice
!=====================================================================
do J = 1, NIG

!- Number of particle in the slice
 NP_SLC = 0
 NP_SLC1 = 0
 do I = 1,NPART_LOC(J)

   ISLC = NSLICE
   IDP = PART(I,J)%IDP
   DSLC1 = DEL_X(J,IDP,ISLC)
   DSLC2 = DEL_X(J,IDP,1)

   if(PART(I,J)%XP>=XSTAT(J,IDP,ISLC+1)-DSLC1.and.PART(I,J)%XP<=XSTAT(J,IDP,ISLC+1)+DSLC1) NP_SLC = NP_SLC + 1
   if(PART(I,J)%XP>=XSTAT(J,IDP,1)-DSLC2.and.PART(I,J)%XP<=XSTAT(J,IDP,1)+DSLC2) NP_SLC1 = NP_SLC1 + 1

 end do

 !-Print filename
 write(FILENAME,10103)'postprocessing/2dpart',trim(FILE_EXT),'_p',J,'_t',NFILEOUT,'.dat'

 !- Open text file
 open(unit=300,file=trim(FILENAME))

 write(300,2004)TIME, NP_SLC 
 do I = 1,NPART_LOC(J)

   ISLC = NSLICE
   IDP = PART(I,J)%IDP
   DSLC1 = DEL_X(J,IDP,ISLC)
   DSLC2 = DEL_X(J,IDP,ISLC)

   if((PART(I,J)%XP>=XSTAT(J,IDP,ISLC+1)-DSLC1.and.PART(I,J)%XP<=XSTAT(J,IDP,ISLC+1)+DSLC1) .or. &
      (PART(I,J)%XP>=XSTAT(J,IDP,1)-DSLC2.and.PART(I,J)%XP<=XSTAT(J,IDP,1)+DSLC2)) &
         write(300,10000)PART(I,J)%XP, &
         PART(I,J)%YP, &
         PART(I,J)%ZP, &
         PART(I,J)%UP, &
         PART(I,J)%VP, &
         PART(I,J)%WP, &
         DPART(J,IDP)
 
 end do

 !- close file
 close(300)

end do !- end do J = 1, NIG

end if

!===========================================================================
2000 format ('VARIABLES = "x", "y", "z", "u", "v", "w", "omega"')
2001 format ('ZONE F=POINT I=',i6,' J=',i4,' K=',i4,' SOLUTIONTIME=',E13.6)

2004 format ('# Time =', e17.7,'Particles=',i7)

10000 format (15(e17.7))

10101 format(A,I3.3,A,A)
10102 format(A,A,A,I3.3,A)
10103 format(A,A,A,I3.3,A,I3.3,A)
10205 format(A,I4.4)

end subroutine PRINT_2D
