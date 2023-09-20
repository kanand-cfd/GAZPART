
subroutine STAT_COLL_PART_CHAN(NCYCLE, TIME)

!!--------------------------------------------------------------
!!   		   (1): <up>
!!                  2 : <vp>
!!                  3 : <wp>
!!                  4 : <up.up>
!!                  5 : <vp.vp>
!!                  6 : <wp.wp>
!!                  7 : <up.vp>
!!                  8 : <up.wp>
!!                  9 : <vp.wp>
!!                 10 : Number of Collisions
!!                 11 : <uf@p>
!!                 12 : <vf@p>
!!                 13 : <wf@p>
!!                 14 : <uf@p.uf@p>
!!                 15 : <vf@p.vf@p>
!!                 16 : <wf@p.wf@p>
!!                 17 : <uf@p.vf@p>
!!                 18 : <uf@p.wf@p>
!!                 19 : <vf@p.wf@p>
!!                 20 : <up.uf@p>
!!                 21 : <vp.vf@p>
!!                 22 : <wp.wf@p>
!!                 23 : <up.vf@p>
!!                 24 : <up.wf@p>
!!                 25 : <vp.wf@p>
!!--------------------------------------------------------------


use DNS_DIM            !- Dimension
use PARAM_PHYS         !- Physical & numerical parameters
use PARTICLE_PARALLEL
use GEOMETRIC_VARIABLE !- Mesh parameters
use STATISTICS         !- Statistics
use CHECK_CPU 
use FLUID_VARIABLE        

implicit none

!---------------------------------------------------------------------
!- Global Variables
!---------------------------------------------------------------------
!- Curent time
real(kind=8), intent(in) :: NCYCLE, TIME
!--------------------------------------------------------------------
!	ARRAYS STATEMENT
!--------------------------------------------------------------------
!- Local arrays
!---------------------------------------------------------------------
!- Arrays containing statistics
real(kind=8), dimension(NSTAT,NIG,POLYDISPMAX,POLYDISPMAX, NSLICE) :: MEAN_COLL_CHAN_LOC
real(kind=8), dimension(NSTAT,NIG,POLYDISPMAX,POLYDISPMAX, NSLICE) :: MEAN_COLL_CHAN

!- Index of particle class for polydisperse stat
integer :: IDP1, IDP2 

!- Time control variable
real(kind=8) :: TIME_START, TIME_END, n_col

integer :: I,J,K, IFLAG1, IPART, JPART, KPART, NS

!- file number
integer :: NUMFILE

real(kind=8), dimension(NIG, POLYDISPMAX, NSLICE) :: Number_Particles
!---------------------------------------------------------------------

!!- CPU check
if(MYID == 0) then
 TIME_START=MPI_WTIME()
end if


MEAN_COLL_CHAN_LOC(:,:,:,:,:) = ZERO
MEAN_COLL_CHAN(:,:,:,:,:) = ZERO

Number_Particles(:,:,:) = ZERO


do J = 1, NIG

    do I = 1, NPART_LOC(J)

        IDP1 = PART(I,J)%IDP
!--------------------------------------------------------------------
! If particle lies in the vicinity (+-delx) of slice
!--------------------------------------------------------------------
        call LOCATE_PART(J, PART(I,J), 2.0*EMAJ_PART(J,IDP1)/APR_PART(J), NS)

        Number_Particles(J, IDP1, NS) = Number_Particles(J, IDP1, NS) + 1.0

!--------------------------------------------------------------------
        do IDP2 = 1, POLYDISP(J)
            
!- <up>
            MEAN_COLL_CHAN_LOC(1,J,IDP1,IDP2,NS) = MEAN_COLL_CHAN_LOC(1,J,IDP1,IDP2,NS) + (A_COLL(1, J, I, IDP1, IDP2) - B_COLL(1, J, I, IDP1, IDP2))

!- <vp>
            MEAN_COLL_CHAN_LOC(2,J,IDP1,IDP2,NS) = MEAN_COLL_CHAN_LOC(2,J,IDP1,IDP2,NS) + (A_COLL(2, J, I, IDP1, IDP2) - B_COLL(2, J, I, IDP1, IDP2))

!- <wp>
            MEAN_COLL_CHAN_LOC(3,J,IDP1,IDP2,NS) = MEAN_COLL_CHAN_LOC(3,J,IDP1,IDP2,NS) + (A_COLL(3, J, I, IDP1, IDP2) - B_COLL(3, J, I, IDP1, IDP2))

!- <up.up>
            MEAN_COLL_CHAN_LOC(4,J,IDP1,IDP2,NS) = MEAN_COLL_CHAN_LOC(4,J,IDP1,IDP2,NS) + (A_COLL(4, J, I, IDP1, IDP2) - B_COLL(4, J, I, IDP1, IDP2))

!- <vp.vp>
            MEAN_COLL_CHAN_LOC(5,J,IDP1,IDP2,NS) = MEAN_COLL_CHAN_LOC(5,J,IDP1,IDP2,NS) + (A_COLL(5, J, I, IDP1, IDP2) - B_COLL(5, J, I, IDP1, IDP2))

!- <wp.wp>
            MEAN_COLL_CHAN_LOC(6,J,IDP1,IDP2,NS) = MEAN_COLL_CHAN_LOC(6,J,IDP1,IDP2,NS) + (A_COLL(6, J, I, IDP1, IDP2) - B_COLL(6, J, I, IDP1, IDP2))

!- <up.vp>
            MEAN_COLL_CHAN_LOC(7,J,IDP1,IDP2,NS) = MEAN_COLL_CHAN_LOC(7,J,IDP1,IDP2,NS) + (A_COLL(7, J, I, IDP1, IDP2) - B_COLL(7, J, I, IDP1, IDP2))

!- <up.wp>
            MEAN_COLL_CHAN_LOC(8,J,IDP1,IDP2,NS) = MEAN_COLL_CHAN_LOC(8,J,IDP1,IDP2,NS) + (A_COLL(8, J, I, IDP1, IDP2) - B_COLL(8, J, I, IDP1, IDP2))

!- <vp.wp>
            MEAN_COLL_CHAN_LOC(9,J,IDP1,IDP2,NS) = MEAN_COLL_CHAN_LOC(9,J,IDP1,IDP2,NS) + (A_COLL(9, J, I, IDP1, IDP2) - B_COLL(9, J, I, IDP1, IDP2))

!- Number of Collisions
            MEAN_COLL_CHAN_LOC(10,J,IDP1,IDP2,NS) = MEAN_COLL_CHAN_LOC(10,J,IDP1,IDP2,NS) + (A_COLL(10, J, I, IDP1, IDP2))


        end do


    end do 

end do 



!!====================================================================
!! Summation overall domain and normalization
!!====================================================================
!if (NPROC>1) then
do J = 1, NIG

	do IDP1 = 1, POLYDISP(J)

        do IDP2 = 1, POLYDISP(J)


                do NS = 1, NSLICE

                        MEAN_COLL_CHAN_LOC(10, J, IDP1, IDP2, NS) = MEAN_COLL_CHAN_LOC(10, J, IDP1, IDP2, NS)/DTIME

                        MEAN_COLL_CHAN_LOC(11, J, IDP1, IDP2, NS) = MEAN_COLL_CHAN_LOC(10, J, IDP1, IDP2, NS)/Number_Particles(J, IDP1, NS)

                        call MPI_ALLREDUCE(MEAN_COLL_CHAN_LOC(:,J,IDP1,IDP2,NS),MEAN_COLL_CHAN(:,J,IDP1,IDP2,NS),NSTAT,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

                end do

        end do

    end do

end do


!!======================================================================
!! Time averaging if so
!!======================================================================
if(STAT_TIME) then
 MEAN_TIME_COLL_CHAN = MEAN_TIME_COLL_CHAN + MEAN_COLL_CHAN
end if 



!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()
 CPU_PART(6) = CPU_PART(6) + TIME_END - TIME_START
end if


10000 format (30(e17.7))


end subroutine STAT_COLL_PART_CHAN
