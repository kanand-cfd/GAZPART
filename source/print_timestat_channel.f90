subroutine PRINT_TIMESTAT_CHANNEL
!!===================================================================!!
!!            Finalize the channel statistics and print              !!
!!===================================================================!!

use DNS_DIM            !- Dimension
use PARAM_PHYS         !- Physical & numerical parameters
use PARTICLE_PARALLEL
use GEOMETRIC_VARIABLE !- Mesh parameters
use STATISTICS         !- Statistics
use CHECK_CPU 
use FLUID_VARIABLE        

implicit none

!--------------------------------------------------------------------
!integer, parameter :: NPDF = 1000

integer :: I, J, IDP, NS, IDP1, IDP2

real(kind=8) :: DELX, NP_SLICE, XM1

! File name
character(len=40) :: FILENAME

! File Number
integer :: NUMFILE1, NUMFILE2, NUMFILE3
!--------------------------------------------------------------------

!!-------------------------------------------------------------------
!! Time-averaged Statistics
!!-------------------------------------------------------------------
MEAN_TIME_PART_CHAN = MEAN_TIME_PART_CHAN / real(NEVEN_CHAN)
!MEAN_TIME_PDF_CHAN = MEAN_TIME_PDF_CHAN / real(NEVEN_PDF)
!!-------------------------------------------------------------------


!! Averaging over Slices
do J = 1, NIG

    do IDP = 1, POLYDISP(J) 

        do NS = 1, NSLICE


            NP_SLICE = MEAN_TIME_PART_CHAN(44, J, IDP, NS)

            MEAN_TIME_PART_CHAN(1:43,J,IDP,NS) = MEAN_TIME_PART_CHAN(1:43,J,IDP,NS) / NP_SLICE

        end do

    end do

end do


!!-------------------------------------------------------------------
!! Calculation of Collision Terms
!!-------------------------------------------------------------------
MEAN_TIME_COLL_CHAN(:9, :, :, :, :) = MEAN_TIME_COLL_CHAN(:9, :, :, :, :) / (DTIME)
MEAN_TIME_COLL_CHAN = MEAN_TIME_COLL_CHAN / real(NEVEN_COLL_CHAN)
!!-------------------------------------------------------------------

do J = 1, NIG
    
    do IDP = 1, POLYDISP(J)

        NUMFILE1 = 8000 + J*100 + IDP*10

        do NS = 1, NSLICE

            write(NUMFILE1,10000)   &
            XSTAT(J,IDP,NS)+0.5*DEL_X(J,IDP,NS), DEL_X(J,IDP,NS), &  
            MEAN_TIME_PART_CHAN( 1,J,IDP,NS),   & !- <up>
            MEAN_TIME_PART_CHAN( 2,J,IDP,NS),   & !- <vp>
            MEAN_TIME_PART_CHAN( 3,J,IDP,NS),   & !- <wp>

            MEAN_TIME_PART_CHAN(10,J,IDP,NS),   & !- <uf@p>
            MEAN_TIME_PART_CHAN(11,J,IDP,NS),   & !- <vf@p>
            MEAN_TIME_PART_CHAN(12,J,IDP,NS),   & !- <wf@p>

            MEAN_TIME_PART_CHAN( 4,J,IDP,NS),   & !- <up.up>
            MEAN_TIME_PART_CHAN( 5,J,IDP,NS),   & !- <vp.vp>
            MEAN_TIME_PART_CHAN( 6,J,IDP,NS),   & !- <wp.wp>

            MEAN_TIME_PART_CHAN( 7,J,IDP,NS),   & !- <up.vp>
            MEAN_TIME_PART_CHAN( 8,J,IDP,NS),   & !- <up.wp>
            MEAN_TIME_PART_CHAN( 9,J,IDP,NS),   & !- <vp.wp>

            MEAN_TIME_PART_CHAN(13,J,IDP,NS),   & !- <uf@p.uf@p>
            MEAN_TIME_PART_CHAN(14,J,IDP,NS),   & !- <vf@p.vf@p>
            MEAN_TIME_PART_CHAN(15,J,IDP,NS),   & !- <wf@p.wf@p>

            MEAN_TIME_PART_CHAN(16,J,IDP,NS),   & !- <omega_px>
            MEAN_TIME_PART_CHAN(17,J,IDP,NS),   & !- <omega_py>
            MEAN_TIME_PART_CHAN(18,J,IDP,NS),   & !- <omega_pz>

            MEAN_TIME_PART_CHAN(19,J,IDP,NS),   & !- <omega_px.omega_px>
            MEAN_TIME_PART_CHAN(20,J,IDP,NS),   & !- <omega_py.omega_py>
            MEAN_TIME_PART_CHAN(21,J,IDP,NS),   & !- <omega_pz.omega_pz>

            MEAN_TIME_PART_CHAN(22,J,IDP,NS),   & !- <omega_px.omega_py>
            MEAN_TIME_PART_CHAN(23,J,IDP,NS),   & !- <omega_px.omega_pz>
            MEAN_TIME_PART_CHAN(24,J,IDP,NS),   & !- <omega_py.omega_pz>

            MEAN_TIME_PART_CHAN(25,J,IDP,NS),   & !- F-Drag-X
            MEAN_TIME_PART_CHAN(26,J,IDP,NS),   & !- F-Drag-Y
            MEAN_TIME_PART_CHAN(27,J,IDP,NS),   & !- F-Drag-Z

            MEAN_TIME_PART_CHAN(28,J,IDP,NS),   & !- F-Lift-X
            MEAN_TIME_PART_CHAN(29,J,IDP,NS),   & !- F-Lift-Y
            MEAN_TIME_PART_CHAN(30,J,IDP,NS),   & !- F-Lift-Z

            MEAN_TIME_PART_CHAN(31,J,IDP,NS),   & !- Pitching Torque-X
            MEAN_TIME_PART_CHAN(32,J,IDP,NS),   & !- Pitching Torque-Y
            MEAN_TIME_PART_CHAN(33,J,IDP,NS),   & !- Pitching Torque-Z

            MEAN_TIME_PART_CHAN(34,J,IDP,NS),   & !- Rotation Torque-X
            MEAN_TIME_PART_CHAN(35,J,IDP,NS),   & !- Rotation Torque-Y
            MEAN_TIME_PART_CHAN(36,J,IDP,NS),   & !- Rotation Torque-Z

            MEAN_TIME_PART_CHAN(37,J,IDP,NS),   & !- Rep
            MEAN_TIME_PART_CHAN(38,J,IDP,NS),   & !- phi
            MEAN_TIME_PART_CHAN(39,J,IDP,NS),   & !- CDRAG
            MEAN_TIME_PART_CHAN(40,J,IDP,NS),   & !- CLIFT

            MEAN_TIME_PART_CHAN(41,J,IDP,NS),   & !- CTPITCH
            MEAN_TIME_PART_CHAN(42,J,IDP,NS),   & !- CTROTATION
            MEAN_TIME_PART_CHAN(43,J,IDP,NS),   & !- CTROTATION_2

            MEAN_TIME_PART_CHAN(44,J,IDP,NS)      !- Np

        end do



        do IDP2 = 1, POLYDISP(J)

            NUMFILE2 = 5000 + J*200 + 10*IDP + IDP2

            do NS = 1, NSLICE

                write(NUMFILE2,10000)   &
                XSTAT(J,IDP,NS)+0.5*DEL_X(J,IDP,NS),       & !- x
                DEL_X(J,IDP,NS),                           & !- del x 
                MEAN_TIME_COLL_CHAN( 1,J, IDP, IDP2,NS),   & !- <up>
                MEAN_TIME_COLL_CHAN( 2,J, IDP, IDP2,NS),   & !- <vp>
                MEAN_TIME_COLL_CHAN( 3,J, IDP, IDP2,NS),   & !- <wp>
                MEAN_TIME_COLL_CHAN( 4,J, IDP, IDP2,NS),   & !- <up.up>
                MEAN_TIME_COLL_CHAN( 5,J, IDP, IDP2,NS),   & !- <vp.vp>
                MEAN_TIME_COLL_CHAN( 6,J, IDP, IDP2,NS),   & !- <wp.wp>
                MEAN_TIME_COLL_CHAN( 7,J, IDP, IDP2,NS),   & !- <up.vp>
                MEAN_TIME_COLL_CHAN( 8,J, IDP, IDP2,NS),   & !- <up.wp>
                MEAN_TIME_COLL_CHAN( 9,J, IDP, IDP2,NS),   & !- <vp.wp>
                MEAN_TIME_COLL_CHAN(10,J, IDP, IDP2,NS),   & !- f_col
                MEAN_TIME_COLL_CHAN(11,J, IDP, IDP2,NS)      !- 1/tau_col
                
            end do

        end do

    end do

end do

!if(LEVELX_PDF) then

    !! Printing one point pdf statistics
!do J = 1, NIG

!    do IDP = 1, POLYDISP(J)

!         NUMFILE1 = 10000 + J*1000 + IDP*100
!         NUMFILE2 = 20000 + J*1000 + IDP*100
!         NUMFILE3 = 30000 + J*1000 + IDP*100
!
!
!         do I = 1, NPDF
!
!            write(NUMFILE1,10000) real(I), &
!            MEAN_TIME_PDF_CHAN(1, J, IDP, NSLICE/2, I), &
!            MEAN_TIME_PDF_CHAN(2, J, IDP, NSLICE/2, I), &
!            MEAN_TIME_PDF_CHAN(3, J, IDP, NSLICE/2, I)
!
!
!            write(NUMFILE2,10000) real(I), &
!            MEAN_TIME_PDF_CHAN(1, J, IDP, NSLICE/4, I), &
!            MEAN_TIME_PDF_CHAN(2, J, IDP, NSLICE/4, I), &
!            MEAN_TIME_PDF_CHAN(3, J, IDP, NSLICE/4, I)
!
!            write(NUMFILE3,10000) real(I), &
!            MEAN_TIME_PDF_CHAN(1, J, IDP, NSLICE-2, I), &
!            MEAN_TIME_PDF_CHAN(2, J, IDP, NSLICE-2, I), &
!            MEAN_TIME_PDF_CHAN(3, J, IDP, NSLICE-2, I), &
!            MEAN_TIME_PDF_CHAN(1, J, IDP, 3, I),        &
!            MEAN_TIME_PDF_CHAN(2, J, IDP, 3, I),        &
!            MEAN_TIME_PDF_CHAN(3, J, IDP, 3, I)
!
!
!         end do

!    end do

!end do



!do J = 1, NIG
!
!    do IDP = 1, POLYDISP(J)
!
!        do IDP2 = 1, POLYDISP(J)
!
!            NUMFILE1 = 40000 + J*1000 + IDP*100 + IDP2
!            NUMFILE2 = 50000 + J*1000 + IDP*100 + IDP2
!            NUMFILE3 = 60000 + J*1000 + IDP*100 + IDP2
!
!            do I = 1, NPDF
!
!                write(NUMFILE1, 10000) real(I), &
!                MEAN_TIME_PDF_RELV_CHAN(1, J, IDP, IDP2, NSLICE/2, I), &
!                MEAN_TIME_PDF_RELV_CHAN(2, J, IDP, IDP2, NSLICE/2, I), &
!                MEAN_TIME_PDF_RELV_CHAN(3, J, IDP, IDP2, NSLICE/2, I), &
!                MEAN_TIME_PDF_RELV_CHAN(4, J, IDP, IDP2, NSLICE/2, I), &
!                MEAN_TIME_PDF_RELV_CHAN(5, J, IDP, IDP2, NSLICE/2, I), &
!                MEAN_TIME_PDF_RELV_CHAN(6, J, IDP, IDP2, NSLICE/2, I), &
!                MEAN_TIME_PDF_RELV_CHAN(7, J, IDP, IDP2, NSLICE/2, I), &
!                MEAN_TIME_PDF_RELV_CHAN(8, J, IDP, IDP2, NSLICE/2, I), &
!                MEAN_TIME_PDF_RELV_CHAN(9, J, IDP, IDP2, NSLICE/2, I)
!
!
!                write(NUMFILE2, 10000) real(I), &
!                MEAN_TIME_PDF_RELV_CHAN(1, J, IDP, IDP2, NSLICE/4, I), &
!                MEAN_TIME_PDF_RELV_CHAN(2, J, IDP, IDP2, NSLICE/4, I), &
!                MEAN_TIME_PDF_RELV_CHAN(3, J, IDP, IDP2, NSLICE/4, I), &
!                MEAN_TIME_PDF_RELV_CHAN(4, J, IDP, IDP2, NSLICE/4, I), &
!                MEAN_TIME_PDF_RELV_CHAN(5, J, IDP, IDP2, NSLICE/4, I), &
!                MEAN_TIME_PDF_RELV_CHAN(6, J, IDP, IDP2, NSLICE/4, I), &
!                MEAN_TIME_PDF_RELV_CHAN(7, J, IDP, IDP2, NSLICE/4, I), &
!                MEAN_TIME_PDF_RELV_CHAN(8, J, IDP, IDP2, NSLICE/4, I), &
!                MEAN_TIME_PDF_RELV_CHAN(9, J, IDP, IDP2, NSLICE/4, I)
!
!
!                write(NUMFILE3, 10000) real(I), &
!                MEAN_TIME_PDF_RELV_CHAN(1, J, IDP, IDP2, NSLICE-2, I), &
!                MEAN_TIME_PDF_RELV_CHAN(2, J, IDP, IDP2, NSLICE-2, I), &
!                MEAN_TIME_PDF_RELV_CHAN(3, J, IDP, IDP2, NSLICE-2, I), &
!                MEAN_TIME_PDF_RELV_CHAN(4, J, IDP, IDP2, NSLICE-2, I), &
!                MEAN_TIME_PDF_RELV_CHAN(5, J, IDP, IDP2, NSLICE-2, I), &
!                MEAN_TIME_PDF_RELV_CHAN(6, J, IDP, IDP2, NSLICE-2, I), &
!                MEAN_TIME_PDF_RELV_CHAN(7, J, IDP, IDP2, NSLICE-2, I), &
!                MEAN_TIME_PDF_RELV_CHAN(8, J, IDP, IDP2, NSLICE-2, I), &
!                MEAN_TIME_PDF_RELV_CHAN(9, J, IDP, IDP2, NSLICE-2, I)
!
!
!            end do
!
!        end do
!
!    end do
!
!end do
!
!end if



10000 format (60(e17.7))

end subroutine PRINT_TIMESTAT_CHANNEL
