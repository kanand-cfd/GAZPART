subroutine TIME_AVERAGED_STATISTICS
    

use DNS_DIM            !- Dimension
use STATISTICS         !- Statistics
use PARAM_PHYS
use GEOMETRIC_VARIABLE
use COLLISION_VARIABLE
use PARTICLE_PARALLEL

!0- Np
!1- alpha_p
!2- Major axis
!3- Aspect Ratio
!4- Density

!5- <up>
!6- <vp>
!7- <wp>

!8- <omega_x>
!9- <omega_y>
!10- <omega_z>

!11- global <omega_x>
!12- global <omega_y>
!13- global <omega_z>

!14- <qp> - Translational
!15- <qp> - Rotational
!16- <qp> - Total

!17- sqrt(<oxp.oxp>)
!18- sqrt(<oyp.oyp>)
!19- sqrt(<ozp.ozp>)

!20- sqrt(Global-<oxp.oxp>)
!21- sqrt(Global-<oyp.oyp>)
!22- sqrt(Global-<ozp.ozp>)

!23- fcol
!24- tau-col
!25- N-col
!26- Nover
!27- Overlap
!28- Overlap/Nover

!29- <up.up>
!30- <vp.vp>
!31- <wp.wp>

!32- <up.vp>
!33- <up.wp>
!34- <vp.wp>

!35- <oxp.oxp>
!36- <oyp.oyp>
!37- <ozp.ozp>

!38- <oxp.oyp>
!39- <oxp.ozp>
!40- <oyp.ozp>

!41- Ixx<oxp.oxp>
!42- Iyy<oyp.oyp>
!43- Izz<ozp.ozp>

!44- Global-<oxp.oxp>
!45- Global-<oyp.oyp>
!46- Global-<ozp.ozp>

!47- Global-<oxp.oyp>
!48- Global-<oyp.ozp>
!49- Global-<oxp.ozp>

!50- Global-<Ixx.oxp.oxp>
!51- Global-<Iyy.oyp.oyp>
!52- Global-<Izz.ozp.ozp>

!53- Global-<Ixy.oxp.oyp>
!54- Global-<Ixz.oxp.ozp>
!55- Global-<Iyz.oyp.ozp>

! Normalized by qp-Total
!56- <upup>p/(2/3qp)
!57- <vpvp>p/(2/3qp)
!58- <wpwp>p/(2/3qp)

!59- <upvp>p/(2/3qp)
!60- <upwp>p/(2/3qp)
!61- <vpwp>p/(2/3qp)

!62- Ixx<oxp.oxp>/(2/3qp)
!63- Iyy<oyp.oyp>/(2/3qp)
!64- Izz<ozp.ozp>/(2/3qp)

!65- Global-<Ixx.oxp.oxp>/(2/3qp)
!66- Global-<Ixx.oxp.oxp>/(2/3qp)
!67- Global-<Ixx.oxp.oxp>/(2/3qp)

! Normalized by qp-Translation
!68- <upup>p/(2/3qp)-Translation
!69- <vpvp>p/(2/3qp)-Translation
!70- <wpwp>p/(2/3qp)-Translation

!71- <upvp>p/(2/3qp)-Translation
!72- <upwp>p/(2/3qp)-Translation
!73- <vpwp>p/(2/3qp)-Translation

! Normalized by qp-Rotation
!74- Ixx<oxp.oxp>/(2/3qp)-Rotation
!75- Iyy<oyp.oyp>/(2/3qp)-Rotation
!76- Izz<ozp.ozp>/(2/3qp)-Rotation

!77- Global-<Ixx.oxp.oxp>/(2/3qp)-Rotation
!78- Global-<Iyy.oyp.oyp>/(2/3qp)-Rotation
!79- Global-<Izz.ozp.ozp>/(2/3qp)-Rotation
!=========================================================!


implicit none

!- File name
character(len=40) :: FILENAME
! File number
integer :: NUMFILE, NUMFILE1
!- index
integer :: IDP1, IDP2
!- Index
integer :: I, J, LGR, M ,N, K

!---------------------------------------------------------------------!
!---------------------------------------------------------------------!

do I = 1, NIG
   
    do IDP1 = 1, POLYDISP(I)

        NUMFILE = 800 + 10*I + IDP1

        write(NUMFILE, 10000) MEAN_TIME_PART(23, I, IDP1) !- Np
        write(NUMFILE, 10000) MEAN_TIME_PART(23, I, IDP1) * &
                              (4.0*PPI/3.0)*(EMAJ_PART(I,IDP1)**3)/(LXMAX*LYMAX*LZMAX * APR_PART(I)**2) !- alpha_p
        
        write(NUMFILE, 10000) EMAJ_PART(I,IDP1) !- Major axis

        write(NUMFILE, 10000) APR_PART(I)      !- Aspect Ratio
        write(NUMFILE, 10000) RHOP(I,IDP1)     !- Density


        write(NUMFILE, 10000) MEAN_TIME_PART(1, I, IDP1)  !- <up>
        write(NUMFILE, 10000) MEAN_TIME_PART(2, I, IDP1)  !- <vp>
        write(NUMFILE, 10000) MEAN_TIME_PART(3, I, IDP1)  !- <wp>

        write(NUMFILE, 10000) MEAN_TIME_PART(11, I, IDP1) !- <omega_x>
        write(NUMFILE, 10000) MEAN_TIME_PART(12, I, IDP1) !- <omega_y>
        write(NUMFILE, 10000) MEAN_TIME_PART(13, I, IDP1) !- <omega_z>

        write(NUMFILE, 10000) MEAN_TIME_PART(14, I, IDP1) !- global <omega_x>
        write(NUMFILE, 10000) MEAN_TIME_PART(15, I, IDP1) !- global <omega_y>
        write(NUMFILE, 10000) MEAN_TIME_PART(16, I, IDP1) !- global <omega_z>

        write(NUMFILE, 10000) MEAN_TIME_PART(10, I, IDP1) !- <qp> - Translational
        write(NUMFILE, 10000) MEAN_TIME_PART(40, I, IDP1) !- <qp> - Rotational
        write(NUMFILE, 10000) MEAN_TIME_PART(41, I, IDP1) !- <qp> - Total

        write(NUMFILE, 10000) sqrt(MEAN_TIME_PART(17, I, IDP1)) !- sqrt(<oxp.oxp>)
        write(NUMFILE, 10000) sqrt(MEAN_TIME_PART(18, I, IDP1)) !- sqrt(<oyp.oyp>)
        write(NUMFILE, 10000) sqrt(MEAN_TIME_PART(19, I, IDP1)) !- sqrt(<ozp.ozp>)

        write(NUMFILE, 10000) sqrt(MEAN_TIME_PART(27, I, IDP1)) !- sqrt(Global-<oxp.oxp>)
        write(NUMFILE, 10000) sqrt(MEAN_TIME_PART(28, I, IDP1)) !- sqrt(Global-<oyp.oyp>)
        write(NUMFILE, 10000) sqrt(MEAN_TIME_PART(29, I, IDP1)) !- sqrt(Global-<ozp.ozp>)



        if((PARTDEF(I)>1) .and. SOLVE_COLLISION) then

            do IDP2 = 1, POLYDISP(I)

                write(NUMFILE, 10000) MEAN_TIME_PART_COL(3,IDP1,IDP2)/MEAN_TIME_PART(23,I,IDP1) !- fcol
                write(NUMFILE, 10000) MEAN_TIME_PART(23,I,IDP1)/MEAN_TIME_PART_COL(3,IDP1,IDP2) !- tau-col
                write(NUMFILE, 10000) MEAN_TIME_PART_COL(2,IDP1,IDP2)                           !- N-col 
                write(NUMFILE, 10000) MEAN_TIME_PART_COL(1,IDP1,IDP2)                           !- Nover
                write(NUMFILE, 10000) MEAN_TIME_PART_COL(4,IDP1,IDP2)                           !- Overlap
                write(NUMFILE, 10000) MEAN_TIME_PART_COL(5,IDP1,IDP2)                           !- Overlap/Nover

            end do

        end if


        write(NUMFILE, 10000) MEAN_TIME_PART(4, I, IDP1)  !- <up.up>
        write(NUMFILE, 10000) MEAN_TIME_PART(5, I, IDP1)  !- <vp.vp>
        write(NUMFILE, 10000) MEAN_TIME_PART(6, I, IDP1)  !- <wp.wp>

        write(NUMFILE, 10000) MEAN_TIME_PART(7, I, IDP1)  !- <up.vp>
        write(NUMFILE, 10000) MEAN_TIME_PART(8, I, IDP1)  !- <up.wp>
        write(NUMFILE, 10000) MEAN_TIME_PART(9, I, IDP1)  !- <vp.wp>

        write(NUMFILE, 10000) MEAN_TIME_PART(17, I, IDP1) !- <oxp.oxp>
        write(NUMFILE, 10000) MEAN_TIME_PART(18, I, IDP1) !- <oyp.oyp>
        write(NUMFILE, 10000) MEAN_TIME_PART(19, I, IDP1) !- <ozp.ozp>

        write(NUMFILE, 10000) MEAN_TIME_PART(20, I, IDP1) !- <oxp.oyp>
        write(NUMFILE, 10000) MEAN_TIME_PART(21, I, IDP1) !- <oxp.ozp>
        write(NUMFILE, 10000) MEAN_TIME_PART(22, I, IDP1) !- <oyp.ozp>

        write(NUMFILE, 10000) MEAN_TIME_PART(24, I, IDP1) !- Ixx<oxp.oxp>
        write(NUMFILE, 10000) MEAN_TIME_PART(25, I, IDP1) !- Iyy<oyp.oyp>
        write(NUMFILE, 10000) MEAN_TIME_PART(26, I, IDP1) !- Izz<ozp.ozp>

        write(NUMFILE, 10000) MEAN_TIME_PART(27, I, IDP1) !- Global-<oxp.oxp>
        write(NUMFILE, 10000) MEAN_TIME_PART(28, I, IDP1) !- Global-<oyp.oyp>
        write(NUMFILE, 10000) MEAN_TIME_PART(29, I, IDP1) !- Global-<ozp.ozp>

        write(NUMFILE, 10000) MEAN_TIME_PART(30, I, IDP1) !- Global-<oxp.oyp>
        write(NUMFILE, 10000) MEAN_TIME_PART(31, I, IDP1) !- Global-<oyp.ozp>
        write(NUMFILE, 10000) MEAN_TIME_PART(32, I, IDP1) !- Global-<oxp.ozp>

        write(NUMFILE, 10000) MEAN_TIME_PART(33, I, IDP1) !- Global-<Ixx.oxp.oxp>
        write(NUMFILE, 10000) MEAN_TIME_PART(34, I, IDP1) !- Global-<Iyy.oyp.oyp>
        write(NUMFILE, 10000) MEAN_TIME_PART(35, I, IDP1) !- Global-<Izz.ozp.ozp>

        write(NUMFILE, 10000) MEAN_TIME_PART(36, I, IDP1) !- Global-<Ixy.oxp.oyp>
        write(NUMFILE, 10000) MEAN_TIME_PART(37, I, IDP1) !- Global-<Ixz.oxp.ozp>
        write(NUMFILE, 10000) MEAN_TIME_PART(38, I, IDP1) !- Global-<Iyz.oyp.ozp>

        ! Normalized by qp-Total
!        write(NUMFILE, 10000) MEAN_TIME_PART(4, I, IDP1)/MEAN_TIME_PART(41, I, IDP1)  !- <upup>p/(2/3qp)
!        write(NUMFILE, 10000) MEAN_TIME_PART(5, I, IDP1)/MEAN_TIME_PART(41, I, IDP1)  !- <vpvp>p/(2/3qp)
!        write(NUMFILE, 10000) MEAN_TIME_PART(6, I, IDP1)/MEAN_TIME_PART(41, I, IDP1)  !- <wpwp>p/(2/3qp)

!        write(NUMFILE, 10000) MEAN_TIME_PART(7, I, IDP1)/MEAN_TIME_PART(41, I, IDP1)  !- <upvp>p/(2/3qp)
!        write(NUMFILE, 10000) MEAN_TIME_PART(8, I, IDP1)/MEAN_TIME_PART(41, I, IDP1)  !- <upwp>p/(2/3qp)
!        write(NUMFILE, 10000) MEAN_TIME_PART(9, I, IDP1)/MEAN_TIME_PART(41, I, IDP1)  !- <vpwp>p/(2/3qp)

!        write(NUMFILE, 10000) MEAN_TIME_PART(24, I, IDP1)/MEAN_TIME_PART(41, I, IDP1) !- Ixx<oxp.oxp>/(2/3qp)
!        write(NUMFILE, 10000) MEAN_TIME_PART(25, I, IDP1)/MEAN_TIME_PART(41, I, IDP1) !- Iyy<oyp.oyp>/(2/3qp)
!        write(NUMFILE, 10000) MEAN_TIME_PART(26, I, IDP1)/MEAN_TIME_PART(41, I, IDP1) !- Izz<ozp.ozp>/(2/3qp)

!        write(NUMFILE, 10000) MEAN_TIME_PART(33, I, IDP1)/MEAN_TIME_PART(41, I, IDP1) !- Global-<Ixx.oxp.oxp>/(2/3qp)
!        write(NUMFILE, 10000) MEAN_TIME_PART(34, I, IDP1)/MEAN_TIME_PART(41, I, IDP1) !- Global-<Ixx.oxp.oxp>/(2/3qp)
!        write(NUMFILE, 10000) MEAN_TIME_PART(35, I, IDP1)/MEAN_TIME_PART(41, I, IDP1) !- Global-<Ixx.oxp.oxp>/(2/3qp)


        ! Normalized by qp-Translation
!        write(NUMFILE, 10000) MEAN_TIME_PART(4, I, IDP1)/MEAN_TIME_PART(10, I, IDP1)  !- <upup>p/(2/3qp)-Translation
!        write(NUMFILE, 10000) MEAN_TIME_PART(5, I, IDP1)/MEAN_TIME_PART(10, I, IDP1)  !- <vpvp>p/(2/3qp)-Translation
!        write(NUMFILE, 10000) MEAN_TIME_PART(6, I, IDP1)/MEAN_TIME_PART(10, I, IDP1)  !- <wpwp>p/(2/3qp)-Translation

!        write(NUMFILE, 10000) MEAN_TIME_PART(7, I, IDP1)/MEAN_TIME_PART(10, I, IDP1)  !- <upvp>p/(2/3qp)-Translation
!        write(NUMFILE, 10000) MEAN_TIME_PART(8, I, IDP1)/MEAN_TIME_PART(10, I, IDP1)  !- <upwp>p/(2/3qp)-Translation
!        write(NUMFILE, 10000) MEAN_TIME_PART(9, I, IDP1)/MEAN_TIME_PART(10, I, IDP1)  !- <vpwp>p/(2/3qp)-Translation


        ! Normalized by qp-Rotation
!        write(NUMFILE, 10000) MEAN_TIME_PART(24, I, IDP1)/MEAN_TIME_PART(40, I, IDP1)  !- Ixx<oxp.oxp>/(2/3qp)-Rotation
!        write(NUMFILE, 10000) MEAN_TIME_PART(25, I, IDP1)/MEAN_TIME_PART(40, I, IDP1)  !- Iyy<oyp.oyp>/(2/3qp)-Rotation
!        write(NUMFILE, 10000) MEAN_TIME_PART(26, I, IDP1)/MEAN_TIME_PART(40, I, IDP1)  !- Izz<ozp.ozp>/(2/3qp)-Rotation

!        write(NUMFILE, 10000) MEAN_TIME_PART(33, I, IDP1)/MEAN_TIME_PART(40, I, IDP1)  !- Global-<Ixx.oxp.oxp>/(2/3qp)-Rotation
!        write(NUMFILE, 10000) MEAN_TIME_PART(34, I, IDP1)/MEAN_TIME_PART(40, I, IDP1)  !- Global-<Iyy.oyp.oyp>/(2/3qp)-Rotation
!        write(NUMFILE, 10000) MEAN_TIME_PART(35, I, IDP1)/MEAN_TIME_PART(40, I, IDP1)  !- Global-<Izz.ozp.ozp>/(2/3qp)-Rotation


        if(FROZEN_FLOW == 1) then


            NUMFILE1 = 810+ 10*I + IDP1

            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID( 1,I,IDP1)   !- Rep
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID( 2,I,IDP1)   !- CDRAG
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID( 3,I,IDP1)   !- CLIFT
  
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID( 4,I,IDP1)   !- FDRAG(1,1)
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID( 5,I,IDP1)   !- FDRAG(2,1)
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID( 6,I,IDP1)   !- FDRAG(3,1)

            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID( 7,I,IDP1)   !- FLIFT(1,1)
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID( 8,I,IDP1)   !- FLIFT(2,1)
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID( 9,I,IDP1)   !- FLIFT(3,1)

            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(10,I,IDP1)   !- !! FDRAG !!
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(11,I,IDP1)   !- !! FLIFT !!

            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(12,I,IDP1)   !- Re_R
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(13,I,IDP1)   !- CTPITCH
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(14,I,IDP1)   !- CTROTATION
              
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(15,I,IDP1)   !- TORQUE_PITCHING(1,1)
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(16,I,IDP1)   !- TORQUE_PITCHING(2,1)
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(17,I,IDP1)   !- TORQUE_PITCHING(3,1)

            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(18,I,IDP1)   !- TORQUE_ROTATION(1,1) 
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(19,I,IDP1)   !- TORQUE_ROTATION(2,1)
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(20,I,IDP1)   !- TORQUE_ROTATION(3,1)

            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(21,I,IDP1)   !- !! TORQUE_PITCHING !! 
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(22,I,IDP1)   !- !! TORQUE_ROTATION !!

            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(23,I,IDP1)   !- TAUPF
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(24,I,IDP1)   !- Gravity
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(25,I,IDP1)   !-phi

            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(26,I,IDP1)   !-principal_axis(1,1)    
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(27,I,IDP1)   !-principal_axis(2,1)
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(28,I,IDP1)   !-principal_axis(3,1)

            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(29,I,IDP1)   !-INVTAUP 

            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(30,I,IDP1)   !-UP 
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(31,I,IDP1)   !-VP 
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(32,I,IDP1)   !-WP 

            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(33,I,IDP1)   !-OXP 
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(34,I,IDP1)   !-OYP 
            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(35,I,IDP1)   !-OZP  

            write(NUMFILE1, 10000) MEAN_TIME_PARTFLUID(36,I,IDP1)   !- <Np>


        end if 


    end do

end do


10000 format(2(e17.7))
    
end subroutine TIME_AVERAGED_STATISTICS