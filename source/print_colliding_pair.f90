subroutine print_colliding_pair(IG, NP1, NP2, NPRINT, &
								RX1, RX2, cnormal)

use PARTICLE_PARALLEL
use PARAM_PHYS
use mod_quaternion
use DNS_DIM

implicit none

!!=====================================================================!!
!!=====================================================================!!
integer, intent(in) :: IG, NP1, NP2, NPRINT
real(kind=8), dimension(ndim, 1), intent(in) :: RX1, RX2, cnormal

!- File name
character(len=50) :: FILENAME
character(len=8) :: NUMFILE

!!=====================================================================!!
!!=====================================================================!!
if(NPRINT == 1) then

    write(FILENAME,10102)'postprocessing/part_coll_pair_p',IG,'_t',NPRINT,'.bin'
    open(unit=160,file=trim(FILENAME),status='replace',form='unformatted')

    write(160)  cnormal(1,1), cnormal(2,1), cnormal(3,1), &
                PART(NP1,IG)%ELLQUAT%a, PART(NP1,IG)%ELLQUAT%b, PART(NP1,IG)%ELLQUAT%c, PART(NP1,IG)%ELLQUAT%d, &
                PART(NP2,IG)%ELLQUAT%a, PART(NP2,IG)%ELLQUAT%b, PART(NP2,IG)%ELLQUAT%c, PART(NP2,IG)%ELLQUAT%d

else if(NPRINT .gt. 1 .and. NPRINT .le. 200000) then


    write(160)  cnormal(1,1), cnormal(2,1), cnormal(3,1), &
                PART(NP1,IG)%ELLQUAT%a, PART(NP1,IG)%ELLQUAT%b, PART(NP1,IG)%ELLQUAT%c, PART(NP1,IG)%ELLQUAT%d, &
                PART(NP2,IG)%ELLQUAT%a, PART(NP2,IG)%ELLQUAT%b, PART(NP2,IG)%ELLQUAT%c, PART(NP2,IG)%ELLQUAT%d

else

    close(160)
    write(*,*) 'Collisions Sampled'
    return

end if

!write(160)PART(NP1,IG)%XP, PART(NP2,IG)%XP
!write(160)PART(NP1,IG)%YP, PART(NP2,IG)%YP
!write(160)PART(NP1,IG)%ZP, PART(NP2,IG)%ZP
!write(160)PART(NP1,IG)%UP, PART(NP2,IG)%UP
!write(160)PART(NP1,IG)%VP, PART(NP2,IG)%VP
!write(160)PART(NP1,IG)%WP, PART(NP2,IG)%WP
!write(160)PART(NP1,IG)%ELLQUAT%a, PART(NP2,IG)%ELLQUAT%a
!write(160)PART(NP1,IG)%ELLQUAT%b, PART(NP2,IG)%ELLQUAT%b
!write(160)PART(NP1,IG)%ELLQUAT%c, PART(NP2,IG)%ELLQUAT%c
!write(160)PART(NP1,IG)%ELLQUAT%d, PART(NP2,IG)%ELLQUAT%d
!write(160)PART(NP1,IG)%OMEGAX, PART(NP2,IG)%OMEGAX
!write(160)PART(NP1,IG)%OMEGAY, PART(NP2,IG)%OMEGAY
!write(160)PART(NP1,IG)%OMEGAZ, PART(NP2,IG)%OMEGAZ
!write(160)RX1(1,1), RX2(2,1), RX1(3,1)
!write(160)RX1(1,1), RX2(2,1), RX2(3,1)
!write(160)cnormal(1,1), cnormal(2,1), cnormal(3,1)
!close(160)




10102 format (A,I2.2,A,I8.8,A)
    
end subroutine print_colliding_pair