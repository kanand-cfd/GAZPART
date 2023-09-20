subroutine BEFORE_COLLISION
!!--------------------------------------------------------------
!!   			   (1): <up>
!!                  2 : <vp>
!!                  3 : <wp>
!!                  4 : <up.up>
!!                  5 : <vp.vp>
!!                  6 : <wp.wp>
!!                  7 : <up.vp>
!!                  8 : <up.wp>
!!                  9 : <vp.wp>
!!                 10 : <uf@p>
!!                 11 : <vf@p>
!!                 12 : <wf@p>
!!                 13 : <uf@p.uf@p>
!!                 14 : <vf@p.vf@p>
!!                 15 : <wf@p.wf@p>
!!                 16 : <uf@p.vf@p>
!!                 17 : <uf@p.wf@p>
!!                 18 : <vf@p.wf@p>
!!                 19 : <up.uf@p>
!!                 20 : <vp.vf@p>
!!                 21 : <wp.wf@p>
!!                 22 : <up.vf@p>
!!                 23 : <up.wf@p>
!!                 24 : <vp.wf@p>
!!--------------------------------------------------------------
use STATISTICS 
use PARTICLE_PARALLEL
use MPI_STRUCTURES
use DNS_DIM               
use PARAM_PHYS
use CHECK_CPU
use COLLISION_VARIABLE

implicit none

!---------------------------------------------------------------
integer :: I, J, IDP
!---------------------------------------------------------------


do J = 1, NIG

    do I = 1, NPART_LOC(J)

        IDP = PART(I,J)%IDP


        B_COLL(1, J, I, IDP, :) = PART(I,J)%UP
        B_COLL(2, J, I, IDP, :) = PART(I,J)%VP
        B_COLL(3, J, I, IDP, :) = PART(I,J)%WP
        B_COLL(4, J, I, IDP, :) = PART(I,J)%UP*PART(I,J)%UP
        B_COLL(5, J, I, IDP, :) = PART(I,J)%VP*PART(I,J)%VP
        B_COLL(6, J, I, IDP, :) = PART(I,J)%WP*PART(I,J)%WP
        B_COLL(7, J, I, IDP, :) = PART(I,J)%UP*PART(I,J)%VP
        B_COLL(8, J, I, IDP, :) = PART(I,J)%UP*PART(I,J)%WP
        B_COLL(9, J, I, IDP, :) = PART(I,J)%VP*PART(I,J)%WP
 


        PART(I,J)%COLL_BIDISP = .false.
        PART(I,J)%COLL_FLAG = .false.

        A_COLL(1, J, I, IDP, :) = PART(I,J)%UP
        A_COLL(2, J, I, IDP, :) = PART(I,J)%VP
        A_COLL(3, J, I, IDP, :) = PART(I,J)%WP
        A_COLL(4, J, I, IDP, :) = PART(I,J)%UP*PART(I,J)%UP
        A_COLL(5, J, I, IDP, :) = PART(I,J)%VP*PART(I,J)%VP
        A_COLL(6, J, I, IDP, :) = PART(I,J)%WP*PART(I,J)%WP
        A_COLL(7, J, I, IDP, :) = PART(I,J)%UP*PART(I,J)%VP
        A_COLL(8, J, I, IDP, :) = PART(I,J)%UP*PART(I,J)%WP
        A_COLL(9, J, I, IDP, :) = PART(I,J)%VP*PART(I,J)%WP
        A_COLL(10,J, I, IDP, :) = ZERO

    end do 

end do


end subroutine BEFORE_COLLISION
