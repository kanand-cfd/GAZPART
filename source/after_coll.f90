subroutine AFTER_COLLISION
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
!!                 10 : Number of Collision
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
use STATISTICS 
use PARTICLE_PARALLEL
use MPI_STRUCTURES
use DNS_DIM               
use PARAM_PHYS
use CHECK_CPU
use COLLISION_VARIABLE

implicit none

!---------------------------------------------------------------
integer :: I, J, IDP1, IDP2
!---------------------------------------------------------------

!A_COLL(:,:,:,:,:) = ZERO


do J = 1, NIG

    do I = 1, NPART_LOC(J)

        IDP1 = PART(I,J)%IDP
        IDP2 = PART(I,J)%IDP

        if (PART(I,J)%COLL_BIDISP .and. IDP1 == 1) then
            IDP2 = 2
        elseif (PART(I,J)%COLL_BIDISP .and. IDP1 == 2) then
            IDP2 = 1
        end if 

        A_COLL(1, J, I, IDP1, IDP2) = PART(I,J)%UP
        A_COLL(2, J, I, IDP1, IDP2) = PART(I,J)%VP
        A_COLL(3, J, I, IDP1, IDP2) = PART(I,J)%WP
        A_COLL(4, J, I, IDP1, IDP2) = PART(I,J)%UP*PART(I,J)%UP
        A_COLL(5, J, I, IDP1, IDP2) = PART(I,J)%VP*PART(I,J)%VP
        A_COLL(6, J, I, IDP1, IDP2) = PART(I,J)%WP*PART(I,J)%WP
        A_COLL(7, J, I, IDP1, IDP2) = PART(I,J)%UP*PART(I,J)%VP
        A_COLL(8, J, I, IDP1, IDP2) = PART(I,J)%UP*PART(I,J)%WP
        A_COLL(9, J, I, IDP1, IDP2) = PART(I,J)%VP*PART(I,J)%WP

        if(PART(I,J)%COLL_FLAG) A_COLL(10,J,I,IDP1,IDP2) = 1.0


!if (PART(I,J)%COLL_BIDISP) write(*,*) B_COLL(3,I,J,1,2), A_COLL(3,I,J,1,2)


    end do 

end do



end subroutine AFTER_COLLISION
