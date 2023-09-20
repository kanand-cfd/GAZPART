subroutine COLLISION_RESPONSE(IG, IDP, &
                             part1, part2, &
                             R_impact1, R_impact2, &
                             normal)

use PARTICLE_PARALLEL
use DNS_DIM               
use PARAM_PHYS
use GEOMETRIC_VARIABLE
use COLLISION_VARIABLE
use STATISTICS
use mod_quaternion

implicit none

!=================================================!
!=================================================!
integer, intent(in) :: IG, IDP

type(PARTTYPE), intent(inout) :: part1, part2

real(kind=8), dimension(ndim, 1), intent(in) :: R_impact1, R_impact2, normal

!integer, intent(in) :: I, J
!=================================================!
!=================================================!
! Velocity at contact point
real(kind=8), dimension(ndim, 1) :: UPC1, UPC2

! Relative Velocity of Contact Points
real(kind=8), dimension(ndim, 1) :: UREL

real(kind=8) :: URN

! Velocity Difference 
real(kind=8), dimension(ndim, 1) :: DELU
!=================================================!
!=================================================!
! Inverse Mass Identity Matrix
real(kind=8), dimension(ndim, ndim) :: m_I

! Inertia Matrix and its Inverse
real(kind=8), dimension(ndim, ndim) :: I_p1, I_p2, Inv_Ip1, Inv_Ip2

! Skew Symmetric Matrix of the Impact Arm vector and its transpose
real(kind=8), dimension(ndim, ndim) :: R_X1, R_X_Tr1, R_X2, R_X_Tr2

! Final mass-Inetria Symmetric Matrix and its Inverse
real(kind=8), dimension(ndim, ndim) :: K1, K2, K, invK

!=================================================!
!=================================================!

! Impulse vector with Inertia
real(kind=8), dimension(ndim, 1) :: PC, K_N, V_IMP
!real(kind=8) :: PCX, PCY, PCZ

! Tangential and Normal Component of Impulse Vector
real(kind=8) :: PCT, PCN, KN_dot_N, PCT_X, PCT_Y, PCT_Z

real(kind=8) :: NRM_TW

real(kind=8), dimension(ndim, 1) :: TW

! Velcity Vector
real(kind=8), dimension(ndim, 1) :: U_P1, U_P2

! Angular Velocity vector
real(kind=8), dimension(ndim, 1) :: OMG1, OMG2

real(kind=8) :: E_initial1, E_initial_trans1, E_initial_rot1, E_final1, E_final_trans1, E_final_rot1
real(kind=8) :: E_initial2, E_initial_trans2, E_initial_rot2, E_final2, E_final_trans2, E_final_rot2
real(kind=8), dimension(ndim, 1) :: Angular_Momentum1, Angular_Momentum2
real(kind=8) :: E_initial, E_final
!=================================================!
!=================================================!

U_P1(1,1) = part1%UP
U_P1(2,1) = part1%VP
U_P1(3,1) = part1%WP

OMG1(1,1) = part1%OMEGAX
OMG1(2,1) = part1%OMEGAY
OMG1(3,1) = part1%OMEGAZ

!write(*,*) " "
!write(*,*) "E_initial_rot (before tranformation) : ", E_initial_rot
!write(*,*) "Particle 1 Angular Velocity before transformation", OMG1(1,1), OMG1(2,1), OMG1(3,1)
call transform_basis(OMG1, part1%ELLQUAT, shape(OMG1))

U_P2(1,1) = part2%UP
U_P2(2,1) = part2%VP
U_P2(3,1) = part2%WP

OMG2(1,1) = part2%OMEGAX
OMG2(2,1) = part2%OMEGAY
OMG2(3,1) = part2%OMEGAZ

!write(*,*) "Particle 2 Angular Velocity before transformation", OMG2(1,1), OMG2(2,1), OMG2(3,1)
call transform_basis(OMG2, part2%ELLQUAT, shape(OMG2))
!write(*,*) " "


! Calculate the velocity at Contact Point
! Particle 1
UPC1(1,1) = part1%UP + (OMG1(2,1)*R_impact1(3,1) - OMG1(3,1)*R_impact1(2,1))
UPC1(2,1) = part1%VP + (OMG1(3,1)*R_impact1(1,1) - OMG1(1,1)*R_impact1(3,1))
UPC1(3,1) = part1%WP + (OMG1(1,1)*R_impact1(2,1) - OMG1(2,1)*R_impact1(1,1))

! Particle 2
UPC2(1,1) = part2%UP + (OMG2(2,1)*R_impact2(3,1) - OMG2(3,1)*R_impact2(2,1))
UPC2(2,1) = part2%VP + (OMG2(3,1)*R_impact2(1,1) - OMG2(1,1)*R_impact2(3,1))
UPC2(3,1) = part2%WP + (OMG2(1,1)*R_impact2(2,1) - OMG2(2,1)*R_impact2(1,1))


! Relative Velocity of Contact Points
UREL = UPC2 - UPC1

! Dot Product of Relative Velocity and normal vector
URN = UREL(1,1)*normal(1,1) + UREL(2,1)*normal(2,1) + UREL(3,1)*normal(3,1)

! Velocity Difference DELU
DELU(:,1) = UREL(:,1) + ECP(IG)*URN*normal(:,1)

!================================================================!
! Calculation of mass-Inertia Positive Definite Symmetric Matrix !
!================================================================!
! Inverse Mass Identity Matrix
m_I(:,:) = 0.0

m_I(1,1) = 1.0
m_I(2,2) = 1.0
m_I(3,3) = 1.0

! Inertia matrix 
I_p1(:,:) = 0.0
Inv_Ip1(:,:) = 0.0
Inv_Ip2(:,:) = 0.0

I_p1(1,1) = IPXX(IG, IDP)
I_p1(2,2) = IPYY(IG, IDP)
I_p1(3,3) = IPZZ(IG, IDP)

!I_p1 = I_p
I_p2 = I_p1 
!write(*,*) 'I_p' I_p(1,:)

! Transform Inertia Matrix to the global coordinate system
call transform_basis(I_p1, part1%ELLQUAT, shape(I_p1))

call transform_basis(I_p2, part2%ELLQUAT, shape(I_p2))

!write(*,*) 'I_p after transform'

! Invert the Transformed Inertia Matrix
call invert_ndim3_matrix(I_p1, Inv_Ip1)

call invert_ndim3_matrix(I_p2, Inv_Ip2)


!!==== Skew Symmetric Matrix of the Impact Arm vector ====!!
! Particle 1
R_X1(1,1) =  0.0
R_X1(1,2) = -R_impact1(3,1)
R_X1(1,3) =  R_impact1(2,1)

R_X1(2,1) =  R_impact1(3,1) 
R_X1(2,2) =  0.0
R_X1(2,3) = -R_impact1(1,1)

R_X1(3,1) = -R_impact1(2,1)
R_X1(3,2) =  R_impact1(1,1)
R_X1(3,3) =  0.0

R_X_Tr1 = transpose(R_X1)

!write(*,*) 'R_X'

!write(*,*) 'transpose of R_X'

K1 = m_I + matmul(R_X_Tr1, matmul(Inv_Ip1, R_X1))

!write(*,*) 'K1'

!write(*,*) 'Inverse of K'

! Particle 2
R_X2(1,1) =  0.0
R_X2(1,2) = -R_impact2(3,1)
R_X2(1,3) =  R_impact2(2,1)

R_X2(2,1) =  R_impact2(3,1) 
R_X2(2,2) =  0.0
R_X2(2,3) = -R_impact2(1,1)

R_X2(3,1) = -R_impact2(2,1)
R_X2(3,2) =  R_impact2(1,1)
R_X2(3,3) =  0.0

R_X_Tr2 = transpose(R_X2)

!write(*,*) 'R_X'

!write(*,*) 'transpose of R_X'

K2 = m_I + matmul(R_X_Tr2, matmul(Inv_Ip2, R_X2))

!write(*,*) 'K2'

K = K1 + K2

call invert_ndim3_matrix(K, invK)

!write(*,*) 'K'

!write(*,*) 'Inverse of K'

!======================================================================!
!                       Calculation of Impulse                         !
!======================================================================!
!! Impulse for sticking
PC = - matmul(invK, DELU) 

!! Check Friction Cone
! Normal Impulse 
PCN = PC(1,1)*normal(1,1) + PC(2,1)*normal(2,1) + PC(3,1)*normal(3,1)
!write(*,*) PCN

! Tangential Impulse
PCT_X = PC(1,1) - PCN*normal(1,1)
PCT_Y = PC(2,1) - PCN*normal(2,1)
PCT_Z = PC(3,1) - PCN*normal(3,1)

PCT = sqrt(PCT_X*PCT_X + PCT_Y*PCT_Y + PCT_Z*PCT_Z)
!write(*,*) PCT


if (PCT > MUP(IG)*PCN) then

    ! Tangential unit vector
    TW(1,1) = UREL(1,1) - URN*normal(1,1)
    TW(2,1) = UREL(2,1) - URN*normal(2,1)
    TW(3,1) = UREL(3,1) - URN*normal(3,1)

    NRM_TW = sqrt(TW(1,1)*TW(1,1) + TW(2,1)*TW(2,1) + TW(3,1)*TW(3,1))

    if (NRM_TW /= 0.0) then

        TW(1,1) = -TW(1,1)/NRM_TW
        TW(2,1) = -TW(2,1)/NRM_TW
        TW(3,1) = -TW(3,1)/NRM_TW

    end if

    !! Impulse for sliding
    ! Calculate (n + mu*t)
    V_IMP = normal + MUP(IG)*TW 

    ! Calculate K.(n + mu*t)
    K_N = matmul(K, V_IMP)

    ! Dot product of K.(n + mu*t) and n
    KN_dot_N = K_N(1,1)*normal(1,1) + K_N(2,1)*normal(2,1) + K_N(3,1)*normal(3,1)

    ! Dot product of DELU and n
    PCN = -1.0*(DELU(1,1)*normal(1,1) + DELU(2,1)*normal(2,1) + DELU(3,1)*normal(3,1))

    ! Final expression of Impulse
    PC = (PCN/KN_dot_N)*(V_IMP)

    !write(*,*) '       '
    !write(*,*) 'Checked slip condition'
    !write(*,*) 'Norm of Tangential Velocity :', NRM_TW
    !write(*,*) 'Tangential Vector :', TW
    !write(*,*) '       '

end if

!write(*,*) 'Impulse vector dot normal', &
!             (PC(1,1)*normal(1,1) + PC(2,1)*normal(2,1) + PC(3,1)*normal(3,1))/ & 
!             (sqrt(PC(1,1)**2 + PC(2,1)**2 + PC(3,1)**2))
!write(*,*) '       '
!----------------------------------------------------------------!

!############################################################!
!################# ENERGY CALCULATION #######################!
!############################################################!
E_initial_trans1 = 0.5*(U_P1(1,1)*U_P1(1,1) + U_P1(2,1)*U_P1(2,1) + U_P1(3,1)*U_P1(3,1))

Angular_Momentum1 = matmul(I_p1, OMG1)

E_initial_rot1 = 0.5*(OMG1(1,1)*Angular_Momentum1(1,1) &
                    + OMG1(2,1)*Angular_Momentum1(2,1) &
                    + OMG1(3,1)*Angular_Momentum1(3,1))

E_initial1 = E_initial_trans1 + E_initial_rot1

E_initial_trans2 = 0.5*(U_P2(1,1)*U_P2(1,1) + U_P2(2,1)*U_P2(2,1) + U_P2(3,1)*U_P2(3,1))

Angular_Momentum2 = matmul(I_p2, OMG2)

E_initial_rot2 = 0.5*(OMG2(1,1)*Angular_Momentum2(1,1) &
                    + OMG2(2,1)*Angular_Momentum2(2,1) &
                    + OMG2(3,1)*Angular_Momentum2(3,1))

E_initial2 = E_initial_trans2 + E_initial_rot2

E_initial = E_initial1 + E_initial2

!write(*,*) "E_initial_rot : ", E_initial_rot
!write(*,*) "Error in tranformation : ", E_initial_rot - 0.5*(IPXX*omegax*omegax &
!                                                           + IPYY*omegay*omegay &
!                                                           + IPZZ*omegaz*omegaz)
!----------------------------------------------------------------!


!============================================================!
!=================== Update the Velocity ====================!
!============================================================!
! Particle 1
!write(*,*) "Particle 1 "
!write(*,*) "Translational Velocity before Collision : ", part1%UP, part1%VP, part1%WP
!write(*,*) "Angular Velocity before Collision : ", OMG1(1,1), OMG1(2,1), OMG1(3,1)

part1%UP = part1%UP - PC(1,1)
part1%VP = part1%VP - PC(2,1)
part1%WP = part1%WP - PC(3,1)

OMG1 = OMG1 + matmul(Inv_Ip1, matmul(R_X1, -PC))


! Particle 2
!write(*,*) "Particle 2 "
!write(*,*) "Translational Velocity before Collision : ", part2%UP, part2%VP, part2%WP
!write(*,*) "Angular Velocity before Collision : ", OMG2(1,1), OMG2(2,1), OMG2(3,1)

part2%UP = part2%UP + PC(1,1)
part2%VP = part2%VP + PC(2,1)
part2%WP = part2%WP + PC(3,1)

OMG2 = OMG2 + matmul(Inv_Ip2, matmul(R_X2, PC))

!============================================================!
!============================================================!
!============================================================!
!write(*,*) " "
!write(*,*) "Particle 1 "
!write(*,*) "Translational Velocity after Collision : ", part1%UP, part1%VP, part1%WP
!write(*,*) "Angular Velocity after Collision : ", OMG1(1,1), OMG1(2,1), OMG1(3,1)

!write(*,*) "Particle 2 "
!write(*,*) "Translational Velocity after Collision : ", part2%UP, part2%VP, part2%WP
!write(*,*) "Angular Velocity after Collision : ", OMG2(1,1), OMG2(2,1), OMG2(3,1)
!write(*,*) " "
!############################################################!
!################# ENERGY CALCULATION #######################!
!############################################################!

E_final_trans1 = 0.5*(part1%UP*part1%UP + part1%VP*part1%VP + part1%WP*part1%WP)

Angular_Momentum1 = matmul(I_p1, OMG1)

E_final_rot1  = 0.5*(OMG1(1,1)*Angular_Momentum1(1,1) &
                   + OMG1(2,1)*Angular_Momentum1(2,1) &
                   + OMG1(3,1)*Angular_Momentum1(3,1))

E_final1 = E_final_trans1 + E_final_rot1


E_final_trans2 = 0.5*(part2%UP*part2%UP + part2%VP*part2%VP + part2%WP*part2%WP)

Angular_Momentum2 = matmul(I_p2, OMG2)

E_final_rot2 = 0.5*(OMG2(1,1)*Angular_Momentum2(1,1) &
                  + OMG2(2,1)*Angular_Momentum2(2,1) &
                  + OMG2(3,1)*Angular_Momentum2(3,1))

E_final2 = E_final_trans2 + E_final_rot2

E_final = E_final1 + E_final2
!############################################################!

call transform_basis(OMG1, conj_q(part1%ELLQUAT), shape(OMG1))

call transform_basis(OMG2, conj_q(part2%ELLQUAT), shape(OMG2))

!write(*,*) "Particle 1 Angular Velocity after transformation ", OMG1(1,1), OMG1(2,1), OMG1(3,1)

!write(*,*) "Error in tranformation : ", E_final_rot1 - 0.5*(IPXX(IG, IDP)*OMG1(1,1)*OMG1(1,1) &
!                                                          + IPYY(IG, IDP)*OMG1(2,1)*OMG1(2,1) &
!                                                          + IPZZ(IG, IDP)*OMG1(3,1)*OMG1(3,1))

!write(*,*) "Particle 2 Angular Velocity after transformation ", OMG2(1,1), OMG2(2,1), OMG2(3,1)

!write(*,*) "Error in tranformation : ", E_final_rot2 - 0.5*(IPXX(IG, IDP)*OMG2(1,1)*OMG2(1,1) &
!                                                          + IPYY(IG, IDP)*OMG2(2,1)*OMG2(2,1) &
!                                                          + IPZZ(IG, IDP)*OMG2(3,1)*OMG2(3,1))
!write(*,*) " "

part1%OMEGAX = OMG1(1,1)
part1%OMEGAY = OMG1(2,1)
part1%OMEGAZ = OMG1(3,1)

part2%OMEGAX = OMG2(1,1)
part2%OMEGAY = OMG2(2,1)
part2%OMEGAZ = OMG2(3,1)

if(abs((E_initial- E_final)/E_initial) > 1.0e-9) then
    !ERROR_COUNT = ERROR_COUNT + 1

    !AVG_ERROR = AVG_ERROR + abs((E_initial- E_final)/E_initial)
    !write(*,*) '   '
    write(*,*) ' Error in Particle Collision '
    write(*,*) ' BALANCE OF KINETIC ENERGY'
    write(*,*) '-----------------------------------------------------------'    
    print*, ' KE Initial', E_initial
    print*, ' KE Final  ', E_final
    print*, ' Error in KE (Total)%', abs((E_initial- E_final)/E_initial)*100.0
    write(*,*) '-----------------------------------------------------------'
    print*, ' KE Final Translational', E_final_trans1 + E_final_trans2
    print*, ' KE Final Rotational', E_final_rot1 + E_final_rot2
    !print*, ' KE Final Total', E_final_trans + E_final_rot
    write(*,*) '-----------------------------------------------------------'
    !write(*,*) 'Particle 1 information'
    !write(*,*) 'X = ', part1%XP, 'Y = ', part1%YP, ' Z =', part1%ZP
    !write(*,*) 'Particle 2 information'
    !write(*,*) 'X = ', part2%XP, 'Y = ', part2%YP, ' Z =', part2%ZP
    !write(*,*) 'Penetration =', sqrt((point1(1,1) - point2(1,1))**2 + (point1(2,1) - point2(2,1))**2 + (point1(3,1) - point2(3,1))**2)
    !write(*,*) 'Common normal', normal
    !write(*,*) 'Point 1', point1
    !write(*,*) 'Point 2', point2
    !write(*,*) 'Point1 - Point2', point1 - point2

    !write(*, *) 'Collision count =',  COLLISION_COUNT 
    stop
end if 
!stop
end subroutine COLLISION_RESPONSE