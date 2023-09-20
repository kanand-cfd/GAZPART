subroutine QUATERNION_ANGULAR_VELOCITY_INTGERATION(I, IG)

!!
!!
!!

use PARAM_PHYS
use DNS_DIM
use PARTICLE_PARALLEL
use mod_quaternion

implicit none

integer, intent(in) :: I, IG

integer :: IDP

real(kind=8), dimension(ndim, 1) :: FTORQUE

real(kind=8), dimension(ndim,1) :: omg_n_4, omg_n_2, omg_dt_n, omg_dt_n_2

type(quaternion) :: quat_n_2, quat_n, quat_map

real(kind=8) :: q0, q1, q2, q3, norm_omg

!=============================================!
IDP = PART(I,IG)%IDP


!! Compute Torque at time step n 
if(SOLVE_FLUID > 0 .or. FROZEN_FLOW > 0) then

    call FLUID_TORQUE_ELLIPSOID(I, IG, IDP,        &
                                PART(I,IG)%OMEGAX, &
                                PART(I,IG)%OMEGAY, &
                                PART(I,IG)%OMEGAZ, &
                                FTORQUE)

else

    FTORQUE(:,1) = 0.0


end if


omg_dt_n(1,1) = (FTORQUE(1,1) + PART(I,IG)%OMEGAY*PART(I,IG)%OMEGAZ*(IPYY(IG,IDP) - IPZZ(IG,IDP)))/IPXX(IG, IDP)

omg_dt_n(2,1) = (FTORQUE(2,1) + PART(I,IG)%OMEGAZ*PART(I,IG)%OMEGAX*(IPZZ(IG,IDP) - IPXX(IG,IDP)))/IPYY(IG, IDP)

omg_dt_n(3,1) = (FTORQUE(3,1) + PART(I,IG)%OMEGAX*PART(I,IG)%OMEGAY*(IPXX(IG,IDP) - IPYY(IG,IDP)))/IPZZ(IG, IDP)

! Angular velocity at time step n+1/4
omg_n_4(1,1) = PART(I,IG)%OMEGAX + (1.0/4.0)*DTIME*omg_dt_n(1,1)
omg_n_4(2,1) = PART(I,IG)%OMEGAY + (1.0/4.0)*DTIME*omg_dt_n(2,1)
omg_n_4(3,1) = PART(I,IG)%OMEGAZ + (1.0/4.0)*DTIME*omg_dt_n(3,1)
    

! Angular velocity at time step n+1/2
omg_n_2(1,1) = PART(I,IG)%OMEGAX + (1.0/2.0)*DTIME*omg_dt_n(1,1)
omg_n_2(2,1) = PART(I,IG)%OMEGAY + (1.0/2.0)*DTIME*omg_dt_n(2,1)
omg_n_2(3,1) = PART(I,IG)%OMEGAZ + (1.0/2.0)*DTIME*omg_dt_n(3,1)


! Transform Angular velocity at (n+1/4) to fixed frame
call transform_basis(omg_n_4, PART(I,IG)%ELLQUAT, shape(omg_n_4))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine Quaternion at (n+1/2) using fixed Angular velocity at (n+1/4)

! Calculate the norm of Angular velocity
norm_omg = sqrt(omg_n_4(1,1)*omg_n_4(1,1) + omg_n_4(2,1)*omg_n_4(2,1) + omg_n_4(3,1)*omg_n_4(3,1))
!write(*,*) norm_omg

if (norm_omg /= 0) then
    ! Calculation of exponential map converted to quaternion
    q0 = cos((norm_omg*DTIME)/4.0)

    q1 = sin((norm_omg*DTIME)/4.0)*(omg_n_4(1,1)/norm_omg) 

    q2 = sin((norm_omg*DTIME)/4.0)*(omg_n_4(2,1)/norm_omg)

    q3 = sin((norm_omg*DTIME)/4.0)*(omg_n_4(3,1)/norm_omg)

else

    q0 = 1.0
    q1 = 0.0
    q2 = 0.0
    q3 = 0.0

end if ! if (norm_omg /= 0) then


! Update the quaternion map  
quat_map%a = q0
quat_map%b = q1
quat_map%c = q2
quat_map%d = q3

! Update the Quaternion at next time-step by direct multiplication
quat_n_2 = mult_quat(quat_map, PART(I,IG)%ELLQUAT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Compute Torque at time step n + 1/2 
!! Compute Torque at time step n 
if(SOLVE_FLUID > 0 .or. FROZEN_FLOW > 0) then

    call FLUID_TORQUE_ELLIPSOID(I, IG, IDP,        &
                                PART(I,IG)%OMEGAX, &
                                PART(I,IG)%OMEGAY, &
                                PART(I,IG)%OMEGAZ, &
                                FTORQUE)

else

    FTORQUE(:,1) = 0.0


end if


omg_dt_n_2(1,1) = (FTORQUE(1,1) + omg_n_2(2,1)*omg_n_2(3,1)*(IPYY(IG,IDP) - IPZZ(IG,IDP)))/IPXX(IG, IDP)

omg_dt_n_2(2,1) = (FTORQUE(2,1) + omg_n_2(3,1)*omg_n_2(1,1)*(IPZZ(IG,IDP) - IPXX(IG,IDP)))/IPYY(IG, IDP)

omg_dt_n_2(3,1) = (FTORQUE(3,1) + omg_n_2(1,1)*omg_n_2(2,1)*(IPXX(IG,IDP) - IPYY(IG,IDP)))/IPZZ(IG, IDP)


! Transform Angular velocity at (n+1/2) to fixed frame using quat at (n+1/2)
call transform_basis(omg_n_2, quat_n_2, shape(omg_n_2))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine corrected Quaternion at (n+1) using fixed Angular velocity at (n+1/2)

! Calculate the norm of Angular velocity
norm_omg = sqrt(omg_n_2(1,1)*omg_n_2(1,1) + omg_n_2(2,1)*omg_n_2(2,1) + omg_n_2(3,1)*omg_n_2(3,1))
!write(*,*) norm_omg

if (norm_omg /= 0) then
    ! Calculation of exponential map converted to quaternion
    q0 = cos((norm_omg*DTIME)/2.0)

    q1 = sin((norm_omg*DTIME)/2.0)*(omg_n_2(1,1)/norm_omg) 

    q2 = sin((norm_omg*DTIME)/2.0)*(omg_n_2(2,1)/norm_omg)

    q3 = sin((norm_omg*DTIME)/2.0)*(omg_n_2(3,1)/norm_omg)

else

    q0 = 1.0
    q1 = 0.0
    q2 = 0.0
    q3 = 0.0

end if ! if (norm_omg /= 0) then


! Update the quaternion map  
quat_map%a = q0
quat_map%b = q1
quat_map%c = q2
quat_map%d = q3

! Update the Quaternion at next time-step by direct multiplication
quat_n = mult_quat(quat_map, PART(I,IG)%ELLQUAT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PART(I, IG)%ELLQUAT = quat_n

PART(I,IG)%OMEGAX = PART(I,IG)%OMEGAX + DTIME*omg_dt_n_2(1,1)
PART(I,IG)%OMEGAY = PART(I,IG)%OMEGAY + DTIME*omg_dt_n_2(2,1)
PART(I,IG)%OMEGAZ = PART(I,IG)%OMEGAZ + DTIME*omg_dt_n_2(3,1)

!if(I==1) write(*,*) 'OMEGA ', PART(I,IG)%OMEGAX, PART(I,IG)%OMEGAY, PART(I,IG)%OMEGAZ
!if(I==1) write(*,*) 'FTORQUE ',FTORQUE
!if(I==1) write(*,*) 'norm quat ', norm_q(PART(I,IG)%ELLQUAT) 


if((PART(I,IG)%ELLQUAT%a /= PART(I,IG)%ELLQUAT%a) .or. (PART(I,IG)%OMEGAX /= PART(I,IG)%OMEGAX)) then
    
    write(*,*) I, IG
    write(*,*) PART(I,IG)%OMEGAX, PART(I,IG)%OMEGAY, PART(I,IG)%OMEGAZ
    write(*,*) PART(I,IG)%ELLQUAT
    stop 'At the end of time-step'

end if 

return
end subroutine QUATERNION_ANGULAR_VELOCITY_INTGERATION