!!====================================================================
!!
!! 
!! This routine build the Right Hand Side of Navier-Stokes equations
!!
!!
!!====================================================================

subroutine NONLINEAR_TERMS(   UFLU,    VFLU,    WFLU, &
                              UFOU,    VFOU,    WFOU, &
                          RHS_UFOU,RHS_VFOU,RHS_WFOU, &
                                                 IRK  )

!!====================================================================
!!
!!  1- compute the vorticity in Fourier space (VORTX,VORTY,VORTZ)
!!   and FFT-1 to get vorticity in physical space (ZETAX,ZETAY,ZETAZ)
!!
!!  2- Build non-linear terms in physical space (TX,TY,TZ)
!!
!!  3- Transform non-linear terms in Fourier space (VORTX,VORTY,VORTZ)
!!
!!  4- Build RHS term according to soloneidal basis projection 
!!     insuring the divergence free velocity field
!!
!!
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE, only: KX, KY, KZ, FILTER
use PARAM_PHYS,  only: VISC
use CHECK_CPU
use STATISTICS, only:EPS_FLU

use P3DFFT

implicit none


!------------------------------------------------------------------
! ARRAYS STATEMENT
!------------------------------------------------------------------
real(kind=8), dimension(ISTART(1)         :IEND(1)             &
                       ,ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL    &
                       ,ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL ), intent(in) :: UFLU

real(kind=8), dimension(ISTART(1)         :IEND(1)             &
                       ,ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL    &
                       ,ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL ), intent(in) :: VFLU

real(kind=8), dimension(ISTART(1)         :IEND(1)             &
                       ,ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL    &
                       ,ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL ), intent(in) :: WFLU

!- Fluid velocity is Fourier space
double complex, &
   dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)), intent(in) :: UFOU
double complex, &
   dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)), intent(in) :: VFOU
double complex, &
   dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)), intent(in) :: WFOU

!- Right-Hand-Side
double complex, &
   dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)), intent(inout) :: RHS_UFOU
double complex, &
   dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)), intent(inout) :: RHS_VFOU
double complex, &
   dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)), intent(inout) :: RHS_WFOU


!- Runge-Kutta cycle
integer, intent(in) :: IRK
!------------------------------------------------------------------

!- Vorticity
double complex, &
   dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)) :: VORTX
double complex, &
   dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)) :: VORTY
double complex, &
   dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)) :: VORTZ

real(kind=8),   &
   dimension(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)) :: ZETAX
real(kind=8),   &
   dimension(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)) :: ZETAY
real(kind=8),   &
   dimension(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)) :: ZETAZ

real(kind=8),   &
   dimension(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)) :: TX
real(kind=8),   &
   dimension(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)) :: TY
real(kind=8),   &
   dimension(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)) :: TZ
!- squared wavenumber
real(kind=8) :: KAPPA2

!- Time control variable
real(kind=8) :: TIME_START, TIME_END

!- Index
integer :: I, J, K,IJK
!------------------------------------------------------------------

!!- CPU check
if(MYID == 0) then
 TIME_START=MPI_WTIME()
end if


!!=================================================================
!! 1. Vorticity
!!=================================================================

!!-----------------------------------------------------------------
!! 1.1. Vorticity in Fourier space
!!-----------------------------------------------------------------
!! [vort_x]f = i*ky*[wf]fou - i*kz*[vf]fou

do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
   VORTX(I,J,K) = ICMPL*WFOU(I,J,K)*KY(J) - ICMPL*VFOU(I,J,K)*KZ(K)
   VORTY(I,J,K) = ICMPL*UFOU(I,J,K)*KZ(K) - ICMPL*WFOU(I,J,K)*KX(I)
   VORTZ(I,J,K) = ICMPL*VFOU(I,J,K)*KX(I) - ICMPL*UFOU(I,J,K)*KY(J)
  end do
 end do
end do



!!-----------------------------------------------------------------
!! 1.2. Vorticity back in physical space
!!-----------------------------------------------------------------
!- Synchronize all the process
 !call MPI_BARRIER(MPI_COMM_WORLD,IERR)

!- Back in physical space
 call P3DFFT_BTRAN_C2R(VORTX,ZETAX,FFTFLAG)       
 call P3DFFT_BTRAN_C2R(VORTY,ZETAY,FFTFLAG)
 call P3DFFT_BTRAN_C2R(VORTZ,ZETAZ,FFTFLAG)


!!-----------------------------------------------------------------
!! 1.3. Statistics on vorticity 
!!-----------------------------------------------------------------
if(IRK==1) then
!!- Dissipation
 EPS_FLU = 0.
 EPS_FLU = sum(ZETAX**2) + sum(ZETAY**2) + sum(ZETAZ**2)
 EPS_FLU = VISC*EPS_FLU
end if

!!=================================================================
!! 2. Non-linear terms
!!=================================================================
TX(:,:,:) = VFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))*ZETAZ(:,:,:) &
          - WFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))*ZETAY(:,:,:)

TY(:,:,:) = WFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))*ZETAX(:,:,:) &
          - UFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))*ZETAZ(:,:,:)

TZ(:,:,:) = UFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))*ZETAY(:,:,:) &
          - VFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))*ZETAX(:,:,:)

 

!!=================================================================
!! 3. Non-linear terms --> Fourier space
!!=================================================================
!- Synchronize all the process
! call MPI_BARRIER(MPI_COMM_WORLD,IERR)

call P3DFFT_FTRAN_R2C(TX,VORTX,FFTFLAG)
call P3DFFT_FTRAN_R2C(TY,VORTY,FFTFLAG)
call P3DFFT_FTRAN_R2C(TZ,VORTZ,FFTFLAG)


!- FFT Normalization
VORTX(:,:,:) = VORTX(:,:,:)*FACTOR
VORTY(:,:,:) = VORTY(:,:,:)*FACTOR
VORTZ(:,:,:) = VORTZ(:,:,:)*FACTOR


!!- Warning -------------------------------------------------------
!!  The arrays VORTX, VORTY, VORTZ are used to store the nonlinear
!!  terms in Fourier space instead of the vorticity in Fourier space
!!- Warning -------------------------------------------------------


!!=================================================================
!! 4 RHS Non-linear terms
!!=================================================================
do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)

   KAPPA2 = KX(I)**2 + KY(J)**2 + KZ(K)**2
   if (KAPPA2 <= ZERO) KAPPA2 = INFINITY

   !- x-component
   RHS_UFOU(I,J,K) =   RHS_UFOU(I,J,K) &
                     + (KY(J)**2+ KZ(K)**2)/KAPPA2 * VORTX(I,J,K) &
                     -  KX(I)   * KY(J)    /KAPPA2 * VORTY(I,J,K) &
                     -  KX(I)   * KZ(K)    /KAPPA2 * VORTZ(I,J,K)

   !- y-component
   RHS_VFOU(I,J,K) =   RHS_VFOU(I,J,K) &
                     -  KY(J)   * KX(I)    /KAPPA2 * VORTX(I,J,K) &
                     + (KX(I)**2+ KZ(K)**2)/KAPPA2 * VORTY(I,J,K) &
                     -  KY(J)   * KZ(K)    /KAPPA2 * VORTZ(I,J,K)

   !- z-component
   RHS_WFOU(I,J,K) =   RHS_WFOU(I,J,K) &
                     -  KZ(K)   * KX(I)    /KAPPA2 * VORTX(I,J,K) &
                     -  KZ(K)   * KY(J)    /KAPPA2 * VORTY(I,J,K) &
                     + (KX(I)**2+ KY(J)**2)/KAPPA2 * VORTZ(I,J,K) 

  end do
 end do
end do




!!- CPU check
if(MYID == 0) then
 TIME_END=MPI_WTIME()
 CPU_FLUID(2) = CPU_FLUID(2) + TIME_END - TIME_START
end if


end subroutine NONLINEAR_TERMS
