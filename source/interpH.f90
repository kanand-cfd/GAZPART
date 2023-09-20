!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
!                            INTERPOLATION                              !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
! The Eulerian velocity field "UAIN" located on the mesh "XAIN,YAIN,    !
! ZAIN" is interpolated on a Larangian position "XINT,YINT,ZINT" and    !
! the result is "UINT".                                                 !
!                                                                       !
!-----------------------------------------------------------------------!
!                                                                       !
!    UAIN(XDIM,YDIM,ZDIM) \            / UINT(NINT)                     !
!    XAIN(XDIM)            \          /  XINT(NINT)                     !
!    YAIN(YDIM)            /  ======> \  YINT(NINT)                     !
!    ZAIN(ZDIM)           /            \ ZINT(NINT)                     !
!                                                                       !
! 3 interpolation schemes are avaible :                                 !
!    INT_SCHEME  = 1 : 1st order Lagrangian polynomial                  !
!                = 2 : 2nd order Lagrangian polynomial                  !
!                = 3 : 3rd order Lagrangian polynomial                  !
!                = 4 : 4th order Lagrangian polynomial                  !
!                                                                       !
!=======================================================================!
!                                                                       !
subroutine INTERPH(INT_SCHEME,     & !- User interpolation scheme choice!
                   XAIN,YAIN,ZAIN, & !- Mesh for interpolation          !
                   UAIN,           & !- Velocity field for interpolation!
                   NINT,           & !- Interpolated array size         !
                   XINT,YINT,ZINT, & !- Interpolated Positions          !
                   UINT            ) !- Interpolated velocity field     !
!                                                                       !
!=======================================================================!

use DNS_DIM
use CHECK_CPU

implicit none

!=======================================================================!
! Input arrays
!=============
!- User interpolation scheme choice
integer, intent(in) :: INT_SCHEME

!- Mesh of data for interpolation
real(kind=8), dimension(ISTART(1):IEND(1)), intent(in) :: XAIN
real(kind=8), dimension(ISTART(2):IEND(2)), intent(in) :: YAIN
real(kind=8), dimension(ISTART(3):IEND(3)), intent(in) :: ZAIN


!- Data field for interpolation
real(kind=8),                                         &
   dimension(ISTART(1)         :IEND(1)               &
            ,ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL      &
            ,ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL ) ,  &
   intent(in) :: UAIN


!- Interpolated array size
integer, intent(in) :: NINT


!- Interpolated Positions 
real(kind=8), dimension(NINT), intent(in) :: XINT
real(kind=8), dimension(NINT), intent(in) :: YINT
real(kind=8), dimension(NINT), intent(in) :: ZINT


! Ouput array
!=============
!
!- Interpolated velocity field
real(kind=8), dimension(NINT), intent(out) :: UINT


!- Time control variable
real(kind=8) :: TIME_START, TIME_END


external INTERP_LAG1, INTERP_LAG2, INTERP_LAG3
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!- CPU check
if(MYID == 0) then
 TIME_START=MPI_WTIME()
end if



!!====================================================================
!! 1. 1st order interpolation scheme
!!====================================================================
if (INT_SCHEME == 1) then

 call INTERP_LAG1( XAIN,YAIN,ZAIN, &
                   UAIN,           &
                   NINT,           &
                   XINT,YINT,ZINT, &
                   UINT            )

!!====================================================================
!! 2. 2nd order interpolation scheme
!!====================================================================
else if (INT_SCHEME == 2) then

 call INTERP_LAG2( XAIN,YAIN,ZAIN, &
                   UAIN,           &
                   NINT,           &
                   XINT,YINT,ZINT, &
                   UINT            )


!!====================================================================
!! 3. 3rd order interpolation scheme
!!====================================================================
else if (INT_SCHEME == 3) then 

 call INTERP_LAG3( XAIN,YAIN,ZAIN, &
                   UAIN,           &
                   NINT,           &
                   XINT,YINT,ZINT, &
                   UINT            )
!!====================================================================
!! 3. 4th order interpolation scheme
!!====================================================================
else if (INT_SCHEME == 4) then 

 call INTERP_LAG4( XAIN,YAIN,ZAIN, &
                   UAIN,           &
                   NINT,           &
                   XINT,YINT,ZINT, &
                   UINT            )
end if



!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()
 CPU_PART(4) = CPU_PART(4) + TIME_END - TIME_START
end if

end subroutine INTERPH
