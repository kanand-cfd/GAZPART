!!=====================================================================
!! Variables for checking CPU's time
!!=====================================================================
module CHECK_CPU

implicit none

!!- Full CPU time elapsed since the beginning
real(kind=8) :: CPU_ELAPSED


!!- Full CPU time for a time cycle
real(kind=8) :: CPU_CYCLE


!!- CPU time for the fluid resolution
real(kind=8), dimension(8) :: CPU_FLUID
!! CPU_FLUID(1) --> All fluid steps
!! CPU_FLUID(2) --> Non-linear terms
!! CPU_FLUID(3) --> Time-advancing 
!! CPU_FLUID(4) --> Momentum forcing
!! CPU_FLUID(5) --> Projection on soloneidal
!! CPU_FLUID(6) --> Statistics
!! CPU_FLUID(7) --> Save fluid
!! CPU_FLUID(8) --> FFT


!!- CPU time for the particle tracking
real(kind=8), dimension(7) :: CPU_PART
!! CPU_PART(1) --> All particle steps
!! CPU_PART(2) --> Position Time-advancing 
!! CPU_PART(3) --> Velocity Time-advancing 
!! CPU_PART(4) --> Interpolation
!! CPU_PART(5) --> Boundary condition (including particle passing)
!! CPU_PART(6) --> Statistics
!! CPU_PART(7) --> Save particles




!!- CPU time for run initiation
real(kind=8) :: CPU_INITIATION

!!- CPU time for run completion
real(kind=8) :: CPU_FINISH

end module CHECK_CPU


