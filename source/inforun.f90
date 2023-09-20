!!====================================================================
!!
!! This routine print informations about the simulation in file
!!                           'run.info'
!!
!!====================================================================

subroutine INFORUN

!!====================================================================
!! 
!!====================================================================

use DNS_DIM
use PARAM_PHYS 
use FORCING
use GEOMETRIC_VARIABLE
use ENSIGHT_VAR
use STATISTICS
use PARTICLE_PARALLEL

implicit none

!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- index
integer :: I, J, N, IDP
!---------------------------------------------------------------------




!!=====================================================================
!! RUN INFORMATIONS
!!=====================================================================
I = 1

write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'====================================================================='
write(UNIT_INFO(I),*)' EQUATIONS SOLVED DURING THE SIMULATION :'
write(UNIT_INFO(I),*)'====================================================================='
if(SOLVE_FLUID>0)write(UNIT_INFO(I),*)'   --> Fluid momentum'
if(SOLVE_SCALAR) write(UNIT_INFO(I),*)'   --> Passive scalar transport'
if(SOLVE_PART  ) write(UNIT_INFO(I),*)'   --> Particle momentum'
if(SOLVE_PART  ) write(UNIT_INFO(I),*)'   --> Particle position'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'====================================================================='
write(UNIT_INFO(I),*)'NUMERICAL PARAMETERS'
write(UNIT_INFO(I),*)'====================================================================='
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),1603)'CPU Number = ',NPROC
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),1608)'Decomposition, NDIM = ',NDIM 
write(UNIT_INFO(I),1608)'              IPROC = ',IPROC
write(UNIT_INFO(I),1608)'              JPROC = ',JPROC
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'====================================================================='
write(UNIT_INFO(I),*)'GEOMETRICS'
write(UNIT_INFO(I),*)'====================================================================='
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Mesh'
write(UNIT_INFO(I),*)'----'
write(UNIT_INFO(I),*)'Physical space:'
write(UNIT_INFO(I),1600)'Nx = ',NX,',  Ny = ',NY,',  Nz = ',NZ
write(UNIT_INFO(I),1601)'Lx = ',LXMAX,', Lx = ',LYMAX,', Lx = ',LZMAX
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),1603)'For CPU:',MYID
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)' Physical space:'
write(UNIT_INFO(I),1600)'min(I) = ',ISTART(1),', min(J) = ',ISTART(2),', min(K) = ',ISTART(3)
write(UNIT_INFO(I),1600)'max(I) = ',IEND(1),', max(J) = ',IEND(2),', max(K) = ',IEND(3)
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),1601)'min(x) = ',XMESH(ISTART(1)),',  max(x) =',XMESH(IEND(1))
write(UNIT_INFO(I),1601)'min(y) = ',YMESH(ISTART(2)),',  max(y) =',YMESH(IEND(2))
write(UNIT_INFO(I),1601)'min(z) = ',ZMESH(ISTART(3)),',  max(z) =',ZMESH(IEND(3))
write(UNIT_INFO(I),*)
if(SOLVE_FLUID>0) then 
write(UNIT_INFO(I),*)'Fourier space:'
write(UNIT_INFO(I),1600)'min(I) = ',FSTART(1),', min(J) = ',FSTART(2),', min(K) = ',FSTART(3)
write(UNIT_INFO(I),1600)'max(I) = ',  FEND(1),', max(J) = ',  FEND(2),', max(K) = ',  FEND(3)
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),1601)'Kx(1) = ',KX(FSTART(1)),'  max(kx) = ',KX(FEND(1))
write(UNIT_INFO(I),1601)'Ky(1) = ',KY(FSTART(2)),'  max(ky) = ',KY(FEND(2))
write(UNIT_INFO(I),1601)'Kz(1) = ',KZ(FSTART(3)),'  max(kz) = ',KZ(FEND(3))
end if
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Time control'
write(UNIT_INFO(I),*)'------------'
write(UNIT_INFO(I),1602)'              dt = ',DTIME_USER
write(UNIT_INFO(I),1603)'Cycle number max = ',NCYCLEMAX
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)


write(UNIT_INFO(I),*)'====================================================================='
write(UNIT_INFO(I),*)'INITIATION'
write(UNIT_INFO(I),*)'====================================================================='
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Fluid:'
write(UNIT_INFO(I),*)'------'
if(INIT_FLUID_VELOCITY == 0) then
 write(UNIT_INFO(I),*) 'Initiation --> Zero'
elseif(INIT_FLUID_VELOCITY == 1) then
 write(UNIT_INFO(I),*) 'Initiation --> Single Eddy'
elseif(INIT_FLUID_VELOCITY == 2) then
 write(UNIT_INFO(I),*) 'Initiation --> Random velocity field'
elseif(INIT_FLUID_VELOCITY == 3) then
 write(UNIT_INFO(I),*) 'Initiation --> Read in file: FLUID.ini'
 if (STEADY) write(UNIT_INFO(I),*) '                             FORCE.ini'
elseif(INIT_FLUID_VELOCITY == 4) then
 write(UNIT_INFO(I),*) 'Initiation --> Start from a stored fluid velocity field'
else
 write(UNIT_INFO(I),*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 write(UNIT_INFO(I),*)'!!         WRONG PARAMETER         !!'
 write(UNIT_INFO(I),*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 write(UNIT_INFO(I),*)'!! file   : inforun.f90            !!'
 write(UNIT_INFO(I),*)'!! problem: INIT_FLUID_VELOCITY>4  !!'
 write(UNIT_INFO(I),*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 stop
end if

if(SOLVE_SCALAR) then
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Scalar:'
write(UNIT_INFO(I),*)'-------'
if(INIT_SCALAR == 0 .or. INIT_SCALAR == 1) then
 write(UNIT_INFO(I),*) 'Initiation --> Zero'
elseif(INIT_SCALAR == 2) then
 write(UNIT_INFO(I),*) 'Initiation --> Random value'
elseif(INIT_SCALAR == 3) then
 write(UNIT_INFO(I),*) 'Initiation --> Read in file: SCALAR.ini'
end if
end if


if(SOLVE_PART) then
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Particles:'
write(UNIT_INFO(I),*)'----------'
if(INTERP_SCHEME == 1) then
 write(UNIT_INFO(I),*) '1st order Lagrangian interpolation'
elseif(INTERP_SCHEME == 2) then
 write(UNIT_INFO(I),*) '2nd order Lagrangian interpolation'
elseif(INTERP_SCHEME == 3) then
 write(UNIT_INFO(I),*) '3rd order Lagrangian interpolation'
elseif(INTERP_SCHEME == 4) then
 write(UNIT_INFO(I),*) 'SFM interpolation scheme'
else
 write(UNIT_INFO(I),*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 write(UNIT_INFO(I),*)'!!         WRONG PARAMETER         !!'
 write(UNIT_INFO(I),*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 write(UNIT_INFO(I),*)'!! file   : inforun.f90            !!'
 write(UNIT_INFO(I),*)'!! problem: INTERP_SCHEME>4        !!'
 write(UNIT_INFO(I),*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 stop
end if 

if(INIT_PART_POSITION == 0.or.INIT_PART_POSITION == 1) then
 write(UNIT_INFO(I),*) 'Initial particle position --> Uniformaly distributed'
elseif(INIT_PART_POSITION == 2) then
 write(UNIT_INFO(I),*) 'Initial particle position --> Random'
elseif(INIT_PART_POSITION == 3) then
 write(UNIT_INFO(I),*) 'Initial particle position --> Read in file: part.ini'
elseif(INIT_PART_POSITION == 4) then
 write(UNIT_INFO(I),*) 'Initial particle position --> Injection at edge (x,y)'
else
 write(UNIT_INFO(I),*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 write(UNIT_INFO(I),*)'!!         WRONG PARAMETER         !!'
 write(UNIT_INFO(I),*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 write(UNIT_INFO(I),*)'!! file   : inforun.f90            !!'
 write(UNIT_INFO(I),*)'!! problem: INIT_PART_POSITION>3   !!'
 write(UNIT_INFO(I),*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 stop
end if 

if(INIT_PART_VELOCITY == 0) then
 write(UNIT_INFO(I),*) 'Initial particle velocity --> Zero'
elseif(INIT_PART_VELOCITY == 1) then
 write(UNIT_INFO(I),*) 'Initial particle velocity --> Equal to fluid'
elseif(INIT_PART_VELOCITY == 2) then
 write(UNIT_INFO(I),*) 'Initial particle velocity --> Random'
elseif(INIT_PART_VELOCITY == 3) then
 write(UNIT_INFO(I),*) 'Initial particle velocity --> Read in file: part.ini'
elseif(INIT_PART_VELOCITY == 4) then
 write(UNIT_INFO(I),*) 'Initial particle velocity --> Injection at the edge (x,y)'
else
 write(UNIT_INFO(I),*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 write(UNIT_INFO(I),*)'!!         WRONG PARAMETER         !!'
 write(UNIT_INFO(I),*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 write(UNIT_INFO(I),*)'!! file   : inforun.f90            !!'
 write(UNIT_INFO(I),*)'!! problem: INIT_PART_VELOCITY>4   !!'
 write(UNIT_INFO(I),*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 stop
end if 

if(SOLVE_SCALAR) then
if(INIT_PART_SCALAR == 0) then
 write(UNIT_INFO(I),*) 'Initial particle velocity --> Zero'
elseif(INIT_PART_SCALAR == 1) then
 write(UNIT_INFO(I),*) 'Initial particle velocity --> Equal to fluid'
elseif(INIT_PART_SCALAR == 2) then
 write(UNIT_INFO(I),*) 'Initial particle velocity --> Random ... not done'
elseif(INIT_PART_SCALAR == 3) then
 write(UNIT_INFO(I),*) 'Initial particle velocity --> Read in file: traj.ini'
else
 write(UNIT_INFO(I),*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 write(UNIT_INFO(I),*)'!!         WRONG PARAMETER         !!'
 write(UNIT_INFO(I),*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 write(UNIT_INFO(I),*)'!! file   : inforun.f90            !!'
 write(UNIT_INFO(I),*)'!! problem: INIT_PART_VELOCITY>3   !!'
 write(UNIT_INFO(I),*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 stop
end if
end if 

end if

write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'====================================================================='
write(UNIT_INFO(I),*)'PHYSICAL PARAMETERS'
write(UNIT_INFO(I),*)'====================================================================='
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Fluid properties:'
write(UNIT_INFO(I),*)'-----------------'
write(UNIT_INFO(I),*)'     Molecular viscosity,  nu = ',VISC,' m2/s'
write(UNIT_INFO(I),*)'                 Density, rho = ',RHOF,' kg/m3'
write(UNIT_INFO(I),*)'     Dynamical viscosity,     = ',VISC*RHOF,' kg/m/s'
if(STEADY) then
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Turbulence --> Forced'
write(UNIT_INFO(I),1608)' --  Number of forced Wavenumber = ',NFORCE_FULL
write(UNIT_INFO(I),1604)' --      First forced Wavenumber = ',KFORCE_MIN/KFIRST,'*K0'
write(UNIT_INFO(I),1604)' --       Last forced Wavenumber = ',KFORCE_MAX/KFIRST,'*K0'
write(UNIT_INFO(I),1605)' --    Stochastic force variance = ',SIGMA_FORCE
write(UNIT_INFO(I),1605)' --   Stochastic force timescale = ',TIME_FORCE
write(UNIT_INFO(I),1605)' --                        dt/tf = ',DTIME_USER/TIME_FORCE
else
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Turbulence --> decaying'
end if

if(SOLVE_SCALAR) then
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Scalar properties:'
write(UNIT_INFO(I),*)'-----------------'
write(UNIT_INFO(I),1612)'       Diffusivity,     K = ',DIFF_SCL,' [m2/s]'
write(UNIT_INFO(I),1612)'   Specific "heat",    Cp = ',CP_SCL,' [J/kg/K]'
write(UNIT_INFO(I),1612)'  Imposed gradient, dT/dy = ',GRAD_SCL,' [kg/m/s]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),1605)'     Prandtl (or Schmidt) = ',VISC/DIFF_SCL
write(UNIT_INFO(I),*)


end if

if(SOLVE_PART) then
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Particle properties:'
write(UNIT_INFO(I),*)'--------------------'
do J = 1, NIG
 if(PARTDEF(J) == 0) then
  write(UNIT_INFO(I),1606)'Particle:',J,'  ==> Motionless'
  write(UNIT_INFO(I),*)
 elseif(PARTDEF(J) == 1) then
  write(UNIT_INFO(I),1606)'Particle:',J,'  ==> Fluid element'
  write(UNIT_INFO(I),*)
 elseif(PARTDEF(J) == 2) then
  if(POLYDISP(J) ==1) then
   write(UNIT_INFO(I),1606)'Particle:',J,'  ==> Monodisperse particles'
   write(UNIT_INFO(I),1609) '--    Np = ',NPART_FULL
   write(UNIT_INFO(I),1612) '--  ep_a = ',EMAJ_PART_USER(J),' [m]'
   write(UNIT_INFO(I),1612) 'Aspect Ratio = ',APR_PART(J) 
   write(UNIT_INFO(I),1612) '--  rhop = ',RHOP_USER(J),' [kg/m3]'
   write(UNIT_INFO(I),1612) '--   alp = ',(4.0*PPI/3.0)*(NPART_FULL*EMAJ_PART_USER(J)**3)/(LXMAX*LYMAX*LZMAX*APR_PART(J)**2),' [-]'
   !write(UNIT_INFO(I),1612) '--  taup = ',RHOP_USER(J)*DPART_USER(J)**2/18./(VISC*RHOF),' [s]'
  else
   write(UNIT_INFO(I),1606)'Particle:',J,'  ==> Polydisperse/Polysolid particles'
   write(UNIT_INFO(I),1609)'--     Np = ',NPART_FULL
   write(UNIT_INFO(I),1609)'-- Nclass = ',POLYDISP(J)
   write(UNIT_INFO(I),*)
   do IDP = 1, POLYDISP(J)
    write(UNIT_INFO(I),1614)'-- ->   Np(',IDP,') = ',NPCLASS(J,IDP),' [-]'
    !write(UNIT_INFO(I),1613)'-- ->   dp(',IDP,') = ',DPART(J,IDP),' [m]'
    write(UNIT_INFO(I),1615)'-- -> rhop(',IDP,') = ',RHOP(J,IDP),' [kg/m3]'
    !write(UNIT_INFO(I),1613)'-- ->  alp(',IDP,') = ',PPI/6./LXMAX/LYMAX/LZMAX*DPART(J,IDP)**3*NPCLASS(J,IDP),' [-]'
    !write(UNIT_INFO(I),1613)'-- -> taup(',IDP,') = ',RHOP(J,IDP)*DPART(J,IDP)**2/18./VISC/RHOF,' [s]'
    write(UNIT_INFO(I),*)
   end do
  end if !- end if polydisp
  write(UNIT_INFO(I),*)
 end if !- end ifPARTDEF
end do
end if
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)



if(LEVEL_STFLU>0.or.LEVEL0_STPAR) then
write(UNIT_INFO(I),*)'====================================================================='
write(UNIT_INFO(I),*)'STATISTICS'
write(UNIT_INFO(I),*)'====================================================================='
write(UNIT_INFO(I),*)

if(STAT_TIME) then
 if(READSTAT) then
 write(UNIT_INFO(I),*)' -->>> Time-averaged statistics from file: STAT.ini'
 else
 write(UNIT_INFO(I),*)' -->>> Time-averaged statistics'
 end if
 write(UNIT_INFO(I),*)
end if
end if
if(LEVEL_STFLU>0) then
 write(UNIT_INFO(I),*)' Statistics over the fluid'
 write(UNIT_INFO(I),*)' -------------------------'
if(READSTAT) then
 write(UNIT_INFO(I),1611) ' Number of stored event: NEVEN=',NEVEN
end if
 write(UNIT_INFO(I),*)
if(LEVEL_STFLU>=1) then
 write(UNIT_INFO(I),*)' -- Mean and fluctuating motion'
 write(UNIT_INFO(I),*)'      + Mean fluid velocity'
 write(UNIT_INFO(I),*)'      + Fluctuating motion'
 write(UNIT_INFO(I),*)'      + Reynolds stress tensor'
 write(UNIT_INFO(I),*)
 write(UNIT_INFO(I),*)' -- Small scales statistics'
 write(UNIT_INFO(I),*)'      + Dissipation'
 write(UNIT_INFO(I),*)'      + Kolmogorov scales'
 write(UNIT_INFO(I),*)'      + Energy budget'
 write(UNIT_INFO(I),*)
end if
if(LEVEL_STFLU>=2) then
 write(UNIT_INFO(I),*)' -- Divergence'
 write(UNIT_INFO(I),*)
end if
if(LEVEL_STFLU>=3) then
 write(UNIT_INFO(I),*)' -- PDF and gradient based statistics'
 write(UNIT_INFO(I),*)'      + Skewness and flatness coefficients'
 write(UNIT_INFO(I),*)'      + PDF of fluid velocity'
 write(UNIT_INFO(I),*)'      + PDF of fluid velocity gradient'
end if
if(LEVEL_STFLU>=4) then
 write(UNIT_INFO(I),*)' -- Spatial correlation'
 write(UNIT_INFO(I),*)'      + No longer available !!!!!!!!'
 write(UNIT_INFO(I),*)
end if
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)
end if


if(LEVEL0_STSCL) then
 write(UNIT_INFO(I),*)' Statistics on scalar'
 write(UNIT_INFO(I),*)' --------------------'
if(LEVEL1_STPAR) then
 write(UNIT_INFO(I),*)' -- Mean and fluctuating'
 write(UNIT_INFO(I),*)'  + Mean and variance'
 write(UNIT_INFO(I),*)'  + Scalalr-velocity cross correlation'
 write(UNIT_INFO(I),*)
end if
if(LEVEL2_STPAR) then
 write(UNIT_INFO(I),*)' -- Dissipation and gradients'
 write(UNIT_INFO(I),*)
end if

end if


if(LEVEL0_STPAR) then
 write(UNIT_INFO(I),*)' Statistics over the particle'
 write(UNIT_INFO(I),*)' -------------------------'
if(READSTAT) then
 write(UNIT_INFO(I),1611) ' Number of stored event: NEVEN=',NEVEN
end if
 write(UNIT_INFO(I),*)
if(LEVEL1_STPAR) then
 write(UNIT_INFO(I),*)' -- Mean and fluctuating motion'
 write(UNIT_INFO(I),*)'  + Mean velocity'
 write(UNIT_INFO(I),*)'  + Fluctuating motion'
 write(UNIT_INFO(I),*)'  + Particle kinetic stress tensor'
 write(UNIT_INFO(I),*)'  + Mean drag coeffecient and relaxation time'
 write(UNIT_INFO(I),*)
end if
if(LEVEL2_STPAR) then
 write(UNIT_INFO(I),*)' -- Lagrangian correlation function'
 write(UNIT_INFO(I),1608)'  + maximum record of correlation function=',DIMLGR
 write(UNIT_INFO(I),*)
end if
if(LEVEL3_STPAR) then
 write(UNIT_INFO(I),*)' -- Radial Distribution Function'
 write(UNIT_INFO(I),*)
end if

write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)
end if

!!====================================================================
1600 format (1x,A,I4,A,I4,A,I4)
1601 format (1x,A,f10.4,A,f10.4,A,f10.4)
1602 format (1x,A,e10.3)
1603 format (1x,A,I5)
1604 format (1x,A,f5.2,A)
1605 format (1x,A,e9.2)
1606 format (1x,2x,A,I2,A)
1607 format (1x,A,I2)
1608 format (1x,A,I5)
1609 format (1x,A,I7)
1610 format (1x,A,I3,3x,A,I4,3x,A,E13.6)
1611 format (1x,A,f8.2)
1612 format (1x,A,e10.3,A)
1613 format (1x,A,I2.2,A,e10.3,A)
1614 format (1x,A,I2.2,A,I8,A)
1615 format (1x,A,I2.2,A,F8.2,A)

end subroutine INFORUN
