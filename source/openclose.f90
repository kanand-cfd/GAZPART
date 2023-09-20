!!======================================================================
!!
!!
!!======================================================================
subroutine OPENCLOSE(FLAG)

!!======================================================================
!!
!! Units between 300 and 399 are reserved for the fluid
!!               400 and XXX                      particles 
!!
!!======================================================================
use DNS_DIM            !- Dimension
use PARAM_PHYS         !- Physical & numerical parameters

implicit none


!-----------------------------------------------------------------------
! ARRAYS STATEMENT
!-----------------------------------------------------------------------
logical, intent(in) :: FLAG

character(len=40) :: FILENAME
integer :: I, J, IJ, IDP1, IDP2

!- file number
integer :: NUMFILE

!- number of opened files
integer :: NBFILE
!-----------------------------------------------------------------------

NBFILE = 0

if(FLAG) then

!!======================================================================
!! 1. Open files for fluid
!!======================================================================
 if(MYID==0) then
 
  if(LEVEL_STFLU>=1) then
   NBFILE = NBFILE + 1
   FILENAME = 'fluid_l1.stat'
   open(unit=301, file=trim(FILENAME), status='replace')
   write(301,20001)
   if(MYID==0) write(*,10700)trim(FILENAME)
  end if
  if(LEVEL_STFLU>=2) then
   NBFILE = NBFILE + 1
   FILENAME = 'fluid_l2.stat'
   open(unit=302, file=trim(FILENAME), status='replace')
   write(302,20002)
   if(MYID==0) write(*,10700)trim(FILENAME)
  end if
  if(LEVEL_STFLU>=3) then
   NBFILE = NBFILE + 1
   FILENAME = 'fluid_l3.stat'
   open(unit=303, file=trim(FILENAME), status='replace')
   write(303,20003)
   if(MYID==0) write(*,10700)trim(FILENAME)
  end if
  if(LEVEL_STFLU>=4) then
   NBFILE = NBFILE + 1
   FILENAME = 'fluid_l4.stat'
   open(unit=304, file=trim(FILENAME), status='replace')
   write(304,20004)
   if(MYID==0) write(*,10700)trim(FILENAME)
  end if


!!======================================================================
!! 2. Open files for scalar
!!======================================================================
if(LEVEL0_STSCL) then 
  if(LEVEL1_STSCL) then
   NBFILE = NBFILE + 1
   FILENAME = 'scalar_l1.stat'
   open(unit=401, file=trim(FILENAME), status='replace')
   write(401,20201)
   if(MYID==0) write(*,10700)trim(FILENAME)
  end if

  if(LEVEL2_STSCL) then
   NBFILE = NBFILE + 1
   FILENAME = 'scalar_l2.stat'
   open(unit=402, file=trim(FILENAME), status='replace')
   write(402,20202)
   if(MYID==0) write(*,10700)trim(FILENAME)
  end if
end if




!!======================================================================
!! 3. Open files for particles
!!======================================================================
if(LEVEL0_STPAR) then 

 if(LEVEL1_STPAR) then
  do I =1, NIG

   do IDP1 =1, POLYDISP(I)
    !--------------------------------------------------------
    NBFILE = NBFILE + 1
    NUMFILE = 500+(I-1)*POLYDISP(I)+IDP1  
    write(FILENAME,10601)'part_l1_trans_p',I,'_c',IDP1,'.stat'
    open(unit=NUMFILE, file=trim(FILENAME), status='replace')
    write(NUMFILE,20101)
    if(MYID==0) write(*,10700)trim(FILENAME)

    !--------------------------------------------------------
    NBFILE = NBFILE + 1
    NUMFILE = 520+(I-1)*POLYDISP(I)+IDP1  
    write(FILENAME,10601)'part_l1_rot_p',I,'_c',IDP1,'.stat'
    open(unit=NUMFILE, file=trim(FILENAME), status='replace')
    write(NUMFILE,20102)
    if(MYID==0) write(*,10700)trim(FILENAME)


    if(LEVEL1_STPDF) then

      NBFILE = NBFILE + 1
      NUMFILE = 600+(I-1)*POLYDISP(I)+IDP1  
      write(FILENAME,10601)'part_l1_pdf_p',I,'_c',IDP1,'.stat'
      open(unit=NUMFILE, file=trim(FILENAME), status='replace')
      write(NUMFILE,20103)
      if(MYID==0) write(*,10700)trim(FILENAME)

    end if

    if(FROZEN_FLOW>0) then

      NBFILE = NBFILE + 1
      NUMFILE = 540+(I-1)*POLYDISP(I)+IDP1  
      write(FILENAME,10601)'part_l1_fluid_p',I,'_c',IDP1,'.stat'
      open(unit=NUMFILE, file=trim(FILENAME), status='replace')
      write(NUMFILE,20108)
      if(MYID==0) write(*,10700)trim(FILENAME)

    end if


    if(SOLVE_COLLISION) then

      if(LEVEL1_STPDF) then

        NBFILE = NBFILE + 1
        NUMFILE = 650 + (I-1)*POLYDISP(I)+IDP1
        write(FILENAME,10601)'part_coll_pdf_p',I,'_c',IDP1,'.stat'
        open(unit=NUMFILE, file=trim(FILENAME), status='replace')
        write(NUMFILE, 20205)
        if(MYID==0) write(*, 10700)trim(FILENAME)

      end if

     !--------------------------------------------------------
     do IDP2 = 1, POLYDISP(I)
      NBFILE = NBFILE + 1
      NUMFILE = 560 
      NUMFILE = NUMFILE + sum(POLYDISP(1:I)**2)-POLYDISP(I)**2
      NUMFILE = NUMFILE + (IDP1-1)*POLYDISP(I)+IDP2

      write(FILENAME,10602)'part_l1_col_p',I,'_c',IDP1,IDP2,'.stat'
      open(unit=NUMFILE, file=trim(FILENAME), status='replace')
      write(NUMFILE,20107)
      if(MYID==0) write(*,10700)trim(FILENAME)

     end do !- Loop on IDP2
  
    end if

  end do !- Loop on IDP1


  end do

 end if

if(LEVELX_STPAR) then
  
  do J = 1, NIG
    
    do IDP1 = 1, POLYDISP(J)
      
      NBFILE = NBFILE + 1
      NUMFILE = 8000 + J*100 + IDP1*10 
      write(FILENAME,10601)'part_chan_p',J,'_',IDP1,'.stat'
      open(unit=NUMFILE, file=trim(FILENAME), status='replace')
      write(NUMFILE,20203)
      if(MYID==0) write(*,10700)trim(FILENAME)

      do IDP2 = 1, POLYDISP(J)

        NBFILE = NBFILE + 1
        NUMFILE = 5000 + J*200 + 10*IDP1 + IDP2
        write(FILENAME, 10602)'part_chan_coll_',J,'_c',IDP1,IDP2,'.stat'
        open(unit=NUMFILE, file=trim(FILENAME), status='replace')
        write(NUMFILE, 20204)
        if(MYID==0) write(*,10700)trim(FILENAME)

      end do

    end do
  
  end do
  
end if

!if(LEVEL2_STPAR) then
!  do I =1, NIG
!
!    NBFILE = NBFILE + 1
!    write(FILENAME,10600)'part_l2_Rp_p',I,'.stat'
!    open(unit=600+I, file=trim(FILENAME), status='replace')
!    
!    if(SOLVE_SCALAR) then 
!      write(600+I,20106)
!    else
!     write(600+I,20104)
!    end if 
!
!    if(MYID==0) write(*,10700)trim(FILENAME)
!  end do
!end if


end if

end if


if(LEVEL1_STPAR .and. STAT_TIME) then

  do I = 1, NIG
    do IDP1 = 1, POLYDISP(I)

      NBFILE = NBFILE + 1
      NUMFILE = 800 + 10*I + IDP1
      write(FILENAME, 10601)'part_time_avg',I,'_c',IDP1,'.stat'
      open(unit=NUMFILE, file=trim(FILENAME), status='replace')
      if(MYID==0) write(*,10700)trim(FILENAME)


      if(FROZEN_FLOW > 0) then
      NBFILE = NBFILE + 1
      NUMFILE = 810+ 10*I + IDP1
      write(FILENAME,10601)'partfluid_time_avg',I,'_c',IDP1,'.stat'
      open(unit=NUMFILE, file=trim(FILENAME), status='replace')
      if(MYID==0) write(*,10700)trim(FILENAME)
      end if

    end do
  end do

end if


!if(LEVEL1_STPAR .and. SOLVE_COLLISION) then
!  do I = 1, NIG
!    do IDP1 = 1, POLYDISP(I)
!      do IDP2 = 1, POLYDISP(I)
!        NBFILE = NBFILE + 1
!        NUMFILE = 7000 + 100*I + 10*IDP1 + IDP2
!        write(FILENAME, 10602)'part_l1_pdfcol_p',I,'c',IDP1,IDP2,'.stat'
!        open(unit=NUMFILE, file=trim(FILENAME), status='replace')
!        if(MYID==0) write(*,10700)trim(FILENAME)
!      end do
!    end do
!  end do
!end if




if(MYID==0) write(*,10800) 'Files opened --> ',NBFILE

!!====================================================================
!! 4. CLOSE FILES
!!====================================================================
else

if(MYID==0) then

 close(301)
 close(302)
 close(303)
 close(304)


end if


end if


!-----------------------------------------------------------------------
10600 format (A,I2.2,A)
10601 format (A,I2.2,A,I1,A)
10602 format (A,I2.2,A,I1,I1,A)

10700 format (1x,'   +   ',A)
10800 format (1x,A,I3)

20001 format('# t, <uf>, <vf>, <wf>, <kf>, <epsf>, <uf.uf>,<vf.vf>,<wf.wf>,<uf.vf>,<uf.wf>,<vf.wf>,eta_K,tau_K,v_K, Taylor g')
20002 format('# t, max(|dui/dxi|), < dui/dxi >, <(dui/dxi)^2>')
20003 format('# r, Rii, Rij, f, g')
20004 format('# t, Sk, Tk, Taylor f, Taylor g, Tay_f/Tay_g/2^0.5')


20101 format('# t, <up>, <vp>, <wp>, <up.up>, <vp.vp>, <wp.wp>, <up.vp>, <up.wp>, <vp.wp>')
20102 format('# t, <ox>, <oy>, <oz>, g<ox>, g<oy>, g<oz>, <ox.ox>, <oy.oy>, <oz.oz>, <ox.oy>, <ox.oz>, <oy.oz>, Ixx.<ox.ox>, Iyy.<oy.oy>, Izz.<oz.oz>, g<ox.ox>, g<oy.oy>, g<oz.oz>, g<ox.oy>, g<ox.oz>, g<oy.oz>, g<Ixx.ox.ox>, g<Iyy.oy.oy>, g<Izz.oz.oz>, g<Ixy.ox.oy>, g<Ixz.ox.oz>, g<Iyz.oy.oz>, Total<qp>, <qp>Trans, <qp>Rot, Np')
20103 format('# idpf, up, vp, wp, oxp, oyp, ozp, global-oxp, global-oyp, global-ozp, px, py, pz')

!20101 format('# t, <up>, <vp>, <wp>, <uf@p>, <vf@p>, <wf@p>, Rp,xx, Rp,yy, Rp,zz, Rp,xy, Rp,xz, Rp,yz,Rf@p,xx, Rf@p,yy, Rf@p,zz, Rf@p,xy, Rf@p,xz, Rf@p,yz,Rfp,xx, Rfp,yy, Rfp,zz, Rfp,xy, Rfp,xz, Rfp,yz')
!20103 format('# t, qp, qfp, qf@p, <1/tp>, <Vr>, <Cd>, <Rep>, <Np>, <dp>, <dp^2>, <dp^3>')
20104 format('# t, Rpx, Rpy, Rpz, Rp, Rf@px, Rf@py, Rf@pz, Rf@p')
20105 format('# t, <Tp>, <Tf@p>, <Tp^2>, <Tf@p^2>, <Tp.up>, <Tp.vp>, <Tp.wp>,<Tp.uf@p>,<Tp.vf@p>, <Tp.wf@p>, <Tf@p.up>, <Tf@p.vp>, <Tf@p.wp>, <Tf@p.uf@p>,<Tf@p.vf@p>, <Tf@p.wf@p>,<1/tp>, <Nup>, <Tp.Tf@p>')

20106 format('# t, Rpx, Rpy, Rpz, Rp, Rf@px, Rf@py, Rf@pz, Rf@p, RTp, RTf@p, RTfvf, RvfTf')
20107 format('# t, Nclose, Nover, Ncol, 1/tcol, dqp, <overlap>/dp')

20108 format('# t, Rep, CDRAG, CLIFT, FDRAGx, FDRAGy, FDRAGz, FLIFTx, FLIFTy, FLIFTz, !FDRAG!, !FLIFT!, Re_R, CTPITCH, CTROTATION, T_P_x, T_P_y, T_P_z, T_R_x, T_R_y, T_R_z, !T_P!, !T_R!, TAUPF, Gravity, phi, px, py, pz, INVTAUP, <Np>')

20201 format('# t, <theta>, <theta^2>, <uf.theta>, <vf.theta>, <wf.theta>')
20202 format('# t, <theta^2>, diss, prod, Sx, Sy, Sz, Kx, Ky, Kz')

20203 format('# x, delx, <up>,<vp>,<wp>,<uf@p>,<vf@p>,<wf@p>,<up.up>,<vp.vp>,<wp.wp>,<up.vp>,<up.wp>,<vp.wp>,<uf@p.uf@p>,<vf@p.vf@p>,<wf@p.wf@p>,<omega_px>,<omega_py>,<omega_pz>,<omega_px.omega_px>,<omega_py.omega_py>,<omega_pz.omega_pz>,<omega_px.omega_py>,<omega_px.omega_pz>,<omega_py.omega_pz>,FD-X,FD-Y,FD-Z,FL-X,FL-Y,FL-Z,PTq-X,PTq-Y,PTq-Z,RTq-X,RTq-Y,RTq-Z,Rep,phi,Cd,Cl,CTPITCH,CTROTATION,CTROTATION_2,Np')

20204 format('# x, <up>, <vp>, <wp>, <up.up>, <vp.vp>, <wp.wp>, <up.vp>, <up.wp>, <vp.wp>')
20205 format('# i, Wk, !Wk!, Theta, px1, py1, pz1, px2, py2, pz2, px12, py12, pz12')

end subroutine OPENCLOSE
