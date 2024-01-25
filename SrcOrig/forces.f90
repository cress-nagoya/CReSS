!***********************************************************************
      module m_forces
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/02/01
!     Modification: 2004/03/05, 2004/04/15, 2004/05/07, 2004/05/31,
!                   2004/06/10, 2006/01/10, 2006/02/13, 2006/04/03,
!                   2006/05/12, 2006/06/21, 2006/07/21, 2006/11/06,
!                   2007/01/20, 2007/05/07, 2007/07/30, 2008/05/02,
!                   2008/08/25, 2008/12/11, 2009/02/27, 2009/03/23,
!                   2010/05/17, 2010/09/22, 2011/05/16, 2011/08/18,
!                   2011/09/22, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures for the scalar forcing terms
!     calculation and bridge the complicated procedure for eddy
!     diffusivity.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_comkind
      use m_forcept
      use m_forceq
      use m_forceqa
      use m_forceqcg
      use m_forceqt
      use m_forceqv
      use m_getiname
      use m_kh8uv
      use m_ndgsat

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: forces, s_forces

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface forces

        module procedure s_forces

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_forces(fpcphopt,fpaslopt,fptrkopt,fptubopt,fmois,    &
     &                    ksp0,ctime,nggdmp,ngrdmp,gtinc,atinc,rtinc,   &
     &                    ni,nj,nk,nqw,nnw,nqi,nni,nqa,zph,j31,j32,     &
     &                    jcb,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,pbr,ptbr,  &
     &                    qvbr,rbr,rst,w,wp,ppp,ptp,ptpp,qv,qvp,        &
     &                    qwtr,qwtrp,nwtr,nwtrp,qice,qicep,nice,nicep,  &
     &                    qcwtr,qcwtrp,qcice,qcicep,qasl,qaslp,qt,qtp,  &
     &                    rbcxy,rbct,ptpgpv,ptptd,qvgpv,qvtd,qwgpv,qwtd,&
     &                    qigpv,qitd,qagpv,qatd,qwrdr,qwrtd,qirdr,qirtd,&
     &                    ptfrc,qvfrc,qwfrc,nwfrc,qifrc,nifrc,          &
     &                    qcwfrc,qcifrc,qafrc,qtfrc,rstxu,rstxv,rstxwc, &
     &                    rkh,rkv8w,rkh8u,rkh8v,tmp1,tmp2,tmp3,         &
     &                    tmp4,tmp5,tmp6)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fpaslopt
                       ! Formal parameter of unique index of aslopt

      integer, intent(in) :: fptrkopt
                       ! Formal parameter of unique index of trkopt

      integer, intent(in) :: fptubopt
                       ! Formal parameter of unique index of tubopt

      integer, intent(in) :: ksp0(1:2)
                       ! Index of lowest vertical sponge level

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nnw
                       ! Number of categories of water concentrations

      integer, intent(in) :: nqi
                       ! Number of categories of ice hydrometeor

      integer, intent(in) :: nni
                       ! Number of categories of ice concentrations

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

      real, intent(in) :: nggdmp
                       ! Analysis nudging damping coefficient for GPV

      real, intent(in) :: ngrdmp(1:2)
                       ! Analysis nudging damping coefficient for radar

      real, intent(in) :: gtinc
                       ! Lapse of forecast time
                       ! from GPV data reading

      real, intent(in) :: atinc
                       ! Lapse of forecast time
                       ! from aerosol data reading

      real, intent(in) :: rtinc(1:2)
                       ! Lapse of forecast time
                       ! from radar data reading

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(in) :: j31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of Jacobian

      real, intent(in) :: j32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of Jacobian

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: jcb8u(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at u points

      real, intent(in) :: jcb8v(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at v points

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(in) :: rmf8u(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at u points

      real, intent(in) :: rmf8v(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at v points

      real, intent(in) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: qvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state water vapor mixing ratio

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at present

      real, intent(in) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(in) :: ppp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at past

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at present

      real, intent(in) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at present

      real, intent(in) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(in) :: qwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at present

      real, intent(in) :: qwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at past

      real, intent(in) :: nwtr(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at present

      real, intent(in) :: nwtrp(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at past

      real, intent(in) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at present

      real, intent(in) :: qicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at past

      real, intent(in) :: nice(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at present

      real, intent(in) :: nicep(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at past

      real, intent(in) :: qcwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at present

      real, intent(in) :: qcwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at past

      real, intent(in) :: qcice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at present

      real, intent(in) :: qcicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at past

      real, intent(in) :: qasl(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at present

      real, intent(in) :: qaslp(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at past

      real, intent(in) :: qt(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at present

      real, intent(in) :: qtp(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at past

      real, intent(in) :: rbcxy(1:ni,1:nj)
                       ! Relaxed lateral sponge damping coefficients

      real, intent(in) :: rbct(1:ni,1:nj,1:nk,1:2)
                       ! Relaxed top sponge damping coefficients

      real, intent(in) :: ptpgpv(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation of GPV data
                       ! at marked time

      real, intent(in) :: ptptd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! potential temperature perturbation of GPV data

      real, intent(in) :: qvgpv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio of GPV data
                       ! at marked time

      real, intent(in) :: qvtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! water vapor mixing ratio of GPV data

      real, intent(in) :: qwgpv(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor of GPV data at marked time

      real, intent(in) :: qwtd(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Time tendency of water hydrometeor of GPV data

      real, intent(in) :: qigpv(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor of GPV data at marked time

      real, intent(in) :: qitd(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Time tendency of ice hydrometeor of GPV data

      real, intent(in) :: qagpv(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Mixing ratio of aerosol data at marked time

      real, intent(in) :: qatd(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Time tendency of mixing ratio of aerosol data

      real, intent(in) :: qwrdr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor of radar data
                       ! at marked time

      real, intent(in) :: qwrtd(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Time tendency of
                       ! water hydrometeor of radar data

      real, intent(in) :: qirdr(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor of radar data
                       ! at marked time

      real, intent(in) :: qirtd(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Time tendency of
                       ! ice hydrometeor of radar data

      real, intent(in) :: rstxu(0:ni+1,0:nj+1,1:nk)
                       ! u x base state density x Jacobian

      real, intent(in) :: rstxv(0:ni+1,0:nj+1,1:nk)
                       ! v x base state density x Jacobian

      real, intent(in) :: rstxwc(0:ni+1,0:nj+1,1:nk)
                       ! wc x base state density x Jacobian

! Input and output variables

      real, intent(inout) :: rkh(0:ni+1,0:nj+1,1:nk)
                       ! rbr x horizontal eddy diffusivity / Jacobian

      real, intent(inout) :: rkv8w(0:ni+1,0:nj+1,1:nk)
                       ! rbr x vertical eddy diffusivity / Jacobian
                       ! at w points

      real, intent(inout) :: ptfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in potential temperature equation

      real, intent(inout) :: qvfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in water vapor mixing ratio equation

! Output variables

      real, intent(out) :: qwfrc(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Forcing terms in water hydrometeor equations

      real, intent(out) :: nwfrc(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Forcing terms in water concentrations equations

      real, intent(out) :: qifrc(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Forcing terms in ice hydrometeor equations

      real, intent(out) :: nifrc(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Forcing terms in ice concentrations equations

      real, intent(out) :: qcwfrc(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Forcing terms
                       ! in charging distribution for water equations

      real, intent(out) :: qcifrc(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Forcing terms
                       ! in charging distribution for ice equations

      real, intent(out) :: qafrc(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Forcing terms in aerosol equations

      real, intent(out) :: qtfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in tracer equation

! Internal shared variables

      integer cphopt   ! Option for cloud micro physics
      integer aslopt   ! Option for aerosol processes
      integer trkopt   ! Option for mixing ratio tracking
      integer tubopt   ! Option for turbulent mixing

      real, intent(inout) :: rkh8u(0:ni+1,0:nj+1,1:nk)
                       ! 2.0 x rbr x horizontal eddy diffusivity / jcb
                       ! at u points

      real, intent(inout) :: rkh8v(0:ni+1,0:nj+1,1:nk)
                       ! 2.0 x rbr x horizontal eddy diffusivity / jcb
                       ! at v points

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp5(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp6(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Remark

!     rkh: This variable is also temporary, because it is not used
!          again.

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpcphopt,cphopt)
      call getiname(fpaslopt,aslopt)
      call getiname(fptrkopt,trkopt)
      call getiname(fptubopt,tubopt)

! -----

! Set the horizontal eddy diffusivity at the u and v points.

      if(tubopt.ge.1) then

        call kh8uv(idmpopt,idmfcopt,ni,nj,nk,rmf,rkh,rkh8u,rkh8v)

      end if

! -----

! Calculate the forcing term in the potential temperature equation.

      call forcept(idnggvar,idlspvar,idvspvar,idlspopt,idvspopt,        &
     &            idgwmopt,idsmtopt,idadvopt,idtubopt,ksp0,nggdmp,gtinc,&
     &            ni,nj,nk,j31,j32,jcb,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,  &
     &            ptbr,rbr,rst,rstxu,rstxv,rstxwc,w,wp,ptp,ptpp,        &
     &            rkh8u,rkh8v,rkv8w,rbcxy,rbct,ptpgpv,ptptd,ptfrc,      &
     &            tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,rkh)

! -----

!! Calculate the forcing terms for the hydrometeor.

      if(fmois(1:5).eq.'moist') then

! Calculate the forcing term in the water vapor mixing ratio equation.

        call forceqv(idnggvar,idlspvar,idvspvar,idlspopt,idvspopt,      &
     &              idsmtopt,idtubopt,ksp0,nggdmp,gtinc,ni,nj,nk,       &
     &              j31,j32,jcb,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,qvbr,rbr,&
     &              rst,rstxu,rstxv,rstxwc,qv,qvp,rkh8u,rkh8v,rkv8w,    &
     &              rbcxy,rbct,qvgpv,qvtd,qvfrc,                        &
     &              tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)

! -----

! Calculate the forcing terms in the water and ice hydrometeor
! equations.

        call forceq(idnggvar,idlspvar,idvspvar,idngrvar,                &
     &              idlspopt,idvspopt,idngropt,idsmtopt,idcphopt,       &
     &              idhaiopt,idtubopt,ksp0,nggdmp,ngrdmp,gtinc,rtinc,   &
     &              ni,nj,nk,nqw,nnw,nqi,nni,j31,j32,jcb,jcb8u,jcb8v,   &
     &              mf,rmf,rmf8u,rmf8v,rbr,rst,rstxu,rstxv,rstxwc,      &
     &              qwtr,qwtrp,nwtr,nwtrp,qice,qicep,nice,nicep,        &
     &              rkh8u,rkh8v,rkv8w,rbcxy,rbct,qwgpv,qwtd,qigpv,qitd, &
     &              qwrdr,qwrtd,qirdr,qirtd,qwfrc,nwfrc,qifrc,nifrc,    &
     &              tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)

! -----

! Calculate the forcing terms in the charging distribution equation.

        call forceqcg(idlspvar,idvspvar,idlspopt,idvspopt,              &
     &                idsmtopt,idcphopt,idhaiopt,idqcgopt,idtubopt,     &
     &                ksp0,ni,nj,nk,nqw,nqi,j31,j32,jcb,jcb8u,jcb8v,    &
     &                mf,rmf,rmf8u,rmf8v,rbr,rst,rstxu,rstxv,rstxwc,    &
     &                qcwtr,qcwtrp,qcice,qcicep,rkh8u,rkh8v,rkv8w,      &
     &                rbcxy,rbct,qcwfrc,qcifrc,tmp1,tmp2,tmp3,          &
     &                tmp4,tmp5,tmp6)

! -----

      end if

!! -----

! Calculate the forcing terms in the aerosol equations.

      if(aslopt.ge.1) then

        call forceqa(idnggvar,idlspvar,idvspvar,idlspopt,idvspopt,      &
     &               idsmtopt,idtubopt,ksp0,nggdmp,atinc,ni,nj,nk,nqa,  &
     &               j31,j32,jcb,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,rbr,    &
     &               rst,rstxu,rstxv,rstxwc,qasl,qaslp,rkh8u,rkh8v,     &
     &               rkv8w,rbcxy,rbct,qagpv,qatd,qafrc,                 &
     &               tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)

      end if

! -----

! Calculate the forcing term in the tracer equation.

      if(trkopt.ge.1) then

        call forceqt(idlspvar,idvspvar,idlspopt,idvspopt,idsmtopt,      &
     &               idtrkopt,idtubopt,idqtstr,idqtend,ksp0,ctime,      &
     &               ni,nj,nk,zph,j31,j32,jcb,jcb8u,jcb8v,mf,rmf,       &
     &               rmf8u,rmf8v,rbr,rst,rstxu,rstxv,rstxwc,qt,qtp,     &
     &               rkh8u,rkh8v,rkv8w,rbcxy,rbct,qtfrc,                &
     &               tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)

      end if

! -----

! Perform analysis nudging to radar data of water vapor mixing ratio.

      if(fmois(1:5).eq.'moist') then

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.2) then

            call s_ndgsat(idngrvar,idngropt,ngrdmp,rtinc,               &
     &                    ni,nj,nk,nqw,nqi,zph,pbr,ptbr,rst,ppp,ptpp,   &
     &                    qvp,qwrdr,qwrtd,qirdr,qirtd,qvfrc,            &
     &                    tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)

          end if

        end if

      end if

! -----

      end subroutine s_forces

!-----7--------------------------------------------------------------7--

      end module m_forces
