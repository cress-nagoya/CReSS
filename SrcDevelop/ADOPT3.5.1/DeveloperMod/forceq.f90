!***********************************************************************
      module m_forceq
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/11/01
!     Modification: 1999/11/19, 1999/11/24, 2000/01/17, 2000/04/18,
!                   2000/06/01, 2000/12/19, 2001/03/13, 2001/04/15,
!                   2001/05/29, 2001/06/06, 2001/07/13, 2001/08/07,
!                   2002/04/02, 2002/06/18, 2002/07/23, 2002/08/15,
!                   2002/09/09, 2002/10/31, 2002/12/11, 2003/01/04,
!                   2003/03/13, 2003/03/21, 2003/04/30, 2003/05/19,
!                   2003/07/15, 2003/10/10, 2003/11/28, 2003/12/12,
!                   2004/02/01, 2004/03/05, 2004/04/15, 2004/05/31,
!                   2004/06/10, 2004/08/01, 2004/08/20, 2004/09/10,
!                   2005/10/05, 2006/01/10, 2006/02/13, 2006/04/03,
!                   2006/05/12, 2006/06/21, 2006/09/21, 2006/11/06,
!                   2007/05/07, 2007/07/30, 2007/11/26, 2008/05/02,
!                   2008/08/25, 2008/12/11, 2009/01/30, 2009/02/27,
!                   2009/03/23, 2010/05/17, 2011/08/18, 2011/09/22,
!                   2013/02/13, 2013/03/27

!     Author      : Satoki Tsujino
!     Modification: 2024/12/26

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the forcing terms in the water hydrometeor and the ice
!     hydrometeor and concentrations equations.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_advs
      use m_comindx
      use m_getcname
      use m_getiname
      use m_inichar
      use m_lsps
      use m_lsps0
      use m_nlsms
      use m_qp2rdr
      use m_s2gpv
      use m_smoo2s
      use m_smoo4s
      use m_turbflx
      use m_turbs
      use m_vsps
      use m_vsps0

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: forceq, s_forceq

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface forceq

        module procedure s_forceq

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_forceq(fpnggvar,fplspvar,fpvspvar,fpngrvar,          &
     &                    fplspopt,fpvspopt,fpngropt,fpsmtopt,          &
     &                    fpcphopt,fphaiopt,fptubopt,ksp0,nggdmp,ngrdmp,&
     &                    gtinc,rtinc,ni,nj,nk,nqw,nnw,nqi,nni,j31,j32, &
     &                    jcb,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,rbr,rst,   &
     &                    rstxu,rstxv,rstxwc,qwtr,qwtrp,nwtr,nwtrp,     &
     &                    qice,qicep,nice,nicep,rkh8u,rkh8v,rkv8w,      &
     &                    rbcxy,rbct,qwgpv,qwtd,qigpv,qitd,qwrdr,qwrtd, &
     &                    qirdr,qirtd,qwfrc,nwfrc,qifrc,nifrc,h3,       &
     &                    tmp1,tmp2,tmp3,tmp4,tmp5)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpnggvar
                       ! Formal parameter of unique index of nggvar

      integer, intent(in) :: fplspvar
                       ! Formal parameter of unique index of lspvar

      integer, intent(in) :: fpvspvar
                       ! Formal parameter of unique index of vspvar

      integer, intent(in) :: fpngrvar
                       ! Formal parameter of unique index of ngrvar

      integer, intent(in) :: fplspopt
                       ! Formal parameter of unique index of lspopt

      integer, intent(in) :: fpvspopt
                       ! Formal parameter of unique index of vspopt

      integer, intent(in) :: fpngropt
                       ! Formal parameter of unique index of ngropt

      integer, intent(in) :: fpsmtopt
                       ! Formal parameter of unique index of smtopt

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

      integer, intent(in) :: fptubopt
                       ! Formal parameter of unique index of tubopt

      integer, intent(in) :: ksp0(1:2)
                       ! Index of lowest vertical sponge level

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

      real, intent(in) :: nggdmp
                       ! Analysis nudging damping coefficient for GPV

      real, intent(in) :: ngrdmp(1:2)
                       ! Analysis nudging damping coefficient for radar

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

      real, intent(in) :: rtinc(1:2)
                       ! Lapse of forecast time from radar data reading

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

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: rstxu(0:ni+1,0:nj+1,1:nk)
                       ! u x base state density x Jacobian

      real, intent(in) :: rstxv(0:ni+1,0:nj+1,1:nk)
                       ! v x base state density x Jacobian

      real, intent(in) :: rstxwc(0:ni+1,0:nj+1,1:nk)
                       ! wc x base state density x Jacobian

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

      real, intent(in) :: rkh8u(0:ni+1,0:nj+1,1:nk)
                       ! 2.0 x rbr x horizontal eddy diffusivity / jcb
                       ! at u points

      real, intent(in) :: rkh8v(0:ni+1,0:nj+1,1:nk)
                       ! 2.0 x rbr x horizontal eddy diffusivity / jcb
                       ! at v points

      real, intent(in) :: rkv8w(0:ni+1,0:nj+1,1:nk)
                       ! rbr x vertical eddy diffusivity / jcb
                       ! at w points

      real, intent(in) :: rbcxy(1:ni,1:nj)
                       ! Relaxed lateral sponge damping coefficients

      real, intent(in) :: rbct(1:ni,1:nj,1:nk,1:2)
                       ! Relaxed top sponge damping coefficients

      real, intent(in) :: qwgpv(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor of GPV data at marked time

      real, intent(in) :: qwtd(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Time tendency of water hydrometeor of GPV data

      real, intent(in) :: qigpv(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor of GPV data at marked time

      real, intent(in) :: qitd(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Time tendency of ice hydrometeor of GPV data

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

! Output variables

      real, intent(out) :: qwfrc(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Forcing terms in water hydrometeor equations

      real, intent(out) :: nwfrc(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Forcing terms in water concentrations equations

      real, intent(out) :: qifrc(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Forcing terms in ice hydrometeor equations

      real, intent(out) :: nifrc(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Forcing terms in ice concentrations equations

! Internal shared variables

      character(len=108) nggvar
                       ! Control flag of
                       ! analysis nudged variables to GPV

      character(len=108) lspvar
                       ! Control flag of
                       ! lateral sponge damped variables

      character(len=108) vspvar
                       ! Control flag of
                       ! vertical sponge damped variables

      character(len=108) ngrvar
                       ! Control flag of
                       ! analysis nudged variables to radar data

      integer lspopt   ! Option for lateral sponge damping
      integer vspopt   ! Option for vertical sponge damping
      integer ngropt   ! Option for analysis nudging to radar

      integer smtopt   ! Option for numerical smoothing
      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes
      integer tubopt   ! Option for turbulent mixing

      integer n        ! Array index in bin categories

      real, intent(inout) :: h3(0:ni+1,0:nj+1,1:nk)
                       ! z components of turbulent fluxes

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

      character(len=108) dmpvar
                       ! Control flag of dump variables

!-----7--------------------------------------------------------------7--

! Initialize the character variables.

      call inichar(nggvar)
      call inichar(lspvar)
      call inichar(vspvar)
      call inichar(ngrvar)

! -----

! Get the required namelist variables.

      call getcname(fpnggvar,nggvar)
      call getcname(fplspvar,lspvar)
      call getcname(fpvspvar,vspvar)
      call getcname(fpngrvar,ngrvar)
      call getiname(fplspopt,lspopt)
      call getiname(fpvspopt,vspopt)
      call getiname(fpngropt,ngropt)
      call getiname(fpsmtopt,smtopt)
      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)
      call getiname(fptubopt,tubopt)

      call getcname(iddmpvar,dmpvar)

! -----

!!!!! Calculate the forcing terms in the water hydrometeor and the ice
!!!!! hydrometeor and concentrations equations for bulk categories.

      if(abs(cphopt).lt.10) then

!!! Calculate the forcing terms in the water hydrometeor equations.

        if(abs(cphopt).ge.1) then

!! For the cloud water.

! Calculate the advection.

          call s_advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,       &
     &                iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc, &
     &                qwtr(0,0,1,1),qwfrc(0,0,1,1),tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

          if(mod(smtopt,10).eq.1) then

            call s_smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,rbr,qwtrp(0,0,1,1),&
     &                    qwfrc(0,0,1,1),tmp1)

          end if

! -----

! Calculate the 4th order smoothing.

          if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

            call s_smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth,   &
     &                    idsmhcoe,idsmvcoe,ni,nj,nk,rbr,               &
     &                    qwtrp(0,0,1,1),qwfrc(0,0,1,1),                &
     &                    tmp1,tmp2,tmp3,tmp4,tmp5)

          end if

! -----

! Calculate the non linear smoothing.

          if(smtopt.ge.11) then

            call s_nlsms(idnlhcoe,idnlvcoe,.001e0,ni,nj,nk,             &
     &                   rbr,qwtrp(0,0,1,1),qwfrc(0,0,1,1),             &
     &                   tmp1,tmp2,tmp3,tmp4)

          end if

! -----

! Calculate the turbulent mixing.

          if(tubopt.ge.1) then

            call s_turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,      &
     &                     ni,nj,nk,j31,j32,jcb,qwtrp(0,0,1,1),         &
     &                     qwfrc(0,0,1,1),rkh8u,rkh8v,rkv8w,            &
     &                     tmp1,tmp2,h3,tmp3,tmp4,tmp5)

            call s_turbs(idtrnopt,idmpopt,idmfcopt,                     &
     &                   iddxiv,iddyiv,iddziv,iddmpvar,ni,nj,nk,'xx',   &
     &                   j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,        &
     &                   tmp1,tmp2,h3,qwfrc(0,0,1,1),tmp3,tmp4,tmp5)

          end if

! -----

! Perform the analysis nudging to GPV.

          if(nggdmp.gt.0.e0.and.nggvar(7:7).eq.'o') then

            call s_s2gpv(idgpvvar,3,nggdmp,gtinc,                       &
     &                   ni,nj,nk,rst,qwtrp(0,0,1,1),                   &
     &                   qwgpv(0,0,1,1),qwtd(0,0,1,1),qwfrc(0,0,1,1))

          end if

! -----

! Calculate the lateral sponge damping.

          if(lspopt.ge.1.and.lspvar(7:7).eq.'o') then

            call s_lsps(idgpvvar,idlspopt,idwdnews,idlsnews,idlspsmt,   &
     &                  3,gtinc,ni,nj,nk,rst,qwtrp(0,0,1,1),rbcxy,      &
     &                  qwgpv(0,0,1,1),qwtd(0,0,1,1),qwfrc(0,0,1,1),    &
     &                  tmp1)

          end if

! -----

! Calculate the vertical sponge damping.

          if(vspopt.ge.1.and.vspvar(7:7).eq.'o') then

            call s_vsps(idgpvvar,idvspopt,3,ksp0,                       &
     &                  gtinc,ni,nj,nk,rst,qwtrp(0,0,1,1),rbct,         &
     &                  qwgpv(0,0,1,1),qwtd(0,0,1,1),qwfrc(0,0,1,1))

          end if

! -----

!! -----

!! For the rain water.

! Calculate the advection.

          call s_advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,       &
     &                iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc, &
     &                qwtr(0,0,1,2),qwfrc(0,0,1,2),tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

          if(mod(smtopt,10).eq.1) then

            call s_smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,rbr,qwtrp(0,0,1,2),&
     &                    qwfrc(0,0,1,2),tmp1)

          end if

! -----

! Calculate the 4th order smoothing.

          if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

            call s_smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth,   &
     &                    idsmhcoe,idsmvcoe,ni,nj,nk,rbr,               &
     &                    qwtrp(0,0,1,2),qwfrc(0,0,1,2),                &
     &                    tmp1,tmp2,tmp3,tmp4,tmp5)

          end if

! -----

! Calculate the non linear smoothing.

          if(smtopt.ge.11) then

            call s_nlsms(idnlhcoe,idnlvcoe,.001e0,ni,nj,nk,             &
     &                   rbr,qwtrp(0,0,1,2),qwfrc(0,0,1,2),             &
     &                   tmp1,tmp2,tmp3,tmp4)

          end if

! -----

! Calculate the turbulent mixing.

          if(tubopt.ge.1) then

            call s_turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,      &
     &                     ni,nj,nk,j31,j32,jcb,qwtrp(0,0,1,2),         &
     &                     qwfrc(0,0,1,2),rkh8u,rkh8v,rkv8w,            &
     &                     tmp1,tmp2,h3,tmp3,tmp4,tmp5)

            call s_turbs(idtrnopt,idmpopt,idmfcopt,                     &
     &                   iddxiv,iddyiv,iddziv,iddmpvar,ni,nj,nk,'xx',   &
     &                   j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,        &
     &                   tmp1,tmp2,h3,qwfrc(0,0,1,2),tmp3,tmp4,tmp5)

          end if

! -----

! Perform the analysis nudging to GPV.

          if(nggdmp.gt.0.e0.and.nggvar(7:7).eq.'o') then

            call s_s2gpv(idgpvvar,4,nggdmp,gtinc,                       &
     &                   ni,nj,nk,rst,qwtrp(0,0,1,2),                   &
     &                   qwgpv(0,0,1,2),qwtd(0,0,1,2),qwfrc(0,0,1,2))

          end if

! -----

! Calculate the lateral sponge damping.

          if(lspopt.ge.1.and.lspvar(7:7).eq.'o') then

            call s_lsps(idgpvvar,idlspopt,idwdnews,idlsnews,idlspsmt,   &
     &                  4,gtinc,ni,nj,nk,rst,qwtrp(0,0,1,2),rbcxy,      &
     &                  qwgpv(0,0,1,2),qwtd(0,0,1,2),qwfrc(0,0,1,2),    &
     &                  tmp1)

          end if

! -----

! Calculate the vertical sponge damping.

          if(vspopt.ge.1.and.vspvar(7:7).eq.'o') then

            call s_vsps(idgpvvar,idvspopt,4,ksp0,                       &
     &                  gtinc,ni,nj,nk,rst,qwtrp(0,0,1,2),rbct,         &
     &                  qwgpv(0,0,1,2),qwtd(0,0,1,2),qwfrc(0,0,1,2))

          end if

! -----

! Perform the analysis nudging to radar data.

          if(ngropt.ge.1.and.ngrvar(5:5).eq.'o') then

            call s_qp2rdr(idngropt,ngrdmp,rtinc,                        &
     &                    ni,nj,nk,rst,qwtrp(0,0,1,2),                  &
     &                    qwrdr(0,0,1,2),qwrtd(0,0,1,2),qwfrc(0,0,1,2))

          end if

! -----

!! -----

        end if

!!! -----

!!! Calculate the forcing terms in the water concentrations equations.

        if(abs(cphopt).eq.4) then

!! For the cloud water concentrations.

! Calculate the advection.

          call s_advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,       &
     &                iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc, &
     &                nwtr(0,0,1,1),nwfrc(0,0,1,1),tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

          if(mod(smtopt,10).eq.1) then

            call s_smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,rbr,nwtrp(0,0,1,1),&
     &                    nwfrc(0,0,1,1),tmp1)

          end if

! -----

! Calculate the 4th order smoothing.

          if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

            call s_smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth,   &
     &                    idsmhcoe,idsmvcoe,ni,nj,nk,rbr,               &
     &                    nwtrp(0,0,1,1),nwfrc(0,0,1,1),                &
     &                    tmp1,tmp2,tmp3,tmp4,tmp5)

          end if

! -----

! Calculate the non linear smoothing.

          if(smtopt.ge.11) then

            call s_nlsms(idnlhcoe,idnlvcoe,1.e20,ni,nj,nk,              &
     &                   rbr,nwtrp(0,0,1,1),nwfrc(0,0,1,1),             &
     &                   tmp1,tmp2,tmp3,tmp4)

          end if

! -----

! Calculate the turbulent mixing.

          if(tubopt.ge.1) then

            call s_turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,      &
     &                     ni,nj,nk,j31,j32,jcb,nwtrp(0,0,1,1),         &
     &                     nwfrc(0,0,1,1),rkh8u,rkh8v,rkv8w,            &
     &                     tmp1,tmp2,h3,tmp3,tmp4,tmp5)

            call s_turbs(idtrnopt,idmpopt,idmfcopt,                     &
     &                   iddxiv,iddyiv,iddziv,iddmpvar,ni,nj,nk,'xx',   &
     &                   j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,        &
     &                   tmp1,tmp2,h3,nwfrc(0,0,1,1),tmp3,tmp4,tmp5)

          end if

! -----

! Calculate the lateral sponge damping.

          if(lspopt.ge.1.and.lspvar(7:7).eq.'o') then

            call s_lsps0(idlspopt,idwdnews,idlsnews,idlspsmt,           &
     &                   ni,nj,nk,rst,nwtrp(0,0,1,1),rbcxy,             &
     &                   nwfrc(0,0,1,1),tmp1)

          end if

! -----

! Calculate the vertical sponge damping.

          if(vspopt.ge.1.and.vspvar(7:7).eq.'o') then

            call s_vsps0(ksp0,ni,nj,nk,rst,nwtrp(0,0,1,1),rbct,         &
     &                   nwfrc(0,0,1,1))

          end if

! -----

!! -----

!! For the rain water concentrations.

! Calculate the advection.

          call s_advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,       &
     &                iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc, &
     &                nwtr(0,0,1,2),nwfrc(0,0,1,2),tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

          if(mod(smtopt,10).eq.1) then

            call s_smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,rbr,nwtrp(0,0,1,2),&
     &                    nwfrc(0,0,1,2),tmp1)

          end if

! -----

! Calculate the 4th order smoothing.

          if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

            call s_smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth,   &
     &                    idsmhcoe,idsmvcoe,ni,nj,nk,rbr,               &
     &                    nwtrp(0,0,1,2),nwfrc(0,0,1,2),                &
     &                    tmp1,tmp2,tmp3,tmp4,tmp5)

          end if

! -----

! Calculate the non linear smoothing.

          if(smtopt.ge.11) then

            call s_nlsms(idnlhcoe,idnlvcoe,1.e9,ni,nj,nk,               &
     &                   rbr,nwtrp(0,0,1,2),nwfrc(0,0,1,2),             &
     &                   tmp1,tmp2,tmp3,tmp4)

          end if

! -----

! Calculate the turbulent mixing.

          if(tubopt.ge.1) then

            call s_turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,      &
     &                     ni,nj,nk,j31,j32,jcb,nwtrp(0,0,1,2),         &
     &                     nwfrc(0,0,1,2),rkh8u,rkh8v,rkv8w,            &
     &                     tmp1,tmp2,h3,tmp3,tmp4,tmp5)

            call s_turbs(idtrnopt,idmpopt,idmfcopt,                     &
     &                   iddxiv,iddyiv,iddziv,iddmpvar,ni,nj,nk,'xx',   &
     &                   j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,        &
     &                   tmp1,tmp2,h3,nwfrc(0,0,1,2),tmp3,tmp4,tmp5)

          end if

! -----

! Calculate the lateral sponge damping.

          if(lspopt.ge.1.and.lspvar(7:7).eq.'o') then

            call s_lsps0(idlspopt,idwdnews,idlsnews,idlspsmt,           &
     &                   ni,nj,nk,rst,nwtrp(0,0,1,2),rbcxy,             &
     &                   nwfrc(0,0,1,2),tmp1)

          end if

! -----

! Calculate the vertical sponge damping.

          if(vspopt.ge.1.and.vspvar(7:7).eq.'o') then

            call s_vsps0(ksp0,ni,nj,nk,rst,nwtrp(0,0,1,2),rbct,         &
     &                   nwfrc(0,0,1,2))

          end if

! -----

!! -----

        end if

!!! -----

!!! Calculate the forcing terms in the ice hydrometeor equations.

        if(abs(cphopt).ge.2) then

!! For the cloud ice.

! Calculate the advection.

          call s_advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,       &
     &                iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc, &
     &                qice(0,0,1,1),qifrc(0,0,1,1),tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

          if(mod(smtopt,10).eq.1) then

            call s_smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,rbr,qicep(0,0,1,1),&
     &                    qifrc(0,0,1,1),tmp1)

          end if

! -----

! Calculate the 4th order smoothing.

          if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

            call s_smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth,   &
     &                    idsmhcoe,idsmvcoe,ni,nj,nk,rbr,               &
     &                    qicep(0,0,1,1),qifrc(0,0,1,1),                &
     &                    tmp1,tmp2,tmp3,tmp4,tmp5)

          end if

! -----

! Calculate the non linear smoothing.

          if(smtopt.ge.11) then

            call s_nlsms(idnlhcoe,idnlvcoe,.0001e0,ni,nj,nk,            &
     &                   rbr,qicep(0,0,1,1),qifrc(0,0,1,1),             &
     &                   tmp1,tmp2,tmp3,tmp4)

          end if

! -----

! Calculate the turbulent mixing.

          if(tubopt.ge.1) then

            call s_turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,      &
     &                     ni,nj,nk,j31,j32,jcb,qicep(0,0,1,1),         &
     &                     qifrc(0,0,1,1),rkh8u,rkh8v,rkv8w,            &
     &                     tmp1,tmp2,h3,tmp3,tmp4,tmp5)

            call s_turbs(idtrnopt,idmpopt,idmfcopt,                     &
     &                   iddxiv,iddyiv,iddziv,iddmpvar,ni,nj,nk,'xx',   &
     &                   j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,        &
     &                   tmp1,tmp2,h3,qifrc(0,0,1,1),tmp3,tmp4,tmp5)

          end if

! -----

! Perform the analysis nudging to GPV.

          if(nggdmp.gt.0.e0.and.nggvar(7:7).eq.'o') then

            call s_s2gpv(idgpvvar,5,nggdmp,gtinc,                       &
     &                   ni,nj,nk,rst,qicep(0,0,1,1),                   &
     &                   qigpv(0,0,1,1),qitd(0,0,1,1),qifrc(0,0,1,1))

          end if

! -----

! Calculate the lateral sponge damping.

          if(lspopt.ge.1.and.lspvar(7:7).eq.'o') then

            call s_lsps(idgpvvar,idlspopt,idwdnews,idlsnews,idlspsmt,   &
     &                  5,gtinc,ni,nj,nk,rst,qicep(0,0,1,1),rbcxy,      &
     &                  qigpv(0,0,1,1),qitd(0,0,1,1),qifrc(0,0,1,1),    &
     &                  tmp1)

          end if

! -----

! Calculate the vertical sponge damping.

          if(vspopt.ge.1.and.vspvar(7:7).eq.'o') then

            call s_vsps(idgpvvar,idvspopt,5,ksp0,                       &
     &                  gtinc,ni,nj,nk,rst,qicep(0,0,1,1),rbct,         &
     &                  qigpv(0,0,1,1),qitd(0,0,1,1),qifrc(0,0,1,1))

          end if

! -----

!! -----

!! For the snow.

! Calculate the advection.

          call s_advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,       &
     &                iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc, &
     &                qice(0,0,1,2),qifrc(0,0,1,2),tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

          if(mod(smtopt,10).eq.1) then

            call s_smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,rbr,qicep(0,0,1,2),&
     &                    qifrc(0,0,1,2),tmp1)

          end if

! -----

! Calculate the 4th order smoothing.

          if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

            call s_smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth,   &
     &                    idsmhcoe,idsmvcoe,ni,nj,nk,rbr,               &
     &                    qicep(0,0,1,2),qifrc(0,0,1,2),                &
     &                    tmp1,tmp2,tmp3,tmp4,tmp5)

          end if

! -----

! Calculate the non linear smoothing.

          if(smtopt.ge.11) then

            call s_nlsms(idnlhcoe,idnlvcoe,.001e0,ni,nj,nk,             &
     &                   rbr,qicep(0,0,1,2),qifrc(0,0,1,2),             &
     &                   tmp1,tmp2,tmp3,tmp4)

          end if

! -----

! Calculate the turbulent mixing.

          if(tubopt.ge.1) then

            call s_turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,      &
     &                     ni,nj,nk,j31,j32,jcb,qicep(0,0,1,2),         &
     &                     qifrc(0,0,1,2),rkh8u,rkh8v,rkv8w,            &
     &                     tmp1,tmp2,h3,tmp3,tmp4,tmp5)

            call s_turbs(idtrnopt,idmpopt,idmfcopt,                     &
     &                   iddxiv,iddyiv,iddziv,iddmpvar,ni,nj,nk,'xx',   &
     &                   j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,        &
     &                   tmp1,tmp2,h3,qifrc(0,0,1,2),tmp3,tmp4,tmp5)

          end if

! -----

! Perform the analysis nudging to GPV.

          if(nggdmp.gt.0.e0.and.nggvar(7:7).eq.'o') then

            call s_s2gpv(idgpvvar,6,nggdmp,gtinc,                       &
     &                   ni,nj,nk,rst,qicep(0,0,1,2),                   &
     &                   qigpv(0,0,1,2),qitd(0,0,1,2),qifrc(0,0,1,2))

          end if

! -----

! Calculate the lateral sponge damping.

          if(lspopt.ge.1.and.lspvar(7:7).eq.'o') then

            call s_lsps(idgpvvar,idlspopt,idwdnews,idlsnews,idlspsmt,   &
     &                  6,gtinc,ni,nj,nk,rst,qicep(0,0,1,2),rbcxy,      &
     &                  qigpv(0,0,1,2),qitd(0,0,1,2),qifrc(0,0,1,2),    &
     &                  tmp1)

          end if

! -----

! Calculate the vertical sponge damping.

          if(vspopt.ge.1.and.vspvar(7:7).eq.'o') then

            call s_vsps(idgpvvar,idvspopt,6,ksp0,                       &
     &                  gtinc,ni,nj,nk,rst,qicep(0,0,1,2),rbct,         &
     &                  qigpv(0,0,1,2),qitd(0,0,1,2),qifrc(0,0,1,2))

          end if

! -----

! Perform the analysis nudging to radar data.

          if(ngropt.ge.1.and.ngrvar(5:5).eq.'o') then

            call s_qp2rdr(idngropt,ngrdmp,rtinc,                        &
     &                    ni,nj,nk,rst,qicep(0,0,1,2),                  &
     &                    qirdr(0,0,1,2),qirtd(0,0,1,2),qifrc(0,0,1,2))

          end if

! -----

!! -----

!! For the graupel.

! Calculate the advection.

          call s_advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,       &
     &                iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc, &
     &                qice(0,0,1,3),qifrc(0,0,1,3),tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

          if(mod(smtopt,10).eq.1) then

            call s_smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,rbr,qicep(0,0,1,3),&
     &                    qifrc(0,0,1,3),tmp1)

          end if

! -----

! Calculate the 4th order smoothing.

          if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

            call s_smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth,   &
     &                    idsmhcoe,idsmvcoe,ni,nj,nk,rbr,               &
     &                    qicep(0,0,1,3),qifrc(0,0,1,3),                &
     &                    tmp1,tmp2,tmp3,tmp4,tmp5)

          end if

! -----

! Calculate the non linear smoothing.

          if(smtopt.ge.11) then

            call s_nlsms(idnlhcoe,idnlvcoe,.001e0,ni,nj,nk,             &
     &                   rbr,qicep(0,0,1,3),qifrc(0,0,1,3),             &
     &                   tmp1,tmp2,tmp3,tmp4)

          end if

! -----

! Calculate the turbulent mixing.

          if(tubopt.ge.1) then

            call s_turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,      &
     &                     ni,nj,nk,j31,j32,jcb,qicep(0,0,1,3),         &
     &                     qifrc(0,0,1,3),rkh8u,rkh8v,rkv8w,            &
     &                     tmp1,tmp2,h3,tmp3,tmp4,tmp5)

            call s_turbs(idtrnopt,idmpopt,idmfcopt,                     &
     &                   iddxiv,iddyiv,iddziv,iddmpvar,ni,nj,nk,'xx',   &
     &                   j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,        &
     &                   tmp1,tmp2,h3,qifrc(0,0,1,3),tmp3,tmp4,tmp5)

          end if

! -----

! Perform the analysis nudging to GPV.

          if(nggdmp.gt.0.e0.and.nggvar(7:7).eq.'o') then

            call s_s2gpv(idgpvvar,7,nggdmp,gtinc,                       &
     &                   ni,nj,nk,rst,qicep(0,0,1,3),                   &
     &                   qigpv(0,0,1,3),qitd(0,0,1,3),qifrc(0,0,1,3))

          end if

! -----

! Calculate the lateral sponge damping.

          if(lspopt.ge.1.and.lspvar(7:7).eq.'o') then

            call s_lsps(idgpvvar,idlspopt,idwdnews,idlsnews,idlspsmt,   &
     &                  7,gtinc,ni,nj,nk,rst,qicep(0,0,1,3),rbcxy,      &
     &                  qigpv(0,0,1,3),qitd(0,0,1,3),qifrc(0,0,1,3),    &
     &                  tmp1)

          end if

! -----

! Calculate the vertical sponge damping.

          if(vspopt.ge.1.and.vspvar(7:7).eq.'o') then

            call s_vsps(idgpvvar,idvspopt,7,ksp0,                       &
     &                  gtinc,ni,nj,nk,rst,qicep(0,0,1,3),rbct,         &
     &                  qigpv(0,0,1,3),qitd(0,0,1,3),qifrc(0,0,1,3))

          end if

! -----

! Perform the analysis nudging to radar data.

          if(ngropt.ge.1.and.ngrvar(5:5).eq.'o') then

            call s_qp2rdr(idngropt,ngrdmp,rtinc,                        &
     &                    ni,nj,nk,rst,qicep(0,0,1,3),                  &
     &                    qirdr(0,0,1,3),qirtd(0,0,1,3),qifrc(0,0,1,3))

          end if

! -----

!! -----

!! For the hail.

          if(haiopt.eq.1) then

! Calculate the advection.

            call s_advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,     &
     &                 iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc,&
     &                 qice(0,0,1,4),qifrc(0,0,1,4),tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

            if(mod(smtopt,10).eq.1) then

              call s_smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,                 &
     &                      rbr,qicep(0,0,1,4),qifrc(0,0,1,4),tmp1)

            end if

! -----

! Calculate the 4th order smoothing.

            if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

              call s_smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth, &
     &                      idsmhcoe,idsmvcoe,ni,nj,nk,rbr,             &
     &                      qicep(0,0,1,4),qifrc(0,0,1,4),              &
     &                      tmp1,tmp2,tmp3,tmp4,tmp5)

            end if

! -----

! Calculate the non linear smoothing.

            if(smtopt.ge.11) then

              call s_nlsms(idnlhcoe,idnlvcoe,.001e0,ni,nj,nk,           &
     &                     rbr,qicep(0,0,1,4),qifrc(0,0,1,4),           &
     &                     tmp1,tmp2,tmp3,tmp4)

            end if

! -----

! Calculate the turbulent mixing.

            if(tubopt.ge.1) then

              call s_turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,    &
     &                       ni,nj,nk,j31,j32,jcb,qicep(0,0,1,4),       &
     &                       qifrc(0,0,1,4),rkh8u,rkh8v,rkv8w,          &
     &                       tmp1,tmp2,h3,tmp3,tmp4,tmp5)

              call s_turbs(idtrnopt,idmpopt,idmfcopt,                   &
     &                     iddxiv,iddyiv,iddziv,iddmpvar,ni,nj,nk,'xx', &
     &                     j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,      &
     &                     tmp1,tmp2,h3,qifrc(0,0,1,4),tmp3,tmp4,tmp5)

            end if

! -----

! Perform the analysis nudging to GPV.

            if(nggdmp.gt.0.e0.and.nggvar(7:7).eq.'o') then

              call s_s2gpv(idgpvvar,8,nggdmp,gtinc,                     &
     &                     ni,nj,nk,rst,qicep(0,0,1,4),                 &
     &                     qigpv(0,0,1,4),qitd(0,0,1,4),qifrc(0,0,1,4))

            end if

! -----

! Calculate the lateral sponge damping.

            if(lspopt.ge.1.and.lspvar(7:7).eq.'o') then

              call s_lsps(idgpvvar,idlspopt,idwdnews,idlsnews,idlspsmt, &
     &                    8,gtinc,ni,nj,nk,rst,qicep(0,0,1,4),rbcxy,    &
     &                    qigpv(0,0,1,4),qitd(0,0,1,4),qifrc(0,0,1,4),  &
     &                    tmp1)

            end if

! -----

! Calculate the vertical sponge damping.

            if(vspopt.ge.1.and.vspvar(7:7).eq.'o') then

              call s_vsps(idgpvvar,idvspopt,8,ksp0,                     &
     &                    gtinc,ni,nj,nk,rst,qicep(0,0,1,4),rbct,       &
     &                    qigpv(0,0,1,4),qitd(0,0,1,4),qifrc(0,0,1,4))

            end if

! -----

! Perform the analysis nudging to radar data.

            if(ngropt.ge.1.and.ngrvar(5:5).eq.'o') then

              call s_qp2rdr(idngropt,ngrdmp,rtinc,                      &
     &                     ni,nj,nk,rst,qicep(0,0,1,4),                 &
     &                     qirdr(0,0,1,4),qirtd(0,0,1,4),qifrc(0,0,1,4))

            end if

! -----

          end if

!! -----

        end if

!!! -----

!!!! Calculate the forcing terms in the ice concentrations equations.

!! For the cloud ice concentrations.

        if(abs(cphopt).ge.2) then

! Calculate the advection.

          call s_advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,       &
     &                iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc, &
     &                nice(0,0,1,1),nifrc(0,0,1,1),tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

          if(mod(smtopt,10).eq.1) then

            call s_smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,rbr,nicep(0,0,1,1),&
     &                    nifrc(0,0,1,1),tmp1)

          end if

! -----

! Calculate the 4th order smoothing.

          if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

            call s_smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth,   &
     &                    idsmhcoe,idsmvcoe,ni,nj,nk,rbr,               &
     &                    nicep(0,0,1,1),nifrc(0,0,1,1),                &
     &                    tmp1,tmp2,tmp3,tmp4,tmp5)

          end if

! -----

! Calculate the non linear smoothing.

          if(smtopt.ge.11) then

            call s_nlsms(idnlhcoe,idnlvcoe,1.e20,ni,nj,nk,              &
     &                   rbr,nicep(0,0,1,1),nifrc(0,0,1,1),             &
     &                   tmp1,tmp2,tmp3,tmp4)

          end if

! -----

! Calculate the turbulent mixing.

          if(tubopt.ge.1) then

            call s_turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,      &
     &                     ni,nj,nk,j31,j32,jcb,nicep(0,0,1,1),         &
     &                     nifrc(0,0,1,1),rkh8u,rkh8v,rkv8w,            &
     &                     tmp1,tmp2,h3,tmp3,tmp4,tmp5)

            call s_turbs(idtrnopt,idmpopt,idmfcopt,                     &
     &                   iddxiv,iddyiv,iddziv,iddmpvar,ni,nj,nk,'xx',   &
     &                   j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,        &
     &                   tmp1,tmp2,h3,nifrc(0,0,1,1),tmp3,tmp4,tmp5)

          end if

! -----

! Calculate the lateral sponge damping.

          if(lspopt.ge.1.and.lspvar(7:7).eq.'o') then

            call s_lsps0(idlspopt,idwdnews,idlsnews,idlspsmt,           &
     &                   ni,nj,nk,rst,nicep(0,0,1,1),rbcxy,             &
     &                   nifrc(0,0,1,1),tmp1)

          end if

! -----

! Calculate the vertical sponge damping.

          if(vspopt.ge.1.and.vspvar(7:7).eq.'o') then

            call s_vsps0(ksp0,ni,nj,nk,rst,nicep(0,0,1,1),rbct,         &
     &                   nifrc(0,0,1,1))

          end if

! -----

        end if

!! -----

!!! For the snow, graupel and hail concentrations.

        if(abs(cphopt).ge.3) then

!! For the snow concentrations.

! Calculate the advection.

          call s_advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,       &
     &                iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc, &
     &                nice(0,0,1,2),nifrc(0,0,1,2),tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

          if(mod(smtopt,10).eq.1) then

            call s_smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,rbr,nicep(0,0,1,2),&
     &                    nifrc(0,0,1,2),tmp1)

          end if

! -----

! Calculate the 4th order smoothing.

          if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

            call s_smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth,   &
     &                    idsmhcoe,idsmvcoe,ni,nj,nk,rbr,               &
     &                    nicep(0,0,1,2),nifrc(0,0,1,2),                &
     &                    tmp1,tmp2,tmp3,tmp4,tmp5)

          end if

! -----

! Calculate the non linear smoothing.

          if(smtopt.ge.11) then

            call s_nlsms(idnlhcoe,idnlvcoe,1.e9,ni,nj,nk,               &
     &                   rbr,nicep(0,0,1,2),nifrc(0,0,1,2),             &
     &                   tmp1,tmp2,tmp3,tmp4)

          end if

! -----

! Calculate the turbulent mixing.

          if(tubopt.ge.1) then

            call s_turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,      &
     &                     ni,nj,nk,j31,j32,jcb,nicep(0,0,1,2),         &
     &                     nifrc(0,0,1,2),rkh8u,rkh8v,rkv8w,            &
     &                     tmp1,tmp2,h3,tmp3,tmp4,tmp5)

            call s_turbs(idtrnopt,idmpopt,idmfcopt,                     &
     &                   iddxiv,iddyiv,iddziv,iddmpvar,ni,nj,nk,'xx',   &
     &                   j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,        &
     &                   tmp1,tmp2,h3,nifrc(0,0,1,2),tmp3,tmp4,tmp5)

          end if

! -----

! Calculate the lateral sponge damping.

          if(lspopt.ge.1.and.lspvar(7:7).eq.'o') then

            call s_lsps0(idlspopt,idwdnews,idlsnews,idlspsmt,           &
     &                   ni,nj,nk,rst,nicep(0,0,1,2),rbcxy,             &
     &                   nifrc(0,0,1,2),tmp1)

          end if

! -----

! Calculate the vertical sponge damping.

          if(vspopt.ge.1.and.vspvar(7:7).eq.'o') then

            call s_vsps0(ksp0,ni,nj,nk,rst,nicep(0,0,1,2),rbct,         &
     &                   nifrc(0,0,1,2))

          end if

! -----

!! -----

!! For the graupel concentrations.

! Calculate the advection.

          call s_advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,       &
     &                iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc, &
     &                nice(0,0,1,3),nifrc(0,0,1,3),tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

          if(mod(smtopt,10).eq.1) then

            call s_smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,rbr,nicep(0,0,1,3),&
     &                    nifrc(0,0,1,3),tmp1)

          end if

! -----

! Calculate the 4th order smoothing.

          if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

            call s_smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth,   &
     &                    idsmhcoe,idsmvcoe,ni,nj,nk,rbr,               &
     &                    nicep(0,0,1,3),nifrc(0,0,1,3),                &
     &                    tmp1,tmp2,tmp3,tmp4,tmp5)

          end if

! -----

! Calculate the non linear smoothing.

          if(smtopt.ge.11) then

            call s_nlsms(idnlhcoe,idnlvcoe,1.e9,ni,nj,nk,               &
     &                   rbr,nicep(0,0,1,3),nifrc(0,0,1,3),             &
     &                   tmp1,tmp2,tmp3,tmp4)

          end if

! -----

! Calculate the turbulent mixing.

          if(tubopt.ge.1) then

            call s_turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,      &
     &                     ni,nj,nk,j31,j32,jcb,nicep(0,0,1,3),         &
     &                     nifrc(0,0,1,3),rkh8u,rkh8v,rkv8w,            &
     &                     tmp1,tmp2,h3,tmp3,tmp4,tmp5)

            call s_turbs(idtrnopt,idmpopt,idmfcopt,                     &
     &                   iddxiv,iddyiv,iddziv,iddmpvar,ni,nj,nk,'xx',   &
     &                   j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,        &
     &                   tmp1,tmp2,h3,nifrc(0,0,1,3),tmp3,tmp4,tmp5)

          end if

! -----

! Calculate the lateral sponge damping.

          if(lspopt.ge.1.and.lspvar(7:7).eq.'o') then

            call s_lsps0(idlspopt,idwdnews,idlsnews,idlspsmt,           &
     &                   ni,nj,nk,rst,nicep(0,0,1,3),rbcxy,             &
     &                   nifrc(0,0,1,3),tmp1)

          end if

! -----

! Calculate the vertical sponge damping.

          if(vspopt.ge.1.and.vspvar(7:7).eq.'o') then

            call s_vsps0(ksp0,ni,nj,nk,rst,nicep(0,0,1,3),rbct,         &
     &                   nifrc(0,0,1,3))

          end if

! -----

!! -----

!! For the hail concentrations.

          if(haiopt.eq.1) then

! Calculate the advection.

            call s_advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,     &
     &                 iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc,&
     &                 nice(0,0,1,4),nifrc(0,0,1,4),tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

            if(mod(smtopt,10).eq.1) then

              call s_smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,                 &
     &                      rbr,nicep(0,0,1,4),nifrc(0,0,1,4),tmp1)

            end if

! -----

! Calculate the 4th order smoothing.

            if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

              call s_smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth, &
     &                      idsmhcoe,idsmvcoe,ni,nj,nk,rbr,             &
     &                      nicep(0,0,1,4),nifrc(0,0,1,4),              &
     &                      tmp1,tmp2,tmp3,tmp4,tmp5)

            end if

! -----

! Calculate the non linear smoothing.

            if(smtopt.ge.11) then

              call s_nlsms(idnlhcoe,idnlvcoe,1.e9,ni,nj,nk,             &
     &                     rbr,nicep(0,0,1,4),nifrc(0,0,1,4),           &
     &                     tmp1,tmp2,tmp3,tmp4)

            end if

! -----

! Calculate the turbulent mixing.

            if(tubopt.ge.1) then

              call s_turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,    &
     &                       ni,nj,nk,j31,j32,jcb,nicep(0,0,1,4),       &
     &                       nifrc(0,0,1,4),rkh8u,rkh8v,rkv8w,          &
     &                       tmp1,tmp2,h3,tmp3,tmp4,tmp5)

              call s_turbs(idtrnopt,idmpopt,idmfcopt,                   &
     &                     iddxiv,iddyiv,iddziv,iddmpvar,ni,nj,nk,'xx', &
     &                     j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,      &
     &                     tmp1,tmp2,h3,nifrc(0,0,1,4),tmp3,tmp4,tmp5)

            end if

! -----

! Calculate the lateral sponge damping.

            if(lspopt.ge.1.and.lspvar(7:7).eq.'o') then

              call s_lsps0(idlspopt,idwdnews,idlsnews,idlspsmt,         &
     &                     ni,nj,nk,rst,nicep(0,0,1,4),rbcxy,           &
     &                     nifrc(0,0,1,4),tmp1)

            end if

! -----

! Calculate the vertical sponge damping.

            if(vspopt.ge.1.and.vspvar(7:7).eq.'o') then

              call s_vsps0(ksp0,ni,nj,nk,rst,nicep(0,0,1,4),rbct,       &
     &                     nifrc(0,0,1,4))

            end if

! -----

          end if

!! -----

        end if

!!! -----

!!!! -----

!!!!! -----

!!!! Calculate the forcing terms in the water hydrometeor and the ice
!!!! hydrometeor and concentrations equations for bin categories.

      else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

!!! Calculate the forcing terms in the water hydrometeor equations.

        if(abs(cphopt).ge.11) then

!! For the water mixing ratio.

          do n=1,nqw

! Calculate the advection.

            call s_advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,     &
     &                 iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc,&
     &                 qwtr(0,0,1,n),qwfrc(0,0,1,n),tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

            if(mod(smtopt,10).eq.1) then

              call s_smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,                 &
     &                      rbr,qwtrp(0,0,1,n),qwfrc(0,0,1,n),tmp1)

            end if

! -----

! Calculate the 4th order smoothing.

            if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

              call s_smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth, &
     &                      idsmhcoe,idsmvcoe,ni,nj,nk,rbr,             &
     &                      qwtrp(0,0,1,n),qwfrc(0,0,1,n),              &
     &                      tmp1,tmp2,tmp3,tmp4,tmp5)

            end if

! -----

! Calculate the turbulent mixing.

            if(tubopt.ge.1) then

              call s_turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,    &
     &                       ni,nj,nk,j31,j32,jcb,qwtrp(0,0,1,n),       &
     &                       qwfrc(0,0,1,n),rkh8u,rkh8v,rkv8w,          &
     &                       tmp1,tmp2,h3,tmp3,tmp4,tmp5)

              call s_turbs(idtrnopt,idmpopt,idmfcopt,                   &
     &                     iddxiv,iddyiv,iddziv,iddmpvar,ni,nj,nk,'xx', &
     &                     j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,      &
     &                     tmp1,tmp2,h3,qwfrc(0,0,1,n),tmp3,tmp4,tmp5)

            end if

! -----

! Calculate the lateral sponge damping.

            if(lspopt.ge.1.and.lspvar(7:7).eq.'o') then

              call s_lsps0(idlspopt,idwdnews,idlsnews,idlspsmt,         &
     &                     ni,nj,nk,rst,qwtrp(0,0,1,n),rbcxy,           &
     &                     qwfrc(0,0,1,n),tmp1)

            end if

! -----

! Calculate the vertical sponge damping.

            if(vspopt.ge.1.and.vspvar(7:7).eq.'o') then

              call s_vsps0(ksp0,ni,nj,nk,rst,qwtrp(0,0,1,n),rbct,       &
     &                     qwfrc(0,0,1,n))

            end if

! -----

          end do

!! -----

!! For the water concentrations.

          do n=1,nnw

! Calculate the advection.

            call s_advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,     &
     &                 iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc,&
     &                 nwtr(0,0,1,n),nwfrc(0,0,1,n),tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

            if(mod(smtopt,10).eq.1) then

              call s_smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,                 &
     &                      rbr,nwtrp(0,0,1,n),nwfrc(0,0,1,n),tmp1)

            end if

! -----

! Calculate the 4th order smoothing.

            if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

              call s_smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth, &
     &                      idsmhcoe,idsmvcoe,ni,nj,nk,rbr,             &
     &                      nwtrp(0,0,1,n),nwfrc(0,0,1,n),              &
     &                      tmp1,tmp2,tmp3,tmp4,tmp5)

            end if

! -----

! Calculate the turbulent mixing.

            if(tubopt.ge.1) then

              call s_turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,    &
     &                       ni,nj,nk,j31,j32,jcb,nwtrp(0,0,1,n),       &
     &                       nwfrc(0,0,1,n),rkh8u,rkh8v,rkv8w,          &
     &                       tmp1,tmp2,h3,tmp3,tmp4,tmp5)

              call s_turbs(idtrnopt,idmpopt,idmfcopt,                   &
     &                     iddxiv,iddyiv,iddziv,iddmpvar,ni,nj,nk,'xx', &
     &                     j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,      &
     &                     tmp1,tmp2,h3,nwfrc(0,0,1,n),tmp3,tmp4,tmp5)

            end if

! -----

! Calculate the lateral sponge damping.

            if(lspopt.ge.1.and.lspvar(7:7).eq.'o') then

              call s_lsps0(idlspopt,idwdnews,idlsnews,idlspsmt,         &
     &                     ni,nj,nk,rst,nwtrp(0,0,1,n),rbcxy,           &
     &                     nwfrc(0,0,1,n),tmp1)

            end if

! -----

! Calculate the vertical sponge damping.

            if(vspopt.ge.1.and.vspvar(7:7).eq.'o') then

              call s_vsps0(ksp0,ni,nj,nk,rst,nwtrp(0,0,1,n),rbct,       &
     &                     nwfrc(0,0,1,n))

            end if

! -----

          end do

!! -----

        end if

!!! -----

!!! Calculate the forcing terms in the ice hydrometeor equations.

        if(abs(cphopt).eq.12) then

!! For the ice mixing ratio.

          do n=1,nqi

! Calculate the advection.

            call s_advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,     &
     &                 iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc,&
     &                 qice(0,0,1,n),qifrc(0,0,1,n),tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

            if(mod(smtopt,10).eq.1) then

              call s_smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,                 &
     &                      rbr,qicep(0,0,1,n),qifrc(0,0,1,n),tmp1)

            end if

! -----

! Calculate the 4th order smoothing.

            if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

              call s_smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth, &
     &                      idsmhcoe,idsmvcoe,ni,nj,nk,rbr,             &
     &                      qicep(0,0,1,n),qifrc(0,0,1,n),              &
     &                      tmp1,tmp2,tmp3,tmp4,tmp5)

            end if

! -----

! Calculate the turbulent mixing.

            if(tubopt.ge.1) then

              call s_turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,    &
     &                       ni,nj,nk,j31,j32,jcb,qicep(0,0,1,n),       &
     &                       qifrc(0,0,1,n),rkh8u,rkh8v,rkv8w,          &
     &                       tmp1,tmp2,h3,tmp3,tmp4,tmp5)

              call s_turbs(idtrnopt,idmpopt,idmfcopt,                   &
     &                     iddxiv,iddyiv,iddziv,iddmpvar,ni,nj,nk,'xx', &
     &                     j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,      &
     &                     tmp1,tmp2,h3,qifrc(0,0,1,n),tmp3,tmp4,tmp5)

            end if

! -----

! Calculate the lateral sponge damping.

            if(lspopt.ge.1.and.lspvar(7:7).eq.'o') then

              call s_lsps0(idlspopt,idwdnews,idlsnews,idlspsmt,         &
     &                     ni,nj,nk,rst,qicep(0,0,1,n),rbcxy,           &
     &                     qifrc(0,0,1,n),tmp1)

            end if

! -----

! Calculate the vertical sponge damping.

            if(vspopt.ge.1.and.vspvar(7:7).eq.'o') then

              call s_vsps0(ksp0,ni,nj,nk,rst,qicep(0,0,1,n),rbct,       &
     &                     qifrc(0,0,1,n))

            end if

! -----

          end do

!! -----

!! For the ice concentrations.

          do n=1,nni

! Calculate the advection.

            call s_advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,     &
     &                 iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc,&
     &                 nice(0,0,1,n),nifrc(0,0,1,n),tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

            if(mod(smtopt,10).eq.1) then

              call s_smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,                 &
     &                      rbr,nicep(0,0,1,n),nifrc(0,0,1,n),tmp1)

            end if

! -----

! Calculate the 4th order smoothing.

            if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

              call s_smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth, &
     &                      idsmhcoe,idsmvcoe,ni,nj,nk,rbr,             &
     &                      nicep(0,0,1,n),nifrc(0,0,1,n),              &
     &                      tmp1,tmp2,tmp3,tmp4,tmp5)

            end if

! -----

! Calculate the turbulent mixing.

            if(tubopt.ge.1) then

              call s_turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,    &
     &                       ni,nj,nk,j31,j32,jcb,nicep(0,0,1,n),       &
     &                       nifrc(0,0,1,n),rkh8u,rkh8v,rkv8w,          &
     &                       tmp1,tmp2,h3,tmp3,tmp4,tmp5)

              call s_turbs(idtrnopt,idmpopt,idmfcopt,                   &
     &                     iddxiv,iddyiv,iddziv,iddmpvar,ni,nj,nk,'xx', &
     &                     j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,      &
     &                     tmp1,tmp2,h3,nifrc(0,0,1,n),tmp3,tmp4,tmp5)

            end if

! -----

! Calculate the lateral sponge damping.

            if(lspopt.ge.1.and.lspvar(7:7).eq.'o') then

              call s_lsps0(idlspopt,idwdnews,idlsnews,idlspsmt,         &
     &                     ni,nj,nk,rst,nicep(0,0,1,n),rbcxy,           &
     &                     nifrc(0,0,1,n),tmp1)

            end if

! -----

! Calculate the vertical sponge damping.

            if(vspopt.ge.1.and.vspvar(7:7).eq.'o') then

              call s_vsps0(ksp0,ni,nj,nk,rst,nicep(0,0,1,n),rbct,       &
     &                     nifrc(0,0,1,n))

            end if

! -----

          end do

!! -----

        end if

!!! -----

      end if

!!!! -----

      end subroutine s_forceq

!-----7--------------------------------------------------------------7--

      end module m_forceq
