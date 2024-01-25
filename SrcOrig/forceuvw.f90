!***********************************************************************
      module m_forceuvw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/05/20,
!                   1999/06/07, 1999/07/05, 1999/07/19, 1999/08/03,
!                   1999/08/09, 1999/08/18, 1999/08/23, 1999/09/16,
!                   1999/09/30, 1999/10/07, 1999/10/12, 1999/10/21,
!                   1999/11/01, 1999/11/19, 1999/11/30, 2000/01/17,
!                   2000/02/07, 2000/04/18, 2000/12/19, 2001/01/15,
!                   2001/03/13, 2001/04/15, 2001/05/29, 2001/06/06,
!                   2001/07/13, 2001/08/07, 2001/11/20, 2002/01/15,
!                   2002/04/02, 2002/06/18, 2002/07/23, 2002/08/15,
!                   2002/09/02, 2002/09/09, 2002/10/31, 2002/11/11,
!                   2002/12/11, 2003/01/04, 2003/01/20, 2003/02/13,
!                   2003/03/13, 2003/03/21, 2003/04/30, 2003/05/19,
!                   2003/07/15, 2003/10/10, 2003/11/28, 2003/12/12,
!                   2004/02/01, 2004/03/05, 2004/04/15, 2004/05/31,
!                   2004/06/10, 2004/08/01, 2004/08/20, 2005/08/05,
!                   2006/01/10, 2006/02/13, 2006/04/03, 2006/05/12,
!                   2006/06/21, 2006/10/20, 2006/11/06, 2007/01/20,
!                   2007/05/07, 2007/07/30, 2008/05/02, 2008/06/09,
!                   2008/08/25, 2008/12/11, 2009/02/27, 2009/03/23,
!                   2011/09/22, 2011/11/10, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the forcing terms in the velocity equations.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_advuvw
      use m_buoywb
      use m_comindx
      use m_coriuv
      use m_coriuvw
      use m_curveuvw
      use m_getiname
      use m_lspuvw
      use m_nlsmuvw
      use m_phy2cnt
      use m_rstuvwc
      use m_setcst3d
      use m_smoo2uvw
      use m_smoo4uvw
      use m_turbdrv
      use m_uvw2gpv
      use m_uvw2rdr
      use m_vspuvw

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: forceuvw, s_forceuvw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface forceuvw

        module procedure s_forceuvw

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_forceuvw(fplspopt,fpvspopt,fpcoropt,fpcrvopt,        &
     &                      fpbuyopt,fpsmtopt,fpadvopt,fptubopt,        &
     &                      fmois,ksp0,dtb,nggdmp,ngrdmp,gtinc,         &
     &                      ni,nj,nk,zph,j31,j32,jcb,jcb8u,jcb8v,jcb8w, &
     &                      mf,mf8u,mf8v,rmf,rmf8u,rmf8v,fc,ubr,vbr,    &
     &                      pbr,ptbr,qvbr,rbr,rst,rst8u,rst8v,rst8w,    &
     &                      u,up,v,vp,w,wp,ppp,ptp,ptpp,qv,qvp,tkep,    &
     &                      rbcx,rbcy,rbcxy,rbct,ugpv,utd,vgpv,vtd,     &
     &                      wgpv,wtd,urdr,vrdr,wrdr,qall,qallp,         &
     &                      ufrc,vfrc,wfrc,wc,rstxu,rstxv,rstxwc,       &
     &                      rkh,rkv,priv,ssq,nsq8w,t13,t23,t33,         &
     &                      tmp1,tmp2)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: fplspopt
                       ! Formal parameter of unique index of lspopt

      integer, intent(in) :: fpvspopt
                       ! Formal parameter of unique index of vspopt

      integer, intent(in) :: fpcoropt
                       ! Formal parameter of unique index of coropt

      integer, intent(in) :: fpcrvopt
                       ! Formal parameter of unique index of crvopt

      integer, intent(in) :: fpbuyopt
                       ! Formal parameter of unique index of buyopt

      integer, intent(in) :: fpsmtopt
                       ! Formal parameter of unique index of smtopt

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

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

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: nggdmp
                       ! Analysis nudging damping coefficient for GPV

      real, intent(in) :: ngrdmp(1:2)
                       ! Analysis nudging damping coefficient for radar

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

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

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: mf8u(0:ni+1,0:nj+1)
                       ! Map scale factors at u points

      real, intent(in) :: mf8v(0:ni+1,0:nj+1)
                       ! Map scale factors at v points

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(in) :: rmf8u(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at u points

      real, intent(in) :: rmf8v(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at v points

      real, intent(in) :: fc(0:ni+1,0:nj+1,1:2)
                       ! 0.25 x Coriolis parameters

      real, intent(in) :: ubr(0:ni+1,0:nj+1,1:nk)
                       ! Base state x components of velocity

      real, intent(in) :: vbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state y components of velocity

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

      real, intent(in) :: rst8u(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at u points

      real, intent(in) :: rst8v(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at v points

      real, intent(in) :: rst8w(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at w points

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at present

      real, intent(in) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at present

      real, intent(in) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

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

      real, intent(in) :: tkep(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at past

      real, intent(in) :: rbcx(1:ni)
                       ! Relaxed lateral sponge damping coefficients
                       ! in x direction

      real, intent(in) :: rbcy(1:nj)
                       ! Relaxed lateral sponge damping coefficients
                       ! in y direction

      real, intent(in) :: rbcxy(1:ni,1:nj)
                       ! Relaxed lateral sponge damping coefficients

      real, intent(in) :: rbct(1:ni,1:nj,1:nk,1:2)
                       ! Relaxed top sponge damping coefficients

      real, intent(in) :: ugpv(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: utd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! x components of velocity of GPV data

      real, intent(in) :: vgpv(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: vtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! y components of velocity of GPV data

      real, intent(in) :: wgpv(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: wtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! z components of velocity of GPV data

      real, intent(in) :: urdr(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity of radar data
                       ! at marked time

      real, intent(in) :: vrdr(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity of radar data
                       ! at marked time

      real, intent(in) :: wrdr(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity of radar data
                       ! at marked time

      real, intent(in) :: qall(0:ni+1,0:nj+1,1:nk)
                       ! Total water and ice mixing ratio at present

      real, intent(in) :: qallp(0:ni+1,0:nj+1,1:nk)
                       ! Total water and ice mixing ratio at past

! Input and output variables

      real, intent(inout) :: ufrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in u equation

      real, intent(inout) :: vfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in v equation

! Output variables

      real, intent(out) :: wfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in w equation

      real, intent(out) :: wc(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity

      real, intent(out) :: rstxu(0:ni+1,0:nj+1,1:nk)
                       ! u x base state density x Jacobian

      real, intent(out) :: rstxv(0:ni+1,0:nj+1,1:nk)
                       ! v x base state density x Jacobian

      real, intent(out) :: rstxwc(0:ni+1,0:nj+1,1:nk)
                       ! wc x base state density x Jacobian

      real, intent(out) :: rkh(0:ni+1,0:nj+1,1:nk)
                       ! rbr x horizontal eddy diffusivity

      real, intent(out) :: rkv(0:ni+1,0:nj+1,1:nk)
                       ! rbr x vertical eddy diffusivity

      real, intent(out) :: priv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of turbulent Prandtl number

      real, intent(out) :: ssq(0:ni+1,0:nj+1,1:nk)
                       ! Magnitude of deformation squared

      real, intent(out) :: nsq8w(0:ni+1,0:nj+1,1:nk)
                       ! Half value of Brunt Vaisala frequency squared
                       ! at w points

! Internal shared variables

      integer lspopt   ! Option for lateral sponge damping
      integer vspopt   ! Option for vertical sponge damping

      integer coropt   ! Option for coriolis force
      integer crvopt   ! Option for earth curvature
      integer buyopt   ! Option for buoyancy calculation
      integer smtopt   ! Option for numerical smoothing
      integer advopt   ! Option for advection scheme
      integer tubopt   ! Option for turbulent mixing

      real, intent(inout) :: t13(0:ni+1,0:nj+1,1:nk)
                       ! x-z components of stress tensor

      real, intent(inout) :: t23(0:ni+1,0:nj+1,1:nk)
                       ! y-z components of stress tensor

      real, intent(inout) :: t33(0:ni+1,0:nj+1,1:nk)
                       ! z-z components of stress tensor

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Remark

!     wc,rstxu,rstxv,rstxwc: These variables are also temporary.

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fplspopt,lspopt)
      call getiname(fpvspopt,vspopt)
      call getiname(fpcoropt,coropt)
      call getiname(fpcrvopt,crvopt)
      call getiname(fpbuyopt,buyopt)
      call getiname(fpsmtopt,smtopt)
      call getiname(fpadvopt,advopt)
      call getiname(fptubopt,tubopt)

! -----

! Calculate the turbulent mixing.

      if(tubopt.ge.1) then

        call turbdrv(fmois,dtb,ni,nj,nk,zph,j31,j32,                    &
     &               jcb,jcb8u,jcb8v,jcb8w,mf,mf8u,mf8v,rmf,rmf8u,rmf8v,&
     &               pbr,ptbr,rbr,up,vp,wp,ppp,ptpp,qvp,tkep,qallp,     &
     &               ufrc,vfrc,wfrc,rkh,rkv,priv,ssq,nsq8w,             &
     &               t13,t23,t33,rstxu,rstxv,rstxwc,                    &
     &               tmp1,tmp2,wc)

      else

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ufrc)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,vfrc)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,wfrc)

      end if

! -----

! Calculate the zeta components of contravariant velocity.

      if(advopt.le.3) then

        call s_phy2cnt(idsthopt,idtrnopt,idmpopt,idmfcopt,idoneopt,     &
     &                 ni,nj,nk,j31,j32,jcb8w,mf,u,v,w,wc,              &
     &                 rstxu,rstxv,rstxwc)

      end if

! -----

! Calculate the 2nd order smoothing.

      if(mod(smtopt,10).eq.1) then

        call smoo2uvw(idsmhcoe,idsmvcoe,ni,nj,nk,jcb8w,ubr,vbr,rbr,     &
     &                rst8w,up,vp,wp,ufrc,vfrc,wfrc,tmp1)

      end if

! -----

! Calculate the 4th order smoothing.

      if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

        call smoo4uvw(idsmtopt,idiwest,idieast,idjsouth,idjnorth,       &
     &                idsmhcoe,idsmvcoe,ni,nj,nk,jcb8u,jcb8v,jcb8w,     &
     &                ubr,vbr,rst8u,rst8v,rst8w,up,vp,wp,ufrc,vfrc,wfrc,&
     &                rstxu,rstxv,rstxwc,tmp1,tmp2)

      end if

! -----

! Calculate the non linear smoothing.

      if(smtopt.ge.11) then

        call nlsmuvw(idnlhcoe,idnlvcoe,2.e0,ni,nj,nk,jcb8w,ubr,vbr,rbr, &
     &               rst8w,up,vp,wp,ufrc,vfrc,wfrc,rstxu,rstxv,rstxwc,  &
     &               tmp1)

      end if

! -----

! Perform the analysis nudging to GPV.

      if(nggdmp.gt.0.e0) then

        call uvw2gpv(idnggvar,nggdmp,gtinc,ni,nj,nk,rst8u,rst8v,rst8w,  &
     &               up,vp,wp,ugpv,utd,vgpv,vtd,wgpv,wtd,ufrc,vfrc,wfrc)

      end if

! -----

! Calculate the lateral sponge damping.

      if(lspopt.ge.1) then

        call s_lspuvw(idgpvvar,idlspvar,idlspopt,                       &
     &                idwdnews,idwdnorm,idlsnews,idlsnorm,idlspsmt,     &
     &                gtinc,ni,nj,nk,ubr,vbr,rst8u,rst8v,rst8w,up,vp,wp,&
     &                rbcx,rbcy,rbcxy,ugpv,utd,vgpv,vtd,wgpv,wtd,       &
     &                ufrc,vfrc,wfrc,tmp1,tmp2)

      end if

! -----

! Calculate the vertical sponge damping.

      if(vspopt.ge.1) then

        call vspuvw(idgpvvar,idvspvar,idvspopt,ksp0,gtinc,              &
     &              ni,nj,nk,ubr,vbr,rst8u,rst8v,rst8w,up,vp,wp,        &
     &              rbct,ugpv,utd,vgpv,vtd,wgpv,wtd,ufrc,vfrc,wfrc,tmp1)

      end if

! -----

! Perform the analysis nudging to radar data.

      if(ngrdmp(2).gt.0.e0) then

        call uvw2rdr(idngrvar,ngrdmp,ni,nj,nk,rst8u,rst8v,rst8w,        &
     &               up,vp,wp,urdr,vrdr,wrdr,ufrc,vfrc,wfrc)

      end if

! -----

! Calculate the Coriolis force.

      if(coropt.eq.1) then

        if(advopt.le.3) then

         call coriuv(ni,nj,nk,fc,rst,u,v,ufrc,vfrc,tmp1)

        else

         call coriuv(ni,nj,nk,fc,rst,up,vp,ufrc,vfrc,tmp1)

        end if

      else if(coropt.eq.2) then

        if(advopt.le.3) then

         call coriuvw(ni,nj,nk,fc,rst,u,v,w,ufrc,vfrc,wfrc,tmp1,tmp2)

        else

         call coriuvw(ni,nj,nk,fc,rst,up,vp,wp,ufrc,vfrc,wfrc,tmp1,tmp2)

        end if

      end if

! -----

! Calculate the curvature of earth.

      if(crvopt.eq.1) then

        if(advopt.le.3) then

          call curveuvw(idmpopt,idmfcopt,ni,nj,nk,rmf8u,rmf8v,rst,      &
     &                  u,v,w,ufrc,vfrc,wfrc,rstxu,rstxv,rstxwc,        &
     &                  tmp1,tmp2)

        else

          call curveuvw(idmpopt,idmfcopt,ni,nj,nk,rmf8u,rmf8v,rst,      &
     &                 up,vp,wp,ufrc,vfrc,wfrc,rstxu,rstxv,rstxwc,      &
     &                 tmp1,tmp2)

        end if

      end if

! -----

! Calculate the buoyancy in the large time steps.

      if(buyopt.eq.1) then

        if(advopt.le.3) then

          call buoywb(idgwmopt,idcphopt,fmois,ni,nj,nk,ptbr,qvbr,rst,   &
     &                ptp,qv,qall,wfrc,tmp1,tmp2)

        else

          call buoywb(idgwmopt,idcphopt,fmois,ni,nj,nk,ptbr,qvbr,rst,   &
     &                ptpp,qvp,qallp,wfrc,tmp1,tmp2)

        end if

      end if

! -----

! The base state density x the Jacobian is multiplyed by u, v and w.

      if(advopt.le.3) then

        call rstuvwc(idmpopt,idmfcopt,idiwest,idieast,idjsouth,idjnorth,&
     &               ni,nj,nk,mf8u,mf8v,rst8u,rst8v,rst8w,              &
     &               u,v,wc,rstxu,rstxv,rstxwc)

      end if

! -----

! Calculate the advection.

      if(advopt.le.3) then

        call advuvw(idadvopt,idiwest,idieast,idjsouth,idjnorth,         &
     &              iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc,   &
     &              u,v,w,ufrc,vfrc,wfrc,t13,t23,t33,tmp1,tmp2)

      end if

! -----

      end subroutine s_forceuvw

!-----7--------------------------------------------------------------7--

      end module m_forceuvw
