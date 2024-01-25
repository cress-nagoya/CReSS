!***********************************************************************
      module m_forcetke
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/10/12
!     Modification: 1999/11/01, 1999/11/19, 1999/11/24, 1999/11/30,
!                   2000/01/17, 2000/04/18, 2000/12/19, 2001/04/15,
!                   2001/06/06, 2001/11/20, 2002/01/15, 2002/04/02,
!                   2002/08/15, 2002/10/31, 2003/01/04, 2003/02/13,
!                   2003/03/13, 2003/03/21, 2003/04/30, 2003/05/19,
!                   2003/07/15, 2003/11/28, 2003/12/12, 2004/02/01,
!                   2004/03/05, 2004/05/31, 2004/06/10, 2004/08/01,
!                   2004/08/20, 2004/09/01, 2005/06/10, 2006/01/10,
!                   2006/02/13, 2006/04/03, 2006/05/12, 2006/06/21,
!                   2006/09/21, 2006/11/06, 2008/05/02, 2008/08/25,
!                   2009/01/30, 2009/02/27, 2011/09/22, 2011/11/10,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the forcing terms in the turbulent kinetic energy
!     equation.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_advs
      use m_buoytke
      use m_comindx
      use m_disptke
      use m_eddydif
      use m_eddyvisj
      use m_getcname
      use m_getiname
      use m_inichar
      use m_lsps0
      use m_sheartke
      use m_smoo2s
      use m_smoo4s
      use m_tkeflx
      use m_turbtke
      use m_vsps0

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: forcetke, s_forcetke

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface forcetke

        module procedure s_forcetke

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
      subroutine s_forcetke(fplspvar,fpvspvar,                          &
     &                      fplspopt,fpvspopt,fpsmtopt,ksp0,            &
     &                      ni,nj,nk,j31,j32,jcb,jcb8u,jcb8v,jcb8w,mf,  &
     &                      rmf,rmf8u,rmf8v,rbr,rst,rstxu,rstxv,rstxwc, &
     &                      tke,tkep,rbcxy,rbct,priv,nsq8w,ssq,rkh,rkv, &
     &                      tkefrc,h3,rkv8s,tmp1,tmp2,tmp3)
!***********************************************************************

! Input variables

      integer, intent(in) :: fplspvar
                       ! Formal parameter of unique index of lspvar

      integer, intent(in) :: fpvspvar
                       ! Formal parameter of unique index of vspvar

      integer, intent(in) :: fplspopt
                       ! Formal parameter of unique index of lspopt

      integer, intent(in) :: fpvspopt
                       ! Formal parameter of unique index of vspopt

      integer, intent(in) :: fpsmtopt
                       ! Formal parameter of unique index of smtopt

      integer, intent(in) :: ksp0(1:2)
                       ! Index of lowest vertical sponge level

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

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

      real, intent(in) :: tke(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at present

      real, intent(in) :: tkep(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at past

      real, intent(in) :: rbcxy(1:ni,1:nj)
                       ! Relaxed lateral sponge damping coefficients

      real, intent(in) :: rbct(1:ni,1:nj,1:nk,1:2)
                       ! Relaxed top sponge damping coefficients

      real, intent(in) :: priv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of turbulent Prandtl number

      real, intent(in) :: nsq8w(0:ni+1,0:nj+1,1:nk)
                       ! Half value of Brunt Vaisala frequency squared
                       ! at w points

! Input and output variables

      real, intent(inout) :: ssq(0:ni+1,0:nj+1,1:nk)
                       ! Magnitude of deformation squared

      real, intent(inout) :: rkh(0:ni+1,0:nj+1,1:nk)
                       ! rbr x horizontal eddy viscosity or diffusivity

      real, intent(inout) :: rkv(0:ni+1,0:nj+1,1:nk)
                       ! rbr x vertical eddy viscosity or diffusivity

! Output variable

      real, intent(out) :: tkefrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in turbulent kinetic energy equation

! Internal shared variables

      character(len=108) lspvar
                       ! Control flag of
                       ! lateral sponge damped variables

      character(len=108) vspvar
                       ! Control flag of
                       ! vertical sponge damped variables

      integer lspopt   ! Option for lateral sponge damping
      integer vspopt   ! Option for vertical sponge damping

      integer smtopt   ! Option for numerical smoothing

      real, intent(inout) :: h3(0:ni+1,0:nj+1,1:nk)
                       ! z components of turbulent fluxes

      real, intent(inout) :: rkv8s(0:ni+1,0:nj+1,1:nk)
                       ! Half value of rbr x vertical eddy diffusivity
                       ! at scalar points

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Remark

!     ssq: This variable is also temporary, because it is not used
!          again.

!-----7--------------------------------------------------------------7--

! Initialize the character variables.

      call inichar(lspvar)
      call inichar(vspvar)

! -----

! Get the required namelist variables.

      call getcname(fplspvar,lspvar)
      call getcname(fpvspvar,vspvar)
      call getiname(fplspopt,lspopt)
      call getiname(fpvspopt,vspopt)
      call getiname(fpsmtopt,smtopt)

! -----

! Calculate the advection.

      call advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,             &
     &          iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc,tke,   &
     &          tkefrc,tmp1,tmp2,tmp3,rkv8s)

! -----

! Calculate the shear production.

      call sheartke(ni,nj,nk,jcb,ssq,rkv,tkefrc)

! -----

! Calculate the 2nd order smoothing.

      if(mod(smtopt,10).eq.1) then

        call smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,rbr,tkep,tkefrc,tmp1)

      end if

! -----

! Calculate the 4th order smoothing.

      if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

        call smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth,         &
     &              idsmhcoe,idsmvcoe,ni,nj,nk,rbr,tkep,tkefrc,         &
     &              tmp1,tmp2,tmp3,ssq,rkv8s)

      end if

! -----

! Calculate the dissipation.

      call disptke(idmpopt,idmfcopt,idisoopt,iddx,iddy,iddz,ni,nj,nk,   &
     &             jcb,rmf,rst,priv,tkep,tkefrc)

! -----

! Calculate the turbulent mixing.

      call eddyvisj(idmfcopt,ni,nj,nk,jcb,mf,rkh,rkv)

      call tkeflx(idtrnopt,idmpopt,idmfcopt,iddxiv,iddyiv,iddziv,       &
     &            ni,nj,nk,j31,j32,jcb,rmf,tkep,rkh,rkv,                &
     &            tmp1,tmp2,h3,tmp3,ssq,rkv8s)

      call turbtke(idtrnopt,idmpopt,idmfcopt,iddxiv,iddyiv,iddziv,      &
     &             ni,nj,nk,j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,     &
     &             tmp1,tmp2,h3,tkefrc,tmp3,ssq,rkv8s)

! -----

! Calculate the lateral sponge damping.

      if(lspopt.ge.1.and.lspvar(10:10).eq.'o') then

        call lsps0(idlspopt,idwdnews,idlsnews,idlspsmt,ni,nj,nk,rst,    &
     &             tkep,rbcxy,tkefrc,tmp1)

      end if

! -----

! Calculate the vertical sponge damping.

      if(vspopt.ge.1.and.vspvar(10:10).eq.'o') then

        call vsps0(ksp0,ni,nj,nk,rst,tkep,rbct,tkefrc)

      end if

! -----

! Calculate the eddy diffusivity.

      call eddydif(idmfcopt,idtubopt,idisoopt,ni,nj,nk,jcb,mf,priv,     &
     &             rkh,rkv,rkv8s)

! -----

! Calculate the buoyancy production.

      call buoytke(ni,nj,nk,jcb8w,nsq8w,rkv8s,tkefrc,tmp1)

! -----

      end subroutine s_forcetke

!-----7--------------------------------------------------------------7--

      end module m_forcetke
