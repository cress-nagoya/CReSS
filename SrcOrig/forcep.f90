!***********************************************************************
      module m_forcep
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/05/20,
!                   1999/06/07, 1999/07/05, 1999/08/03, 1999/08/09,
!                   1999/08/18, 1999/08/23, 1999/09/30, 1999/10/07,
!                   1999/10/12, 1999/11/01, 1999/11/19, 1999/11/24,
!                   2000/01/17, 2000/02/02, 2000/04/18, 2000/12/18,
!                   2001/01/15, 2001/03/13, 2001/04/15, 2001/05/29,
!                   2001/06/06, 2001/07/13, 2001/08/07, 2001/11/20,
!                   2002/04/02, 2002/06/18, 2002/07/23, 2002/08/15,
!                   2002/09/09, 2002/10/31, 2002/12/11, 2003/01/04,
!                   2003/03/21, 2003/04/30, 2003/05/19, 2003/07/15,
!                   2003/09/01, 2003/10/10, 2003/12/12, 2004/03/05,
!                   2004/04/15, 2004/05/31, 2004/06/10, 2004/08/01,
!                   2004/08/20, 2006/01/10, 2006/02/13, 2006/04/03,
!                   2006/05/12, 2006/06/21, 2006/09/21, 2006/11/06,
!                   2007/05/07, 2007/07/30, 2008/05/02, 2008/08/25,
!                   2008/12/11, 2009/02/27, 2009/03/23, 2011/09/22,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the forcing term in the pressure equation.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_advp
      use m_comindx
      use m_getcname
      use m_getiname
      use m_inichar
      use m_lspp
      use m_p2gpv
      use m_smoo2p
      use m_smoo4p
      use m_vspp

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: forcep, s_forcep

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface forcep

        module procedure s_forcep

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
      subroutine s_forcep(fpnggvar,fplspvar,fpvspvar,                   &
     &                    fplspopt,fpvspopt,fpsmtopt,ksp0,nggdmp,gtinc, &
     &                    ni,nj,nk,jcb,jcb8u,jcb8v,jcb8w,mf8u,mf8v,     &
     &                    u,v,wc,pp,ppp,rbcxy,rbct,ppgpv,pptd,pfrc,     &
     &                    tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpnggvar
                       ! Formal parameter of unique index of nggvar

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

      real, intent(in) :: nggdmp
                       ! Analysis nudging damping coefficient for GPV

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: jcb8u(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at u points

      real, intent(in) :: jcb8v(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at v points

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: mf8u(0:ni+1,0:nj+1)
                       ! Map scale factors at u points

      real, intent(in) :: mf8v(0:ni+1,0:nj+1)
                       ! Map scale factors at v points

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at present

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at present

      real, intent(in) :: wc(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity
                       ! at present

      real, intent(in) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at present

      real, intent(in) :: ppp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at past

      real, intent(in) :: rbcxy(1:ni,1:nj)
                       ! Relaxed lateral sponge damping coefficients

      real, intent(in) :: rbct(1:ni,1:nj,1:nk,1:2)
                       ! Relaxed top sponge damping coefficients

      real, intent(in) :: ppgpv(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation of GPV data
                       ! at marked time

      real, intent(in) :: pptd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! pressure perturbation of GPV data

! Input and output variable

      real, intent(inout) :: pfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in pressure equation

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

      integer lspopt   ! Option for lateral sponge damping
      integer vspopt   ! Option for vertical sponge damping

      integer smtopt   ! Option for numerical smoothing

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

      real, intent(inout) :: tmp7(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp8(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! Initialize the character variables.

      call inichar(nggvar)
      call inichar(lspvar)
      call inichar(vspvar)

! -----

! Get the required namelist variables.

      call getcname(fpnggvar,nggvar)
      call getcname(fplspvar,lspvar)
      call getcname(fpvspvar,vspvar)
      call getiname(fplspopt,lspopt)
      call getiname(fpvspopt,vspopt)
      call getiname(fpsmtopt,smtopt)

! -----

! Calculate the pressure advection.

      call advp(idadvopt,idmpopt,idmfcopt,iddiaopt,                     &
     &          idiwest,idieast,idjsouth,idjnorth,iddxiv,iddyiv,        &
     &          iddziv,ni,nj,nk,mf8u,mf8v,jcb8u,jcb8v,jcb8w,u,v,wc,     &
     &          pp,pfrc,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8)

! -----

! Calculate the 2nd order smoothing.

      if(mod(smtopt,10).eq.1) then

        call smoo2p(idsmhcoe,idsmvcoe,ni,nj,nk,ppp,pfrc)

      end if

! -----

! Calculate the 4th order smoothing.

      if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

        call smoo4p(idsmtopt,idiwest,idieast,idjsouth,idjnorth,         &
     &              idsmhcoe,idsmvcoe,ni,nj,nk,ppp,pfrc,                &
     &              tmp1,tmp2,tmp3,tmp4)

      end if

! -----

! Perform the analysis nudging to GPV.

      if(nggdmp.gt.0.e0.and.nggvar(4:4).eq.'o') then

        call p2gpv(nggdmp,gtinc,ni,nj,nk,jcb,ppp,ppgpv,pptd,pfrc)

      end if

! -----

! Calculate the letaral sponge damping.

      if(lspopt.ge.1.and.lspvar(4:4).eq.'o') then

        call lspp(idgpvvar,idlspopt,idwdnews,idlsnews,idlspsmt,gtinc,   &
     &            ni,nj,nk,jcb,ppp,rbcxy,ppgpv,pptd,pfrc,tmp1)

      end if

! -----

! Calculate the vertical sponge damping.

      if(vspopt.ge.1.and.vspvar(4:4).eq.'o') then

        call vspp(idgpvvar,idvspopt,ksp0,gtinc,ni,nj,nk,jcb,ppp,rbct,   &
     &            ppgpv,pptd,pfrc)

      end if

! -----

      end subroutine s_forcep

!-----7--------------------------------------------------------------7--

      end module m_forcep
