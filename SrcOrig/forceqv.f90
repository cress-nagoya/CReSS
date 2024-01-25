!***********************************************************************
      module m_forceqv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/08/03
!     Modification: 1999/08/09, 1999/08/18, 1999/08/23, 1999/09/30,
!                   1999/10/07, 1999/10/12, 1999/11/01, 1999/11/19,
!                   1999/11/24, 2000/01/17, 2000/02/02, 2000/02/07,
!                   2000/04/11, 2000/04/18, 2000/12/19, 2001/01/15,
!                   2001/03/13, 2001/04/15, 2001/05/29, 2001/06/06,
!                   2001/07/13, 2001/08/07, 2001/11/20, 2002/04/02,
!                   2002/06/18, 2002/07/23, 2002/08/15, 2002/09/09,
!                   2002/10/31, 2002/12/11, 2003/01/04, 2003/01/20,
!                   2003/03/13, 2003/03/21, 2003/04/30, 2003/05/19,
!                   2003/07/15, 2003/10/10, 2003/11/28, 2003/12/12,
!                   2004/02/01, 2004/03/05, 2004/04/15, 2004/05/31,
!                   2004/06/10, 2004/08/01, 2004/08/20, 2006/01/10,
!                   2006/02/13, 2006/04/03, 2006/05/12, 2006/06/21,
!                   2006/09/21, 2006/11/06, 2007/05/07, 2007/07/30,
!                   2008/05/02, 2008/08/25, 2008/12/11, 2009/02/27,
!                   2009/03/23, 2011/09/22, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the forcing term in the water vapor mixing ratio
!     equation.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_advs
      use m_comindx
      use m_getcname
      use m_getiname
      use m_inichar
      use m_lspqv
      use m_nlsmqv
      use m_s2gpv
      use m_smoo2qv
      use m_smoo4qv
      use m_turbflx
      use m_turbs
      use m_vspqv

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: forceqv, s_forceqv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface forceqv

        module procedure s_forceqv

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
      subroutine s_forceqv(fpnggvar,fplspvar,fpvspvar,                  &
     &                     fplspopt,fpvspopt,fpsmtopt,fptubopt,         &
     &                     ksp0,nggdmp,gtinc,ni,nj,nk,j31,j32,jcb,      &
     &                     jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,qvbr,rbr,rst, &
     &                     rstxu,rstxv,rstxwc,qv,qvp,rkh8u,rkh8v,rkv8w, &
     &                     rbcxy,rbct,qvgpv,qvtd,qvfrc,h3,              &
     &                     tmp1,tmp2,tmp3,tmp4,tmp5)
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

      real, intent(in) :: nggdmp
                       ! Analysis nudging damping coefficient for GPV

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

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

      real, intent(in) :: qvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state water vapor mixing ratio

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

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at present

      real, intent(in) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

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

      real, intent(in) :: qvgpv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio of GPV data
                       ! at marked time

      real, intent(in) :: qvtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! water vapor mixing ratio of GPV data

! Input and output variable

      real, intent(inout) :: qvfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in water vapor mixing ratio equation

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
      integer tubopt   ! Option for turbulent mixing

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
      call getiname(fptubopt,tubopt)

! -----

! Calculate the advection.

      call advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,             &
     &          iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc,       &
     &          qv,qvfrc,tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

      if(mod(smtopt,10).eq.1) then

       call smoo2qv(idsmhcoe,idsmvcoe,ni,nj,nk,qvbr,rbr,qvp,qvfrc,tmp1)

      end if

! -----

! Calculate the 4th order smoothing.

      if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

        call smoo4qv(idsmtopt,idiwest,idieast,idjsouth,idjnorth,        &
     &               idsmhcoe,idsmvcoe,ni,nj,nk,qvbr,rbr,qvp,qvfrc,     &
     &               tmp1,tmp2,tmp3,tmp4,tmp5)

      end if

! -----

! Calculate the non linear smoothing.

      if(smtopt.ge.11) then

       call nlsmqv(idnlhcoe,idnlvcoe,.001e0,ni,nj,nk,qvbr,rbr,qvp,qvfrc,&
     &             tmp1,tmp2,tmp3,tmp4)

      end if

! -----

! Calculate the turbulent mixing.

      if(tubopt.ge.1) then

        call turbflx(idtrnopt,idsfcopt,iddxiv,iddyiv,iddziv,ni,nj,nk,   &
     &               j31,j32,jcb,qvp,qvfrc,rkh8u,rkh8v,rkv8w,           &
     &               tmp1,tmp2,h3,tmp3,tmp4,tmp5)

        call turbs(idtrnopt,idmpopt,idmfcopt,iddxiv,iddyiv,iddziv,      &
     &             ni,nj,nk,j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,     &
     &             tmp1,tmp2,h3,qvfrc,tmp3,tmp4,tmp5)

      end if

! -----

! Perform the analysis nudging to GPV.

      if(nggdmp.gt.0.e0.and.nggvar(6:6).eq.'o') then

        call s2gpv(idgpvvar,2,nggdmp,gtinc,                             &
     &             ni,nj,nk,rst,qvp,qvgpv,qvtd,qvfrc)

      end if

! -----

! Calculate the lateral sponge damping.

      if(lspopt.ge.1.and.lspvar(6:6).eq.'o') then

        call lspqv(idgpvvar,idlspopt,idwdnews,idlsnews,idlspsmt,gtinc,  &
     &             ni,nj,nk,rst,qvbr,qvp,rbcxy,qvgpv,qvtd,qvfrc,tmp1)

      end if

! -----

! Calculate the vertical sponge damping.

      if(vspopt.ge.1.and.vspvar(6:6).eq.'o') then

        call vspqv(idgpvvar,idvspopt,ksp0,gtinc,ni,nj,nk,rst,qvbr,qvp,  &
     &             rbct,qvgpv,qvtd,qvfrc)

      end if

! -----

      end subroutine s_forceqv

!-----7--------------------------------------------------------------7--

      end module m_forceqv
