!***********************************************************************
      module m_forceqt
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/05/31
!     Modification: 2004/06/10, 2004/08/01, 2004/08/20, 2006/01/10,
!                   2006/02/13, 2006/04/03, 2006/05/12, 2006/06/21,
!                   2006/09/21, 2006/11/06, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2011/05/16, 2011/08/18, 2011/09/22,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the forcing term in the tracer equation.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_advs
      use m_comindx
      use m_comkind
      use m_emitqt
      use m_getcname
      use m_getiname
      use m_getrname
      use m_inichar
      use m_lsps0
      use m_smoo2s
      use m_smoo4s
      use m_turbflx
      use m_turbs
      use m_vsps0

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: forceqt, s_forceqt

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface forceqt

        module procedure s_forceqt

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic int
      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_forceqt(fplspvar,fpvspvar,fplspopt,fpvspopt,fpsmtopt,&
     &                     fptrkopt,fptubopt,fpqtstr,fpqtend,ksp0,ctime,&
     &                     ni,nj,nk,zph,j31,j32,jcb,jcb8u,jcb8v,mf,rmf, &
     &                     rmf8u,rmf8v,rbr,rst,rstxu,rstxv,rstxwc,      &
     &                     qt,qtp,rkh8u,rkh8v,rkv8w,rbcxy,rbct,         &
     &                     qtfrc,h3,tmp1,tmp2,tmp3,tmp4,tmp5)
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

      integer, intent(in) :: fptrkopt
                       ! Formal parameter of unique index of trkopt

      integer, intent(in) :: fptubopt
                       ! Formal parameter of unique index of tubopt

      integer, intent(in) :: fpqtstr
                       ! Formal parameter of unique index of qtstr

      integer, intent(in) :: fpqtend
                       ! Formal parameter of unique index of qtend

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

      real, intent(in) :: qt(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at present

      real, intent(in) :: qtp(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at past

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

! Output variable

      real, intent(out) :: qtfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in tracer equation

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
      integer trkopt   ! Option for mixing ratio tracking
      integer tubopt   ! Option for turbulent mixing

      real qtstr       ! User specified emitting start time
      real qtend       ! User specified emitting end time

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

      call inichar(lspvar)
      call inichar(vspvar)

! -----

! Get the required namelist variables.

      call getcname(fplspvar,lspvar)
      call getcname(fpvspvar,vspvar)
      call getiname(fplspopt,lspopt)
      call getiname(fpvspopt,vspopt)
      call getiname(fpsmtopt,smtopt)
      call getiname(fptrkopt,trkopt)
      call getiname(fptubopt,tubopt)
      call getrname(fpqtstr,qtstr)
      call getrname(fpqtend,qtend)

! -----

! Calculate the advection.

      call advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,             &
     &          iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc,       &
     &          qt,qtfrc,tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

      if(mod(smtopt,10).eq.1) then

        call smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,rbr,qtp,qtfrc,tmp1)

      end if

! -----

! Calculate the 4th order smoothing.

      if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

        call smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth,         &
     &              idsmhcoe,idsmvcoe,ni,nj,nk,rbr,qtp,qtfrc,           &
     &              tmp1,tmp2,tmp3,tmp4,tmp5)

      end if

! -----

! Calculate the turbulent mixing.

      if(tubopt.ge.1) then

        call turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,ni,nj,nk,   &
     &               j31,j32,jcb,qtp,qtfrc,rkh8u,rkh8v,rkv8w,           &
     &               tmp1,tmp2,h3,tmp3,tmp4,tmp5)

        call turbs(idtrnopt,idmpopt,idmfcopt,iddxiv,iddyiv,iddziv,      &
     &             ni,nj,nk,j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,     &
     &             tmp1,tmp2,h3,qtfrc,tmp3,tmp4,tmp5)

      end if

! -----

! Calculate the lateral sponge damping.

      if(lspopt.ge.1.and.lspvar(9:9).eq.'o') then

        call lsps0(idlspopt,idwdnews,idlsnews,idlspsmt,ni,nj,nk,rst,    &
     &             qtp,rbcxy,qtfrc,tmp1)

      end if

! -----

! Calculate the vertical sponge damping.

      if(vspopt.ge.1.and.vspvar(9:9).eq.'o') then

        call vsps0(ksp0,ni,nj,nk,rst,qtp,rbct,qtfrc)

      end if

! -----

! Emit the tracer mixing ratio from user specified location.

      if(trkopt.eq.2) then

        if(ctime.ge.1000_i8*int(qtstr+.1e0,i8).and.                     &
     &     ctime.le.1000_i8*int(qtend+.1e0,i8)) then

          call s_emitqt(idqt0opt,idqt0num,idqt0rx,idqt0ry,idqt0rz,      &
     &                  idqt0cx,idqt0cy,idqt0cz,idqt0ds,idqtdt,         &
     &                  ni,nj,nk,zph,qtfrc,tmp1,tmp2)

        end if

      end if

! -----

      end subroutine s_forceqt

!-----7--------------------------------------------------------------7--

      end module m_forceqt
