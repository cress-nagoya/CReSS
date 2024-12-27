!***********************************************************************
      module m_forceqa
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/08/18
!     Modification: 2011/09/22, 2013/02/13, 2013/03/27

!     Author      : Satoki Tsujino
!     Modification: 2024/12/26

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the forcing terms in the aerosol equations.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_advs
      use m_comindx
      use m_getcname
      use m_getiname
      use m_inichar
      use m_lsps
      use m_s2gpv
      use m_smoo2s
      use m_smoo4s
      use m_turbflx
      use m_turbs
      use m_vsps

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: forceqa, s_forceqa

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface forceqa

        module procedure s_forceqa

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
      subroutine s_forceqa(fpnggvar,fplspvar,fpvspvar,                  &
     &                     fplspopt,fpvspopt,fpsmtopt,fptubopt,         &
     &                     ksp0,nggdmp,atinc,ni,nj,nk,nqa,j31,j32,      &
     &                     jcb,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,rbr,rst,  &
     &                     rstxu,rstxv,rstxwc,qasl,qaslp,rkh8u,rkh8v,   &
     &                     rkv8w,rbcxy,rbct,qagpv,qatd,qafrc,           &
     &                     h3,tmp1,tmp2,tmp3,tmp4,tmp5)
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

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

      real, intent(in) :: nggdmp
                       ! Analysis nudging damping coefficient for GPV

      real, intent(in) :: atinc
                       ! Lapse of forecast time
                       ! from aerosol data reading

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

      real, intent(in) :: qasl(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at present

      real, intent(in) :: qaslp(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at past

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

      real, intent(in) :: qagpv(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Mixing ratio of aerosol data at marked time

      real, intent(in) :: qatd(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Time tendency of mixing ratio of aerosol data

! Output variable

      real, intent(out) :: qafrc(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Forcing terms in aerosol equations

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

      integer n        ! Array index in aerosol categories

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

! -----

! Get the required namelist variables.

      call getcname(fpnggvar,nggvar)
      call getcname(fplspvar,lspvar)
      call getcname(fpvspvar,vspvar)
      call getiname(fplspopt,lspopt)
      call getiname(fpvspopt,vspopt)
      call getiname(fpsmtopt,smtopt)
      call getiname(fptubopt,tubopt)

      call getcname(iddmpvar,dmpvar)

! -----

!! Repeat processes for all types of aerosol.

      do n=1,nqa(0)

! Calculate the advection.

        call s_advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,         &
     &              iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc,   &
     &              qasl(0,0,1,n),qafrc(0,0,1,n),tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

        if(mod(smtopt,10).eq.1) then

          call s_smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,rbr,qaslp(0,0,1,n),  &
     &                  qafrc(0,0,1,n),tmp1)

        end if

! -----

! Calculate the 4th order smoothing.

        if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

          call s_smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth,     &
     &                  idsmhcoe,idsmvcoe,ni,nj,nk,rbr,qaslp(0,0,1,n),  &
     &                  qafrc(0,0,1,n),tmp1,tmp2,tmp3,tmp4,tmp5)

        end if

! -----

! Calculate the turbulent mixing.

        if(tubopt.ge.1) then

          call s_turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,        &
     &                   ni,nj,nk,j31,j32,jcb,qaslp(0,0,1,n),           &
     &                   qafrc(0,0,1,n),rkh8u,rkh8v,rkv8w,              &
     &                   tmp1,tmp2,h3,tmp3,tmp4,tmp5)

          call s_turbs(idtrnopt,idmpopt,idmfcopt,iddxiv,iddyiv,iddziv,  &
     &                 iddmpvar,ni,nj,nk,'xx',                          &
     &                 j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,          &
     &                 tmp1,tmp2,h3,qafrc(0,0,1,n),tmp3,tmp4,tmp5)

        end if

! -----

! Perform the analysis nudging to aerosol data.

        if(nggdmp.gt.0.e0.and.nggvar(8:8).eq.'o') then

          call s_s2gpv(idgpvvar,9,nggdmp,atinc,                         &
     &                 ni,nj,nk,rst,qaslp(0,0,1,n),                     &
     &                 qagpv(0,0,1,n),qatd(0,0,1,n),qafrc(0,0,1,n))

        end if

! -----

! Calculate the lateral sponge damping.

        if(lspopt.ge.1.and.lspvar(8:8).eq.'o') then

          call s_lsps(idgpvvar,idlspopt,idwdnews,idlsnews,idlspsmt,     &
     &                9,atinc,ni,nj,nk,rst,qaslp(0,0,1,n),rbcxy,        &
     &                qagpv(0,0,1,n),qatd(0,0,1,n),qafrc(0,0,1,n),tmp1)

        end if

! -----

! Calculate the vertical sponge damping.

        if(vspopt.ge.1.and.vspvar(8:8).eq.'o') then

          call s_vsps(idgpvvar,idvspopt,9,ksp0,                         &
     &                atinc,ni,nj,nk,rst,qaslp(0,0,1,n),rbct,           &
     &                qagpv(0,0,1,n),qatd(0,0,1,n),qafrc(0,0,1,n))

        end if

! -----

      end do

!! -----

      end subroutine s_forceqa

!-----7--------------------------------------------------------------7--

      end module m_forceqa
