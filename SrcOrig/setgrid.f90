!***********************************************************************
      module m_setgrid
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/01/25, 1999/03/25, 1999/04/06,
!                   1999/05/10, 1999/05/20, 1999/07/05, 1999/08/18,
!                   1999/08/23, 1999/09/01, 1999/09/30, 1999/11/01,
!                   1999/11/30, 1999/12/06, 2000/01/17, 2000/04/18,
!                   2000/12/18, 2001/03/13, 2001/04/15, 2001/05/29,
!                   2001/06/06, 2001/10/18, 2002/04/02, 2002/07/03,
!                   2002/07/15, 2002/10/31, 2003/01/04, 2003/03/21,
!                   2003/05/19, 2003/09/01, 2003/11/05, 2004/06/10,
!                   2004/08/01, 2004/08/20, 2005/02/10, 2005/04/04,
!                   2006/01/10, 2006/04/03, 2006/09/30, 2006/11/06,
!                   2006/12/04, 2007/01/05, 2007/01/20, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2011/08/09, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the model grid.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bc2d
      use m_bcycle
      use m_bcyclex
      use m_combuf
      use m_comindx
      use m_getbufgx
      use m_getbufgy
      use m_getbufsx
      use m_getbufsy
      use m_getiname
      use m_gettrn
      use m_getxy
      use m_getz
      use m_jacobian
      use m_mapfct
      use m_phycood
      use m_putbufgx
      use m_putbufgy
      use m_putbufsx
      use m_putbufsy
      use m_setcst2d
      use m_setcst3d
      use m_setproj
      use m_shiftgx
      use m_shiftgy
      use m_shiftsx
      use m_shiftsy
      use m_trilat
      use m_xy2ll

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: setgrid, s_setgrid

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface setgrid

        module procedure s_setgrid

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
      subroutine s_setgrid(fpwbc,fpebc,fptrnopt,fpexbopt,               &
     &                     fpmfcopt,fpcoropt,area,ni,nj,nk,             &
     &                     zph,zsth,lat,lon,j31,j32,jcb,jcb8u,jcb8v,    &
     &                     jcb8w,mf,mf8u,mf8v,rmf,rmf8u,rmf8v,fc,x,y,z)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fptrnopt
                       ! Formal parameter of unique index of trnopt

      integer, intent(in) :: fpexbopt
                       ! Formal parameter of unique index of exbopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: fpcoropt
                       ! Formal parameter of unique index of coropt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Output variables

      real, intent(out) :: area(0:4)
                       ! Area of each boundary plane

      real, intent(out) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(out) :: zsth(1:nk)
                       ! 1 dimensional stretched z coordinates

      real, intent(out) :: lat(0:ni+1,0:nj+1)
                       ! Latitude

      real, intent(out) :: lon(0:ni+1,0:nj+1)
                       ! Longitude

      real, intent(out) :: j31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of Jacobian

      real, intent(out) :: j32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of Jacobian

      real, intent(out) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(out) :: jcb8u(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at u points

      real, intent(out) :: jcb8v(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at v points

      real, intent(out) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(out) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(out) :: mf8u(0:ni+1,0:nj+1)
                       ! Map scale factors at u points

      real, intent(out) :: mf8v(0:ni+1,0:nj+1)
                       ! Map scale factors at v points

      real, intent(out) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(out) :: rmf8u(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at u points

      real, intent(out) :: rmf8v(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at v points

      real, intent(out) :: fc(0:ni+1,0:nj+1,1:2)
                       ! 0.25 x Coriolis parameters

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions

      integer trnopt   ! Option for terrain height setting
      integer exbopt   ! Option for external boundary forcing
      integer mfcopt   ! Option for map scale factor
      integer coropt   ! Option for Coriolis force

      integer ni_sub   ! Substitute for ni
      integer nj_sub   ! Substitute for nj

      real x0          ! x origin
      real y0          ! y origin
      real cpj(1:7)    ! Map projection parameters

      real, intent(inout) :: x(0:ni+1)
                       ! x coordinates

      real, intent(inout) :: y(0:nj+1)
                       ! y coordinates

      real, intent(inout) :: z(1:nk)
                       ! zeta coordinates

! Remark

!     jcb,jcb8w: These variables are also temporary.

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fptrnopt,trnopt)
      call getiname(fpexbopt,exbopt)
      call getiname(fpmfcopt,mfcopt)
      call getiname(fpcoropt,coropt)

! -----

! Set the substituted variables.

      ni_sub=ni
      nj_sub=nj

! -----

! Get the x and the y coordinets at the model grid points.

      call getxy(iddx,iddy,'xx',0,ni+1,0,nj+1,x,y)

! -----

! Get the map projection parameters and the latitude and longitude.

      call setproj(idmpopt,idnspol,iddx,iddy,idulat,idulon,idriu,idrju, &
     &             idtlat1,idtlat2,idtlon,'solver  ',6,x0,y0,cpj)

      call xy2ll(idmpopt,idnspol,idtlon,'solver  ',6,'xx',x0,y0,cpj,    &
     &           0,ni+1,0,nj+1,x,y,lat,lon)

! -----

! Fill in the array with 0 in the case the map scale factors are not
! applied.

      if(mfcopt.eq.0) then

        call setcst2d(0,ni+1,0,nj+1,1.e0,mf)
        call setcst2d(0,ni+1,0,nj+1,1.e0,mf8u)
        call setcst2d(0,ni+1,0,nj+1,1.e0,mf8v)

        call setcst3d(0,ni+1,0,nj+1,1,3,1.e0,rmf)
        call setcst3d(0,ni+1,0,nj+1,1,3,1.e0,rmf8u)
        call setcst3d(0,ni+1,0,nj+1,1,3,1.e0,rmf8v)

! -----

! Calculate the map scale factors.

      else if(mfcopt.eq.1) then

        call s_mapfct(idmpopt,idnspol,idadvopt,                         &
     &                idtubopt,iddisr,iddxiv,iddyiv,cpj,                &
     &                ni,nj,x,lat,mf,mf8u,mf8v,rmf,rmf8u,rmf8v,jcb)

      end if

! -----

! Fill in the array fc with 0 in the case the Coriolis force is not
! applied.

      if(coropt.eq.0) then

        call setcst3d(0,ni+1,0,nj+1,1,2,0.e0,fc)

! -----

! Get the Coriolis parameters.

      else if(coropt.ge.1) then

        call trilat(idcoropt,ni,nj,lat,fc)

      end if

! -----

! Get the terrain height.

      if(trnopt.eq.2.and.mod(exbopt,10).eq.2) then

        call s_gettrn(idtrnopt,idzsfc,idmnthgh,idmntwx,idmntwy,         &
     &                idmntcx,idmntcy,'terrain.damp',12,0,ni,nj,x,y,jcb)

      else

        call s_gettrn(idtrnopt,idzsfc,idmnthgh,idmntwx,idmntwy,         &
     &                idmntcx,idmntcy,'terrain     ',7,0,ni,nj,x,y,jcb)

      end if

! -----

!! Set the lateral boundary conditions.

      if(exbopt.eq.0) then

! Exchange the value between first halo regions.

        call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,jcb,1,1,sbuf)

        call s_shiftsx(idwbc,idebc,'bnd',nj,1,1,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,jcb,1,1,rbuf)

        call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,1,jcb,1,1,sbuf)

        call s_shiftsy(idsbc,idnbc,'bnd',ni,1,1,sbuf,rbuf)

        call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,1,jcb,1,1,rbuf)

        call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,jcb,1,1,sbuf)

        call s_shiftgx(idwbc,idebc,'bnd',nj,1,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,jcb,1,1,rbuf)

        call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,1,jcb,1,1,sbuf)

        call s_shiftgy(idsbc,idnbc,'bnd',ni,1,1,sbuf,rbuf)

        call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,1,jcb,1,1,rbuf)

        call bcycle(idwbc,idebc,idsbc,idnbc,                            &
     &              2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,1,jcb)

! -----

! Exchange the value between second halo regions.

        call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,1,jcb,1,1,sbuf)

        call s_shiftsx(idwbc,idebc,'all',nj,1,1,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,1,jcb,1,1,rbuf)

        call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,1,jcb,1,1,sbuf)

        call s_shiftsy(idsbc,idnbc,'all',ni,1,1,sbuf,rbuf)

        call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,1,jcb,1,1,rbuf)

        call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,1,jcb,1,1,sbuf)

        call s_shiftgx(idwbc,idebc,'all',nj,1,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,1,jcb,1,1,rbuf)

        call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,1,jcb,1,1,sbuf)

        call s_shiftgy(idsbc,idnbc,'all',ni,1,1,sbuf,rbuf)

        call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,1,jcb,1,1,rbuf)

        call bcycle(idwbc,idebc,idsbc,idnbc,                            &
     &              3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,ni,nj,1,jcb)

! -----

! Set the boundary conditions.

        call s_bc2d(idwbc,idebc,idsbc,idnbc,ni,nj,jcb)

! ------

      end if

!! -----

! Set the periodic boundary conditions.

      if(exbopt.ge.1.and.abs(wbc).eq.1.and.abs(ebc).eq.1) then

        call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,jcb,1,1,sbuf)

        call s_shiftsx(idwbc,idebc,'bnd',nj,1,1,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,jcb,1,1,rbuf)

        call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,jcb,1,1,sbuf)

        call s_shiftgx(idwbc,idebc,'bnd',nj,1,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,jcb,1,1,rbuf)

        call s_bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,1,jcb)

        call s_putbufsx(idwbc,idebc,'bnd',3,ni-3,ni,nj,1,jcb,1,1,sbuf)

        call s_shiftsx(idwbc,idebc,'bnd',nj,1,1,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,1,jcb,1,1,rbuf)

        call s_putbufgx(idwbc,idebc,'bnd',3,ni-3,ni,nj,1,jcb,1,1,sbuf)

        call s_shiftgx(idwbc,idebc,'bnd',nj,1,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,1,jcb,1,1,rbuf)

        call bcyclex(idwbc,idebc,3,0,ni-3,ni_sub,ni,nj,1,jcb)

      end if

! -----

! Get the zeta coordinates at the model grid points.

      call getz(iddz,idzsfc,nk,z)

! -----

! Calculate the z physical coordinates.

      call s_phycood(idsthopt,idzsfc,idzflat,ni,nj,nk,z,jcb,zph,zsth,   &
     &               jcb8w)

! -----

! Get the x and the y coordinates at the model grid points.

      call getxy(iddx,iddy,'oo',0,ni+1,0,nj+1,x,y)

! -----

! Calculate the Jacobian.

      call jacobian(idwbc,idebc,idexbopt,idadvopt,idsmtopt,idtubopt,    &
     &              area,ni,nj,nk,x,y,z,zph,rmf,rmf8u,rmf8v,j31,j32,    &
     &              jcb,jcb8u,jcb8v,jcb8w)

! -----

      end subroutine s_setgrid

!-----7--------------------------------------------------------------7--

      end module m_setgrid
