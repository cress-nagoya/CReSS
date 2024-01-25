!***********************************************************************
      module m_getllsw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/04/11
!     Modification: 2007/06/27, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2009/11/13, 2011/08/09, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the latitude and the longitude at south-west corner.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_commath
      use m_commpi
      use m_getiname
      use m_getrname
      use m_setproj

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getllsw, s_getllsw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getllsw

        module procedure s_getllsw

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic atan
      intrinsic cos
      intrinsic exp
      intrinsic log
      intrinsic max
      intrinsic min
      intrinsic real
      intrinsic sqrt
      intrinsic tan

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_getllsw(fpmpopt,fpnspol,fptlon,fpdx,fpdy,            &
     &                     ni,nj,latsw,lonsw)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpnspol
                       ! Formal parameter of unique index of nspol

      integer, intent(in) :: fptlon
                       ! Formal parameter of unique index of tlon

      integer, intent(in) :: fpdx
                       ! Formal parameter of unique index of dx

      integer, intent(in) :: fpdy
                       ! Formal parameter of unique index of dy

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

! Output variables

      real, intent(out) :: latsw
                       ! Latitude at south-west corner

      real, intent(out) :: lonsw
                       ! Longitude at south-west corner

! Internal shared variables

      integer mpopt    ! Option for map projection
      integer nspol    ! Option for projected region

      real tlon        ! True longitude

      real dx          ! Grid distance in x direction
      real dy          ! Grid distance in y direction

      real x0          ! x origin
      real y0          ! y origin
      real cpj(1:7)    ! Map projection parameters

      real xsw         ! x coordinates at south-west corner
      real ysw         ! y coordinates at south-west corner

      real xx          ! x coordinate on map coordinates system
      real yy          ! y coordinate on map coordinates system

      real rr          ! Radius on map coordinates system

      real rpol        ! real(nspol)

      real r2d2        ! 2.0 x r2d
      real r2d3        ! cpj(3) x r2d
      real r2d5        ! cpj(5) x r2d

      real tlonw       ! tlon + 180.0

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpmpopt,mpopt)
      call getiname(fpnspol,nspol)
      call getrname(fptlon,tlon)
      call getrname(fpdx,dx)
      call getrname(fpdy,dy)

! -----

! Get the map projection parameters.

      call setproj(idmpopt,idnspol,iddx,iddy,idulat,idulon,idriu,idrju, &
     &             idtlat1,idtlat2,idtlon,'solver  ',6,x0,y0,cpj)

! -----

! Calculate the x and the y coordinates at south-west corner.

      xsw=.5e0*real(2*(ni-3)*nisub*igrp+1)*dx
      ysw=.5e0*real(2*(nj-3)*njsub*jgrp+1)*dy

! -----

! Set the common used variables.

      rpol=real(nspol)

      r2d2=2.e0*r2d
      r2d3=cpj(3)*r2d
      r2d5=cpj(5)*r2d

      tlonw=tlon+180.e0

! -----

!! Calculate the latitude and the longitude at south-west corner.

! Calculate the latitude and the longitude with latitude and longitude
! coordinates.

      if(mpopt.eq.0.or.mpopt.eq.10) then

        xx=xsw+x0
        yy=ysw+y0

        latsw=max(min(yy*r2d3,90.e0),-90.e0)

        lonsw=xx*r2d3

        if(lonsw.gt.180.e0) then
          lonsw=lonsw-360.e0
        end if

        if(lonsw.lt.-180.e0) then
          lonsw=lonsw+360.e0
        end if

! -----

! Calculate the latitude and the longitude with the Polar Stereographic
! projection method.

      else if(mpopt.eq.1) then

        xx=xsw+x0
        yy=rpol*ysw+y0

        rr=sqrt(xx*xx+yy*yy)*cpj(3)

        latsw=max(min(rpol*(90.e0-r2d2*atan(rr)),90.e0),-90.e0)

        if(yy.gt.0.e0) then
          lonsw=tlonw+atan(-xx/(yy+eps))*r2d
        else
          lonsw=tlon+atan(-xx/(yy-eps))*r2d
        end if

        if(lonsw.gt.180.e0) then
          lonsw=lonsw-360.e0
        end if

        if(lonsw.lt.-180.e0) then
          lonsw=lonsw+360.e0
        end if

! -----

! Calculate the latitude and the longitude with the Lambert Conformal
! Conic projection method.

      else if(mpopt.eq.2) then

        xx=xsw+x0
        yy=rpol*ysw+y0

        rr=cpj(2)*exp(cpj(5)*log(sqrt(xx*xx+yy*yy)*cpj(7)+eps))

        latsw=max(min(rpol*(90.e0-r2d2*atan(rr)),90.e0),-90.e0)

        if(yy.gt.0.e0) then
          lonsw=tlonw+atan(-xx/(yy+eps))*r2d5
        else
          lonsw=tlon+atan(-xx/(yy-eps))*r2d5
        end if

        if(lonsw.gt.180.e0) then
          lonsw=lonsw-360.e0
        end if

        if(lonsw.lt.-180.e0) then
          lonsw=lonsw+360.e0
        end if

! -----

! Calculate the latitude and the longitude with the Mercator projection
! method.

      else if(mpopt.eq.3.or.mpopt.eq.13) then

        xx=xsw+x0
        yy=ysw+y0

        latsw                                                           &
     &    =max(min(90.e0-2.e0*atan(exp(-yy*cpj(3)))*r2d,90.e0),-90.e0)

        lonsw=xx*r2d3

        if(lonsw.gt.180.e0) then
          lonsw=lonsw-360.e0
        end if

        if(lonsw.lt.-180.e0) then
          lonsw=lonsw+360.e0
        end if

! -----

! Calculate the latitude and the longitude directly without any
! projection method.

      else if(mpopt.eq.4) then

        xx=xsw+x0
        yy=ysw+y0

        latsw=max(min(cpj(2)*yy,90.e0),-90.e0)

        lonsw=cpj(2)*xx/cos(latsw*d2r)+tlon

        if(lonsw.gt.180.e0) then
          lonsw=lonsw-360.e0
        end if

        if(lonsw.lt.-180.e0) then
          lonsw=lonsw+360.e0
        end if

! -----

! Calculate the latitude and the longitude with the circular cylinder
! coordinates system.

      else if(mpopt.eq.5) then

        xx=xsw+x0

        latsw=max(min(cpj(2)*y0,90.e0),-90.e0)

        lonsw=cpj(2)*xx/cos(latsw*d2r)

        if(lonsw.gt.180.e0) then
          lonsw=lonsw-360.e0
        end if

        if(lonsw.lt.-180.e0) then
          lonsw=lonsw+360.e0
        end if

      end if

! -----

!! -----

      end subroutine s_getllsw

!-----7--------------------------------------------------------------7--

      end module m_getllsw
