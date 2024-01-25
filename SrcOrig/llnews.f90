!***********************************************************************
      module m_llnews
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/04/11
!     Modification: 2007/06/27, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2009/11/13, 2011/08/09, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the latitude and the longitude at model corner.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: llnews, s_llnews

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface llnews

        module procedure s_llnews

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
      subroutine s_llnews(fpmpopt,fpnspol,fptlon,x0,y0,cpj,             &
     &                    x4,y4,lat4,lon4)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpnspol
                       ! Formal parameter of unique index of nspol

      integer, intent(in) :: fptlon
                       ! Formal parameter of unique index of tlon

      real, intent(in) :: x0
                       ! x origin

      real, intent(in) :: y0
                       ! y origin

      real, intent(in) :: cpj(1:7)
                       ! Map projection parameters

      real, intent(in) :: x4
                       ! x coordinates at model corners

      real, intent(in) :: y4
                       ! y coordinates at model corners

! Output variables

      real, intent(out) :: lat4
                       ! Latitude at model corners

      real, intent(out) :: lon4
                       ! Longitude at model corners

! Internal shared variables

      integer mpopt    ! Option for map projection
      integer nspol    ! Option for projected region

      real tlon        ! True longitude

      real rr          ! Radius on map coordinates system

      real xx          ! x coordinate on map coordinates system
      real yy          ! y coordinate on map coordinates system

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

! -----

! Set the common used variables.

      rpol=real(nspol)

      r2d2=2.e0*r2d
      r2d3=cpj(3)*r2d
      r2d5=cpj(5)*r2d

      tlonw=tlon+180.e0

! -----

!! Calculate the latitude and the longitude at model corner.

! Calculate the latitude and the longitude with latitude and longitude
! coordinates.

      if(mpopt.eq.0.or.mpopt.eq.10) then

        xx=x4+x0
        yy=y4+y0

        lat4=max(min(yy*r2d3,90.e0),-90.e0)

        lon4=xx*r2d3

        if(lon4.gt.180.e0) then
          lon4=lon4-360.e0
        end if

        if(lon4.lt.-180.e0) then
          lon4=lon4+360.e0
        end if

! -----

! Calculate the latitude and the longitude with the Polar Stereographic
! projection method.

      else if(mpopt.eq.1) then

        xx=x4+x0
        yy=rpol*y4+y0

        rr=sqrt(xx*xx+yy*yy)*cpj(3)

        lat4=max(min(rpol*(90.e0-r2d2*atan(rr)),90.e0),-90.e0)

        if(yy.gt.0.e0) then
          lon4=tlonw+atan(-xx/(yy+eps))*r2d
        else
          lon4=tlon+atan(-xx/(yy-eps))*r2d
        end if

        if(lon4.gt.180.e0) then
          lon4=lon4-360.e0
        end if

        if(lon4.lt.-180.e0) then
          lon4=lon4+360.e0
        end if

! -----

! Calculate the latitude and the longitude with the Lambert Conformal
! Conic projection method.

      else if(mpopt.eq.2) then

        xx=x4+x0
        yy=rpol*y4+y0

        rr=cpj(2)*exp(cpj(5)*log(sqrt(xx*xx+yy*yy)*cpj(7)+eps))

        lat4=max(min(rpol*(90.e0-r2d2*atan(rr)),90.e0),-90.e0)

        if(yy.gt.0.e0) then
          lon4=tlonw+atan(-xx/(yy+eps))*r2d5
        else
          lon4=tlon+atan(-xx/(yy-eps))*r2d5
        end if

        if(lon4.gt.180.e0) then
          lon4=lon4-360.e0
        end if

        if(lon4.lt.-180.e0) then
          lon4=lon4+360.e0
        end if

! -----

! Calculate the latitude and the longitude with the Mercator projection
! method.

      else if(mpopt.eq.3.or.mpopt.eq.13) then

        xx=x4+x0
        yy=y4+y0

        lat4=max(min(90.e0-2.e0*atan(exp(-yy*cpj(3)))*r2d,90.e0),-90.e0)

        lon4=xx*r2d3

        if(lon4.gt.180.e0) then
          lon4=lon4-360.e0
        end if

        if(lon4.lt.-180.e0) then
          lon4=lon4+360.e0
        end if

! -----

! Calculate the latitude and the longitude directly without any
! projection method.

      else if(mpopt.eq.4) then

        xx=x4+x0
        yy=y4+y0

        lat4=max(min(cpj(2)*yy,90.e0),-90.e0)

        lon4=cpj(2)*xx/cos(lat4*d2r)+tlon

        if(lon4.gt.180.e0) then
          lon4=lon4-360.e0
        end if

        if(lon4.lt.-180.e0) then
          lon4=lon4+360.e0
        end if

! -----

! Calculate the latitude and the longitude with the circular cylinder
! coordinates system.

      else if(mpopt.eq.5) then

        xx=x4+x0

        lat4=max(min(cpj(2)*y0,90.e0),-90.e0)

        lon4=cpj(2)*xx/cos(lat4*d2r)

        if(lon4.gt.180.e0) then
          lon4=lon4-360.e0
        end if

        if(lon4.lt.-180.e0) then
          lon4=lon4+360.e0
        end if

      end if

! -----

!! -----

      end subroutine s_llnews

!-----7--------------------------------------------------------------7--

      end module m_llnews
