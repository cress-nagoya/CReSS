!***********************************************************************
      module m_xy2ll
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/03/17, 1999/03/25, 1999/04/06, 1999/05/10,
!                   1999/06/07, 1999/09/30, 1999/10/12, 1999/11/01,
!                   1999/11/19, 2000/01/17, 2000/07/05, 2001/01/09,
!                   2001/04/15, 2001/05/29, 2001/08/07, 2001/11/20,
!                   2002/04/02, 2002/11/20, 2003/04/30, 2003/05/19,
!                   2003/11/28, 2004/05/31, 2004/07/10, 2004/09/01,
!                   2005/02/10, 2006/09/21, 2006/11/06, 2007/06/27,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2008/10/10,
!                   2009/02/27, 2009/11/13, 2011/08/09, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the latitude and the longitude from the x and the y
!     coordinates.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_getindx
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: xy2ll, s_xy2ll

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface xy2ll

        module procedure s_xy2ll

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
      subroutine s_xy2ll(fpmpopt,fpnspol,fptlon,pname,ncpn,xo,          &
     &                   x0,y0,cpj,imin,imax,jmin,jmax,x,y,lat,lon)
!***********************************************************************

! Input variables

      character(len=8), intent(in) :: pname
                       ! Running program name

      character(len=2), intent(in) :: xo
                       ! Control flag of variable arrangement

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpnspol
                       ! Formal parameter of unique index of nspol

      integer, intent(in) :: fptlon
                       ! Formal parameter of unique index of tlon

      integer, intent(in) :: ncpn
                       ! Number of character of pname

      integer, intent(in) :: imin
                       ! Minimum array index in x direction

      integer, intent(in) :: imax
                       ! Maximum array index in x direction

      integer, intent(in) :: jmin
                       ! Minimum array index in y direction

      integer, intent(in) :: jmax
                       ! Maximum array index in y direction

      real, intent(in) :: x0
                       ! x origin

      real, intent(in) :: y0
                       ! y origin

      real, intent(in) :: cpj(1:7)
                       ! Map projection parameters

      real, intent(in) :: x(imin:imax)
                       ! x coordinates

      real, intent(in) :: y(jmin:jmax)
                       ! y coordinates

! Output variables

      real, intent(out) :: lat(imin:imax,jmin:jmax)
                       ! Latitude

      real, intent(out) :: lon(imin:imax,jmin:jmax)
                       ! Longitude

! Internal shared variables

      integer mpopt    ! Option for map projection
      integer nspol    ! Option for projected region

      integer istr     ! Minimum do loops index in x direction
      integer iend     ! Maximum do loops index in x direction
      integer jstr     ! Minimum do loops index in y direction
      integer jend     ! Maximum do loops index in y direction

      real tlon        ! True longitude

      real rpol        ! real(nspol)

      real r2d2        ! 2.0 x r2d
      real r2d3        ! cpj(3) x r2d
      real r2d5        ! cpj(5) x r2d

      real tlonw       ! tlon + 180.0

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      real xx          ! x coordinate on map coordinates system
      real yy          ! y coordinate on map coordinates system

      real rr          ! Radius on map coordinates system

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpmpopt,mpopt)
      call getiname(fpnspol,nspol)
      call getrname(fptlon,tlon)

! -----

! Get the maximum and minimim indices of do loops.

      call getindx(xo,imin,imax,jmin,jmax,istr,iend,jstr,jend)

! -----

! Set the common used variables.

      rpol=real(nspol)

      r2d2=2.e0*r2d
      r2d3=cpj(3)*r2d
      r2d5=cpj(5)*r2d

      tlonw=tlon+180.e0

! -----

!!! Calculate the latitude and the longitude from the x and the y
!!! coordinates.

!$omp parallel default(shared)

!! Calculate the latitude and the longitude with latitude and longitude
!! coordinates.

      if(mpopt.eq.0.or.mpopt.eq.10) then

! For the program, solver.

        if(pname(1:ncpn).eq.'solver') then

!$omp do schedule(runtime) private(i,j,xx,yy)

          do j=jstr,jend
          do i=istr,iend
            xx=x(i)+x0
            yy=y(j)+y0

            lat(i,j)=max(min(yy*r2d3,90.e0),-90.e0)

            lon(i,j)=xx*r2d3

            if(lon(i,j).gt.180.e0) then
              lon(i,j)=lon(i,j)-360.e0
            end if

            if(lon(i,j).lt.-180.e0) then
              lon(i,j)=lon(i,j)+360.e0
            end if

          end do
          end do

!$omp end do

! -----

! For the pre and post processors.

        else

!$omp do schedule(runtime) private(i,j)

          do j=jstr,jend
          do i=istr,iend
            lat(i,j)=max(min(y(j)+y0,90.e0),-90.e0)

            lon(i,j)=x(i)+x0

            if(lon(i,j).gt.180.e0) then
              lon(i,j)=lon(i,j)-360.e0
            end if

            if(lon(i,j).lt.-180.e0) then
              lon(i,j)=lon(i,j)+360.e0
            end if

          end do
          end do

!$omp end do

        end if

! -----

!! -----

! Calculate the latitude and the longitude with the Polar Stereographic
! projection method.

      else if(mpopt.eq.1) then

!$omp do schedule(runtime) private(i,j,xx,yy,rr)

        do j=jstr,jend
        do i=istr,iend
          xx=x(i)+x0
          yy=rpol*y(j)+y0

          rr=sqrt(xx*xx+yy*yy)*cpj(3)

          lat(i,j)=max(min(rpol*(90.e0-r2d2*atan(rr)),90.e0),-90.e0)

          if(yy.gt.0.e0) then
            lon(i,j)=tlonw+atan(-xx/(yy+eps))*r2d
          else
            lon(i,j)=tlon+atan(-xx/(yy-eps))*r2d
          end if

          if(lon(i,j).gt.180.e0) then
            lon(i,j)=lon(i,j)-360.e0
          end if

          if(lon(i,j).lt.-180.e0) then
            lon(i,j)=lon(i,j)+360.e0
          end if

        end do
        end do

!$omp end do

! -----

! Calculate the latitude and the longitude with the Lambert Conformal
! Conic projection method.

      else if(mpopt.eq.2) then

!$omp do schedule(runtime) private(i,j,xx,yy,rr)

        do j=jstr,jend
        do i=istr,iend
          xx=x(i)+x0
          yy=rpol*y(j)+y0

          rr=cpj(2)*exp(cpj(5)*log(sqrt(xx*xx+yy*yy)*cpj(7)+eps))

          lat(i,j)=max(min(rpol*(90.e0-r2d2*atan(rr)),90.e0),-90.e0)

          if(yy.gt.0.e0) then
            lon(i,j)=tlonw+atan(-xx/(yy+eps))*r2d5
          else
            lon(i,j)=tlon+atan(-xx/(yy-eps))*r2d5
          end if

          if(lon(i,j).gt.180.e0) then
            lon(i,j)=lon(i,j)-360.e0
          end if

          if(lon(i,j).lt.-180.e0) then
            lon(i,j)=lon(i,j)+360.e0
          end if

        end do
        end do

!$omp end do

! -----

! Calculate the latitude and the longitude with the Mercator projection
! method.

      else if(mpopt.eq.3.or.mpopt.eq.13) then

!$omp do schedule(runtime) private(i,j,xx,yy)

        do j=jstr,jend
        do i=istr,iend
          xx=x(i)+x0
          yy=y(j)+y0

          lat(i,j)                                                      &
     &      =max(min(90.e0-2.e0*atan(exp(-yy*cpj(3)))*r2d,90.e0),-90.e0)

          lon(i,j)=xx*r2d3

          if(lon(i,j).gt.180.e0) then
            lon(i,j)=lon(i,j)-360.e0
          end if

          if(lon(i,j).lt.-180.e0) then
            lon(i,j)=lon(i,j)+360.e0
          end if

        end do
        end do

!$omp end do

! -----

! Calculate the latitude and the longitude directly without any
! projection method.

      else if(mpopt.eq.4) then

!$omp do schedule(runtime) private(i,j,xx,yy)

        do j=jstr,jend
        do i=istr,iend
          xx=x(i)+x0
          yy=y(j)+y0

          lat(i,j)=max(min(cpj(2)*yy,90.e0),-90.e0)

          lon(i,j)=cpj(2)*xx/cos(lat(i,j)*d2r)+tlon

          if(lon(i,j).gt.180.e0) then
            lon(i,j)=lon(i,j)-360.e0
          end if

          if(lon(i,j).lt.-180.e0) then
            lon(i,j)=lon(i,j)+360.e0
          end if

        end do
        end do

!$omp end do

! -----

! Calculate the latitude and the longitude with the circular cylinder
! coordinates system.

      else if(mpopt.eq.5) then

!$omp do schedule(runtime) private(i,j,xx)

        do j=jstr,jend
        do i=istr,iend
          xx=x(i)+x0

          lat(i,j)=max(min(cpj(2)*y0,90.e0),-90.e0)

          lon(i,j)=cpj(2)*xx/cos(lat(i,j)*d2r)

          if(lon(i,j).gt.180.e0) then
            lon(i,j)=lon(i,j)-360.e0
          end if

          if(lon(i,j).lt.-180.e0) then
            lon(i,j)=lon(i,j)+360.e0
          end if

        end do
        end do

!$omp end do

      end if

! -----

!$omp end parallel

!!! -----

      end subroutine s_xy2ll

!-----7--------------------------------------------------------------7--

      end module m_xy2ll
