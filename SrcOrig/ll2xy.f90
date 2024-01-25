!***********************************************************************
      module m_ll2xy
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/03/17, 1999/03/25, 1999/04/06, 1999/05/10,
!                   1999/06/07, 1999/07/05, 1999/09/30, 1999/10/12,
!                   1999/11/01, 2000/01/17, 2000/07/05, 2001/01/09,
!                   2001/02/13, 2001/05/29, 2001/06/29, 2001/08/07,
!                   2001/11/20, 2002/04/02, 2002/07/03, 2003/04/30,
!                   2003/05/19, 2003/11/28, 2004/05/31, 2004/09/01,
!                   2005/02/10, 2006/09/21, 2006/11/06, 2007/01/31,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2008/10/10,
!                   2009/02/27, 2011/08/09, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the x and the y coordinates from the latitude and the
!     longitude.

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

      public :: ll2xy, s_ll2xy

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface ll2xy

        module procedure s_ll2xy

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic cos
      intrinsic sin
      intrinsic tan
      intrinsic exp
      intrinsic log
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_ll2xy(fpmpopt,fpnspol,fptlon,pname,ncpn,xo,          &
     &                   x0,y0,cpj,imin,imax,jmin,jmax,lat,lon,x2d,y2d)
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

      real, intent(in) :: lat(imin:imax,jmin:jmax)
                       ! Latitude

      real, intent(in) :: lon(imin:imax,jmin:jmax)
                       ! Longitude

! Output variables

      real, intent(out) :: x2d(imin:imax,jmin:jmax)
                       ! 2 dimensional x coordinates

      real, intent(out) :: y2d(imin:imax,jmin:jmax)
                       ! 2 dimensional y coordinates

! Internal shared variables

      integer mpopt    ! Option for map projection
      integer nspol    ! Option for projected region

      integer istr     ! Minimum do loops index in x direction
      integer iend     ! Maximum do loops index in x direction
      integer jstr     ! Minimum do loops index in y direction
      integer jend     ! Maximum do loops index in y direction

      real tlon        ! True longitude

      real rpol        ! real(nspol)
      real rpol2       ! real(nspol) x cpj(2)

      real pol05       ! 0.5 x rpol
      real pold2r      ! rpol x d2r

      real d2r2        ! cpj(2) x d2r
      real d2r4        ! cpj(4) x d2r

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      real dlon        ! Distance of longitude

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

      rpol2=rpol*cpj(2)

      pol05=.5e0*rpol

      pold2r=rpol*d2r

      d2r2=cpj(2)*d2r
      d2r4=cpj(4)*d2r

! -----

!!! Calculate the x and the y coordinates from the latitude and the
!!! longitude.

!$omp parallel default(shared)

!! Calculate the reference latitude and longitude with the latitude and
!! longitude coordinates.

      if(mpopt.eq.0.or.mpopt.eq.10) then

! For the program, solver.

        if(pname(1:ncpn).eq.'solver') then

!$omp do schedule(runtime) private(i,j)

          do j=jstr,jend
          do i=istr,iend
            x2d(i,j)=d2r2*lon(i,j)-x0
            y2d(i,j)=d2r2*lat(i,j)-y0
          end do
          end do

!$omp end do

! -----

! For the pre and post processors.

        else

!$omp do schedule(runtime) private(i,j)

          do j=jstr,jend
          do i=istr,iend
            x2d(i,j)=lon(i,j)-x0
            y2d(i,j)=lat(i,j)-y0
          end do
          end do

!$omp end do

        end if

! -----

!! -----

! Calculate the x and the y coordinates from the latitude and the
! longitude with the Polar Stereographic projection method.

      else if(mpopt.eq.1) then

!$omp do schedule(runtime) private(i,j,rr,dlon)

        do j=jstr,jend
        do i=istr,iend
          rr=cpj(2)*cos(lat(i,j)*d2r)/((1.e0+sin(pold2r*lat(i,j)))+eps)

          dlon=(lon(i,j)-tlon)*d2r

          x2d(i,j)=rr*sin(dlon)-x0
          y2d(i,j)=-rpol*(rr*cos(dlon)+y0)

        end do
        end do

!$omp end do

! -----

! Calculate the x and the y coordinates from the latitude and the
! longitude with the Lambert Conformal Conic projection method.

      else if(mpopt.eq.2) then

!$omp do schedule(runtime) private(i,j,rr,dlon)

        do j=jstr,jend
        do i=istr,iend
          rr=cpj(6)*exp(cpj(4)                                          &
     &      *log(tan((45.e0-pol05*lat(i,j))*d2r)*cpj(3)+eps))

          dlon=d2r4*(lon(i,j)-tlon)

          x2d(i,j)=rr*sin(dlon)-x0
          y2d(i,j)=-rpol*(rr*cos(dlon)+y0)

        end do
        end do

!$omp end do

! -----

! Calculate the x and the y coordinates from the latitude and the
! longitude with the Mercator projection method.

      else if(mpopt.eq.3.or.mpopt.eq.13) then

!$omp do schedule(runtime) private(i,j)

        do j=jstr,jend
        do i=istr,iend
          x2d(i,j)=d2r2*lon(i,j)-x0
          y2d(i,j)=-rpol2*log(tan((45.e0-pol05*lat(i,j))*d2r)+eps)-y0
        end do
        end do

!$omp end do

! -----

! Calculate the x and the y coordinates from the latitude and the
! longitude without any projection method.

      else if(mpopt.eq.4) then

!$omp do schedule(runtime) private(i,j,dlon)

        do j=jstr,jend
        do i=istr,iend
          dlon=lon(i,j)-tlon

          if(dlon.lt.-180.e0) then
            dlon=dlon+360.e0
          end if

          if(dlon.gt.180.e0) then
            dlon=dlon-360.e0
          end if

          x2d(i,j)=cpj(1)*dlon*cos(lat(i,j)*d2r)-x0
          y2d(i,j)=cpj(1)*lat(i,j)-y0

        end do
        end do

!$omp end do

! -----

! Calculate the x and the y coordinates from the latitude and the
! longitude with the circular cylinder coordinates system.

      else if(mpopt.eq.5) then

!$omp do schedule(runtime) private(i,j)

        do j=jstr,jend
        do i=istr,iend
          x2d(i,j)=cpj(1)*lon(i,j)*cos(lat(i,j)*d2r)-x0
          y2d(i,j)=cpj(1)*lat(i,j)-y0
        end do
        end do

!$omp end do

      end if

! -----

!$omp end parallel

!!! -----

      end subroutine s_ll2xy

!-----7--------------------------------------------------------------7--

      end module m_ll2xy
