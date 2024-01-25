!***********************************************************************
      module m_getrij
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/04/06
!     Modification: 1999/05/10, 1999/05/20, 1999/10/12, 1999/11/01,
!                   2000/01/17, 2000/12/18, 2002/06/18, 2002/07/03,
!                   2003/04/30, 2003/05/19, 2003/12/12, 2006/09/21,
!                   2006/11/06, 2008/05/02, 2008/08/25, 2008/10/10,
!                   2009/02/27, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the real indices at the model grid points in the data
!     region.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_getxy
      use m_ll2xy
      use m_xy2ij
      use m_xy2ll

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getrij, s_getrij

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getrij

        module procedure s_getrij

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!     none

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_getrij(fpmpopt_dat,fpnspol_dat,fptlon_dat,           &
     &                    fpdxiv_dat,fpdyiv_dat,pname,ncpn,             &
     &                    xo,x0,y0,cpj,x0dat,y0dat,cpjdat,              &
     &                    ni,nj,ri,rj,lat,lon,x,y)
!***********************************************************************

! Input variables

      character(len=8), intent(in) :: pname
                       ! Running program name

      character(len=2), intent(in) :: xo
                       ! Control flag of variable arrangement

      integer, intent(in) :: fpmpopt_dat
                       ! Formal parameter of unique index of mpopt_XXX
                       ! in namelist table

      integer, intent(in) :: fpnspol_dat
                       ! Formal parameter of unique index of nspol_XXX
                       ! in namelist table

      integer, intent(in) :: fptlon_dat
                       ! Formal parameter of unique index of tlon_XXX
                       ! in namelist table

      integer, intent(in) :: fpdxiv_dat
                       ! Formal parameter of unique index of dxiv_XXX
                       ! in namelist table

      integer, intent(in) :: fpdyiv_dat
                       ! Formal parameter of unique index of dyiv_XXX
                       ! in namelist table

      integer, intent(in) :: ncpn
                       ! Number of character of pname

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      real, intent(in) :: x0
                       ! x origin of model grid

      real, intent(in) :: y0
                       ! y origin of model grid

      real, intent(in) :: cpj(1:7)
                       ! Map projection parameters of model grid

      real, intent(in) :: x0dat
                       ! x origin of data grid

      real, intent(in) :: y0dat
                       ! y origin of data grid

      real, intent(in) :: cpjdat(1:7)
                       ! Map projection parameters of data grid

! Output variables

      real, intent(out) :: ri(0:ni+1,0:nj+1)
                       ! Real indices in data region in x direction

      real, intent(out) :: rj(0:ni+1,0:nj+1)
                       ! Real indices in data region in y direction

! Internal shared variables

      real, intent(inout) :: lat(0:ni+1,0:nj+1)
                       ! Latitude at model grid points

      real, intent(inout) :: lon(0:ni+1,0:nj+1)
                       ! Longitude at model grid points

      real, intent(inout) :: x(0:ni+1,0:nj+1)
                       ! x coordinates at model grid points

      real, intent(inout) :: y(0:ni+1,0:nj+1)
                       ! y coordinates at model grid points

!-----7--------------------------------------------------------------7--

! Calculate the x and the y coodinates at the model grid points.

      call s_getxy(iddx,iddy,xo,0,ni+1,0,nj+1,x,y)

! -----

! Get the latitude and the longitude from the x and the y coodinates at
! the model grid points.

      call s_xy2ll(idmpopt,idnspol,idtlon,'solver  ',6,xo,x0,y0,cpj,    &
     &             0,ni+1,0,nj+1,x,y,lat,lon)

! -----

! Calculate the x and the y coodinates at the model grid points from the
! latitude and the longitude on the data grid.

      call ll2xy(fpmpopt_dat,fpnspol_dat,fptlon_dat,pname,ncpn,xo,      &
     &           x0dat,y0dat,cpjdat,0,ni+1,0,nj+1,lat,lon,x,y)

! -----

! Calculate the real indices from the x and the y coodinates at the
! model grid points in the data region.

      call xy2ij(fpdxiv_dat,fpdyiv_dat,xo,0,ni+1,0,nj+1,x,y,ri,rj)

! -----

      end subroutine s_getrij

!-----7--------------------------------------------------------------7--

      end module m_getrij
