!***********************************************************************
      module m_trndamp
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/07/15
!     Modification: 2002/08/27, 2002/09/02, 2002/10/31, 2003/05/19,
!                   2003/10/10, 2003/12/12, 2004/05/07, 2004/05/31,
!                   2005/02/10, 2005/04/04, 2006/11/06, 2007/01/20,
!                   2008/05/02, 2008/08/19, 2008/08/25, 2008/10/10,
!                   2008/12/11, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the damped terrain height.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_cpondtrn
      use m_getrij
      use m_gettrn
      use m_getxy
      use m_hint2d
      use m_outtrn
      use m_rdsfcdmp

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: trndamp, s_trndamp

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface trndamp

        module procedure s_trndamp

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
      subroutine s_trndamp(x0,y0,cpj,x0gpv,y0gpv,cpjgpv,ni,nj,land,ht,  &
     &                     tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,nid,njd,htdat)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

      real, intent(in) :: x0
                       ! x origin of model grid

      real, intent(in) :: y0
                       ! y origin of model grid

      real, intent(in) :: cpj(1:7)
                       ! Map projection parameters of model grid

      real, intent(in) :: x0gpv
                       ! x origin of GPV data grid

      real, intent(in) :: y0gpv
                       ! y origin of GPV data grid

      real, intent(in) :: cpjgpv(1:7)
                       ! Map projection parameters of GPV data grid

      real, intent(in) :: htdat(1:nid,1:njd)
                       ! Terrain height in data

! Internal shared variables

      integer, intent(inout) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      real, intent(inout) :: ht(0:ni+1,0:nj+1)
                       ! Terrain height

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp5(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp6(0:ni+1,0:nj+1)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

!! Set the terrain height.

! Get the x and the y coordinates at the model grid points.

      call s_getxy(iddx,iddy,'xx',0,ni+1,0,nj+1,tmp1,tmp2)

! -----

! Get the terrain height.

      call s_gettrn(idtrnopt,idzsfc,idmnthgh,idmntwx,idmntwy,           &
     &              idmntcx,idmntcy,'terrain     ',7,1,ni,nj,           &
     &              tmp1,tmp2,ht)

! -----

! Read out the land use categories from the interpolated surface file.

      call rdsfcdmp(idexprim,idcrsdir,idsfcdat,idncexp,idnccrs,         &
     &              idwlngth,idsfcopt,idzsfc,ni,nj,ht,land)

! -----

! Correspond and damping the model height to the external data height
! and read in to the reestimated terrain height.

      call getrij(idmpopt_gpv,idnspol_gpv,idtlon_gpv,                   &
     &            iddxiv_gpv,iddyiv_gpv,'gridata ',7,'xx',              &
     &            x0,y0,cpj,x0gpv,y0gpv,cpjgpv,ni,nj,                   &
     &            tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)

      call hint2d(idmpopt_gpv,idintopt_gpv,'xx','cal ',ni,nj,tmp1,tmp2, &
     &            tmp3,tmp4,tmp5,nid,njd,htdat)

      call s_cpondtrn(idexbwid,ni,nj,land,tmp5,ht,tmp1,tmp2,tmp3)

      call outtrn(idexprim,idcrsdir,idncexp,idnccrs,idwlngth,           &
     &            'terrain.damp',12,1,ni,nj,ht)

! -----

!! -----

      end subroutine s_trndamp

!-----7--------------------------------------------------------------7--

      end module m_trndamp
