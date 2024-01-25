!***********************************************************************
      module m_getnews
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/07/01
!     Modification: 2004/08/01, 2005/02/10, 2005/04/04, 2006/11/06,
!                   2007/01/20, 2008/05/02, 2008/08/25, 2008/12/11,
!                   2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the mininum and maximum data indices to create the base state
!     variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_commpi
      use m_getrij
      use m_newsindx

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getnews, s_getnews

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getnews

        module procedure s_getnews

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
      subroutine s_getnews(nid,njd,idstr,idend,jdstr,jdend,             &
     &                     x0,y0,cpj,x0gpv,y0gpv,cpjgpv,ni,nj,          &
     &                     ri,rj,tmp1,tmp2,tmp3,tmp4)
!***********************************************************************

! Input variables

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

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

      real, intent(in) :: x0gpv
                       ! x origin of GPV data grid

      real, intent(in) :: y0gpv
                       ! y origin of GPV data grid

      real, intent(in) :: cpjgpv(1:7)
                       ! Map projection parameters of GPV data grid

! Input and output variables

      integer, intent(inout) :: idstr
                       ! Minimum index of model grid in data region
                       ! in x direction

      integer, intent(inout) :: idend
                       ! Maximum index of model grid in data region
                       ! in x direction

      integer, intent(inout) :: jdstr
                       ! Minimum index of model grid in data region
                       ! in y direction

      integer, intent(inout) :: jdend
                       ! Maximum index of model grid in data region
                       ! in y direction

! Internal shared variables

      real, intent(inout) :: ri(0:ni+1,0:nj+1)
                       ! Real indices in data region in x direction

      real, intent(inout) :: rj(0:ni+1,0:nj+1)
                       ! Real indices in data region in y direction

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! Initialize the output variables.

      if(mype.eq.root) then

        idstr=2*nid
        idend=-nid

        jdstr=njd
        jdend=1

      end if

! -----

! Calculate the real indices at the model grid points in the data
! region.

      call getrij(idmpopt_gpv,idnspol_gpv,idtlon_gpv,                   &
     &            iddxiv_gpv,iddyiv_gpv,'gridata ',7,'xx',              &
     &            x0,y0,cpj,x0gpv,y0gpv,cpjgpv,ni,nj,ri,rj,             &
     &            tmp1,tmp2,tmp3,tmp4)

! -----

! Get the mininum and maximum data indices to create the base state
! variables.

      call newsindx(idmpopt,idmpopt_gpv,                                &
     &              nid,njd,idstr,idend,jdstr,jdend,ni,nj,ri,rj)

! -----

      end subroutine s_getnews

!-----7--------------------------------------------------------------7--

      end module m_getnews
